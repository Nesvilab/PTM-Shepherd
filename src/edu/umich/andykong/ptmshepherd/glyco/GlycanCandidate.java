/*
 *    Copyright 2022 University of Michigan
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

package edu.umich.andykong.ptmshepherd.glyco;

import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanParser;
import umich.ms.glyco.GlycanResidue;

import java.util.*;

/**
 Container for theoretical glycan compositions (supplied by user or database) to be
 searched against experimental results.
 */
public class GlycanCandidate extends Glycan {
    Map<GlycanResidue, Integer> glycanComposition;     // map of residue type: count of residue to describe the composition
    double decoyMassShift;
    boolean isDecoy;
    public TreeMap<String, GlycanFragment> Yfragments;
    public TreeMap<String, GlycanFragment> oxoniumFragments;
    public String name;

    /**
     * Construct a search candidate from a glycan composition and fragment ion lists
     */
    public GlycanCandidate(Map<GlycanResidue, Integer> glycanComposition,
                           double decoyMassShift,
                           boolean isDecoy,
                           HashMap<String, GlycanResidue> glycanResiduesMap,
                           TreeMap<String, GlycanFragment> Yfragments,
                           TreeMap<String, GlycanFragment> oxoniumFragments) {
        super(glycanComposition);
        this.decoyMassShift = decoyMassShift;
        if (isDecoy) {
            mass += decoyMassShift;
        }
        // make sure that all residue types are accounted for (add Residue with 0 counts for any not included in the file)
        for (GlycanResidue residue : glycanResiduesMap.values()){
            if (!this.glycanComposition.containsKey(residue)) {
                this.glycanComposition.put(residue, 0);
            }
        }
        this.isDecoy = isDecoy;
        this.Yfragments = Yfragments;
        this.oxoniumFragments = oxoniumFragments;
        name = this.toString();
    }

    /**
     * Initialize a GlycanCandidate and its associated GlycanFragment ions from an input Glycan composition.
     */
    public static GlycanCandidate initGlycanCandidate(Map<GlycanResidue, Integer> inputGlycanComp, double decoyMassShift, boolean isDecoy,
                           HashMap<String, GlycanResidue> glycanResiduesMap, boolean nGlycan, Random randomGenerator,
                           HashMap<GlycanResidue, ArrayList<GlycanFragment>> glycoOxoniumDatabase) {
        TreeMap<String, GlycanFragment> Yfragments = initializeYFragments(inputGlycanComp, isDecoy, glycanResiduesMap, nGlycan, randomGenerator);
        TreeMap<String, GlycanFragment> oxoniumFragments = initializeOxoniumFragments(inputGlycanComp, isDecoy, glycoOxoniumDatabase);

        return new GlycanCandidate(inputGlycanComp, decoyMassShift, isDecoy, glycanResiduesMap, Yfragments, oxoniumFragments);
    }

    /**
     * Constructor for a new glycan candidate for a second search using fragment propensities, initialized from a
     * candidate from the first search. Making a new container and fragments to avoid threading issues, but passing
     * original probabilities where propensities not found.
     * @param oldCandidate original search candidate to use as model
     */
    public static GlycanCandidate initCandidateFromProps(GlycanCandidate oldCandidate,
                                                         HashMap<String, GlycanResidue> glycanResiduesMap,
                                                         HashMap<String, Double> yFragmentPropensities,
                                                         HashMap<String, Double> yFragmentIntensities,
                                                         HashMap<String, Double> oxFragmentPropensities,
                                                         HashMap<String, Double> oxFragmentIntensities) {
        // initialize fragments for this candidate
        TreeMap<String, GlycanFragment> Yfragments = initializeYFragmentsFromProps(oldCandidate.Yfragments, yFragmentPropensities, yFragmentIntensities);
        TreeMap<String, GlycanFragment> oxoniumFragments = initializeOxoniumFragmentsFromProps(oldCandidate.oxoniumFragments, oxFragmentPropensities, oxFragmentIntensities);

        return new GlycanCandidate(oldCandidate.glycanComposition, oldCandidate.decoyMassShift, oldCandidate.isDecoy, glycanResiduesMap, Yfragments, oxoniumFragments);

    }

    /**
     * Constructor for copying existing glycan candidate to new object to avoid concurrent access in multi-threading.
     * Take all information from previous candidate, just initialize as a new object. Also re-initialize Fragment
     * objects for same reason.
     * @param otherCandidate candidate to copy
     */
    public static GlycanCandidate copyCandidate(GlycanCandidate otherCandidate, HashMap<String, GlycanResidue> glycanResiduesMap) {
        TreeMap<String, GlycanFragment> Yfragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> entry : otherCandidate.Yfragments.entrySet()) {
            Yfragments.put(entry.getKey(), GlycanFragment.copyFragment(entry.getValue()));
        }
        TreeMap<String, GlycanFragment> oxoniumFragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> entry : otherCandidate.oxoniumFragments.entrySet()) {
            oxoniumFragments.put(entry.getKey(), GlycanFragment.copyFragment(entry.getValue()));
        }
        return new GlycanCandidate(otherCandidate.glycanComposition, otherCandidate.decoyMassShift, otherCandidate.isDecoy, glycanResiduesMap, Yfragments, oxoniumFragments);
    }

    /**
     * Empty candidate
     */
    public static GlycanCandidate emptyCandidate() {
        return new GlycanCandidate(new TreeMap<>(), 0,false, new HashMap<>(), new TreeMap<>(), new TreeMap<>());
    }

    /**
     * Initialize array of all fragment ions to search for this candidate using the original candidate's fragments as a template.
     * Decoy fragments are generated for decoy candidates.
     * Fragment propensities are used for all fragments found in the bootstrap/input data (specified in the input map)
     * and fragments lacking any input are assumed to have 0 input propensity.
     */
    public static TreeMap<String, GlycanFragment> initializeYFragmentsFromProps(TreeMap<String, GlycanFragment> originalYs, HashMap<String, Double> fragmentPropensities, HashMap<String, Double> fragmentIntensities) {
        // Initialize a new Y fragment for each in the input map, adding propensity/intensity from the fragmentInfo container
        TreeMap<String, GlycanFragment> Yfragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> originalFragEntry : originalYs.entrySet()) {
            double expectedIntensity;
            double propensity;
            GlycanFragment origFrag = originalFragEntry.getValue();
            if (fragmentPropensities.containsKey(originalFragEntry.getKey())) {
                // have propensity/intensity info for this fragment - read from input fragmentInfo
                expectedIntensity = fragmentIntensities.get(origFrag.hash);
                propensity = fragmentPropensities.get(origFrag.hash);
            } else {
                // no added info - copy the original
                expectedIntensity = origFrag.expectedIntensity;
                propensity = origFrag.propensity;
            }
            GlycanFragment newFragment = GlycanFragment.copyFragmentWithPropensity(origFrag, expectedIntensity, propensity);
            Yfragments.put(originalFragEntry.getKey(), newFragment);
        }
        return Yfragments;
    }

    /**
     * Init new oxonium ions based on the original candidate's ions, updating expected intensity/propensity if found
     * in the provided fragmentInfo container. Same logic as for Y ions
     * @param originalOxos original candidate's oxonium fragment map
     */
    public static TreeMap<String, GlycanFragment> initializeOxoniumFragmentsFromProps(TreeMap<String, GlycanFragment> originalOxos, HashMap<String, Double> fragmentPropensities, HashMap<String, Double> fragmentIntensities) {
        // Initialize a new oxonium fragment for each in the input map, adding propensity/intensity from the fragmentInfo container
        TreeMap<String, GlycanFragment> oxoniumFragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> originalFragEntry : originalOxos.entrySet()) {
            double expectedIntensity;
            double propensity;
            GlycanFragment origFrag = originalFragEntry.getValue();
            if (fragmentPropensities.containsKey(originalFragEntry.getKey())) {
                // have propensity/intensity info for this fragment - read from input fragmentInfo
                expectedIntensity = fragmentIntensities.get(origFrag.hash);
                propensity = fragmentPropensities.get(origFrag.hash);
            } else {
                // no added info - copy the original
                expectedIntensity = origFrag.expectedIntensity;
                propensity = origFrag.propensity;
            }
            GlycanFragment newFragment = GlycanFragment.copyFragmentWithPropensity(origFrag, expectedIntensity, propensity);
            oxoniumFragments.put(originalFragEntry.getKey(), newFragment);
        }
        return oxoniumFragments;
    }

    /**
     * Initialize array of all fragment ions to search for this candidate. Decoy fragments are generated for decoy candidates.
     */
    public static TreeMap<String, GlycanFragment> initializeYFragments(Map<GlycanResidue, Integer> glycanComposition,
                                                                       boolean isDecoy,
                                                                       HashMap<String, GlycanResidue> glycanResiduesMap,
                                                                       boolean nGlycan,
                                                                       Random randomGenerator) {
        TreeMap<String, GlycanFragment> Yfragments = new TreeMap<>();
        GlycanResidue hexnac = GlycanParser.findResidueName("HexNAc", glycanResiduesMap);
        List<Map.Entry<GlycanResidue, Integer>> remainingComp = new ArrayList<>();
        for (Map.Entry<GlycanResidue, Integer> compEntry : glycanComposition.entrySet()) {
            // Do not include labile residues in Y ions
            if (compEntry.getValue() > 0 && !compEntry.getKey().isLabile) {
                remainingComp.add(compEntry);
            }
        }

        ArrayList<TreeMap<GlycanResidue, Integer>> previousYs = new ArrayList<>();
        while (!remainingComp.isEmpty()) {
            Map.Entry<GlycanResidue, Integer> currentEntry = remainingComp.remove(0);

            // add this residue to all previous Y ions (while also keeping previous Ys without it)
            ArrayList<TreeMap<GlycanResidue, Integer>> newYs = new ArrayList<>();
            for (TreeMap<GlycanResidue, Integer> previousY: previousYs) {
                for (int i=1; i <= currentEntry.getValue(); i++) {
                    TreeMap<GlycanResidue, Integer> newY = (TreeMap<GlycanResidue, Integer>) previousY.clone();
                    newY.put(currentEntry.getKey(), i);
                    newYs.add(newY);
                }
            }
            previousYs.addAll(newYs);

            // add entries for this residue only to the list of Y ions
            for (int i=1; i <= currentEntry.getValue(); i++) {
                TreeMap<GlycanResidue, Integer> currentY = new TreeMap<>();
                currentY.put(currentEntry.getKey(), i);
                previousYs.add(currentY);
            }
        }

        for (TreeMap<GlycanResidue, Integer> Ycomp: previousYs) {
            // require hexnac for N-glycan Ys. todo: replace hard-coded rule with user params
            if (nGlycan && glycanComposition.getOrDefault(hexnac, 0) > 0 && !Ycomp.containsKey(hexnac)) {
                continue;
            }
            if (nGlycan && glycanComposition.getOrDefault(hexnac, 0) > 0) {
                if (Ycomp.containsKey(GlycanParser.findResidueName("Hex", glycanResiduesMap)) && Ycomp.get(hexnac) < 2) {
                    continue;
                }
            }
            GlycanFragment newFragment = GlycanFragment.initializeYFragment(Ycomp, isDecoy, randomGenerator);
            Yfragments.put(newFragment.hash, newFragment);
        }
        return Yfragments;
    }

    /**
     * Helper method to initialize oxonium fragment rules. Only initializes fragments for a
     * residue type if at least one candidate contains that residue type (no need to consider if not).
     * Decoys generated for all residue types that have at least one decoy candidate containing that type.
     */
    public static TreeMap<String, GlycanFragment> initializeOxoniumFragments(Map<GlycanResidue, Integer> glycanComposition, boolean isDecoy, HashMap<GlycanResidue, ArrayList<GlycanFragment>> glycoOxoniumDatabase) {
        TreeMap<String, GlycanFragment> oxoniumFragments = new TreeMap<>();
        for (GlycanResidue residue : glycanComposition.keySet()) {
            if (glycanComposition.get(residue) > 0) {
                if (isDecoy) {
                    oxoniumFragments.putAll(makeOxoniums(residue, true, glycoOxoniumDatabase));
                } else {
                    oxoniumFragments.putAll(makeOxoniums(residue, false, glycoOxoniumDatabase));
                }
            }
        }
        return oxoniumFragments;
    }

    /**
     * Helper method to add fragment ions to the oxonium list
     * @param residue residue type
     * @param isDecoy decoy or not
     * @return updated list
     */
    private static TreeMap<String, GlycanFragment> makeOxoniums(GlycanResidue residue, boolean isDecoy, HashMap<GlycanResidue, ArrayList<GlycanFragment>> glycoOxoniumDatabase) {
        TreeMap<String, GlycanFragment> newFragments = new TreeMap<>();
        ArrayList<GlycanFragment> oxoniumIonDescriptors = glycoOxoniumDatabase.getOrDefault(residue, new ArrayList<>());
        for (GlycanFragment fragmentDescriptor : oxoniumIonDescriptors) {
            GlycanFragment newFragment = GlycanFragment.copyFragment(fragmentDescriptor);
            newFragment.isDecoy = isDecoy;
            newFragments.put(newFragment.hash, newFragment);
        }
        return newFragments;
    }

    /**
     * Return string representation of this glycan for writing to output tables
     * @return string
     */
    public String toString() {
        return toCandidateString(glycanComposition, mass, isDecoy);
    }

    // Format: (Decoy_)<GlycanComposition> % <Mass>
    public static String toCandidateString(Map<GlycanResidue, Integer> composition, double mass, boolean isDecoy) {
        String decoy = isDecoy ? "Decoy_" : "";
        return String.format("%s%s %% %.4f", decoy, Glycan.toGlycanString(composition), mass);
    }
}
