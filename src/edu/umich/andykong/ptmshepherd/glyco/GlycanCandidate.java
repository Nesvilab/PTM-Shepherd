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

import edu.umich.andykong.ptmshepherd.core.AAMasses;

import java.util.*;

/**
 Container for theoretical glycan compositions (supplied by user or database) to be
 searched against experimental results.
 */
public class GlycanCandidate {
    double monoisotopicMass;
    Map<GlycanResidue, Integer> glycanComposition;     // map of residue type: count of residue to describe the composition
    boolean isDecoy;
    public static final double MAX_CANDIDATE_DECOY_SHIFT_DA = 3;
    public static final double DEFAULT_PEPTIDE_MASS = 1500;
    public TreeMap<String, GlycanFragment> Yfragments;
    public TreeMap<String, GlycanFragment> oxoniumFragments;
    public boolean hasFragmentProps;    // if this candidate has fragment propensity info or default values

    /**
     * Constructor for reading fragment probability information from the glycoFrags file
     * @param glycanStr glycan composition string to read
     * @param parsedFragmentInfo list of strings containing fragment ion info
     */
    public GlycanCandidate(String glycanStr, String[] parsedFragmentInfo, GlycoParams glycoParams){
        glycanComposition = glycoParams.parseGlycanString(glycanStr);
        monoisotopicMass = computeMonoisotopicMass(glycanComposition);

        Yfragments = new TreeMap<>();
        oxoniumFragments = new TreeMap<>();
        for (String fragment : parsedFragmentInfo) {
            String[] typeSplits = fragment.split("~");
            // string format is [type]~[composition]~[intensity]
            if (typeSplits[0].matches("Y")) {
                // Y ion
                GlycanFragment newFragment = new GlycanFragment(typeSplits[1], Double.parseDouble(typeSplits[2]), GlycanFragment.FragType.Y, glycoParams);
                Yfragments.put(newFragment.hash, newFragment);
            } else if (typeSplits[0].matches("Ox")) {
                GlycanFragment newFragment = new GlycanFragment(typeSplits[1], Double.parseDouble(typeSplits[2]), GlycanFragment.FragType.Ox, glycoParams);
                oxoniumFragments.put(newFragment.hash, newFragment);
            } else {
                // invalid
            }
        }
    }

    /**
     * Base constructor for a new glycan candidate for initial search (not using fragment propensities)
     * @param inputGlycanComp composition map
     * @param isDecoy bool
     * @param glycoParams parameters container
     */
    public GlycanCandidate(Map<GlycanResidue, Integer> inputGlycanComp, boolean isDecoy, GlycoParams glycoParams) {
        this.glycanComposition = inputGlycanComp;
        this.isDecoy = isDecoy;
        // make sure that all residue types are accounted for (add Residue with 0 counts for any not included in the file)
        for (GlycanResidue residue : glycoParams.glycanResidues){
            if (!this.glycanComposition.containsKey(residue)) {
                this.glycanComposition.put(residue, 0);
            }
        }
        this.monoisotopicMass = setMassHelper(glycoParams);

        // initialize fragments for this candidate
        initializeYFragments(glycoParams.randomGenerator);
        initializeOxoniumFragments(glycoParams.glycoOxoniumDatabase, glycoParams.randomGenerator);
        this.hasFragmentProps = false;
    }

    /**
     * Constructor for a new glycan candidate for a second search using fragment propensities, initialized from a
     * candidate from the first search. Making a new container and fragments to avoid threading issues, but passing
     * original probabilities where propensities not found.
     * @param oldCandidate original search candidate to use as model
     * @param fragmentInfo input propensities and intensities for fragment ions
     * @param glycoParams parameters container
     */
    public GlycanCandidate(GlycanCandidate oldCandidate, GlycanCandidateFragments fragmentInfo, GlycoParams glycoParams) {
        this.glycanComposition = oldCandidate.glycanComposition;
        this.isDecoy = oldCandidate.isDecoy;
        // make sure that all residue types are accounted for (add Residue with 0 counts for any not included in the file)
        for (GlycanResidue residue : glycoParams.glycanResidues){
            if (!this.glycanComposition.containsKey(residue)) {
                this.glycanComposition.put(residue, 0);
            }
        }
        this.monoisotopicMass = setMassHelper(glycoParams);

        // initialize fragments for this candidate
        initializeYFragmentsFromProps(oldCandidate.Yfragments, fragmentInfo, glycoParams.randomGenerator);
        initializeOxoniumFragmentsFromProps(oldCandidate.oxoniumFragments, fragmentInfo, glycoParams.randomGenerator);
        this.hasFragmentProps = true;
    }

    /**
     * Constructor for copying existing glycan candidate to new object to avoid concurrent access in multi-threading.
     * Take all information from previous candidate, just initialize as a new object. Also re-initialize Fragment
     * objects for same reason.
     * @param inputGlycanComp composition
     * @param isDecoy decoy bool
     * @param monoisotopicMass intact mass
     * @param yfragments list of Y fragments
     * @param oxoniumFragments list of oxonium fragments
     */
    public GlycanCandidate(Map<GlycanResidue, Integer> inputGlycanComp, boolean isDecoy, double monoisotopicMass, TreeMap<String, GlycanFragment> yfragments, TreeMap<String, GlycanFragment> oxoniumFragments) {
        this.glycanComposition = inputGlycanComp;
        this.isDecoy = isDecoy;
        this.monoisotopicMass = monoisotopicMass;
        this.Yfragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> entry : yfragments.entrySet()) {
            this.Yfragments.put(entry.getKey(), new GlycanFragment(entry.getValue()));
        }
        this.oxoniumFragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> entry : oxoniumFragments.entrySet()) {
            this.oxoniumFragments.put(entry.getKey(), new GlycanFragment(entry.getValue()));
        }
        this.hasFragmentProps = false;
    }

    /**
     * Empty candidate to avoid null pointers
     */
    public GlycanCandidate() {
        this.Yfragments = new TreeMap<>();
        this.oxoniumFragments = new TreeMap<>();
        this.glycanComposition = new TreeMap<>();
        this.isDecoy = false;
        this.monoisotopicMass = 0;
        this.hasFragmentProps = false;
    }

    // Helper method for determining decoy masses for various decoy mass generation settings
    private double setMassHelper(GlycoParams glycoParams) {
        double mass;
        if (! isDecoy) {
            mass = computeMonoisotopicMass(glycanComposition);
        } else {
            double baseMonoistopicMass = computeMonoisotopicMass(glycanComposition);
            double randomShift = 0;
            switch (glycoParams.decoyType) {
                case 0:
                    // simple mass window
                    randomShift = GlycanFragment.randomMassShift(MAX_CANDIDATE_DECOY_SHIFT_DA, glycoParams.randomGenerator);
                    break;
                case 1:
                    // random isotope and mass error
                    randomShift = getRandomShiftIsotopes(baseMonoistopicMass, glycoParams.glycoIsotopes, glycoParams.glycoPPMtol, glycoParams.randomGenerator);
                    break;
                case 2:
                    // random mass error, no isotope error
                    Integer[] noIsotopes = {0};
                    randomShift = getRandomShiftIsotopes(baseMonoistopicMass, noIsotopes, glycoParams.glycoPPMtol, glycoParams.randomGenerator);
                case 3:
                    // exact target mass - random shift left at 0
                    break;
            }
            mass = baseMonoistopicMass + randomShift;
        }
        return mass;
    }

    /**
     * Initialize array of all fragment ions to search for this candidate using the original candidate's fragments as a template.
     * Decoy fragments are generated for decoy candidates.
     * Fragment propensities are used for all fragments found in the bootstrap/input data (specified in the input map)
     * and fragments lacking any input are assumed to have 0 input propensity.
     */
    public void initializeYFragmentsFromProps(TreeMap<String, GlycanFragment> originalYs, GlycanCandidateFragments fragmentInfo, Random randomGenerator) {
        // Initialize a new Y fragment for each in the input map, adding propensity/intensity from the fragmentInfo container
        this.Yfragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> originalFragEntry : originalYs.entrySet()) {
            double expectedIntensity;
            double propensity;
            GlycanFragment origFrag = originalFragEntry.getValue();
            if (fragmentInfo.yFragmentProps.containsKey(originalFragEntry.getKey())) {
                // have propensity/intensity info for this fragment - read from input fragmentInfo
                expectedIntensity = fragmentInfo.yFragmentIntensities.get(origFrag.hash);
                propensity = fragmentInfo.yFragmentProps.get(origFrag.hash);
            } else {
                // no added info - copy the original
                expectedIntensity = origFrag.expectedIntensity;
                propensity = origFrag.propensity;
            }
            GlycanFragment newFragment = new GlycanFragment(origFrag, expectedIntensity, propensity);
            this.Yfragments.put(originalFragEntry.getKey(), newFragment);
        }
    }

    /**
     * Init new oxonium ions based on the original candidate's ions, updating expected intensity/propensity if found
     * in the provided fragmentInfo container. Same logic as for Y ions
     * @param originalOxos original candidate's oxonium fragment map
     * @param fragmentInfo fragmt info container
     * @param randomGenerator the run's random generator
     */
    public void initializeOxoniumFragmentsFromProps(TreeMap<String, GlycanFragment> originalOxos, GlycanCandidateFragments fragmentInfo, Random randomGenerator) {
        // Initialize a new oxonium fragment for each in the input map, adding propensity/intensity from the fragmentInfo container
        this.oxoniumFragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> originalFragEntry : originalOxos.entrySet()) {
            double expectedIntensity;
            double propensity;
            GlycanFragment origFrag = originalFragEntry.getValue();
            if (fragmentInfo.OxFragmentProps.containsKey(originalFragEntry.getKey())) {
                // have propensity/intensity info for this fragment - read from input fragmentInfo
                expectedIntensity = fragmentInfo.OxFragmentIntensities.get(origFrag.hash);
                propensity = fragmentInfo.OxFragmentProps.get(origFrag.hash);
            } else {
                // no added info - copy the original
                expectedIntensity = origFrag.expectedIntensity;
                propensity = origFrag.propensity;
            }
            GlycanFragment newFragment = new GlycanFragment(origFrag, expectedIntensity, propensity);
            this.oxoniumFragments.put(originalFragEntry.getKey(), newFragment);
        }
    }

    private static <T> List<List<T>> generateCombinations(List<T> source, List<T> comb, int targetSize) {
        if (comb.size() == targetSize) {
            List<List<T>> result = new ArrayList<>();
            result.add(comb);
            return result;
        }

        List<List<T>> result = new ArrayList<>();
        Iterator<T> iterator = source.iterator();
        while (iterator.hasNext()) {
            T item = iterator.next();
            iterator.remove();

            // creating a new combination using existing as a base
            List<T> newComb = new ArrayList<T>(comb);
            newComb.add(item); // adding the element that was removed from the source
            result.addAll(generateCombinations(new ArrayList<>(source), newComb, targetSize)); // adding all the combinations generated
        }
        return result;
    }

    private String hashComp(TreeMap<GlycanResidue, Integer> composition) {
        StringBuilder output = new StringBuilder();
        for (Map.Entry<GlycanResidue, Integer> entry: composition.entrySet()) {
            output.append(entry.getKey().name);
            output.append(entry.getValue());
        }
        return output.toString();
    }

    /**
     * Initialize array of all fragment ions to search for this candidate. Candidate has
     * fragment Ys up to the max HexNAc, Hex, and dHex present. Decoy fragments are generated for decoy candidates.
     */
    public void initializeYFragments(Random randomGenerator) {
        // Initialize list of all Y fragments to consider. Currently using only HexNAc, Hex, and dHex in Y ions
        Yfragments = new TreeMap<>();
        List<GlycanResidue> allGlycans = new ArrayList<>();
        for (Map.Entry<GlycanResidue, Integer> compEntry: glycanComposition.entrySet()) {
            for (int i=0; i < compEntry.getValue(); i++) {
                allGlycans.add(compEntry.getKey());
            }
        }

        // compute all combinations of length k
        HashSet<String> foundCombos = new HashSet<>();
        for (int k=0; k <= allGlycans.size(); k++) {
            ArrayList<GlycanResidue> inputGlycans = new ArrayList<>(allGlycans);
            List<List<GlycanResidue>> combos = generateCombinations(inputGlycans, new ArrayList<>(), k);
            for (List<GlycanResidue> combo: combos) {
                if (combo.size() > 0) {
                    TreeMap<GlycanResidue, Integer> composition = new TreeMap<>();
                    for (GlycanResidue residue : combo) {
                        composition.put(residue, composition.getOrDefault(residue, 0) + 1);
                    }
                    String hash = hashComp(composition);
                    if (!foundCombos.contains(hash)) {
                        foundCombos.add(hash);
                        GlycanFragment fragment = new GlycanFragment(composition, this.isDecoy, randomGenerator, GlycanFragment.FragType.Y);
                        Yfragments.put(hash, fragment);
                    }
                }
            }
        }

        // todo: add filtering rules

    }

    /**
     * Helper method to initialize hard-coded oxonium fragment rules. Only initializes fragments for a
     * residue type if at least one candidate contains that residue type (no need to consider if not).
     * Decoys generated for all residue types that have at least one decoy candidate containing that type.
     */
    public void initializeOxoniumFragments(HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase, Random randomGenerator) {
        this.oxoniumFragments = new TreeMap<>();
        for (GlycanResidue residue : glycanComposition.keySet()) {
            if (this.glycanComposition.get(residue) > 0) {
                if (this.isDecoy) {
                    oxoniumFragments.putAll(makeOxoniums(residue, true, glycoOxoniumDatabase, randomGenerator));
                } else {
                    oxoniumFragments.putAll(makeOxoniums(residue, false, glycoOxoniumDatabase, randomGenerator));

                }
            }
        }
    }

    /**
     * Helper method to add fragment ions to the oxonium list
     * @param residue residue type
     * @param isDecoy decoy or not
     * @return updated list
     */
    private TreeMap<String, GlycanFragment> makeOxoniums(GlycanResidue residue, boolean isDecoy, HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase, Random randomGenerator) {
        TreeMap<String, GlycanFragment> newFragments = new TreeMap<>();
        ArrayList<GlycanFragmentDescriptor> oxoniumIonDescriptors = glycoOxoniumDatabase.getOrDefault(residue, new ArrayList<>());
        for (GlycanFragmentDescriptor fragmentDescriptor : oxoniumIonDescriptors) {
            GlycanFragment newFragment = new GlycanFragment(fragmentDescriptor.requiredComposition, fragmentDescriptor.ruleProbabilies, fragmentDescriptor.massShift, isDecoy, randomGenerator, fragmentDescriptor.comment, GlycanFragment.FragType.Ox);
            newFragments.put(newFragment.hash, newFragment);
        }
        return newFragments;
    }


    /**
     * Take initialized oxonium ions (from initializeOxoniumFragments) and update them with propensity and intensity
     * information from the provided fragment info database
     * @param fragmentInfo fragment info container
     */
    private void updateOxoniums(GlycanCandidateFragments fragmentInfo) {
        for (Map.Entry<String, GlycanFragment> fragmentEntry : oxoniumFragments.entrySet()) {
            fragmentEntry.getValue().propensity = fragmentInfo.OxFragmentProps.getOrDefault(fragmentEntry.getKey(), 0.0);
            fragmentEntry.getValue().expectedIntensity = fragmentInfo.OxFragmentProps.getOrDefault(fragmentEntry.getKey(), 0.0);
        }
    }

    /**
     * Compute exact mass of a given composition
     * @return monoisotopic mass
     */
    public static double computeMonoisotopicMass(Map<GlycanResidue, Integer> glycanComposition) {
        double mass = 0;
        for (Map.Entry<GlycanResidue, Integer> glycanEntry : glycanComposition.entrySet()) {
            // mass = residue mass * residue count
            mass += (glycanEntry.getKey().mass * glycanEntry.getValue());
        }
        return mass;
    }

    /**
     * Generate a random mass shift within tolerancePPM about a randomly selected isotope peak in the
     * provided isotopes list.
     * @param isotopes list of isotopes
     * @param tolerancePPM Match tolerance (ppm) for glycan matching (from input parameter)
     * @param randomGenerator single random generator instance for whole glycan analysis
     * @return random mass shift within specified ranges
     */
    public static double getRandomShiftIsotopes(double glycanMass, Integer[] isotopes, double tolerancePPM, Random randomGenerator) {
        // randomly select isotope (must be sorted in ascending order)
        int minIso = isotopes[0];
        int maxIso = isotopes[isotopes.length - 1];
        // randomInt(0, max - min) + min yields correct range of min : max (including if min < 0)
        int isotope = randomGenerator.nextInt(maxIso + 1 - minIso) + minIso;  // upper bound is not inclusive, need to add 1 to get to max isotope

        // randomly generate mass shift within tolerance and add to chosen isotope
        double random = randomGenerator.nextDouble();       // between 0 and 1
        double baseMassEstimate = glycanMass + DEFAULT_PEPTIDE_MASS + isotope;
        double toleranceDa = baseMassEstimate * 1e-6 * tolerancePPM;
        double randomShift = -toleranceDa + random * (2 * toleranceDa);     // shift to range (min - random * (max - min)), where min = -toleranceDa and max = +toleranceDa
        return isotope * AAMasses.averagineIsotopeMass + randomShift;
    }

    /**
     * Return string representation of this glycan for writing to output tables
     * @return string
     */
    public String toString() {
        return GlycanFragment.toGlycanString(glycanComposition, monoisotopicMass, isDecoy);
    }

    public boolean containsResidueType(GlycanResidue residue) {
        return glycanComposition.containsKey(residue);
    }

}
