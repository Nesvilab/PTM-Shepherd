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
    public GlycanCandidate(String glycanStr, String[] parsedFragmentInfo){
        glycanComposition = StaticGlycoUtilities.parseGlycanString(glycanStr);
        Yfragments = new TreeMap<>();
        oxoniumFragments = new TreeMap<>();

        for (String fragment : parsedFragmentInfo) {
            String[] typeSplits = fragment.split("~");
            // string format is [type]~[composition]~[intensity]
            if (typeSplits[0].matches("Y")) {
                // Y ion
                GlycanFragment newFragment = new GlycanFragment(typeSplits[1], Double.parseDouble(typeSplits[2]));
                Yfragments.put(newFragment.toHashString(), newFragment);
            } else if (typeSplits[0].matches("Ox")) {
                GlycanFragment newFragment = new GlycanFragment(typeSplits[1], Double.parseDouble(typeSplits[2]));
                oxoniumFragments.put(newFragment.toHashString(), newFragment);
            } else {
                // invalid
            }
        }
    }

    /**
     * Base constructor for a new glycan candidate for initial search (not using fragment propensities)
     * @param inputGlycanComp composition map
     * @param isDecoy bool
     * @param decoyType type of decoy (0 - 3)
     * @param glycoPPMtol MS1 tolerance
     * @param glycoIsotopes possible isotope errors
     * @param probabilityTable input prob ratios
     * @param glycoOxoniumDatabase oxonium ion input info
     * @param randomGenerator random object for randomizing masses if needed
     */
    public GlycanCandidate(Map<GlycanResidue, Integer> inputGlycanComp, boolean isDecoy, int decoyType, double glycoPPMtol, Integer[] glycoIsotopes, ProbabilityTables probabilityTable, HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase, Random randomGenerator) {
        this.glycanComposition = inputGlycanComp;
        this.isDecoy = isDecoy;
        // make sure that all residue types are accounted for (add Residue with 0 counts for any not included in the file)
        for (GlycanResidue residue : GlycanResidue.values()){
            if (!this.glycanComposition.containsKey(residue)) {
                this.glycanComposition.put(residue, 0);
            }
        }
        this.monoisotopicMass = setMassHelper(decoyType, glycoPPMtol, glycoIsotopes, randomGenerator);

        // initialize fragments for this candidate
        initializeYFragments(probabilityTable, randomGenerator);
        initializeOxoniumFragments(glycoOxoniumDatabase, randomGenerator);
        this.hasFragmentProps = false;
    }

    /**
     * Constructor for a new glycan candidate for a second search using fragment propensities
     * @param inputGlycanComp composition map
     * @param isDecoy bool
     * @param decoyType type of decoy (0 - 3)
     * @param glycoPPMtol MS1 tolerance
     * @param glycoIsotopes possible isotope errors
     * @param yFragmentProps input propensities for Y ions
     * @param OxFragmentProps input propensities for oxonium ions
     * @param randomGenerator random object for randomizing masses if needed
     */
    public GlycanCandidate(Map<GlycanResidue, Integer> inputGlycanComp, HashMap<String, Double> yFragmentProps, HashMap<String, Double> OxFragmentProps, boolean isDecoy, int decoyType, double glycoPPMtol, Integer[] glycoIsotopes, Random randomGenerator, HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase) {
        this.glycanComposition = inputGlycanComp;
        this.isDecoy = isDecoy;
        // make sure that all residue types are accounted for (add Residue with 0 counts for any not included in the file)
        for (GlycanResidue residue : GlycanResidue.values()){
            if (!this.glycanComposition.containsKey(residue)) {
                this.glycanComposition.put(residue, 0);
            }
        }
        this.monoisotopicMass = setMassHelper(decoyType, glycoPPMtol, glycoIsotopes, randomGenerator);

        // initialize fragments for this candidate
        initializeYFragmentsFromProps(yFragmentProps, randomGenerator);
        initializeOxoniumFragments(glycoOxoniumDatabase, randomGenerator);
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
            this.Yfragments.put(entry.getKey(), new GlycanFragment(entry.getValue().requiredComposition, entry.getValue().ruleProbabilities, entry.getValue().isDecoy, entry.getValue().neutralMass, entry.getValue().expectedIntensity, entry.getValue().propensity));
        }
        this.oxoniumFragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> entry : oxoniumFragments.entrySet()) {
            this.oxoniumFragments.put(entry.getKey(), new GlycanFragment(entry.getValue().requiredComposition, entry.getValue().ruleProbabilities, entry.getValue().isDecoy, entry.getValue().neutralMass, entry.getValue().expectedIntensity, entry.getValue().propensity));        }
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
    private double setMassHelper(int decoyType, double glycoPPMtol, Integer[] glycoIsotopes, Random randomGenerator) {
        double mass;
        if (! isDecoy) {
            mass = computeMonoisotopicMass(glycanComposition);
        } else {
            double baseMonoistopicMass = computeMonoisotopicMass(glycanComposition);
            double randomShift = 0;
            switch (decoyType) {
                case 0:
                    // simple mass window
                    randomShift = GlycanFragment.randomMassShift(MAX_CANDIDATE_DECOY_SHIFT_DA, randomGenerator);
                    break;
                case 1:
                    // random isotope and mass error
                    randomShift = getRandomShiftIsotopes(baseMonoistopicMass, glycoIsotopes, glycoPPMtol, randomGenerator);
                    break;
                case 2:
                    // random mass error, no isotope error
                    Integer[] noIsotopes = {0};
                    randomShift = getRandomShiftIsotopes(baseMonoistopicMass, noIsotopes, glycoPPMtol, randomGenerator);
                case 3:
                    // exact target mass - random shift left at 0
                    break;
            }
            mass = baseMonoistopicMass + randomShift;
        }
        return mass;
    }

    /**
     * Initialize array of all fragment ions to search for this candidate. Candidate has
     * fragment Ys up to the max HexNAc, Hex, and dHex present. Decoy fragments are generated for decoy candidates.
     * Fragment propensities are used for all fragments found in the bootstrap/input data (specified in the input map)
     * and fragments lacking any input are assumed to have 0 input propensity.
     */
    public void initializeYFragmentsFromProps(HashMap<String, Double> fragmentProps, Random randomGenerator) {
        // Initialize list of all Y fragments to consider. Currently using only HexNAc, Hex, and dHex in Y ions
        this.Yfragments = new TreeMap<>();
        for (int hexnac = 0; hexnac <= this.glycanComposition.get(GlycanResidue.HexNAc); hexnac++) {
            for (int hex = 0; hex <= this.glycanComposition.get(GlycanResidue.Hex); hex++) {
                if (!(hexnac == 0 && hex == 0)) {
                    // todo: option or remove
                    if (hexnac == 0 && this.glycanComposition.get(GlycanResidue.HexNAc) > 0) {
                        // do not generate non-HexNAc containing Y ions for compositions that have HexNAc
                        continue;
                    }
                    // add "regular" (no dHex) Y fragment for this HexNAc/Hex combination
                    Map<GlycanResidue, Integer> composition = new HashMap<>();
                    composition.put(GlycanResidue.HexNAc, hexnac);
                    composition.put(GlycanResidue.Hex, hex);
                    GlycanFragment fragment = new GlycanFragment(composition, fragmentProps, this.isDecoy, randomGenerator);
                    Yfragments.put(fragment.toHashString(), fragment);
                }
                for (int dHex = 1; dHex <= this.glycanComposition.get(GlycanResidue.dHex); dHex++) {
                    // add dHex fragments (if allowed)
                    Map<GlycanResidue, Integer> dHexcomposition = new HashMap<>();
                    dHexcomposition.put(GlycanResidue.HexNAc, hexnac);
                    dHexcomposition.put(GlycanResidue.Hex, hex);
                    dHexcomposition.put(GlycanResidue.dHex, dHex);
                    GlycanFragment dHexfragment = new GlycanFragment(dHexcomposition, fragmentProps, this.isDecoy, randomGenerator);
                    Yfragments.put(dHexfragment.toHashString(), dHexfragment);
                }
            }
        }
    }

    /**
     * Initialize array of all fragment ions to search for this candidate. Candidate has
     * fragment Ys up to the max HexNAc, Hex, and dHex present. Decoy fragments are generated for decoy candidates.
     */
    public void initializeYFragments(ProbabilityTables probabilityTable, Random randomGenerator) {
        // Initialize list of all Y fragments to consider. Currently using only HexNAc, Hex, and dHex in Y ions
        this.Yfragments = new TreeMap<>();
        for (int hexnac = 0; hexnac <= this.glycanComposition.get(GlycanResidue.HexNAc); hexnac++) {
            for (int hex = 0; hex <= this.glycanComposition.get(GlycanResidue.Hex); hex++) {
                if (!(hexnac == 0 && hex == 0)) {
                    // todo: option or remove
                    if (hexnac == 0 && this.glycanComposition.get(GlycanResidue.HexNAc) > 0) {
                        // do not generate non-HexNAc containing Y ions for compositions that have HexNAc
                        continue;
                    }
                    // add "regular" (no dHex) Y fragment for this HexNAc/Hex combination
                    Map<GlycanResidue, Integer> composition = new HashMap<>();
                    composition.put(GlycanResidue.HexNAc, hexnac);
                    composition.put(GlycanResidue.Hex, hex);
                    GlycanFragment fragment = new GlycanFragment(composition, probabilityTable.regularYrules, this.isDecoy, randomGenerator);
                    Yfragments.put(fragment.toHashString(), fragment);
                }
                for (int dHex = 1; dHex <= this.glycanComposition.get(GlycanResidue.dHex); dHex++) {
                    // add dHex fragments (if allowed)
                    Map<GlycanResidue, Integer> dHexcomposition = new HashMap<>();
                    dHexcomposition.put(GlycanResidue.HexNAc, hexnac);
                    dHexcomposition.put(GlycanResidue.Hex, hex);
                    dHexcomposition.put(GlycanResidue.dHex, dHex);
                    GlycanFragment dHexfragment = new GlycanFragment(dHexcomposition, probabilityTable.dHexYrules, this.isDecoy, randomGenerator);
                    Yfragments.put(dHexfragment.toHashString(), dHexfragment);
                }
            }
        }
    }

    /**
     * Helper method to initialize hard-coded oxonium fragment rules. Only initializes fragments for a
     * residue type if at least one candidate contains that residue type (no need to consider if not).
     * Decoys generated for all residue types that have at least one decoy candidate containing that type.
     */
    public void initializeOxoniumFragments(HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase, Random randomGenerator) {
        this.oxoniumFragments = new TreeMap<>();
        // HexNAc, Hex oxoniums

        // NeuAc
        if (this.glycanComposition.get(GlycanResidue.NeuAc) > 0) {
            if (this.isDecoy) {
                oxoniumFragments.putAll(makeOxoniums(GlycanResidue.NeuAc, true, glycoOxoniumDatabase, randomGenerator));
            } else {
                oxoniumFragments.putAll(makeOxoniums(GlycanResidue.NeuAc, false, glycoOxoniumDatabase, randomGenerator));
            }
        }
        // NeuGc
        if (this.glycanComposition.get(GlycanResidue.NeuGc) > 0) {
            if (this.isDecoy) {
                oxoniumFragments.putAll(makeOxoniums(GlycanResidue.NeuGc, true, glycoOxoniumDatabase, randomGenerator));
            } else {
                oxoniumFragments.putAll(makeOxoniums(GlycanResidue.NeuGc, false, glycoOxoniumDatabase, randomGenerator));
            }
        }
        // Phospho-Hex
        if (this.glycanComposition.get(GlycanResidue.Phospho) > 0) {
            if (this.isDecoy) {
                oxoniumFragments.putAll(makeOxoniums(GlycanResidue.Phospho, true, glycoOxoniumDatabase, randomGenerator));
            } else {
                oxoniumFragments.putAll(makeOxoniums(GlycanResidue.Phospho, false, glycoOxoniumDatabase, randomGenerator));
            }
        }
        // Sulfo
        if (this.glycanComposition.get(GlycanResidue.Sulfo) > 0) {
            if (this.isDecoy) {
                oxoniumFragments.putAll(makeOxoniums(GlycanResidue.Sulfo, true, glycoOxoniumDatabase, randomGenerator));
            } else {
                oxoniumFragments.putAll(makeOxoniums(GlycanResidue.Sulfo, false, glycoOxoniumDatabase, randomGenerator));
            }
        }
        // dHex
        if (this.glycanComposition.get(GlycanResidue.dHex) > 0) {
            if (this.isDecoy) {
                oxoniumFragments.putAll(makeOxoniums(GlycanResidue.dHex, true, glycoOxoniumDatabase, randomGenerator));
            } else {
                oxoniumFragments.putAll(makeOxoniums(GlycanResidue.dHex, false, glycoOxoniumDatabase, randomGenerator));
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
        ArrayList<GlycanFragmentDescriptor> oxoniumIonDescriptors = glycoOxoniumDatabase.get(residue);
        for (GlycanFragmentDescriptor fragmentDescriptor : oxoniumIonDescriptors) {
            GlycanFragment newFragment = new GlycanFragment(fragmentDescriptor.requiredComposition, fragmentDescriptor.ruleProbabilies, fragmentDescriptor.massShift, isDecoy, randomGenerator);
            newFragments.put(newFragment.toHashString(), newFragment);
        }
        return newFragments;
    }

    /**
     * Compute exact mass of a given composition
     * @return monoisotopic mass
     */
    public static double computeMonoisotopicMass(Map<GlycanResidue, Integer> glycanComposition) {
        double mass = 0;
        for (Map.Entry<GlycanResidue, Integer> glycanEntry : glycanComposition.entrySet()) {
            // mass = residue mass * residue count
            mass += (GlycanMasses.glycoMasses.get(glycanEntry.getKey()) * glycanEntry.getValue());
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
        return toGlycanHash(glycanComposition, isDecoy);
    }

    public boolean containsResidueType(GlycanResidue residue) {
        return glycanComposition.containsKey(residue);
    }

    // Static helper for both glycan fragment and candidate
    public static String toGlycanHash(Map<GlycanResidue, Integer> glycanComposition, boolean isDecoy) {
        StringBuilder stringBuilder = new StringBuilder();
        if (isDecoy) {
            stringBuilder.append("Decoy_");
        }
        ArrayList<String> residues = new ArrayList<>();
        for (GlycanResidue residueKey : GlycanResidue.values()) {
//        for (Map.Entry<GlycanResidue, Integer> residue : glycanComposition.entrySet()) {
            if (glycanComposition.getOrDefault(residueKey, 0) > 0) {
                residues.add(String.format("%s-%d", GlycanMasses.outputGlycoNames.get(residueKey), glycanComposition.get(residueKey)));
            }
        }
        stringBuilder.append(String.join("_", residues));
        return stringBuilder.toString();
    }

    // Static helper for both glycan fragment and candidate
    public static String toGlycanString(Map<GlycanResidue, Integer> glycanComposition, boolean isDecoy, double intensity) {
        return String.format("%s~%.4f", toGlycanHash(glycanComposition, isDecoy), intensity);
    }
}
