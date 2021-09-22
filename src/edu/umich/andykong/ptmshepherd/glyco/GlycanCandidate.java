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
    public GlycanFragment[] Yfragments;
    public GlycanFragment[] oxoniumFragments;

    public GlycanCandidate(Map<GlycanResidue, Integer> inputGlycanComp, boolean isDecoy, int decoyType, double glycoPPMtol, Integer[] glycoIsotopes, ProbabilityTables probabilityTable, HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase, Random randomGenerator) {
        this.glycanComposition = inputGlycanComp;
        this.isDecoy = isDecoy;
        // make sure that all residue types are accounted for (add Residue with 0 counts for any not included in the file)
        for (GlycanResidue residue : GlycanResidue.values()){
            if (!this.glycanComposition.containsKey(residue)) {
                this.glycanComposition.put(residue, 0);
            }
        }
        // set mass
        if (! isDecoy) {
            this.monoisotopicMass = computeMonoisotopicMass(inputGlycanComp);
        } else {
            double baseMonoistopicMass = computeMonoisotopicMass(inputGlycanComp);
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
            this.monoisotopicMass = baseMonoistopicMass + randomShift;
        }

        // initialize fragments for this candidate
        initializeYFragments(probabilityTable, glycoOxoniumDatabase, randomGenerator);
        initializeOxoniumFragments(glycoOxoniumDatabase, randomGenerator);
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
    public GlycanCandidate(Map<GlycanResidue, Integer> inputGlycanComp, boolean isDecoy, double monoisotopicMass, GlycanFragment[] yfragments, GlycanFragment[] oxoniumFragments) {
        this.glycanComposition = inputGlycanComp;
        this.isDecoy = isDecoy;
        this.monoisotopicMass = monoisotopicMass;
        this.Yfragments = new GlycanFragment[yfragments.length];
        for (int i=0; i < yfragments.length; i++) {
            this.Yfragments[i] = new GlycanFragment(yfragments[i].requiredComposition, yfragments[i].ruleProbabilities, yfragments[i].isDecoy, yfragments[i].neutralMass, yfragments[i].expectedIntensity);
        }
        this.oxoniumFragments = new GlycanFragment[oxoniumFragments.length];
        for (int i=0; i < oxoniumFragments.length; i++) {
            this.oxoniumFragments[i] = new GlycanFragment(oxoniumFragments[i].requiredComposition, oxoniumFragments[i].ruleProbabilities, oxoniumFragments[i].isDecoy, oxoniumFragments[i].neutralMass, oxoniumFragments[i].expectedIntensity);
        }
    }

    /**
     * Initialize array of all fragment ions to search for this candidate. Candidate has
     * fragment Ys up to the max HexNAc, Hex, and dHex present. Decoy fragments are generated for decoy candidates.
     */
    public void initializeYFragments(ProbabilityTables probabilityTable, HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase, Random randomGenerator) {
        // Initialize list of all Y fragments to consider. Currently using only HexNAc, Hex, and dHex in Y ions
        ArrayList<GlycanFragment> yFragments = new ArrayList<>();
        for (int hexnac = 0; hexnac <= this.glycanComposition.get(GlycanResidue.HexNAc); hexnac++) {
            for (int hex = 0; hex <= this.glycanComposition.get(GlycanResidue.Hex); hex++) {
                if (!(hexnac == 0 && hex == 0)) {
                    if (hexnac == 0 && this.glycanComposition.get(GlycanResidue.HexNAc) > 0) {
                        // do not generate non-HexNAc containing Y ions for compositions that have HexNAc
                        continue;
                    }
                    // add "regular" (no dHex) Y fragment for this HexNAc/Hex combination
                    Map<GlycanResidue, Integer> composition = new HashMap<>();
                    composition.put(GlycanResidue.HexNAc, hexnac);
                    composition.put(GlycanResidue.Hex, hex);
                    GlycanFragment fragment = new GlycanFragment(composition, probabilityTable.regularYrules, this.isDecoy, randomGenerator);
                    yFragments.add(fragment);
                }
                for (int dHex = 1; dHex <= this.glycanComposition.get(GlycanResidue.dHex); dHex++) {
                    // add dHex fragments (if allowed)
                    Map<GlycanResidue, Integer> dHexcomposition = new HashMap<>();
                    dHexcomposition.put(GlycanResidue.HexNAc, hexnac);
                    dHexcomposition.put(GlycanResidue.Hex, hex);
                    dHexcomposition.put(GlycanResidue.dHex, dHex);
                    GlycanFragment dHexfragment = new GlycanFragment(dHexcomposition, probabilityTable.dHexYrules, this.isDecoy, randomGenerator);
                    yFragments.add(dHexfragment);
                }
            }
        }
        this.Yfragments = yFragments.toArray(new GlycanFragment[0]);
    }

    /**
     * Helper method to initialize hard-coded oxonium fragment rules. Only initializes fragments for a
     * residue type if at least one candidate contains that residue type (no need to consider if not).
     * Decoys generated for all residue types that have at least one decoy candidate containing that type.
     * @return list of GlycanFragments for oxonium ions
     */
    public void initializeOxoniumFragments(HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase, Random randomGenerator) {
        ArrayList<GlycanFragment> oxoniumList = new ArrayList<>();
        // HexNAc, Hex oxoniums

        // NeuAc
        if (this.glycanComposition.get(GlycanResidue.NeuAc) > 0) {
            if (this.isDecoy) {
                oxoniumList.addAll(makeOxoniums(GlycanResidue.NeuAc, true, glycoOxoniumDatabase, randomGenerator));
            } else {
                oxoniumList.addAll(makeOxoniums(GlycanResidue.NeuAc, false, glycoOxoniumDatabase, randomGenerator));
            }
        }
        // NeuGc
        if (this.glycanComposition.get(GlycanResidue.NeuGc) > 0) {
            if (this.isDecoy) {
                oxoniumList.addAll(makeOxoniums(GlycanResidue.NeuGc, true, glycoOxoniumDatabase, randomGenerator));
            } else {
                oxoniumList.addAll(makeOxoniums(GlycanResidue.NeuGc, false, glycoOxoniumDatabase, randomGenerator));
            }
        }
        // Phospho-Hex
        if (this.glycanComposition.get(GlycanResidue.Phospho) > 0) {
            if (this.isDecoy) {
                oxoniumList.addAll(makeOxoniums(GlycanResidue.Phospho, true, glycoOxoniumDatabase, randomGenerator));
            } else {
                oxoniumList.addAll(makeOxoniums(GlycanResidue.Phospho, false, glycoOxoniumDatabase, randomGenerator));
            }
        }
        // Sulfo
        if (this.glycanComposition.get(GlycanResidue.Sulfo) > 0) {
            if (this.isDecoy) {
                oxoniumList.addAll(makeOxoniums(GlycanResidue.Sulfo, true, glycoOxoniumDatabase, randomGenerator));
            } else {
                oxoniumList.addAll(makeOxoniums(GlycanResidue.Sulfo, false, glycoOxoniumDatabase, randomGenerator));
            }
        }
        // dHex
        if (this.glycanComposition.get(GlycanResidue.dHex) > 0) {
            if (this.isDecoy) {
                oxoniumList.addAll(makeOxoniums(GlycanResidue.dHex, true, glycoOxoniumDatabase, randomGenerator));
            } else {
                oxoniumList.addAll(makeOxoniums(GlycanResidue.dHex, false, glycoOxoniumDatabase, randomGenerator));
            }
        }
        this.oxoniumFragments = oxoniumList.toArray(new GlycanFragment[0]);
    }

    /**
     * Helper method to add fragment ions to the oxonium list
     * @param residue residue type
     * @param isDecoy decoy or not
     * @return updated list
     */
    private ArrayList<GlycanFragment> makeOxoniums(GlycanResidue residue, boolean isDecoy, HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase, Random randomGenerator) {
        ArrayList<GlycanFragment> newFragments = new ArrayList<>();
        ArrayList<GlycanFragmentDescriptor> oxoniumIonDescriptors = glycoOxoniumDatabase.get(residue);
        for (GlycanFragmentDescriptor fragmentDescriptor : oxoniumIonDescriptors) {
            newFragments.add(new GlycanFragment(fragmentDescriptor.requiredComposition, fragmentDescriptor.ruleProbabilies, fragmentDescriptor.massShift, isDecoy, randomGenerator));
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
        StringBuilder stringBuilder = new StringBuilder();
        int i=0;
        if (isDecoy) {
            stringBuilder.append("Decoy_");
        }
        for (Map.Entry<GlycanResidue, Integer> residue : glycanComposition.entrySet()) {
            if (residue.getValue() == 0) {
                continue;
            }
            if (i > 0) {
                stringBuilder.append("_");
            }
            i++;
            stringBuilder.append(String.format("%s-%d", GlycanMasses.outputGlycoNames.get(residue.getKey()), residue.getValue()));
        }
        return stringBuilder.toString();
    }

    /**
     * Generate a string representation of this composition in guaranteed order. Uses declaration
     * order of Residue types in the GlycanResidue enum as the output order
     * @return string
     */
    public String toHashString() {
        StringBuilder stringBuilder = new StringBuilder();
        for (GlycanResidue residueKey : GlycanResidue.values()) {
            stringBuilder.append(String.format("%s-%d-", GlycanMasses.outputGlycoNames.get(residueKey), glycanComposition.get(residueKey)));
        }
        return stringBuilder.toString();
    }

    public boolean containsResidueType(GlycanResidue residue) {
        return glycanComposition.containsKey(residue);
    }

}