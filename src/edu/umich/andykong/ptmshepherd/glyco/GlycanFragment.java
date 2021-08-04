package edu.umich.andykong.ptmshepherd.glyco;

import java.util.Map;
import java.util.Random;

/**
 * Container for glycan fragment ion information. Holds composition requirements, mass, and
 * pairwise scoring rules.
 */
public class GlycanFragment {
    double neutralMass;
    Map<GlycanResidue, Integer> requiredComposition;
    double[] ruleProbabilities;
    double foundIntensity;
    boolean isDecoy;
    public static final double MAX_DECOY_FRAGMENT_SHIFT_DA = 20;

    /**
     * Constructor for case that neutral mass is the mass of requiredComposition exactly
     * @param ruleProbabilities probabilities to use
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     * @param randomGenerator the single random number generator instance
     */
    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilities, boolean isDecoy, Random randomGenerator) {
        this.requiredComposition = requiredComposition;
        this.ruleProbabilities = ruleProbabilities;
        this.foundIntensity = 0;
        this.isDecoy = isDecoy;
        if (isDecoy) {
            this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition) + randomMassShift(MAX_DECOY_FRAGMENT_SHIFT_DA, randomGenerator);
        } else {
            this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition);
        }
    }

    /**
     * Constructor for cases where neutral mass needs to be supplied (e.g. some oxonium ions)
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     * @param ruleProbabilities probabilities to use
     * @param neutralMass neutral mass of the fragment
     * @param randomGenerator the single random number generator instance
     */
    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilities, double neutralMass, boolean isDecoy, Random randomGenerator) {
        this.requiredComposition = requiredComposition;
        this.ruleProbabilities = ruleProbabilities;
        this.foundIntensity = 0;
        this.isDecoy = isDecoy;
        if (isDecoy) {
            this.neutralMass = neutralMass + randomMassShift(MAX_DECOY_FRAGMENT_SHIFT_DA, randomGenerator);
        } else {
            this.neutralMass = neutralMass;
        }
    }

    /**
     * Determine if this fragment is an allowed fragment of the provided candidate composition.
     * For comparing candidates, if a given fragment ion is found, there are 4 possibilities for whether that
     * fragment is allowed/expected in both candidates, candidate 1 only, 2 only, or neither. This method
     * determines if the fragment is allowed for a single candidate, and thus should be run twice for each comparison
     * to figure out which of the 4 probabilities is applicable.
     * Note: only target fragments are allowed to match target candidates, and same for decoys
     *
     * Standard logic: if there are at least the required number of each required residue type in the candidate,
     * this is allowed.
     *
     * Special cases:
     * 1)   If Hex is required but not HexNAc, do NOT allow if the composition contains both. This
     *      is to distinguish glycation (Hex linkage) from glycosylation (HexNAc linkage)
     *
     *
     * @param candidate glycan candidate to consider
     * @return true if allowed/expected for this candidate, false if not
     */
    public boolean isAllowedFragment(GlycanCandidate candidate){
        // target fragments can only match target candidates and decoy fragments can only match decoy candidates
        if (this.isDecoy) {
            if (!candidate.isDecoy) {
                return false;
            }
        } else {
            if (candidate.isDecoy) {
                return false;
            }
        }

        // standard logic
        for (Map.Entry<GlycanResidue, Integer> requirementEntry : this.requiredComposition.entrySet()) {
            if (requirementEntry.getValue() > candidate.glycanComposition.get(requirementEntry.getKey())) {
                // more residues of this type required than found, return false
                return false;
            }
        }

        // special cases
        if (this.requiredComposition.containsKey(GlycanResidue.Hex) && this.requiredComposition.containsKey(GlycanResidue.HexNAc)) {
            if (this.requiredComposition.get(GlycanResidue.Hex) > 0 && !(this.requiredComposition.get(GlycanResidue.HexNAc) > 0)) {
                // Hex required but NOT HexNAc. Return false if candidate contains HexNAc
                if (candidate.glycanComposition.get(GlycanResidue.HexNAc) > 0) {
                    return false;
                }
            }
        }
        // all checks pass - return true
        return true;
    }

    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        int i=0;
        if (isDecoy) {
            stringBuilder.append("Decoy_");
        }
        for (Map.Entry<GlycanResidue, Integer> residue : requiredComposition.entrySet()) {
            if (residue.getValue() == 0) {
                continue;
            }
            if (i > 0) {
                stringBuilder.append("_");
            }
            i++;
            stringBuilder.append(String.format("%s-%d", GlycanMasses.outputGlycoNames.get(residue.getKey()), residue.getValue()));
        }
        stringBuilder.append(String.format("_%.0f", foundIntensity));
        return stringBuilder.toString();
    }

    /**
     * Generate a string suitable for detecting duplicate fragments (same required composition and decoy status).
     * Basically same as toString, but without intensity.
     * @return string of composition + decoy
     */
    public String toHashString() {
        StringBuilder stringBuilder = new StringBuilder();
        int i=0;
        if (isDecoy) {
            stringBuilder.append("Decoy_");
        }
        for (Map.Entry<GlycanResidue, Integer> residue : requiredComposition.entrySet()) {
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
     * Generate a random shift in mass, used for shifting decoy fragment ion masses and intact mass.
     * @param maxShift maximum size of shift (+/-) in Da
     * @param randomGenerator single random generator instance for whole glycan analysis
     * @return random shift
     */
    public static double randomMassShift(double maxShift, Random randomGenerator) {
        double random = randomGenerator.nextDouble();       // between 0 and 1
        return 1 + random * (maxShift - 1);                 // between 1 and maxShift
    }

}
