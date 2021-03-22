package edu.umich.andykong.ptmshepherd.glyco;

import java.util.Map;

/**
 * Container for glycan fragment ion information. Holds composition requirements, mass, and
 * pairwise scoring rules.
 */
public class GlycanFragment {
    double neutralMass;
    Map<GlycanResidue, Integer> requiredComposition;
    double[] ruleProbabilities;
    double foundIntensity;

    /**
     * Constructor for case that neutral mass is the mass of requiredComposition exactly
     * @param ruleProbabilities probabilities to use
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     */
    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilities) {
        this.requiredComposition = requiredComposition;
        this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition);
        this.ruleProbabilities = ruleProbabilities;
        this.foundIntensity = 0;
    }

    /**
     * Constructor for cases where neutral mass needs to be supplied (e.g. some oxonium ions)
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     * @param ruleProbabilities probabilities to use
     * @param neutralMass neutral mass of the fragment
     */
    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilities, double neutralMass) {
        this.requiredComposition = requiredComposition;
        this.neutralMass = neutralMass;
        this.ruleProbabilities = ruleProbabilities;
        this.foundIntensity = 0;
    }

    /**
     * Determine if this fragment is an allowed fragment of the provided candidate composition.
     * For comparing candidates, if a given fragment ion is found, there are 4 possibilities for whether that
     * fragment is allowed/expected in both candidates, candidate 1 only, 2 only, or neither. This method
     * determines if the fragment is allowed for a single candidate, and thus should be run twice for each comparison
     * to figure out which of the 4 probabilities is applicable.
     * Standard logic: if there are at least the required number of each required residue type in the candidate,
     * this is allowed.
     *
     * Special cases:
     * 1)   If Hex is required but not HexNAc, do NOT allow if the composition contains both. This
     *      is to distinguish glycation (Hex linkage) from glycosylation (HexNAc linkage)
     *
     *
     * @param candidateComposition glycan candidate's composition to consider
     * @return true if allowed/expected for this candidate, false if not
     */
    public boolean isAllowedFragment(Map<GlycanResidue, Integer> candidateComposition){
        // standard logic
        for (Map.Entry<GlycanResidue, Integer> requirementEntry : this.requiredComposition.entrySet()) {
            if (requirementEntry.getValue() > candidateComposition.get(requirementEntry.getKey())) {
                // more residues of this type required than found, return false
                return false;
            }
        }

        // special cases
        if (this.requiredComposition.containsKey(GlycanResidue.Hex) && !this.requiredComposition.containsKey(GlycanResidue.HexNAc)){
            // Hex required but NOT HexNAc. Return false if candidate contains HexNAc
            if (candidateComposition.get(GlycanResidue.HexNAc) > 0){
                return false;
            }
        }
        // all checks pass - return true
        return true;
    }

}
