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
    double expectedIntensity;
    boolean isDecoy;
    public static final double MAX_DECOY_FRAGMENT_SHIFT_DA = 20;
    double propensity;      // for fragment-specific probability calculations only

    /**
     * Constructor for propensity/bootstrapped analysis, replaces rule probs with propensity
     * @param fragmentPropensities fragment propensity map (fragment hashstring : propensity) from previous/bootstrap data
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     * @param randomGenerator the single random number generator instance
     */
    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, Map<String, Double> fragmentPropensities, boolean isDecoy, Random randomGenerator) {
        this.requiredComposition = requiredComposition;
        this.foundIntensity = 0;
        this.isDecoy = isDecoy;
        this.propensity = fragmentPropensities.getOrDefault(this.toHashString(), 0.0);
        if (isDecoy) {
            this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition) + randomMassShift(MAX_DECOY_FRAGMENT_SHIFT_DA, randomGenerator);
        } else {
            this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition);
        }
    }

    /**
     * Fragment info constructor for fragments from the glycofrags file
     * @param glycanStr string to parse for composition
     * @param expectedIntensity observed relative intensity todo: might be able to remove this if not using later
     */
    public GlycanFragment(String glycanStr, double expectedIntensity) {
        this.requiredComposition = StaticGlycoUtilities.parseGlycanString(glycanStr);
        this.expectedIntensity = expectedIntensity;
        this.isDecoy = false;
        this.foundIntensity = 0;
        this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition);
    }

    /**
     * Constructor for case that neutral mass is the mass of requiredComposition exactly
     * @param ruleProbabilities probabilities to use
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     * @param randomGenerator the single random number generator instance
     */
    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilities, boolean isDecoy, Random randomGenerator) {
        this.requiredComposition = requiredComposition;
        this.ruleProbabilities = ruleProbabilities;
        // propensity default is 1 - miss probRatio to generate equivalent scoring in new system
        this.propensity = 1 - ruleProbabilities[1];
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
     * @param neutralMassShift neutral mass shift of the fragment relative to its composition (e.g., for H2O losses from oxonium ions)
     * @param randomGenerator the single random number generator instance
     */
    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilities, double neutralMassShift, boolean isDecoy, Random randomGenerator) {
        this.requiredComposition = requiredComposition;
        this.ruleProbabilities = ruleProbabilities;
        // propensity default is 1 - miss probRatio to generate equivalent scoring in new system
        this.propensity = 1 - ruleProbabilities[1];
        this.foundIntensity = 0;
        if (ruleProbabilities.length == 2) {
            this.expectedIntensity = -1;    // not provided, set to negative value to ignore
        } else {
            this.expectedIntensity = ruleProbabilities[2];   // 3rd value is expected intensity
        }
        this.isDecoy = isDecoy;
        if (isDecoy) {
            this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition) + neutralMassShift + randomMassShift(MAX_DECOY_FRAGMENT_SHIFT_DA, randomGenerator);
        } else {
            this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition) + neutralMassShift;
        }
    }

    /**
     * Constructor for generating new Fragments for searching. Needed to avoid carrying intensity over from one
     * object to another when searching in multiple threads. Intended to copy all basic info from another
     * fragment and sets intensity to 0.
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     * @param ruleProbabilities probabilities to use
     * @param isDecoy decoy status
     * @param exactMass mass of a fragment (i.e., fragment.neutralMass)
     */
    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilities, boolean isDecoy, double exactMass, double expectedIntensity, double propensity) {
        this.requiredComposition = requiredComposition;
        this.ruleProbabilities = ruleProbabilities;
        this.isDecoy = isDecoy;
        this.neutralMass = exactMass;
        this.foundIntensity = 0;
        this.expectedIntensity = expectedIntensity;
        this.propensity = propensity;
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
        return GlycanCandidate.toGlycanString(requiredComposition, isDecoy, foundIntensity);
    }

    /**
     * Generate a string suitable for detecting duplicate fragments (same required composition and decoy status).
     * Basically same as toString, but without intensity.
     * @return string of composition + decoy
     */
    public String toHashString() {
        return GlycanCandidate.toGlycanHash(requiredComposition, isDecoy);
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
