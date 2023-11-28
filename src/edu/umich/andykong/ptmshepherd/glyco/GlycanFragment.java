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

import java.util.ArrayList;
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
    String hash;
    String compositionComment;      // needed for cases with duplicate compositions (like oxonium ions with non-standard masses from fragmentation)
    FragType fragType;

    /**
     * Fragment info constructor for fragments from the glycofrags file
     * @param glycanStr string to parse for composition
     * @param expectedIntensity observed relative intensity todo: might be able to remove this if not using later
     */
    public GlycanFragment(String glycanStr, double expectedIntensity, FragType fragType, GlycoParams glycoParams) {
        this.requiredComposition = glycoParams.parseGlycanString(glycanStr);
        this.expectedIntensity = expectedIntensity;
        this.isDecoy = false;
        this.foundIntensity = 0;
        this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition);
        this.compositionComment = "";
        this.hash = toFragmentHash();
        this.fragType = fragType;
    }

    /**
     * Constructor Y ions (neutral mass is the mass of requiredComposition exactly)
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     * @param randomGenerator the single random number generator instance
     */
    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, boolean isDecoy, Random randomGenerator, FragType fragType) {
        this.requiredComposition = requiredComposition;
        this.ruleProbabilities = computeYRuleProbs();
        this.propensity = 0;
        this.expectedIntensity = 0;
        this.foundIntensity = 0;
        this.isDecoy = isDecoy;
        if (isDecoy) {
            this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition) + randomMassShift(MAX_DECOY_FRAGMENT_SHIFT_DA, randomGenerator);
        } else {
            this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition);
        }
        this.compositionComment = "";
        this.hash = toFragmentHash();
        this.fragType = fragType;
    }

    /**
     * Constructor for cases where neutral mass needs to be supplied (e.g. some oxonium ions)
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     * @param ruleProbabilities probabilities to use
     * @param neutralMassShift neutral mass shift of the fragment relative to its composition (e.g., for H2O losses from oxonium ions)
     * @param randomGenerator the single random number generator instance
     */
    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilities, double neutralMassShift, boolean isDecoy, Random randomGenerator, String compComment, FragType fragType) {
        this.requiredComposition = requiredComposition;
        this.ruleProbabilities = ruleProbabilities;
        this.propensity = 0;
        this.foundIntensity = 0;
        if (ruleProbabilities.length == 2) {
            this.expectedIntensity = 0;    // not provided, set to negative value to ignore
        } else {
            this.expectedIntensity = ruleProbabilities[2];   // 3rd value is expected intensity
        }
        this.isDecoy = isDecoy;
        if (isDecoy) {
            this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition) + neutralMassShift + randomMassShift(MAX_DECOY_FRAGMENT_SHIFT_DA, randomGenerator);
        } else {
            this.neutralMass = GlycanCandidate.computeMonoisotopicMass(requiredComposition) + neutralMassShift;
        }
        this.compositionComment = compComment;
        this.hash = toFragmentHash();
        this.fragType = fragType;
    }

    /**
     * Constructor for generating new Fragments for searching. Needed to avoid carrying intensity over from one
     * object to another when searching in multiple threads. Intended to copy all basic info from another
     * fragment and sets found intensity to 0.
     * @param baseFragment Fragment to copy from
     */
    public GlycanFragment(GlycanFragment baseFragment) {
        this.requiredComposition = baseFragment.requiredComposition;
        this.ruleProbabilities = baseFragment.ruleProbabilities;
        this.isDecoy = baseFragment.isDecoy;
        this.neutralMass = baseFragment.neutralMass;
        this.foundIntensity = 0;
        this.expectedIntensity = baseFragment.expectedIntensity;
        this.propensity = baseFragment.propensity;
        this.compositionComment = baseFragment.compositionComment;
        this.hash = baseFragment.hash;
        this.fragType = baseFragment.fragType;
    }
    /**
     * Constructor for 2nd pass search for new Fragment with provided intensity and propensity
     * @param baseFragment Fragment to copy from
     */
    public GlycanFragment(GlycanFragment baseFragment, double expectedIntensity, double propensity) {
        this.requiredComposition = baseFragment.requiredComposition;
        this.ruleProbabilities = baseFragment.ruleProbabilities;
        this.isDecoy = baseFragment.isDecoy;
        this.neutralMass = baseFragment.neutralMass;
        this.foundIntensity = 0;
        this.expectedIntensity = expectedIntensity;
        this.propensity = propensity;
        this.compositionComment = baseFragment.compositionComment;
        this.hash = baseFragment.hash;
        this.fragType = baseFragment.fragType;
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
    public boolean isAllowedFragment(GlycanCandidate candidate, GlycoParams glycoParams){
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

        // special cases todo: fix/generalize
        GlycanResidue hexNAc = glycoParams.findResidueName("HexNAc");
        GlycanResidue hexRes = glycoParams.findResidueName("Hex");
        if (this.requiredComposition.containsKey(hexRes) && this.requiredComposition.containsKey(hexNAc)) {
            if (this.requiredComposition.get(hexRes) > 0 && !(this.requiredComposition.get(hexNAc) > 0)) {
                // Hex required but NOT HexNAc. Return false if candidate contains HexNAc
                if (candidate.glycanComposition.get(hexNAc) > 0) {
                    return false;
                }
            }
        }
        // all checks pass - return true
        return true;
    }

    /**
     * Compute the generic scores for a Y ion from the provided scores of the individual residues. The score resulting
     * in the least change is taken if there are disagreements between residue scores.
     * @return
     */
    private double[] computeYRuleProbs() {
        double minProbPlus = 100000;
        double maxProbMinus = -1;
        for (GlycanResidue residue: requiredComposition.keySet()) {
            if (residue.yProbs[0] < minProbPlus && residue.yProbs[0] > 0) {
                minProbPlus = residue.yProbs[0];
            }
            if (residue.yProbs[1] > maxProbMinus) {
                maxProbMinus = residue.yProbs[1];
            }
        }
        return new double[]{minProbPlus, maxProbMinus};
    }

    /**
     * Output format for printing to .rawglyco file
     */
    public String toString() {
        return String.format("%s~%.4f", toGlycanString(requiredComposition, neutralMass, isDecoy), foundIntensity);
    }

    /**
     * Composition-only identifier method, to allow decoys with same comp but different masses to be matched (e.g.,
     * neutral loss oxonium ions)
     * @return string of composition + decoy
     */
    public String toFragmentHash() {
        return toGlycanCompString(requiredComposition, isDecoy, compositionComment);
    }

    /**
     * Hash string is dependent on comp and comments only, not mass. NOTE: should be checked for duplicates in
     * the case of neutral losses (e.g. oxonium ions)
     * @param glycanComposition composition
     * @param isDecoy decoy status
     * @param comment optional string to append to distinguish between same composition but different fragment (e.g. NLs)
     * @return string
     */
    public static String toGlycanCompString(Map<GlycanResidue, Integer> glycanComposition, boolean isDecoy, String comment) {
        StringBuilder stringBuilder = new StringBuilder();
        if (isDecoy) {
            stringBuilder.append("Decoy_");
        }
        ArrayList<String> residues = new ArrayList<>();
        for (GlycanResidue residueKey : glycanComposition.keySet()) {
            if (glycanComposition.getOrDefault(residueKey, 0) > 0) {
                residues.add(String.format("%s(%d)", residueKey.name, glycanComposition.get(residueKey)));
            }
        }
        stringBuilder.append(String.join("", residues));
        stringBuilder.append(comment);
        return stringBuilder.toString();
    }

    /**
     * Generate a string suitable for detecting duplicate fragments (same required composition and decoy status).
     * Basically same as toString, but without intensity.
     * @return string of composition + decoy + mass
     */
    public static String toGlycanString(Map<GlycanResidue, Integer> glycanComposition, double mass, boolean isDecoy) {
        return GlycanFragment.toGlycanCompString(glycanComposition, isDecoy, "") + String.format(" %% %.4f", mass);
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

    // Type of GlycanFragment
    public enum FragType {
        Y,
        Ox
    }
}
