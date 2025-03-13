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

import java.util.HashMap;
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

    public GlycanFragment(Map<GlycanResidue, Integer> requiredComposition, FragType fragType, double[] ruleProbabilities, boolean isDecoy, double neutralMass, double expectedIntensity, double foundIntensity, double propensity, String compositionComment) {
        this.requiredComposition = requiredComposition;
        this.ruleProbabilities = ruleProbabilities;
        this.foundIntensity = foundIntensity;
        this.expectedIntensity = expectedIntensity;
        this.isDecoy = isDecoy;
        this.neutralMass = neutralMass;
        this.propensity = propensity;
        this.compositionComment = compositionComment;
        this.hash = toFragmentHash(requiredComposition, isDecoy, compositionComment);
        this.fragType = fragType;
    }

    /**
     * Initialize a new Fragment from a glycan string, type, and expected intensity.
     * @param glycanStr string to parse for composition
     * @param expectedIntensity observed relative intensity
     */
    public static GlycanFragment parseGlycanFragment(String glycanStr, double expectedIntensity, FragType fragType, HashMap<String, GlycanResidue> glycanResiduesMap) {
        Glycan glycan = GlycanParser.parseGlycanString(glycanStr, glycanResiduesMap);
        double neutralMass = Glycan.computeCompositionMass(glycan.composition);
        return new GlycanFragment(glycan.composition, fragType, new double[]{}, false, neutralMass, expectedIntensity, 0, 0, "");
    }

    /**
     * Initialize Y ion (neutral mass is the mass of requiredComposition exactly)
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     * @param randomGenerator the single random number generator instance
     */
    public static GlycanFragment initializeYFragment(Map<GlycanResidue, Integer> requiredComposition, boolean isDecoy, Random randomGenerator) {
        double[] ruleProbs = computeYRuleProbs(requiredComposition);
        double neutralMass;
        if (isDecoy) {
            neutralMass = Glycan.computeCompositionMass(requiredComposition) + randomMassShift(MAX_DECOY_FRAGMENT_SHIFT_DA, randomGenerator);
        } else {
            neutralMass = Glycan.computeCompositionMass(requiredComposition);
        }
        return new GlycanFragment(requiredComposition, FragType.Y, ruleProbs, isDecoy, neutralMass, 0, 0, 0, "");
    }

    /**
     * Initialize oxonium ion. In some cases, neutral mass needs to be adjusted by the provided mass shift and composition
     * comment (e.g., H2O loss).
     * @param requiredComposition map of residues and counts required to be in the candidate to match this fragment
     * @param ruleProbabilities probabilities to use
     * @param neutralMassShift neutral mass shift of the fragment relative to its composition (e.g., for H2O losses from oxonium ions)
     * @param randomGenerator the single random number generator instance
     */
    public static GlycanFragment initializeOxoniumFragment(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilities, double neutralMassShift, boolean isDecoy, Random randomGenerator, String compComment) {
        double expectedIntensity;
        if (ruleProbabilities.length == 2) {
            expectedIntensity = 0;    // not provided, set to negative value to ignore
        } else {
            expectedIntensity = ruleProbabilities[2];   // 3rd value is expected intensity
        }
        double neutralMass;
        if (isDecoy) {
            neutralMass = Glycan.computeCompositionMass(requiredComposition) + neutralMassShift + randomMassShift(MAX_DECOY_FRAGMENT_SHIFT_DA, randomGenerator);
        } else {
            neutralMass = Glycan.computeCompositionMass(requiredComposition) + neutralMassShift;
        }
        return new GlycanFragment(requiredComposition, FragType.Ox, ruleProbabilities, isDecoy, neutralMass, expectedIntensity, 0, 0, compComment);
    }

    /**
     * Constructor for generating new Fragments for searching. Needed to avoid carrying intensity over from one
     * object to another when searching in multiple threads. Intended to copy all basic info from another
     * fragment and sets found intensity to 0.
     * @param baseFragment Fragment to copy from
     */
    public static GlycanFragment copyFragment(GlycanFragment baseFragment) {
        return new GlycanFragment(baseFragment.requiredComposition,
                baseFragment.fragType,
                baseFragment.ruleProbabilities,
                baseFragment.isDecoy,
                baseFragment.neutralMass,
                baseFragment.expectedIntensity,
                0,
                baseFragment.propensity,
                baseFragment.compositionComment);
    }
    /**
     * Constructor for 2nd pass search for new Fragment with provided intensity and propensity
     * @param baseFragment Fragment to copy from
     */
    public static GlycanFragment copyFragmentWithPropensity(GlycanFragment baseFragment, double expectedIntensity, double propensity) {
        return new GlycanFragment(baseFragment.requiredComposition,
                baseFragment.fragType,
                baseFragment.ruleProbabilities,
                baseFragment.isDecoy,
                baseFragment.neutralMass,
                expectedIntensity,
                0,
                propensity,
                baseFragment.compositionComment);
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
    public boolean isAllowedFragment(GlycanCandidate candidate, HashMap<String, GlycanResidue> glycanResiduesMap){
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
        GlycanResidue hexNAc = GlycanParser.findResidueName("HexNAc", glycanResiduesMap);
        GlycanResidue hexRes = GlycanParser.findResidueName("Hex", glycanResiduesMap);
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
    private static double[] computeYRuleProbs(Map<GlycanResidue, Integer> requiredComposition) {
        double minProbPlus = 100000;
        double maxProbMinus = -1;
        for (GlycanResidue residue: requiredComposition.keySet()) {
            if (residue.yProbPlus < minProbPlus && residue.yProbPlus > 0) {
                minProbPlus = residue.yProbPlus;
            }
            if (residue.yProbMinus > maxProbMinus) {
                maxProbMinus = residue.yProbMinus;
            }
        }
        return new double[]{minProbPlus, maxProbMinus};
    }

    /**
     * Output format for printing to .rawglyco file
     */
    public String toString() {
        return toFragmentIntensityString(requiredComposition, neutralMass, isDecoy, compositionComment, foundIntensity);
    }

    /**
     * Hash string is dependent on comp and comments only, not mass. NOTE: should be checked for duplicates in
     * the case of neutral losses (e.g. oxonium ions)
     * format: (Decoy_)<Composition><Comment>
     * @param glycanComposition composition
     * @param isDecoy decoy status
     * @param comment optional string to append to distinguish between same composition but different fragment (e.g. NLs)
     * @return string
     */
    public static String toFragmentHash(Map<GlycanResidue, Integer> glycanComposition, boolean isDecoy, String comment) {
        StringBuilder stringBuilder = new StringBuilder();
        if (isDecoy) {
            stringBuilder.append("Decoy_");
        }
        stringBuilder.append(Glycan.toGlycanString(glycanComposition));
        stringBuilder.append(comment);
        return stringBuilder.toString();
    }

    /**
     * Fragment hash with mass. Same as GlycanCandidate string but with (optional) composition comment
     * Format: (Decoy_)<Composition><Comment> % <Mass>
     */
    public static String toFragmentMassString(Map<GlycanResidue, Integer> glycanComposition, double mass, boolean isDecoy, String comment) {
        return GlycanFragment.toFragmentHash(glycanComposition, isDecoy, comment) + String.format(" %% %.4f", mass);
    }

    /**
     * Fragment hash with mass and intensity.
     * Format: (Decoy_)<Composition><Comment> % <Mass>~<Intensity>
     */
    public static String toFragmentIntensityString(Map<GlycanResidue, Integer> glycanComposition, double mass, boolean isDecoy, String comment, double foundIntensity) {
        return String.format("%s~%.4f", GlycanFragment.toFragmentMassString(glycanComposition, mass, isDecoy, comment), foundIntensity);
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
