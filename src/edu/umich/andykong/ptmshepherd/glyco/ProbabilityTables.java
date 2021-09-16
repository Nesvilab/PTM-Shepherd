package edu.umich.andykong.ptmshepherd.glyco;

import java.util.HashMap;

/**
 * Central location for holding probability tables for glycan assignment. Initializes with default values, but
 * can be updated by parameters (e.g. for testing)
 */
public class ProbabilityTables {
    public HashMap<Integer, Double> isotopeProbTable;
    public double massProbScaling;

    /* Y ion scores to generate probability ratios
     * Entries 0-3 are if ion is found in spectrum, 4-7 are if ion is not found. In both sets of 4, probs are
     * 0=allowed in both candidates
     * 1=allowed in candidate 1, not in candidate 2
     * 2=allowed in candidate 2, not in candidate 1
     * 3=allowed in neither candidate
     */
    public double[] regularYrules;
    public double[] dHexYrules;

    // Oxonium ion probs. Same format as Y ion probabilities
    public double[] neuacRules;
    public double[] neugcRules;
    public double[] phosphoRules;
    public double[] sulfoRules;
    public HashMap<GlycanResidue, double[]> rulesByResidue;

    // parameter names
    public static String[] probabilityParams = {"prob_neuacOx",
            "prob_neugcOx",
            "prob_phosphoOx",
            "prob_sulfoOx",
            "prob_regY",
            "prob_dhexY",
            "prob_mass",
            "prob_isotope"
    };

    // Constructor with default probability values
    public ProbabilityTables() {
        // Isotope error probability table
        isotopeProbTable = new HashMap<>();
        isotopeProbTable.put(-2, 0.125);
        isotopeProbTable.put(-1, 0.25);
        isotopeProbTable.put(0, 1.0);
        isotopeProbTable.put(1, 0.95);
        isotopeProbTable.put(2, 0.5);
        isotopeProbTable.put(3, 0.25);
        isotopeProbTable.put(4, 0.125);

        massProbScaling = 1.0;

        /* Y ion scores to generate probability ratios
         * Entries 0 is for
         */
        regularYrules = new double[]{5, 0.5};
        dHexYrules = new double[]{5, 0.5};

        // Oxonium ion probs. Same format as Y ion probabilities
        neuacRules = new double[]{10, 0.05};
        neugcRules = new double[]{10, 0.05};
        phosphoRules = new double[]{10, 0.05};
        sulfoRules = new double[]{10, 0.05};

        // initialize lookup table for rule types by residue
        rulesByResidue = new HashMap<>();
    }

    /**
     * Update the list of probabilities after parameters have been parsed
     */
    public void updateRulesByResidue() {
        rulesByResidue.put(GlycanResidue.NeuAc, neuacRules);
        rulesByResidue.put(GlycanResidue.NeuGc, neugcRules);
        rulesByResidue.put(GlycanResidue.Phospho, phosphoRules);
        rulesByResidue.put(GlycanResidue.Sulfo, sulfoRules);
    }



}
