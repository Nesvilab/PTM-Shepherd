package edu.umich.andykong.ptmshepherd.glyco;

import java.util.HashMap;

/**
 * Central location for holding probability tables for glycan assignment
 */
public class ProbabilityTables {

    // Isotope error probability table
    static HashMap<Integer, Double> isotopeProbTable;
    static {
        isotopeProbTable = new HashMap<>();
        isotopeProbTable.put(-2, 0.125);
        isotopeProbTable.put(-1, 0.25);
        isotopeProbTable.put(0, 1.0);
        isotopeProbTable.put(1, 0.95);
        isotopeProbTable.put(2, 0.5);
        isotopeProbTable.put(3, 0.25);
        isotopeProbTable.put(4, 0.125);
    }

    // mass error probability table
    static HashMap<Integer, Double> massProbTable;
    static {
        massProbTable = new HashMap<>();
        massProbTable.put(0, 1.0);
        massProbTable.put(1, 0.9);
        massProbTable.put(2, 0.75);
        massProbTable.put(3, 0.5);
    }

    /* Y ion scores to generate probability ratios
     * Entries 0-3 are if ion is found in spectrum, 4-7 are if ion is not found. In both sets of 4, probs are
     * 0=allowed in both candidates
     * 1=allowed in candidate 1, not in candidate 2
     * 2=allowed in candidate 2, not in candidate 1
     * 3=allowed in neither candidate
     */
    static double[] regularYrules = {1, 1.5, 0.666666, 1, 1, 0.8, 1.25, 1};   // todo: read from params
    static double[] dHexYrules = {1, 2, 0.5, 1, 1, 1, 1, 1};            // todo: read from params
    static double[] twodHexYrules = {1, 2, 0.5, 1, 1, 1, 1, 1};          // todo: read from params

    // Oxonium ion probs. Same format as Y ion probabilities
    static double[] neuacRules = {1, 5, 0.2, 1, 1, 0.25, 4, 1};
    static double[] neugcRules = {1, 5, 0.2, 1, 1, 0.25, 4, 1};
    static double[] phosphoRules = {1, 2, 0.5, 1, 1, 0.5, 2, 1};
    static double[] sulfoRules = {1, 2, 0.5, 1, 1, 0.5, 2, 1};

}
