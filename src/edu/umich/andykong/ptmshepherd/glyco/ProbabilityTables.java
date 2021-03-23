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
}
