package edu.umich.andykong.ptmshepherd.glyco;

import java.util.Arrays;
import java.util.List;

public class GlycanMod {
    public String name;
    public double mass;
    public boolean isLabile;
    public int maxAllowed;

    public GlycanMod(String name, double mass, boolean isLabile, int maxAllowed, boolean isFixed, List<GlycanResidue> requiredResidues, ) {
        this.name = name;
        this.mass = mass;
        this.isLabile = isLabile;
        this.maxAllowed = maxAllowed;
    }

    public static GlycanMod parseMod(String line) {
        String[] splits = line.split("\t");
        double mass = Double.parseDouble(splits[1]);
        double yProbPlus = getOrDefault(splits[2]);
        double yProbMinus = getOrDefault(splits[3]);
        double oxoProbPlus = getOrDefault(splits[4]);
        double oxoProbMinus = getOrDefault(splits[5]);
        double oxoInt = getOrDefault(splits[6]);
        double[] diagnosticIons = new double[0];
        if (splits[7].length() > 0) {
            diagnosticIons = Arrays.stream(splits[7].split(",")).mapToDouble(Double::parseDouble).toArray();
        }
        String[] altNames = splits[8].split(",");
        return new GlycanMod(splits[0],
                mass,
                new double[]{yProbPlus, yProbMinus},
                new double[]{oxoProbPlus, oxoProbMinus, oxoInt},
                diagnosticIons,
                altNames
        );
    }

    // return -1 for any values not provided (not all rules are needed for all residues).
    private static double getOrDefault(String input) {
        if (input.length() > 0) {
            return Double.parseDouble(input);
        } else {
            return -1;
        }
    }
}
