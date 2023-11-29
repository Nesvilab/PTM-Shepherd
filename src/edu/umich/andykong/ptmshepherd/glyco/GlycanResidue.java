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


/**
 * Allowed names for glyco residues. Matched to masses in GlycanMasses
 */

public class GlycanResidue implements Comparable<GlycanResidue> {
    public String name;
    public double mass;
    public double[] yProbs;
    public double[] oxoProbs;
    public String[] alternateNames;
    public boolean islabile;

    public GlycanResidue(String name, double mass, double[] yProbs, double[] oxoProbs, String[] altNames, boolean isLabile) {
        this.name = name;
        this.mass = mass;
        this.yProbs = yProbs;
        this.oxoProbs = oxoProbs;
        this.alternateNames = altNames;
        this.islabile = isLabile;
    }

    public static GlycanResidue parseResidue(String line) {
        String[] splits = line.replace("\"", "").split("\t");
        double mass = Double.parseDouble(splits[1]);
        double yProbPlus = getOrDefault(splits[2]);
        double yProbMinus = getOrDefault(splits[3]);
        double oxoProbPlus = getOrDefault(splits[4]);
        double oxoProbMinus = getOrDefault(splits[5]);
        double oxoInt = getOrDefault(splits[6]);
        String[] altNames = splits[7].split(",");
        boolean isLabile = Boolean.parseBoolean(splits[8]);
        return new GlycanResidue(splits[0],
                mass,
                new double[]{yProbPlus, yProbMinus},
                new double[]{oxoProbPlus, oxoProbMinus, oxoInt},
                altNames,
                isLabile
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

    public String toString() {
        return name;
    }

    @Override
    public int compareTo(GlycanResidue o) {
        return name.compareTo(o.name);
    }
}
