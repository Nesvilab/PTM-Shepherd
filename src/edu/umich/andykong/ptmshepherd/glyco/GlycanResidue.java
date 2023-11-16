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
import java.util.Arrays;
import java.util.Comparator;

/**
 * Allowed names for glyco residues. Matched to masses in GlycanMasses
 */

public class GlycanResidue implements Comparator<GlycanResidue> {
    public String name;
    public double mass;
    public double[] yProbs;
    public double[] oxoProbs;
    public String[] alternateNames;
    public double[] diagnosticIons;

    public GlycanResidue(String name, double mass, double[] yProbs, double[] oxoProbs, double[] diagnosticIons, String[] altNames) {
        this.name = name;
        this.mass = mass;
        this.yProbs = yProbs;
        this.oxoProbs = oxoProbs;
        this.diagnosticIons = diagnosticIons;
        this.alternateNames = altNames;
    }

    public static GlycanResidue parseResidue(String line) {
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
        return new GlycanResidue(splits[0],
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

    @Override
    public int compare(GlycanResidue o1, GlycanResidue o2) {
        return o1.name.compareTo(o2.name);
    }
}
