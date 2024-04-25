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

import umich.ms.glyco.GlycanResidue;

import java.util.Map;

/**
 * Container for parsed fragment ion information parsed from the Fragment Database prior to initializing specific
 * GlycanFragment objects for the actual search.
 */
public class GlycanFragmentDescriptor {
    Map<GlycanResidue, Integer> requiredComposition;    // composition of this fragment ion
    double[] ruleProbabilies;                           // set of rule probabilities to use for this fragment
    double massShift;                                   // mass shift between composition mass and fragment mass (e.g., for loss of H2O in oxonium ions)
    String comment;                                     // comment describing mass shift

    public GlycanFragmentDescriptor(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilies, double massShift, String comment) {
        this.requiredComposition = requiredComposition;
        this.ruleProbabilies = ruleProbabilies;
        this.massShift = massShift;
        this.comment = comment;
    }
}
