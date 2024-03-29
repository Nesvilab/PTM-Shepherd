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

import java.util.HashMap;

/**
 * Container for holding fragment propensities for a given glycan
 */
public class GlycanCandidateFragments {
    HashMap<String, Double> yFragmentProps;
    HashMap<String, Double> yFragmentIntensities;
    HashMap<String, Double> OxFragmentProps;
    HashMap<String, Double> OxFragmentIntensities;


    public GlycanCandidateFragments(HashMap<String, Double> yFragmentProps, HashMap<String, Double> OxFragmentProps, HashMap<String, Double> yFragmentIntensities, HashMap<String, Double> OxFragmentIntensities) {
        this.yFragmentProps = yFragmentProps;
        this.OxFragmentProps = OxFragmentProps;
        this.yFragmentIntensities = yFragmentIntensities;
        this.OxFragmentIntensities = OxFragmentIntensities;
    }

    // empty constructor for candidates without fragment info
    public GlycanCandidateFragments() {
        this.yFragmentProps = new HashMap<>();
        this.OxFragmentProps = new HashMap<>();
        this.yFragmentIntensities = new HashMap<>();
        this.OxFragmentIntensities = new HashMap<>();
    }
}
