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

import java.util.LinkedHashMap;

/**
 * Container for holding fragment propensities for a given glycan
 */
public class GlycanCandidateFragments {
    LinkedHashMap<String, Double> yFragmentIntensities;
    LinkedHashMap<String, Double> OxFragmentIntensities;
    LinkedHashMap<String, Double> generalOxFragmentIntensities;

    // Standard deviations for error bars
    LinkedHashMap<String, Double> yFragmentStdDevs;
    LinkedHashMap<String, Double> OxFragmentStdDevs;
    LinkedHashMap<String, Double> generalOxFragmentStdDevs;


    public GlycanCandidateFragments(LinkedHashMap<String, Double> yFragmentIntensities, LinkedHashMap<String, Double> OxFragmentIntensities, LinkedHashMap<String, Double> generalOxFragmentIntensities) {
        this.yFragmentIntensities = yFragmentIntensities;
        this.OxFragmentIntensities = OxFragmentIntensities;
        this.generalOxFragmentIntensities = generalOxFragmentIntensities;
        this.yFragmentStdDevs = new LinkedHashMap<>();
        this.OxFragmentStdDevs = new LinkedHashMap<>();
        this.generalOxFragmentStdDevs = new LinkedHashMap<>();
    }

    public GlycanCandidateFragments(LinkedHashMap<String, Double> yFragmentIntensities, LinkedHashMap<String, Double> OxFragmentIntensities, LinkedHashMap<String, Double> generalOxFragmentIntensities,
                                   LinkedHashMap<String, Double> yFragmentStdDevs, LinkedHashMap<String, Double> OxFragmentStdDevs, LinkedHashMap<String, Double> generalOxFragmentStdDevs) {
        this.yFragmentIntensities = yFragmentIntensities;
        this.OxFragmentIntensities = OxFragmentIntensities;
        this.generalOxFragmentIntensities = generalOxFragmentIntensities;
        this.yFragmentStdDevs = yFragmentStdDevs;
        this.OxFragmentStdDevs = OxFragmentStdDevs;
        this.generalOxFragmentStdDevs = generalOxFragmentStdDevs;
    }

    // empty constructor for candidates without fragment info
    public GlycanCandidateFragments() {
        this.yFragmentIntensities = new LinkedHashMap<>();
        this.OxFragmentIntensities = new LinkedHashMap<>();
        this.generalOxFragmentIntensities = new LinkedHashMap<>();
        this.yFragmentStdDevs = new LinkedHashMap<>();
        this.OxFragmentStdDevs = new LinkedHashMap<>();
        this.generalOxFragmentStdDevs = new LinkedHashMap<>();
    }
}
