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

import umich.ms.glyco.GlycanCandidate;
import umich.ms.glyco.GlycanFragment;

import java.util.ArrayList;

public class GlycanAssignmentResult {
    // Glycan Assignment results
    public GlycanCandidate bestCandidate;
    public GlycanCandidate bestTarget;
    public GlycanCandidate bestDecoy;
    public double bestDecoyScore;
    public boolean isDecoyGlycan;
    public double glycanScore;
    public double bestTargetScore;
    public double glycanQval;

    public ArrayList<GlycanCandidate> allCandidates = new ArrayList<>(); // all glycan candidates for this PSM
    public ArrayList<Double> allScores = new ArrayList<>(); // all glycan candidates for this PSM

    // Basic PSM info (prior to PTM-S)
    public int psmLineIndex; // index of the PSM line in the input file
    String peptide;
    float deltaMass;
    float pepMass;
    String assignedMods;
    String specName;

    // old-style results strings for printing diagnostic, glycan outputs
    String glycanAssignmentString;
    public static final String NO_GLYCAN_RESULT_STR = "No Glycan Matched";
    public boolean foundGlycan = false;     // true if glycan assignment was successful

    // LDA features
    double YFragmentScore;
    double OxFragmentScore;
    double massErrorScore;
    double isotopeScore;
    double YproportionScore;
    double KLscore;
    int yCount;
    float yHyper;
    public double[] featureVec;

    public GlycanAssignmentResult(int psmLineIndex, String peptide, float deltaMass, float pepMass, String assignedMods, String specName) {
        this.psmLineIndex = psmLineIndex;
        this.peptide = peptide;
        this.deltaMass = deltaMass;
        this.pepMass = pepMass;
        this.assignedMods = assignedMods;
        this.specName = specName;

        // initialize placeholder values
        this.bestTargetScore = Double.NaN;
        this.bestTarget = GlycanCandidate.emptyCandidate();
        this.bestCandidate = GlycanCandidate.emptyCandidate();
        this.bestDecoy = GlycanCandidate.emptyCandidate();
        allCandidates = new ArrayList<>();
    }

    /**
     * Format glycan fragment info for output to file for later fragment probability boostrapping.
     * Note: does not print glycan q-value because FDR performed after the initial run analysis
     * @return formatted string to print one line from this result
     */
    public String printGlycoFragmentInfo() {
        StringBuilder sb = new StringBuilder();
        // initial spectrum data
        sb.append(String.format("%s\t%s\t%s\t%.4f\t%.4f", specName, peptide, assignedMods, pepMass, deltaMass));
        if (glycanAssignmentString.matches("ERROR")) {
            // spectrum not found - print ERROR
            sb.append("\tERROR\n");
            return sb.toString();
        }
        if (glycanAssignmentString.contains(NO_GLYCAN_RESULT_STR)) {
            // no matching glycans found - leave result empty
            sb.append(glycanAssignmentString);
            sb.append("\n");
            return sb.toString();
        }

        if (deltaMass > 3.5 || deltaMass < -1.5) {
            printBestGlycan(sb);
        }
        sb.append("\n");
        return sb.toString();
    }

    private void printBestGlycan(StringBuilder sb) {
        if (!isDecoyGlycan) {
            // for target glycans, append best decoy as well
            if (bestDecoy != null && !Double.isNaN(bestDecoyScore)) {
                sb.append(String.format("\t%s\t%.2f\t%.4f\t%s\t%.2f", bestCandidate, glycanScore, glycanQval,bestDecoy, bestDecoyScore));
            } else {
                sb.append(String.format("\t%s\t%.2f\t%.4f\t%s\t", bestCandidate, glycanScore, glycanQval, "no decoy matches"));
            }
        } else {
            // for decoy glycans, append best target as well
            if (bestTarget != null && !Double.isNaN(bestTargetScore)) {
                sb.append(String.format("\t%s\t%.2f\t%.4f\t%s\t%.2f", bestCandidate, glycanScore, glycanQval, bestTarget, bestTargetScore));
            } else {
                sb.append(String.format("\t%s\t%.2f\t%.4f\t%s\t", bestCandidate, glycanScore, glycanQval, "no target matches"));
            }
        }

        if (featureVec != null) {
            for (double feature : featureVec) {
                sb.append(String.format("\t%.4f", feature)); // append each feature value
            }
        }

        // glycan fragment info for target glycans
        if (!isDecoyGlycan) {
            // Y ions
            for (GlycanFragment ion : bestCandidate.Yfragments.values()) {
                if (ion.foundIntensity > 0) {
                    sb.append(String.format("\tY~%s", ion));      // format is [ion type] [ion comp] [found intensity]
                }
            }
            // oxonium ions
            for (GlycanFragment ion : bestCandidate.oxoniumFragments.values()) {
                if (ion.foundIntensity > 0) {
                    sb.append(String.format("\tOx~%s", ion));      // format is [ion type] [ion comp] [found intensity]
                }
            }
        }
    }

    /**
     * Print all candidates for this PSM (intended for debugging)
     */
    public String printAllCandidates() {
        StringBuilder sb = new StringBuilder();
        // initial spectrum data
        sb.append(String.format("%s\t%s\t%s\t%.4f\t%.4f", specName, peptide, assignedMods, pepMass, deltaMass));
        if (glycanAssignmentString.matches("ERROR")) {
            // spectrum not found - print ERROR
            sb.append("\tERROR\n");
            return sb.toString();
        }
        if (glycanAssignmentString.contains(NO_GLYCAN_RESULT_STR)) {
            // no matching glycans found - leave result empty
            sb.append(glycanAssignmentString);
            sb.append("\n");
            return sb.toString();
        }

        if (foundGlycan) {
            printBestGlycan(sb);
            sb.append("\n");
            for (int i = 1; i < allCandidates.size(); i++) {
                sb.append("\t\t\t\t");
                sb.append(String.format("\t%s\t%.2f", allCandidates.get(i), allScores.get(i)));

                // glycan fragment info for target glycans
                if (!isDecoyGlycan) {
                    // Y ions
                    for (GlycanFragment ion : allCandidates.get(i).Yfragments.values()) {
                        if (ion.foundIntensity > 0) {
                            sb.append(String.format("\tY~%s", ion));      // format is [ion type] [ion comp] [found intensity]
                        }
                    }
                    // oxonium ions
                    for (GlycanFragment ion : allCandidates.get(i).oxoniumFragments.values()) {
                        if (ion.foundIntensity > 0) {
                            sb.append(String.format("\tOx~%s", ion));      // format is [ion type] [ion comp] [found intensity]
                        }
                    }
                }
                sb.append("\n");
            }
        } else {
            sb.append("\n");
        }
        return sb.toString();
    }

}
