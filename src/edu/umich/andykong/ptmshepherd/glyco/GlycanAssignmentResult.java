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

import edu.umich.andykong.ptmshepherd.Mod;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import umich.ms.glyco.GlycanFragment;
import java.util.ArrayList;

public class GlycanAssignmentResult {
    // Glycan Assignment results
    public GlycanCandidateResult bestCandidate;
    public GlycanCandidateResult bestTarget;
    public GlycanCandidateResult bestDecoy;
    public boolean isDecoyGlycan;
    public double glycanScore;
    public double summedScore;      // for comparing to LDA score (debugging)
    public double glycanQval;

    public ArrayList<GlycanCandidateResult> allCandidates = new ArrayList<>(); // all glycan candidates for this PSM
    public boolean foundGlycan = false;

    // Variable mod check results
    public String modChangeDescription = "";  // e.g., "removed 3S(+79.9663)" or empty if original was best
    public Mod removedMod = null;  // the mod that was removed, null if original was best

    // Basic PSM info (prior to PTM-S)
    public int psmLineIndex; // index of the PSM line in the input file
    String peptide;
    float deltaMass;
    float pepMass;
    String assignedMods;
    String specName;
    ArrayList<Mod> assignedModsList;

    public GlycanAssignmentResult(int psmLineIndex, String peptide, float deltaMass, float pepMass, String assignedMods, ArrayList<Mod> assignedModsList, String specName) {
        this.psmLineIndex = psmLineIndex;
        this.peptide = peptide;
        this.deltaMass = deltaMass;
        this.pepMass = pepMass;
        this.assignedMods = assignedMods;
        this.assignedModsList = assignedModsList;
        this.specName = specName;

        this.bestTarget = null;
        this.bestCandidate = null;
        this.bestDecoy = null;
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
        if (deltaMass > 3.5 || deltaMass < -1.5) {
            if (foundGlycan) {
                printBestGlycan(sb);
                printNextBestGlycan(sb);
                printFeatureVector(sb, bestCandidate);
                printFragments(sb, bestCandidate);
            } else {
                sb.append("\tno glycan matched");
            }
        }
        sb.append("\n");
        return sb.toString();
    }

    private void printBestGlycan(StringBuilder sb) {
        sb.append(String.format("\t%s\t%.2f\t%.4f", bestCandidate, glycanScore, glycanQval));
    }

    private void printFeatureVector(StringBuilder sb, GlycanCandidateResult candidate) {
        if (candidate.featureVec != null) {
            for (double feature : candidate.featureVec) {
                sb.append(String.format("\t%.4f", feature)); // append each feature value
            }
        }
    }

    // Print best opposite glycan (target or decoy)
    private void printNextBestGlycan(StringBuilder sb) {
        if (!isDecoyGlycan) {
            // for target glycans, append best decoy as well
            if (bestDecoy != null) {
                sb.append(String.format("\t%s\t%.2f", bestDecoy, bestDecoy.glycanScore));
            } else {
                sb.append(String.format("\t%s\t", "no decoy matches"));
            }
        } else {
            // for decoy glycans, append best target as well
            if (bestTarget != null) {
                sb.append(String.format("\t%s\t%.2f", bestTarget, bestTarget.glycanScore));
            } else {
                sb.append(String.format("\t%s\t", "no target matches"));
            }
        }
    }

    private void printFragments(StringBuilder sb, GlycanCandidateResult candidate) {
        // glycan fragment info for target glycans
        if (!isDecoyGlycan) {
            // Y ions
            for (GlycanFragment ion : candidate.Yfragments.values()) {
                if (ion.foundIntensity > 0) {
                    sb.append(String.format("\tY~%s", ion));      // format is [ion type] [ion comp] [found intensity]
                }
            }
            // oxonium ions
            for (GlycanFragment ion : candidate.oxoniumFragments.values()) {
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

        if (foundGlycan) {
            printBestGlycan(sb);
            printFeatureVector(sb, bestCandidate);
            printFragments(sb, bestCandidate);
            sb.append("\n");
            for (int i = 1; i < allCandidates.size(); i++) {
                sb.append("\t\t\t\t");
                sb.append(String.format("\t%s\t%.2f\t", allCandidates.get(i), allCandidates.get(i).glycanScore));
                printFeatureVector(sb, allCandidates.get(i));
                printFragments(sb, allCandidates.get(i));
                sb.append("\n");
            }
        } else {
            sb.append("\n");
        }
        return sb.toString();
    }

}
