package edu.umich.andykong.ptmshepherd.glyco;

import java.util.TreeMap;

public class GlycanAssignmentResult {
    // Glycan Assignment results
    GlycanCandidate bestCandidate;
    GlycanCandidate bestTarget;
    GlycanCandidate bestDecoy;
    double bestDecoyScore;
    boolean isDecoyGlycan;
    double glycanScore;
    double bestTargetScore;
    double glycanQval;
    TreeMap<GlycanFragment, Float> glycanFragments;     // might not be needed since candidate remembers its fragment intensities...
    
    // Basic PSM info (prior to PTM-S)
    String peptide;
    float deltaMass;
    float pepMass;
    String assignedMods;
    String specName;

    // old-style results strings for printing diagnostic, glycan outputs
    String glycanAssignmentString;
    String diagnosticResultString;


    public GlycanAssignmentResult(String peptide, float deltaMass, float pepMass, String assignedMods, String specName) {
        this.peptide = peptide;
        this.deltaMass = deltaMass;
        this.pepMass = pepMass;
        this.assignedMods = assignedMods;
        this.specName = specName;

        // initialize placeholder values
        this.bestTargetScore = Double.NaN;
        this.bestTarget = new GlycanCandidate();
        this.bestCandidate = new GlycanCandidate();
        this.bestDecoy = new GlycanCandidate();
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

        if (deltaMass > 3.5 || deltaMass < -1.5) {
            // main glycan results
            if (!isDecoyGlycan) {
                // for target glycans, append best decoy as well
                if (bestDecoy != null && !Double.isNaN(bestDecoyScore)) {
                    sb.append(String.format("\t%s\t%.2f\t\t%s\t%.2f", bestCandidate, glycanScore, bestDecoy, bestDecoyScore));
                } else {
                    sb.append(String.format("\t%s\t%.2f\t\t%s\t", bestCandidate, glycanScore, "no decoy matches"));
                }
            } else {
                // for decoy glycans, append best target as well
                if (bestTarget != null && !Double.isNaN(bestTargetScore)) {
                    sb.append(String.format("\t%s\t%.2f\t\t%s\t%.2f", bestCandidate, glycanScore, bestTarget, bestTargetScore));
                } else {
                    sb.append(String.format("\t%s\t%.2f\t\t%s\t", bestCandidate, glycanScore, "no target matches"));
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
        sb.append("\n");
        return sb.toString();
    }

}
