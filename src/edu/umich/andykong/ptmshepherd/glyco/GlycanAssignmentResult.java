package edu.umich.andykong.ptmshepherd.glyco;

import java.util.TreeMap;

public class GlycanAssignmentResult {
    // Glycan Assignment results
    GlycanCandidate bestCandidate;
    GlycanCandidate bestTarget;
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

    // old-style results strings for preserving old outputs
    String glycanAssignmentString;
    String fullRawGlycoString;


    public GlycanAssignmentResult(GlycanCandidate bestCandidate, boolean isDecoyGlycan, double glycanScore, double glycanQval) {
        this.bestCandidate = bestCandidate;
        this.isDecoyGlycan = isDecoyGlycan;
        this.glycanScore = glycanScore;
        this.glycanQval = glycanQval;
    }

    public GlycanAssignmentResult(String peptide, float deltaMass, float pepMass, String assignedMods, String specName) {
        this.peptide = peptide;
        this.deltaMass = deltaMass;
        this.pepMass = pepMass;
        this.assignedMods = assignedMods;
        this.specName = specName;
    }


}
