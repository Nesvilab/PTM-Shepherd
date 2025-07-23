package edu.umich.andykong.ptmshepherd.glyco;

import umich.ms.glyco.GlycanCandidate;
import umich.ms.glyco.GlycanFragment;
import umich.ms.glyco.GlycanResidue;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 * Container for a GlycanCandidate and all its associated results and scores.
 * A GlycanAssignmentResult will contain a list of GlycanCandidateResults for a given PSM.
 */
public class GlycanCandidateResult extends GlycanCandidate {
    public double glycanScore;
    public double ldaScore;
    public double summedScore;
    public double massError;
    public int isotope;
    double YFragmentScore;
    double OxFragmentScore;
    double massErrorScore;
    double isotopeScore;
    double YproportionScore;
    double ms1Score;
    double frequencyPrior;
    double[] featureVec;

    // initialize with the base Candidate and add scores as they are computed
    public GlycanCandidateResult(GlycanCandidate candidate, HashMap<String, GlycanResidue> glycanResiduesMap) {
        super(candidate.composition, candidate.decoyMassShift, candidate.isDecoy, glycanResiduesMap, candidate.Yfragments, candidate.oxoniumFragments);

        // deep copy the fragments
        TreeMap<String, GlycanFragment> copiedYfragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> entry : candidate.Yfragments.entrySet()) {
            copiedYfragments.put(entry.getKey(), GlycanFragment.copyFragment(entry.getValue()));
        }
        TreeMap<String, GlycanFragment> copiedOxoniumFragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> entry : candidate.oxoniumFragments.entrySet()) {
            copiedOxoniumFragments.put(entry.getKey(), GlycanFragment.copyFragment(entry.getValue()));
        }
        this.Yfragments = copiedYfragments;
        this.oxoniumFragments = copiedOxoniumFragments;
    }
}
