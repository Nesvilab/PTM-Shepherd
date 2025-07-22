package edu.umich.andykong.ptmshepherd.glyco;

import umich.ms.glyco.GlycanCandidate;
import umich.ms.glyco.GlycanResidue;

import java.util.HashMap;

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
    double[] featureVec;

    // initialize with the base Candidate and add scores as they are computed
    public GlycanCandidateResult(GlycanCandidate candidate, HashMap<String, GlycanResidue> glycanResiduesMap) {
        super(candidate.composition, candidate.decoyMassShift, candidate.isDecoy, glycanResiduesMap,candidate.Yfragments, candidate.oxoniumFragments);
    }
}
