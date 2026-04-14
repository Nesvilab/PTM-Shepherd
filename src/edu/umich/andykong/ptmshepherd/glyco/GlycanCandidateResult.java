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
    public double nnScore;
    public double summedScore;
    public double massError;
    public int isotope;
    double YFragmentScore;
    double OxFragmentScore;
    double massErrorScore;
    double ms1Score;
    double ms1DeltaScore;
    double klScore;
    double frequencyPrior;
    double ySpecSim;
    double oxSpecSim;
    double[] featureVec;

    // initialize with the base Candidate and add scores as they are computed
    public GlycanCandidateResult(GlycanCandidate candidate, HashMap<String, GlycanResidue> glycanResiduesMap) {
        super(candidate.composition, candidate.decoyMassShift, candidate.isDecoy, glycanResiduesMap, candidate.Yfragments, candidate.oxoniumFragments, candidate.generalOxoniumFragments);

        // deep copy the fragments
        TreeMap<String, GlycanFragment> copiedYfragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> entry : candidate.Yfragments.entrySet()) {
            copiedYfragments.put(entry.getKey(), GlycanFragment.copyFragment(entry.getValue()));
        }
        TreeMap<String, GlycanFragment> copiedOxoniumFragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> entry : candidate.oxoniumFragments.entrySet()) {
            copiedOxoniumFragments.put(entry.getKey(), GlycanFragment.copyFragment(entry.getValue()));
        }
        TreeMap<String, GlycanFragment> copiedGeneralOxoniumFragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> entry : candidate.generalOxoniumFragments.entrySet()) {
            copiedGeneralOxoniumFragments.put(entry.getKey(), GlycanFragment.copyFragment(entry.getValue()));
        }
        this.Yfragments = copiedYfragments;
        this.oxoniumFragments = copiedOxoniumFragments;
        this.generalOxoniumFragments = copiedGeneralOxoniumFragments;
    }


    /**
     * Normalize a map of fragment intensities within itself (e.g., Ys and oxoniums separately)
     *
     * @param fragmentIntensities input map
     */
    public static void normalizeIntensities(TreeMap<String, GlycanFragment> fragmentIntensities) {
        double maxIntensity = 0;
        for (Map.Entry<String, GlycanFragment> entry : fragmentIntensities.entrySet()) {
            if (entry.getValue().foundIntensity > maxIntensity) {
                maxIntensity = entry.getValue().foundIntensity;
            }
        }
        if (maxIntensity > 0) {
            for (Map.Entry<String, GlycanFragment> entry : fragmentIntensities.entrySet()) {
                entry.getValue().foundIntensity = entry.getValue().foundIntensity / maxIntensity;
            }
        }
    }
}
