package edu.umich.andykong.ptmshepherd.glyco;

import java.util.HashMap;

/**
 * Container for holding fragment propensities for a given glycan
 */
public class GlycanCandidateFragments {
    HashMap<String, Double> yFragmentProps;
    HashMap<String, Double> OxFragmentProps;

    public GlycanCandidateFragments(HashMap<String, Double> yFragmentProps, HashMap<String, Double> OxFragmentProps) {
        this.yFragmentProps = yFragmentProps;
        this.OxFragmentProps = OxFragmentProps;
    }
}
