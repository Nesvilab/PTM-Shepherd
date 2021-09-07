package edu.umich.andykong.ptmshepherd.glyco;

import java.util.Map;

/**
 * Container for parsed fragment ion information parsed from the Fragment Database prior to initializing specific
 * GlycanFragment objects for the actual search.
 */
public class GlycanFragmentDescriptor {
    Map<GlycanResidue, Integer> requiredComposition;    // composition of this fragment ion
    double[] ruleProbabilies;                           // set of rule probabilities to use for this fragment
    double massShift;                                   // mass shift between composition mass and fragment mass (e.g., for loss of H2O in oxonium ions)

    public GlycanFragmentDescriptor(Map<GlycanResidue, Integer> requiredComposition, double[] ruleProbabilies, double massShift) {
        this.requiredComposition = requiredComposition;
        this.ruleProbabilies = ruleProbabilies;
        this.massShift = massShift;
    }
}
