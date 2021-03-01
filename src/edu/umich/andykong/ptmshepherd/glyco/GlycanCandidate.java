package edu.umich.andykong.ptmshepherd.glyco;

import java.util.Map;

/**
 Container for theoretical glycan compositions (supplied by user or database) to be
 searched against experimental results.
 */
public class GlycanCandidate {
    double monoisotopicMass;
    Map<GlycanResidue, Integer> glycanComposition;     // map of residue type: count of residue to describe the composition
    double[][] exactIsotopeCluster;

    public GlycanCandidate(Map<GlycanResidue, Integer> inputGlycanComp) {
        this.glycanComposition = inputGlycanComp;
        this.monoisotopicMass = computeMonoisotopicMass();
        this.exactIsotopeCluster = computeIsotopeCluster();
    }

    /**
     * Compute exact mass of this composition
     * @return monoisotopic mass
     */
    public double computeMonoisotopicMass() {
        double mass = 0;
        for (Map.Entry<GlycanResidue, Integer> glycanEntry : glycanComposition.entrySet()) {
            // mass = residue mass * residue count
            mass += (GlycanMasses.glycoMasses.get(glycanEntry.getKey()) * glycanEntry.getValue());
        }
        return mass;
    }

    public boolean containsResidueType(GlycanResidue residue) {
        return glycanComposition.containsKey(residue);
    }

    public double[][] computeIsotopeCluster() {
        // not yet implemented
        double[][] isoCluster = null;
        return isoCluster;
    }

}
