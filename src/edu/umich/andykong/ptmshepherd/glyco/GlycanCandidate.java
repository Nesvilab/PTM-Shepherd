package edu.umich.andykong.ptmshepherd.glyco;

import edu.umich.andykong.ptmshepherd.core.Spectrum;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;

/**
 Container for theoretical glycan compositions (supplied by user or database) to be
 searched against experimental results.
 */
public class GlycanCandidate {
    double monoisotopicMass;
    Map<GlycanResidue, Integer> glycanComposition;     // map of residue type: count of residue to describe the composition
    ArrayList<Float> expectedOxoniumIons;
    ArrayList<Float> disallowedOxoniumIons;
    ArrayList<Float> expectedYIons;
    ArrayList<Float> disallowedYIons;
    boolean isDecoy;
    public static final double MAX_CANDIDATE_DECOY_SHIFT_DA = 3;

    public GlycanCandidate(Map<GlycanResidue, Integer> inputGlycanComp, boolean isDecoy) {
        this.glycanComposition = inputGlycanComp;
        this.isDecoy = isDecoy;
        // make sure that all residue types are accounted for (add Residue with 0 counts for any not included in the file)
        for (GlycanResidue residue : GlycanResidue.values()){
            if (!this.glycanComposition.containsKey(residue)) {
                this.glycanComposition.put(residue, 0);
            }
        }
        // set mass
        if (! isDecoy) {
            this.monoisotopicMass = computeMonoisotopicMass(inputGlycanComp);
        } else {
            // randomly shift monoisotopic mass within range +/- 20 Da
            double randomShift = GlycanFragment.randomMassShift(MAX_CANDIDATE_DECOY_SHIFT_DA);
            this.monoisotopicMass = computeMonoisotopicMass(inputGlycanComp) + randomShift;
        }
    }

    /**
     * Compute exact mass of a given composition
     * @return monoisotopic mass
     */
    public static double computeMonoisotopicMass(Map<GlycanResidue, Integer> glycanComposition) {
        double mass = 0;
        for (Map.Entry<GlycanResidue, Integer> glycanEntry : glycanComposition.entrySet()) {
            // mass = residue mass * residue count
            mass += (GlycanMasses.glycoMasses.get(glycanEntry.getKey()) * glycanEntry.getValue());
        }
        return mass;
    }

    /**
     * Return string representation of this glycan for writing to output tables
     * @return string
     */
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        int i=0;
        if (isDecoy) {
            stringBuilder.append("Decoy_");
        }
        for (Map.Entry<GlycanResidue, Integer> residue : glycanComposition.entrySet()) {
            if (residue.getValue() == 0) {
                continue;
            }
            if (i > 0) {
                stringBuilder.append("_");
            }
            i++;
            stringBuilder.append(String.format("%s-%d", GlycanMasses.outputGlycoNames.get(residue.getKey()), residue.getValue()));
        }
        return stringBuilder.toString();
    }

    /**
     * Generate a string representation of this composition in guaranteed order. Uses declaration
     * order of Residue types in the GlycanResidue enum as the output order
     * @return string
     */
    public String toHashString() {
        StringBuilder stringBuilder = new StringBuilder();
        for (GlycanResidue residueKey : GlycanResidue.values()) {
            stringBuilder.append(String.format("%s-%d-", GlycanMasses.outputGlycoNames.get(residueKey), glycanComposition.get(residueKey)));
        }
        return stringBuilder.toString();
    }

    public boolean containsResidueType(GlycanResidue residue) {
        return glycanComposition.containsKey(residue);
    }

}
