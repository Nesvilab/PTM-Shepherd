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
    double[][] exactIsotopeCluster;
    ArrayList<Float> expectedOxoniumIons;
    ArrayList<Float> disallowedOxoniumIons;
    ArrayList<Float> expectedYIons;
    ArrayList<Float> disallowedYIons;

    public GlycanCandidate(Map<GlycanResidue, Integer> inputGlycanComp) {
        this.glycanComposition = inputGlycanComp;
        this.monoisotopicMass = computeMonoisotopicMass();
        initializeAllowedOxoniums();
        initializeAllowedYs();
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

    /**
     * Initialize mass list for allowed Y ions based on the number of Hex, HexNAc and presence of Fuc
     * in this composition. Also initialize disallowed Ys with (currently hardcoded) number of extra
     * HexNAc and Hex to check for evidence that this isn't the right composition.
     *
     */
    private void initializeAllowedYs() {
        expectedYIons = new ArrayList<>();
        disallowedYIons = new ArrayList<>();

        // HexNAc, Hex: simple (composition only) rule expects up to the total number of HexNAcs in the composition and disallows more than that
        int maxHexNAc = 0;
        int maxHex = 0;
        if (containsResidueType(GlycanResidue.HexNAc)) {
            maxHexNAc = glycanComposition.get(GlycanResidue.HexNAc);
        }
        if (containsResidueType(GlycanResidue.Hex)) {
            maxHex = glycanComposition.get(GlycanResidue.Hex);
        }
        // generate all allowed combinations
        for (int hexnac=0; hexnac <= maxHexNAc; hexnac++) {
            for (int hex=0; hex <= maxHex; hex++) {
                float yMass = hexnac * GlycanMasses.hexnacMass + hex * GlycanMasses.hexMass;
                boolean skipHexOnly = maxHexNAc > 0 && hexnac == 0;    // skip Ys that have Hex but no HexNAc, as long as HexNAc is part of the composition
                if (yMass > 0 && !skipHexOnly) {
                    expectedYIons.add(yMass);
                    // Fuc: allow Y ions with and without Fuc if present
                    if (containsResidueType(GlycanResidue.dHex)) {
                        expectedYIons.add(yMass + GlycanMasses.dhexMass);
                    }
                }
            }
        }
        /* generate disallowed Ys: more Hex, HexNAc than in this composition.
         * The purpose of this is to check for alternative explanations (i.e. if there is a clear series of Y ions
         * that extend beyond this the max in this candidate, it is more likely incorrect than if the observed Y ion
         * series terminates where expected based on this composition).
         * Currently hardcoded how many extras to check. todo: might need tuning/improvement in future
         */
        int maxExtraHexNAc = 3;
        int maxExtraHex = 5;
        for (int hexnac = maxHexNAc; hexnac <= maxHexNAc + maxExtraHexNAc; hexnac++) {
            for (int hex = maxHex; hex <= maxHex + maxExtraHex; hex++) {
                // don't want to add an allowed combination, but do need to start loops at 0 extra to cover all possibilities - skip this case
                if (hexnac == maxHexNAc && hex == maxHex)
                    continue;
                float yMass = hexnac * GlycanMasses.hexnacMass + hex * GlycanMasses.hexMass;
                disallowedYIons.add(yMass);
            }
        }
    }

    /**
     * Set expected/disallowed Oxonium ions for this candidate based on a set of heuristic
     * filtering rules. Expected and disallowed ions are listed as 1+ m/z.
     */
    private void initializeAllowedOxoniums() {
        expectedOxoniumIons = new ArrayList<>();
        disallowedOxoniumIons = new ArrayList<>();

        // NeuAc: expect NeuAc Oxo ions if present and vice versa
        if (containsResidueType(GlycanResidue.NeuAc)) {
            expectedOxoniumIons.addAll(Arrays.asList(GlycanMasses.NeuAcOxoniums));
        } else {
            disallowedOxoniumIons.addAll(Arrays.asList(GlycanMasses.NeuAcOxoniums));
        }

        // NeuGc: expect NeuGc Oxo ions if present and vice versa
        if (containsResidueType(GlycanResidue.NeuGc)) {
            expectedOxoniumIons.addAll(Arrays.asList(GlycanMasses.NeuGcOxoniums));
        } else {
            disallowedOxoniumIons.addAll(Arrays.asList(GlycanMasses.NeuGcOxoniums));
        }

        // Phospho: expect if present and vice versa (todo: not always expect if present?)
        if (containsResidueType(GlycanResidue.Phospho)) {
            expectedOxoniumIons.addAll(Arrays.asList(GlycanMasses.PhosphoOxoniums));
        } else {
            disallowedOxoniumIons.addAll(Arrays.asList(GlycanMasses.PhosphoOxoniums));
        }

        // Sulfo: expect if present and vice versa (todo: not always expect if present?)
        if (containsResidueType(GlycanResidue.Sulfo)) {
            expectedOxoniumIons.addAll(Arrays.asList(GlycanMasses.SulfoOxoniums));
        } else {
            disallowedOxoniumIons.addAll(Arrays.asList(GlycanMasses.SulfoOxoniums));
        }
    }

    /**
     * Generate an array of all m/z values expected for the Y ions of this candidate given the max
     * charge (typically precursor charge - 1). Boolean "expected" means return expected Ys if true,
     * disallowed Ys if false.
     * @param maxCharge max Y ion charge
     * @param expected expected or disallowed ions
     * @return float[] of m/z values to search
     */
    public float[] generateYmzs(int maxCharge, boolean expected) {
        ArrayList<Float> ionsToUse;
        if (expected) {
            ionsToUse = expectedYIons;
        } else {
            ionsToUse = disallowedYIons;
        }
        float[] yMZs = new float[ionsToUse.size() * maxCharge];
        for (int i=0; i < ionsToUse.size(); i++) {
            for (int z=1; z <= maxCharge; z++) {
                yMZs[i+z] = Spectrum.neutralMassToMZ(ionsToUse.get(i), z);
            }
        }
        return yMZs;
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
