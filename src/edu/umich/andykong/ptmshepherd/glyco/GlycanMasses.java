package edu.umich.andykong.ptmshepherd.glyco;

import java.util.TreeMap;

public class GlycanMasses {

    // Building blocks for glycan compositions
    public static final float hexnacMass = (float) 203.079037;
    public static final float hexMass = (float) 162.05282;
    public static final float dhexMass = (float) 146.057909;
    public static final float neuacMass = (float) 291.095417;
    public static final float neugcMass = (float) 307.090334;
    public static final float phosphoMass = (float) 79.96633;
    public static final float sulfoMass = (float) 79.95682;

    // map of residue -> mass
    public static TreeMap<GlycanResidue, Float> glycoMasses;
    static
    {
        glycoMasses = new TreeMap<GlycanResidue, Float>();
        glycoMasses.put(GlycanResidue.Hex, hexMass);
        glycoMasses.put(GlycanResidue.HexNAc, hexnacMass);
        glycoMasses.put(GlycanResidue.dHex, dhexMass);
        glycoMasses.put(GlycanResidue.NeuAc, neuacMass);
        glycoMasses.put(GlycanResidue.NeuGc, neugcMass);
        glycoMasses.put(GlycanResidue.Phospho, phosphoMass);
        glycoMasses.put(GlycanResidue.Sulfo, sulfoMass);
    }

}

