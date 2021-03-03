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

    // lists of expected oxonium ions generated by given residue types
    public static final Float[] NeuAcOxoniums = {274.0921325f, 292.1026925f, 657.2349f};
    public static final Float[] NeuGcOxoniums = {308.09761f};
    public static final Float[] PhosphoOxoniums = {243.026426f, 405.079246f, 485.045576f};
    public static final Float[] SulfoOxoniums = {284.044f, 446.098f, 811.212f};

}

