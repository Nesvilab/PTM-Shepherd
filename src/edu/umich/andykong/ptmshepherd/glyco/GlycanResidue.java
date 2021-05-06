package edu.umich.andykong.ptmshepherd.glyco;

/**
 * Allowed names for glyco residues. Matched to masses in GlycanMasses
 */

public enum GlycanResidue {
    HexNAc,     // N-acetyl hexosamine (e.g. GlcNAc)
    Hex,        // Hexose (e.g. Glucose)
    dHex,       // deoxy Hexose (e.g. Fucose)
    NeuAc,      // N-acetyl Neuraminic acid (sialic acid)
    NeuGc,      // N-Glycolyl Neuraminic acid (sialic acid)
    Phospho,    // Phosphate
    Sulfo,      // Sulfate
    NH3,        // ammonium adduct
    Na,         // sodium adduct
    Fe3,        // iron-3 adduct
    Fe2,        // iron-2 adduct
    Ca,         // calcium adduct
    Al          // aluminum adduct
}
