/*
 *    Copyright 2022 University of Michigan
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

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
