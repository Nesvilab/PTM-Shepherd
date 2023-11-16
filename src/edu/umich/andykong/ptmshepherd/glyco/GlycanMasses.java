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

import java.util.HashMap;
import java.util.TreeMap;

public class GlycanMasses {

    // Building blocks for glycan compositions
    public static final float hexnacMass = (float) 203.07937;
    public static final float hexMass = (float) 162.05282;
    public static final float dhexMass = (float) 146.057909;
    public static final float neuacMass = (float) 291.095417;
    public static final float neugcMass = (float) 307.090334;
    public static final float phosphoMass = (float) 79.96633;
    public static final float sulfoMass = (float) 79.95682;

    // Adduct masses
    public static final float NH3adductMass = (float) 17.02655;
    public static final float NaAdductMass = (float) 21.981944;     // replaced 1 proton by Na+
    public static final float Fe3AdductMass = (float) 52.911462;     // replaced 3 protons by Fe+3
    public static final float Fe2AdductMass = (float) 53.919288;     // replaced 2 protons by Fe+2
    public static final float CaAdductMass = (float) 37.946938;     // replaced 2 protons by Ca+2
    public static final float AlAdductMass = (float) 23.958062;     // replaced 3 protons by Al+3


    // Map of names to residues for parsing database files. NOTE: assumed all lower case
    public static TreeMap<String, GlycanResidue> glycoNames;
    static
    {
        glycoNames = new TreeMap<>();
        glycoNames.put("phospho", GlycanResidue.Phospho);
        glycoNames.put("phosphorylation", GlycanResidue.Phospho);
        glycoNames.put("phosphate", GlycanResidue.Phospho);
        glycoNames.put("sulf", GlycanResidue.Sulfo);
        glycoNames.put("sulfo", GlycanResidue.Sulfo);
        glycoNames.put("sulfation", GlycanResidue.Sulfo);
        glycoNames.put("sulfate", GlycanResidue.Sulfo);
        glycoNames.put("NH3", GlycanResidue.NH3);
        glycoNames.put("nh3", GlycanResidue.NH3);
        glycoNames.put("Nh3", GlycanResidue.NH3);
        glycoNames.put("ammonia", GlycanResidue.NH3);
        glycoNames.put("Na", GlycanResidue.Na);
        glycoNames.put("Na+", GlycanResidue.Na);
        glycoNames.put("na", GlycanResidue.Na);
        glycoNames.put("sodium", GlycanResidue.Na);
        glycoNames.put("Fe3", GlycanResidue.Fe3);
        glycoNames.put("Fe3+", GlycanResidue.Fe3);
        glycoNames.put("fe3", GlycanResidue.Fe3);
        glycoNames.put("fe3+", GlycanResidue.Fe3);
        glycoNames.put("Fe2", GlycanResidue.Fe2);
        glycoNames.put("Fe2+", GlycanResidue.Fe2);
        glycoNames.put("fe2", GlycanResidue.Fe2);
        glycoNames.put("fe2+", GlycanResidue.Fe2);
        glycoNames.put("Ca", GlycanResidue.Ca);
        glycoNames.put("ca", GlycanResidue.Ca);
        glycoNames.put("Al", GlycanResidue.Al);
        glycoNames.put("al", GlycanResidue.Al);
    }
    public static HashMap<GlycanResidue, String> outputGlycoNames;
    static
    {
        outputGlycoNames = new HashMap<>();
        outputGlycoNames.put(GlycanResidue.Hex, "Hex");
        outputGlycoNames.put(GlycanResidue.HexNAc, "HexNAc");
        outputGlycoNames.put(GlycanResidue.dHex, "Fuc");
        outputGlycoNames.put(GlycanResidue.NeuAc, "NeuAc");
        outputGlycoNames.put(GlycanResidue.NeuGc, "NeuGc");
        outputGlycoNames.put(GlycanResidue.Phospho, "Phospho");
        outputGlycoNames.put(GlycanResidue.Sulfo, "Sulfo");
        outputGlycoNames.put(GlycanResidue.NH3, "NH3");
        outputGlycoNames.put(GlycanResidue.Na, "Na");
        outputGlycoNames.put(GlycanResidue.Fe3, "Fe3");
        outputGlycoNames.put(GlycanResidue.Fe2, "Fe2");
        outputGlycoNames.put(GlycanResidue.Ca, "Ca");
        outputGlycoNames.put(GlycanResidue.Al, "Al");
    }

}

