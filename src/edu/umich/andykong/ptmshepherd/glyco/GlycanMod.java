package edu.umich.andykong.ptmshepherd.glyco;

import edu.umich.andykong.ptmshepherd.PTMShepherd;

import java.util.ArrayList;

public class GlycanMod extends GlycanResidue {
    public boolean isFixed;
    public int maxAllowed;
    public ArrayList<GlycanResidue> requiredResidues;

    /**
     * Parse the basic information about a mod from the glycan_mods.tsv file. Mod settings (max allowed, fixed/variable)
     * come from user params and are set later.
     * @param modTsvLine
     */
    public GlycanMod(String modTsvLine, GlycoParams glycoParams) {
        super(modTsvLine);
        String[] splits = modTsvLine.replace("\"", "").split("\t");
        maxAllowed = Integer.parseInt(splits[6]);
        isFixed = Boolean.parseBoolean(splits[7]);

        // parse required residues
        requiredResidues = new ArrayList<>();
        if (!splits[5].matches("")) {
            String[] resSplits = splits[5].split(",");
            for (String residueStr : resSplits) {
                GlycanResidue residue = glycoParams.findResidueName(residueStr);
                if (residue != null) {
                    requiredResidues.add(residue);
                } else {
                    PTMShepherd.die(String.format("Error parsing glycan_mods.tsv: required residue %s in line %s was not recognized. Make sure this glycan residue is in the glycan_residues.tsv and that it is spelled correctly in glycan_mods.tsv and try again.", residueStr, modTsvLine));
                }
            }
        }
    }

    public String printParam() {
        StringBuilder requiredRes = new StringBuilder();
        if (requiredResidues.size() > 0) {
            requiredRes.append(", requires:");
            for (GlycanResidue residue: requiredResidues) {
                requiredRes.append(" ");
                requiredRes.append(residue.name);
            }
        }

        if (isFixed) {
            return String.format("fixed: %s%s", super.printParam(), requiredRes);
        } else if (maxAllowed > 0) {
            return String.format("variable: %s, max %d%s", super.printParam(), maxAllowed, requiredRes);
        } else {
            return "";
        }
    }
}
