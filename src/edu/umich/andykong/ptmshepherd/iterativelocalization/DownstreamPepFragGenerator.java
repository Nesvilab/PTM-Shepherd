package edu.umich.andykong.ptmshepherd.iterativelocalization;

import edu.umich.andykong.ptmshepherd.utils.Peptide;

import java.util.ArrayList;

public class DownstreamPepFragGenerator {
    public static ArrayList<Float> calculatePeptideFragments(Peptide pep, String ionTypes, int fromResidue, int maxCharge) {
        ArrayList<Float> pepFrags = new ArrayList<>();

        ArrayList<Float> tmpPepFrags;
        int startIndex = fromResidue;
        for (Character it : ionTypes.toCharArray()) {
            tmpPepFrags = pep.calculatePeptideFragments(Character.toString(it), maxCharge);
            if (it == 'a' || it == 'b' || it == 'c') { //if n-term ion series {
                if (fromResidue < pep.pepSeq.length()-1) {
                    if (fromResidue != 0) // fragmentation doesn't produce a/b/c_1 ion
                        startIndex--;
                    pepFrags.addAll(tmpPepFrags.subList(startIndex, tmpPepFrags.size()));
                }
            } else if (it == 'x' || it == 'y' || it == 'z') { // c-term ion series
                if (fromResidue > 0) {
                    pepFrags.addAll(tmpPepFrags.subList(pep.pepSeq.length()-(fromResidue+1), tmpPepFrags.size()));
                }
            }
        }

        return pepFrags;
    }
}