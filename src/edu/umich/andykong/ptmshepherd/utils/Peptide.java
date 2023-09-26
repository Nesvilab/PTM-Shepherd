package edu.umich.andykong.ptmshepherd.utils;

import edu.umich.andykong.ptmshepherd.core.AAMasses;

import java.util.ArrayList;

public class Peptide { //TODO theoretical peptide fragments, should this not start at 0? Skip b1/y1/a1?
    public static ArrayList<Float> calculatePeptideFragments(String seq, float[] mods, String ionTypes, int maxCharge) {
        ArrayList<Float> knownFrags = new ArrayList<>(seq.length() * ionTypes.length());

        ArrayList<Character> nIonTypes = new ArrayList<>();
        ArrayList<Character> cIonTypes = new ArrayList<>();
        for (int i = 0; i < ionTypes.length(); i++) {
            char curIonType = ionTypes.charAt(i);
            if (curIonType == 'a' || curIonType == 'b' || curIonType == 'c')
                nIonTypes.add(curIonType);
            else if (curIonType == 'x' || curIonType == 'y' || curIonType == 'z')
                cIonTypes.add(curIonType);
        }

        float [] aaMasses = AAMasses.monoisotopic_masses;
        float [] fragTypeShifts = AAMasses.ionTypeShifts;
        int cLen = seq.length();

        float nTermMass;
        for (Character iType : nIonTypes) {
            nTermMass = fragTypeShifts[iType - 'a'];
            for (int ccharge = 1; ccharge <= maxCharge; ccharge++) { //loop through charge states
                float cmass = AAMasses.monoisotopic_nterm_mass + nTermMass;
                for (int i = 0; i < cLen - 1; i++) { //loop through positions on the peptide
                    cmass += (aaMasses[seq.charAt(i) - 'A'] + mods[i]) / ccharge;
                    knownFrags.add(cmass);
                }
            }
        }
        float cTermMass;
        for (Character iType : cIonTypes) {
            cTermMass = fragTypeShifts[iType - 'x' + 3];
            for (int ccharge = 1; ccharge <= maxCharge; ccharge++) {
                float cmass = (cTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
                for (int i = 0; i < cLen - 1; i++) {
                    cmass += (aaMasses[seq.charAt(cLen - 1 - i) - 'A'] + mods[cLen - 1 - i]) / ccharge;
                    knownFrags.add(cmass);
                }
            }
        }

        return knownFrags;
    }
}
