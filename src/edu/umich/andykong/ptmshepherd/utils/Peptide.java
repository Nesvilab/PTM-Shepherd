package edu.umich.andykong.ptmshepherd.utils;

import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.AAMutationRules;
import org.jetbrains.annotations.NotNull;

import java.nio.charset.StandardCharsets;
import java.util.*;

public class Peptide { //TODO theoretical peptide fragments, should this not start at 0? Skip b1/y1/a1?
    public String pepSeq;
    public float[] mods;
    public Peptide(String pepSeq, float[] mods) {
        this.pepSeq = pepSeq;
        this.mods = mods;
    }

    public static ArrayList<Float> calculatePeptideFragments(String seq, float[] mods, String ionTypes, int maxCharge) {
        ArrayList<Float> knownFrags = new ArrayList<>(seq.length() * ionTypes.length());

        ArrayList<Character> nIonTypes = new ArrayList<>();
        ArrayList<Character> cIonTypes = new ArrayList<>();
        // Todo this should start at b2
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
                    if (i != 0) // todo skip a1/b1 ion, check if this needs to be skipped for c too
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

    public static Peptide generateDecoy(String pep, float[] mods, Random rng, String method) {
        if (method.equals("shuffled"))
            return generateShuffledDecoy(pep, mods, rng);
        else if (method.equals("mutated"))
            return generateMutatedDecoy(pep, mods, rng);
        else
            return null; //TODO
    }

    public static Peptide generateMutatedDecoy(String pep, float[] mods, Random rng) {
        ArrayList<Site> sites = new ArrayList<>(pep.length());
        for (int i = 0; i < pep.length(); i++)
            sites.add(new Site(pep.charAt(i), mods[i]));
        sites.get(1).aa = AAMutationRules.mutateFrom.get(sites.get(1).aa);
        sites.get(pep.length()-2).aa = AAMutationRules.mutateFrom.get(sites.get(pep.length()-2).aa);

        StringBuilder newPep = new StringBuilder();

        for (int i = 0; i < sites.size(); i++)
            newPep.append(sites.get(i).aa);

        return new Peptide(newPep.toString(), mods);
    }

    public static Peptide generateShuffledDecoy(String pep, float[] mods, Random rng) {
        ArrayList<Site> sites = new ArrayList<>(pep.length());

        // Shuffle core
        for (int i = 1; i < pep.length()-1; i++)
            sites.add(new Site(pep.charAt(i), mods[i]));
        Collections.shuffle(sites, rng);

        StringBuilder newPep = new StringBuilder();
        float[] newMods = new float[mods.length];
        // N-term AA
        newPep.append(pep.charAt(0));
        newMods[0] = mods[0];
        for (int i = 0; i < sites.size(); i++) {
            newPep.append(sites.get(i).aa);
            newMods[i+1] = sites.get(i).mod;
        }
        // C-term AA
        newPep.append(pep.charAt(pep.length()-1));
        newMods[newMods.length-1] = mods[mods.length-1];

        return new Peptide(newPep.toString(), newMods);
    }

    static class Site {
        char aa;
        float mod;

        Site(char aa, float mod) {
            this.aa = aa;
            this.mod = mod;
        }
    }
}
