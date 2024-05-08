package edu.umich.andykong.ptmshepherd.utils;

import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.AAMutationRules;
import org.jetbrains.annotations.NotNull;

import java.nio.charset.StandardCharsets;
import java.util.*;

public class Peptide { //TODO theoretical peptide fragments, should this not start at 0? Skip b1/y1/a1?
    public String pepSeq;
    public float[] mods;
    public int mutatedResidue; // for mono-mutated decoys, 0 indexed

    public Peptide(String pepSeq, float[] mods) {
        this.pepSeq = pepSeq;
        this.mods = mods;
    }

    private Peptide(String pepSeq, float[] mods, int mutatedResidue) {
        this(pepSeq, mods);
        this.mutatedResidue = mutatedResidue;
    }

    // Add mod, 0-index
    public void addMod(float dmass, int residue) {
        this.mods[residue] += dmass;
    }

    public ArrayList<Float> calculatePeptideFragments(String ionTypes, int maxCharge) {
       return  calculatePeptideFragments(this.pepSeq, this.mods, ionTypes, maxCharge);
    }

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
                float cmass = AAMasses.protMass + nTermMass;
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
                float cmass = (cTermMass + ccharge * AAMasses.protMass) / ccharge;
                for (int i = 0; i < cLen - 1; i++) {
                    cmass += (aaMasses[seq.charAt(cLen - 1 - i) - 'A'] + mods[cLen - 1 - i]) / ccharge;
                    knownFrags.add(cmass);
                }
            }
        }

        return knownFrags;
    }

    public Peptide generateDecoy(Random rng, String method) {
        if (method.equals("shuffled"))
            return generateShuffledDecoy(this.pepSeq, this.mods, rng);
        else if (method.equals("full-shuffled"))
            return generateFullShuffledDecoy(this.pepSeq, this.mods, rng);
        else if (method.equals("swapped"))
            return generateSwappedDecoy(this.pepSeq, this.mods, rng);
        else if (method.equals("mutated"))
            return generateMutatedDecoy(this.pepSeq, this.mods);
        else if (method.equals("mono-mutated"))
            return generateMonoMutatedDecoy(this.pepSeq, mods, rng);
        else
            return null; //TODO
    }

    public static Peptide generateDecoy(String pep, float[] mods, Random rng, String method) {
        if (method.equals("shuffled"))
            return generateShuffledDecoy(pep, mods, rng);
        else if (method.equals("full-shuffled"))
            return generateFullShuffledDecoy(pep, mods, rng);
        else if (method.equals("swapped"))
            return generateSwappedDecoy(pep, mods, rng);
        else if (method.equals("mutated"))
            return generateMutatedDecoy(pep, mods);
        else if (method.equals("mono-mutated"))
            return generateMonoMutatedDecoy(pep, mods, rng);
        else
            return null; //TODO
    }

    public static Peptide generateDecoy(String pep, float[] mods, int pos, Random rng, String method) {
        if (method.equals("mono-swapped"))
            return generateMonoSwappedDecoy(pep, mods, pos, rng);
        else
            return null; //TODO
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

    public static Peptide generateFullShuffledDecoy(String pep, float[] mods, Random rng) {
        ArrayList<Site> sites = new ArrayList<>(pep.length());

        // Shuffle pep
        for (int i = 0; i < pep.length(); i++)
            sites.add(new Site(pep.charAt(i), mods[i]));
        Collections.shuffle(sites, rng);

        StringBuilder newPep = new StringBuilder();
        float[] newMods = new float[mods.length];

        for (int i = 0; i < sites.size(); i++) {
            newPep.append(sites.get(i).aa);
            newMods[i] = sites.get(i).mod;
        }

        return new Peptide(newPep.toString(), newMods);
    }

    public static Peptide generateMutatedDecoy(String pep, float[] mods) {
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

    public static Peptide generateMonoMutatedDecoy(String pep, float[] mods, Random rng) {
        int randomSite = rng.nextInt(pep.length());
        return generateMonoMutatedDecoy(pep, mods, randomSite);
    }

    public static Peptide generateMonoMutatedDecoy(String pep, float[] mods, int mutSite) {
        StringBuilder newPep = new StringBuilder();
        for (int i = 0; i < pep.length(); i++) {
            if (i != mutSite)
                newPep.append(pep.charAt(i));
            else
                newPep.append(AAMutationRules.mutateFrom.get(pep.charAt(i)));
        }
        return new Peptide(newPep.toString(), mods, mutSite);
    }

    public static Peptide generateSwappedDecoy(String pepSeq, float[] mods, Random rng) {
        int pos = rng.nextInt(pepSeq.length());
        int swapPos = pos;
        while (swapPos == pos)
            swapPos = rng.nextInt(pepSeq.length());

        char swapAA = pepSeq.charAt(swapPos);
        float swapMod = mods[swapPos];

        char[] pepSeqArray = pepSeq.toCharArray();
        pepSeqArray[swapPos] = pepSeqArray[pos];
        mods[swapPos] = mods[pos];

        pepSeqArray[pos] = swapAA;
        mods[pos] = swapMod;

        return new Peptide(pepSeqArray.toString(), mods);
    }

    public static Peptide generateMonoSwappedDecoy(String pepSeq, float[] mods, int pos, Random rng) {
        int swapPos = pos;
        while (Math.abs(swapPos - pos) < (pepSeq.length() / 2))
            swapPos = rng.nextInt(pepSeq.length());

        char swapAA = pepSeq.charAt(swapPos);
        float swapMod = mods[swapPos];

        char[] pepSeqArray = pepSeq.toCharArray();
        pepSeqArray[swapPos] = pepSeqArray[pos];
        mods[swapPos] = mods[pos];

        pepSeqArray[pos] = swapAA;
        mods[pos] = swapMod;

        return new Peptide(String.valueOf(pepSeqArray), mods);
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
