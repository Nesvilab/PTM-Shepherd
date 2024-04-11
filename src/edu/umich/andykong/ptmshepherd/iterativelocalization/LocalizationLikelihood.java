package edu.umich.andykong.ptmshepherd.iterativelocalization;

import java.util.ArrayList;

public class LocalizationLikelihood {
    ArrayList<Mod> mods;

    LocalizationLikelihood() {
        this.mods = new ArrayList<>();
    }

    LocalizationLikelihood(float dMass, double[] siteLikelihoods) {
        this.mods = new ArrayList<>();
        this.mods.add(new Mod(dMass, siteLikelihoods));
    }

    public Mod getMod() { //todo variable mod searches will require multiple modifications
        return mods.get(0);
    }

    public class Mod {
        float dMass;
        double[] siteLikelihoods;

        Mod(float dMass, double[] siteLikelihoods) {
            this.dMass = dMass;
            this.siteLikelihoods = siteLikelihoods;
        }

        public double[] getSiteLikelihoods() {
            return this.siteLikelihoods;
        }
    }
}
