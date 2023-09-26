package edu.umich.andykong.ptmshepherd.iterativelocalization;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;

public class MatchedIonDistribution {
    int[] pdf;
    int[] pdfUnmatched;
    double[] cdf;
    float resolution;
    int binMult;

    public final String mode = "ptminer-unmatched-theoretical";

    public MatchedIonDistribution(float resolution) {
        this.resolution = resolution;
        this.binMult = (int) (1.0f/resolution);
        this.pdf = new int[((int) (100.0 * this.binMult + 1))];
        this.pdfUnmatched = new int[((int) (100.0 * this.binMult + 1))]; //TODO test
    }

    // Original was ptminer mode
    public void addIon(float intensity) {
        if (this.mode.equals("ptminer")) {
            if (intensity < 0.0f)
                return;
            int idx = (int) (intensity * this.binMult);
            this.pdf[idx]++;
        } else if (this.mode.equals("ptminer-unmatched-theoretical")) {
            if (intensity < 0.0f) {//include unmatched ions in the distribution
                int idx = (int) (0.0);
                this.pdf[idx]++;
            } else {
                int idx = (int) (intensity * this.binMult);
                this.pdf[idx]++;
            }
        } else if (this.mode.equals("ptminer-unmatched-theoretical-experimental")) {
            if (intensity < 0.0f) {
                int idx = (int) (-1 * intensity * this.binMult);
                this.pdfUnmatched[idx]++;
            } else {
                int idx = (int) (-1 * intensity * this.binMult);
                this.pdf[idx]++;
            }
        }
    }

    public synchronized void addIons(float [] intensities) {
        for(int i = 0; i < intensities.length; i++) {
            addIon(intensities[i]);
        }
    }

    public void addIon_Original(float intensity) {
        int idx = (int) (intensity * this.binMult);
        this.pdf[idx]++;
    }

    public synchronized void addIons_Original(float [] intensities) {
        for(int i = 0; i < intensities.length; i++) {
             if (intensities[i] < 0)
                continue;
             addIon(intensities[i]);
        }
    }

    public void calculateCdf() {
        this.cdf = new double[((int) (100.0 * this.binMult + 1))];
        int sum = 0;
        for (int i = 0; i < this.pdf.length; i++) {
            sum += this.pdf[i];
            this.cdf[i] = sum;
        }
        for (int i = 0; i < cdf.length; i++)
            this.cdf[i] /= (double) sum;

        // Prevent bad left tail behavior by extending (min value/2) to 0
        int minI = 0;
        while (this.cdf[minI] == 0.0)
            minI++;
        double minVal = this.cdf[minI] / 2.0;
        for (int i = 0; i < minI; i++)
            this.cdf[i] = minVal;

        // Prevent bad right tail behavior by setting final two bins equal to their weighted average
        double lastBinsAverage = (this.cdf[this.cdf.length-2] + this.cdf[this.cdf.length-1]) / 2.0;
        this.cdf[this.cdf.length-1] = this.cdf[this.cdf.length-2] = lastBinsAverage;
    }

    public void calculateCdf_UnmatchedExperimentalTheoretical() { //TODO test
        this.cdf = new double[((int) (100.0 * this.binMult + 1))];
        int sum = 0;
        for (int i = 0; i < this.pdf.length; i++) {
            sum += this.pdf[i];
            this.cdf[i] = sum;
        }
        for (int i = 0; i < cdf.length; i++)
            this.cdf[i] /= (double) sum;

        // Prevent bad left tail behavior by extending (min value/2) to 0
        int minI = 0;
        while (this.cdf[minI] == 0.0)
            minI++;
        double minVal = this.cdf[minI] / 2.0;
        for (int i = 0; i < minI; i++)
            this.cdf[i] = minVal;

        // Prevent bad right tail behavior by setting final two bins equal to their weighted average
        double lastBinsAverage = (this.cdf[this.cdf.length-2] + this.cdf[this.cdf.length-1]) / 2.0;
        this.cdf[this.cdf.length-1] = this.cdf[this.cdf.length-2] = lastBinsAverage;
    }

    public double calcIonProbability(float intensity) {
        double prob;
        if (intensity < 0)
             prob = -1;
        else
            prob = this.cdf[(int) intensity * this.binMult];
        return prob;
    }

    public double[] calcIonProbabilities(float [] intensities) {
        double[] probs = new double[intensities.length];
        for (int i = 0; i < intensities.length; i++)
            probs[i] = calcIonProbability(intensities[i]);
        return probs;
    }
}
