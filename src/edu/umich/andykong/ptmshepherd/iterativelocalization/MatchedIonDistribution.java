package edu.umich.andykong.ptmshepherd.iterativelocalization;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;

public class MatchedIonDistribution {
    int[] pdf;
    int[] pdfDecoy;
    double[] cdf;

    int[] pdfMassError;
    int[] pdfMassErrorDecoy;
    double[] ionPosterior;
    double[] ionPosteriorMassError;
    float resolution;
    int binMult;

    public String mode = "ptminer-unmatched-theoretical";

    public MatchedIonDistribution(float resolution, boolean poissonBinomial) {
        this.resolution = resolution;
        this.binMult = (int) (1.0f/resolution);
        this.pdf = new int[((int) (100.0 * this.binMult + 1))];
        this.pdfDecoy = new int[((int) (100.0 * this.binMult + 1))];
        this.pdfMassError = new int[100]; // Default 100 ppm max, can make it dynamic TODO
        this.pdfMassErrorDecoy = new int[100]; // Default 100 ppm max, can make it dynamic TODO
        if (poissonBinomial)
            this.mode = "poissonbinomial-matched-theoretical-decoy";
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
        } else if (this.mode.equals("poissonbinomial-matched-theoretical-decoy")) {
            if (intensity < 0.0f) {
                int idx = (int) (-1 * intensity * this.binMult);
                this.pdfDecoy[idx]++;
            } else {
                int idx = (int) (-1 * intensity * this.binMult);
                this.pdf[idx]++;
            }
        }
    }

    /**
     * Adds ion to ion posterior probability distribution. Currently only using MATCHED ions.
     * @param intensity
     * @param isDecoy
     */
    public void addIon(float intensity, boolean isDecoy) {
        if (!isDecoy) {
            if (intensity > 0.0f) {
                int idx = (int) (intensity * this.binMult);
                this.pdf[idx]++;
            }
        } else {
            if (intensity > 0.0f) {
                int idx = (int) (intensity * this.binMult);
                this.pdfDecoy[idx]++;
            }
        }
    }

    /**
     * Adds ion to ion posterior probability distribution. Currently only using MATCHED ions.
     * @param intensity
     * @param massError
     * @param isDecoy
     */
    public void addIon(float intensity, float massError, boolean isDecoy) {
        if (!isDecoy) {
            if (intensity > 0.0f) { // Only include matched ions, unmatched ions have negative intensity
                int intIndx = (int) (intensity * this.binMult);
                this.pdf[intIndx]++;
                int massErrorIndx = (int) massError;
                this.pdfMassError[massErrorIndx]++;
            }
        } else {
            if (intensity > 0.0f) { // Only include matched ions, unmatched ions have negative intensity
                int idx = (int) (intensity * this.binMult);
                this.pdfDecoy[idx]++;
                int massErrorIndx = (int) massError;
                this.pdfMassErrorDecoy[massErrorIndx]++;
            }
        }
    }


    public void addIons(float [] intensities) {
        for(int i = 0; i < intensities.length; i++) {
            addIon(intensities[i]);
        }
    }

    public void addIons(float[] intensities, float[] massErrors, boolean isDecoy) {
        for(int i = 0; i < intensities.length; i++) {
            addIon(intensities[i], massErrors[i], isDecoy);
        }
    }

    public void addIons(float [] intensities, boolean isDecoy) {
        for(int i = 0; i < intensities.length; i++) {
            addIon(intensities[i], isDecoy);
        }
    }

    public void addIon_Original(float intensity) {
        int idx = (int) (intensity * this.binMult);
        this.pdf[idx]++;
    }

    public void addIons_Original(float [] intensities) {
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

    public void calculateIonPosterior() {
        calculateIonIntensityPosterior();
        calculateIonMassErrorPosterior();
    }

    private void calculateIonIntensityPosterior() {
        // Calculate local q-val estimate //todo terminology
        this.ionPosterior = new double[((int) (100.0 * this.binMult + 1))];
        int sumTarget = 0;
        int sumDecoy = 1;
        for (int i = this.pdf.length - 1; i >= 0; i--) {
            sumTarget += this.pdf[i];
            sumDecoy += this.pdfDecoy[i];
            this.ionPosterior[i] = 1.0 - ((double) sumDecoy / Math.max(sumTarget, 1)); //todo try decoy+1/target
        }

        // Find q-val monotonic increasing
        // leftStop is always 0
        int rightStop = this.ionPosterior.length;
        while (rightStop >= 0) {
            double max = -1;
            int maxIndx = -1;
            // Search within window for max
            for (int i = 0; i < rightStop; i++) {
                if (this.ionPosterior[i] > max) {
                    max = this.ionPosterior[i];
                    maxIndx = i;
                }
            }
            System.out.println(rightStop);
            if (maxIndx == -1) //Why??? TODO
                break;
            for (int i = maxIndx; i < rightStop; i++) {
                this.ionPosterior[i] = max;
            }
            rightStop = maxIndx - 1;
        }
    }

    private void calculateIonMassErrorPosterior() {
        // Calculate local q-val estimates ?? //todo terminology
        this.ionPosteriorMassError = new double[100];
        int sumTarget = 0;
        int sumDecoy = 1;
        for (int i = 0; i < this.pdfMassError.length; i++) {
            sumTarget += this.pdfMassError[i];
            sumDecoy += this.pdfMassErrorDecoy[i];
            this.ionPosteriorMassError[i] = 1.0 - ((double) sumDecoy / Math.max(sumTarget, 1)); //todo try decoy+1/target
        }

        // Find q-val monotonic increasing
        // rightStop is always this.ionPosteriorMassError.length;
        int leftStop = 0;
        while (leftStop <= this.ionPosteriorMassError.length - 1) {
            double max = -1;
            int maxIndx = -1;
            // Search within window for max
            for (int i = this.ionPosteriorMassError.length - 1; i > leftStop; i--) {
                if (this.ionPosteriorMassError[i] > max) {
                    max = this.ionPosteriorMassError[i];
                    maxIndx = i;
                }
            }
            for (int i = maxIndx; i > leftStop; i--) {
                this.ionPosteriorMassError[i] = max;
            }
            leftStop = maxIndx + 1;
            System.out.println(leftStop);
        }
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
            if (!this.mode.equals("poissonbinomial-matched-theoretical-decoy")) {
                prob = this.cdf[(int) intensity * this.binMult];
            } else {
                prob = this.ionPosterior[(int) intensity * this.binMult];
            }
        return prob;
    }

    public double[] calcIonProbabilities(float [] intensities) {
        double[] probs = new double[intensities.length];
        for (int i = 0; i < intensities.length; i++)
            probs[i] = calcIonProbability(intensities[i]);
        return probs;
    }

    public String intensityToString() {
        StringBuffer sb = new StringBuffer();

        if (!this.mode.equals("poissonbinomial-matched-theoretical-decoy")) {
            sb.append("intensity\tcount\n");
            for (int i = 0; i < this.pdf.length; i++) {
                sb.append(new DecimalFormat("0.00").format((double) i / (double) this.binMult)
                        + "\t" + this.pdf[i] + "\n");
            }
        } else {
            sb.append("intensity\tcount_target\tcount_decoy\tion_PEP\n");
            for (int i = 0; i < this.pdf.length; i++) {
                sb.append(new DecimalFormat("0.00").format((double) i / (double) this.binMult)
                        + "\t" + this.pdf[i] + "\t" + this.pdfDecoy[i] + "\t" + this.ionPosterior[i] + "\n");
            }
        }
        return sb.toString();
    }

    public void printIntensityHisto(String fname) throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(fname));
        out.print(this.intensityToString());
        out.close();
    }

    public String massErrorToString() {
        StringBuffer sb = new StringBuffer();

        sb.append("mass_error\tcount_target\tcount_decoy\tion_PEP\n");
        for (int i = 0; i < this.pdfMassError.length; i++) {sb.append(i + "\t" + this.pdfMassError[i] + "\t" +
                this.pdfMassErrorDecoy[i] + "\t" + this.ionPosteriorMassError[i] + "\n");
        }
        return sb.toString();
    }

    public void printMassErrorHisto(String fname) throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(fname));
        out.print(this.massErrorToString());
        out.close();
    }

}
