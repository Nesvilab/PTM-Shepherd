package edu.umich.andykong.ptmshepherd.iterativelocalization;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import static edu.umich.andykong.ptmshepherd.PTMShepherd.executorService;


public class MatchedIonDistribution {
    int[] pdf;
    int[] pdfDecoy;
    double[] cdf;
    double[] qVals;

    int dNegCount;
    int tNegCount;
    double negPredictiveValue;

    int[] projTargetCounts;
    int[] projDecoyCounts;

    int[] pdfMassError;
    int[] pdfMassErrorDecoy;
    double[] ionPosterior;
    double[] ionPosteriorMassError;
    float resolution;
    int binMult;
    MatchedIonTable datapoints;





    TwoDimJointPMF targetPMF;
    TwoDimJointPMF decoyPMF;
    TwoDimJointQValue qValues;
    LDAProcessor ldaProcessor;

    public String mode = "ptminer-unmatched-theoretical";

    public MatchedIonDistribution(float resolution, boolean poissonBinomial) {
        this.resolution = resolution;
        this.binMult = (int) (1.0f/resolution);
        this.pdf = new int[((int) (100.0 * this.binMult + 1))];
        this.pdfDecoy = new int[((int) (100.0 * this.binMult + 1))];
        this.pdfMassError = new int[100]; // Default 100 ppm max, can make it dynamic TODO
        this.pdfMassErrorDecoy = new int[100]; // Default 100 ppm max, can make it dynamic TODO
        if (poissonBinomial) {
            this.mode = "poissonbinomial-matched-theoretical-decoy";
            this.datapoints = new MatchedIonTable();
            this.targetPMF = new TwoDimJointPMF(101, 3001, true, true);
            this.decoyPMF = new TwoDimJointPMF(101, 3001, true, true);
        }
        this.tNegCount = 0;
        this.dNegCount = 0;
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
    /**
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
     **/
    public void addIon(float intensity, float massError, boolean isDecoy) {
        if (!isDecoy) {
            if (intensity > 0.0f) { // Only include matched ions, unmatched ions have negative intensity
                this.datapoints.addRow(massError, intensity, isDecoy);
            } else {
                int intensityIndex = (int) (-1 * intensity);
                this.pdf[intensityIndex]++;
                this.tNegCount++;
            }
        } else {
            if (intensity > 0.0f) { // Only include matched ions, unmatched ions have negative intensity
                this.datapoints.addRow(massError, intensity, isDecoy);
            } else {
                int intensityIndex = (int) (-1 * intensity);
                this.pdfDecoy[intensityIndex]++;
                this.dNegCount++;
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
        // Sort table and find min and max
        this.datapoints.sortRowsByProjVal();
        // Set up arrays to hold ion counts for targets and decoys and q-values
        int arraySize = (int) (10.0 * (Math.abs(this.datapoints.getMaxProjVal() -
                this.datapoints.getMinProjVal()) + 1));
        this.projTargetCounts = new int[arraySize];
        this.projDecoyCounts = new int[arraySize];
        this.qVals = new double[arraySize];

        // Fill arrays with projected values
        for (MatchedIonTable.MatchedIonRow row : this.datapoints.rows) {
            int valIndex = translateLdaValToIndex(row.projVal);
            if (row.isDecoy)
                projDecoyCounts[valIndex]++;
            else
                projTargetCounts[valIndex]++;
        }

        boolean leftToRight = checkLdaProjectionLeftToRight();

        // Fill q-value array with P(ion == true_positive | LDA score)
        int nTargets = 1;
        int nDecoys = 1;
        double cMin = 10000.0;
        double cmax = -10000.0;
        if (leftToRight) {
            // Fill array
            for (int i = 0; i < arraySize; i++) {
                nTargets = this.projTargetCounts[i];
                nDecoys = this.projDecoyCounts[i];
                double q = (double) (nDecoys + 1) / (double) (nDecoys + nTargets + 2);
                this.qVals[i] = q;
                System.out.println(this.qVals[i]);
            }
            // Make array monotonic
            /**
            for (int i = arraySize-1; i >= 0; i--) {
                if (this.qVals[i] < cMin) {
                    cMin = this.qVals[i];
                } else {
                    while ((i >= 0) && (this.qVals[i] >= cMin)) {
                        this.qVals[i] = cMin;
                        i--;
                    }
                    i++;
                }
            }
             **/
        } else {
            for (int i = arraySize-1; i >= 0; i--) {
                nTargets = this.projTargetCounts[i];
                nDecoys = this.projDecoyCounts[i];
                double q = (double) (nDecoys + 1) / (double) (nDecoys + nTargets + 2);
                this.qVals[i] = q;
                System.out.println(this.qVals[i]);
            }
            // Make array monotonic
            /**
            for (int i = 0; i < arraySize; i++) {
                if (this.qVals[i] < cMin) {
                    cMin = this.qVals[i];
                } else {
                    while ((i < arraySize) && this.qVals[i] >= cMin) {
                        this.qVals[i] = cMin;
                        i++;
                    }
                    i--;
                }
            }
             **/
        }
    }

    // Check which direction LDA projected by checking proportion of decoys in first and last half
    public boolean checkLdaProjectionLeftToRight() {
        int arraySize = this.projDecoyCounts.length;
        boolean leftToRight = true;
        int leftToRightDecoyCount = 0;
        int leftToRightTargetCount = 0;
        int rightToLeftDecoyCount = 0;
        int rightToLeftTargetCount = 0;
        for (int i = 0; i < ((arraySize / 2) + 1); i++) {
            leftToRightDecoyCount += this.projDecoyCounts[i];
            leftToRightTargetCount += this.projTargetCounts[i];
            rightToLeftDecoyCount += this.projDecoyCounts[arraySize - (i + 1)];
            rightToLeftTargetCount += this.projTargetCounts[arraySize - (i + 1)];
        }

        if ((rightToLeftDecoyCount / rightToLeftTargetCount) > (leftToRightDecoyCount / leftToRightTargetCount))
            leftToRight = false;

        return leftToRight;
    }


    public void calculateLdaWeights() throws Exception {
        this.ldaProcessor = new LDAProcessor(this.datapoints);
        this.ldaProcessor.solveLDA(executorService);
        this.datapoints.mergeProjectedData(ldaProcessor.projectedData);
    }

    public void calculateNegativePredictiveValue() { // todo see if adding 1 here helps
        this.negPredictiveValue = (double) (this.tNegCount + 1) / (double) (this.dNegCount);

        this.ionPosterior = new double[((int) (100.0 * this.binMult + 1))];
        for (int i = 0; i < this.ionPosterior.length; i++) {
            this.ionPosterior[i] = (double) ((this.pdf[i] + 1) / (double) (this.pdf[i] + this.pdfDecoy[i] + 2));
            System.out.println(this.ionPosterior[i]);
        }
        //this.negPredictiveValue = Math.max(this.negPredictiveValue, 0.01);
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

    private int translateLdaValToIndex(double projVal) {
        int valIndex = (int) ((projVal - this.datapoints.getMinProjVal()) * 10.0);
        if (valIndex < 0)
            return 0;
        else if (valIndex >= this.qVals.length)
            return (this.qVals.length - 1);
        else
            return valIndex;
    }

    public double[] calcIonProbabilities(float [] intensities) {
        double[] probs = new double[intensities.length];
        for (int i = 0; i < intensities.length; i++)
            probs[i] = calcIonProbability(intensities[i]);
        return probs;
    }

    public double calcIonProbability(float intensity, float massError) {
        double prob;
        if (intensity < 0)
            prob = -1;
        else {
            if (!this.mode.equals("poissonbinomial-matched-theoretical-decoy")) {
                prob = this.cdf[(int) intensity * this.binMult];
            } else {
                int intIndx = (int) intensity;
                int massErrorIndx =  (int) massError;
                prob = 1 - this.qValues.getQVal(intIndx, massErrorIndx);
            }
        }
        return prob;
    }

    public double calculateIonProbabilityLda(float intensity, float massError) {
        double projVal;
        double prob;
        if (intensity < 0)
            prob = this.ionPosterior[(int) (intensity * -1 * this.binMult)];
        else {
            projVal = this.ldaProcessor.projectData(massError, intensity).getEntry(0,0);
            int projValIndex = translateLdaValToIndex(projVal);
            prob = 1.0 - this.qVals[projValIndex];
        }
        return prob;
    }

    public double[] calcIonProbabilities(float [] intensities, float[] massErrors) {
        double[] probs = new double[intensities.length];
        for (int i = 0; i < intensities.length; i++) {
            //probs[i] = calcIonProbability(intensities[i], massErrors[i]); todo original function
            probs[i] = calculateIonProbabilityLda(intensities[i], massErrors[i]);
        }
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
            sb.append(this.qValues.toString());
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

    public void printQVals(String fname) throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(fname));
        out.println("lda_val\tq_val");
        for (int i = 0; i < qVals.length; i++)
            out.println((i / 100.0) + "\t" +  qVals[i]);
        out.close();
    }

    public void printFeatureTable(String fname) {
        try {
            PrintWriter out = new PrintWriter(new FileWriter(fname));
            out.println("mass_error\tintensity\tprojected_value\tis_decoy");
            for (MatchedIonTable.MatchedIonRow row : datapoints.rows)
                out.println(row);
            out.close();
        } catch (IOException e) {
          e.printStackTrace();
          System.exit(1);
        }
    }

    public void printPriorProbabilities(String fname) throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(fname));
    }

}
