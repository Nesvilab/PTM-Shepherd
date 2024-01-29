package edu.umich.andykong.ptmshepherd.iterativelocalization;

import java.util.ArrayList;
import java.util.Arrays;

class BinPriorProbabilities {
    //Holds the values to compute prior probability Pr(Pep_{ij})
    //Note: the original manuscript has 104 values, we use 62 because we do not consider protein n and c termini
    //Note: the length of these is 26 to make character indexing possible, but only the 20 AAs will be valid
    double [] probs;
    double [] nProbs;
    double [] cProbs;
    double n;
    double c;

    // Previous epoch's values
    double [] probsT1;
    double [] nProbsT1;
    double [] cProbsT1;
    double nT1;
    double cT1;
    int [] probsT1Count;
    int [] nProbsT1Count;
    int [] cProbsT1Count;
    int nT1Count;
    int cT1Count;
    int nPsms;
    double priorVectorNorm;
    boolean isConverged;
    int epoch;
    final double INITVAL = 1.0 / 62.0;
    final char[] AAs = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};

    BinPriorProbabilities() {
        this.probs = new double[26];
        this.nProbs = new double[26];
        this.cProbs = new double[26];
        Arrays.fill(this.probs, this.INITVAL);
        Arrays.fill(this.nProbs, this.INITVAL);
        Arrays.fill(this.cProbs, this.INITVAL);
        this.n = this.INITVAL;
        this.c = this.INITVAL;

        this.probsT1Count = new int[26];
        this.nProbsT1Count = new int[26];
        this.cProbsT1Count = new int[26];
        this.nT1Count = 0;
        this.cT1Count = 0;
        this.nPsms = 0;

        this.probsT1 = new double[26];
        this.nProbsT1 = new double[26];
        this.cProbsT1 = new double[26];
        this.nT1 = 0.0;
        this.cT1 = 0.0;

        this.epoch = 0;
        this.priorVectorNorm = 0.0;
        this.isConverged = false;

    }

    double[] computePriorProbs(String pep, boolean[] allowedPoses) {
        double[] priorProbs = new double[pep.length()+2];
        double sumProbs = 0.0;

        // Build peptide prior probability distribution
        // The prior probability function includes two extra values for N- and C-term modifications
        for (int i = 0; i < pep.length(); i++) { // Iterate over middle n-2 indices
            char cAA = pep.charAt(i);
            if (allowedPoses[i+1] == true) {
                if (i == 0) { // Pick max probability for N-terminal options
                    double curProb = Math.max(this.nProbs[cAA-'A'], this.probs[cAA-'A']);
                    priorProbs[i+1] = curProb;
                    sumProbs += curProb;
                } else if (i == pep.length()-1) { // Pick max probability for C-terminal options
                    double curProb = Math.max(this.cProbs[cAA-'A'], this.probs[cAA-'A']);
                    priorProbs[i+1] = curProb;
                    sumProbs += curProb;
                } else {
                    double curProb = this.probs[cAA-'A'];
                    priorProbs[i+1] = curProb;
                    sumProbs += curProb;
                }
            }
        }
        // Add termini //TODO this will not handle the case where termini are allowed but n* is not
        if (allowedPoses[0]) {
            priorProbs[0] = this.n;
            sumProbs += this.n;
        }
        if (allowedPoses[allowedPoses.length-1]) {
            priorProbs[priorProbs.length - 1] = this.c;
            sumProbs += this.c;
        }

        // Normalize prior probability distribution
        for (int i = 0; i < priorProbs.length; i++)
            priorProbs[i] /= sumProbs;

        return priorProbs;
    }

    static double[] computeUniformPriorProbs(String pep, boolean[] allowedPoses) {
        double[] priorProbs = new double[pep.length()+2];
        double sumProbs = 0.0;

        // Build peptide prior probability distribution
        for (int i = 0; i < pep.length()+2; i++) {
            if (allowedPoses[i] == true) {
                priorProbs[i] = 1.0;
                sumProbs += 1.0;
            }
        }

        // Normalize prior probability distribution
        for (int i = 0; i < priorProbs.length; i++)
            priorProbs[i] /= sumProbs;

        return priorProbs;
    }

    public void update(String pep, double[] siteProbs, boolean [] allowedPoses) {
        // Update termini probabilities
        if (allowedPoses[0]) {
            this.nT1 += siteProbs[0];
            this.nT1Count++;
        }
        if (allowedPoses[allowedPoses.length-1]) {
            this.cT1 += siteProbs[siteProbs.length - 1];
            this.cT1Count++;
        }

        // Update AA probabilities
        for (int i = 0; i < pep.length(); i++) {
            char cAA = pep.charAt(i);

            if (i == 0 && allowedPoses[1]) { // Update N-terminal AA probabilities
                this.nProbsT1[cAA-'A'] += siteProbs[1];
                this.nProbsT1Count[cAA-'A']++;

            } else if (i == pep.length()-1 && allowedPoses[pep.length()]) { // Update C-termini AA probabilities
                this.cProbsT1[cAA-'A'] += siteProbs[pep.length()];
                this.cProbsT1Count[cAA-'A']++;
            }

            // Update base AA probabilities
            if (allowedPoses[i+1]) {
                this.probsT1[cAA - 'A'] += siteProbs[i + 1];
                this.probsT1Count[cAA - 'A']++;
            }
        }
        this.nPsms++;
    }

    boolean calcConvergence(double convCriterion) {
        if (this.isConverged == true)
            return true;
        // Normalize T1 probabilities values
        //System.out.println(nT1);
        //System.out.println(cT1);
        //System.out.println(Arrays.toString(nProbsT1));
        //System.out.println(Arrays.toString(cProbsT1));
        //System.out.println(Arrays.toString(probsT1));
        //System.out.println(Arrays.toString(probsT1Count));

        this.nT1 = safeDivide(this.nT1, this.nT1Count);
        this.cT1 = safeDivide(this.cT1, this.cT1Count);
        for (int i = 0; i < this.nProbsT1.length; i++)
            this.nProbsT1[i] = safeDivide(this.nProbsT1[i], this.nProbsT1Count[i]);
        for (int i = 0; i < this.cProbsT1.length; i++)
            this.cProbsT1[i] = safeDivide(this.cProbsT1[i], this.cProbsT1Count[i]);
        for (int i = 0; i < this.probsT1.length; i++)
            this.probsT1[i] = safeDivide(this.probsT1[i], this.probsT1Count[i]);

        // Calculate L1 norm of t1 - t0
        double norm = 0.0;
        for (Character c : this.AAs) {
            norm += Math.abs(this.nProbsT1[c-'A'] - this.nProbs[c-'A']);
            norm += Math.abs(this.cProbsT1[c-'A'] - this.cProbs[c-'A']);
            norm += Math.abs(this.probsT1[c-'A'] - this.probs[c-'A']);
        }
        norm += Math.abs(this.nT1 - this.n);
        norm += Math.abs(this.cT1 - this.c);

        // Reset values
        this.n = this.nT1;
        this.c = this.cT1;
        this.nProbs = Arrays.copyOf(this.nProbsT1, this.nProbs.length);
        this.cProbs = Arrays.copyOf(this.cProbsT1, this.cProbs.length);
        this.probs = Arrays.copyOf(this.probsT1, this.probs.length);

        this.nT1 = 0.0;
        this.cT1 = 0.0;
        Arrays.fill(this.nProbsT1, 0.0);
        Arrays.fill(this.cProbsT1, 0.0);
        Arrays.fill(this.probsT1, 0.0);
        this.nT1Count = 0;
        this.cT1Count = 0;
        Arrays.fill(this.nProbsT1Count, 0);
        Arrays.fill(this.cProbsT1Count, 0);
        Arrays.fill(this.probsT1Count, 0);

        this.priorVectorNorm = norm;

        // Check convergence
        if (norm <= convCriterion)
            this.isConverged = true;

        return this.isConverged;
    }

    private double safeDivide(double x, int y) {
        if (y == 0)
            return 0;
        else
            return x / (double) y;
    }

    public double getPriorVectorNorm() {
        return this.priorVectorNorm;
    }

    public boolean getIsConverged() {
        return this.isConverged;
    }


    public String toString() {
        StringBuffer sb = new StringBuffer("\nprobs\t");
        for (Character c : this.AAs) {
            sb.append(c + ":" + Double.toString(this.probs[c-'A'])+"  ");
        }
        sb.append("\nnprobs\t");
        for (Character c : this.AAs) {
            sb.append(c + ":" + Double.toString(this.nProbs[c-'A'])+"  ");
        }
        sb.append("\ncprobs\t");
        for (Character c : this.AAs) {
            sb.append(c + ":" + Double.toString(this.cProbs[c-'A'])+"  ");
        }
        sb.append("\nn\t");
        sb.append(this.n);
        sb.append("\nc\t");
        sb.append(this.c);
        sb.append("\n");

        return sb.toString();
    }


}

