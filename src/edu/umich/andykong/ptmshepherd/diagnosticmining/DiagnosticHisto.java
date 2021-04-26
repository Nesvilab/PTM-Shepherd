package edu.umich.andykong.ptmshepherd.diagnosticmining;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class DiagnosticHisto {
    Bin [] bins;
    Bin [] smoothBins;
    double min;
    double max;
    int buffersize = 5; //pads endings to mirror internal bins
    double binsPerDa;
    double total;
    double minSignal;
    public ArrayList<Double> filteredPeaks;

    public DiagnosticHisto (double mn, double mx, double binWidth, double minSignal) {
        this.min = (int) mn - buffersize;
        this.max = (int) mx + buffersize;
        this.total = 0;
        this.binsPerDa = (int) (1.0 / binWidth);
        this.minSignal = minSignal;
        this.bins = new Bin [(int)(this.binsPerDa * (this.max - this.min))];
        for (int i = 0; i < this.bins.length; i++)
            this.bins[i] = new Bin();
    }

    private int locateBin(double val) {
        int bini = -1;
        bini = (int) ((val - this.min) * this.binsPerDa);
        return bini;
    }

    public void placeIon(double mz, double intensity) {
        this.bins[locateBin(mz)].bumpVal(intensity);
    }

    public void placeIons(float[][] peaks) {
        HashMap<Integer, ArrayList<Double>> peakToBin = new HashMap<>();
        for (int i = 0; i < peaks.length; i++) {
            int bin = locateBin(peaks[i][0]);
            if (!peakToBin.containsKey(bin))
                peakToBin.put(bin, new ArrayList<>());
            peakToBin.get(bin).add((double) peaks[i][1]);
        }
        for (Integer binIndx : peakToBin.keySet()) {
            double max = -1000000;
            for (Double intensity : peakToBin.get(binIndx)) {
                if (intensity > max)
                    max = intensity;
            }
            this.bins[binIndx].bumpVal(max);
        }
    }

    public double binToMass(int bini) {
        return (((double) bini / (double) this.binsPerDa) + this.min);
    }

    public int findMax() {
        double max = -1;
        int maxi = -1;
        for (int i = 0; i < this.bins.length; i++) {
            if (this.bins[i].val > max) {
                max = this.bins[i].val;
                maxi = i;
            }
        }
        return maxi;
    }

    public void printHisto(String fname) throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(new File(fname)));
        out.printf("mass\theight\tsmooth\n");
        for (int i = 0; i < this.bins.length; i++) {
            StringBuffer sb = new StringBuffer();
            for (int j = 0; j < 10000; j++) {
                if (i + j < this.bins.length) {
                    sb.append(String.format("%.04f\t%.04f\t%.04f\n", binToMass(i), this.bins[i].val, this.smoothBins[i].val));
                    i++;
                } else {
                    continue;
                }
            }
            out.printf(sb.toString());
        }
        out.close();
    }

    public void smoothify() {
        /* Initialize smoothed bins */
        this.smoothBins = new Bin [this.bins.length];
        for (int i = 0; i < this.smoothBins.length; i++)
            this.smoothBins[i] = new Bin();

        /* Calculate bin scaling metrics */
        double ppmTol = 20.0; //todo cast param
        // todo this approach won't work for histos with neg vals...
        // int nBinsBegin = (int) (((ppmTol / 1000000.0) * binToMass(0)) * this.binsPerDa);
        //int nBinsBegin = 3; //not the issue
        //System.out.println(nBinsBegin);
        //int nBinsEnd = (int) ((ppmTol / 1000000.0) * binToMass(this.bins.length - 1) * this.binsPerDa); //todo
        //int nBinsEnd = 3; //not the issue
        //System.out.println(nBinsEnd);
        //int nBinsIncrement = this.bins.length / (nBinsEnd - nBinsBegin); //todo
        //System.out.println(nBinsIncrement);

        /* Initialize rolling average metrics */
        int beginIndx = locateBin(this.min + 1);
        //int beginIndx = 3;
        int endIndx = locateBin(this.max - 1); // Ending a few bins before the end
        int cBins = 1; //todo
        //int cBins = 7;
        int cBinsPerSide = calculateCBinsSide(cBins);

        /* Initialize vals */
        double [] cVals;// = getCVals(beginIndx, cBinsPerSide);
        //int cEldestIndx = 0; // Stores oldest value in cVals
        double cMean;// = calculateCMean(cVals);

        /* Iterate through bins and place new vals in smoothBins */
        for (int i = beginIndx; i < endIndx; i++) {
            /* If this a bin where we need to increase size of cVals */ //todo figure out ppm peakpicking because neg values will break it
            // if (i % nBinsIncrement == 0) {
                /* Reinitialize all variables and recalc mean to reduce floating point error */
            //    cBins++;
            //    cBinsPerSide = calculateCBinsSide(cBins);
            //    cVals = getCVals(i, cBinsPerSide);
            //    cMean = calculateCMean(cVals);
            //    cEldestIndx = 0;
            //}
            /* If this is not a bin where we increase size of cVals */
            //else {
                /* Recalculate mean */
            /* todo
            cVals = getCVals(i, cBinsPerSide);
            cMean = calculateCMean(cVals);
            */
            cMean = this.bins[i].val;

                //cMean -= cVals[cEldestIndx] / cVals.length;
                //cVals[cEldestIndx] = this.bins[i + cBinsPerSide].val;
                //cMean += cVals[cEldestIndx] / cVals.length;
                /* Find where the new oldest value is */
                //if (cEldestIndx == cVals.length - 1)
                //    cEldestIndx = 0;
                //else
                //    cEldestIndx++;
            //}
            this.smoothBins[i].val = cMean;
            this.total += cMean;
            //System.out.printf("%.04f\t%.04f\t%.04f\t%.04f\t%.04f\t%.04f\n", this.bins[i].val, this.smoothBins[i].val, cVals[0], cVals[1], cVals[2], cMean);
        }
    }

    //todo this function currently penalizes wider peaks, aka those from larger mass shifts
    public void findPeaks() {
        this.filteredPeaks = new ArrayList<>();

        double minVal = this.total * this.minSignal; // Peak must be least 0.01 of total / cbins todo??
        System.out.println("Minimum:" + minVal);
        ArrayList<Peak> peaks = new ArrayList<>();
        for (int i = 0; i < this.smoothBins.length; i++) {
            if (this.smoothBins[i].val <= minVal)
                continue;
            int peakIndx = i;
            if (this.smoothBins[i+1].val > this.smoothBins[i].val) { //this can be while instead of if todo
                /* Find local maximum */
                //System.out.println("Entered:"+binToMass(i));
                while (this.smoothBins[i+1].val > this.smoothBins[i].val)
                    i++;
                /* If local max is shared with other peaks, find source peak in unsmoothed histo */
                int nEqualBins = 1;
                while (this.smoothBins[i+nEqualBins-1].val == this.smoothBins[i+nEqualBins].val)
                    nEqualBins++;
                double cUnsmoothedMax = 0;
                for (int j = 0; j < nEqualBins; j++) {
                    if (this.bins[i+j].val > cUnsmoothedMax) {
                        cUnsmoothedMax = this.bins[i+j].val;
                        peakIndx = i+j;
                    }
                }
                /* Store peak and val */ // This will have to move if we decide to measure peak area
                peaks.add(new Peak(binToMass(peakIndx), this.smoothBins[peakIndx].val));
                /* Collect number of downslope bins */
                i += nEqualBins - 1;
                //System.out.println("Equal:"+binToMass(i));
                while (this.smoothBins[i+1].val < this.smoothBins[i].val)
                    i++;
                //System.out.println("Downslope:"+binToMass(i));
                //TODO increment i on the downhill slope as well, otherwise it'll duplicate peaks on the right side
                //TODO am i getting lucky? I think I should be adding i rather than peakIndx...
                // TODO peakindx might be redundant....
            } else {
                //System.out.println("Entered2:" + binToMass(i));
                /* If local max is shared with other peaks, find source peak in unsmoothed histo */
                int nEqualBins = 1;
                while (this.smoothBins[i + nEqualBins - 1].val == this.smoothBins[i + nEqualBins].val)
                    nEqualBins++;
                double cUnsmoothedMax = 0;
                for (int j = 0; j < nEqualBins; j++) {
                    if (this.bins[i + j].val > cUnsmoothedMax) {
                        cUnsmoothedMax = this.bins[i + j].val;
                        peakIndx = i + j;
                    }
                }
                peaks.add(new Peak(binToMass(peakIndx), this.smoothBins[peakIndx].val));
                /* Collect number of downslope bins */
                i += nEqualBins - 1;
                //System.out.println("Equal2:" + binToMass(i));
                while (this.smoothBins[i + 1].val < this.smoothBins[i].val)
                    i++;
                //System.out.println("Downslope2:" + binToMass(i));
            }
        }
        int maxPeaks = Math.min(100, peaks.size());
        Collections.sort(peaks);
        for (int i = 0; i < maxPeaks; i++)
            this.filteredPeaks.add(peaks.get(i).MZ);
        //for (Peak p : peaks) {
        //    System.out.println(p.MZ + "\t" + p.Int);
        //}
    }

    private int calculateCBinsSide(int cBins) {
        int nEachSide;
        if (cBins % 2 == 0)
            nEachSide = (cBins / 2) - 1;
        else
            nEachSide = cBins / 2;
        return nEachSide;
    }

    private double[] getCVals(int cInd, int cBinsSide) {
        double [] cVals = new double[2 * cBinsSide + 1];
        for (int i = cInd - cBinsSide; i <= cInd + cBinsSide; i++)
            cVals[i - cInd + cBinsSide] = this.bins[i].val;
        return cVals;
    }

    private double calculateCMean(double[] cVals) {
        double cMean = 0;
        for (int i = 0; i < cVals.length; i++)
            cMean += cVals[i];
        cMean /= cVals.length;
        return cMean;
    }

    public void clearMemory() {
        this.bins = null;
        this.smoothBins = null;
    }

}

class Bin {
    public double val;
    public int n;
    public Bin() {
        this.val = 0;
    }
    public synchronized void bumpVal(double v) {
        this.val += v;
        this.n++;
    }
}

class Peak implements Comparable<Peak> {
    double MZ, Int;
    public Peak(double MZ, double Int) {
        this.MZ = MZ;
        this.Int = Int;
    }
    public int compareTo(Peak arg0) {
        return -1* Double.compare(this.Int, arg0.Int);
    }
}
