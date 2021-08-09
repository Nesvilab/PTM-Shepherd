package edu.umich.andykong.ptmshepherd.diagnosticmining;

import com.google.common.util.concurrent.AtomicDoubleArray;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.DoubleAdder;
import com.google.common.util.concurrent.AtomicDouble;

public class DiagnosticHisto {
    double peakApex;
    Bin [] bins;
    Bin [] smoothBins;
    AtomicDoubleArray bins2;
    AtomicDoubleArray smoothBins2;
    double min;
    double max;
    int buffersize = 5; //pads endings to mirror internal bins

    int smoothFactor;
    double smoothFactorScale; //multiplying by inverse speeds up math
    double[] gweights;
    double gfrac = 0.95;
    int nBinsPerSide;

    double binsPerDa;
    AtomicDouble total;
    double minSignal;
    double ppmTol;
    double normMass;
    public ArrayList<Double> filteredPeaks;
    public ArrayList<HistoPeak> filteredHistoPeaks;

    public DiagnosticHisto (double peakApex, double mn, double mx, double binWidth, double minSignal, double ppmTol, double normMass) {
        this.peakApex = peakApex;
        this.min = (int) mn - buffersize;
        this.max = (int) mx + buffersize;
        this.total = new AtomicDouble();
        this.binsPerDa = (int) (1.0 / binWidth);
        this.ppmTol = ppmTol;
        this.normMass = normMass;
        this.minSignal = minSignal;
        this.bins2 = new AtomicDoubleArray((int)(this.binsPerDa * (this.max - this.min)));
        this.smoothBins2 = new AtomicDoubleArray((int)(this.binsPerDa * (this.max - this.min)));

        /* Calc optimal smoothing factor */
        this.smoothFactor = calcSmoothFactor(ppmTol, normMass);
        this.smoothFactorScale = 1.0 / smoothFactor;
        this.gweights = calcWeights(this.smoothFactor);
        this.nBinsPerSide = (int) this.smoothFactor / 2;

    }

    private int locateBin(double val) {
        int bini = -1;
        bini = (int) ((val - this.min) * this.binsPerDa);
        return bini;
    }

    public void placeIon(double mz, double intensity) {
        int cBin = locateBin(mz);
        //this.bins[cBin].bumpVal(intensity);
        this.bins2.addAndGet(cBin, intensity);
        this.total.addAndGet(intensity);
        double intPerBin = intensity * this.smoothFactorScale;
        for (int i = -this.nBinsPerSide; i < this.nBinsPerSide; i++)
            //this.smoothBins[cBin + i].bumpVal(intPerBin);
            this.smoothBins2.addAndGet(cBin + i, intPerBin);
            //this.smoothBins[cBin + i].bumpVal(intensity * this.gweights[i + this.nBinsPerSide]);
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

    private int calcSmoothFactor(double ppmTol, double normMass) {
        double daTol = (normMass * ppmTol) / 1000000.0;
        int sf;
        if ((int)(this.binsPerDa * daTol) > 1)
            sf = (int)(this.binsPerDa * daTol);
        else
            sf = 1;
        if (sf % 2 == 0 && sf > 1)
            sf--;
        return sf;
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
        // double ppmTol = 20.0; //todo cast param
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
        int cBins = 7; //todo
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

            cVals = getCVals(i, cBinsPerSide);
            cMean = calculateCMean(cVals);

            //System.out.print(i+"\t" + this.bins[i].val+"\t");
            //for (int j = 0; j < cVals.length; j++)
            //    System.out.print(cVals[j]+" ");
            //System.out.println(cMean);

            //cMean = this.bins[i].val;

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
            //this.total += cMean;
            //System.out.printf("%.04f\t%.04f\t%.04f\t%.04f\t%.04f\t%.04f\n", this.bins[i].val, this.smoothBins[i].val, cVals[0], cVals[1], cVals[2], cMean);
        }
    }

    public void smoothify(ExecutorService executorService, int nThreads) throws ExecutionException, InterruptedException {
        /* Calc optimal smoothing factor */
        this.smoothFactor = calcSmoothFactor(ppmTol, normMass);
        /* Initialize smoothed bins */
        this.smoothBins = new Bin[this.bins.length];
        for (int i = 0; i < this.smoothBins.length; i++)
            this.smoothBins[i] = new Bin();

        /* Initialize rolling average metrics */
        int beginIndx = locateBin(this.min + smoothFactor * 2);
        int endIndx = locateBin(this.max - smoothFactor * 2); // Ending a few bins before the end
        int totalBins = endIndx - beginIndx + 1;

        /* Split histo info blocks */
        final int BLOCKSIZE = totalBins / nThreads; //number of scans to be parsed per thread (to cut down on thread creation overhead)

        /* Set up containers for threads */
        ArrayList<Future> futureList = new ArrayList<>(nThreads);
        for (int i = 0; i < nThreads; i++) {
            int istart = i * BLOCKSIZE + beginIndx;
            int iend = (i + 1) * BLOCKSIZE + beginIndx;
            if (i == nThreads - 1)
                iend = endIndx;
            int finalIend = iend;
            futureList.add(executorService.submit(() -> smoothifyBlock(istart, finalIend)));
        }
        for (Future future : futureList) { //checks to make sure all threads are done
            future.get();
        }

    }

    private void smoothifyBlock(int istart, int iend) {
        int cBins = this.smoothFactor; //todo
        int cBinsPerSide = calculateCBinsSide(cBins);

        /* Initialize vals */
        double [] cVals = getCVals(istart, cBinsPerSide);
        int cEldestIndx = 0; // Stores oldest value in cVals
        double cMean;
        double locTotal = 0;

        /* Iterate through bins and place new vals in smoothBins */
        for (int i = istart; i < iend; i++) {
            //cVals = getCVals(i, cBinsPerSide);
            //cMean = calculateCMean(cVals);
            //double q = cVals[cEldestIndx];
            //double t = this.bins[i+cBinsPerSide].val;
            cVals[cEldestIndx] = this.bins[i+cBinsPerSide].val;
            cMean = calculateCMean(cVals);

            if (cEldestIndx == cBins - 1)
                cEldestIndx = 0;
            else
                cEldestIndx++;

            this.smoothBins[i].val = cMean;
            locTotal += cMean;
        }

        bumpTotal(locTotal);
    }

    private synchronized void bumpTotal(double v) {} //this.total += v; }

    public void findPeaks() {
        this.filteredPeaks = new ArrayList<>();

        double minVal = this.total.doubleValue() * this.minSignal; // Peak must be least 0.01 of total / cbins todo??
        double minEntryVal = 2e-100;
        double prominence = 0.99;

        ArrayList<HistoPeak> peaks = new ArrayList<>();
        for (int i = 0; i < this.smoothBins2.length(); i++) {
            /* Skip bins not worth looking at */
            if (this.smoothBins2.get(i) <= minEntryVal)
                continue;
            /* Once we find a bin worth looking at, set up vals */
            int peakApex = i;
            double maxRaw = this.bins2.get(i);
            double peakArea = this.smoothBins2.get(i);
            /* Integrate uphill slope and find apex */
            while (this.smoothBins2.get(i+1) >= this.smoothBins2.get(i) * prominence) {
                if (this.smoothBins2.get(i+1) > this.smoothBins2.get(i)) {
                    maxRaw = this.bins2.get(i+1);
                    peakApex = i+1;
                }
                else if (this.smoothBins2.get(i+1) == this.smoothBins2.get(i)) {
                    if (this.bins2.get(i+1) > maxRaw) {
                        maxRaw = this.bins2.get(i + 1);
                        peakApex = i + 1;
                    }
                }
                peakArea += this.smoothBins2.get(i);
                i++;
            }
            /* Integrate downhill slope */
            while ((this.smoothBins2.get(i+1) * prominence <= this.smoothBins2.get(i)) &&
                    (this.smoothBins2.get(i) >= minEntryVal)) {
                i++;
                peakArea += this.smoothBins2.get(i);
            }
            /* Check if peak meets criteria and store it */
            if (peakArea >= minVal)
                peaks.add(new HistoPeak(binToMass(peakApex), peakArea));
        }
        /* Select top peaks */
        int maxPeaks = Math.min(1000, peaks.size());
        Collections.sort(peaks);
        for (int i = 0; i < maxPeaks; i++) {
            this.filteredPeaks.add(peaks.get(i).MZ);
            //this.filteredHistoPeaks.add(peaks.get(i));
            //System.out.println(peaks.get(i).MZ + "\t" + (peaks.get(i).Int / this.total));
        }
        ArrayList<Integer> removePeaks = new ArrayList<>();
        for (int i = 0; i < this.filteredPeaks.size(); i++) {
            for (int j = i+1; j < this.filteredPeaks.size(); j++) {
                if (removePeaks.contains(j))
                    continue;
                if (Math.abs(this.filteredPeaks.get(i) - this.filteredPeaks.get(j)) <= 0.0021)
                    removePeaks.add(j);
            }
        }
        Collections.sort(removePeaks, Collections.reverseOrder());
        for (Integer j : removePeaks)
            this.filteredPeaks.remove(j);
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
        this.bins2 = null;
        this.smoothBins = null;
        this.smoothBins2 = null;
    }


    private double [] calcWeights(int Nbins) {
        NormalDistribution nd = new NormalDistribution();
        double a = (1 - this.gfrac) / 2;
        double lower = nd.inverseCumulativeProbability(a);
        double upper = nd.inverseCumulativeProbability(1-a);
        double width = (upper - lower) / Nbins;
        double[] gweights = new double[Nbins];
        for(int i = 0; i < Nbins; i++) {
            double lo = nd.cumulativeProbability(lower+i*width);
            double hi = nd.cumulativeProbability(lower +(i+1)*width);
            gweights[i] = hi-lo;
        }

        return gweights;
    }

}

class Bin {
    public double val;
    public int n;
    public Bin() {
        this.val = 0;
    }
    public synchronized void setVal(double v) { this.val = v; }
    //public synchronized void bumpVal(double v) {
    public synchronized void bumpVal(double v) {
        this.val += v;
        this.n++;
    }
}

class HistoPeak implements Comparable<HistoPeak> {
    double MZ, Int;
    public HistoPeak(double MZ, double Int) {
        this.MZ = MZ;
        this.Int = Int;
    }
    public int compareTo(HistoPeak arg0) {
        return -1* Double.compare(this.Int, arg0.Int);
    }
}
