package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.core.Spectrum;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class BinDiagMetric {
    double peakApex;
    double leftBound;
    double rightBound;
    //public ArrayList<Spectrum> spectra;
    //public ArrayList<String> specNames;
    public int nSpecs;
    public DiagnosticHisto immoniumIons;
    public DiagnosticHisto capYIons;
    public ArrayList<DiagnosticHisto> tildeIons;
    public String ionTypes;
    public double binMinMax[][]; /* [n][2], iontype -> {min, max} */
    public double ppmTol;
    public HashMap<String, ArrayList<DiagnosticRecord>> peptideMap;
    public PeakCompareTester testResults;

    public BinDiagMetric(double[] peakBounds, String ions, double ppmTol) {
        this.peakApex = peakBounds[0];
        this.leftBound = peakBounds[1];
        this.rightBound = peakBounds[2];
        this.ionTypes = ions;
        this.binMinMax = new double[2+ions.length()][2];
        this.peptideMap = new HashMap<>();
        this.ppmTol = ppmTol;
    }

    /* Build map of peptides so that each peptide (not PSM) can be weighted equally */
    public void addPSMToPeptideMap(DiagnosticRecord dr) {
        String pepKey = dr.pepSeq + dr.modifications + dr.charge; //peptidekey includes charge ///todo is that ideal?
        /* Check if pepKey exists and add line to pepKey */
        if (!this.peptideMap.containsKey(pepKey))
            this.peptideMap.put(pepKey, new ArrayList<>());
        this.peptideMap.get(pepKey).add(dr);
    }

    /* Sends the peptides to the histogram */
    public void processPeptideMap(ExecutorService executorService, int nThreads, double minSignal, int maxCharge) throws Exception {

        /* Prepopulate histo min/max for each ion type */
        for (int i = 0; i < this.binMinMax.length; i++) {
            this.binMinMax[i][0] = 1000000.0; // Big minimum
            this.binMinMax[i][1] = -1000000.0; // Small maximum
        }

        double avgImm = 150;
        double avgPepPrec = 0;
        double avgFrag[] = new double[this.ionTypes.length()];
        double avgPepLen = 0;

        /* Find the mins and maxes of each histo */
        int nPepKeys = this.peptideMap.keySet().size();
        for (String pepKey : this.peptideMap.keySet()) {
            int nDrs = this.peptideMap.get(pepKey).size();
            int pepLen = 0;
            for (DiagnosticRecord dr : this.peptideMap.get(pepKey)) {
                pepLen = dr.pepSeq.length();
                int binMinMaxIndx = 0;
                for (int i = 0; i < dr.immoniumPeaks.length; i++) {
                    float mz = dr.immoniumPeaks[i][0]; // Cycle through mz vals
                    //avgImm += (mz / nDrs) / nPepKeys;
                    if (mz < this.binMinMax[binMinMaxIndx][0])
                        this.binMinMax[binMinMaxIndx][0] = mz;
                    if (mz > this.binMinMax[binMinMaxIndx][1])
                        this.binMinMax[binMinMaxIndx][1] = mz;
                }
                binMinMaxIndx++;
                for (int i = 0; i < dr.capYPeaks.length; i++) {
                    float mz = dr.capYPeaks[i][0]; // Cycle through mz vals
                    if (mz < this.binMinMax[binMinMaxIndx][0])
                        this.binMinMax[binMinMaxIndx][0] = mz;
                    if (mz > this.binMinMax[binMinMaxIndx][1])
                        this.binMinMax[binMinMaxIndx][1] = mz;
                }
                avgPepPrec += (dr.calcAvgCapYTol(1, maxCharge) / nDrs) / nPepKeys;
                for (int h = 0; h < dr.ionTypes.length(); h++) {
                    binMinMaxIndx++;
                    char ionType = dr.ionTypes.charAt(h);
                    for (int i = 0; i < dr.squigglePeaks.get(ionType).length; i++) {
                        float mz = dr.squigglePeaks.get(ionType)[i][0]; // Cycle through mz vals
                        if (mz < this.binMinMax[binMinMaxIndx][0])
                            this.binMinMax[binMinMaxIndx][0] = mz;
                        if (mz > this.binMinMax[binMinMaxIndx][1])
                            this.binMinMax[binMinMaxIndx][1] = mz;
                    }
                    avgFrag[h] += (dr.calcAvgFragTol(ionType, 1) / nDrs) / nPepKeys; //todo charge states
                }
            }
            avgPepLen += (double) pepLen / nPepKeys;
        }

        /* Initialize histograms */
        this.immoniumIons = new DiagnosticHisto(this.peakApex, this.binMinMax[0][0], this.binMinMax[0][1], 0.001, minSignal, this.ppmTol, avgImm);
        this.capYIons = new DiagnosticHisto(this.peakApex, this.binMinMax[1][0], this.binMinMax[1][1], 0.001, minSignal, this.ppmTol, avgPepPrec);
        this.tildeIons = new ArrayList<>();
        for (int i = 0; i < this.ionTypes.length(); i++) { /* +2 because of immonium and Y ions in first 2 i's */
            tildeIons.add(new DiagnosticHisto(this.peakApex, this.binMinMax[i + 2][0], this.binMinMax[i + 2][1], 0.001, minSignal / avgPepLen, this.ppmTol, avgFrag[i]));
        }

        /* Assign data to histograms */
        List<Future> futureList = new ArrayList<>(this.peptideMap.size());
        for (String pepKey : this.peptideMap.keySet()) {
            futureList.add(executorService.submit(() -> {
                double nPsms = this.peptideMap.get(pepKey).size();
                for (DiagnosticRecord dr : this.peptideMap.get(pepKey)) {
                    this.nSpecs++;
                    for (int i = 0; i < dr.immoniumPeaks.length; i++)
                        this.immoniumIons.placeIon(dr.immoniumPeaks[i][0], (double) dr.immoniumPeaks[i][1] / nPsms);
                    for (int i = 0; i < dr.capYPeaks.length; i++)
                        this.capYIons.placeIon(dr.capYPeaks[i][0], (double) dr.capYPeaks[i][1] / nPsms);
                    for (int h = 0; h < this.ionTypes.length(); h++) {
                        char ionType = this.ionTypes.charAt(h);
                        for (int i = 0; i < dr.squigglePeaks.get(ionType).length; i++) {
                            double mz = dr.squigglePeaks.get(ionType)[i][0];
                            if (!(-3.5 < mz && mz < 3.5)) // Filter out small mzs here //todo should be a histo param
                                this.tildeIons.get(h).placeIon(mz, ((dr.squigglePeaks.get(ionType)[i][1] / (double) dr.pepSeq.length()) / nPsms));
                        }
                    }
                }
            }));
        }

        /* Wait for histograms to be populated */
        for (Future future : futureList)
            future.get();

        /* Free pepMap memory */
        this.peptideMap = null;

        System.out.println("\nPeakApex:"+peakApex);

        this.immoniumIons.smoothify(executorService, nThreads);
        this.capYIons.smoothify(executorService, nThreads);
        //System.out.println("ImmoniumAbund");
        this.immoniumIons.findPeaks();
        //System.out.println("CapYAbund");
        this.capYIons.findPeaks();


        this.immoniumIons.clearMemory();
        this.capYIons.clearMemory();

        for (int i = 0; i < this.tildeIons.size(); i++) {
            long t1 = System.currentTimeMillis();
            //System.out.println("IonType:"+ionTypes.charAt(i));
            this.tildeIons.get(i).smoothify(executorService, nThreads);
            long t2 = System.currentTimeMillis();
            //System.out.println("SquiggleAbund");
            if (Math.abs(this.peakApex - 226.08) < 1)
                this.tildeIons.get(i).printHisto(this.peakApex + "_" + this.ionTypes.charAt(i) + ".tsv");
            this.tildeIons.get(i).findPeaks();
            long t3 = System.currentTimeMillis();
            this.tildeIons.get(i).clearMemory();
            //System.out.printf("Processing time (%d ms smoothing - %d ms peakpicking)\n", t2-t1, t3-t2);
        }
    }

    public void setTestResults(PeakCompareTester pct) {
        this.testResults = pct;
    }

    public String toString() {
        StringBuffer newLines = new StringBuffer();

        /* Format immonium tests */
        for (int i = 0; i < this.testResults.immoniumTests.size(); i++) {
            Test t = this.testResults.immoniumTests.get(i);
            String newLine = String.format("%.04f\tdiagnostic\t%.04f\t%.04f\t%e\t%f\t%.04f\t%d\t%d\n",
                    this.peakApex, t.mass, t.adjustedMass, t.q, t.rbc, t.u, t.n1, t.n2);
            newLines.append(newLine);
        }

        /* Format capY tests */
        for (int i = 0; i < this.testResults.capYTests.size(); i++) {
            Test t = this.testResults.capYTests.get(i);
            String newLine = String.format("%.04f\tY\t%.04f\t%.04f\t%e\t%f\t%.04f\t%d\t%d\n",
                    this.peakApex, t.mass, t.adjustedMass, t.q, t.rbc, t.u, t.n1, t.n2);
            newLines.append(newLine);
        }

        /* Format squiggle ions */
        for (Character cIon : this.testResults.squigglesTests.keySet()) {
            for (int i = 0; i < this.testResults.squigglesTests.get(cIon).size(); i++) {
                Test t = this.testResults.squigglesTests.get(cIon).get(i);
                String newLine = String.format("%.04f\t%c\t%.04f\t%.04f\t%e\t%f\t%.04f\t%d\t%d\n",
                        this.peakApex, cIon, t.mass, t.adjustedMass, t.q, t.rbc, t.u, t.n1, t.n2);
                newLines.append(newLine);
            }
        }

        return newLines.toString();
    }
}

class SquigglePeak implements Comparable<SquigglePeak> {
    double MZ, Int;
    public SquigglePeak(double MZ, double Int) {
        this.MZ = MZ;
        this.Int = Int;
    }
    public int compareTo(SquigglePeak arg0) {
        return -1* Double.compare(this.MZ, arg0.MZ);
    }
}
