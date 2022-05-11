/*
 *    Copyright 2022 University of Michigan
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

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
    public String modName;

    ArrayList<DiagnosticProfileRecord> diagProfRecs;

    public BinDiagMetric(double[] peakBounds, String ions, double ppmTol, String modName) {
        this.peakApex = peakBounds[0];
        this.leftBound = peakBounds[1];
        this.rightBound = peakBounds[2];
        this.ionTypes = ions;
        this.binMinMax = new double[2+ions.length()][2];
        this.peptideMap = new HashMap<>();
        this.ppmTol = ppmTol;
        this.modName = modName;
    }

    /* Build map of peptides so that each peptide (not PSM) can be weighted equally */
    public void addPSMToPeptideMap(DiagnosticRecord dr) {
        if (dr.isMangled == true)
            return;
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
            this.binMinMax[i][0] = 0; // Big minimum
            this.binMinMax[i][1] = 0; // Small maximum
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
                pepLen = dr.pepSeq.length() - 1;
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
                    avgFrag[h] += (dr.calcAvgFragTol(ionType, 1) / nDrs) / nPepKeys;
                }
            }
            avgPepLen += (double) pepLen / nPepKeys;
        }

        /* Initialize histograms */
        this.immoniumIons = new DiagnosticHisto(this.peakApex, this.binMinMax[0][0], this.binMinMax[0][1], 0.0002, minSignal, this.ppmTol, avgImm, nPepKeys);
        this.capYIons = new DiagnosticHisto(this.peakApex, this.binMinMax[1][0], this.binMinMax[1][1], 0.0002, minSignal, this.ppmTol, avgPepPrec, nPepKeys);
        this.tildeIons = new ArrayList<>();
        for (int i = 0; i < this.ionTypes.length(); i++) { /* +2 because of immonium and Y ions in first 2 i's */
            tildeIons.add(new DiagnosticHisto(this.peakApex, this.binMinMax[i + 2][0], this.binMinMax[i + 2][1], 0.0002, minSignal, this.ppmTol, avgFrag[i], nPepKeys));
        }

        /* Assign data to histograms */
        List<Future> futureList = new ArrayList<>(this.peptideMap.size());
        for (String pepKey : this.peptideMap.keySet()) {
            futureList.add(executorService.submit(() -> {
                int nPsms = this.peptideMap.get(pepKey).size();
                for (DiagnosticRecord dr : this.peptideMap.get(pepKey)) {
                    this.nSpecs++;
                    for (int i = 0; i < dr.immoniumPeaks.length; i++)
                        this.immoniumIons.placeIon(dr.immoniumPeaks[i][0], (double) dr.immoniumPeaks[i][1] / nPsms, nPsms);
                    for (int i = 0; i < dr.capYPeaks.length; i++)
                        this.capYIons.placeIon(dr.capYPeaks[i][0], (double) dr.capYPeaks[i][1] / nPsms, nPsms);
                    for (int h = 0; h < this.ionTypes.length(); h++) {
                        char ionType = this.ionTypes.charAt(h);
                        for (int i = 0; i < dr.squigglePeaks.get(ionType).length; i++) {
                            double mz = dr.squigglePeaks.get(ionType)[i][0];
                            this.tildeIons.get(h).placeIon(mz, (double) (dr.squigglePeaks.get(ionType)[i][1] / nPsms), nPsms);
                        }
                    }
                }
                //System.out.println(nSpecs + "\t" + nPsms);
            }));
        }

        /* Wait for histograms to be populated */
        for (Future future : futureList)
            future.get();

        /* Free pepMap memory */
        this.peptideMap = null;

        //System.out.println("\nPeakApex:"+peakApex);

        //this.immoniumIons.smoothify(executorService, nThreads);
        //this.capYIons.smoothify(executorService, nThreads);
        //System.out.println("ImmoniumAbund");
        this.immoniumIons.findPeaks();
        //this.immoniumIons.printHisto(peakApex+"_diag.tsv");
        //System.out.println("CapYAbund");
        this.capYIons.findPeaks();
        //this.capYIons.printHisto(peakApex+"_peptide.tsv");


        this.immoniumIons.clearMemory();
        this.capYIons.clearMemory();

        for (int i = 0; i < this.tildeIons.size(); i++) {
            long t1 = System.currentTimeMillis();
            //System.out.println("IonType:"+ionTypes.charAt(i));
            //this.tildeIons.get(i).smoothify(executorService, nThreads);
            long t2 = System.currentTimeMillis();
            //System.out.println("SquiggleAbund");
            this.tildeIons.get(i).findPeaks();
            //this.tildeIons.get(i).printHisto(this.peakApex + "_" + this.ionTypes.charAt(i) + ".tsv");
            long t3 = System.currentTimeMillis();
            this.tildeIons.get(i).clearMemory();
            //System.out.printf("Processing time (%d ms smoothing - %d ms peakpicking)\n", t2-t1, t3-t2);
        }
    }

    public void setTestResults(PeakCompareTester pct) {
        this.testResults = pct;
    }

    public String toString(boolean printNonReps) {
        StringBuffer newLines = new StringBuffer();

        /* Format immonium tests */
        for (int i = 0; i < this.testResults.immoniumTests.size(); i++) {
            Test t = this.testResults.immoniumTests.get(i);
            if (t.isDecoy == true)
                continue;
            if (!printNonReps && !t.isIsotopeRep)
                continue;
            String newLine = String.format("%.04f\tdiagnostic\t%.04f\t%.04f\t%e\t%f\t%.02f\t%.02f\t%.02f\t%.02f\t%d\t%d\n",
                    this.peakApex, t.mass, t.adjustedMass, t.q, t.rbc,
                    t.propWIonTreat, Double.isNaN(t.propWIonCont) ? 0 : t.propWIonCont,
                    t.propWIonIntensity, t.propWIonIntensityCont,
                    t.n1, t.n2);
            newLines.append(newLine);
        }

        /* Format capY tests */
        for (int i = 0; i < this.testResults.capYTests.size(); i++) {
            Test t = this.testResults.capYTests.get(i);
            if (t.isDecoy == true)
                continue;
            if (!printNonReps && !t.isIsotopeRep)
                continue;
            String newLine = String.format("%.04f\tY\t%.04f\t%.04f\t%e\t%f\t%.02f\t%.02f\t%.02f\t%.02f\t%d\t%d\n",
                    this.peakApex, t.mass, t.adjustedMass, t.q, t.rbc,
                    t.propWIonTreat, Double.isNaN(t.propWIonCont) ? 0 : t.propWIonCont,
                    t.propWIonIntensity, t.propWIonIntensityCont,
                    t.n1, t.n2);
            newLines.append(newLine);
        }

        /* Format squiggle ions */
        for (Character cIon : this.testResults.squigglesTests.keySet()) {
            for (int i = 0; i < this.testResults.squigglesTests.get(cIon).size(); i++) {
                Test t = this.testResults.squigglesTests.get(cIon).get(i);
                if (t.isDecoy == true)
                    continue;
                if (!printNonReps && !t.isIsotopeRep)
                    continue;
                if (!printNonReps && !t.isValid)
                    continue;
                String newLine = String.format("%.04f\t%c\t%.04f\t%.04f\t%e\t%f\t%.02f\t%.02f\t%.02f\t%.02f\t%d\t%d\n",
                        this.peakApex, cIon, t.mass, t.adjustedMass, t.q, t.rbc,
                        t.propWIonTreat, Double.isNaN(t.propWIonCont) ? 0 : t.propWIonCont,
                        t.propWIonIntensity, t.propWIonIntensityCont,
                        t.n1, t.n2);
                newLines.append(newLine);
            }
        }

        return newLines.toString();
    }

    public String toString() {
        StringBuffer newLines = new StringBuffer();
        /* Format tests */
        for (int i = 0; i < this.diagProfRecs.size(); i++) {
            String newLine = this.diagProfRecs.get(i).toString();
            if (!newLine.equals(""))
                newLines.append(newLine);
        }
        return newLines.toString();
    }

    public void initializeDiagProfRecs() {
        this.diagProfRecs =  new ArrayList<>();
        for (int i = 0; i < this.testResults.immoniumTests.size(); i++) {
            Test t = this.testResults.immoniumTests.get(i);
            if (!t.isIsotopeRep)
                continue;
            this.diagProfRecs.add(new DiagnosticProfileRecord(this.peakApex, this.modName, "diagnostic", t.mass, t.adjustedMass, t.q, t.rbc, t.propWIonTreat, t.propWIonCont, t.propWIonIntensity, t.propWIonIntensityCont));
        }
        for (int i = 0; i < this.testResults.capYTests.size(); i++) {
            Test t = this.testResults.capYTests.get(i);
            if (!t.isIsotopeRep)
                continue;
            this.diagProfRecs.add(new DiagnosticProfileRecord(this.peakApex, this.modName,"peptide", t.mass, t.adjustedMass, t.q, t.rbc, t.propWIonTreat, t.propWIonCont, t.propWIonIntensity, t.propWIonIntensityCont));
        }
        for (Character cIon : this.testResults.squigglesTests.keySet()) {
            for (int i = 0; i < this.testResults.squigglesTests.get(cIon).size(); i++) {
                Test t = this.testResults.squigglesTests.get(cIon).get(i);
                if (!t.isIsotopeRep)
                    continue;
                if (!t.isValid)
                    continue;
                this.diagProfRecs.add(new DiagnosticProfileRecord(this.peakApex, this.modName, Character.toString(cIon), t.mass, t.adjustedMass, t.q, t.rbc, t.propWIonTreat, t.propWIonCont, t.propWIonIntensity, t.propWIonIntensityCont));
            }
        }
    }
}
