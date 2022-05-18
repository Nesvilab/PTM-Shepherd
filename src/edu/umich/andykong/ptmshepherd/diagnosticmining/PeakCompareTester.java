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

import com.google.common.util.concurrent.AtomicDouble;
import com.google.gson.internal.$Gson$Types;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
//import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import umich.ms.datatypes.lcmsrun.Hash;
import umich.ms.fileio.filetypes.agilent.cef.jaxb.P;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.math3.stat.inference.TTest;

public class PeakCompareTester {
    double peakApex;
    HashMap<Double, ArrayList<Double>> immoniumX;
    HashMap<Double, ArrayList<Double>> immoniumY;
    HashMap<Double, ArrayList<Double>> immoniumDecoy;
    HashMap<Double, ArrayList<Double>> capYX;
    HashMap<Double, ArrayList<Double>> capYY;
    HashMap<Double, ArrayList<Double>> capYDecoy;
    HashMap<Character, HashMap<Double,ArrayList<Double>>> squigglesX;
    HashMap<Character, HashMap<Double,ArrayList<Double>>> squigglesY;
    HashMap<Character, HashMap<Double,ArrayList<Double>>> squigglesDecoy;

    HashMap<String, ArrayList<DiagnosticRecord>> contPeptideMap;
    int nControlPsms;
    HashMap<String, ArrayList<DiagnosticRecord>> treatPeptideMap;
    int nTreatPsms;
    HashMap<String, ArrayList<DiagnosticRecord>> decoyPeptideMap;

    double [] x;
    double [] y;
    ArrayList<Double> xNums;
    ArrayList<Double> yNums;

    public ArrayList<Test> immoniumTests;
    public ArrayList<Test> capYTests;
    public HashMap<Character, ArrayList<Test>> squigglesTests;
    public ArrayList<Test> immoniumDecoyTests;
    public ArrayList<Test> capYDecoyTests;
    public HashMap<Character, ArrayList<Test>> squigglesDecoyTests;

    double maxP;
    double minRbc;
    double minSpecDiff;
    double minFoldChange;
    boolean twoTailedTests;
    int minPeps;

    double diagMinFoldChange;
    double diagMinSpecDiff;
    double pepMinFoldChange;
    double pepMinSpecDiff;
    double fragMinFoldChange;
    double fragMinSpecDiff;
    double fragMinPropensity;

    double specTol;

    boolean debug;

    public PeakCompareTester(double peakApex, ArrayList<Double> unifImm, ArrayList<Double> unifCapY, HashMap<Character, ArrayList<Double>> unifSquig, double maxP, double minRbc, double minSpecDiff, double minFoldChange, boolean twoTailedTests, double specTol) {
        this.peakApex = peakApex;
        this.immoniumX = new HashMap<>();
        this.immoniumY = new HashMap<>();
        this.immoniumDecoy = new HashMap<>();
        for (Double peak : unifImm) {
            this.immoniumX.put(peak, new ArrayList<>());
            this.immoniumY.put(peak, new ArrayList<>());
            this.immoniumDecoy.put(peak, new ArrayList<>());
        }
        this.capYX = new HashMap<>();
        this.capYY = new HashMap<>();
        this.capYDecoy = new HashMap<>();
        for (Double peak : unifCapY) {
            this.capYX.put(peak, new ArrayList<>());
            this.capYY.put(peak, new ArrayList<>());
            this.capYDecoy.put(peak, new ArrayList<>());
        }
        this.squigglesX = new HashMap<>();
        this.squigglesY = new HashMap<>();
        this.squigglesDecoy = new HashMap<>();
        for (Character ion : unifSquig.keySet()) {
            this.squigglesX.put(ion, new HashMap<>());
            this.squigglesY.put(ion, new HashMap<>());
            this.squigglesDecoy.put(ion, new HashMap<>());
            for (Double peak : unifSquig.get(ion)) {
                this.squigglesX.get(ion).put(peak, new ArrayList<>());
                this.squigglesY.get(ion).put(peak, new ArrayList<>());
                this.squigglesDecoy.get(ion).put(peak, new ArrayList<>());
            }
        }
        this.contPeptideMap = new HashMap<>();
        this.treatPeptideMap = new HashMap<>();
        this.decoyPeptideMap = new HashMap<>();
        this.nControlPsms = 0; //actually pepkey count
        this.nTreatPsms = 0; //actually pepkey count

        this.maxP = maxP;
        this.minRbc = minRbc;
        this.minSpecDiff = minSpecDiff;
        this.twoTailedTests = twoTailedTests;
        this.minPeps = Integer.parseInt(PTMShepherd.getParam("diagmine_minPeps"));
        this.specTol = specTol;

        this.diagMinFoldChange = Double.parseDouble(PTMShepherd.getParam("diagmine_diagMinFoldChange"));
        this.diagMinSpecDiff = Double.parseDouble(PTMShepherd.getParam("diagmine_diagMinSpecDiff")) / 100.0;
        this.pepMinFoldChange = Double.parseDouble(PTMShepherd.getParam("diagmine_pepMinFoldChange"));
        this.pepMinSpecDiff = Double.parseDouble(PTMShepherd.getParam("diagmine_pepMinSpecDiff")) / 100.0;
        this.fragMinFoldChange = Double.parseDouble(PTMShepherd.getParam("diagmine_fragMinFoldChange"));
        this.fragMinSpecDiff = Double.parseDouble(PTMShepherd.getParam("diagmine_fragMinSpecDiff")) / 100.0;
        this.fragMinPropensity = Double.parseDouble(PTMShepherd.getParam("diagmine_fragMinFoldChange"));
        this.minFoldChange = minFoldChange;

        this.debug = Boolean.parseBoolean(PTMShepherd.getParam("diagmine_printDebug"));
    }

    public synchronized void addDrs(ArrayList<DiagnosticRecord> drs, boolean isControl, boolean isDecoy) {
        for (DiagnosticRecord dr : drs)
            addVals(dr, isControl, isDecoy);
    }

    public void addVals(DiagnosticRecord dr, boolean isControl, boolean isDecoy) {
        String pepKey = dr.pepSeq + dr.modifications + dr.charge; //peptidekey includes charge
        if (isControl) {
            this.nControlPsms++;
            if (!this.contPeptideMap.containsKey(pepKey))
                this.contPeptideMap.put(pepKey, new ArrayList<>());
            this.contPeptideMap.get(pepKey).add(dr);
        } else if (isDecoy) {
            if (!this.decoyPeptideMap.containsKey(pepKey))
                this.decoyPeptideMap.put(pepKey, new ArrayList<>());
            this.decoyPeptideMap.get(pepKey).add(dr);
        } else {
            this.nTreatPsms++;
            if (!this.treatPeptideMap.containsKey(pepKey))
                this.treatPeptideMap.put(pepKey, new ArrayList<>());
            this.treatPeptideMap.get(pepKey).add(dr);
        }
    }

    /* This function will rotate the HashMap containing pepKey -> {ion intensity list}
        into one that will be ionKey -> {pepKey intensity list}
     */
    private void collapseHashMap(String target) {
        /* declare variables that will be used */
        HashMap<String, ArrayList<DiagnosticRecord>> pepMap;
        HashMap<Double, ArrayList<Double>> imm;
        HashMap<Double, ArrayList<Double>> capY;
        HashMap<Character, HashMap<Double,ArrayList<Double>>> squiggles;
        HashMap<Double, ArrayList<Double>> immN;
        HashMap<Double, ArrayList<Double>> capYN;
        HashMap<Character, HashMap<Double, ArrayList<Double>>> squigglesN;

        if (target.equals("control")) {
            pepMap = contPeptideMap;
            this.nControlPsms = pepMap.keySet().size();
            imm = immoniumX;
            capY = capYX;
            squiggles = squigglesX;
        } else if (target.equals("treatment")) {
            pepMap = treatPeptideMap;
            this.nTreatPsms = pepMap.keySet().size();
            imm = immoniumY;
            capY = capYY;
            squiggles = squigglesY;
        } else { // "decoy"
            pepMap = decoyPeptideMap;
            imm = immoniumDecoy;
            capY = capYDecoy;
            squiggles = squigglesDecoy;
        }

        //immoniumXN = {peak : <values>}
        for (String pepKey : pepMap.keySet()) {
            for (DiagnosticRecord dr : pepMap.get(pepKey)) { // This is only n = 1, can be fixed
                for (Double peak : dr.selectedImmoniumPeaks.keySet())
                    imm.get(peak).add(dr.selectedImmoniumPeaks.get(peak));
                for (Double peak : dr.selectedCapYPeaks.keySet())
                    capY.get(peak).add(dr.selectedCapYPeaks.get(peak));
                for (Character c : dr.selectedSquigglePeaks.keySet()) {
                    for (Double peak : dr.selectedSquigglePeaks.get(c).keySet())
                        squiggles.get(c).get(peak).add(dr.selectedSquigglePeaks.get(c).get(peak));
                }
            }
        }

    }



    private void collapseHashMaps() {
        collapseHashMap("control");
        collapseHashMap("treatment");
        collapseHashMap("decoy");
    }

    private double calcProportionWIon(ArrayList<Double> vals, int nPsms) {
        int n = 0;
        for (Double val : vals) {
            if (val > 0.0)
                n++;
        }
        if (n > 0 && nPsms > 0)
            return (double) n / (double) nPsms;
        else
            return 0;
    }

    private double calcWIonIntensity(ArrayList<Double> vals, int nPsms) {
        int n = 0;
        double intensity = 0.0;
        for (Double val: vals) {
            if (val > 0.0) {
                n++;
                intensity += val;
            }
        }
        if (n > 0)
            return intensity / (double) n;
        else
            return 0.0;
    }

    private double calcFoldChange(double propWMod, double wModInt, double propWModCont, double wModIntCont) {
        if (propWModCont == 0.0 || wModIntCont == 0.0)
            return 100;
        return (propWMod * wModInt) / (propWModCont * wModIntCont);
    }


    public void performTests() throws IOException {
        collapseHashMaps();
        MannWhitneyUTest mwu = new MannWhitneyUTest();

        this.immoniumTests = new ArrayList<>();
        for (Double peak : this.immoniumX.keySet()) {
            if (this.immoniumX.get(peak).size() < this.minPeps || this.immoniumY.get(peak).size() < this.minPeps)
                continue;
            double[] x = this.immoniumX.get(peak).stream().mapToDouble(i -> i).toArray();//.filter(i -> i > 0.0).toArray();
            double[] y = this.immoniumY.get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
            double p = mwu.mannWhitneyUTest(x, y);
            double u2 = mwu.mannWhitneyU2(x, y);
            long n1n2 = (long) this.immoniumX.get(peak).size() * (long) this.immoniumY.get(peak).size();
            double rankBiserCorr = u2 / (n1n2);
            //double rankBiserCorr = 1;
            double u1 = n1n2 - u2;
            boolean greaterThan = u2 > u1 ? true : false;
            if (this.twoTailedTests == false)
                p = convertP(p, greaterThan);
            p *= this.immoniumX.size();
            double propWIon = calcProportionWIon(this.immoniumY.get(peak), this.nTreatPsms);
            double wIonIntensity = calcWIonIntensity(this.immoniumY.get(peak), this.nTreatPsms);
            double propWIonCont = calcProportionWIon(this.immoniumX.get(peak), this.nControlPsms);
            double wIonIntensityCont = calcWIonIntensity(this.immoniumX.get(peak), this.nControlPsms);
            double foldChange = calcFoldChange(propWIon, wIonIntensity, propWIonCont, wIonIntensityCont);
            double maxAbsFoldChange = Math.max(foldChange, Math.pow(foldChange, -1));

            if (debug)
                System.out.println(peakApex + "\t" + peak + "\t" + foldChange + "\t" + propWIon);
            if (this.twoTailedTests == false) {
                if (p <= this.maxP && foldChange >= this.diagMinFoldChange && propWIon >= this.diagMinSpecDiff) {
                    this.immoniumTests.add(new Test(peakApex, peak, p, rankBiserCorr, false, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont,
                            u2, this.immoniumX.get(peak).size(), this.immoniumY.get(peak).size()));
                }
            } else {
                if (p <= this.maxP && maxAbsFoldChange >= this.diagMinFoldChange && propWIon >= this.diagMinSpecDiff ) {
                    this.immoniumTests.add(new Test(peakApex, peak, p, rankBiserCorr, false, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont,
                            u2, this.immoniumX.get(peak).size(), this.immoniumY.get(peak).size()));
                }
            }
        }
//        for (Double peak : this.immoniumX.keySet()) {
//            if (this.immoniumX.get(peak).size() < this.minPeps || this.immoniumDecoy.get(peak).size() < this.minPeps)
//                continue;
//            //double[] x = this.immoniumX.get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
//            //double[] y = this.immoniumDecoy.get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
//            double[] x = this.immoniumX.get(peak).stream().mapToDouble(i -> i).toArray();//.filter(i -> i > 0.0).toArray();
//            double[] y = this.immoniumDecoy.get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
//            double p = mwu.mannWhitneyUTest(x, y);
//            double u2 = mwu.mannWhitneyU2(x, y);
//            /*
//            double p;
//            double u2;
//            if (x.length > 1 && y.length > 1) {
//                p = tTest.tTest(x, y);
//                u2 = tTest.t(x, y);
//            } else {
//                p = 1;
//                u2 = 0;
//            }
//            */
//            //double rankBiserCorr = Math.abs((2.0 * uStat / (this.immoniumX.get(peak).size() * this.immoniumY.get(peak).size())) - 1);
//            long n1n2 = (long)x.length * (long)y.length;
//            double rankBiserCorr = u2 / (n1n2);
//            //double rankBiserCorr = 1;
//            double u1 = n1n2 - u2;
//            boolean greaterThan = u2 > u1 ? true : false;
//            if (this.twoTailedTests == false)
//                p = convertP(p, greaterThan);
//            p *= this.immoniumX.size();
//            double propWIon = calcProportionWIon(this.immoniumDecoy.get(peak), this.nTreatPsms);
//            double wIonIntensity = calcWIonIntensity(this.immoniumDecoy.get(peak), this.nTreatPsms);
//            double propWIonCont = calcProportionWIon(this.immoniumX.get(peak), this.nControlPsms);
//            double wIonIntensityCont = calcWIonIntensity(this.immoniumX.get(peak), this.nControlPsms);
//            double specDiff = Math.max(propWIon - propWIonCont, (propWIon*wIonIntensityCont - propWIonCont*wIonIntensityCont) / 100.0);
//            if (this.twoTailedTests == false) {
//                if (p <= this.maxP && (rankBiserCorr - 0.5) >= this.minRbc) {
//                    this.immoniumTests.add(new Test(peakApex, peak, p, rankBiserCorr, true, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
//                            this.immoniumX.get(peak).size(), this.immoniumDecoy.get(peak).size()));
//                }
//            } else {
//                if (p <= this.maxP && Math.abs(rankBiserCorr - 0.5) >= this.minRbc) {
//                    this.immoniumTests.add(new Test(peakApex, peak, p, rankBiserCorr, true, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
//                            this.immoniumX.get(peak).size(), this.immoniumDecoy.get(peak).size()));
//                }
//            }
//
//            //System.out.printf("Immonium %.04f\t%e\t%f\n", peak, p, rankBiserCorr);
//        }

        this.capYTests = new ArrayList<>();
        for (Double peak : this.capYX.keySet()) {
            if (this.capYX.get(peak).size() < this.minPeps || this.capYY.get(peak).size() < this.minPeps)
                continue;
            double[] x = this.capYX.get(peak).stream().mapToDouble(i -> i).toArray();
            double[] y = this.capYY.get(peak).stream().mapToDouble(i -> i).toArray();
            double p = mwu.mannWhitneyUTest(x, y);
            double u2 = mwu.mannWhitneyU2(x, y);

            long n1n2 = (long)x.length * (long)y.length;
            double rankBiserCorr = u2 / (n1n2);
            double u1 = n1n2 - u2;
            boolean greaterThan = u2 > u1 ? true : false;
            if (this.twoTailedTests == false)
                p = convertP(p, greaterThan);
            p *= this.capYX.size();
            double propWIon = calcProportionWIon(this.capYY.get(peak), this.nTreatPsms);
            double wIonIntensity = calcWIonIntensity(this.capYY.get(peak), this.nTreatPsms);
            double propWIonCont = calcProportionWIon(this.capYX.get(peak), this.nControlPsms);
            double wIonIntensityCont = calcWIonIntensity(this.capYX.get(peak), this.nControlPsms);
            double specDiff = propWIon;
            double foldChange = calcFoldChange(propWIon, wIonIntensity, propWIonCont, wIonIntensityCont);
            double maxAbsFoldChange = Math.max(foldChange, Math.pow(foldChange, -1));
            if (this.twoTailedTests == false) {
                if (p <= this.maxP && propWIon > this.pepMinSpecDiff && foldChange >= this.pepMinFoldChange) {
                    this.capYTests.add(new Test(peakApex, peak, p, rankBiserCorr, false, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                            this.capYX.get(peak).size(), this.capYY.get(peak).size()));
                }
            } else {
                if (p <= this.maxP && propWIon > this.pepMinSpecDiff && maxAbsFoldChange >= this.pepMinFoldChange) {
                    this.capYTests.add(new Test(peakApex, peak, p, rankBiserCorr, false, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                            this.capYX.get(peak).size(), this.capYY.get(peak).size()));
                }
            }
        }
//        for (Double peak : this.capYX.keySet()) {
//            if (this.capYX.get(peak).size() < this.minPeps || this.capYDecoy.get(peak).size() < this.minPeps)
//                continue;
//            //double[] x = this.capYX.get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
//            //double[] y = this.capYDecoy.get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
//            double[] x = this.capYX.get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
//            double[] y = this.capYDecoy.get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
//            double p = mwu.mannWhitneyUTest(x, y);
//            double u2 = mwu.mannWhitneyU2(x, y);
//            /*
//            double p;
//            double u2;
//            if (x.length > 1 && y.length > 1) {
//                p = tTest.tTest(x, y);
//                u2 = tTest.t(x, y);
//            } else {
//                p = 1;
//                u2 = 0;
//            }
//            */
//            long n1n2 = (long) x.length * (long) y.length;
//            double rankBiserCorr = u2 / (n1n2);
//            //double rankBiserCorr = 1;
//            double u1 = n1n2 - u2;
//            boolean greaterThan = u2 > u1 ? true : false;
//            if (this.twoTailedTests == false)
//                p = convertP(p, greaterThan);
//            //double rankBiserCorr = Math.abs((2.0 * uStat / (this.capYX.get(peak).size() * this.capYY.get(peak).size())) - 1);
//            p *= this.capYX.size();
//            double propWIon = calcProportionWIon(this.capYDecoy.get(peak), this.nTreatPsms);
//            double wIonIntensity = calcWIonIntensity(this.capYDecoy.get(peak), this.nTreatPsms);
//            double propWIonCont = calcProportionWIon(this.capYX.get(peak), this.nControlPsms);
//            double wIonIntensityCont = calcWIonIntensity(this.capYX.get(peak), this.nControlPsms);
//            double specDiff = Math.max(propWIon - propWIonCont, (propWIon*wIonIntensityCont - propWIonCont*wIonIntensityCont) / 100.0);
//            if (this.twoTailedTests == false) {
//                if (p <= this.maxP && (rankBiserCorr - 0.5) >= this.minRbc) {
//                    this.capYTests.add(new Test(peakApex, peak, p, rankBiserCorr, true, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
//                            this.capYX.get(peak).size(), this.capYDecoy.get(peak).size()));
//                }
//            } else {
//                if (p <= this.maxP && Math.abs(rankBiserCorr - 0.5) >= this.minRbc) {
//                    this.capYTests.add(new Test(peakApex, peak, p, rankBiserCorr, true, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
//                            this.capYX.get(peak).size(), this.capYDecoy.get(peak).size()));
//                }
//            }
//            //System.out.printf("CapY %.04f\t%e\t%f\n", peak, p, rankBiserCorr);
//        }

        this.squigglesTests = new HashMap<>();
        for (Character c : this.squigglesX.keySet()) {
            this.squigglesTests.put(c, new ArrayList<>());
            for (Double peak : this.squigglesX.get(c).keySet()) {
                if (this.squigglesX.get(c).get(peak).size() < this.minPeps || this.squigglesY.get(c).get(peak).size() < this.minPeps)
                    continue;
                double[] x = this.squigglesX.get(c).get(peak).stream().mapToDouble(i -> i).toArray();
                double[] y = this.squigglesY.get(c).get(peak).stream().mapToDouble(i -> i).toArray();
                double p = mwu.mannWhitneyUTest(x, y);
                double u2 = mwu.mannWhitneyU2(x, y);
                long n1n2 = (long)x.length * (long)y.length;
                double rankBiserCorr = u2 / (n1n2);
                double u1 = n1n2 - u2;
                boolean greaterThan = u2 > u1 ? true : false;
                if (this.twoTailedTests == false)
                    p = convertP(p, greaterThan);
                p *= this.squigglesX.get(c).size();

                double propWIon = calcProportionWIon(this.squigglesY.get(c).get(peak), this.nTreatPsms);
                double wIonIntensity = calcWIonIntensity(this.squigglesY.get(c).get(peak), this.nTreatPsms);
                double propWIonCont = calcProportionWIon(this.squigglesX.get(c).get(peak), this.nControlPsms);
                double wIonIntensityCont = calcWIonIntensity(this.squigglesX.get(c).get(peak), this.nControlPsms);

                double specDiff = propWIon;
                double foldChange = calcFoldChange(propWIon, wIonIntensity, propWIonCont, wIonIntensityCont);
                double maxAbsFoldChange = Math.max(foldChange, Math.pow(foldChange, -1));

                if (this.twoTailedTests == false) {
                    if (p <= this.maxP && propWIon > this.fragMinSpecDiff && foldChange >= this.fragMinFoldChange) {
                        this.squigglesTests.get(c).add(new Test(peakApex, peak, p, rankBiserCorr, false, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                                this.squigglesX.get(c).get(peak).size(), this.squigglesY.get(c).get(peak).size()));
                    }
                } else {
                    if (p <= this.maxP && propWIon > this.fragMinSpecDiff && maxAbsFoldChange >= this.fragMinFoldChange) {
                        this.squigglesTests.get(c).add(new Test(peakApex, peak, p, rankBiserCorr, false, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                                this.squigglesX.get(c).get(peak).size(), this.squigglesY.get(c).get(peak).size()));
                    }
                }
            }
        }
        for (Character c : this.squigglesX.keySet()) {

            for (Double peak : this.squigglesX.get(c).keySet()) {
                if (this.squigglesX.get(c).get(peak).size() < this.minPeps || this.squigglesDecoy.get(c).get(peak).size() < this.minPeps)
                    continue;
                //double[] x = this.squigglesX.get(c).get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
                //double[] y = this.squigglesDecoy.get(c).get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
                double[] x = this.squigglesX.get(c).get(peak).stream().mapToDouble(i -> i).toArray();//.filter(i -> i > 0.0).toArray();
                double[] y = this.squigglesDecoy.get(c).get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
                double p = mwu.mannWhitneyUTest(x, y);
                double u2 = mwu.mannWhitneyU2(x, y);
                /*
                double p;
                double u2;
                if (x.length > 1 && y.length > 1) {
                    p = tTest.tTest(x, y);
                    u2 = tTest.t(x, y);
                } else {
                    p = 1;
                    u2 = 0;
                }
                 */
                long n1n2 = (long)x.length * (long)y.length;
                double rankBiserCorr = u2 / (n1n2);
                //double rankBiserCorr = 1;
                double u1 = n1n2 - u2;
                boolean greaterThan = u2 > u1 ? true : false;
                if (this.twoTailedTests == false)
                    p = convertP(p, greaterThan);
                p *= this.squigglesX.get(c).keySet().size();
                //double rankBiserCorr = Math.abs((2.0 * uStat / (this.squigglesX.get(c).get(peak).size() * this.squigglesY.get(c).get(peak).size())) - 1);
                //out.printf("%.04f\t%e\t%f\n", peak, p, rankBiserCorr);
                p *= this.squigglesX.get(c).size();
                double propWIon = calcProportionWIon(this.squigglesDecoy.get(c).get(peak), this.nTreatPsms);
                double wIonIntensity = calcWIonIntensity(this.squigglesDecoy.get(c).get(peak), this.nTreatPsms);
                double propWIonCont = calcProportionWIon(this.squigglesX.get(c).get(peak), this.nControlPsms);
                double wIonIntensityCont = calcWIonIntensity(this.squigglesX.get(c).get(peak), this.nControlPsms);
                if (this.twoTailedTests == false) {
                    if (p <= this.maxP && (rankBiserCorr - 0.5) >= this.minRbc) {
                        this.squigglesTests.get(c).add(new Test(peakApex, peak, p, rankBiserCorr, true, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                                this.squigglesX.get(c).get(peak).size(), this.squigglesDecoy.get(c).get(peak).size()));
                    }
                } else {
                    if (p <= this.maxP && Math.abs(rankBiserCorr - 0.5) >= this.minRbc) {
                        this.squigglesTests.get(c).add(new Test(peakApex, peak, p, rankBiserCorr, true, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                                this.squigglesX.get(c).get(peak).size(), this.squigglesDecoy.get(c).get(peak).size()));
                    }
                }
                //if (p * this.squigglesX.get(c).get(peak).size() < 0.05 && rankBiserCorr > 0.5)
                //System.out.printf("%.04f\t%e\t%f\n", peak, p, rankBiserCorr);
            }
            //out.close();
        }

        relocalizeDeltas(1, this.specTol); //todo tol and nAdjacentAAs??
        //System.out.println("immonium");
        collapseTests(this.immoniumTests, 0.01, false); //todo tol
        //System.out.println("capY");
        collapseTests(this.capYTests, 0.05, false); //todo tol
        calibrateTests(this.capYTests, 0.05, false);
        //System.out.println("squiggle");
        for (Character c : this.squigglesTests.keySet()) {
            collapseTests(this.squigglesTests.get(c), 0.05, true); //todo tol
            //findRemainderMassesCutoff(this.squigglesTests.get(c), 0.05);
            calibrateTests(this.squigglesTests.get(c), 0.05, true);
            findRemainderMassesCutoff(this.squigglesTests.get(c), 0.05);
        }

        /* clear out most of memory */
        clearMemory();

        //Print all tests

        if (debug) {
            Collections.sort(this.immoniumTests);
            for (Test t : this.immoniumTests)
                System.out.printf("Immonium %.04f\t%e\t%.04f\t%b\t%d\t%d\n", t.mass, t.q, t.u, t.isDecoy, t.n1, t.n2);
            Collections.sort(this.capYTests);
            for (Test t : this.capYTests)
                System.out.printf("CapY %.04f\t%e\t%.04f\t%b\t%d\t%d\n", t.mass, t.q, t.u, t.isDecoy, t.n1, t.n2);
            for (Character c : this.squigglesTests.keySet()) {
                Collections.sort(this.squigglesTests.get(c));
                for (Test t : this.squigglesTests.get(c))
                    System.out.printf("Squiggle %s %.04f\t%e\t%.04f\t%b\t%d\t%d\n", Character.toString(c), t.mass, t.q, t.u, t.isDecoy, t.n1, t.n2);
            }
        }
    }

    /* Converts two-tailed p-value to one-tailed p-value
    * Hipparchus only provides a two-tailed version of MWU
    * Currently using two-tailed tests for everything
    * */
    private double convertP(double oldPVal, boolean greaterThan) {
        double p;
        /*
        if (greaterThan)
            p = oldPVal / 2.0;
        else
            p = 1.0 - (oldPVal / 2.0);
         */
        return oldPVal;
    }

    private void clearMemory() {
        this.contPeptideMap = null;
        this.treatPeptideMap = null;
        this.decoyPeptideMap = null;
        this.immoniumX = null;
        this.immoniumY = null;
        this.immoniumDecoy = null;
        this.capYX = null;
        this.capYY = null;
        this.capYDecoy = null;
        this.squigglesX = null;
        this.squigglesY = null;
        this.squigglesDecoy = null;
    }

    private void _collapseTests(ArrayList<Test> tests, double tol, boolean adjustedMass) {
        Collections.sort(tests, new Comparator<Test>() {
            @Override
            public int compare(Test o1, Test o2) {
                return Double.compare(o1.q, o2.q);
            }
        });

        // First, group tests based on isotopic status (0-3 isotpe error, mass 1.00235)
        int group = 1;
        for (int i = 0; i < tests.size(); i++) {
            if (tests.get(i).group > 0)
                continue;
            tests.get(i).group = group;
            for (int j = i + 1; j < tests.size(); j++) {
                if (tests.get(j).group > 0)
                    continue;
                for (int c = -1; c < 4; c++) {
                    double c13 = 1.00235;
                    double shiftMass = c13 * c;
                    if (adjustedMass) {
                        if (Math.abs(tests.get(i).adjustedMass - (tests.get(j).adjustedMass - shiftMass)) < tol) {
                            tests.get(j).group = group;
                            if (c == 0)
                                tests.get(j).isRedundant = true;
                            else
                                tests.get(j).isIsotopeRep = false;
                        }
                    } else {
                        if (Math.abs(tests.get(i).mass - (tests.get(j).mass - shiftMass)) < tol) {
                            tests.get(j).group = group;
                            tests.get(j).isGroupRep = false;
                        }
                    }
                }
            }
            group++;
        }
    }

    private void collapseTests(ArrayList<Test> tests, double tol, boolean adjustedMass) {
        // First, filter out reudndant tests by selecting most sigfnicant ions
        Collections.sort(tests, new Comparator<Test>() {
            @Override
            public int compare(Test o1, Test o2) {
                return Double.compare(o1.q, o2.q);
            } // Sorting by most to least significant
        });

        int group = 1;
        for (int i = 0; i < tests.size(); i++) {
            if (tests.get(i).group > 0)
                continue;
            tests.get(i).group = group;
            //System.out.println(tests.get(i).mass);
            //todo this can give a suboptimal result if the monoisotopic peak is not the most significant
            //todo search tests recursively instead of in sequence
            for (int j = i + 1; j < tests.size(); j++) {
                if (tests.get(j).group > 0)
                    continue;
                for (int c = -1; c < 4; c++) {
                    double c13 = 1.00235;
                    double shiftMass = c13 * c;
                    if (adjustedMass) {
                        if (Math.abs(tests.get(i).adjustedMass - (tests.get(j).adjustedMass - shiftMass)) < tol) {
                            tests.get(j).group = group;
                            if (c == 0)
                                tests.get(j).isRedundant = true;
                        }
                    } else {
                        if (Math.abs(tests.get(i).mass - (tests.get(j).mass - shiftMass)) < tol) {
                            tests.get(j).group = group;
                            if (c == 0)
                                tests.get(j).isRedundant = true;
                        }
                    }
                }
            }
            group++;
        }

        // Second, find monoisotopc peak by looking at most instense ion
        Collections.sort(tests, new Comparator<Test>() {
            @Override
            public int compare(Test o1, Test o2) {
                int g1 = o1.group;
                int g2 = o2.group;
                int firstComp = Integer.compare(g1, g2);

                if (firstComp != 0)
                    return firstComp;

                double int1 = o1.propWIonTreat * o1.propWIonIntensity;
                double int2 = o2.propWIonTreat * o2.propWIonIntensity;
                return -1 * Double.compare(int1, int2);
            }
        });

        int cGroup = 0;
        for (int i = 0; i < tests.size(); i++) {
            if (tests.get(i).isRedundant)
                continue;
            if (tests.get(i).group == cGroup)
                continue;
            else {
                tests.get(i).isIsotopeRep = true;
                cGroup = tests.get(i).group;
            }
        }

    }

    private void calibrateTests(ArrayList<Test> tests, double tol, boolean adjustedMass) {
        for (Test t: tests)
            t.calibrateMass(adjustedMass, tol);
    }

    /* Checks for enrichment of adjacent AAs for each squiggle peak, then adjusts mass appropriately */
    private void relocalizeDeltas(int nAdjacentAas, double tol) {
        //System.out.println("Relocalizing deltas");
        for (Character c : this.squigglesTests.keySet()) {
            //System.out.println("**" + c);
            int[][] shiftedCounts = new int[nAdjacentAas * 2 + 1][26];
            for (Test t : this.squigglesTests.get(c)) {
                //System.out.println(t.mass);
                int[] resCounts2 = new int[26];
                int[][] resCounts = new int[nAdjacentAas * 2 + 1][26]; //todo tmp
                int[][] bgResCounts = new int[nAdjacentAas * 2 + 1][26];
                int[] totalChances = new int[nAdjacentAas*2+1];
                boolean reassigned = false;
                int totalMatchedIons = 0;
                for (String pepKey : this.treatPeptideMap.keySet()) {
                    ArrayList<int[]> ionCounts = new ArrayList<>();
                    String pepSeq = "";
                    for (DiagnosticRecord dr : this.treatPeptideMap.get(pepKey)) {
                        /* Get and normalize scores */
                        int[] scores = dr.localizeRemainderMass((float) t.mass, c);
                        pepSeq = dr.pepSeq;
                        ionCounts.add(scores);
                    }
                    /* Add pepKey's modal ion counts to total relative positions */
                    int[] ionModes = calculateIonModes(ionCounts);
                    for (int i = 0; i < ionModes.length; i++) {
                        for (int cRelPos = -1 * nAdjacentAas; cRelPos <= nAdjacentAas; cRelPos++) {
                            int cAbsPos = i + cRelPos;
                            if (cAbsPos < 0 || cAbsPos > ionModes.length - 1)
                                continue;
                            char cres = Character.toLowerCase(pepSeq.charAt(cAbsPos));
                            //System.out.println(cres);
                            resCounts[cRelPos+nAdjacentAas][cres - 'a'] += ionModes[i];
                        }
                        totalMatchedIons += ionModes[i];
                    }

                    /* Add background counts to total */
                    int[][] backgroundCounts = calculateResidueBackground(pepSeq, c, nAdjacentAas);
                    for (int i = 0; i < backgroundCounts.length; i++) {
                        for (int j = 0; j < backgroundCounts[i].length; j++) {
                            bgResCounts[i][j] += backgroundCounts[i][j];
                            totalChances[i] += backgroundCounts[i][j];
                        }
                    }
                }
                double newRemainderMass = t.mass;
                double bestResP = -1;
                /* Merge leucine and isoleucine onto leucine*/
                for (int i = 0; i < resCounts.length; i++) {
                    resCounts[i]['l' - 'a'] += resCounts[i]['i' - 'a'];
                    resCounts[i]['i' - 'a'] = 0;
                }

                //Adjust mass by res mass if it is very prevalent
                /*
                for (int i  = 0; i < resCounts2.length; i++) {
                    if ((double) resCounts2[i] / resCounts2Sum > 0.75) {
                        newRemainderMass = calcNewRemainderMass(t.mass, c, 0, i);
                    }
                }
                */

                for (int i = 0; i < resCounts.length; i++) {
                    for (int j = 0; j < resCounts[i].length; j++) {
                        if (resCounts[i][j] > 0.0) {
                            char cRes = (char) (j + 'a');
                            int cResLoc = resCounts[i][j];
                            int cResBg = bgResCounts[i][j];
                            int cResNotLoc = cResBg - cResLoc;
                            int notCResLoc = totalMatchedIons - cResLoc;
                            int notCResBg = totalChances[i] - cResBg;
                            int notCResNotLoc = notCResBg - notCResLoc;
                            double resP = (double) resCounts[i][j] / (double) totalMatchedIons;

                            //System.out.println(cRes + "\t" + cResLoc + "\t" + resP);

                            if (resP > 0.5 && resP > bestResP) {
                                newRemainderMass = calcNewRemainderMass(t.mass, c, i-nAdjacentAas, j);
                                //System.out.println(t.mass + "\t" + newRemainderMass + "\t" + c + "\t" + (i-nAdjacentAas) + "\t" + j);
                                bestResP = resP;
                                shiftedCounts[i][j]++;
                            }
                        }
                    }
                }

                t.adjustedMass = newRemainderMass;
            }
        }
    }

    private int[] calculateIonModes(ArrayList<int[]> ionCounters) {
        int pepLen = ionCounters.get(0).length;

        /* First pass finds max value at each peptide position */
        int[] maxIonsByPepPos = new int[pepLen];
        for (int i = 0; i < ionCounters.size(); i++) {
            for (int j = 0; j < pepLen; j++) {
                if (ionCounters.get(i)[j] > maxIonsByPepPos[j])
                    maxIonsByPepPos[j] = ionCounters.get(i)[j];
            }
        }

        /* Second pass calculates modal ionCount for each position */
        int[] ionModes = new int[pepLen];
        for (int i = 0; i < maxIonsByPepPos.length; i++) {
            /* For position i on peptide, count appearance for each nFrags */
            int[] ionCount = new int[maxIonsByPepPos[i]+1]; // Each position contains # of appearances for each nFrags
            for (int j = 0; j < ionCounters.size(); j++)
                ionCount[ionCounters.get(j)[i]]++;

            /* Find modal nFrags for peptide position i */
            int cMaxFragCounts = -1;
            int cMode = -1;
            for (int k = 0; k < ionCount.length; k++) {
                if (ionCount[k] > cMaxFragCounts) {
                    cMaxFragCounts = ionCount[k];
                    cMode = k;
                }
            }
            ionModes[i] = cMode;
        }

        return ionModes;
    }

    private int[][] calculateResidueBackground(String pepSeq, char ionType, int nAdjacentAas) {
        int pepLen = pepSeq.length();
        int[][] resCounts = new int[nAdjacentAas * 2 + 1][26];
        if (ionType == 'a' || ionType == 'b' || ionType == 'c') {
            for (int i = 0; i < pepLen - 1; i++) {
                for (int cRelPos = -1 * nAdjacentAas; cRelPos <= nAdjacentAas; cRelPos++) {
                    int cAbsPos = i + cRelPos;
                    if (cAbsPos < 0 || cAbsPos > pepLen - 1)
                        continue;
                    char cres = Character.toLowerCase(pepSeq.charAt(cAbsPos));
                    resCounts[cRelPos + nAdjacentAas][cres - 'a']++;
                }
            }
        }
        if (ionType == 'x' || ionType == 'y' || ionType == 'z') {
            for (int i = 0; i < pepLen - 1; i++) {
                for (int cRelPos = -1 * nAdjacentAas; cRelPos <= nAdjacentAas; cRelPos++) {
                    int cAbsPos = (pepLen - 1 - i) + cRelPos;
                    if (cAbsPos < 0 || cAbsPos > pepLen - 1)
                        continue;
                    char cres = Character.toLowerCase(pepSeq.charAt(cAbsPos));
                    resCounts[cRelPos + nAdjacentAas][cres - 'a']++;
                }
            }
        }

        return resCounts;
    }

    private double calcNewRemainderMass(double dmass, char ionType, int relativePosition, int res) {
        //double newDmass = 100000000; //todo some errors on edge cases (doesn't change results)
        double newDmass = dmass;

        if (ionType > 'n') { // c-term
            if (relativePosition < 0) // to the left
                newDmass = dmass - AAMasses.monoisotopic_masses[res];
            else if (relativePosition == 0) // to the right
                newDmass = dmass + AAMasses.monoisotopic_masses[res];
        } else if (ionType < 'n') {
            if (relativePosition == 0)
                newDmass = dmass + AAMasses.monoisotopic_masses[res];
            else if (relativePosition > 0)
                newDmass = dmass - AAMasses.monoisotopic_masses[res];
        }

        //if (newDmass > 10000) { //todo some errors on edge cases (doesn't change results)
        //    System.out.println(dmass + "\t" + ionType + "\t" + relativePosition + "\t" + res);
        //}
        return newDmass;
    }

    /*
    Find the first mass that was adjusted because an adjacent residue in enriched
     */
    private void findRemainderMassesCutoff(ArrayList<Test> tests, double tol) {
        int breakPoint = -1;
        for (int i  = 0; i < Math.min(tests.size(), 10); i++) {
            for (int j  = i + 1; j < Math.min(tests.size(), 10); j++) {
                if (breakPoint > -1)
                    break;
                if (Math.abs(tests.get(i).mass - tests.get(j).mass) > 5 && // Make sure it isn't an isotopic collapse
                        Math.abs(tests.get(i).adjustedMass - tests.get(j).adjustedMass) < 0.5) { // If it is the same mass
                    breakPoint = j - 1;
                }
            }
            if (breakPoint > -1)
                break;
        }
        for (int i  = 0; i < tests.size(); i++) {
            if (i >= 10)
                tests.get(i).isValid = false;
            else if (breakPoint == -1)
                tests.get(i).isValid = true;
            else
                if (i <= breakPoint)
                    tests.get(i).isValid = true;
                else
                    tests.get(i).isValid = false;
        }
    }
}

class Test implements Comparable<Test> {
    public double peakApex;
    public double mass;
    public double adjustedMass;
    public double q;
    public double rbc;
    public boolean isDecoy;
    public boolean isValid;
    public int group;
    public boolean isGroupRep;
    public double u;
    public long n1;
    public long n2;
    public double propWIonTreat;
    public double propWIonCont;
    public double propWIonIntensity;
    public double propWIonIntensityCont;

    public boolean isRedundant;
    public boolean isIsotopeRep;

    Test(double peakApex, double mass, double q, double rbc, boolean isDecoy, double propWIonTreat, double propWIonCont,
         double propWIonIntensity, double propWIonIntensityCont, double u, long n1, long n2) {
        this.peakApex = peakApex;
        this.mass = mass;
        this.adjustedMass = mass;
        this.q = q;
        this.rbc = rbc;
        this.isDecoy = isDecoy;
        this.propWIonTreat = propWIonTreat;
        this.propWIonCont = propWIonCont;
        this.propWIonIntensity = propWIonIntensity;
        this.propWIonIntensityCont = propWIonIntensityCont;
        this.u = u;
        this.group = 0;
        this.isGroupRep = false;
        this.isIsotopeRep = false;
        this.isRedundant = false;
        this.n1 = n1;
        this.n2 = n2;
    }

    public void calibrateMass(boolean adjustedMass, double tol) { //TODO is this only correcting some Y ions?
        if (adjustedMass) {
            //if (Math.abs(this.peakApex - 2529.8778) < 0.1)
            //    System.out.println(this.adjustedMass + "\t" + tol);
            if (Math.abs(this.adjustedMass) < tol) {//shift pep remainder to 0 if present
                //if (Math.abs(this.peakApex - 2529.8778) < 0.1)
                //    System.out.println("**" + this.adjustedMass + "\t" + tol);
                this.adjustedMass = 0;
                //if (Math.abs(this.peakApex - 2529.8778) < 0.1)
                //    System.out.println("**" + this.adjustedMass + "\t" + tol);
            }
            else if (Math.abs(this.adjustedMass - this.peakApex) < tol)
                this.adjustedMass = this.peakApex;
        } else {
            if (Math.abs(this.mass - this.peakApex) < tol)
                this.adjustedMass = this.peakApex;
        }

    }

    public String toString() {
        return String.format("%.04f\t%e\t%.04f\t%b\t%d\t%d", this.mass, this.q, this.u, this.isDecoy, this.n1, this.n2);
    }

    @Override
    public int compareTo(Test arg0) {
        return Double.compare(this.adjustedMass, arg0.adjustedMass);
    }

}
