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
    HashMap<Double, ArrayList<Double>> immoniumXN;
    HashMap<Double, ArrayList<Double>> immoniumYN;
    HashMap<Double, ArrayList<Double>> immoniumDecoyN;
    HashMap<Double, ArrayList<Double>> capYXN;
    HashMap<Double, ArrayList<Double>> capYYN;
    HashMap<Double, ArrayList<Double>> capYDecoyN;
    HashMap<Character, HashMap<Double, ArrayList<Double>>> squigglesXN;
    HashMap<Character, HashMap<Double, ArrayList<Double>>> squigglesYN;
    HashMap<Character, HashMap<Double, ArrayList<Double>>> squigglesDecoyN;

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
    boolean twoTailedTests;
    int minPeps;

    double specTol;

    public PeakCompareTester(ArrayList<Double> xVals, ArrayList<Double> yVals) {
        this.x = xVals.stream().mapToDouble(i -> i).toArray();
        this.y = yVals.stream().mapToDouble(i -> i).toArray();
    }

    public PeakCompareTester(double peakApex, ArrayList<Double> unifImm, ArrayList<Double> unifCapY, HashMap<Character, ArrayList<Double>> unifSquig, double maxP, double minRbc, boolean twoTailedTests, double specTol) {
        this.peakApex = peakApex;
        this.immoniumX = new HashMap<>();
        this.immoniumY = new HashMap<>();
        this.immoniumDecoy = new HashMap<>();
        this.immoniumXN = new HashMap<>();
        this.immoniumYN = new HashMap<>();
        this.immoniumDecoyN = new HashMap<>();
        for (Double peak : unifImm) {
            this.immoniumX.put(peak, new ArrayList<>());
            this.immoniumY.put(peak, new ArrayList<>());
            this.immoniumDecoy.put(peak, new ArrayList<>());
            this.immoniumXN.put(peak, new ArrayList<>());
            this.immoniumYN.put(peak, new ArrayList<>());
            this.immoniumDecoyN.put(peak, new ArrayList<>());
        }
        this.capYX = new HashMap<>();
        this.capYY = new HashMap<>();
        this.capYDecoy = new HashMap<>();
        this.capYXN = new HashMap<>();
        this.capYYN = new HashMap<>();
        this.capYDecoyN = new HashMap<>();
        for (Double peak : unifCapY) {
            this.capYX.put(peak, new ArrayList<>());
            this.capYY.put(peak, new ArrayList<>());
            this.capYDecoy.put(peak, new ArrayList<>());
            this.capYXN.put(peak, new ArrayList<>());
            this.capYYN.put(peak, new ArrayList<>());
            this.capYDecoyN.put(peak, new ArrayList<>());
        }
        this.squigglesX = new HashMap<>();
        this.squigglesY = new HashMap<>();
        this.squigglesDecoy = new HashMap<>();
        this.squigglesXN = new HashMap<>();
        this.squigglesYN = new HashMap<>();
        this.squigglesDecoyN = new HashMap<>();
        for (Character ion : unifSquig.keySet()) {
            this.squigglesX.put(ion, new HashMap<>());
            this.squigglesY.put(ion, new HashMap<>());
            this.squigglesDecoy.put(ion, new HashMap<>());
            this.squigglesXN.put(ion, new HashMap<>());
            this.squigglesYN.put(ion, new HashMap<>());
            this.squigglesDecoyN.put(ion, new HashMap<>());
            for (Double peak : unifSquig.get(ion)) {
                this.squigglesX.get(ion).put(peak, new ArrayList<>());
                this.squigglesY.get(ion).put(peak, new ArrayList<>());
                this.squigglesDecoy.get(ion).put(peak, new ArrayList<>());
                this.squigglesXN.get(ion).put(peak, new ArrayList<>());
                this.squigglesYN.get(ion).put(peak, new ArrayList<>());
                this.squigglesDecoyN.get(ion).put(peak, new ArrayList<>());
            }
        }
        this.contPeptideMap = new HashMap<>();
        this.treatPeptideMap = new HashMap<>();
        this.decoyPeptideMap = new HashMap<>();
        this.nControlPsms = 0; //actually pepkey count
        this.nTreatPsms = 0; //actually pepkey count

        this.maxP = maxP;
        this.minRbc = minRbc;
        this.twoTailedTests = twoTailedTests;
        this.minPeps = Integer.parseInt(PTMShepherd.getParam("diagmine_minPeps"));

        this.specTol = specTol;
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


    private void collapseHashMap(String target) {

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
            immN = immoniumXN;
            capYN = capYXN;
            squigglesN = squigglesXN;
        } else if (target.equals("treatment")) {
            pepMap = treatPeptideMap;
            this.nTreatPsms = pepMap.keySet().size();
            imm = immoniumY;
            capY = capYY;
            squiggles = squigglesY;
            immN = immoniumYN;
            capYN = capYYN;
            squigglesN = squigglesYN;
        } else { // "decoy"
            pepMap = decoyPeptideMap;
            imm = immoniumDecoy;
            capY = capYDecoy;
            squiggles = squigglesDecoy;
            immN = immoniumDecoyN;
            capYN = capYDecoyN;
            squigglesN = squigglesDecoyN;
        }

        /* Collapse pepkeys in pep map */
        for (String pepKey : pepMap.keySet()) {
            /* Init normalizing params and temp objects, this is for a single pepKey */
            int nPsms = pepMap.get(pepKey).size();
            HashMap<Double, Double> immPeaks = new HashMap<>(); // MS2 Peak Apex, Avg intensity
            HashMap<Double, AtomicInteger> immNPeaks = new HashMap<>(); // MS2 Peak Apex, n spec nonzero intensity
            HashMap<Double, AtomicDouble> immNSumPeaks = new HashMap<>(); // MS2 Peak Apex, sum of specs with nonzero intensity
            HashMap<Double, Double> capYPeaks = new HashMap<>();
            HashMap<Double, AtomicInteger> capYNPeaks = new HashMap<>();
            HashMap<Double, AtomicDouble> capYNSumPeaks = new HashMap<>();
            HashMap<Character, HashMap<Double, Double>> squigglePeaks = new HashMap<>();
            HashMap<Character, HashMap<Double, AtomicInteger>> squiggleNPeaks = new HashMap<>();
            HashMap<Character, HashMap<Double, AtomicDouble>> squiggleNSumPeaks = new HashMap<>();
            for (Double peak : imm.keySet()) {
                immPeaks.put(peak, 0.0);
                immNPeaks.put(peak, new AtomicInteger());
                immNSumPeaks.put(peak, new AtomicDouble());
            }
            for (Double peak : capY.keySet()) {
                capYPeaks.put(peak, 0.0);
                capYNPeaks.put(peak, new AtomicInteger());
                capYNSumPeaks.put(peak, new AtomicDouble());
            }
            for (Character c : squiggles.keySet()) {
                squigglePeaks.put(c, new HashMap<>());
                squiggleNPeaks.put(c, new HashMap<>());
                squiggleNSumPeaks.put(c, new HashMap<>());
                for (Double peak : squiggles.get(c).keySet()) {
                    squigglePeaks.get(c).put(peak, 0.0);
                    squiggleNPeaks.get(c).put(peak, new AtomicInteger());
                    squiggleNSumPeaks.get(c).put(peak, new AtomicDouble());
                }
            }
            /* Collapse into mean of respective vals */
            for (DiagnosticRecord dr : pepMap.get(pepKey)) {
                for (Double peak : dr.selectedImmoniumPeaks.keySet()) {
                    double cVal = immPeaks.get(peak);
                    cVal += dr.selectedImmoniumPeaks.get(peak) / nPsms;
                    immPeaks.put(peak, cVal);
                    if (dr.selectedImmoniumPeaks.get(peak) > 0.00000000001) { //count PSMs with ions
                        immNPeaks.get(peak).getAndIncrement();
                        immNSumPeaks.get(peak).getAndAdd(dr.selectedImmoniumPeaks.get(peak));
                    }
                }
                for (Double peak : dr.selectedCapYPeaks.keySet()) {
                    double cVal = capYPeaks.get(peak);
                    cVal += dr.selectedCapYPeaks.get(peak) / nPsms;
                    capYPeaks.put(peak, cVal);
                    if (dr.selectedCapYPeaks.get(peak) > 0.00000000001) {//count PSMs with ions
                        capYNPeaks.get(peak).getAndIncrement();
                        capYNSumPeaks.get(peak).getAndAdd(dr.selectedCapYPeaks.get(peak));
                    }
                }
                for (Character c : dr.selectedSquigglePeaks.keySet()) {
                    for (Double peak : dr.selectedSquigglePeaks.get(c).keySet()) {
                        double cVal = squigglePeaks.get(c).get(peak);
                        cVal += dr.selectedSquigglePeaks.get(c).get(peak) / nPsms;
                        squigglePeaks.get(c).put(peak, cVal);
                        if (dr.selectedSquigglePeaks.get(c).get(peak) > 0.00000000001) {//count PSMs with ions
                            squiggleNPeaks.get(c).get(peak).getAndIncrement();
                            squiggleNSumPeaks.get(c).get(peak).getAndAdd(dr.selectedSquigglePeaks.get(c).get(peak));
                        }
                    }
                }
            }
            /* Add vals to class lists */
            for (Double peak : immPeaks.keySet()) {
                imm.get(peak).add(immPeaks.get(peak));
                if (immNPeaks.get(peak).get() >= (double) nPsms / 1.999) { // Slight round down for floating point error
                    double cMean = immNSumPeaks.get(peak).get() / (double) immNPeaks.get(peak).get();
                    immN.get(peak).add(cMean);
                } else {
                    immN.get(peak).add(0.0);
                }
            }
            for (Double peak : capYPeaks.keySet()) {
                capY.get(peak).add(capYPeaks.get(peak));
                if (capYNPeaks.get(peak).get() >= (double) nPsms / 1.999) { // Slight round down for floating point error
                    double cMean = capYNSumPeaks.get(peak).get() / (double) capYNPeaks.get(peak).get();
                    capYN.get(peak).add(cMean);
                } else {
                    capYN.get(peak).add(0.0);
                }
            }
            for (Character c : squigglePeaks.keySet()) {
                for (Double peak : squigglePeaks.get(c).keySet()) {
                    squiggles.get(c).get(peak).add(squigglePeaks.get(c).get(peak));
                    if (squiggleNPeaks.get(c).get(peak).get() >= (double) nPsms / 1.999) {
                        double cMean = squiggleNSumPeaks.get(c).get(peak).get() / (double) squiggleNPeaks.get(c).get(peak).get();
                        squigglesN.get(c).get(peak).add(cMean);
                    } else {
                        squigglesN.get(c).get(peak).add(0.0);
                    }
                }
            }
        }
    }

    private void collapseHashMaps() {
        collapseHashMap("control");
        collapseHashMap("treatment");
        collapseHashMap("decoy");
    }

    private void __collapseHashMaps() {
        /* Collapse pepkeys in control spectra */
        for (String pepKey : this.contPeptideMap.keySet()) {
            /* Init normalizing params and temp objects */
            int nPsms = this.contPeptideMap.get(pepKey).size();
            HashMap<Double, Double> immPeaks = new HashMap<>();
            HashMap<Double, Double> capYPeaks = new HashMap<>();
            HashMap<Character, HashMap<Double, Double>> squigglePeaks = new HashMap<>();
            for (Double peak : this.immoniumX.keySet())
                immPeaks.put(peak, 0.0);
            for (Double peak : this.capYX.keySet())
                capYPeaks.put(peak, 0.0);
            for (Character c : this.squigglesX.keySet()) {
                squigglePeaks.put(c, new HashMap<>());
                for (Double peak : this.squigglesX.get(c).keySet())
                    squigglePeaks.get(c).put(peak, 0.0);
            }
            /* Collapse into mean of respective vals */
            for (DiagnosticRecord dr : this.contPeptideMap.get(pepKey)) {
                for (Double peak : dr.selectedImmoniumPeaks.keySet()) {
                    double cVal = immPeaks.get(peak);
                    cVal += dr.selectedImmoniumPeaks.get(peak) / nPsms;
                    immPeaks.put(peak, cVal);
                   //if (dr.selectedImmoniumPeaks.get(peak) > 0.00000000001) { //count PSMs with ions
                    //    this.immoniumXN.put(peak, (this.immoniumXN.get(peak) + cVal));
                    //}
                }
                for (Double peak : dr.selectedCapYPeaks.keySet()) {
                    double cVal = capYPeaks.get(peak);
                    cVal += dr.selectedCapYPeaks.get(peak) / nPsms;
                    capYPeaks.put(peak, cVal);
                    //if (dr.selectedCapYPeaks.get(peak) > 0.00000000001) {//count PSMs with ions
                    //    this.capYXN.put(peak, (this.capYXN.get(peak) + cVal));
                    //}
                }
                for (Character c : dr.selectedSquigglePeaks.keySet()) {
                    for (Double peak : dr.selectedSquigglePeaks.get(c).keySet()) {
                        double cVal = squigglePeaks.get(c).get(peak);
                        cVal += dr.selectedSquigglePeaks.get(c).get(peak) / nPsms;
                        squigglePeaks.get(c).put(peak, cVal);
                        //if (dr.selectedSquigglePeaks.get(c).get(peak) > 0.00000000001) //count PSMs with ions
                        //    this.squigglesXN.get(c).put(peak, (this.squigglesXN.get(c).get(peak) + cVal));
                    }
                }
            }
            /* Add vals to class lists */
            for (Double peak : immPeaks.keySet())
                this.immoniumX.get(peak).add(immPeaks.get(peak));
            for (Double peak : capYPeaks.keySet())
                this.capYX.get(peak).add(capYPeaks.get(peak));
            for (Character c : squigglePeaks.keySet()) {
                for (Double peak : squigglePeaks.get(c).keySet()) {
                    this.squigglesX.get(c).get(peak).add(squigglePeaks.get(c).get(peak));
                }
            }
        }
        /* Collapse pepkeys in spec of interest */
        for (String pepKey : this.treatPeptideMap.keySet()) {
            /* Init normalizing params and temp objects */
            int nPsms = this.treatPeptideMap.get(pepKey).size();

            HashMap<Double, Double> immPeaks = new HashMap<>();
            HashMap<Double, Double> capYPeaks = new HashMap<>();
            HashMap<Character, HashMap<Double, Double>> squigglePeaks = new HashMap<>();
            for (Double peak : this.immoniumY.keySet())
                immPeaks.put(peak, 0.0);
            for (Double peak : this.capYY.keySet())
                capYPeaks.put(peak, 0.0);
            for (Character c : this.squigglesY.keySet()) {
                squigglePeaks.put(c, new HashMap<>());
                for (Double peak : this.squigglesY.get(c).keySet())
                    squigglePeaks.get(c).put(peak, 0.0);
            }
            /* Collapse into mean of respective vals */
            for (DiagnosticRecord dr : this.treatPeptideMap.get(pepKey)) {
                for (Double peak : dr.selectedImmoniumPeaks.keySet()) {
                    double cVal = immPeaks.get(peak);
                    cVal += dr.selectedImmoniumPeaks.get(peak) / nPsms;
                    immPeaks.put(peak, cVal);
                    //if (dr.selectedImmoniumPeaks.get(peak) > 1e-6) { //count PSMs with ions
                    //    this.immoniumYN.put(peak, (this.immoniumYN.get(peak) + cVal));
                    //}
                }
                for (Double peak : dr.selectedCapYPeaks.keySet()) {
                    double cVal = capYPeaks.get(peak);
                    cVal += dr.selectedCapYPeaks.get(peak) / nPsms;
                    capYPeaks.put(peak, cVal);
                    //if (dr.selectedCapYPeaks.get(peak) > 1e-6) //count PSMs with ions
                    //    this.capYYN.put(peak, (this.capYYN.get(peak) + cVal));
                }
                for (Character c : dr.selectedSquigglePeaks.keySet()) {
                    for (Double peak : dr.selectedSquigglePeaks.get(c).keySet()) {
                        double cVal = squigglePeaks.get(c).get(peak);
                        cVal += dr.selectedSquigglePeaks.get(c).get(peak) / nPsms;
                        squigglePeaks.get(c).put(peak, cVal);
                        //if (dr.selectedSquigglePeaks.get(c).get(peak) > 1e-6) //count PSMs with ions
                        //    this.squigglesYN.get(c).put(peak, (this.squigglesYN.get(c).get(peak) + cVal));
                    }
                }
            }

            /* Add vals to class lists */
            for (Double peak : immPeaks.keySet())
                this.immoniumY.get(peak).add(immPeaks.get(peak));
            for (Double peak : capYPeaks.keySet())
                this.capYY.get(peak).add(capYPeaks.get(peak));
            for (Character c : squigglePeaks.keySet()) {
                for (Double peak : squigglePeaks.get(c).keySet()) {
                    this.squigglesY.get(c).get(peak).add(squigglePeaks.get(c).get(peak));
                }
            }
        }
        /* Collapse pepkey spectra for decoys */
        for (String pepKey : this.decoyPeptideMap.keySet()) {
            /* Init normalizing params and temp objects */
            int nPsms = this.decoyPeptideMap.get(pepKey).size();
            HashMap<Double, Double> immPeaks = new HashMap<>();
            HashMap<Double, Double> capYPeaks = new HashMap<>();
            HashMap<Character, HashMap<Double, Double>> squigglePeaks = new HashMap<>();
            for (Double peak : this.immoniumDecoy.keySet())
                immPeaks.put(peak, 0.0);
            for (Double peak : this.capYDecoy.keySet())
                capYPeaks.put(peak, 0.0);
            for (Character c : this.squigglesDecoy.keySet()) {
                squigglePeaks.put(c, new HashMap<>());
                for (Double peak : this.squigglesDecoy.get(c).keySet())
                    squigglePeaks.get(c).put(peak, 0.0);
            }
            /* Collapse into mean of respective vals */
            for (DiagnosticRecord dr : this.decoyPeptideMap.get(pepKey)) {
                for (Double peak : dr.selectedImmoniumPeaks.keySet()) {
                    double cVal = immPeaks.get(peak);
                    cVal += dr.selectedImmoniumPeaks.get(peak) / nPsms;
                    immPeaks.put(peak, cVal);
                    //if (dr.selectedImmoniumPeaks.get(peak) > 1e-6) {//count PSMs with ions
                    //    this.immoniumDecoyN.put(peak, (this.immoniumDecoyN.get(peak) + cVal));
                    //}
                }
                for (Double peak : dr.selectedCapYPeaks.keySet()) {
                    double cVal = capYPeaks.get(peak);
                    cVal += dr.selectedCapYPeaks.get(peak) / nPsms;
                    capYPeaks.put(peak, cVal);
                    //if (dr.selectedCapYPeaks.get(peak) > 1e-6) //count PSMs with ions
                    //    this.capYDecoyN.put(peak, (this.capYDecoyN.get(peak) + cVal));
                }
                for (Character c : dr.selectedSquigglePeaks.keySet()) {
                    for (Double peak : dr.selectedSquigglePeaks.get(c).keySet()) {
                        double cVal = squigglePeaks.get(c).get(peak);
                        cVal += dr.selectedSquigglePeaks.get(c).get(peak) / nPsms;
                        squigglePeaks.get(c).put(peak, cVal);
                        //if (dr.selectedSquigglePeaks.get(c).get(peak) > 1e-6) //count PSMs with ions
                        //    this.squigglesDecoyN.get(c).put(peak, (this.squigglesDecoyN.get(c).get(peak) + cVal));
                    }
                }
            }

            /* Add vals to class lists */
            for (Double peak : immPeaks.keySet())
                this.immoniumDecoy.get(peak).add(immPeaks.get(peak));
            for (Double peak : capYPeaks.keySet())
                this.capYDecoy.get(peak).add(capYPeaks.get(peak));
            for (Character c : squigglePeaks.keySet()) {
                for (Double peak : squigglePeaks.get(c).keySet()) {
                    this.squigglesDecoy.get(c).get(peak).add(squigglePeaks.get(c).get(peak));
                }
            }
        }
    }

    public void _addVals(DiagnosticRecord dr, boolean isControl) {
        if (isControl) {
            for (Double peak : dr.selectedImmoniumPeaks.keySet())
                this.immoniumX.get(peak).add(dr.selectedImmoniumPeaks.get(peak));
            for (Double peak : dr.selectedCapYPeaks.keySet())
                this.capYX.get(peak).add(dr.selectedCapYPeaks.get(peak));
            for (Character c : dr.selectedSquigglePeaks.keySet()) {
                for (Double peak : dr.selectedSquigglePeaks.get(c).keySet())
                    this.squigglesX.get(c).get(peak).add(dr.selectedSquigglePeaks.get(c).get(peak));
            }
        } else {
            for (Double peak : dr.selectedImmoniumPeaks.keySet())
                this.immoniumY.get(peak).add(dr.selectedImmoniumPeaks.get(peak));
            for (Double peak : dr.selectedCapYPeaks.keySet())
                this.capYY.get(peak).add(dr.selectedCapYPeaks.get(peak));
            for (Character c : dr.selectedSquigglePeaks.keySet()) {
                for (Double peak : dr.selectedSquigglePeaks.get(c).keySet())
                    this.squigglesY.get(c).get(peak).add(dr.selectedSquigglePeaks.get(c).get(peak));
            }
        }
    }

    private double calcProportionWIon(ArrayList<Double> vals, int nPsms) {
        int n = 0;
        for (Double val : vals) {
            if (val > 0.0)
                n++;
        }
        double prop = (double) n / (double) nPsms;
        return prop;
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
        intensity /= (double) n;
        return intensity;
    }


    public void performTests() throws IOException {
        collapseHashMaps();
        MannWhitneyUTest mwu = new MannWhitneyUTest();
        TTest tTest = new TTest();

        this.immoniumTests = new ArrayList<>();
        for (Double peak : this.immoniumX.keySet()) {
            if (this.immoniumX.get(peak).size() < this.minPeps || this.immoniumY.get(peak).size() < this.minPeps)
                continue;
            double[] x = this.immoniumX.get(peak).stream().mapToDouble(i -> i).toArray();//.filter(i -> i > 0.0).toArray();
            double[] y = this.immoniumY.get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
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
            //double rankBiserCorr = Math.abs((2.0 * uStat / (this.immoniumX.get(peak).size() * this.immoniumY.get(peak).size())) - 1);
            long n1n2 = (long) this.immoniumX.get(peak).size() * (long) this.immoniumY.get(peak).size();
            //long n1n2 = (long)x.length * (long)y.length;
            double rankBiserCorr = u2 / (n1n2);
            //double rankBiserCorr = 1;
            double u1 = n1n2 - u2;
            boolean greaterThan = u2 > u1 ? true : false;
            if (this.twoTailedTests == false)
                p = convertP(p, greaterThan);
            //p *= this.immoniumX.size();
            double propWIon = calcProportionWIon(this.immoniumY.get(peak), this.nTreatPsms);
            double wIonIntensity = calcWIonIntensity(this.immoniumY.get(peak), this.nTreatPsms);
            double propWIonCont = calcProportionWIon(this.immoniumX.get(peak), this.nControlPsms);
            double wIonIntensityCont = calcWIonIntensity(this.immoniumX.get(peak), this.nControlPsms);
            //System.out.println(p + "\t" + rankBiserCorr);
            if (p <= this.maxP && Math.abs(rankBiserCorr - 0.5) >= this.minRbc)
                this.immoniumTests.add(new Test(peak, p, rankBiserCorr, false, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont,
                        u2, this.immoniumX.get(peak).size(), this.immoniumY.get(peak).size()));
            //System.out.printf("Immonium %.04f\t%e\t%f\n", peak, p, rankBiserCorr);
        }
        for (Double peak : this.immoniumX.keySet()) {
            if (this.immoniumX.get(peak).size() < this.minPeps || this.immoniumDecoyN.get(peak).size() < this.minPeps)
                continue;
            //double[] x = this.immoniumX.get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
            //double[] y = this.immoniumDecoy.get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
            double[] x = this.immoniumX.get(peak).stream().mapToDouble(i -> i).toArray();//.filter(i -> i > 0.0).toArray();
            double[] y = this.immoniumDecoy.get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
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
            //double rankBiserCorr = Math.abs((2.0 * uStat / (this.immoniumX.get(peak).size() * this.immoniumY.get(peak).size())) - 1);
            long n1n2 = (long)x.length * (long)y.length;
            double rankBiserCorr = u2 / (n1n2);
            //double rankBiserCorr = 1;
            double u1 = n1n2 - u2;
            boolean greaterThan = u2 > u1 ? true : false;
            if (this.twoTailedTests == false)
                p = convertP(p, greaterThan);
            //p *= this.immoniumX.size();
            double propWIon = calcProportionWIon(this.immoniumDecoy.get(peak), this.nTreatPsms);
            double wIonIntensity = calcWIonIntensity(this.immoniumDecoy.get(peak), this.nTreatPsms);
            double propWIonCont = calcProportionWIon(this.immoniumX.get(peak), this.nControlPsms);
            double wIonIntensityCont = calcWIonIntensity(this.immoniumX.get(peak), this.nControlPsms);
            if (p <= this.maxP && Math.abs(rankBiserCorr - 0.5) >= this.minRbc)
                this.immoniumTests.add(new Test(peak, p, rankBiserCorr, true,propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                        this.immoniumX.get(peak).size(), this.immoniumDecoy.get(peak).size()));
            //System.out.printf("Immonium %.04f\t%e\t%f\n", peak, p, rankBiserCorr);
        }

        this.capYTests = new ArrayList<>();
        for (Double peak : this.capYX.keySet()) {
            if (this.capYX.get(peak).size() < this.minPeps || this.capYYN.get(peak).size() < this.minPeps)
                continue;
            double[] x = this.capYX.get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
            double[] y = this.capYY.get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
            //double[] x = this.capYX.get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
            //double[] y = this.capYY.get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
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
            //double rankBiserCorr = Math.abs((2.0 * uStat / (this.capYX.get(peak).size() * this.capYY.get(peak).size())) - 1);
            //p *= this.capYX.size();
            double propWIon = calcProportionWIon(this.capYY.get(peak), this.nTreatPsms);
            double wIonIntensity = calcWIonIntensity(this.capYY.get(peak), this.nTreatPsms);
            double propWIonCont = calcProportionWIon(this.capYX.get(peak), this.nControlPsms);
            double wIonIntensityCont = calcWIonIntensity(this.capYX.get(peak), this.nControlPsms);
            if (p <= this.maxP && Math.abs(rankBiserCorr - 0.5) >= this.minRbc)
                this.capYTests.add(new Test(peak, p, rankBiserCorr, false, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                        this.capYX.get(peak).size(), this.capYY.get(peak).size()));
            //System.out.printf("CapY %.04f\t%e\t%f\n", peak, p, rankBiserCorr);
        }
        for (Double peak : this.capYX.keySet()) {
            if (this.capYX.get(peak).size() < this.minPeps || this.capYDecoyN.get(peak).size() < this.minPeps)
                continue;
            //double[] x = this.capYX.get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
            //double[] y = this.capYDecoy.get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
            double[] x = this.capYX.get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
            double[] y = this.capYDecoy.get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
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
            //double rankBiserCorr = Math.abs((2.0 * uStat / (this.capYX.get(peak).size() * this.capYY.get(peak).size())) - 1);
            //p *= this.capYX.size();
            double propWIon = calcProportionWIon(this.capYDecoy.get(peak), this.nTreatPsms);
            double wIonIntensity = calcWIonIntensity(this.capYDecoy.get(peak), this.nTreatPsms);
            double propWIonCont = calcProportionWIon(this.capYX.get(peak), this.nControlPsms);
            double wIonIntensityCont = calcWIonIntensity(this.capYX.get(peak), this.nControlPsms);
            if (p <= this.maxP && Math.abs(rankBiserCorr - 0.5) >= this.minRbc)
                this.capYTests.add(new Test(peak, p, rankBiserCorr, true,propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                        this.capYX.get(peak).size(), this.capYDecoy.get(peak).size()));
            //System.out.printf("CapY %.04f\t%e\t%f\n", peak, p, rankBiserCorr);
        }

        this.squigglesTests = new HashMap<>();
        for (Character c : this.squigglesX.keySet()) {
            this.squigglesTests.put(c, new ArrayList<>());
            //String fname = this.peakApex + "_" + c + "sig_cles.tsv";
            //PrintWriter out = new PrintWriter(new FileWriter(new File(fname)));
            //out.printf("mass\tp\trbc\n");
            //System.out.println("Squiggle"+c);
            for (Double peak : this.squigglesX.get(c).keySet()) {
                if (this.squigglesX.get(c).get(peak).size() < this.minPeps || this.squigglesY.get(c).get(peak).size() < this.minPeps)
                    continue;
                //double[] x = this.squigglesX.get(c).get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
                //double[] y = this.squigglesY.get(c).get(peak).stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
                double[] x = this.squigglesX.get(c).get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
                double[] y = this.squigglesY.get(c).get(peak).stream().mapToDouble(i -> i).toArray();//filter(i -> i > 0.0).toArray();
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
                //double rankBiserCorr = Math.abs((2.0 * uStat / (this.squigglesX.get(c).get(peak).size() * this.squigglesY.get(c).get(peak).size())) - 1);
                //out.printf("%.04f\t%e\t%f\n", peak, p, rankBiserCorr);
                //p *= this.squigglesX.get(c).size();
                double propWIon = calcProportionWIon(this.squigglesY.get(c).get(peak), this.nTreatPsms);
                double wIonIntensity = calcWIonIntensity(this.squigglesY.get(c).get(peak), this.nTreatPsms);
                double propWIonCont = calcProportionWIon(this.squigglesX.get(c).get(peak), this.nControlPsms);
                double wIonIntensityCont = calcWIonIntensity(this.squigglesX.get(c).get(peak), this.nControlPsms);
                //if (p <= this.maxP && Math.abs(rankBiserCorr - 0.5) >= this.minRbc) {
                if (p <= this.maxP && Math.abs(rankBiserCorr - 0.5) >= this.minRbc) {
                    this.squigglesTests.get(c).add(new Test(peak, p, rankBiserCorr, false, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                            this.squigglesX.get(c).get(peak).size(), this.squigglesY.get(c).get(peak).size()));
                }
                //if (p * this.squigglesX.get(c).get(peak).size() < 0.05 && rankBiserCorr > 0.5)
                //System.out.printf("%.04f\t%e\t%f\n", peak, p, rankBiserCorr);
            }
            //out.close();
        }
        for (Character c : this.squigglesX.keySet()) {
            //this.squigglesTests.put(c, new ArrayList<>());
            //String fname = this.peakApex + "_" + c + "sig_cles.tsv";
            //PrintWriter out = new PrintWriter(new FileWriter(new File(fname)));
            //out.printf("mass\tp\trbc\n");
            //System.out.println("Squiggle"+c);
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
                //double rankBiserCorr = Math.abs((2.0 * uStat / (this.squigglesX.get(c).get(peak).size() * this.squigglesY.get(c).get(peak).size())) - 1);
                //out.printf("%.04f\t%e\t%f\n", peak, p, rankBiserCorr);
                //p *= this.squigglesX.get(c).size();
                double propWIon = calcProportionWIon(this.squigglesDecoy.get(c).get(peak), this.nTreatPsms);
                double wIonIntensity = calcWIonIntensity(this.squigglesDecoy.get(c).get(peak), this.nTreatPsms);
                double propWIonCont = calcProportionWIon(this.squigglesX.get(c).get(peak), this.nControlPsms);
                double wIonIntensityCont = calcWIonIntensity(this.squigglesX.get(c).get(peak), this.nControlPsms);
                if (p <= this.maxP && Math.abs(rankBiserCorr - 0.5) >= this.minRbc) {
                    this.squigglesTests.get(c).add(new Test(peak, p, rankBiserCorr, true, propWIon, propWIonCont, wIonIntensity, wIonIntensityCont, u2,
                            this.squigglesX.get(c).get(peak).size(), this.squigglesDecoy.get(c).get(peak).size()));
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
        collapseTests(this.capYTests, 0.01, false); //todo tol
        //System.out.println("squiggle");
        for (Character c : this.squigglesTests.keySet())
            collapseTests(this.squigglesTests.get(c), 0.01, true); //todo tol


        /* clear out most of memory */
        clearMemory();

        //Print all tests
        boolean debug = false;
        if (debug) {
            Collections.sort(this.immoniumTests);
            for (Test t : this.immoniumTests)
                System.out.printf("Immonium %.04f\t%e\t%f\t%.04f\t%b\t%d\t%d\n", t.mass, t.q, t.rbc, t.u, t.isDecoy, t.n1, t.n2);
            Collections.sort(this.capYTests);
            for (Test t : this.capYTests)
                System.out.printf("CapY %.04f\t%e\t%f\t%.04f\t%b\t%d\t%d\n", t.mass, t.q, t.rbc, t.u, t.isDecoy, t.n1, t.n2);
            for (Character c : this.squigglesTests.keySet()) {
                Collections.sort(this.squigglesTests.get(c));
                for (Test t : this.squigglesTests.get(c))
                    System.out.printf("Squiggle %.04f\t%e\t%f\t%.04f\t%b\t%d\t%d\n", t.mass, t.q, t.rbc, t.u, t.isDecoy, t.n1, t.n2);
            }
        }
    }

    /* Converts two-tailed p-value to one-tailed p-value
    * Hipparchus only provides a two-tailed version of MWU
    * */
    private double convertP(double oldPVal, boolean greaterThan) {
        double p;
        /*
        if (greaterThan)
            p = oldPVal / 2.0;
        else
            p = 1.0 - oldPVal / 2.0;
         */
        return oldPVal;
    }

    private void clearMemory() {
        this.contPeptideMap = null;
        this.treatPeptideMap = null;
        this.decoyPeptideMap = null;
        this.immoniumX = null;
        this.immoniumXN = null;
        this.immoniumY = null;
        this.immoniumYN = null;
        this.immoniumDecoy = null;
        this.immoniumDecoyN = null;
        this.capYX = null;
        this.capYXN = null;
        this.capYY = null;
        this.capYYN = null;
        this.capYDecoy = null;
        this.capYDecoyN = null;
        this.squigglesX = null;
        this.squigglesY = null;
        this.squigglesDecoy = null;
        this.squigglesXN = null;
        this.squigglesYN = null;
        this.squigglesDecoyN = null;
    }

    private void collapseTests(ArrayList<Test> tests, double tol, boolean adjustedMass) {
        /* Collapse immonium tests */
        int group = 1;
        Collections.sort(tests, new Comparator<Test>() {
            @Override
            public int compare(Test o1, Test o2) {
                return -1*Double.compare(o1.rbc, o2.rbc);
            }
        });
        for (int i = 0; i < tests.size(); i++) {
            if (tests.get(i).group > 0)
                continue;
            tests.get(i).group = group;
            //System.out.println(tests.get(i).mass);
            for (int j = 0; j < tests.size(); j++) {
                if (i == j)
                    continue;
                if (tests.get(j).group > 0)
                    continue;
                for (int c = -1; c < 4; c++) {
                    if (c == 0)
                        continue;
                    double c13 = 1.003355;
                    double shiftMass = c13 * c;
                    if (Math.abs(tests.get(i).mass - (tests.get(j).mass - shiftMass)) < tol) {
                        tests.get(j).group = group;
                        //System.out.println("* " + tests.get(j).mass);
                    }
                }
            }
            group++;
        }
    }

    /* Checks for enrichment of adjacent AAs for each squiggle peak, then adjusts mass appropriately */
    private void relocalizeDeltas(int nAdjacentAas, double tol) {
        //System.out.println("Relocalizing deltas");
        for (Character c : this.squigglesTests.keySet()) {
            //System.out.println("**" + c);
            int[][] shiftedCounts = new int[nAdjacentAas * 2 + 1][26];
            for (Test t : this.squigglesTests.get(c)) {
                //System.out.println(t.mass);
                int[][] resCounts = new int[nAdjacentAas * 2 + 1][26]; //todo tmp
                int[][] bgResCounts = new int[nAdjacentAas * 2 + 1][26];
                int[] totalChances = new int[nAdjacentAas*2+1];
                int totalMatchedIons = 0;
                for (String pepKey : this.treatPeptideMap.keySet()) {
                    //todo currently collapsing on pepKeys
                    ArrayList<int[]> ionCounts = new ArrayList<>();
                    String pepSeq = "";
                    for (DiagnosticRecord dr : this.treatPeptideMap.get(pepKey)) {
                        /* Get and normalize scores */
                        int[] scores = dr.localizeRemainderMass((float) t.mass, c);
                        pepSeq = dr.pepSeq;
                        ionCounts.add(scores);
                    }
                    /* Add pepKey's modal ion counts to total */
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

                            if (resP > 0.75 && resP > bestResP) {
                                newRemainderMass = calcNewRemainderMass(t.mass, c, i-nAdjacentAas, j);
                                //System.out.println(t.mass + "\t" + newRemainderMass + "\t" + c + "\t" + (i-nAdjacentAas) + "\t" + j);
                                bestResP = resP;
                                shiftedCounts[i][j]++;
                            }
                        }
                    }
                }
                //System.out.println(newRemainderMass + "\t" + t.mass);
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

    private void performTests(ArrayList<Double> xNums, ArrayList<Double> yNums) {
        /* Get nonzero values from array */
        double[] x = xNums.stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
        double[] y = yNums.stream().mapToDouble(i -> i).filter(i -> i > 0.0).toArray();
        if (x.length < this.minPeps || y.length < this.minPeps) {

        }
    }

}

class Test implements Comparable<Test> {
    public double mass;
    public double adjustedMass;
    public double q;
    public double rbc;
    public boolean isDecoy;
    public int group;
    public double u;
    public long n1;
    public long n2;
    public double propWIonTreat;
    public double propWIonCont;
    public double propWIonIntensity;
    public double propWIonIntensityCont;

    Test(double mass, double q, double rbc, boolean isDecoy, double propWIonTreat, double propWIonCont,
         double propWIonIntensity, double propWIonIntensityCont, double u, long n1, long n2) {
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
        this.n1 = n1;
        this.n2 = n2;
    }

    public String toString() {
        return String.format("%.04f\t%e\t%f\t%.04f\t%b\t%d\t%d", this.mass, this.q, this.rbc, this.u, this.isDecoy, this.n1, this.n2);
    }

    @Override
    public int compareTo(Test arg0) {
        return Double.compare(this.adjustedMass, arg0.adjustedMass);
    }

}
