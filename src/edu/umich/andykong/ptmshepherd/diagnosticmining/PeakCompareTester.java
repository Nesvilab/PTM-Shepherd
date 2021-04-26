package edu.umich.andykong.ptmshepherd.diagnosticmining;

import com.google.gson.internal.$Gson$Types;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import umich.ms.datatypes.lcmsrun.Hash;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;

public class PeakCompareTester {
    double peakApex;
    HashMap<Double, ArrayList<Double>> immoniumX;
    HashMap<Double,ArrayList<Double>> immoniumY;
    HashMap<Double,ArrayList<Double>> capYX;
    HashMap<Double,ArrayList<Double>> capYY;
    HashMap<Character, HashMap<Double,ArrayList<Double>>> squigglesX;
    HashMap<Character, HashMap<Double,ArrayList<Double>>> squigglesY;
    HashMap<String, ArrayList<DiagnosticRecord>> contPeptideMap;
    HashMap<String, ArrayList<DiagnosticRecord>> treatPeptideMap;
    double [] x;
    double [] y;
    ArrayList<Double> xNums;
    ArrayList<Double> yNums;

    public ArrayList<Test> immoniumTests;
    public ArrayList<Test> capYTests;
    public HashMap<Character, ArrayList<Test>> squigglesTests;

    public PeakCompareTester(ArrayList<Double> xVals, ArrayList<Double> yVals) {
        this.x = xVals.stream().mapToDouble(i -> i).toArray();
        this.y = yVals.stream().mapToDouble(i -> i).toArray();
    }

    public PeakCompareTester(double peakApex, ArrayList<Double> unifImm, ArrayList<Double> unifCapY, HashMap<Character, ArrayList<Double>> unifSquig) {
        this.peakApex = peakApex;
        this.immoniumX = new HashMap<>();
        this.immoniumY = new HashMap<>();
        for (Double peak : unifImm) {
            this.immoniumX.put(peak, new ArrayList<>());
            this.immoniumY.put(peak, new ArrayList<>());
        }
        this.capYX = new HashMap<>();
        this.capYY = new HashMap<>();
        for (Double peak : unifCapY) {
            this.capYX.put(peak, new ArrayList<>());
            this.capYY.put(peak, new ArrayList<>());
        }
        this.squigglesX = new HashMap<>();
        this.squigglesY = new HashMap<>();
        for (Character ion : unifSquig.keySet()) {
            this.squigglesX.put(ion, new HashMap<>());
            this.squigglesY.put(ion, new HashMap<>());
            for (Double peak : unifSquig.get(ion)) {
                this.squigglesX.get(ion).put(peak, new ArrayList<>());
                this.squigglesY.get(ion).put(peak, new ArrayList<>());
            }
        }
        this.contPeptideMap = new HashMap<>();
        this.treatPeptideMap = new HashMap<>();
    }

    public void addVals(DiagnosticRecord dr, boolean isControl) {
        String pepKey = dr.pepSeq + dr.modifications + dr.charge; //peptidekey includes charge ///todo is that ideal?
        if (isControl) {
            if (!this.contPeptideMap.containsKey(pepKey))
                this.contPeptideMap.put(pepKey, new ArrayList<>());
            this.contPeptideMap.get(pepKey).add(dr);
        } else {
            if (!this.treatPeptideMap.containsKey(pepKey))
                this.treatPeptideMap.put(pepKey, new ArrayList<>());
            this.treatPeptideMap.get(pepKey).add(dr);
        }
    }

    private void collapseHashMaps() {
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
                }
                for (Double peak : dr.selectedCapYPeaks.keySet()) {
                    double cVal = capYPeaks.get(peak);
                    cVal += dr.selectedCapYPeaks.get(peak) / nPsms;
                    capYPeaks.put(peak, cVal);
                }
                for (Character c : dr.selectedSquigglePeaks.keySet()) {
                    for (Double peak : dr.selectedSquigglePeaks.get(c).keySet()) {
                        double cVal = squigglePeaks.get(c).get(peak);
                        cVal += dr.selectedSquigglePeaks.get(c).get(peak) / nPsms;
                        squigglePeaks.get(c).put(peak, cVal);
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
                }
                for (Double peak : dr.selectedCapYPeaks.keySet()) {
                    double cVal = capYPeaks.get(peak);
                    cVal += dr.selectedCapYPeaks.get(peak) / nPsms;
                    capYPeaks.put(peak, cVal);
                }
                for (Character c : dr.selectedSquigglePeaks.keySet()) {
                    for (Double peak : dr.selectedSquigglePeaks.get(c).keySet()) {
                        double cVal = squigglePeaks.get(c).get(peak);
                        cVal += dr.selectedSquigglePeaks.get(c).get(peak) / nPsms;
                        squigglePeaks.get(c).put(peak, cVal);
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

    public void performTests() throws IOException {
        collapseHashMaps();
        MannWhitneyUTest mwu = new MannWhitneyUTest();

        this.immoniumTests = new ArrayList<>();
        for (Double peak : this.immoniumX.keySet()) {
            double p = mwu.mannWhitneyUTest(this.immoniumX.get(peak).stream().mapToDouble(i -> i).toArray(),
                    this.immoniumY.get(peak).stream().mapToDouble(i -> i).toArray());
            double uStat = mwu.mannWhitneyU(this.immoniumX.get(peak).stream().mapToDouble(i -> i).toArray(),
                    this.immoniumY.get(peak).stream().mapToDouble(i -> i).toArray());
            double rankBiserCorr = (2.0 * uStat / (this.immoniumX.get(peak).size() * this.immoniumY.get(peak).size())) - 1;
            //if (p * this.immoniumY.get(peak).size() < 0.05 && rankBiserCorr > 0.5)
            //    System.out.printf("Immonium %.04f\t%e\t%f\n", peak, p, rankBiserCorr);
            this.immoniumTests.add(new Test(peak, p, rankBiserCorr));
        }

        this.capYTests = new ArrayList<>();
        for (Double peak : this.capYX.keySet()) {
            double p = mwu.mannWhitneyUTest(this.capYX.get(peak).stream().mapToDouble(i -> i).toArray(),
                    this.capYY.get(peak).stream().mapToDouble(i -> i).toArray());
            double uStat = mwu.mannWhitneyU(this.capYX.get(peak).stream().mapToDouble(i -> i).toArray(),
                    this.capYY.get(peak).stream().mapToDouble(i -> i).toArray());
            double rankBiserCorr = (2.0 * uStat / (this.capYX.get(peak).size() * this.capYY.get(peak).size())) - 1;
            if (p * this.capYY.get(peak).size() < 0.05 && rankBiserCorr > 0.5)
                System.out.printf("CapY %.04f\t%e\t%f\n", peak, p, rankBiserCorr);
            this.capYTests.add(new Test(peak, p, rankBiserCorr));
        }

        this.squigglesTests = new HashMap<>();
        for (Character c : this.squigglesX.keySet()) {
            this.squigglesTests.put(c, new ArrayList<>());
            //String fname = this.peakApex + "_" + c + "sig_cles.tsv";
            //PrintWriter out = new PrintWriter(new FileWriter(new File(fname)));
            //out.printf("mass\tp\trbc\n");
            //System.out.println("Squiggle"+c);
            for (Double peak : this.squigglesX.get(c).keySet()) {
                double p = mwu.mannWhitneyUTest(this.squigglesX.get(c).get(peak).stream().mapToDouble(i -> i).toArray(),
                        this.squigglesY.get(c).get(peak).stream().mapToDouble(i -> i).toArray());
                double uStat = mwu.mannWhitneyU(this.squigglesX.get(c).get(peak).stream().mapToDouble(i -> i).toArray(),
                        this.squigglesY.get(c).get(peak).stream().mapToDouble(i -> i).toArray());
                double rankBiserCorr = (2.0 * uStat / (this.squigglesX.get(c).get(peak).size() * this.squigglesY.get(c).get(peak).size())) - 1;
                //out.printf("%.04f\t%e\t%f\n", peak, p, rankBiserCorr);
                if (p * this.squigglesX.get(c).get(peak).size() < 0.05 && rankBiserCorr > 0.5)
                    this.squigglesTests.get(c).add(new Test(peak, p * this.squigglesX.get(c).get(peak).size(), rankBiserCorr));
                //if (p * this.squigglesX.get(c).get(peak).size() < 0.05 && rankBiserCorr > 0.5)
                //    System.out.printf("%.04f\t%e\t%f\n", peak, p, rankBiserCorr);
            }
            //out.close();
        }
        relocalizeDeltas(0.05, 0.5, 1, 20);
        this.contPeptideMap = null;
        this.treatPeptideMap = null;
    }

    private void relocalizeDeltas(double minP, double minRoc, int nAdjacentAas, double tol) {
        for (Character c : this.squigglesTests.keySet()) {
            int[][] shiftedCounts = new int[nAdjacentAas * 2 + 1][26];
            for (Test t : this.squigglesTests.get(c)) {
                //ArrayList<double[][]> allResCounts = new ArrayList<>();
                if (t.q < minP) {
                    int[][] resCounts = new int[nAdjacentAas * 2 + 1][26]; //todo tmp
                    int[][] bgResCounts = new int[nAdjacentAas * 2 + 1][26];
                    int[] totalChances = new int[nAdjacentAas*2+1];
                    int totalMatchedIons = 0;
                    //System.out.println("***************"+t.mass + "\t"+c);
                    for (String pepKey : this.treatPeptideMap.keySet()) {
                        //double[][] resCounts = new double[nAdjacentAas * 2 + 1][26];
                        double nPsms = this.treatPeptideMap.keySet().size();
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
                        //System.out.println(pepSeq);
                        //for (int i = 0; i < ionModes.length; i++)
                        //    System.out.print(ionModes[i]);
                        //System.out.println();
                        //allResCounts.add(resCounts); todo reinstall, tmp comment
                    }
                    double newRemainderMass = t.mass;
                    double bestResP = -1;
                    /* Merge leucine and isoleucine onto leucine*/
                    for (int i = 0; i < resCounts.length; i++) {
                        resCounts[i]['l' - 'a'] += resCounts[i]['i' - 'a'];
                        resCounts[i]['i' - 'a'] = 0;
                    }
                    for (int i = 0; i < resCounts.length; i++) {
                        //double maxOddsR = -1;
                        //int maxIndx = -1;
                        //if (i - nAdjacentAas == 0) // Skip relative position of 0
                        //    continue;
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

                                if (resP > 0.75 && resP > bestResP) {
                                    newRemainderMass = calcNewRemainderMass(t.mass, c, i-nAdjacentAas, j);
                                    bestResP = resP;
                                    shiftedCounts[i][j]++;
                                }

                                /*
                                double bgP = (double) bgResCounts[i][j] / (double) totalChances[i];
                                double oddsRatio = resP / bgP;
                                */
                                //System.out.println(cRes + "\t" + i + "\t" + cResLoc + "\t" + totalMatchedIons + "\t" + resP);
                                /*
                                System.out.println(cRes + "\t" + cResBg + "\t" + totalChances[i] + "\t" + bgP +  "\t" + 0);
                                System.out.println(i+"\t"+j+"\t"+(double)resCounts[i][j]/(double)bgResCounts[i][j]);
                                FisherExact fe = new FisherExact(cResLoc + notCResLoc + cResNotLoc + notCResNotLoc);
                                double p = fe.getRightTailedP(cResLoc , notCResLoc, cResNotLoc, notCResNotLoc);
                                System.out.println(cRes + "\t" + i + "\t" + j + "\t" + oddsRatio + "\t" + p);
                                */
                            }
                        }
                    }
                    //System.out.printf("Shifted apexes: %f\t%c\t%s\t%f\t%f%n", this.peakApex, c, t.mass, newRemainderMass, t.rbc);
                    t.adjustedMass = newRemainderMass;
                }
            }
            //for(int i = 0; i < shiftedCounts.length; i++) {
            //    System.out.println(i);
            //    for(int j = 0; j < shiftedCounts[i].length; j++)
            //        System.out.print((char)(j+'a')+" "+shiftedCounts[i][j]+" ");
             //   System.out.println();
            //}
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
        double newDmass = 100000000;
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
        return newDmass;
    }
}

class Test {
    public double mass;
    public double adjustedMass;
    public double q;
    public double rbc;

    Test(double mass, double q, double rbc) {
        this.mass = mass;
        this.q = q;
        this.rbc = rbc;
    }
}
