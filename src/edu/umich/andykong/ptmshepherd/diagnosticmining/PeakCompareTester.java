package edu.umich.andykong.ptmshepherd.diagnosticmining;

import com.google.gson.internal.$Gson$Types;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
//import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import umich.ms.datatypes.lcmsrun.Hash;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.*;
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

    double maxP;
    double minRbc;
    boolean twoTailedTests;

    double specTol;

    public PeakCompareTester(ArrayList<Double> xVals, ArrayList<Double> yVals) {
        this.x = xVals.stream().mapToDouble(i -> i).toArray();
        this.y = yVals.stream().mapToDouble(i -> i).toArray();
    }

    public PeakCompareTester(double peakApex, ArrayList<Double> unifImm, ArrayList<Double> unifCapY, HashMap<Character, ArrayList<Double>> unifSquig, double maxP, double minRbc, boolean twoTailedTests, double specTol) {
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

        this.maxP = maxP;
        this.minRbc = minRbc;
        this.twoTailedTests = twoTailedTests;

        this.specTol = specTol;
    }

    public synchronized void addDrs(ArrayList<DiagnosticRecord> drs, boolean isControl) {
        for (DiagnosticRecord dr : drs)
            addVals(dr, isControl);
    }

    public void addVals(DiagnosticRecord dr, boolean isControl) {
        String pepKey = dr.pepSeq + dr.modifications + dr.charge; //peptidekey includes charge
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
            double u2 = mwu.mannWhitneyU2(this.immoniumX.get(peak).stream().mapToDouble(i -> i).toArray(),
                    this.immoniumY.get(peak).stream().mapToDouble(i -> i).toArray());
            //double rankBiserCorr = Math.abs((2.0 * uStat / (this.immoniumX.get(peak).size() * this.immoniumY.get(peak).size())) - 1);
            long n1n2 = (long) this.immoniumX.get(peak).size() * (long) this.immoniumY.get(peak).size();
            double rankBiserCorr = u2 / (n1n2);
            double u1 = n1n2 - u2;
            boolean greaterThan = u2 > u1 ? true : false;
            if (this.twoTailedTests == false)
                p = convertP(p, greaterThan);
            if (p < this.maxP / this.immoniumY.get(peak).size() && rankBiserCorr > this.minRbc)
                this.immoniumTests.add(new Test(peak, p, rankBiserCorr, u2, this.immoniumX.get(peak).size(), this.immoniumY.get(peak).size()));
            //System.out.printf("Immonium %.04f\t%e\t%f\n", peak, p, rankBiserCorr);
        }

        this.capYTests = new ArrayList<>();
        for (Double peak : this.capYX.keySet()) {
            double p = mwu.mannWhitneyUTest(this.capYX.get(peak).stream().mapToDouble(i -> i).toArray(),
                    this.capYY.get(peak).stream().mapToDouble(i -> i).toArray());
            double u2 = mwu.mannWhitneyU2(this.capYX.get(peak).stream().mapToDouble(i -> i).toArray(),
                    this.capYY.get(peak).stream().mapToDouble(i -> i).toArray());
            long n1n2 = (long) this.capYX.get(peak).size() * (long) this.capYY.get(peak).size();
            double rankBiserCorr = u2 / n1n2;
            double u1 = n1n2 - u2;
            boolean greaterThan = u2 > u1 ? true : false;
            if (this.twoTailedTests == false)
                p = convertP(p, greaterThan);
            //double rankBiserCorr = Math.abs((2.0 * uStat / (this.capYX.get(peak).size() * this.capYY.get(peak).size())) - 1);
            if (p < this.maxP / this.capYY.get(peak).size() && rankBiserCorr > this.minRbc)
                this.capYTests.add(new Test(peak, p, rankBiserCorr, u2, this.capYX.get(peak).size(), this.capYY.get(peak).size()));
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
                double p = mwu.mannWhitneyUTest(this.squigglesX.get(c).get(peak).stream().mapToDouble(i -> i).toArray(),
                        this.squigglesY.get(c).get(peak).stream().mapToDouble(i -> i).toArray());
                double u2 = mwu.mannWhitneyU2(this.squigglesX.get(c).get(peak).stream().mapToDouble(i -> i).toArray(),
                        this.squigglesY.get(c).get(peak).stream().mapToDouble(i -> i).toArray());
                long n1n2 = (long) this.squigglesX.get(c).get(peak).size() * (long) this.squigglesY.get(c).get(peak).size();
                double rankBiserCorr = 2 * u2 / (n1n2) - 1;
                double u1 = n1n2 - u2;
                boolean greaterThan = u2 > u1 ? true : false;
                if (this.twoTailedTests == false)
                    p = convertP(p, greaterThan);
                //double rankBiserCorr = Math.abs((2.0 * uStat / (this.squigglesX.get(c).get(peak).size() * this.squigglesY.get(c).get(peak).size())) - 1);
                //out.printf("%.04f\t%e\t%f\n", peak, p, rankBiserCorr);
                if (p < this.maxP / this.squigglesX.get(c).get(peak).size() && rankBiserCorr > this.minRbc)
                    this.squigglesTests.get(c).add(new Test(peak, p, rankBiserCorr, u2, this.squigglesX.get(c).get(peak).size(), this.squigglesY.get(c).get(peak).size()));
                //if (p * this.squigglesX.get(c).get(peak).size() < 0.05 && rankBiserCorr > 0.5)
                //System.out.printf("%.04f\t%e\t%f\n", peak, p, rankBiserCorr);
            }
            //out.close();
        }

        relocalizeDeltas(1, this.specTol); //todo tol and nAdjacentAAs??
        System.out.println("immonium");
        collapseTests(this.immoniumTests, 0.01, false); //todo tol
        System.out.println("capY");
        collapseTests(this.capYTests, 0.01, false); //todo tol
        System.out.println("squiggle");
        for (Character c : this.squigglesTests.keySet())
            collapseTests(this.squigglesTests.get(c), 0.01, true); //todo tol


        this.contPeptideMap = null;
        this.treatPeptideMap = null;

        //Print all tests
        Collections.sort(this.immoniumTests);
        for (Test t : this.immoniumTests)
            System.out.printf("Immonium %.04f\t%e\t%f\t%.04f\t%d\t%d\n", t.mass, t.q, t.rbc, t.u, t.n1, t.n2);
        Collections.sort(this.capYTests);
        for (Test t : this.capYTests)
            System.out.printf("CapY %.04f\t%e\t%f\t%.04f\t%d\t%d\n", t.mass, t.q, t.rbc, t.u, t.n1, t.n2);
        for (Character c : this.squigglesTests.keySet()) {
            Collections.sort(this.squigglesTests.get(c));
            for (Test t : this.squigglesTests.get(c))
                System.out.printf("Squiggle %.04f\t%e\t%f\t%.04f\t%d\t%d\n", t.mass, t.q, t.rbc, t.u, t.n1, t.n2);
        }
    }

    /* Converts two tailed p-value to one tailed p-value
    * Hipparchus only provides a two-tailed version of MWU
    * */
    private double convertP(double oldPVal, boolean greaterThan) {
        double p;
        if (greaterThan)
            p = oldPVal / 2;
        else
            p = 1- oldPVal / 2;
        return p;
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
            System.out.println(tests.get(i).mass);
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
                        System.out.println("* " + tests.get(j).mass);
                    }
                }
            }
            group++;
        }
    }

    /* Checks for enrichment of adjacent AAs for each squiggle peak, then adjusts mass appropriately */
    private void relocalizeDeltas(int nAdjacentAas, double tol) {
        for (Character c : this.squigglesTests.keySet()) {
            int[][] shiftedCounts = new int[nAdjacentAas * 2 + 1][26];
            for (Test t : this.squigglesTests.get(c)) {
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

                            if (resP > 0.75 && resP > bestResP) {
                                newRemainderMass = calcNewRemainderMass(t.mass, c, i-nAdjacentAas, j);
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

class Test implements Comparable<Test> {
    public double mass;
    public double adjustedMass;
    public double q;
    public double rbc;
    public int group;
    public double u;
    public long n1;
    public long n2;

    Test(double mass, double q, double rbc, double u, long n1, long n2) {
        this.mass = mass;
        this.q = q;
        this.rbc = rbc;
        this.u = u;
        this.group = 0;
        this.n1 = n1;
        this.n2 = n2;
    }

    @Override
    public int compareTo(Test arg0) {
        return Double.compare(this.mass, arg0.mass);
    }

}
