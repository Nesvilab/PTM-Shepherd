package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import umich.ms.datatypes.lcmsrun.Hash;

import java.lang.reflect.Array;
import java.nio.CharBuffer;
import java.util.*;

/* This class holds the transformed peaklists of a single PSM
*  Each instance of this will be entered into a binary diagnostic output file */
public class DiagnosticRecord implements Comparable<DiagnosticRecord>  {
    public int scanNum;
    public float[][] immoniumPeaks;
    public float[][] capYPeaks;
    public HashMap<Character, float[][]> squigglePeaks;
    public String ionTypes;
    public String pepSeq;
    public int charge;
    public float[] modificationsArray;
    public HashMap<Integer, Float> modifications;
    public float dmass;
    public HashMap<Double, Double> selectedImmoniumPeaks;
    public HashMap<Double, Double> selectedCapYPeaks;
    public HashMap<Character, HashMap<Double, Double>> selectedSquigglePeaks;
    public int nPeaks;
    static double [] fact;

    public int compareTo(DiagnosticRecord dr) {
        return Integer.compare(this.scanNum, dr.scanNum);
    };

    public DiagnosticRecord(Spectrum spec, String ionTypes, String pepSeq, float[] mods, float dmass, int charge) {
        this.scanNum = spec.scanNum;
        this.ionTypes = ionTypes;
        this.pepSeq = pepSeq;
        this.modificationsArray = mods;
        this.modifications = formatMods(mods);
        this.dmass = dmass;
        this.charge = charge;
    }

    public DiagnosticRecord(int scanNum, ArrayList<Character> ionTypes, String pepSeq, float[] mods, float dmass,
                            int charge, float[][] immoniumPeaks, float[][] capYPeaks,
                            HashMap<Character, float[][]> squigglePeaks) {
        this.scanNum = scanNum;
        this.ionTypes = formatIonTypes(ionTypes);
        this.pepSeq = pepSeq;
        this.modificationsArray = mods;
        this.modifications = formatMods(mods);
        this.dmass = dmass;
        this.charge = charge;
        this.immoniumPeaks = immoniumPeaks;
        this.capYPeaks = capYPeaks;
        this.squigglePeaks = squigglePeaks;

    }

    public void setImmoniumPeaks(float[][] immoniumPeaks) {
        this.immoniumPeaks = immoniumPeaks;
    }

    public void setCapYPeaks(float[][] capYPeaks) {
        this.capYPeaks = capYPeaks;
    }

    public void setSquigglePeaks(HashMap<Character, float[][]> squigglePeaks) {
        this.squigglePeaks = squigglePeaks;
    }

    private HashMap<Integer, Float> formatMods(float[] mods) {
        HashMap<Integer, Float> hashMods = new HashMap<>();
        for (int i = 0; i < mods.length; i++) {
            if (mods[i] > 0.0)
                hashMods.put(i, mods[i]);
        }
        return hashMods;
    }

    private String formatIonTypes(ArrayList<Character> ionTypes) {
        StringBuffer sb = new StringBuffer();
        for (Character c : ionTypes)
            sb.append(Character.toLowerCase(c));
        return sb.toString();
    }

    public void filterIons_old(ArrayList<Double> immMasses, ArrayList<Double> capYMasses, HashMap<Character, ArrayList<Double>> squigglePeaksMasses, double tol) {
        /* Initialize variables that hold data for Mann-Whitney-U */
        this.selectedImmoniumPeaks = new HashMap<>();
        this.selectedCapYPeaks = new HashMap<>();
        this.selectedSquigglePeaks = new HashMap<>();

        /* Collect intensities in range of each immonium peak */
        for (Double peak : immMasses) {
            double cInt = 0;
            for (int i = 0; i < this.immoniumPeaks.length; i++) {
                double cDist = Math.abs(this.immoniumPeaks[i][0] - peak);
                if (cDist < tol)
                    cInt += this.immoniumPeaks[i][1];
            }
            this.selectedImmoniumPeaks.put(peak, cInt);
        }

        /* Collect intensities in range of each capY peak */
        for (Double peak : capYMasses) {
            double cInt = 0;
            for (int i = 0; i < this.capYPeaks.length; i++) {
                double cDist = Math.abs(this.capYPeaks[i][0] - peak);
                if (cDist < tol)
                    cInt += this.capYPeaks[i][1];
            }
            this.selectedCapYPeaks.put(peak, cInt);
        }

        for (Character c : squigglePeaksMasses.keySet()) {
            this.selectedSquigglePeaks.put(c, new HashMap<>());
            for (Double peak : squigglePeaksMasses.get(c)) {
                double cInt = 0;
                for (int i = 0; i < this.squigglePeaks.get(c).length; i++) {
                    double cDist = Math.abs(this.squigglePeaks.get(c)[i][0] - peak);
                    if (cDist < tol)
                        cInt += this.squigglePeaks.get(c)[i][1];
                }
                cInt /= this.pepSeq.length();
                this.selectedSquigglePeaks.get(c).put(peak, cInt);
            }
        }
    }

    public void filterIons(ArrayList<Double> immMasses, ArrayList<Double> capYMasses, HashMap<Character, ArrayList<Double>> squigglePeaksMasses, double tol) {
        /* Initialize variables that hold data for Mann-Whitney-U */
        this.selectedImmoniumPeaks = new HashMap<>();
        this.selectedCapYPeaks = new HashMap<>();
        this.selectedSquigglePeaks = new HashMap<>();

        double[][] selectImmoniumPeaks = collectImmIonIntsensities(immMasses, this.immoniumPeaks, tol);
        for (double[] peak : selectImmoniumPeaks)
            this.selectedImmoniumPeaks.put(peak[0], peak[1]);

        double[][] selectCapYPeaks = collectIonIntsensities(capYMasses, this.capYPeaks, tol);
        for (double[] peak : selectCapYPeaks)
            this.selectedCapYPeaks.put(peak[0], peak[1]);

        for (Character c : squigglePeaksMasses.keySet()) {
            this.selectedSquigglePeaks.put(c, new HashMap<>());
            double[][] selectSquigglePeaks = collectIonIntsensities(squigglePeaksMasses.get(c), this.squigglePeaks.get(c), tol);
            for (int i = 0 ; i < selectSquigglePeaks.length; i++)
                selectSquigglePeaks[i][1] /= this.pepSeq.length();
            for (double[] peak : selectSquigglePeaks)
                this.selectedSquigglePeaks.get(c).put(peak[0], peak[1]);
        }
    }

    private double[][] collectImmIonIntsensities(ArrayList<Double> searchKeyIonList, float[][] expPeakList, double tol) {

        double[][] selectedPeaks = new double[searchKeyIonList.size()][2];
        for (int i = 0; i < searchKeyIonList.size(); i++)
            selectedPeaks[i][0] = searchKeyIonList.get(i);
        Arrays.sort(selectedPeaks, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return -1*Double.compare(o2[0], o1[0]);
            }
        });

        Arrays.sort(expPeakList, new Comparator<float[]>() {
            @Override
            public int compare(float[] o1, float[] o2) {
                return -1*Double.compare(o2[0], o1[0]);
            }
        });

        int expPeaksIndx = 0;
        int maxExpPeaksIndx = expPeakList.length;
        //System.out.println("***");
        for (int i = 0; i < selectedPeaks.length; i++) {
            //System.out.println("SelectPeak\t"+i + "\t" + selectedPeaks[i][0]);
            double minVal = selectedPeaks[i][0] - selectedPeaks[i][0] * tol / 1000000;
            double maxVal = selectedPeaks[i][0] + selectedPeaks[i][0] * tol / 1000000;

            while (expPeaksIndx < maxExpPeaksIndx) {
                //System.out.println("ExpPeakPass\t"+expPeaksIndx + "\t" + expPeakList[expPeaksIndx][0]);
                if (expPeakList[expPeaksIndx][0] < minVal)
                    expPeaksIndx++;
                else
                    break;
            }
            int cIndx = expPeaksIndx;
            while (cIndx < maxExpPeaksIndx) {
                //System.out.println("ExpPeakAdd\t"+cIndx + "\t" + expPeakList[cIndx][0]);
                if (expPeakList[cIndx][0] <= maxVal) {
                    selectedPeaks[i][1] += expPeakList[cIndx][1];
                    cIndx++;
                } else {
                    break;
                }
            }
        }

        /*
        for(int i = 0; i < expPeakList.length; i++)
            System.out.print(expPeakList[i][0] + " " + expPeakList[i][1]);
        System.out.println("exp");
        for(int i = 0; i < selectedPeaks.length; i++)
            System.out.print(selectedPeaks[i][0] + " " + selectedPeaks[i][1] + " ");
        System.out.println("sel");
        */

        return selectedPeaks;

    }

    private double[][] collectCapYIonIntsensities(ArrayList<Double> searchKeyIonList, float[][] expPeakList, double tol) {

        double[][] selectedPeaks = new double[searchKeyIonList.size()][2];
        for (int i = 0; i < searchKeyIonList.size(); i++)
            selectedPeaks[i][0] = searchKeyIonList.get(i);
        Arrays.sort(selectedPeaks, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return -1*Double.compare(o2[0], o1[0]);
            }
        });

        Arrays.sort(expPeakList, new Comparator<float[]>() {
            @Override
            public int compare(float[] o1, float[] o2) {
                return -1*Double.compare(o2[0], o1[0]);
            }
        });

        int expPeaksIndx = 0;
        int maxExpPeaksIndx = expPeakList.length;
        //System.out.println("***"); //todo tol over multiple charge states
        for (int i = 0; i < selectedPeaks.length; i++) {
            //System.out.println("SelectPeak\t"+i + "\t" + selectedPeaks[i][0]);
            double precMass = calcPrecursorMass(1) + selectedPeaks[i][0]; //todo charge states
            double minVal = precMass - precMass * tol / 1000000;
            double maxVal = precMass + precMass * tol / 1000000;

            while (expPeaksIndx < maxExpPeaksIndx) {
                //System.out.println("ExpPeakPass\t"+expPeaksIndx + "\t" + expPeakList[expPeaksIndx][0]);
                if (expPeakList[expPeaksIndx][0] < minVal)
                    expPeaksIndx++;
                else
                    break;
            }
            int cIndx = expPeaksIndx;
            while (cIndx < maxExpPeaksIndx) {
                //System.out.println("ExpPeakAdd\t"+cIndx + "\t" + expPeakList[cIndx][0]);
                if (expPeakList[cIndx][0] <= maxVal) {
                    selectedPeaks[i][1] += expPeakList[cIndx][1];
                    cIndx++;
                } else {
                    break;
                }
            }
        }

        /*
        for(int i = 0; i < expPeakList.length; i++)
            System.out.print(expPeakList[i][0] + " " + expPeakList[i][1]);
        System.out.println("exp");
        for(int i = 0; i < selectedPeaks.length; i++)
            System.out.print(selectedPeaks[i][0] + " " + selectedPeaks[i][1] + " ");
        System.out.println("sel");
        */

        return selectedPeaks;

    }

    private double[][] collectIonIntsensities(ArrayList<Double> searchKeyIonList, float[][] expPeakList, double tol) {

        double[][] selectedPeaks = new double[searchKeyIonList.size()][2];
        for (int i = 0; i < searchKeyIonList.size(); i++)
            selectedPeaks[i][0] = searchKeyIonList.get(i);
        Arrays.sort(selectedPeaks, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return -1*Double.compare(o2[0], o1[0]);
            }
        });

        Arrays.sort(expPeakList, new Comparator<float[]>() {
            @Override
            public int compare(float[] o1, float[] o2) {
                return -1*Double.compare(o2[0], o1[0]);
            }
        });

        int expPeaksIndx = 0;
        int maxExpPeaksIndx = expPeakList.length;
        //System.out.println("***");
        for (int i = 0; i < selectedPeaks.length; i++) {
            //System.out.println("SelectPeak\t"+i + "\t" + selectedPeaks[i][0]);
            //todo calculate tol
            double minVal = selectedPeaks[i][0] - tol;
            double maxVal = selectedPeaks[i][0] + tol;

            while (expPeaksIndx < maxExpPeaksIndx) {
                //System.out.println("ExpPeakPass\t"+expPeaksIndx + "\t" + expPeakList[expPeaksIndx][0]);
                if (expPeakList[expPeaksIndx][0] < minVal)
                    expPeaksIndx++;
                else
                    break;
            }
            int cIndx = expPeaksIndx;
            while (cIndx < maxExpPeaksIndx) {
                //System.out.println("ExpPeakAdd\t"+cIndx + "\t" + expPeakList[cIndx][0]);
                if (expPeakList[cIndx][0] <= maxVal) {
                    selectedPeaks[i][1] += expPeakList[cIndx][1];
                    cIndx++;
                } else {
                    break;
                }
            }
        }

        /*
        for(int i = 0; i < expPeakList.length; i++)
            System.out.print(expPeakList[i][0] + " " + expPeakList[i][1]);
        System.out.println("exp");
        for(int i = 0; i < selectedPeaks.length; i++)
            System.out.print(selectedPeaks[i][0] + " " + selectedPeaks[i][1] + " ");
        System.out.println("sel");
        */

        return selectedPeaks;

    }

    //public double[] localizeRemainderMass(float dmass, char ionType) {
    public int[] localizeRemainderMass(float dmass, char ionType) {
        /* Preinitialize factorials if not done yet */
        if(fact == null) {
            fact = new double[64];
            fact[0] = 1;
            for(int i = 1; i < 64; i++)
                fact[i] = (float)(fact[i-1] * i);
        }

        double ppmTol = Double.parseDouble(PTMShepherd.getParam("spectra_ppmtol")); //todo
        /* Calc deltaHypers */
        ///double[] scores = getDeltaHyper(ppmTol, ionType, dmass);
        int[] scores = getDeltaHyper(ppmTol, ionType, dmass);
        //for (int i = 0; i < scores.length; i++) {
        //    System.out.println(scores[i]);
        //}

        return scores;
    }

    //private double[] getDeltaHyper(double ppmTol, char iType, double dmass) {
    private int[] getDeltaHyper(double ppmTol, char iType, double dmass) {
        int maxCharge = 1; //Math.min(Integer.parseInt(PTMShepherd.getParam("spectra_maxfragcharge")), charge);

        /* Locally stores important masses */
        float [] aaMasses = AAMasses.monoisotopic_masses;
        float [] fragTypeShifts = AAMasses.ionTypeShifts;

        /* init scoring variables for each pos in peptide */
        int cLen = this.pepSeq.length();
        int[] matchedIons = new int[cLen];
        double[] matchedInts = new double[cLen];
        double[] ntols = new double[cLen];
        double[] ctols = new double[cLen];

        /* calculate the Da tolerances for each position on the peptide */
        double tol;
        double nTermMass;
        if (iType == 'a' || iType == 'b' || iType =='c') {
            nTermMass = fragTypeShifts[iType - 'a'] + dmass;
            for (int ccharge = 1; ccharge <= maxCharge; ccharge++) { //todo this will break once maxcharge >1
                double cmass = (nTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
                for (int i = 0; i < cLen - 1; i++) {
                    cmass += (aaMasses[this.pepSeq.charAt(i) - 'A'] + this.modificationsArray[i]) / ccharge;
                    tol = cmass * (ppmTol / 1000000.0);
                    ntols[i] = tol;
                }
            }
        }
        double cTermMass;
        if (iType == 'x' || iType == 'y' || iType =='z') {
            cTermMass = fragTypeShifts[iType - 'x' + 3] + dmass;
            for (int ccharge = 1; ccharge <= maxCharge; ccharge++) { //todo this will break once maxcharge >1
                double cmass = (cTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
                for (int i = 0; i < cLen - 1; i++) {
                    cmass += (aaMasses[this.pepSeq.charAt(cLen - 1 - i) - 'A'] + this.modificationsArray[cLen - 1 - i]) / ccharge;
                    tol = cmass * (ppmTol / 1000000.0);
                    ctols[cLen - 1 - i] = tol;
                }
            }
        }

        /* add ions and intensities to each peptide pos */
        if (iType == 'a' || iType == 'b' || iType =='c') {
            for (int i = 0; i < this.squigglePeaks.get(iType).length; i++) {
                float[] peak = this.squigglePeaks.get(iType)[i];
                if (peak[1] > 0.0) {
                    int cPos = i % (cLen - 1);
                    if (Math.abs(peak[0] - dmass) < ntols[cPos]) {
                        matchedIons[cPos]++;
                        matchedInts[cPos] += peak[1];
                    }
                }
            }
        }
        if (iType == 'x' || iType == 'y' || iType =='z') {
            for (int i = 0; i < this.squigglePeaks.get(iType).length; i++) {
                float[] peak = this.squigglePeaks.get(iType)[i];
                if (peak[1] > 0.0) {
                    int cPos = i % (cLen - 1); // enumerated position
                    cPos = cLen - cPos - 1; // actual position in pepSeq
                    if (Math.abs(peak[0] - dmass) < ctols[cPos]) {
                        matchedIons[cPos]++;
                        matchedInts[cPos] += peak[1];
                    }
                }
            }
        }

        /* calculate scores for each peptide pos and subtract off lowest*/
        /*
        double[] scores = new double[cLen];
        double minScore = 0.0;
        for (int i = 0; i < cLen; i++) {
            scores[i] = fact[matchedIons[i]];
            scores[i] *= matchedInts[i];
            //if (scores[i] < minScore)
            //    minScore = scores[i];
        }
        //for (int i = 0; i < cLen; i++)
        //    scores[i] -= minScore;
        //todo maybe try this later
        */
        return matchedIons;
    }

    public double calcAvgFragTol(char iType, int maxCharge) {
        double avgFragTol = 0;
        if (iType == 'a' || iType == 'b' || iType =='c') {
            double[] fragTols = calcNMasses(iType, maxCharge, this.pepSeq.length());
            for (int i = 0; i < fragTols.length; i++) {
                //System.out.printf("%.4f\t", fragTols[i]);
                avgFragTol += fragTols[i] / this.pepSeq.length();
            }
            //System.out.println();
        } else if (iType == 'x' || iType == 'y' || iType =='z') {
            double[] fragTols = calcCMasses(iType, maxCharge, this.pepSeq.length());
            for (int i = 0; i < fragTols.length; i++) {
                //System.out.printf("%.4f\t", fragTols[i]);
                avgFragTol += fragTols[i] / this.pepSeq.length();
            }
            //System.out.println();
        }
        //System.out.println(iType + "\t" + avgFragTol);
        return avgFragTol;
    }

    public double calcAvgCapYTol(int minCharge, int maxCharge) {
        maxCharge = Math.min(maxCharge, this.charge);
        float [] aaMasses = AAMasses.monoisotopic_masses;
        double nTermMass = AAMasses.monoisotopic_nterm_mass;
        double averagePrecursorMass = 0;
        for (int ccharge = minCharge; ccharge <= maxCharge; ccharge++) {
            double cmass = (nTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
            for (int i = 0; i < this.pepSeq.length() - 1; i++) {
                cmass += this.modificationsArray[i] / ccharge;
                cmass += (aaMasses[this.pepSeq.charAt(i) - 'A'] + this.modificationsArray[i]) / ccharge;
            }
            averagePrecursorMass += cmass / (maxCharge - minCharge + 1);
        }

        return averagePrecursorMass;
    }

    private double calcPrecursorMass(int charge) {
        float [] aaMasses = AAMasses.monoisotopic_masses;
        double nTermMass = AAMasses.monoisotopic_nterm_mass;
        double cmass = (nTermMass + charge * AAMasses.monoisotopic_nterm_mass) / charge;
        for (int i = 0; i < this.pepSeq.length() - 1; i++) {
            cmass += this.modificationsArray[i] / charge;
            cmass += (aaMasses[this.pepSeq.charAt(i) - 'A'] + this.modificationsArray[i]) / charge;
        }

        return cmass;
    }

    public double[] calcNMasses(char iType, int maxCharge, int cLen) {
        /* calculate the Da tolerances for each position on the peptide */
        float [] aaMasses = AAMasses.monoisotopic_masses;
        float [] fragTypeShifts = AAMasses.ionTypeShifts;
        double[] nmasses = new double[cLen];
        double nTermMass;
        nTermMass = fragTypeShifts[iType - 'a'];
        for (int ccharge = 1; ccharge <= maxCharge; ccharge++) { //todo this will break once maxcharge >1
            double cmass = (nTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
            for (int i = 0; i < cLen - 1; i++) {
                cmass += (aaMasses[this.pepSeq.charAt(i) - 'A'] + this.modificationsArray[i]) / ccharge;
                nmasses[i] = cmass;
            }
        }
        return nmasses;
    }

    public double[] calcCMasses(char iType, int maxCharge, int cLen) {
        /* calculate the Da tolerances for each position on the peptide */
        float [] aaMasses = AAMasses.monoisotopic_masses;
        float [] fragTypeShifts = AAMasses.ionTypeShifts;
        double[] cmasses = new double[cLen];
        double cTermMass;
        cTermMass = fragTypeShifts[iType - 'x' + 3];
        for (int ccharge = 1; ccharge <= maxCharge; ccharge++) { //todo this will break once maxcharge >1
            double cmass = (cTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
            for (int i = 0; i < cLen - 1; i++) {
                cmass += (aaMasses[this.pepSeq.charAt(cLen - 1 - i) - 'A'] + this.modificationsArray[cLen - 1 - i]) / ccharge;
                cmasses[cLen - 1 - i] = cmass;
            }
        }
        return cmasses;
    }

}
