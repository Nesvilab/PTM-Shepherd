package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.core.Spectrum;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class BinDiagMetric { //todo figure out bounds in a smarter way
    double peakApex;
    double leftBound;
    double rightBound;
    public ArrayList<Spectrum> spectra;
    public ArrayList<String> specNames;
    public DiagnosticHisto immoniumIons;
    public DiagnosticHisto capYIons;
    public ArrayList<DiagnosticHisto> tildeIons;
    public String ionTypes;
    public double binMinMax[][]; /* [n][2], iontype -> {min, max} */
    public HashMap<String, ArrayList<String>> peptideMap;

    public BinDiagMetric(double[] peakBounds, String ions) {
        this.peakApex = peakBounds[0];
        this.leftBound = peakBounds[1];
        this.rightBound = peakBounds[2];
        this.ionTypes = ions;
        this.binMinMax = new double[2+ions.length()][2];
        this.peptideMap = new HashMap<>();
    }

    /* Build map of peptides so that each peptide (not PSM) can be weighted equally */
    public void addPSMToPeptideMap(String cline) {
        /* Get PSM info from table */
        String[] sp = cline.split("\t");
        String specName = sp[0];
        String charge = Character.toString(specName.charAt(specName.lastIndexOf(".")+1));
        String pep = sp[1];
        String mods = sp[2];
        String pepKey = pep + mods + charge; //peptidekey includes charge ///todo is that ideal?

        /* Check if peptkey exists and add line to pepkey */
        if (!this.peptideMap.containsKey(pepKey))
            this.peptideMap.put(pepKey, new ArrayList<>());
        this.peptideMap.get(pepKey).add(cline);
    }

    /* Sends the peptides to the histogram */
    public void processPeptideMap() {

        /* Prepopulate histo min/max for each ion type */
        for (int i = 0; i < this.binMinMax.length; i++) {
            this.binMinMax[i][0] = 1000000.0; // Big minimum
            this.binMinMax[i][1] = -1000000.0; // Small maximum
        }

        /* Find the mins and maxes of each histo */
        for (String pepKey : this.peptideMap.keySet()) {
            for (String cline : this.peptideMap.get(pepKey)) {
                String[] sp = cline.substring(0, cline.length() - 1).split("\t");
                for (int i = 4; i < sp.length; i++) { /* First 4 cols in preprocessed file are spec metadata */
                    String[] ions = sp[i].split(",");
                    for (String ion : ions) {
                        String[] peak = ion.split(":");
                        double mz = Double.parseDouble(peak[0]);
                        if (mz < this.binMinMax[i-4][0])
                            this.binMinMax[i-4][0] = mz;
                        if (mz > this.binMinMax[i-4][1])
                            this.binMinMax[i-4][1] = mz;
                    }
                }
            }
        }

        //for (int i = 0; i < this.binMinMax.length; i++) {
        //    System.out.printf("%.4f\t%.4f\n",this.binMinMax[i][0], this.binMinMax[i][1]);
        //}

        /* Initialize histograms */
        this.immoniumIons = new DiagnosticHisto(this.binMinMax[0][0], this.binMinMax[0][1], 0.01);
        this.capYIons = new DiagnosticHisto(this.binMinMax[1][0], this.binMinMax[1][1], 0.01);
        this.tildeIons = new ArrayList<>();
        for (int i = 0; i < this.ionTypes.length(); i++) { /* +2 because of immonium and Y ions in first 2 i's */
            tildeIons.add(new DiagnosticHisto(this.binMinMax[i+2][0], this.binMinMax[i+2][1], 0.01));
        }

        /* Assign data to histograms */
        for (String pepKey : this.peptideMap.keySet()) {
            double nPsms = this.peptideMap.get(pepKey).size();
            for (String cline : this.peptideMap.get(pepKey)) {
                String[] sp = cline.substring(0, cline.length() - 1).split("\t");
                for (int i = 4; i < sp.length; i++) { /* First 4 cols in preprocessed file are spec metadata */
                    String[] ions = sp[i].split(",");
                    for (String ion : ions) {
                        String[] peak = ion.split(":");
                        double mz = Double.parseDouble(peak[0]);
                        double intensity = Double.parseDouble(peak[1]) / nPsms;
                        if (i == 4)
                            this.immoniumIons.placeIon(mz, intensity);
                        else if (i == 5)
                            this.capYIons.placeIon(mz, intensity);
                        else if (i > 5) {
                            if (!(-3.5 < mz && mz < 3.5)) {
                                //System.out.println(mz);
                                this.tildeIons.get(i - 6).placeIon(mz, intensity); /* 6 = 4 metadata cols + immon + Y */
                            }
                        }
                    }
                }
            }
        }
        System.out.println(peakApex);
        System.out.println(this.tildeIons.get(0).findMax());
        System.out.println(this.tildeIons.get(0).binToMass(this.tildeIons.get(0).findMax()));



    }
}