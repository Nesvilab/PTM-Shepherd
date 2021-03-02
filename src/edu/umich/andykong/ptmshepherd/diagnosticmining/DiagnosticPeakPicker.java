package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

public class DiagnosticPeakPicker {
    String dsName;
    File diagFile;
    MXMLReader mr;
    String ionTypes;
    float precursorTol, spectraTol;
    int condPeaks;
    int precursorMassUnits;
    double condRatio;
    FastLocator locate;
    float noiseLevel;

    TreeMap<Integer, TreeMap<String, ArrayList<Integer>>> peakToFileToLine;

    double [][] peaks; //[3][n] apex, left, right
    BinDiagMetric [] binDiagMetrics;

    public DiagnosticPeakPicker(float noiseLevel, double [][] peakApexBounds, double peakTol, int precursorMassUnits, String ions) {
        this.peaks = peakApexBounds;
        this.noiseLevel = noiseLevel;
        this.peakToFileToLine = new TreeMap<>();
        this.ionTypes = ions;
        this.locate = new FastLocator(peakApexBounds, peakTol, precursorMassUnits);
    }

    /* Construct mass shift peak -> preprocessed file -> line mapping datastructure */
    public void addFileToIndex(String dataset) {
        String fname = dataset + ".diagions";
        long t1 = System.currentTimeMillis();
        System.out.printf("\t\tIndexing data from %s\n", dataset);

        try {
            BufferedReader in = new BufferedReader(new FileReader(fname));
            String cline;
            in.readLine();

            int lineIndx = 1;
            while ((cline = in.readLine()) != null) {
                if (cline.equals("COMPLETE")) //eof condition
                    break;
                String[] sp = cline.split("\t");

                float dmass = Float.parseFloat(sp[3]);
                int peakIndx = this.locate.getIndex(dmass);
                if (peakIndx == -1)
                    continue;

                if (!this.peakToFileToLine.containsKey(peakIndx)) {
                    TreeMap<String, ArrayList<Integer>> fileToLine = new TreeMap<>();
                    ArrayList<Integer> lineIndxs = new ArrayList<>();
                    lineIndxs.add(lineIndx);
                    fileToLine.put(fname, lineIndxs);
                    this.peakToFileToLine.put(peakIndx, fileToLine);
                } else if (!this.peakToFileToLine.get(peakIndx).containsKey(fname)) {
                    ArrayList<Integer> lineIndxs = new ArrayList<>();
                    lineIndxs.add(lineIndx);
                    this.peakToFileToLine.get(peakIndx).put(fname, lineIndxs);
                } else {
                    this.peakToFileToLine.get(peakIndx).get(fname).add(lineIndx);
                }

                lineIndx++;
            }

        } catch (Exception e) {
            System.out.println(e);
        }

        long t2 = System.currentTimeMillis();
        System.out.printf("\t\tDone indexing data from %s (%d ms)\n", dataset, t2-t1);
    }

    /* Send ions to BinDiagnosticMetric containers */
    public void process() throws IOException {
        for (Integer peakIndx : this.peakToFileToLine.keySet()) {
            double[] peakVals = new double[] {this.peaks[0][peakIndx], this.peaks[1][peakIndx], this.peaks[0][peakIndx]};
            BinDiagMetric bdMetrics  = new BinDiagMetric(peakVals, this.ionTypes);

            /* Collect file data */
            for (String fname : this.peakToFileToLine.get(peakIndx).keySet()) {
                ArrayList<String> linesData = new ArrayList<>();
                try {
                    BufferedReader in = new BufferedReader(new FileReader(fname));
                    String cline;
                    while ((cline = in.readLine()) != null) {
                        if (cline.trim().length() > 0)
                            linesData.add(cline);
                    }
                    in.close();
                } catch (Exception e) { //todo this should be handled better
                    System.out.println(e);
                    System.exit(1);
                }

                for (Integer lineIndx : this.peakToFileToLine.get(peakIndx).get(fname)) {
                    bdMetrics.addPSMToPeptideMap(linesData.get(lineIndx));
                }
            }

            bdMetrics.processPeptideMap();
        }
    }

}
