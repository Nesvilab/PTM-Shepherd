package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import umich.ms.datatypes.lcmsrun.Hash;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class DiagnosticPeakPicker {
    String dsName;
    File diagFile;
    MXMLReader mr;
    String ionTypes;
    char[] cIonTypes;
    float precursorTol, spectraTol;
    int condPeaks;
    int precursorMassUnits;
    double condRatio;
    FastLocator locate;
    float noiseLevel;

    TreeMap<Integer, TreeMap<String, ArrayList<Integer>>> peakToFileToScan;

    double [][] peaks; //[3][n] apex, left, right
    BinDiagMetric [] binDiagMetrics;

    public DiagnosticPeakPicker(float noiseLevel, double [][] peakApexBounds, double peakTol, int precursorMassUnits, String ions) {
        this.peaks = peakApexBounds;
        this.noiseLevel = noiseLevel;
        this.peakToFileToScan = new TreeMap<>();
        this.ionTypes = ions; //redundant
        Arrays.sort(this.cIonTypes = ions.toCharArray());
        this.binDiagMetrics = new BinDiagMetric[this.peaks[0].length];
        this.locate = new FastLocator(peakApexBounds, peakTol, precursorMassUnits);
    }

    /* Construct mass shift peak -> preprocessed file -> scan nums datastructure */
    public void addFilesToIndex(String dataset, HashMap<String, File> mzMappings, ExecutorService executorService, int nThreads) {
        long t1 = System.currentTimeMillis();
        System.out.printf("\t\tIndexing data from %s\n", dataset);

        for (String cf : mzMappings.keySet()) {
            try {
                System.out.println(cf);
                String dgbinFname = cf + ".diagBIN";
                DiagBINFile dbf = new DiagBINFile(executorService, nThreads, dgbinFname, false);
                LinkedHashMap<Integer, Float> scanToDmass = dbf.getDmasses();
                for (Integer scan : scanToDmass.keySet()) {
                    float dmass = scanToDmass.get(scan);
                    int peakIndx = this.locate.getIndex(dmass);
                    //System.out.println(peakIndx + "**" + this.peaks.length);
                    if (peakIndx == -1)
                        continue;
                    if (!this.peakToFileToScan.containsKey(peakIndx)) {
                        TreeMap<String, ArrayList<Integer>> fileToScan = new TreeMap<>();
                        ArrayList<Integer> scans = new ArrayList<>();
                        scans.add(scan);
                        fileToScan.put(dgbinFname, scans);
                        this.peakToFileToScan.put(peakIndx, fileToScan);
                    } else if (!this.peakToFileToScan.get(peakIndx).containsKey(dgbinFname)) {
                        ArrayList<Integer> lineIndxs = new ArrayList<>();
                        lineIndxs.add(scan);
                        this.peakToFileToScan.get(peakIndx).put(dgbinFname, lineIndxs);
                    } else {
                        this.peakToFileToScan.get(peakIndx).get(dgbinFname).add(scan);
                    }

                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        long t2 = System.currentTimeMillis();
        System.out.printf("\t\tDone indexing data from %s (%d ms)\n", dataset, t2-t1);
    }

    /* Send ions to BinDiagnosticMetric containers */
    public void process(ExecutorService executorService, int nThreads) throws Exception {
        /* Finds the peaks for each BinDiagnosticMetric container */
        int peakCount = 0;
        for (Integer peakIndx : this.peakToFileToScan.keySet()) {
            if (peakCount++ > 25)
                break;
            //System.out.println("**" + peakIndx + "**" + this.peaks.length);
            double[] peakVals = new double[] {this.peaks[0][peakIndx], this.peaks[1][peakIndx], this.peaks[2][peakIndx]};
            BinDiagMetric bdMetrics  = new BinDiagMetric(peakVals, this.ionTypes);

            for (String fname : this.peakToFileToScan.get(peakIndx).keySet()) {
                try {
                    DiagBINFile dgBin = new DiagBINFile(executorService, nThreads, fname, true);
                    for (Integer scan : this.peakToFileToScan.get(peakIndx).get(fname)) {
                        DiagnosticRecord dr = dgBin.getScan(scan); //todo if this is slow, only load some scans
                        if (dr != null)
                            bdMetrics.addPSMToPeptideMap(dr);
                        else
                            System.out.printf("%d not found in %s\n", scan, fname);
                        }
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
            }
            bdMetrics.processPeptideMap(executorService, nThreads);
            this.binDiagMetrics[peakIndx] = bdMetrics;
        }

        /* Tests the peaks in each BinDiagnosticMetric container against zero bin */
        int zeroBin = this.locate.getIndex(0.0); //todo, better define zero bin index
        //todo this should be functions
        peakCount = 0;
        for (Integer peakIndx : this.peakToFileToScan.keySet()) {
            if (peakCount++ > 25) //todo remove this when done
                break;
            /* Construct unified peaklists that will be used in downstream processes */
            System.out.println("\nApex\t" + this.peaks[0][peakIndx]);
            ArrayList<Double> unifiedImmPeakList = new ArrayList<>(binDiagMetrics[peakIndx].immoniumIons.filteredPeaks);
            ArrayList<Double> unifiedCapYPeakList = new ArrayList<>(binDiagMetrics[peakIndx].capYIons.filteredPeaks);
            HashMap<Character, ArrayList<Double>> unifiedSquiggleIonPeakLists = new HashMap<>();
            for (Character it : this.cIonTypes)
                unifiedSquiggleIonPeakLists.put(it, new ArrayList<>(binDiagMetrics[peakIndx].tildeIons.get(this.ionTypes.indexOf(it)).filteredPeaks));
            //todo check if including zero-bin peaks is useful
            unifiedImmPeakList.addAll(binDiagMetrics[zeroBin].immoniumIons.filteredPeaks);
            /*
            unifiedCapYPeakList.addAll(binDiagMetrics[zeroBin].capYIons.filteredPeaks);
            for (Character it : this.cIonTypes)
                unifiedSquiggleIonPeakLists.get(it).addAll(binDiagMetrics[zeroBin].tildeIons.get(this.ionTypes.indexOf(it)).filteredPeaks);
             */

            /* Remove duplicate peaks */
            Collections.sort(unifiedImmPeakList);
            for (int i = unifiedImmPeakList.size()-1; i > 0; i--) {
                if (Math.abs(unifiedImmPeakList.get(i) - unifiedImmPeakList.get(i - 1)) < 0.001)
                    unifiedImmPeakList.remove(i);
            }
            Collections.sort(unifiedCapYPeakList);
            for (int i = unifiedCapYPeakList.size()-1; i > 0; i--) {
                if (Math.abs(unifiedCapYPeakList.get(i) - unifiedCapYPeakList.get(i - 1)) < 0.001)
                    unifiedCapYPeakList.remove(i);
            }
            for (Character it : this.cIonTypes) {
                Collections.sort(unifiedSquiggleIonPeakLists.get(it));
                for (int i = unifiedSquiggleIonPeakLists.get(it).size() - 1; i > 0; i--) {
                    if (Math.abs(unifiedSquiggleIonPeakLists.get(it).get(i) - unifiedSquiggleIonPeakLists.get(it).get(i - 1)) < 0.001)
                        unifiedSquiggleIonPeakLists.get(it).remove(i);
                }
            }

            /* Initialize peak testers */
            PeakCompareTester pct = new PeakCompareTester(this.peaks[0][peakIndx], unifiedImmPeakList, unifiedCapYPeakList, unifiedSquiggleIonPeakLists); // x is baseline, y is deltamasses
            /* Add peaks to peak testers */
            Set<String> unifiedFileList = new HashSet();
            for (String fname : this.peakToFileToScan.get(peakIndx).keySet())
                unifiedFileList.add(fname);
            for (String fname : this.peakToFileToScan.get(zeroBin).keySet())
                unifiedFileList.add(fname);

            /* Add all bin-wise values to PeakCompareTester, then test */
            for (String fname : unifiedFileList) {
                try {
                    DiagBINFile dgBin = new DiagBINFile(executorService, nThreads, fname, true);
                    if (this.peakToFileToScan.get(zeroBin).containsKey(fname)) {
                        for (Integer scan : this.peakToFileToScan.get(zeroBin).get(fname)) {
                            DiagnosticRecord dr = dgBin.getScan(scan); //todo if this is slow, only load some scans
                            dr.filterIons(unifiedImmPeakList, unifiedCapYPeakList, unifiedSquiggleIonPeakLists, 7, 0.01); //todo cast last param to min pep size and other param to tolerance
                            pct.addVals(dr, true);
                        }
                    }
                    if (this.peakToFileToScan.get(peakIndx).containsKey(fname)) {
                        for (Integer scan : this.peakToFileToScan.get(peakIndx).get(fname)) {
                            DiagnosticRecord dr = dgBin.getScan(scan); //todo if this is slow, only load some scans
                            dr.filterIons(unifiedImmPeakList, unifiedCapYPeakList, unifiedSquiggleIonPeakLists, 7, 0.01); //todo cast last param to min pep size and other param to tolerance
                            pct.addVals(dr, false);
                        }
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
            }
            pct.performTests();
            this.binDiagMetrics[peakIndx].setTestResults(pct);
        }

    }

    public String getBasename(String f) {
        String baseName;
        if (f.contains("_calibrated"))
            baseName = f.substring(0, f.indexOf("_calibrated"));
        else if (f.contains("_uncalibrated"))
            baseName = f.substring(0, f.indexOf("_uncalibrated"));
        else
            baseName = f.substring(0, f.lastIndexOf("."));
        return baseName;
    }

}
