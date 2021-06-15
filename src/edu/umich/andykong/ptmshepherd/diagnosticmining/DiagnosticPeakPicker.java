package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import sun.reflect.generics.tree.Tree;
import umich.ms.datatypes.lcmsrun.Hash;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
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
    double minSignal;
    int maxPrecursorCharge;
    int MAXSCANS = 1000;

    double maxP;
    double minRbc;
    boolean twoTailedTests;

    TreeMap<Integer, TreeMap<String, ArrayList<Integer>>> peakToFileToScan;

    double[][] peaks; //[3][n] apex, left, right
    BinDiagMetric[] binDiagMetrics;

    public DiagnosticPeakPicker(double minSignal, double[][] peakApexBounds, double peakTol, int precursorMassUnits, String ions, float specTol, int maxPrecursorCharge, double maxP, double minRbc, int twoTailedTests) {
        this.peaks = peakApexBounds;
        this.minSignal = minSignal;
        this.peakToFileToScan = new TreeMap<>();
        this.ionTypes = ions; //redundant
        Arrays.sort(this.cIonTypes = ions.toCharArray());
        this.binDiagMetrics = new BinDiagMetric[this.peaks[0].length];
        this.locate = new FastLocator(peakApexBounds, peakTol, precursorMassUnits);
        this.spectraTol = specTol;
        this.maxPrecursorCharge = maxPrecursorCharge;

        this.maxP = maxP;
        this.minRbc = minRbc;
        this.twoTailedTests = twoTailedTests == 1 ? true : false;
    }

    /* Construct mass shift peak -> preprocessed file -> scan nums datastructure */
    public void addFilesToIndex(String dataset, HashMap<String, File> mzMappings, ExecutorService executorService, int nThreads) {
        long t1 = System.currentTimeMillis();
        System.out.printf("\t\tIndexing data from %s\n", dataset);

        for (String cf : mzMappings.keySet()) {
            try {
                //System.out.println(cf);
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
        System.out.printf("\t\tDone indexing data from %s (%d ms)\n", dataset, t2 - t1);
    }

    /* Send ions to BinDiagnosticMetric containers */
    public void process(ExecutorService executorService, int nThreads) throws Exception {
        /* Finds the peaks for each BinDiagnosticMetric container */
        for (Integer peakIndx : this.peakToFileToScan.keySet()) {
            double[] peakVals = new double[]{this.peaks[0][peakIndx], this.peaks[1][peakIndx], this.peaks[2][peakIndx]};
            BinDiagMetric bdMetrics = new BinDiagMetric(peakVals, this.ionTypes, this.spectraTol);

            for (String fname : this.peakToFileToScan.get(peakIndx).keySet()) {
                try {
                    ArrayList<Integer> scanNums = new ArrayList<>();
                    for (Integer scan : this.peakToFileToScan.get(peakIndx).get(fname))
                        scanNums.add(scan);
                    DiagBINFile dgBin = new DiagBINFile(executorService, nThreads, fname, false);
                    dgBin.loadDiagBinSpectra(executorService, nThreads, scanNums);
                    for (Integer scan : scanNums) {
                        DiagnosticRecord dr = dgBin.getScan(scan);
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
            bdMetrics.processPeptideMap(executorService, nThreads, this.minSignal, this.maxPrecursorCharge);
            this.binDiagMetrics[peakIndx] = bdMetrics;
        }

        /* Tests the peaks in each BinDiagnosticMetric container against zero bin */
        int zeroBin = this.locate.getIndex(0.0);
        TreeMap<Integer, TreeMap<String, ArrayList<Integer>>> peakToFileToScan = filterScanNums(this.MAXSCANS);

        //todo this should be functions
        for (Integer peakIndx : peakToFileToScan.keySet()) {
            /* Construct unified peaklists that will be used in downstream processes */
            System.out.println("\nApex\t" + this.peaks[0][peakIndx]);

            ArrayList<Double> unifiedImmPeakList = new ArrayList<>(binDiagMetrics[peakIndx].immoniumIons.filteredPeaks);
            ArrayList<Double> unifiedCapYPeakList = new ArrayList<>(binDiagMetrics[peakIndx].capYIons.filteredPeaks);
            HashMap<Character, ArrayList<Double>> unifiedSquiggleIonPeakLists = new HashMap<>();
            for (Character it : this.cIonTypes)
                unifiedSquiggleIonPeakLists.put(it, new ArrayList<>(binDiagMetrics[peakIndx].tildeIons.get(this.ionTypes.indexOf(it)).filteredPeaks));

            if (this.twoTailedTests == true) {//todo check if including zero-bin peaks is useful
                unifiedImmPeakList.addAll(binDiagMetrics[zeroBin].immoniumIons.filteredPeaks);
                unifiedCapYPeakList.addAll(binDiagMetrics[zeroBin].capYIons.filteredPeaks);
                for (Character it : this.cIonTypes)
                    unifiedSquiggleIonPeakLists.get(it).addAll(binDiagMetrics[zeroBin].tildeIons.get(this.ionTypes.indexOf(it)).filteredPeaks);

                /* Remove duplicate peaks */
                Collections.sort(unifiedImmPeakList);
                for (int i = unifiedImmPeakList.size() - 1; i > 0; i--) {
                    if (Math.abs(unifiedImmPeakList.get(i) - unifiedImmPeakList.get(i - 1)) <= 0.002)
                        unifiedImmPeakList.remove(i);
                }
                Collections.sort(unifiedCapYPeakList);
                for (int i = unifiedCapYPeakList.size() - 1; i > 0; i--) {
                    if (Math.abs(unifiedCapYPeakList.get(i) - unifiedCapYPeakList.get(i - 1)) <= 0.002)
                        unifiedCapYPeakList.remove(i);
                }
                for (Character it : this.cIonTypes) {
                    Collections.sort(unifiedSquiggleIonPeakLists.get(it));
                    for (int i = unifiedSquiggleIonPeakLists.get(it).size() - 1; i > 0; i--) {
                        if (Math.abs(unifiedSquiggleIonPeakLists.get(it).get(i) - unifiedSquiggleIonPeakLists.get(it).get(i - 1)) <= 0.005)
                            unifiedSquiggleIonPeakLists.get(it).remove(i);
                    }
                }
            }

            /* Initialize peak testers */
            PeakCompareTester pct = new PeakCompareTester(this.peaks[0][peakIndx], unifiedImmPeakList, unifiedCapYPeakList, unifiedSquiggleIonPeakLists, this.maxP, this.minRbc, this.twoTailedTests, this.spectraTol); // x is baseline, y is deltamasses
            /* Add peaks to peak testers */
            Set<String> unifiedFileList = new HashSet();
            for (String fname : peakToFileToScan.get(peakIndx).keySet())
                unifiedFileList.add(fname);
            for (String fname : peakToFileToScan.get(zeroBin).keySet())
                unifiedFileList.add(fname);

            /* Add all bin-wise values to PeakCompareTester, then test */
            for (String fname : unifiedFileList) {
                try {
                    /* Collect scan nums to parse */
                    ArrayList<Integer> scanNums = new ArrayList<>();
                    ArrayList<Integer> peakScanNums = new ArrayList<>();
                    ArrayList<Integer> zeroScanNums = new ArrayList<>();
                    if (peakToFileToScan.get(peakIndx).containsKey(fname)) {
                        for (Integer scan : peakToFileToScan.get(peakIndx).get(fname)) {
                            scanNums.add(scan);
                            peakScanNums.add(scan);
                        }
                    }
                    if (peakToFileToScan.get(zeroBin).containsKey(fname)) {
                        for (Integer scan : peakToFileToScan.get(zeroBin).get(fname)) {
                            scanNums.add(scan);
                            zeroScanNums.add(scan);
                        }
                    }
                    DiagBINFile dgBin = new DiagBINFile(executorService, nThreads, fname, false);
                    dgBin.loadDiagBinSpectra(executorService, nThreads, scanNums);

                    /* Set up multithreading structures for peak bin, get diagnosticRecords, send to threads */
                    int factor = Math.max(4, 64 / nThreads);
                    ArrayList<Future> futureList = new ArrayList<>(factor * nThreads);
                    for (int i = 0; i < factor * nThreads; i++) {
                        ArrayList<DiagnosticRecord> localDrs = new ArrayList<>();
                        int start = (peakScanNums.size() * i) / (factor * nThreads);
                        int end = (peakScanNums.size() * (i + 1)) / (factor * nThreads);
                        for (int j = start; j < end; j++)
                            localDrs.add(dgBin.getScan(peakScanNums.get(j)));
                        futureList.add(executorService.submit(() ->
                                filterSpecBlock(localDrs, unifiedImmPeakList, unifiedCapYPeakList,
                                        unifiedSquiggleIonPeakLists, spectraTol, pct, false)));  //todo tol
                    }
                    /* Get results */
                    for (Future future : futureList)
                        future.get();

                    /* Set up multithreading structures for zero bin, get diagnosticRecords, send to threads */
                    futureList = new ArrayList<>(factor * nThreads);
                    for (int i = 0; i < factor * nThreads; i++) {
                        ArrayList<DiagnosticRecord> localDrs = new ArrayList<>();
                        int start = (zeroScanNums.size() * i) / (factor * nThreads);
                        int end = (zeroScanNums.size() * (i + 1)) / (factor * nThreads);
                        for (int j = start; j < end; j++)
                            localDrs.add(dgBin.getScan(zeroScanNums.get(j)));
                        futureList.add(executorService.submit(() ->
                                filterSpecBlock(localDrs, unifiedImmPeakList, unifiedCapYPeakList,
                                        unifiedSquiggleIonPeakLists, spectraTol, pct, true)));  //todo tol
                    }
                    /* Get results */
                    for (Future future : futureList)
                        future.get();

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

    private void filterSpecBlock(ArrayList<DiagnosticRecord> localDrs, ArrayList<Double> unifiedImmPeakList, ArrayList<Double> unifiedCapYPeakList,
                                 HashMap<Character, ArrayList<Double>> unifiedSquiggleIonPeakLists, double tol, PeakCompareTester pct, boolean isControl) {
        for (DiagnosticRecord dr : localDrs)
            dr.filterIons(unifiedImmPeakList, unifiedCapYPeakList, unifiedSquiggleIonPeakLists, tol); //todo tol
        pct.addDrs(localDrs, isControl);
    }

    /* Selects MAXSCANS random scans to be included in testing for each bin */
    private TreeMap<Integer, TreeMap<String, ArrayList<Integer>>> filterScanNums(int maxScans) {
        /* Set up structure to be sorted */
        class FileScanTuple {
            public String f;
            public int s;

            FileScanTuple(String file, int scan) {
                this.f = file;
                this.s = scan;
            }
        }

        TreeMap<Integer, TreeMap<String, ArrayList<Integer>>> filteredPeakToFileToScan = new TreeMap<>();

        for (Integer peak : this.peakToFileToScan.keySet()) {
            ArrayList<FileScanTuple> fileToScan = new ArrayList<>();
            for (String file : this.peakToFileToScan.get(peak).keySet()) {
                for (Integer scan : this.peakToFileToScan.get(peak).get(file)) {
                    fileToScan.add(new FileScanTuple(file, scan));
                }
            }
            Collections.shuffle(fileToScan);

            TreeMap<String, ArrayList<Integer>> filteredFileToScan = new TreeMap<>();
            for (int i = 0; i < Math.min(maxScans, fileToScan.size()); i++) {
                String cFile = fileToScan.get(i).f;
                int cScan = fileToScan.get(i).s;
                if (!filteredFileToScan.containsKey(cFile))
                    filteredFileToScan.put(cFile, new ArrayList<>());
                filteredFileToScan.get(cFile).add(cScan);
            }
            filteredPeakToFileToScan.put(peak, filteredFileToScan);
        }

        return filteredPeakToFileToScan;
    }

    public void print(String fout) throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(fout,false));

        out.print("peak_apex\tion_type\tdiagnostic_mass\tadjusted_mass\tp_value\tauc\tprop_mod_spectra\tprop_unmod_spectra\tu_stat\tn_control\tn_test\n");
        for (int i = 0; i < this.binDiagMetrics.length; i++)
            out.print(this.binDiagMetrics[i].toString());

        out.close();
    }

    private ArrayList<Integer> filterScanNums(ArrayList<Integer> scanNums, int maxScans) {
        Collections.shuffle(scanNums);
        ArrayList<Integer> scans = new ArrayList<>(scanNums.subList(0, maxScans));
        return scans;
    }
}
