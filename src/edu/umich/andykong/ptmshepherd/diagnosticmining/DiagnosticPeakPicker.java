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

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
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
    double minSpecDiff;
    double minFoldChange;
    int minIonEvidence;
    boolean twoTailedTests;

    String[] annotations;

    HashMap<Integer, HashMap<String, Pepkey>> peakToPepkeys;
    TreeMap<Integer, TreeMap<String, ArrayList<Integer>>> peakToFileToScan;

    double[][] peaks; //[3][n] apex, left, right
    BinDiagMetric[] binDiagMetrics;



    public DiagnosticPeakPicker(double minSignal, double[][] peakApexBounds, double peakTol, int precursorMassUnits, String ions, float specTol, int maxPrecursorCharge, double maxP, double minRbc, double minSpecDiff, double minFoldChange, int minIonEvidence, int twoTailedTests, int condPeaks, double condRatio, String[] annotations) {
        this.peaks = peakApexBounds;
        this.minSignal = minSignal;
        this.minFoldChange = minFoldChange;
        this.minIonEvidence = minIonEvidence;
        this.peakToFileToScan = new TreeMap<>();
        this.peakToPepkeys = new HashMap<>();
        this.ionTypes = ions; //redundant
        Arrays.sort(this.cIonTypes = ions.toCharArray());
        this.binDiagMetrics = new BinDiagMetric[this.peaks[0].length];
        this.locate = new FastLocator(peakApexBounds, peakTol, precursorMassUnits);
        this.spectraTol = specTol;
        this.maxPrecursorCharge = maxPrecursorCharge;

        this.condPeaks = condPeaks;
        this.condRatio = condRatio;

        this.maxP = maxP;
        this.minRbc = minRbc;
        this.minSpecDiff = minSpecDiff;
        this.twoTailedTests = twoTailedTests == 1 ? true : false;

        this.annotations = annotations;
    }

    /* Construct mass shift peak -> preprocessed file -> scan nums datastructure */
    public void addFilesToIndex(String dataset, HashMap<String, File> mzMappings, ExecutorService executorService, int nThreads) {
        long t1 = System.currentTimeMillis();
        System.out.printf("\t\tIndexing data from %s\n", dataset);

        for (String cf : mzMappings.keySet()) {
            try {
                String dgbinFname = PTMShepherd.normFName(cf + ".diagBIN");
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

    public void addPepkeysToIndex(PSMFile pf) {
        int dmassCol = pf.dMassCol;
        int eValCol = pf.getColumn("Expectation");
        int pepSeqCol = pf.getColumn("Peptide");
        int chargeCol = pf.getColumn("Charge");
        int modCol = pf.getColumn("Assigned Modifications");
        int specCol = pf.getColumn("Spectrum");

        for (int i = 0; i < pf.data.size(); i++) {
            String[] sp = pf.data.get(i).split("\t");
            String charge = sp[chargeCol];
            String pepSeq = sp[pepSeqCol];
            String[] specSp = sp[specCol].split("\\.");
            String mzFile = specSp[0];
            int scanNum = Integer.parseInt(specSp[specSp.length-2]);
            String mods = sp[modCol];
            float eVal = Float.parseFloat(sp[eValCol]);
            float dmass = Float.parseFloat(sp[dmassCol]);
            int peakIndx = this.locate.getIndex(dmass);
            String pepKey = pepSeq + mods + charge;

            if (peakIndx == -1)
                continue;

            if (!this.peakToPepkeys.containsKey(peakIndx))
                this.peakToPepkeys.put(peakIndx, new HashMap<>());
            if (!this.peakToPepkeys.get(peakIndx).containsKey(pepKey))
                this.peakToPepkeys.get(peakIndx).put(pepKey, new Pepkey(mzFile, pepKey, eVal, scanNum));
            else {
                if (this.peakToPepkeys.get(peakIndx).get(pepKey).eVal > eVal)
                    this.peakToPepkeys.get(peakIndx).put(pepKey, new Pepkey(mzFile, pepKey, eVal, scanNum));
            }
        }
    }

    public void filterPepkeys() {
        this.peakToFileToScan = new TreeMap<>();

        for (Integer peakIndx : this.peakToPepkeys.keySet()) {
            this.peakToFileToScan.put(peakIndx, new TreeMap<>());
            for (String pk : this.peakToPepkeys.get(peakIndx).keySet()) {
                Pepkey tpk = this.peakToPepkeys.get(peakIndx).get(pk);
                if (!this.peakToFileToScan.get(peakIndx).containsKey(tpk.mzFile + ".diagBIN"))
                    this.peakToFileToScan.get(peakIndx).put(tpk.mzFile + ".diagBIN", new ArrayList<>());
                this.peakToFileToScan.get(peakIndx).get(tpk.mzFile + ".diagBIN").add(tpk.scan);
            }
        }
    }

    /* Send ions to BinDiagnosticMetric containers */
    public void process(ExecutorService executorService, int nThreads) throws Exception {
        /* Finds the peaks for each BinDiagnosticMetric container */
        for (Integer peakIndx : this.peakToFileToScan.keySet()) {
            double[] peakVals = new double[]{this.peaks[0][peakIndx], this.peaks[1][peakIndx], this.peaks[2][peakIndx]};
            BinDiagMetric bdMetrics = new BinDiagMetric(peakVals, this.ionTypes, this.spectraTol, this.annotations[peakIndx]);

            for (String fname : this.peakToFileToScan.get(peakIndx).keySet()) {
                try {
                    ArrayList<Integer> scanNums = new ArrayList<>();
                    for (Integer scan : this.peakToFileToScan.get(peakIndx).get(fname))
                        scanNums.add(scan);
                    DiagBINFile dgBin = new DiagBINFile(executorService, nThreads, PTMShepherd.normFName(fname), false);
                    dgBin.loadDiagBinSpectra(executorService, nThreads, scanNums);
                    for (Integer scan : scanNums) {
                        DiagnosticRecord dr = dgBin.getScan(scan);
                        if (dr == null)
                            continue;
                        if (!dr.isMangled)
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
            //System.out.println("\nApex\t" + this.peaks[0][peakIndx]);

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
                        if (Math.abs(unifiedSquiggleIonPeakLists.get(it).get(i) - unifiedSquiggleIonPeakLists.get(it).get(i - 1)) <= 0.002)
                            unifiedSquiggleIonPeakLists.get(it).remove(i);
                    }
                }
            }

            /* Initialize peak testers */
            PeakCompareTester pct = new PeakCompareTester(this.peaks[0][peakIndx], unifiedImmPeakList, unifiedCapYPeakList, unifiedSquiggleIonPeakLists, this.maxP, this.minRbc, this.minSpecDiff, this.minFoldChange, this.twoTailedTests, this.spectraTol); // x is baseline, y is deltamasses
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
                    DiagBINFile dgBin = new DiagBINFile(executorService, nThreads, PTMShepherd.normFName(fname), false);
                    dgBin.loadDiagBinSpectra(executorService, nThreads, scanNums);

                    /* Set up multithreading structures for peak bin, get diagnosticRecords, send to threads */
                    int factor = Math.max(4, 64 / nThreads);
                    ArrayList<Future> futureList = new ArrayList<>(factor * nThreads);
                    for (int i = 0; i < factor * nThreads; i++) {
                        ArrayList<DiagnosticRecord> localDrs = new ArrayList<>();
                        int start = (peakScanNums.size() * i) / (factor * nThreads);
                        int end = (peakScanNums.size() * (i + 1)) / (factor * nThreads);
                        for (int j = start; j < end; j++) {
                            DiagnosticRecord dr = dgBin.getScan(peakScanNums.get(j));
                            if (dr == null)
                                continue;
                            localDrs.add(dgBin.getScan(peakScanNums.get(j)));
                        }
                        futureList.add(executorService.submit(() ->
                                filterSpecBlock(localDrs, unifiedImmPeakList, unifiedCapYPeakList,
                                        unifiedSquiggleIonPeakLists, spectraTol, pct, false, false)));  //todo tol
                    }
                    /* Get results */
                    for (Future future : futureList)
                        future.get();

                    /* Set up multithreading structures for decoy instance of peak bin, get DRs, send to threads */
                    /*
                    futureList = new ArrayList<>(factor * nThreads);
                    for (int i = 0; i < factor * nThreads; i++) {
                        ArrayList<DiagnosticRecord> localDrs = new ArrayList<>();
                        int start = (peakScanNums.size() * i) / (factor * nThreads);
                        int end = (peakScanNums.size() * (i + 1)) / (factor * nThreads);
                        for (int j = start; j < end; j++) {
                            DiagnosticRecord dr = dgBin.getScan(peakScanNums.get(j));
                            if (dr == null)
                                continue;
                            localDrs.add(dr.clone());
                        }
                        futureList.add(executorService.submit(() ->
                                filterSpecBlock(localDrs, unifiedImmPeakList, unifiedCapYPeakList,
                                        unifiedSquiggleIonPeakLists, spectraTol, pct, false, true)));  //todo tol
                    }
                    */
                    /* Get results */
                    //for (Future future : futureList)
                    //    future.get();

                    /* Set up multithreading structures for zero bin, get diagnosticRecords, send to threads */
                    futureList = new ArrayList<>(factor * nThreads);
                    for (int i = 0; i < factor * nThreads; i++) {
                        ArrayList<DiagnosticRecord> localDrs = new ArrayList<>();
                        int start = (zeroScanNums.size() * i) / (factor * nThreads);
                        int end = (zeroScanNums.size() * (i + 1)) / (factor * nThreads);
                        for (int j = start; j < end; j++) {
                            DiagnosticRecord dr = dgBin.getScan(zeroScanNums.get(j));
                            if (dr == null)
                                continue;
                            localDrs.add(dgBin.getScan(zeroScanNums.get(j)));
                        }
                        futureList.add(executorService.submit(() ->
                                filterSpecBlock(localDrs, unifiedImmPeakList, unifiedCapYPeakList,
                                        unifiedSquiggleIonPeakLists, spectraTol, pct, true, false)));  //todo tol
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
                                 HashMap<Character, ArrayList<Double>> unifiedSquiggleIonPeakLists, double tol, PeakCompareTester pct, boolean isControl, boolean isDecoy) {
        for (DiagnosticRecord dr : localDrs) {
            if (dr == null)
                continue;
            //if (isDecoy)
            //    dr.makeDecoy();
            dr.filterIons(unifiedImmPeakList, unifiedCapYPeakList, unifiedSquiggleIonPeakLists, tol, this.minIonEvidence);
        }
        pct.addDrs(localDrs, isControl, isDecoy);
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
        /*
        PrintWriter out = new PrintWriter(new FileWriter(fout,false));

        out.print("peak_apex\tion_type\tdiagnostic_mass\tadjusted_mass\te_value\tauc\tprop_mod_spectra\tprop_unmod_spectra\tmod_spectra_int\tunmod_spectra_int\tn_control\tn_test\n");
        for (int i = 0; i < this.binDiagMetrics.length; i++) {
            out.print(this.binDiagMetrics[i].toString(false));
        }

        out.close();
        */
        PrintWriter out2 = new PrintWriter(new FileWriter(fout, false));

        out2.print("peak_apex\tmod_annotation\tion_type\t" +
                "diagnostic_mass\t" +
                "remainder_propensity\t" + "delta_diagnostic_mass\t" +
                "percent_mod\tpercent_unmod\t" +
                "avg_intensity_mod\tavg_intensity_unmod\t" +
                "auc\n");
        for (int i = 0; i < this.binDiagMetrics.length; i++) {
            out2.print(this.binDiagMetrics[i].toString());
        }

        out2.close();
    }

    /* Initialize diagnostic profile */
    public void initDiagProfRecs() {
        for (int i = 0; i < this.binDiagMetrics.length; i++)
            this.binDiagMetrics[i].initializeDiagProfRecs();
    }

    /* This function goes back into the scans to get PSM-level info on diagnostic ion propensity */
    public void diagIonsPSMs(PSMFile pf, HashMap<String, File> mzMappings, ExecutorService executorService) throws Exception {
        /* Get PSM table headers for parsing */
        int specCol = pf.getColumn("Spectrum");
        int pepCol = pf.getColumn("Peptide");
        int deltaCol = pf.dMassCol;
        int modCol = pf.getColumn("Assigned Modifications");
        int pepMassCol = pf.getColumn("Calculated Peptide Mass");

        /* Map PSM lines to each fraction */
        HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
        for (int i = 0; i < pf.data.size(); i++) {
            String[] sp = pf.data.get(i).split("\t");
            String bn = sp[specCol].substring(0, sp[specCol].indexOf(".")); //fraction
            if (!mappings.containsKey(bn))
                mappings.put(bn, new ArrayList<>());
            mappings.get(bn).add(i);
        }

        /* Process spectral files one at a time */
        for (String cf : mappings.keySet()) {
            long t1 = System.currentTimeMillis();
            mr = new MXMLReader(mzMappings.get(cf), Integer.parseInt(PTMShepherd.getParam("threads")));
            mr.readFully();
            long t2 = System.currentTimeMillis();

            ArrayList<Integer> clines = mappings.get(cf); //lines corr to curr spec file
            /* set up parallelization blocks */
            final int BLOCKSIZE = 100; //number of scans to be parsed per thread (to cut down on thread creation overhead)
            int nBlocks = clines.size() / (BLOCKSIZE); //number of jobs submitted to queue
            if (clines.size() % BLOCKSIZE != 0) //if there are missing scans, add one more block
                nBlocks++;
            ArrayList<Future> futureList = new ArrayList<>(nBlocks);

            /* Process PSM chunks and add them to diagnosticRecords*/
            for (int i = 0; i < nBlocks; i++) {
                int startInd = i * BLOCKSIZE;
                int endInd = Math.min((i + 1) * BLOCKSIZE, clines.size());
                ArrayList<String> cBlock = new ArrayList<>();
                for (int j = startInd; j < endInd; j++)
                    cBlock.add(pf.data.get(clines.get(j)));
                futureList.add(executorService.submit(() -> extractIonsBlock(cBlock, specCol, pepCol, deltaCol, modCol, pepMassCol)));
            }

            /* Wait for all processes to finish */
            for (Future future : futureList)
                future.get();

            long t3 = System.currentTimeMillis();
            PTMShepherd.print(String.format("\t\t%s - %d lines (%d ms reading, %d ms processing)", cf, clines.size(), t2-t1,t3-t2));
        }

    }

    //HERE TODO HERE
    public void extractIonsBlock(ArrayList<String> cBlock, int specCol, int pepCol, int deltaCol, int modCol, int pepMassCol) {
        for (int i = 0; i < cBlock.size(); i++)
            extractIonsLine(cBlock.get(i), specCol, pepCol, deltaCol, modCol, pepMassCol);
    }

    // Gets the intensities for diagnostic ions for a single PSM
    public void extractIonsLine(String in, int specCol, int pepCol, int deltaCol, int modCol, int pepMassCol) {
        // PSM metadata
        String[] sp = in.split("\t");
        String specName = sp[specCol];
        int charge = Integer.parseInt(specName.split("\\.")[specName.split("\\.").length - 1]);
        String pepSeq = sp[pepCol];
        String [] smods = sp[modCol].split(",");
        float dmass = Float.parseFloat(sp[deltaCol]);
        float pepMass = Float.parseFloat(sp[pepMassCol]);

        int dmassIndx = locate.getIndex(dmass);
        if (dmassIndx == -1)
            return;

        // Get spec
        Spectrum spec = mr.getSpectrum(reNormName(specName));
        if (spec == null)
            return;
        spec.condition(this.condPeaks, this.condRatio);

        // Find ions of interest
        for (DiagnosticProfileRecord dpr : binDiagMetrics[dmassIndx].diagProfRecs) {
            if (dpr.type.equals("imm")) {
                double immInt = spec.findIon(dpr.adjustedMass, this.spectraTol);
                if (immInt > 0.01) {
                    dpr.nWithIon.incrementAndGet();
                    dpr.wIonInt.addAndGet(immInt);
                }
                dpr.nTotal.incrementAndGet();
            } else if (dpr.type.equals("Y")) {
                double capYInt = spec.findIonNeutral(pepMass + dpr.adjustedMass, this.spectraTol, 1);
                if (capYInt > 0.01) {
                    dpr.nWithIon.incrementAndGet();
                    dpr.wIonInt.addAndGet(capYInt);
                }
                dpr.nTotal.incrementAndGet();
            } else {
                dpr.nTotal.incrementAndGet();
                int nUnshiftedIons = spec.getFrags(pepSeq, formatMods(smods, pepSeq), this.spectraTol, dpr.type, 0.0f);
                dpr.nUnshiftedIons.getAndAdd(pepSeq.length());
                int nShiftedIons = spec.getFrags(pepSeq, formatMods(smods, pepSeq), this.spectraTol, dpr.type, (float)dpr.adjustedMass);
                dpr.nShiftedIons.addAndGet(nShiftedIons);
                dpr.pctCoverage.addAndGet(nShiftedIons / (double) pepSeq.length());
            }
        }
    }

    public String reNormName(String s) {
        String[] sp = s.split("\\.");
        int sn = Integer.parseInt(sp[1]);
        //with charge state
        //return String.format("%s.%d.%d.%s",sp[0],sn,sn,sp[3]);
        //without charge state
        return String.format("%s.%d.%d", sp[0], sn, sn);
    }

    public float[] formatMods(String[] smods, String seq) {
        float [] mods = new float[seq.length()];
        Arrays.fill(mods, 0f);
        for(int i = 0; i < smods.length; i++) {
            smods[i] = smods[i].trim();
            if (smods[i].length() == 0)
                continue;
            int p = smods[i].indexOf("(");
            int q = smods[i].indexOf(")");
            String spos = smods[i].substring(0, p).trim();
            double mass = Double.parseDouble(smods[i].substring(p + 1, q).trim());
            int pos = -1;
            if (spos.equals("N-term"))
                pos = 0;
            else if (spos.equals("c"))
                pos = mods.length - 1;
            else
                pos = Integer.parseInt(spos.substring(0, spos.length() - 1)) - 1;
            mods[pos] += mass;
        }
        return mods;
    }

    private class Pepkey implements Comparable<Pepkey> {
        public String mzFile;
        public String pepkey;
        public float eVal;
        public int scan;

        Pepkey(String mzFile, String pepkey, float eVal, int scan) {
            this.mzFile = mzFile;
            this.pepkey = pepkey;
            this.eVal = eVal;
            this.scan = scan;
        }

        @Override
        public int compareTo(Pepkey p) {
            return Float.valueOf(p.eVal).compareTo(this.eVal);
        }
    }
}
