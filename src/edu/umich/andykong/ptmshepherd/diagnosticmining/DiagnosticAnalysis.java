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
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class DiagnosticAnalysis {
    String dsName;
    MXMLReader mr;
    String ionTypes;
    String filterIonTypes;
    float precursorTol, spectraTol;
    int condPeaks;
    int precursorMassUnits;
    int specCol, pepCol, modPepCol, deltaCol, rtCol, intCol, pmassCol, modCol, pepMassCol;
    double condRatio;
    int minImmon = 0;
    int maxImmon = 10000;
    FastLocator locate;
    private ArrayList<DiagnosticRecord> diagnosticRecords;
    double [][] peaks; //[3][n] apex, left, right
    BinDiagMetric [] binDiagMetrics;

    public DiagnosticAnalysis(String ds) throws Exception {
        this.dsName = ds;
        //get necessary params
        this.precursorTol = Float.parseFloat(PTMShepherd.getParam("precursor_tol"));
        this.precursorMassUnits = Integer.parseInt(PTMShepherd.getParam("precursor_mass_units"));
        this.spectraTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
        this.condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
        this.condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));
        this.filterIonTypes = PTMShepherd.getParam("diagmine_filterIonTypes");
        this.ionTypes = PTMShepherd.getParam("diagmine_ionTypes");

    }

    public void initializeBinBoundaries(double [][] peakApexBounds) throws Exception {
        this.peaks = peakApexBounds;
    }

    /* This function adds a PSM list file to this DiagnosticAnalysis */
    public void diagIonsPSMs(PSMFile pf, HashMap<String, File> mzMappings, ExecutorService executorService, int nThread) throws Exception {
        /* Map PSM lines to each fraction */
        HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
        for (int i = 0; i < pf.data.size(); i++) {
            String[] sp = pf.data.get(i).split("\t");
            String bn = sp[specCol].substring(0, sp[specCol].indexOf(".")); //fraction
            if (!mappings.containsKey(bn))
                mappings.put(bn, new ArrayList<>());
            mappings.get(bn).add(i);
        }

        /* Get PSM table headers for parsing */
        specCol = pf.getColumn("Spectrum");
        pepCol = pf.getColumn("Peptide");
        deltaCol = pf.dMassCol;
        pepMassCol = pf.getPrecursorCol();
        pmassCol = pf.getColumn("Calculated Peptide Mass");
        intCol = pf.getColumn("Intensity");
        modCol = pf.getColumn("Assigned Modifications");

        /* Process spectral file one at a time */
        for (String cf : mappings.keySet()) {
            if (new File(PTMShepherd.normFName(cf+".diagBIN")).exists()) {
                System.out.println("\tFound existing cached diagnostic data for " + cf);
                continue;
            }
            long t1 = System.currentTimeMillis();
            /* Initialize new DiagnosticRecord list */
            this.diagnosticRecords = new ArrayList<>();

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
                futureList.add(executorService.submit(() -> processLinesBlock(cBlock)));
            }

            /* Wait for all processes to finish */
            for (Future future : futureList)
                future.get();

            /* Write results to DiagBINFile */
            DiagBINFile diagBinFile = new DiagBINFile(this.diagnosticRecords, PTMShepherd.normFName(cf+".diagBIN"), this.ionTypes);
            diagBinFile.writeDiagBinFile();

            long t3 = System.currentTimeMillis();
            PTMShepherd.print(String.format("\t\t%s - %d lines (%d ms reading, %d ms processing)", cf, clines.size(), t2-t1,t3-t2));
        }

    }

    public void processLinesBlock(ArrayList<String> cBlock) {
        /* Temporary arrays that hold blocks of spectra to be added to shared Collections synchronously */
        ArrayList<DiagnosticRecord> diagnosticRecords = new ArrayList<>();

        for (int i = 0; i < cBlock.size(); i++) {
            DiagnosticRecord dr = processLine(cBlock.get(i));
            if (dr == null) // Spectrum not found
                continue;
            if (!dr.isMangled) // redundant?
                diagnosticRecords.add(dr);
        }

        addAllTempDiagnosticRecords(diagnosticRecords);
    }

    /* Processes line in PSM list and turns turns it into a DiagnosticRecord of transformed spectra */
    public DiagnosticRecord processLine(String s) {
        /* Get metadata from PSM list line */
        String[] sp = s.split("\t");
        String specName = sp[specCol];
        int charge = Integer.parseInt(specName.split("\\.")[specName.split("\\.").length - 1]);
        String pepSeq = sp[pepCol];
        String mods = sp[modCol];
        String [] smods = mods.split(",");
        float dmass = Float.parseFloat(sp[deltaCol]);
        float pepMass = Float.parseFloat(sp[pmassCol]);

        /* Prep spec and normalize to base peak */
        Spectrum spec = mr.getSpectrum(reNormName(specName));
        if (spec != null)
            spec.conditionOptNorm(condPeaks, condRatio, true);
        else
            return null;

        /* Initialize DiagnosticRecord and add relevant data */
        DiagnosticRecord diagnosticRecord =  new DiagnosticRecord(spec, this.ionTypes, pepSeq, parseModifications(smods, pepSeq), dmass, charge);
        diagnosticRecord.setImmoniumPeaks(calcImmoniumPeaks(spec, this.minImmon, this.maxImmon, pepSeq, parseModifications(smods, pepSeq), this.filterIonTypes, 1, dmass, this.spectraTol)); //todo max charge
        diagnosticRecord.setCapYPeaks(calcCapYPeaks(spec, pepSeq, parseModifications(smods, pepSeq), this.filterIonTypes, 1, pepMass, dmass, this.spectraTol)); //todo max charge
        diagnosticRecord.setSquigglePeaks(calcSquigglePeaks(spec, this.spectraTol, pepSeq, smods, this.ionTypes, this.filterIonTypes, 1)); //todo max charge

        return diagnosticRecord;
    }

    private synchronized void addAllTempDiagnosticRecords(ArrayList<DiagnosticRecord> diagnosticRecordsBlock) {
        this.diagnosticRecords.addAll(diagnosticRecordsBlock);
    }

    public String reNormName(String s) {
        String[] sp = s.split("\\.");
        int sn = Integer.parseInt(sp[1]);
        //with charge state
        //return String.format("%s.%d.%d.%s",sp[0],sn,sn,sp[3]);
        //without charge state
        return String.format("%s.%d.%d", sp[0], sn, sn);
    }

    public float[][] calcImmoniumPeaks(Spectrum spec, int min, int max, String seq, float[] mods, String filterIonTypes, int maxCharge, float dmass, float specTol) {
        float[][] peaks = spec.calcImmoniumPeaks(min, max, seq, mods, filterIonTypes, maxCharge, dmass, specTol);
        return peaks;
    }

    public String peaksToString(float[][] immoniumPeaks) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < immoniumPeaks.length; i++)
            sb.append(String.format("%.4f:%.4f,", immoniumPeaks[i][0], immoniumPeaks[i][1]));
        if (sb.length() == 0)
            return sb.substring(0,0);
        return sb.substring(0, sb.length()-1);
    }

    public float[][] calcCapYPeaks(Spectrum spec, String seq, float[] mods, String filterIonTypes, int maxCharge, float pepMass, float dmass, float specTol) {
        float[][] peaks = spec.calcCapYPeaks(seq, mods, filterIonTypes, maxCharge, pepMass, dmass, specTol);
        return peaks;
    }

    public HashMap<Character, float[][]> calcSquigglePeaks(Spectrum spec, float specTol, String pepSeq, String[] smods, String ionTypes, String filterIonTypes, int maxCharge) {
        //Get dem mods and put 'em on the peptide
        float [] mods = parseModifications(smods, pepSeq);
        HashMap<Character, float[][]> squigglePeaks = spec.calcSquigglePeaks(specTol, pepSeq, mods, ionTypes, filterIonTypes, maxCharge);
        return squigglePeaks;
    }

    public float[] parseModifications(String[] smods, String pepSeq) {
        float [] mods = new float[pepSeq.length()];
        Arrays.fill(mods, 0f);
        for(int i = 0; i < smods.length; i++) {
            smods[i] = smods[i].trim();
            if(smods[i].length() == 0)
                continue;
            int p = smods[i].indexOf("(");
            int q = smods[i].indexOf(")");
            String spos = smods[i].substring(0, p).trim();
            double mass = Double.parseDouble(smods[i].substring(p+1, q).trim());
            int pos = -1;
            if(spos.equals("N-term"))
                pos = 0;
            else if(spos.equals("c"))
                pos = mods.length - 1;
            else
                pos = Integer.parseInt(spos.substring(0,spos.length()-1)) - 1;
            mods[pos] += mass;
        }
        return mods;
    }

    /* todo make iscomplete function for diagBIN files
    public void complete() throws Exception {
        PrintWriter out = new PrintWriter(new FileWriter(diagFile,true));
        out.println("COMPLETE");
        out.close();
    }

    public boolean isComplete() throws Exception {
        if(diagFile.exists()) {
            RandomAccessFile raf = new RandomAccessFile(diagFile, "r");
            raf.seek(Math.max(0, diagFile.length() - 20));
            String cline;
            while((cline = raf.readLine())!=null)
                if(cline.equals("COMPLETE")) {
                    raf.close();
                    return true;
                }
            raf.close();
            diagFile.delete();
        }
        return false;
    }
     */
}
