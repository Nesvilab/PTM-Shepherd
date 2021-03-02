package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class DiagnosticAnalysis {
    String dsName;
    File diagFile;
    MXMLReader mr;
    String ionTypes;
    float precursorTol, spectraTol;
    int condPeaks;
    int precursorMassUnits;
    int specCol, pepCol, modPepCol, deltaCol, rtCol, intCol, pmassCol, modCol, pepMassCol;
    double condRatio;
    int minImmon = 0;
    int maxImmon = 500;
    FastLocator locate;

    double [][] peaks; //[3][n] apex, left, right
    BinDiagMetric [] binDiagMetrics;

    public DiagnosticAnalysis(String ds) throws Exception {
        this.dsName = ds;
        diagFile = new File(dsName+".diagions");
        //get necessary params
        this.precursorTol = Float.parseFloat(PTMShepherd.getParam("precursor_tol"));
        this.precursorMassUnits = Integer.parseInt(PTMShepherd.getParam("precursor_mass_units"));
        this.spectraTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
        this.condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
        this.condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));
        this.ionTypes = "by"; //todo
    }

    public void initializeBinBoundaries(double [][] peakApexBounds) throws Exception {
        this.peaks = peakApexBounds;
        //this.binDiagMetrics = new BinDiagMetric[this.peaks[0].length];
        //for(int i = 0; i < this.binDiagMetrics.length; i++)
        //    this.binDiagMetrics[i] = new BinDiagMetric(peakApexBounds[i]);
        //locate = new FastLocator(peaks, this.precursorTol, this.precursorMassUnits);
    }



    public void diagIonsPSMs(PSMFile pf, HashMap<String, File> mzMappings, ExecutorService executorService, int nThread) throws Exception {
        HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
        PrintWriter out = new PrintWriter(diagFile);
        StringBuffer headBuff = new StringBuffer();

        headBuff.append("spectrum\tpeptide\tmodifications\tdelta_mass\timmonium_ions\tY_ions");
        for (int i = 0; i < this.ionTypes.length(); i++)
            headBuff.append(String.format("\t%c~_ions", this.ionTypes.charAt(i)));
        out.println(headBuff.toString());

        specCol = pf.getColumn("Spectrum");
        pepCol = pf.getColumn("Peptide");
        deltaCol = pf.dMassCol;
        pepMassCol = pf.getPrecursorCol();
        pmassCol = pf.getColumn("Calculated Peptide Mass");
        intCol = pf.getColumn("Intensity");
        modCol = pf.getColumn("Assigned Modifications");

        //loop through psm file to assign spectra to bins
        for (int i = 0; i < pf.data.size(); i++) {
            String[] sp = pf.data.get(i).split("\t");
            String bn = sp[specCol].substring(0, sp[specCol].indexOf(".")); //fraction
            if (!mappings.containsKey(bn))
                mappings.put(bn, new ArrayList<>());
            mappings.get(bn).add(i);
        }

        for (String cf : mappings.keySet()) { //for file in relevant spectral files // DEBUGGING1
            long t1 = System.currentTimeMillis();
            mr = new MXMLReader(mzMappings.get(cf), Integer.parseInt(PTMShepherd.getParam("threads")));
            mr.readFully();
            long t2 = System.currentTimeMillis();
            ArrayList<Integer> clines = mappings.get(cf); //lines corr to curr spec file

            /* set up parallelization blocks */
            final int BLOCKSIZE = 10; //number of scans to be parsed per thread (to cut down on thread creation overhead)
            int nBlocks = clines.size() / (BLOCKSIZE); //number of jobs submitted to queue
            if (clines.size() % BLOCKSIZE != 0) //if there are missing scans, add one more block
                nBlocks++;
            ArrayList<Future> futureList = new ArrayList<>(nBlocks);

            /* Process PSM chunks */
            for (int i = 0; i < nBlocks; i++) {
                int startInd = i * BLOCKSIZE;
                int endInd = Math.min((i + 1) * BLOCKSIZE, clines.size() - 1);
                ArrayList<String> cBlock = new ArrayList<>();
                for (int j = startInd; j < endInd; j++)
                    cBlock.add(pf.data.get(clines.get(j)));
                futureList.add(executorService.submit(() -> processLinesBlock(cBlock, out)));
            }

            /* Wait for all processes to finish */
            for (Future future : futureList)
                future.get();

            //for (int i = 0; i < clines.size(); i++) { //for relevant line in curr spec file
            //    StringBuffer newLine = processLine(pf.data.get(clines.get(i)));
            //    out.println(newLine.toString());
            //}

            long t3 = System.currentTimeMillis();
            PTMShepherd.print(String.format("\t\t%s - %d lines (%d ms reading, %d ms processing)", cf, clines.size(), t2-t1,t3-t2));
        }
        out.close();
    }

    public void processLinesBlock(ArrayList<String> cBlock, PrintWriter out) {
        StringBuffer newBlock  = new StringBuffer();
        for (int i = 0; i < cBlock.size(); i++) {
            processLine(cBlock.get(i));
            newBlock.append(processLine(cBlock.get(i)) + "\n");
        }
        printLines(out, newBlock.toString());
    }

    private synchronized void printLines(PrintWriter out, String linesBlock) {
        out.print(linesBlock);
    }

    public String reNormName(String s) {
        String[] sp = s.split("\\.");
        int sn = Integer.parseInt(sp[1]);
        //with charge state
        //return String.format("%s.%d.%d.%s",sp[0],sn,sn,sp[3]);
        //without charge state
        return String.format("%s.%d.%d", sp[0], sn, sn);
    }

    public StringBuffer processLine(String s) {
        StringBuffer sb = new StringBuffer();
        String[] sp = s.split("\t");
        String dmass = sp[deltaCol];
        String specName = sp[specCol];
        String pepSeq = sp[pepCol];
        String [] smods = sp[modCol].split(",");
        float pepMass = Float.parseFloat(sp[pepMassCol]);

        sb.append(String.format("%s\t%s\t%s\t%s", specName, pepSeq, sp[modCol], dmass));

        Spectrum spec = mr.getSpectrum(reNormName(specName));

        if (spec == null) {
            System.out.println("SPECISSUE");
            System.out.println();}
        spec.conditionOptNorm(condPeaks, condRatio, true);

        float[][] immoniumPeaks;
        float[][] capYPeaks;
        HashMap<Character, float[][]> squigglePeaks;

        immoniumPeaks = calcImmoniumIons(spec, minImmon, maxImmon);
        sb.append(String.format("\t%s", peaksToString(immoniumPeaks)));
        capYPeaks = calcCapYIons(spec, pepMass);
        sb.append(String.format("\t%s", peaksToString(capYPeaks)));
        squigglePeaks = calcSquigglePeaks(spec, 1.0f, this.spectraTol, pepSeq, smods, ionTypes, 1); //todo is pepmass needed? todo ppmtol for spectra

        for (int i = 0; i < ionTypes.length(); i++)
            sb.append(String.format("\t%s", peaksToString(squigglePeaks.get(ionTypes.charAt(i)))));

        return sb;
    }

    public float[][] calcImmoniumIons(Spectrum spec, int min, int max) {
        float[][] peaks = spec.calcImmoniumIons(min, max);
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

    public float[][] calcCapYIons(Spectrum spec, float precursorMass) {
        float[][] peaks = spec.calcCapYIons(precursorMass);
        return peaks;
    }

    public HashMap<Character, float[][]> calcSquigglePeaks(Spectrum spec, float pepMass, float specTol, String pepSeq, String[] smods, String ionTypes, int maxCharge) { //TODO do i need pepmass?
        //Get dem mods and put 'em on the peptide
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

        HashMap<Character, float[][]> squigglePeaks = spec.calcSquigglePeaks(pepMass, specTol, pepSeq, mods, ionTypes, maxCharge);

        return squigglePeaks;
    }

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
}
