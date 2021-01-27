package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public class DiagnosticAnalysis {
    String dsName;
    File diagFile;
    MXMLReader mr;
    float precursorTol, spectraTol;
    int condPeaks;
    int precursorMassUnits;
    int specCol, pepCol, modPepCol, deltaCol, rtCol, intCol, pmassCol, modCol;
    double condRatio;

    double [][] peakApexBounds; //[n][3] apex, left, right
    BinDiagMetric [] binDiagMetrics;

    public DiagnosticAnalysis(String ds) throws Exception {
        this.dsName = dsName;
        diagFile = new File(dsName+".diagions");
        //get necessary params
        this.precursorTol = Float.parseFloat(PTMShepherd.getParam(""));
        this.precursorMassUnits = Integer.parseInt(PTMShepherd.getParam(""));
        this.spectraTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
        this.condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
        this.condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));
    }

    public void initializeBinBoundaries(double [][] peakApexBounds) throws Exception {
        this.peakApexBounds = peakApexBounds;
    }

    public void diagIonsPSMs(PSMFile pf, HashMap<String, File> mzMappings) throws Exception{
        HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
        specCol = pf.getColumn("Spectrum");
        pepCol = pf.getColumn("Peptide");
        deltaCol = pf.dMassCol;
        intCol = pf.getColumn("Intensity");

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
            //System.out.println(cf);
            mr = new MXMLReader(mzMappings.get(cf), Integer.parseInt(PTMShepherd.getParam("threads")));
            mr.readFully();
            ArrayList<Integer> clines = mappings.get(cf); //lines corr to curr spec file
            //for (int i = 0; i < clines.size(); i++) {//for relevant line in curr spec file
            //    StringBuffer newline = processLine(pf.data.get(clines.get(i)));
            //    out.println(newline.toString());
            }
        }
        }
    //}


//}
