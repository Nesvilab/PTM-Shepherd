package edu.umich.andykong.ptmshepherd.diagnosticanalysis;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;

public class DiagnosticExtractor {
    private File diagFile;
    private float specTol;
    private int condPeaks;
    private double condRatio;

    public DiagnosticExtractor(float specTol, int condPeaks, double condRatio) {
        this.diagFile = new File(PTMShepherd.normFName("global.rawdiag.tsv"));
        this.specTol = specTol;
    }

    public void diagPSMs(PSMFile pf, HashMap<String, File> mzMappings, ExecutorService executorService, int numThreads) throws Exception {
        HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
        PrintWriter out = new PrintWriter(new FileWriter(this.diagFile));
        ArrayList<String> linesWithoutSpectra = new ArrayList<>();

    }
}
