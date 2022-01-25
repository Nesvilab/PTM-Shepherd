package edu.umich.andykong.ptmshepherd.glyco;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class GlycoAnalysis {
    String dsName;
    File glycoFile;             // .rawglyco file
    File glycoFragmentFile;     // fragment boostrapping file
    MXMLReader mr;
    ArrayList<String> lineWithoutSpectra = new ArrayList<>();
    int totalLines;
    HashMap<String, MXMLReader> multiMr;
    float ppmTol, peakTol;
    int condPeaks;
    int specCol, pepCol, modpepCol, chargecol, deltaCol, rtCol, intCol, pmassCol, modCol;
    double condRatio;
    double[] capYShifts;// = new double[]{0,203.07937,406.15874,568.21156,730.26438,892.3172,349.137279};
    double[] oxoniumIons;//= new double[]{204.086646,186.076086,168.065526,366.139466,144.0656,138.055,512.197375,292.1026925,274.0921325,657.2349,243.026426,405.079246,485.045576,308.09761};
    double[] remainderMasses;// = new double[]{203.07937,406.15874,568.21156,730.26438,892.3172,349.137279};
    ArrayList<GlycanCandidate> glycanDatabase;
    double meanMassError;
    double massErrorWidth;
    ProbabilityTables probabilityTable;
    public static final int NUM_ADDED_GLYCO_PSM_COLUMNS = 3;
    public static final int NUM_ADDED_RAWGLYCO_COLUMNS = 5;
    boolean normYions;
    double defaultMassErrorAbsScore;
    Integer[] glycoIsotopes;
    double glycoPPMtol;
    boolean runGlycanAssignment;
    public static final double DEFAULT_GLYCO_PPM_TOL = 50;
    public static final double DEFAULT_GLYCO_FDR = 0.01;
    public static final int DEFAULT_GLYCO_DECOY_TYPE = 1;
    public static final double DEFAULT_GLYCO_ABS_SCORE_BASE = 5;
    public boolean useFragmentSpecificProbs;
    HashMap<Integer, HashMap<String, Integer>> glycanMassBinMap;
    public static final int MIN_GLYCO_PSMS_FOR_BOOTSTRAP = 5;

    // Default constructor
    public GlycoAnalysis(String dsName, boolean runGlycanAssignment, ArrayList<GlycanCandidate> glycoDatabase, ProbabilityTables inputProbabilityTable, boolean normYs, double absMassErrorDefault, Integer[] glycoIsotopes, double glycoPPMtol) {
        this.dsName = dsName;
        this.runGlycanAssignment = runGlycanAssignment;
        this.glycoFile = new File(PTMShepherd.normFName(dsName + ".rawglyco"));
        this.glycoFragmentFile = new File(PTMShepherd.normFName(dsName + ".glycofrags"));
        this.glycanDatabase = glycoDatabase;

        // init with default values, can be changed by params
        this.probabilityTable = inputProbabilityTable;
        this.normYions = normYs;
        this.defaultMassErrorAbsScore = absMassErrorDefault;
        this.glycoPPMtol = glycoPPMtol;
        this.glycoIsotopes = glycoIsotopes;
        this.useFragmentSpecificProbs = false;
        this.glycanMassBinMap = new HashMap<>();
    }

    public void glycoPSMs(PSMFile pf, HashMap<String, File> mzMappings, ExecutorService executorService, int numThreads) throws Exception {
        //open up output file
        HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
        PrintWriter out = new PrintWriter(new FileWriter(glycoFile));
        PrintWriter fragmentOut = new PrintWriter(new FileWriter(glycoFragmentFile));
        ArrayList<String> linesWithoutSpectra = new ArrayList<>();

        //get necessary params
        ppmTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
        condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
        condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));
        //cap y ions
        String[] capYstrs;
        if (PTMShepherd.getParam("cap_y_ions").length() > 0)
            capYstrs = PTMShepherd.getParam("cap_y_ions").split(",| |/");
        else
            capYstrs = new String[0];
        capYShifts = new double[capYstrs.length];
        for (int i = 0; i < capYstrs.length; i++)
            capYShifts[i] = Double.parseDouble(capYstrs[i]);
        //oxonium ions
        String[] oxStrs;
        if (PTMShepherd.getParam("diag_ions").length() > 0)
            oxStrs = PTMShepherd.getParam("diag_ions").split(",| |/");
        else
            oxStrs = new String[0];
        oxoniumIons = new double[oxStrs.length];
        for (int i = 0; i < oxStrs.length; i++)
            oxoniumIons[i] = Double.parseDouble(oxStrs[i]);
        //remainder masses
        String[] remainderStrs;
        if (PTMShepherd.getParam("remainder_masses").length() > 0)
            remainderStrs = PTMShepherd.getParam("remainder_masses").split(",| |/");
        else
            remainderStrs = new String[0];
        remainderMasses = new double[remainderStrs.length];
        for (int i = 0; i < remainderStrs.length; i++)
            remainderMasses[i] = Double.parseDouble(remainderStrs[i]);

        //write header
        StringBuilder headbuff = new StringBuilder(String.format("%s\t%s\t%s\t%s\t%s", "Spectrum", "Peptide", "Mods", "Pep Mass", "Mass Shift"));
        if (runGlycanAssignment) {
            headbuff.append("\tBest Glycan\tGlycan Score\tGlycan q-value\tBest Target Glycan\tBest Target Score");
        }
        fragmentOut.println(headbuff + "\tFragments:");
        for (double capYShift : capYShifts) headbuff.append(String.format("\tY_%.4f_intensity", capYShift));
        for (double oxoniumIon : oxoniumIons) headbuff.append(String.format("\tox_%.4f_intensity", oxoniumIon));
        for (double remainderMass : remainderMasses) headbuff.append(String.format("\tdeltascore_%.4f\tlocalization_%.4f", remainderMass, remainderMass));
        out.println(headbuff.toString());
        //get necessary col indices
        specCol = pf.getColumn("Spectrum");
        pepCol = pf.getColumn("Peptide");
        modpepCol = pf.getColumn("Modified Peptide");
        modCol = pf.getColumn("Assigned Modifications");
        deltaCol = pf.dMassCol;
        pmassCol = pf.getColumn("Calculated Peptide Mass");
        rtCol = pf.getColumn("Retention");
        intCol = pf.getColumn("Intensity");

        //map PSMs to file
        for (int i = 0; i < pf.data.size(); i++) {
            String[] sp = pf.data.get(i).split("\t");
            String bn = sp[specCol].substring(0, sp[specCol].indexOf(".")); //fraction
            if (!mappings.containsKey(bn))
                mappings.put(bn, new ArrayList<>());
            mappings.get(bn).add(i);
        }

        /* Loop through spectral files -> indexed lines in PSM -> process each line */
        for (String cf : mappings.keySet()) { //for file in relevant spectral files
            long t1 = System.currentTimeMillis();
            //System.out.println(cf);
            mr = new MXMLReader(mzMappings.get(cf), Integer.parseInt(PTMShepherd.getParam("threads")));
            mr.readFully();
            long t2 = System.currentTimeMillis();
            ArrayList<Integer> clines = mappings.get(cf); //lines corr to curr spec file

            getMassErrorWidth(pf, clines);

            /* set up parallelization blocks */
            final int BLOCKSIZE = 100; //number of scans to be parsed per thread (to cut down on thread creation overhead)
            int nBlocks = clines.size() / (BLOCKSIZE); //number of jobs submitted to queue
            if (clines.size() % BLOCKSIZE != 0) //if there are missing scans, add one more block
                nBlocks++;

            ArrayList<Future> futureList = new ArrayList<>(nBlocks);
            /* Process PSM chunks */
            for (int i = 0; i < nBlocks; i++) {
                int startInd = i * BLOCKSIZE;
                int endInd = Math.min((i + 1) * BLOCKSIZE, clines.size());
                ArrayList<String> cBlock = new ArrayList<>();
                for (int j = startInd; j < endInd; j++)
                    cBlock.add(pf.data.get(clines.get(j)));
                futureList.add(executorService.submit(() -> processLinesBlock(cBlock, out, fragmentOut)));
            }
            /* Wait for all processes to finish */
            for (Future future : futureList)
                future.get();

            long t3 = System.currentTimeMillis();
            PTMShepherd.print(String.format("\t%s - %d (%d ms, %d ms)", cf, clines.size(), t2 - t1, t3 - t2));
        }
        out.close();

        if (!linesWithoutSpectra.isEmpty()) {
            PTMShepherd.print(String.format("Could not find %d/%d (%.1f%%) spectra.\n", linesWithoutSpectra.size(), this.totalLines,
                    100.0*((double)linesWithoutSpectra.size()/this.totalLines)));
            int previewSize = Math.min(linesWithoutSpectra.size(), 5);
            PTMShepherd.print(String.format("Showing first %d of %d spectra IDs that could not be found: \n\t%s\n", previewSize, linesWithoutSpectra.size(),
                    String.join("\n\t", linesWithoutSpectra.subList(0, previewSize))));
        }
    }

    public void processLinesBlock(ArrayList<String> cBlock, PrintWriter out, PrintWriter fragmentOutWriter) {
        StringBuilder newBlock  = new StringBuilder();
        StringBuilder fragmentBlock = new StringBuilder();
        for (String line : cBlock) {
            GlycanAssignmentResult glycoResult = processLine(line);
            newBlock.append(glycoResult.fullRawGlycoString).append("\n");
            fragmentBlock.append(glycoResult.printGlycoFragmentInfo());
        }
        printLines(out, newBlock.toString());
        printLines(fragmentOutWriter, fragmentBlock.toString());
    }

    private synchronized void printLines(PrintWriter out, String linesBlock) {
        out.print(linesBlock);
    }

    /**
     * Read the generated glycofrags file to determine glycan fragment probabilities for each glycan in the database.
     * Option to save prevalence file for diagnostics/info to be added?
     * @return Map of glycan string : fragment propensities container
     */
    public HashMap<String, GlycanCandidateFragments> computeGlycanFragmentProbs() throws IOException {
        HashMap<String, GlycanCandidateFragments> glycanCandidateFragmentsMap = new HashMap<>();
        HashMap<String, ArrayList<GlycanCandidate>> glycanInputMap = new HashMap<>();    // container for glycan: glycan fragment info (read in from file)

        // read info from glycofrags file
        BufferedReader in = new BufferedReader(new FileReader(glycoFragmentFile), 1 << 22);
        String currentLine;
        String headers = in.readLine();
        // todo: find headers dynamically and merge with rawglyco file header reading to single helper method
        int fragmentStartCol = 10;
        int glycanCol = 5;
        int qValCol = 7;
        int deltaMassCol = 4;

        // read all glycan info in
        while ((currentLine = in.readLine()) != null) {
            String[] splits = currentLine.split("\t", 0);       // limit 0 to discard extra empty cells if present
            // only read lines with glycan info (after column 5)
            if (splits.length > 5) {
                // todo: add FDR check (and add q-val to the glycofrags file)
                String glycanString = splits[glycanCol];
                String[] fragmentInfo;  // String[] fragmentInfo = splits.length >= fragmentStartCol ? Arrays.copyOfRange(splits, fragmentStartCol, splits.length) : new String[]{};
                if (splits.length >= fragmentStartCol) {
                    fragmentInfo = Arrays.copyOfRange(splits, fragmentStartCol, splits.length);
                } else {
                    fragmentInfo = new String[]{};
                }
                GlycanCandidate fragmentInfoContainer = new GlycanCandidate(glycanString, fragmentInfo);
                String glycanHash = fragmentInfoContainer.hash;
                if (glycanInputMap.containsKey(glycanHash)) {
                    glycanInputMap.get(glycanHash).add(fragmentInfoContainer);
                } else {
                    ArrayList<GlycanCandidate> newList = new ArrayList<>();
                    newList.add(fragmentInfoContainer);
                    glycanInputMap.put(glycanHash, newList);
                }

                // add to delta mass map for calculating glycan prevalence priors
                double deltaMass = Double.parseDouble(splits[deltaMassCol]);
                int massBin = (int) Math.floor(deltaMass);
                if (glycanMassBinMap.containsKey(massBin)) {
                    // seen this mass bin before. Get the count-by-glycan dict and increment the count for this glycan
                    HashMap<String, Integer> massBinGlycanCounts = glycanMassBinMap.get(massBin);
                    int glycanCount = massBinGlycanCounts.getOrDefault(glycanHash, 0);
                    glycanCount++;
                    massBinGlycanCounts.put(glycanHash, glycanCount);
                } else {
                    // New mass bin. Create a new count-by-glycan dict
                    HashMap<String, Integer> massBinGlycanCounts = new HashMap<>();
                    massBinGlycanCounts.put(glycanHash, 1);
                    glycanMassBinMap.put(massBin, massBinGlycanCounts);
                }
            }
        }

        // summarize results for each glycan to get final fragment propensities
        for (Map.Entry<String, ArrayList<GlycanCandidate>> glycanEntry : glycanInputMap.entrySet()) {
            // Determine the fragment likelihoods based on all PSMs for this entry
            HashMap<String, Integer> YCounts = new HashMap<>();
            HashMap<String, Integer> OxCounts = new HashMap<>();
            ArrayList<GlycanCandidate> allPSMsWithThisGlycan = glycanEntry.getValue();
            // skip generating fragment information for glycans with too few PSMs to get reasonable values
            if (allPSMsWithThisGlycan.size() < MIN_GLYCO_PSMS_FOR_BOOTSTRAP) {
                continue;
            }

            for (GlycanCandidate inputGlycan : allPSMsWithThisGlycan) {
                // read all fragments from each input glycan into the count database
                for (GlycanFragment fragment : inputGlycan.Yfragments) {
                    String fragmentHash = fragment.toHashString();
                    int count = YCounts.getOrDefault(fragmentHash, 0);
                    count++;
                    YCounts.put(fragmentHash, count);
                }
                for (GlycanFragment fragment : inputGlycan.oxoniumFragments) {
                    String fragmentHash = fragment.toHashString();
                    int count = OxCounts.getOrDefault(fragmentHash, 0);
                    count++;
                    OxCounts.put(fragmentHash, count);
                }
            }

            // now that all fragment info from all PSMs of this glycan is collected, determine propensities for each fragment
            HashMap<String, Double> yFragmentProps = new HashMap<>();
            for (Map.Entry<String, Integer> fragmentEntry : YCounts.entrySet()) {
                // save the proportion of PSMs that had this fragment
                yFragmentProps.put(fragmentEntry.getKey(), fragmentEntry.getValue() / (double) allPSMsWithThisGlycan.size());
            }
            HashMap<String, Double> OxFragmentProps = new HashMap<>();
            for (Map.Entry<String, Integer> fragmentEntry : OxCounts.entrySet()) {
                // save the proportion of PSMs that had this fragment
                OxFragmentProps.put(fragmentEntry.getKey(), fragmentEntry.getValue() / (double) allPSMsWithThisGlycan.size());
            }

            // save determined propensities to the output container
            GlycanCandidateFragments fragmentInfo = new GlycanCandidateFragments(yFragmentProps, OxFragmentProps);
            glycanCandidateFragmentsMap.put(glycanEntry.getKey(), fragmentInfo);
        }
        return glycanCandidateFragmentsMap;
    }

    /**
     * Read rawglyco file during 2nd pass to compute FDR across whole dataset and write updated
     * information back to rawglyco file. Requires that first pass has already been done and
     * Glycans assigned to PSMs.
     *
     * @param desiredRatio: desired FDR (typically 0.01 = 1%)
     */
    public void computeGlycanFDR(double desiredRatio) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(glycoFile), 1 << 22);

        // read rawglyco file into map of spectrum index: full line (string)
        LinkedHashMap<String, String[]> glyLines = new LinkedHashMap<>();   // linkedHashMap to preserve spectrum order
        String cgline;

        // detect headers
        String[] headerSplits = in.readLine().split("\t");
        int gSpecCol = 0;
        int absScoreCol = 0;
        int bestGlycanCol = 0;
        int qValCol = 0;
        for (int i = 0; i < headerSplits.length; i++) {
            switch (headerSplits[i].trim()) {
                case "Spectrum":
                    gSpecCol = i;
                    break;
                case "Glycan Score":
                    absScoreCol = i;
                    break;
                case "Best Glycan":
                    bestGlycanCol = i;
                    break;
                case "Glycan q-value":
                    qValCol = i;
                    break;
            }
        }
        if (absScoreCol == 0 || bestGlycanCol == 0 || qValCol == 0) {
            PTMShepherd.print(String.format("Warning: rawglyco file headers not found! FDR calculation may fail for file %s\n", glycoFile));
        }

        // read file, accumulating scores
        HashMap<String, Double> scoreMap = new HashMap<>();
        int targets = 0;
        int decoys = 0;
        while ((cgline = in.readLine()) != null) {
            if (cgline.equals("COMPLETE")) {
                break;
            }
            if (cgline.startsWith("ERROR"))
                continue;
            String[] splits = cgline.split("\t", -1);
            String spectrumID = splits[gSpecCol];
            glyLines.put(spectrumID, splits);     // save full line for later editing/writing
            // only consider columns with actual glycan info
            if (!splits[bestGlycanCol].matches("") && !splits[bestGlycanCol].toLowerCase(Locale.ROOT).matches("no matches")) {
                if (!splits[qValCol].matches("")) {
                    // glycan FDR already performed on this dataset - skip
                    PTMShepherd.print("\tGlycan FDR calculation already performed, skipping");
                    return;
                }

                // detect if target or decoy and save score
                if (splits[bestGlycanCol].toLowerCase(Locale.ROOT).contains("decoy")) {
                    decoys++;
                } else {
                    targets++;
                }
                scoreMap.put(spectrumID, Double.parseDouble(splits[absScoreCol]));
            }
        }
        in.close();

        PTMShepherd.print("Calculating Glycan FDR");
        // sort scoreMap in order of ascending score
        List<Map.Entry<String, Double>> entries = new ArrayList<>(scoreMap.entrySet());
        entries.sort(Map.Entry.comparingByValue());
        Map<String, Double> sortedScoreMap = new LinkedHashMap<>();
        for (Map.Entry<String, Double> entry : entries) {
            sortedScoreMap.put(entry.getKey(), entry.getValue());
        }
        double targetDecoyRatio = decoys / (double) targets;
        if (targetDecoyRatio < desiredRatio) {
            // not enough decoys to compute FDR - already above desired ratio. Do not update table
            PTMShepherd.print(String.format("\tNot enough decoys to compute FDR at %.1f%%, started at %.2f%%", desiredRatio * 100, targetDecoyRatio * 100));
            if (desiredRatio > targetDecoyRatio * 100) {
                PTMShepherd.print(("\tNot enough decoys to compute FDR at 0.01 * initial ratio. Check data and parameters. No FDR calculation performed!\n"));
                return;
            } else {
                // only missed by a little, try reducing desired FDR to accomodate
                desiredRatio = targetDecoyRatio - (targetDecoyRatio * 0.1);
                PTMShepherd.print(String.format("\tFDR reduced to %.2f pct due to limited decoys", desiredRatio * 100));
            }
        }

        /*
            Find threshold at which target/decoy ratio hits desired value and update rawglyco lines
            NOTE: calculation is FDR <= desired ratio, so target/decoy counts updated AFTER the current PSM. This means the last
            decoy passes FDR, and the last target has the correct FDR instead of 0/0
            NOTE2: q-value is set to min(current FDR, FDR of all PSMs with lower score) to make it a step down rather than sawtooth shape
         */
        boolean foundThreshold = false;
        double currentMinQ = 1;
        for (Map.Entry<String, Double> scoreEntry : sortedScoreMap.entrySet()) {
            String[] rawGlycoLine = glyLines.get(scoreEntry.getKey());

            // update counts
            if (rawGlycoLine[bestGlycanCol].toLowerCase(Locale.ROOT).contains("decoy")) {
                decoys--;
            } else {
                targets--;
            }
            // compute TD ratio and q-val
            targetDecoyRatio = decoys / (double) targets;
            if (decoys > targets) {
                targetDecoyRatio = 1.0;     // cap FDR at 1
            } else if (targets == 0) {
                targetDecoyRatio = 0.0;     // min FDR = 0. Using else-if with the above block so that if decoys are nonzero with 0 targets, FDR = 1
            }
            double qval = Math.min(targetDecoyRatio, currentMinQ);
            if (qval < currentMinQ){
                currentMinQ = qval;
            }
            // Write q-value to output, and write q=1 for decoys
            if (rawGlycoLine[bestGlycanCol].toLowerCase(Locale.ROOT).contains("decoy")) {
                rawGlycoLine[qValCol] = "1";
            } else {
                rawGlycoLine[qValCol] = String.format("%s", qval);
            }
            
            if (!foundThreshold) {
                // still below the threshold: continue checking decoys/targets and appending 'failfdr'
                if (targetDecoyRatio <= desiredRatio) {
                    // stop here, found cutoff
                    foundThreshold = true;
                    PTMShepherd.print(String.format("\tConverged to %.1f%% FDR with %d targets and %d decoys (%d total inputs)", targetDecoyRatio * 100, targets, decoys, sortedScoreMap.size()));
                }
                if (!rawGlycoLine[bestGlycanCol].contains("FailFDR")) {
                    // only add failFDR annotation once (prevents multiple writes on re-analyses)
                    rawGlycoLine[bestGlycanCol] = "FailFDR_" + rawGlycoLine[bestGlycanCol];
                }
            }

            // update the output text with the new info
            glyLines.put(scoreEntry.getKey(), rawGlycoLine);
        }

        // write output back to rawglyco file
        PrintWriter out = new PrintWriter(new FileWriter(glycoFile));
        out.println(String.join("\t", Arrays.asList(headerSplits)));
        for (String[] glyLine : glyLines.values()) {
            out.println(String.join("\t", Arrays.asList(glyLine)));
        }
        out.flush();
        out.close();
    }

    /**
     * Determine the width of mass errors in PSMs without delta mass to use for mass error probability
     * estimation. Returns the sigma of a Gaussian distribution fit to the mass errors of all PSMs without
     * delta masses (isotope corrected)
     * @param psmFile PSM file to analyze
     * @param clines line numbers in the PSM file?
     */
    public void getMassErrorWidth(PSMFile psmFile, ArrayList<Integer> clines) {
        ArrayList<Double> massErrors = new ArrayList<>();

        // Get mass errors for PSMs with delta mass in exclusion range (-1.5 to 3.5)
        double minError = 10;
        double maxError = -10;
        for (Integer cline : clines) {//for relevant line in curr spec file
            String line = psmFile.data.get(cline);
            String[] sp = line.split("\\t");
            float deltaMass = Float.parseFloat(sp[deltaCol]);
//            float pepMass = Float.parseFloat(sp[pmassCol]);

            if (deltaMass > -1.5 && deltaMass < 3.5) {
                int isotopeError = Math.round(deltaMass);
                double massError = deltaMass - (isotopeError * AAMasses.averagineIsotopeMass);
                if (massError > maxError) {
                    maxError = massError;
                }
                if (massError < minError) {
                    minError = massError;
                }
                massErrors.add(massError);
            }
        }

        if (massErrors.size() < 200) {
            // not enough unmodified PSMs to compute mass error stats - use defaults
            PTMShepherd.print("\tNot enough unmodified PSMs to determine mass error distribution, using default values");
            massErrorWidth = 0.005;
            meanMassError = 0;
            return;
        }
        // Bin error values into a histogram
        final int numBins = 200;
        final int[] binCounts = new int[numBins];
        final double binSize = (maxError - minError) / numBins;
        for (double massError : massErrors) {
            int bin = (int) ((massError - minError) / binSize);
            // catch overflow from rounding errors
            if (bin > binCounts.length - 1)
                bin = binCounts.length - 1;
            if (bin < 0)
                bin = 0;
            binCounts[bin] += 1;
        }

        // Fit Gaussian and save center and width
        GaussianCurveFitter fitter = GaussianCurveFitter.create();
        WeightedObservedPoints massErrorObservations = new WeightedObservedPoints();
        for (int i = 0; i < binCounts.length; i++) {
            // x-value is minError + binSize*i + binSize/2 (middle of the bin), y-value is counts
            double xval = minError + i * binSize + binSize / 2.0;
            massErrorObservations.add(xval, binCounts[i]);
        }
        double[] fitParameters = fitter.fit(massErrorObservations.toList());
        meanMassError = fitParameters[1];   // output is amplitude, mean, sigma of fitted curve
        massErrorWidth = fitParameters[2];
    }


    public GlycanAssignmentResult processLine(String line) {
        StringBuilder sb = new StringBuilder();
        String[] sp = line.split("\\t");
        String seq = sp[pepCol];
        float dmass = Float.parseFloat(sp[deltaCol]);
        float pepMass = Float.parseFloat(sp[pmassCol]);
        String[] smods = sp[modCol].split(",");
        String specName = sp[specCol];

        sb.append(String.format("%s\t%s\t%s\t%.4f\t%.4f", specName, seq, sp[modCol], pepMass, dmass));
        GlycanAssignmentResult glycoResult = new GlycanAssignmentResult(seq, dmass, pepMass, sp[modCol], specName);

        Spectrum spec = mr.getSpectrum(reNormName(specName));
        if (spec == null) {
            this.lineWithoutSpectra.add(reNormName(specName));
            glycoResult.fullRawGlycoString = "ERROR";
            return glycoResult;
        }
        spec.conditionOptNorm(condPeaks, condRatio, false);

        if (runGlycanAssignment) {
            glycoResult = assignGlycanToPSM(spec, glycoResult, glycanDatabase, massErrorWidth, meanMassError);
            sb.append(glycoResult.glycanAssignmentString);
        }
        //System.out.println("got spec");
        double[] capYIonIntensities;
        double[] oxoniumIonIntensities;
        capYIonIntensities = findCapitalYIonMasses(spec, pepMass);
        oxoniumIonIntensities = findOxoniumIonMasses(spec, pepMass);

        for (double capYIonIntensity : capYIonIntensities)
            sb.append(String.format("\t%.2f", capYIonIntensity));
        for (double oxoniumIonIntensity : oxoniumIonIntensities)
            sb.append(String.format("\t%.2f", oxoniumIonIntensity));
        float[] deltaScores = new float[remainderMasses.length];
        boolean[][] isMaxScores = localizeRemainderFragments(spec, sp[pepCol], smods, deltaScores);

        for (int i = 0; i < remainderMasses.length; i++) {
            sb.append(String.format("\t%.1f", deltaScores[i]));
            StringBuilder locSb = new StringBuilder("\t");
            for (int j = 0; j < seq.length(); j++) {
                if (isMaxScores[i][j]) {
                    locSb.append(String.format("%d%c", j + 1, seq.charAt(j))); //position (1 indexed), character
                }
            }
            sb.append(locSb);
        }
        glycoResult.fullRawGlycoString = sb.toString();
        return glycoResult;
    }

    /**
     * Main glycan assignment method at PSM level. Searches Y/Oxonium ions (and eventually exact mass/isotope) to compare
     * to possible glycan candidates. Goal is to return best glycan candidate and score.
     * Formats results for writing to rawglyco file, which are later written to PSM table.
     *
     * @param spec           spectrum being searched
     * @param glycoResult    result container with spectrum information. Will have glycan results added
     * @param glycanDatabase possible glycan candidates
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     */
    public GlycanAssignmentResult assignGlycanToPSM(Spectrum spec, GlycanAssignmentResult glycoResult, ArrayList<GlycanCandidate> glycanDatabase, double massErrorWidth, double meanMassError) {
        // skip non-delta mass PSMs - leave added columns empty
        if (glycoResult.deltaMass < 3.5 && glycoResult.deltaMass > -1.5) {
            StringBuilder sb = new StringBuilder();
            for (int i=0; i < NUM_ADDED_RAWGLYCO_COLUMNS; i++){
                sb.append("\t");
            }
            glycoResult.glycanAssignmentString = sb.toString();
            return glycoResult;
        }

        // Determine possible glycan candidates from mass
        ArrayList<GlycanCandidate> searchCandidates = getMatchingGlycansByMass(glycoResult.pepMass, glycoResult.deltaMass, glycanDatabase, glycoIsotopes, glycoPPMtol);
        String output;
        if (searchCandidates.size() > 0) {
            // Search Y and oxonium ions in spectrum for each candidate
            float ppmTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
            for (GlycanCandidate candidate : searchCandidates) {
                for (GlycanFragment yFragment : candidate.Yfragments) {
                    yFragment.foundIntensity = spec.findIonNeutral(yFragment.neutralMass + glycoResult.pepMass, ppmTol, spec.charge) / spec.basePeakInt;  // sum of charge state intensities if >1 found
                }
                for (GlycanFragment oxoniumFragment : candidate.oxoniumFragments) {
                    // save oxonium ion intensity relative to base peak
                    oxoniumFragment.foundIntensity = spec.findIon(oxoniumFragment.neutralMass + AAMasses.protMass, ppmTol) / spec.basePeakInt;
                }
            }

            // score candidates and save results
            // todo: separate methods for fragment-specific analysis
            int bestCandidateIndex = 0;
            double[] scoresVsBestCandidate = new double[searchCandidates.size()];
            for (int i = 0; i < searchCandidates.size(); i++) {
                if (i == bestCandidateIndex) {
                    continue;
                }
                double comparisonScore;
                if (useFragmentSpecificProbs) {
                    comparisonScore = pairwiseCompareFragments(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass);
                } else {
                    comparisonScore = pairwiseCompareGlycans(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass, meanMassError);
                }

                if (comparisonScore < 0) {
                    // new best candidate - reset best candidate position and update scores at all other positions
                    bestCandidateIndex = i;
                    scoresVsBestCandidate[i] = -1 * comparisonScore;    // reverse the score since we have a new best candidate
                } else {
                    scoresVsBestCandidate[i] = comparisonScore;
                }
            }

            // update comparison scores against the final best candidate for those that weren't compared to best in the first pass
            for (int i = 0; i <= bestCandidateIndex; i++) {
                if (useFragmentSpecificProbs) {
                    scoresVsBestCandidate[i] = pairwiseCompareFragments(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass);
                } else {
                    scoresVsBestCandidate[i] = pairwiseCompareGlycans(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass, meanMassError);
                }
            }

            int[] sortedIndicesOfBestScores = new int[scoresVsBestCandidate.length];
            double previousBest = -1000;
            for (int i = 0; i < scoresVsBestCandidate.length; i++) {
                double currentMinScore = 1000;
                int indexOfCurrentMin = 0;
                for (int j = 0; j < scoresVsBestCandidate.length; j++) {
                    if (scoresVsBestCandidate[j] < currentMinScore && scoresVsBestCandidate[j] > previousBest) {
                        // new current max found
                        currentMinScore = scoresVsBestCandidate[j];
                        indexOfCurrentMin = j;
                    }
                }
                sortedIndicesOfBestScores[i] = indexOfCurrentMin;
                previousBest = scoresVsBestCandidate[indexOfCurrentMin];
            }

            // compute absolute score for best glycan
            // todo: add fragment-specific method
            double absoluteScore = computeAbsoluteScore(searchCandidates.get(bestCandidateIndex), glycoResult.deltaMass, massErrorWidth, meanMassError);
            glycoResult.bestCandidate = searchCandidates.get(bestCandidateIndex);
            glycoResult.glycanScore = absoluteScore;

            // output - best glycan, scores, etc back to PSM table
            // write top glycan info
            output = String.format("\t%s\t%.4f\t", searchCandidates.get(bestCandidateIndex).toString(), absoluteScore);

            // if top glycan is a decoy, also write best target and best target score to subsequent columns
            if (searchCandidates.get(bestCandidateIndex).isDecoy) {
                glycoResult.isDecoyGlycan = true;
                boolean foundTarget = false;
                for (int bestIndex : sortedIndicesOfBestScores) {
                    GlycanCandidate nextCandidate = searchCandidates.get(bestIndex);
                    if (!nextCandidate.isDecoy) {
                        // add the best target's information
                        double bestTargetScore = computeAbsoluteScore(nextCandidate, glycoResult.deltaMass, massErrorWidth, meanMassError);
                        output = String.format("%s\t%s\t%.4f", output, nextCandidate.toString(), bestTargetScore);
                        foundTarget = true;
                        glycoResult.bestTarget = nextCandidate;
                        glycoResult.bestTargetScore = bestTargetScore;
                        break;
                    }
                }
                if (!foundTarget) {
                    // no target candidates found (only matched to a decoy)
                    output = output + "\tNo target matches\t";
                }
            } else {
                glycoResult.isDecoyGlycan = false;
                output = output + "\t\t";
            }
        } else {
            output = "\tNo Matches\t\t\t\t";
        }
        glycoResult.glycanAssignmentString = output;
        return glycoResult;
    }

    /**
     * Perform pairwise comparison of two glycans. Uses sum of log probability ratios between candidates for
     * each category (mass/iso error and fragment ion) being considered. Returns a single score of combined
     * probability of first glycan candidate over second.
     *
     * @param glycan1        candidate 1
     * @param glycan2        candidate 2
     * @param deltaMass      observed delta mass
     * @return output probability score (sum of log ratios)
     */
    public double pairwiseCompareGlycans(GlycanCandidate glycan1, GlycanCandidate glycan2, double deltaMass, double meanMassError) {
        double sumLogRatio = 0;
        // Y ions
        sumLogRatio += pairwiseCompareYFragments(glycan1, glycan2, normYions);

        // oxonium ions
        sumLogRatio += pairwiseCompareOxoFragments(glycan1, glycan2);

        // isotope and mass errors
        sumLogRatio += determineIsotopeAndMassErrorProbs(glycan1, glycan2, deltaMass, meanMassError);

        return sumLogRatio;
    }

    public double pairwiseCompareFragments(GlycanCandidate glycan1, GlycanCandidate glycan2, double deltaMass) {
        // determine the overall likelihood priors of these glycans given the observed delta mass
        int glyc1Count = 0;
        int glyc2Count = 0;
        int totalGlycCount = 0;
        HashMap<String, Integer> emptyMap = new HashMap<>();
        // count glycans in this delta mass bin and nearby allowed bins
        for (int isotope : glycoIsotopes) {
            int massBin = (int) Math.floor(deltaMass + isotope);
            HashMap<String, Integer> glycanCountMap = glycanMassBinMap.getOrDefault(massBin, emptyMap);
            if (glycanCountMap.size() > 0) {
                // count instances of glycan 1, glycan 2, and all glycans
                glyc1Count = glyc1Count + glycanCountMap.getOrDefault(glycan1.hash, 0);
                glyc2Count = glyc2Count + glycanCountMap.getOrDefault(glycan2.hash, 0);
                for (int glycanCount : glycanCountMap.values()) {
                    totalGlycCount += glycanCount;
                }
            }
        }
        double propGlycan1 = glyc1Count / (double) totalGlycCount;
        double propGlycan2 = glyc2Count / (double) totalGlycCount;

        // calculate fragment-specific prob estimates based on observed fragment ions
        double sumLogRatio = 0;
        for (GlycanFragment fragment1 : glycan1.Yfragments) {
            double prob;
            if (fragment1.isAllowedFragment(glycan2)) {
                GlycanFragment fragment2 = glycan2.Yfragments[0];   // todo: fix
                if (fragment1.foundIntensity > 0) {
                    // "hit": fragment found in spectrum. Compute prob of glycans given the presence of this ion
                    prob = fragment1.propensity * propGlycan1 / (fragment1.propensity * propGlycan1 + fragment2.propensity * propGlycan2);
                } else {
                    // "miss": fragment not found. Compute prob of glycans given absence of this ion. Miss propensity = 1 - hit propensity
                    prob = (1 - fragment1.propensity) * propGlycan1 / ((1 - fragment1.propensity) * propGlycan1 + (1 - fragment2.propensity) * propGlycan2);
                }
            } else {
                // fragment only possible for glycan 1, use glycan 1 only estimate
                if (fragment1.foundIntensity > 0) {
                    // "hit": fragment found in spectrum. Compute prob of glycans given the presence of this ion
                    prob = fragment1.propensity * propGlycan1;
                } else {
                    // "miss": fragment not found. Compute prob of glycans given absence of this ion. Miss propensity = 1 - hit propensity
                    prob = (1 - fragment1.propensity) * propGlycan1;
                }
            }
            sumLogRatio += Math.log(prob);
        }
        // glycan 2 fragments
        for (GlycanFragment fragment2 : glycan2.Yfragments) {
            double prob;
            // only consider fragments unique to glycan 2 because the shared fragments have already been included from glycan 1
            if (!fragment2.isAllowedFragment(glycan2)) {
                // fragment only possible for glycan 1, use glycan 1 only estimate
                if (fragment2.foundIntensity > 0) {
                    // "hit": fragment found in spectrum. Compute prob of glycans given the presence of this ion
                    prob = fragment2.propensity * propGlycan1;
                } else {
                    // "miss": fragment not found. Compute prob of glycans given absence of this ion. Miss propensity = 1 - hit propensity
                    prob = (1 - fragment2.propensity) * propGlycan1;
                }
                sumLogRatio += Math.log(prob);
            }
        }
        // todo: add oxoniums

        // mass and isotope error
        sumLogRatio += determineIsotopeAndMassErrorProbs(glycan1, glycan2, deltaMass, meanMassError);

        return sumLogRatio;
    }

    /**
     * Compute sum log probability ratios for the compared glycans for a particular fragment type with normalization
     * of misses.
     * NOTE: does NOT allow fragment specific probabilities (for missed ions)
     *
     * @param glycan1       candidate 1
     * @param glycan2       candidate 2
     * @return sum log probability with normalization included
     */
    public double pairwiseCompareYFragments(GlycanCandidate glycan1, GlycanCandidate glycan2, boolean normYions) {
        int cand1Misses = 0;
        int cand2Misses = 0;
        int cand1Hits = 0;
        int cand2Hits = 0;
        double sumLogRatio = 0;

        // Loop over each candidate's fragments, scoring unique (i.e., not in the other candidate) fragments as hit/miss if found/not in spectrum
        for (GlycanFragment fragment1 : glycan1.Yfragments) {
            if (!fragment1.isAllowedFragment(glycan2)) {
                boolean foundInSpectrum = fragment1.foundIntensity > 0;
                if (foundInSpectrum) {
                    cand1Hits++;
                } else {
                    cand1Misses++;
                }
            }
        }
        for (GlycanFragment fragment2 : glycan2.Yfragments) {
            if (!fragment2.isAllowedFragment(glycan1)) {
                boolean foundInSpectrum = fragment2.foundIntensity > 0;
                if (foundInSpectrum) {
                    cand2Hits++;
                } else {
                    cand2Misses++;
                }
            }
        }

        // add prob normalized for number of misses for each candidate
        double cand1MissProb = glycan1.Yfragments[0].ruleProbabilities[1];
        double cand2MissProb = glycan2.Yfragments[0].ruleProbabilities[1];
        double cand1HitProb = glycan1.Yfragments[0].ruleProbabilities[0];
        double cand2HitProb = glycan2.Yfragments[0].ruleProbabilities[0];
        /* score is [log(prob1) * hits1 - log(prob2) * hits2] - [log(prob2miss) * miss2 - log(prob1miss) * miss1], or
         *      log(prob1) * hits1 - log(prob2) * hits2 - log(prob1miss) * miss1 + log(prob2miss) * miss2
         * Because log(prob) of hits is always positive, and log(prob) of miss is always negative,
         * the effect of each count on the score is:
         *      hits1 - hits2 - misses1 + misses2
         * so score is increased by hits for 1 and misses for 2, and decreased by hits for 2 and misses for 1
         */
        if (normYions) {
            sumLogRatio += Math.sqrt(cand1Hits) * Math.log(cand1HitProb);       // candidate 1 hits - added
            sumLogRatio -= Math.sqrt(cand2Hits) * Math.log(cand2HitProb);       // candidate 2 hits - subtracted
            sumLogRatio += Math.sqrt(cand1Misses) * Math.log(cand1MissProb);    // candidate 1 misses - negative value added
            sumLogRatio -= Math.sqrt(cand2Misses) * Math.log(cand2MissProb);    // candidate 2 misses - negative value subtracted
        } else {
            sumLogRatio += cand1Hits * Math.log(cand1HitProb);       // candidate 1 hits - added
            sumLogRatio -= cand2Hits * Math.log(cand2HitProb);       // candidate 2 hits - subtracted
            sumLogRatio += cand1Misses * Math.log(cand1MissProb);    // candidate 1 misses - negative value added
            sumLogRatio -= cand2Misses * Math.log(cand2MissProb);    // candidate 2 misses - negative value subtracted
        }
        return sumLogRatio;
    }

    /**
     * Compute sum log probability ratios for the compared glycans for a particular fragment type. Does NOT
     * normalize fragment miss rate. Allows fragment specific probabilities.
     *
     * @param glycan1       candidate 1
     * @param glycan2       candidate 2
     * @return sum log probability with normalization included
     */
    public double pairwiseCompareOxoFragments(GlycanCandidate glycan1, GlycanCandidate glycan2) {
        double sumLogRatio = 0;

        // Loop over each candidate's fragments, scoring unique (i.e., not in the other candidate) fragments as hit/miss if found/not in spectrum
        for (GlycanFragment fragment1 : glycan1.oxoniumFragments) {
            if (!fragment1.isAllowedFragment(glycan2)) {
                boolean foundInSpectrum = fragment1.foundIntensity > 0;
                if (foundInSpectrum) {
                    double intensityRatio = computeIntensityRatio(fragment1);
                    sumLogRatio += Math.log(fragment1.ruleProbabilities[0] * intensityRatio);    // candidate 1 hit - added
                } else {
                    sumLogRatio += Math.log(fragment1.ruleProbabilities[1]);    // candidate 1 miss - negative value added
                }
            }
        }
        for (GlycanFragment fragment2 : glycan2.oxoniumFragments) {
            if (!fragment2.isAllowedFragment(glycan1)) {
                boolean foundInSpectrum = fragment2.foundIntensity > 0;
                if (foundInSpectrum) {
                    double intensityRatio = computeIntensityRatio(fragment2);
                    sumLogRatio -= Math.log(fragment2.ruleProbabilities[0] * intensityRatio);    // candidate 2 hit - subtracted
                } else {
                    sumLogRatio -= Math.log(fragment2.ruleProbabilities[1]);    // candidate 2 miss - negative value subtracted
                }
            }
        }
        return sumLogRatio;
    }

    /**
     * Determine the probability ratio for this pairwise comparison based on isotope error. Currently
     * uses hard-coded isotope probabilities, but could be updated to get rate from dataset
     *
     * @param glycan1        glycan 1
     * @param glycan2        glycan 2
     * @param deltaMass      observed delta mass
     * @param meanMassError  mean mass error of non-delta mass peptides
     * @return probability ratio (glycan 1 over 2)
     */
    public double determineIsotopeAndMassErrorProbs(GlycanCandidate glycan1, GlycanCandidate glycan2, double deltaMass, double meanMassError) {
        // Determine isotopes
        float iso1 = (float) (deltaMass - glycan1.monoisotopicMass);
        int roundedIso1 = Math.round(iso1);
        float iso2 = (float) (deltaMass - glycan2.monoisotopicMass);
        int roundedIso2 = Math.round(iso2);

        double isotopeProbRatio = probabilityTable.isotopeProbTable.get(roundedIso1) / probabilityTable.isotopeProbTable.get(roundedIso2);

        // mass error calc
        double massProbRatio;
        if (probabilityTable.massProbScaling == 0) {
            // enable skipping mass error calc for testing
            massProbRatio = 1.0;
        } else {
            double minMassError = deltaMass * (glycoPPMtol * 0.01) * 1e-6;  // min mass error is ppmTol / 100
            double massError1 = deltaMass - glycan1.monoisotopicMass - (roundedIso1 * AAMasses.averagineIsotopeMass);
            double massStDevs1 = massError1 - meanMassError;
            if (Math.abs(massStDevs1) < minMassError)
                massStDevs1 = minMassError;
            double massError2 = deltaMass - glycan2.monoisotopicMass - (roundedIso2 * AAMasses.averagineIsotopeMass);
            double massStDevs2 = massError2 - meanMassError;
            if (Math.abs(massStDevs2) < minMassError)
                massStDevs2 = minMassError;
            massProbRatio = Math.abs(massStDevs2 / massStDevs1) * probabilityTable.massProbScaling;     // divide #2 by #1 to get ratio for likelihood of #1 vs #2, adjust by scaling factor
        }
        return Math.log(isotopeProbRatio) + Math.log(massProbRatio);
    }

    /**
     * Compute the "absolute" score of the Y ions of the provided glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     * Normalized to account for different glycan sizes. Not compatible with fragment-specific probabilities
     *
     * @param bestGlycan glycan candidate to calculate score for
     * @return absolute score
     */
    public double computeYAbsoluteScoreNormed(GlycanCandidate bestGlycan) {
        double sumLogRatio = 0;
        // Y ions - check if allowed for this composition and score if so (ignore if not)
        int hitCount = 0;
        int missCount = 0;
        // probabilities MUST be the same for all Y ions for this method to work
        double hitProb = bestGlycan.Yfragments[0].ruleProbabilities[0];
        double missProb = bestGlycan.Yfragments[0].ruleProbabilities[1];

        for (GlycanFragment yFragment : bestGlycan.Yfragments) {
            if (yFragment.foundIntensity > 0) {
                hitCount++;
            } else {
                missCount++;
            }
        }
        // normalize hit/miss counts and return final probability
        sumLogRatio = Math.sqrt(hitCount) * Math.log(hitProb) + Math.sqrt(missCount) * Math.log(missProb);
        return sumLogRatio;
    }

    /**
     * Compute the "absolute" score of the Y ions of the provided glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     * @param bestGlycan glycan candidate to calculate score for
     * @return absolute score
     */
    public double computeYAbsoluteScore(GlycanCandidate bestGlycan){
        double sumLogRatio = 0;
        for (GlycanFragment yFragment : bestGlycan.Yfragments) {
            if (yFragment.foundIntensity > 0) {
                sumLogRatio += Math.log(yFragment.ruleProbabilities[0]);     // found in spectrum - ion supports this glycan
            } else {
                sumLogRatio += Math.log(yFragment.ruleProbabilities[1]);     // not found in spectrum - ion does not support this glycan
            }
        }

        return sumLogRatio;
    }
    /**
     * Compute the "absolute" score of the Y ions of the provided glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     * @param bestGlycan glycan candidate to calculate score for
     * @return absolute score
     */
    public double computeOxoAbsoluteScore(GlycanCandidate bestGlycan){
        double sumLogRatio = 0;
        for (GlycanFragment fragment : bestGlycan.oxoniumFragments) {
            if (fragment.foundIntensity > 0) {
                double intensityRatio = computeIntensityRatio(fragment);
                sumLogRatio += Math.log(fragment.ruleProbabilities[0] * intensityRatio);     // found in spectrum - ion supports this glycan
            } else {
                sumLogRatio += Math.log(fragment.ruleProbabilities[1]);     // not found in spectrum - ion does not support this glycan
            }
        }
        return sumLogRatio;
    }

    /**
     * Compute the "absolute" score of this glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     * @param bestGlycan glycan candidate to calculate score for
     * @param deltaMass spectrum delta mass
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     * @param meanMassError mean mass error of non-delta mass peptides
     * @return absolute score
     */
    public double computeAbsoluteScore(GlycanCandidate bestGlycan, double deltaMass, double massErrorWidth, double meanMassError) {
        double sumLogRatio;
        if (normYions) {
            sumLogRatio = computeYAbsoluteScoreNormed(bestGlycan);
        } else {
            sumLogRatio = computeYAbsoluteScore(bestGlycan);
        }
        sumLogRatio += computeOxoAbsoluteScore(bestGlycan);

        // isotope and mass errors. Isotope is ratio relative to no isotope error (0)
        float iso1 = (float) (deltaMass - bestGlycan.monoisotopicMass);
        int roundedIso1 = Math.round(iso1);
        double isotopeProbRatio = probabilityTable.isotopeProbTable.get(roundedIso1) / probabilityTable.isotopeProbTable.get(0);
        sumLogRatio += Math.log(isotopeProbRatio);
        if (! (probabilityTable.massProbScaling == 0)) {
            // mass error is computed in the absolute sense - the number of std devs from mean is used instead of the ratio of two such numbers
            double massError1 = deltaMass - bestGlycan.monoisotopicMass - (roundedIso1 * AAMasses.averagineIsotopeMass);
            double massStDevs1 = (massError1 - meanMassError) / massErrorWidth;
            double massDist = Math.abs(massStDevs1) * probabilityTable.massProbScaling;
            sumLogRatio += Math.log(defaultMassErrorAbsScore / massDist);      // Compare to "default" mass error, set to 5 std devs since glycopeps tend to have larger error than regular peps, AND we're more concerned about penalizing large misses
        }
        return sumLogRatio;
    }

    /**
     * Helper method to compute the ratio of observed to expected intensity for a fragment.
     * NOTE: If the observed intensity is less than the "critical point" value (i.e., the product
     * of hit probability * intensity ratio < 1), the intensity ratio is set to the critical point value.
     * This makes it so that finding a low-intensity peak will have no effect on score rather than causing a
     * reduction in score.
     * @param fragment fragment of interest with intensity information already stored
     * @return intensity ratio
     */
    public double computeIntensityRatio(GlycanFragment fragment) {
        double intensityRatio;
        if (fragment.expectedIntensity > 0) {
            intensityRatio = fragment.foundIntensity / fragment.expectedIntensity;
            double criticalValue = 1 / fragment.ruleProbabilities[0];
            if (intensityRatio < criticalValue) {
                // beyond critical point - low intensity of hit will cause log to change sign. Cap at no effect rather than allowing sign to go negative
                intensityRatio = criticalValue;
            }
        } else {
            // ignore if parameter not provided
            intensityRatio = 1;
        }
        return intensityRatio;
    }

    /**
     * Get glycan candidates to consider for a given delta mass and isotope errors/mass tolerance.
     * Might optimize for speed at some point by indexing glycan database by mass (if needed)
     * @param pepMass peptide mass - needed for correct PPM error calculation
     * @param deltaMass delta mass being searched
     * @param glycanDatabase list of glycan candidates
     * @param isotopesToSearch list of isotope errors to consider
     * @param ms1TolerancePPM MS1 tolerance to consider around delta mass and isotope errors
     * @return list of glycan candidates with masses within the delta mass + iso errors and tolerance
     */
    public ArrayList<GlycanCandidate> getMatchingGlycansByMass(double pepMass, double deltaMass, ArrayList<GlycanCandidate> glycanDatabase, Integer[] isotopesToSearch, double ms1TolerancePPM) {
        ArrayList<GlycanCandidate> matchingGlycans = new ArrayList<>();
        for (int isotope : isotopesToSearch) {
            // add isotope error, which is recorded as an increase relative to delta mass
            double isotopeCorrMass = deltaMass - (isotope * AAMasses.averagineIsotopeMass) + pepMass;   // add peptide mass to get correct PPM calculation
            double massRangeDa = isotopeCorrMass * 0.000001 * ms1TolerancePPM;
            double massLo = isotopeCorrMass - massRangeDa - pepMass;    // remove pep mass after PPM calc for final calculation
            double massHi = isotopeCorrMass + massRangeDa - pepMass;
            for (GlycanCandidate glycan : glycanDatabase) {
                // see if mass within specified ranges
                if (glycan.monoisotopicMass >= massLo && glycan.monoisotopicMass <= massHi) {
                    // match. todo: check duplicates (could be if user inputs them)
                    // add copy of candidate to allow multi-threading without competing access
                    matchingGlycans.add(new GlycanCandidate(glycan.glycanComposition, glycan.isDecoy, glycan.monoisotopicMass, glycan.Yfragments, glycan.oxoniumFragments));
                }
            }
        }
        return matchingGlycans;
    }

    public double[] findCapitalYIonMasses(Spectrum spec, double pepMass) {
        //implement charge states //todo

        int normToBasePeak = Integer.parseInt(PTMShepherd.getParam("glyco_cap_y_ions_normalize"));
        //System.out.println(normToBasePeak);

        //initialize final capYion masses
        double[] capYIons = new double[capYShifts.length];
        double [] capYIonIntensities = new double[capYShifts.length];
        for (int i = 0; i < capYIons.length; i++)
            capYIons[i] = capYShifts[i] + pepMass;
        //find capital Y ion intensities
        for (int i = 0; i < capYIons.length; i++) {
            //System.out.println(capYIons[i]);
            capYIonIntensities[i] = spec.findIonNeutral(capYIons[i],
                    Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol")),
                    Integer.parseInt(PTMShepherd.getParam("spectra_maxPrecursorCharge"))); //todo simplify parameter calling
            if (normToBasePeak == 1) {
                //System.out.print(capYIonIntensities[i]);
                //System.out.println(" 1");
                //System.out.print(spec.findBasePeakInt());
                //System.out.println(" 2");
                capYIonIntensities[i] /= spec.basePeakInt;
                capYIonIntensities[i] *= 100.0;
                //System.out.print(capYIonIntensities[i]);
                //System.out.println(" 3");
            }
        }
        return capYIonIntensities;
    }

    public double[] findOxoniumIonMasses(Spectrum spec, double pepMass) {
        //initialize oxonium masses //todo
        //initialize capYion masses /todo
        //implement charge states //todo
        //initialize oxonium ion intensities
        int normToBasePeak = Integer.parseInt(PTMShepherd.getParam("glyco_diag_ions_normalize"));
        double[] oxoniumIonIntensities = new double[oxoniumIons.length];
        //for ion in oxonium masses/capYions
        for (int i = 0; i < oxoniumIons.length; i++) {
            oxoniumIonIntensities[i] = spec.findIon(oxoniumIons[i], Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"))); //todo simplify parameter calling
            if (normToBasePeak == 1) {
                //System.out.print(oxoniumIonIntensities[i]);
                //System.out.println(" 1");
                //System.out.print(spec.findBasePeakInt());
                //System.out.println(" 2");
                oxoniumIonIntensities[i] /= spec.basePeakInt;
                oxoniumIonIntensities[i] *= 100.0;
                //System.out.print(oxoniumIonIntensities[i]);
                //System.out.println(" 3");
            }
        }
        return oxoniumIonIntensities;
    }

    public boolean[][] localizeRemainderFragments(Spectrum spec, String seq, String[] smods, float[] deltaScores) {
        //initialize allowed positions
        boolean [] allowedPoses = SiteLocalization.parseAllowedPositions(seq, PTMShepherd.getParam("localization_allowed_res"));
        //initialize remainder delta scores
        //double[] remainderDscores = new double[remainderMasses.length];
        //add variable and fixed mods to frag masses for peptide
        float [] mods = new float[seq.length()];
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
            if(spos.equals("N-term")) {
                pos = 0;
            }
            else if(spos.equals("c")) {
                pos = mods.length - 1;
            }
            else
                pos = Integer.parseInt(spos.substring(0,spos.length()-1)) - 1;
            mods[pos] += mass;
        }
        //iterate through remainder masses
        //these 3 variables store values for each remainder mass
        float [] maxScores = new float[remainderMasses.length];
        int [] maxFrags = new int[remainderMasses.length];
        boolean [][] isMaxScores = new boolean[remainderMasses.length][seq.length()]; //1 if localized AND = max score
        //these 3 variables store values that are constant for the PSM

        float baseScore = spec.getHyper(seq, mods, ppmTol);
        int baseFrags = spec.getFrags(seq, mods, ppmTol);
        //these 3 variables need to be reinitialized every remainder mass
        float [] scores;
        int [] frags;
        double dmass;
        //begin iterate through remainder masses
        for (int i = 0; i < remainderMasses.length; i++){
            //reinit for each remainder mass
            dmass = remainderMasses[i];
            scores = new float[seq.length()];
            frags = new int[seq.length()];
            maxScores[i] = baseScore;
            maxFrags[i] = baseFrags;
            //localize at each position
            for(int j = 0; j < seq.length(); j++) {
                if (allowedPoses[j])
                    mods[j] += dmass;
                scores[j] = spec.getHyper(seq, mods, ppmTol);
                //System.out.println(scores[j] + "score");
                frags[j] = spec.getFrags(seq, mods, ppmTol);
                if(frags[j] > maxFrags[i])
                    maxFrags[i] = frags[j];
                if(scores[j] > maxScores[i])
                    maxScores[i] = scores[j];
                if (allowedPoses[j])
                    mods[j] -= dmass;
            }
            //System.out.println(maxScores[i]+"maxscore");
            //determine if localized and record max positions
            if (maxScores[i] > baseScore) {
                deltaScores[i] = maxScores[i] - baseScore;
                for (int j = 0; j < seq.length(); j++) {
                    if (scores[j] == maxScores[i]) {
                        isMaxScores[i][j] = true;
                    } else {
                        isMaxScores[i][j] = false;
                    }
                }
            } else {
                for (int j = 0; j < seq.length(); j++) {
                    isMaxScores[i][j] = false;
                }
            }
        }
        return isMaxScores;
    }

    public String reNormName(String s) {
        String[] sp = s.split("\\.");
        int sn = Integer.parseInt(sp[1]);
        //with charge state
        //return String.format("%s.%d.%d.%s",sp[0],sn,sn,sp[3]);
        //without charge state
        return String.format("%s.%d.%d", sp[0], sn, sn);
    }

    public boolean isComplete() throws Exception {
        if(glycoFile.exists()) {
            RandomAccessFile raf = new RandomAccessFile(glycoFile, "r");
            raf.seek(Math.max(0, glycoFile.length() - 20));
            String cline;
            while((cline = raf.readLine())!=null)
                if(cline.equals("COMPLETE")) {
                    raf.close();
                    return true;
                }
            raf.close();
            glycoFile.delete();
        }
        return false;
    }

    public void complete() throws Exception {
        PrintWriter out = new PrintWriter(new FileWriter(glycoFile,true));
        out.println("COMPLETE");
        out.close();
    }

    public void updateGlycoProfiles(GlycoProfile[] profiles) throws Exception {
        BufferedReader in = new BufferedReader(new FileReader(glycoFile));
        String cline;
        in.readLine();
        int glycanAssignmentAddedLines = runGlycanAssignment ? NUM_ADDED_RAWGLYCO_COLUMNS : 0;     // number of added lines to rawglyco file - 0 if not running glycan assignment
        while ((cline = in.readLine()) != null) {
            if (cline.equals("COMPLETE"))
                break;
            if (cline.startsWith("Spectrum"))
                continue;
            if (cline.startsWith("ERROR"))
                continue;
            String[] sp = cline.split("\\t");
            double md = Double.parseDouble(sp[4]);
            for (int i = 0; i < profiles.length; i++) {
                int cind = profiles[i].locate.getIndex(md);
                if (cind != -1) {
                    profiles[i].records[cind].updateWithLine(sp, glycanAssignmentAddedLines);
                }
            }
        }
        in.close();
    }


    /**
     * Find the first location of N-glycan sequon (N-X-S/T, X is not P) in the provided peptide sequence.
     * @param pepSeq peptide sequence string to search
     * @return index of N in the sequon. Position is 0-indexed
     */
    public static int findNGlycSequon(String pepSeq) {
        byte[] pepseq_glyco = pepSeq.getBytes();
        boolean xFlag = false;
        boolean asnFlag = false;
        for (int i = 0; i < pepseq_glyco.length; ++i) {
            // Check the sequon: look for next part of the sequence (if nothing, look for N; if N, look for X; etc)
            if (xFlag) {
                // Found N-X, look for S/T
                if (pepseq_glyco[i] == 'S' || pepseq_glyco[i] == 'T') {
                    // update sequon index and keep looking
                    return i - 2;       // sequon index is 2 behind S/T
                }
                xFlag = false;
            }
            if (asnFlag) {
                // look for X (anything other than pro)
                xFlag = pepseq_glyco[i] != 'P';
            }
            // always reset Asn flag in case of multiple Asn in a row
            asnFlag = pepseq_glyco[i] == 'N';
        }
        // if we reach this point, no sequon was found. This can happen if the sequence ends in NX, where X is an enzyme cut point.
        // Assume the second to last residue is the desired index
        return pepSeq.length() - 2;
    }

    /**
     * Parse the internal list of oxonium ion descriptions to generate the internal database of fragment info
     * by residue type
     * @return map of residue type : list of fragment descriptors parsed from the table
     */
    public static HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> parseOxoniumDatabase(ProbabilityTables probabilityTable) {
        HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> oxoniumDB = new HashMap<>();
        BufferedReader in;
        try {
            // no glycan database provided - fall back to default glycan list in PeakAnnotator
            String defaultDB = "oxonium_ion_list.txt";
            in = new BufferedReader(new InputStreamReader(GlycoAnalysis.class.getResourceAsStream(defaultDB)));
            String line;
            while ((line = in.readLine()) != null) {
                if (line.startsWith("#"))
                    continue;
                String[] splits = line.split("\t");
                GlycanResidue residue = GlycanMasses.glycoNames.get(splits[0].trim().toLowerCase(Locale.ROOT));
                TreeMap<GlycanResidue, Integer> ionComposition = StaticGlycoUtilities.parseGlycanString(splits[1]);
                double massShift = Double.parseDouble(splits[2]);

                // Add to existing list if present or create new list if residue type not seen yet
                if (oxoniumDB.containsKey(residue)) {
                    oxoniumDB.get(residue).add(new GlycanFragmentDescriptor(ionComposition, probabilityTable.rulesByResidue.get(residue), massShift));
                } else {
                    ArrayList<GlycanFragmentDescriptor> residueList = new ArrayList<>();
                    residueList.add(new GlycanFragmentDescriptor(ionComposition, probabilityTable.rulesByResidue.get(residue), massShift));
                    oxoniumDB.put(residue, residueList);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            PTMShepherd.die("IO Exception while reading oxonium database file");
        }
        return oxoniumDB;
    }

}
