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
    File glycoFile;
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
    boolean normYions;
    double defaultMassErrorAbsScore;
    Integer[] glycoIsotopes;
    double glycoPPMtol;
    boolean runGlycanAssignment;

    public GlycoAnalysis(String dsName, boolean runGlycanAssignment, ArrayList<GlycanCandidate> glycoDatabase, ProbabilityTables inputProbabilityTable, boolean normYs, double absMassErrorDefault, Integer[] glycoIsotopes, double glycoPPMtol) {
        this.dsName = dsName;
        this.runGlycanAssignment = runGlycanAssignment;
        this.glycoFile = new File(PTMShepherd.normFName(dsName + ".rawglyco"));
        this.glycanDatabase = glycoDatabase;

        // init with default values, can be changed by params
        this.probabilityTable = inputProbabilityTable;
        this.normYions = normYs;
        this.defaultMassErrorAbsScore = absMassErrorDefault;
        this.glycoPPMtol = glycoPPMtol;
        this.glycoIsotopes = glycoIsotopes;
    }

    public void glycoPSMs(PSMFile pf, HashMap<String, File> mzMappings, ExecutorService executorService, int numThreads) throws Exception {
        //open up output file
        HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
        PrintWriter out = new PrintWriter(new FileWriter(glycoFile));
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
            headbuff.append("\tBest Glycan\tGlycan Score\tGlycan q-value");
        }
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
                futureList.add(executorService.submit(() -> processLinesBlock(cBlock, out)));
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

    public void processLinesBlock(ArrayList<String> cBlock, PrintWriter out) {
        StringBuilder newBlock  = new StringBuilder();
        for (String line : cBlock) {
            newBlock.append(processLine(line)).append("\n");
        }
        printLines(out, newBlock.toString());
    }

    private synchronized void printLines(PrintWriter out, String linesBlock) {
        out.print(linesBlock);
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
            if (desiredRatio > targetDecoyRatio * 10) {
                PTMShepherd.print(("\tNot enough decoys to compute FDR at 0.1 * initial ratio. Check data and parameters. No FDR calculation performed!\n"));
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
            rawGlycoLine[qValCol] = String.format("%s", qval);

            if (!foundThreshold) {
                // still below the threshold: continue checking decoys/targets and appending 'failfdr'
                if (targetDecoyRatio <= desiredRatio) {
                    // stop here, found cutoff
                    foundThreshold = true;
                    PTMShepherd.print(String.format("\tConverged to %.1f pct FDR with %d targets and %d decoys", targetDecoyRatio * 100, targets, decoys));
                }
                if (!rawGlycoLine[bestGlycanCol].contains("FailFDR")) {
                    // only add failFDR annotation once (prevents multiple writes on re-analyses)
                    rawGlycoLine[bestGlycanCol] = "FailFDR_" + rawGlycoLine[bestGlycanCol];
                }
            }

            // update counts
            if (rawGlycoLine[bestGlycanCol].toLowerCase(Locale.ROOT).contains("decoy")) {
                decoys--;
            } else {
                targets--;
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


    public String processLine(String line) {
        StringBuilder sb = new StringBuilder();
        String[] sp = line.split("\\t");
        String seq = sp[pepCol];
        float dmass = Float.parseFloat(sp[deltaCol]);
        float pepMass = Float.parseFloat(sp[pmassCol]);
        String[] smods = sp[modCol].split(",");
        String specName = sp[specCol];

        sb.append(String.format("%s\t%s\t%s\t%.4f\t%.4f", specName, seq, sp[modCol], pepMass, dmass));

        Spectrum spec = mr.getSpectrum(reNormName(specName));
        if (spec == null) {
            this.lineWithoutSpectra.add(reNormName(specName));
            return "ERROR";
        }
        spec.conditionOptNorm(condPeaks, condRatio, false);

        if (runGlycanAssignment) {
            sb.append(assignGlycanToPSM(spec, pepMass, dmass, glycanDatabase, massErrorWidth, meanMassError));
        }
        //System.out.println("got spec");
        double[] capYIonIntensities;
        double[] oxoniumIonIntensities;
        capYIonIntensities = findCapitalYIonMasses(spec, pepMass);
        oxoniumIonIntensities = findOxoniumIonMasses(spec, pepMass);

        for (double capYIonIntensity : capYIonIntensities) sb.append(String.format("\t%.2f", capYIonIntensity));
        for (double oxoniumIonIntensity : oxoniumIonIntensities) sb.append(String.format("\t%.2f", oxoniumIonIntensity));
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
        return sb.toString();
    }

    /**
     * Main glycan assignment method at PSM level. Searches Y/Oxonium ions (and eventually exact mass/isotope) to compare
     * to possible glycan candidates. Goal is to return best glycan candidate and score.
     *
     * @param spec           spectrum being searched
     * @param pepMass        peptide mass (without glycan)
     * @param glycanDatabase possible glycan candidates
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     * @param deltaMass      observed delta mass from PSM
     */
    public String assignGlycanToPSM(Spectrum spec, double pepMass, double deltaMass, ArrayList<GlycanCandidate> glycanDatabase, double massErrorWidth, double meanMassError) {
        // skip non-delta mass PSMs - leave added columns empty
        if (deltaMass < 3.5 && deltaMass > -1.5) {
            StringBuilder sb = new StringBuilder();
            for (int i=0; i < NUM_ADDED_GLYCO_PSM_COLUMNS; i++){
                sb.append("\t");
            }
            return sb.toString();
        }

        // Determine possible glycan candidates from mass
        ArrayList<GlycanCandidate> searchCandidates = getMatchingGlycansByMass(pepMass, deltaMass, glycanDatabase, glycoIsotopes, glycoPPMtol);

        // Search Y and oxonium ions in spectrum for each candidate
        float ppmTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
        for (GlycanCandidate candidate : searchCandidates) {
            for (GlycanFragment yFragment : candidate.Yfragments) {
                yFragment.foundIntensity = spec.findIonNeutral(yFragment.neutralMass + pepMass, ppmTol, spec.charge);  // sum of charge state intensities if >1 found
            }
            for (GlycanFragment oxoniumFragment: candidate.oxoniumFragments) {
                oxoniumFragment.foundIntensity = spec.findIon(oxoniumFragment.neutralMass + AAMasses.protMass, ppmTol);
            }
        }

        // score candidates and save results
        int bestCandidateIndex = 0;
        int nextBestCandidateIndex = 1;

        for (int i = 1; i < searchCandidates.size(); i++) {
            if (i == bestCandidateIndex) {
                continue;
            }
            double comparisonScore = pairwiseCompareGlycans(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), deltaMass, meanMassError);
            if (comparisonScore < 0) {
                // new best candidate - reset best candidate position and update scores at all other positions
                nextBestCandidateIndex = bestCandidateIndex;
                bestCandidateIndex = i;
            }
        }
        // it's possible for the 2nd best score to be worse than those for candidates searched AFTER the last change to best candidate. Compare those against 2nd best score to get correct result
        for (int i = bestCandidateIndex + 1; i < searchCandidates.size(); i++) {
            double comparisonScore2 = pairwiseCompareGlycans(searchCandidates.get(nextBestCandidateIndex), searchCandidates.get(i), deltaMass, meanMassError);
            if (comparisonScore2 < 0) {
                // currently listed 2nd best glycan had worse score than this one: change index to this glycan
                nextBestCandidateIndex = i;
                // don't need to recalculate comparison score because we still have it relative to the best glycan from before
            }
        }

        // compute absolute score for best glycan
        GlycanFragment[] foundYions = mergeYFragments(searchCandidates);
        GlycanFragment[] foundOxoniumIons = mergeOxoniumFragments(searchCandidates);
        double absoluteScore = 0;
        if (searchCandidates.size() > 0) {
            absoluteScore = computeAbsoluteScore(searchCandidates.get(bestCandidateIndex), foundYions, foundOxoniumIons, deltaMass, massErrorWidth, meanMassError);
        }
        // output - best glycan, scores, etc back to PSM table
        String output;
        if (searchCandidates.size() == 0) {
            output = "\tNo Matches\t\t";
        } else {
            output = String.format("\t%s\t%.1f\t", searchCandidates.get(bestCandidateIndex).toString(), absoluteScore);
        }
        return output;
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
                    sumLogRatio += Math.log(fragment1.ruleProbabilities[0]);    // candidate 1 hit - added
                } else {
                    sumLogRatio += Math.log(fragment1.ruleProbabilities[1]);    // candidate 1 miss - negative value added
                }
            }
        }
        for (GlycanFragment fragment2 : glycan2.oxoniumFragments) {
            if (!fragment2.isAllowedFragment(glycan1)) {
                boolean foundInSpectrum = fragment2.foundIntensity > 0;
                if (foundInSpectrum) {
                    sumLogRatio -= Math.log(fragment2.ruleProbabilities[0]);    // candidate 2 hit - subtracted
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
            double massError1 = deltaMass - glycan1.monoisotopicMass - (roundedIso1 * AAMasses.averagineIsotopeMass);
            double massStDevs1 = massError1 - meanMassError;
            double massError2 = deltaMass - glycan2.monoisotopicMass - (roundedIso2 * AAMasses.averagineIsotopeMass);
            double massStDevs2 = massError2 - meanMassError;
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
     * @param yFragments array of possible Y fragments
     * @return absolute score
     */
    public double computeYAbsoluteScoreNormed(GlycanCandidate bestGlycan, GlycanFragment[] yFragments) {
        double sumLogRatio = 0;
        // Y ions - check if allowed for this composition and score if so (ignore if not)
        int hitCount = 0;
        int missCount = 0;
        int disallowedHitCount = 0;
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
        // yFragments contains all fragments FOUND in spectrum. Check for those not belonging to the best glycan
        for (GlycanFragment otherFragment : yFragments) {
            if (!otherFragment.isAllowedFragment(bestGlycan)) {
                // only allow target-target and decoy-decoy matches to affect absolute score
                if ((bestGlycan.isDecoy && otherFragment.isDecoy) || (!bestGlycan.isDecoy && !otherFragment.isDecoy)) {
                    disallowedHitCount++;
                }
            }
        }
        // normalize hit/miss counts and return final probability
        sumLogRatio = Math.sqrt(hitCount) * Math.log(hitProb) + Math.sqrt(missCount) * Math.log(missProb) + Math.sqrt(disallowedHitCount) * Math.log(1 / hitProb);
        return sumLogRatio;
    }

    /**
     * Compute the "absolute" score of the Y ions of the provided glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     * @param bestGlycan glycan candidate to calculate score for
     * @param fragments array of possible fragment ions with intensities from spectrum
     * @return absolute score
     */
    public double computeYAbsoluteScore(GlycanCandidate bestGlycan, GlycanFragment[] fragments){
        double sumLogRatio = 0;
        for (GlycanFragment yFragment : bestGlycan.Yfragments) {
            if (yFragment.foundIntensity > 0) {
                sumLogRatio += Math.log(yFragment.ruleProbabilities[0]);     // found in spectrum - ion supports this glycan
            } else {
                sumLogRatio += Math.log(yFragment.ruleProbabilities[1]);     // not found in spectrum - ion does not support this glycan
            }
        }

        // yFragments contains all fragments FOUND in spectrum. Check for those not belonging to the best glycan
        for (GlycanFragment otherFragment : fragments) {
            if (!otherFragment.isAllowedFragment(bestGlycan)) {
                // only allow target-target and decoy-decoy matches to affect absolute score
                if ((bestGlycan.isDecoy && otherFragment.isDecoy) || (!bestGlycan.isDecoy && !otherFragment.isDecoy)) {
                    sumLogRatio += Math.log(1 / otherFragment.ruleProbabilities[0]);
                }
            }
        }
        return sumLogRatio;
    }
    /**
     * Compute the "absolute" score of the Y ions of the provided glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     * @param bestGlycan glycan candidate to calculate score for
     * @param fragments array of possible fragment ions with intensities from spectrum
     * @return absolute score
     */
    public double computeOxoAbsoluteScore(GlycanCandidate bestGlycan, GlycanFragment[] fragments){
        double sumLogRatio = 0;
        for (GlycanFragment fragment : bestGlycan.oxoniumFragments) {
            if (fragment.foundIntensity > 0) {
                sumLogRatio += Math.log(fragment.ruleProbabilities[0]);     // found in spectrum - ion supports this glycan
            } else {
                sumLogRatio += Math.log(fragment.ruleProbabilities[1]);     // not found in spectrum - ion does not support this glycan
            }
        }

        // fragments contains all fragments of given type FOUND in spectrum. Check for those not belonging to the best glycan
        for (GlycanFragment otherFragment : fragments) {
            if (!otherFragment.isAllowedFragment(bestGlycan)) {
                // only allow target-target and decoy-decoy matches to affect absolute score
                if ((bestGlycan.isDecoy && otherFragment.isDecoy) || (!bestGlycan.isDecoy && !otherFragment.isDecoy)) {
                    sumLogRatio += Math.log(1 / otherFragment.ruleProbabilities[0]);
                }
            }
        }
        return sumLogRatio;
    }
    /**
     * Compute the "absolute" score of this glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     * @param bestGlycan glycan candidate to calculate score for
     * @param yFragments array of possible Y fragments
     * @param oxoFragments array of possible oxonium fragments
     * @param deltaMass spectrum delta mass
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     * @param meanMassError mean mass error of non-delta mass peptides
     * @return absolute score
     */
    public double computeAbsoluteScore(GlycanCandidate bestGlycan, GlycanFragment[] yFragments, GlycanFragment[] oxoFragments, double deltaMass, double massErrorWidth, double meanMassError) {
        double sumLogRatio;
        if (normYions) {
            sumLogRatio = computeYAbsoluteScoreNormed(bestGlycan, yFragments);
        } else {
            sumLogRatio = computeYAbsoluteScore(bestGlycan, yFragments);
        }
        sumLogRatio += computeOxoAbsoluteScore(bestGlycan, oxoFragments);

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

    /**
     * Generate a list of all Y ions with intensity found in spectrum from amongst all candidates in the list.
     * @param searchCandidates list of candidates being searched
     * @return array of fragments
     */
    private GlycanFragment[] mergeYFragments(ArrayList<GlycanCandidate> searchCandidates) {
        ArrayList<GlycanFragment> allFragments = new ArrayList<>();
        HashMap<String, Boolean> fragmentInList = new HashMap<>();
        for (GlycanCandidate candidate : searchCandidates) {
            for (GlycanFragment fragment : candidate.Yfragments) {
                if (!fragmentInList.containsKey(fragment.toHashString())) {
                    if (fragment.foundIntensity > 0) {
                        allFragments.add(fragment);
                        fragmentInList.put(fragment.toHashString(), true);
                    }
                }
            }
        }
        return allFragments.toArray(new GlycanFragment[0]);
    }
    /**
     * Generate a list of all oxonium ions with intensity found in spectrum from amongst all candidates in the list.
     * @param searchCandidates list of candidates being searched
     * @return array of fragments
     */
    private GlycanFragment[] mergeOxoniumFragments(ArrayList<GlycanCandidate> searchCandidates) {
        ArrayList<GlycanFragment> allFragments = new ArrayList<>();
        HashMap<String, Boolean> fragmentInList = new HashMap<>();
        for (GlycanCandidate candidate : searchCandidates) {
            for (GlycanFragment fragment : candidate.oxoniumFragments) {
                if (!fragmentInList.containsKey(fragment.toHashString())) {
                    if (fragment.foundIntensity > 0) {
                        allFragments.add(fragment);
                        fragmentInList.put(fragment.toHashString(), true);
                    }
                }
            }
        }
        return allFragments.toArray(new GlycanFragment[0]);
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
                capYIonIntensities[i] /= spec.findBasePeakInt();
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
                oxoniumIonIntensities[i] /= spec.findBasePeakInt();
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
        int glycanAssignmentAddedLines = runGlycanAssignment ? NUM_ADDED_GLYCO_PSM_COLUMNS : 0;     // number of added lines to rawglyco file - 0 if not running glycan assignment
        while ((cline = in.readLine()) != null) {
            if (cline.equals("COMPLETE"))
                break;
            if (cline.startsWith("Spectrum"))
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
                TreeMap<GlycanResidue, Integer> ionComposition = PTMShepherd.parseGlycanString(splits[1]);
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
