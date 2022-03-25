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

package edu.umich.andykong.ptmshepherd.glyco;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import static edu.umich.andykong.ptmshepherd.PTMShepherd.reNormName;

public class GlycoAnalysis {
    String dsName;
    File glycoFile;                     // .rawglyco file
    MXMLReader mr;
    ArrayList<String> lineWithoutSpectra = new ArrayList<>();
    int totalLines;
    HashMap<String, MXMLReader> multiMr;
    float ppmTol, peakTol;
    int condPeaks;
    int specCol, pepCol, modpepCol, chargecol, deltaCol, rtCol, intCol, pmassCol, modCol;
    double condRatio;
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
    public static final double DEFAULT_GLYCO_PPM_TOL = 50;
    public static final double DEFAULT_GLYCO_FDR = 0.01;
    public static final int DEFAULT_GLYCO_DECOY_TYPE = 1;
    public static final double DEFAULT_GLYCO_ABS_SCORE_BASE = 5;
    public boolean useFragmentSpecificProbs;
    public HashMap<Integer, HashMap<String, Integer>> glycanMassBinMap;
    public static final int MIN_GLYCO_PSMS_FOR_BOOTSTRAP = 10;
    public static final int MAX_PROB_RATIO_ABSOLUTE = 10;
    public double finalGlycoFDR;

    // Default constructor
    public GlycoAnalysis(String dsName, ArrayList<GlycanCandidate> glycoDatabase, ProbabilityTables inputProbabilityTable, boolean normYs, double absMassErrorDefault, Integer[] glycoIsotopes, double glycoPPMtol) {
        this.dsName = dsName;
        this.glycoFile = new File(PTMShepherd.normFName(dsName + PTMShepherd.rawGlycoName));
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
        PrintWriter glycoOut = new PrintWriter(new FileWriter(glycoFile));
        ArrayList<String> linesWithoutSpectra = new ArrayList<>();

        //get necessary params
        ppmTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
        condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
        condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));

        //write header
        glycoOut.println(String.format("%s\t%s\t%s\t%s\t%s", "Spectrum", "Peptide", "Mods", "Pep Mass", "Mass Shift") + "\tBest Glycan\tGlycan Score\tGlycan q-value\tBest Target Glycan\tBest Target Score" + "\tFragments:");

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
                futureList.add(executorService.submit(() -> processLinesBlock(cBlock, glycoOut)));
            }
            /* Wait for all processes to finish */
            for (Future future : futureList)
                future.get();

            long t3 = System.currentTimeMillis();
            PTMShepherd.print(String.format("\t%s - %d (%d ms, %d ms)", cf, clines.size(), t2 - t1, t3 - t2));
        }
        glycoOut.close();

        if (!linesWithoutSpectra.isEmpty()) {
            PTMShepherd.print(String.format("Could not find %d/%d (%.1f%%) spectra.\n", linesWithoutSpectra.size(), this.totalLines,
                    100.0*((double)linesWithoutSpectra.size()/this.totalLines)));
            int previewSize = Math.min(linesWithoutSpectra.size(), 5);
            PTMShepherd.print(String.format("Showing first %d of %d spectra IDs that could not be found: \n\t%s\n", previewSize, linesWithoutSpectra.size(),
                    String.join("\n\t", linesWithoutSpectra.subList(0, previewSize))));
        }
    }

    public void processLinesBlock(ArrayList<String> cBlock, PrintWriter fragmentOutWriter) {
        StringBuilder fragmentBlock = new StringBuilder();
        for (String line : cBlock) {
            GlycanAssignmentResult glycoResult = processLine(line);
            fragmentBlock.append(glycoResult.printGlycoFragmentInfo());
        }
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
        BufferedReader in = new BufferedReader(new FileReader(glycoFile), 1 << 22);
        String[] headerSplits = in.readLine().split("\t");
        int glycanCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Best Glycan");
        int qValCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Glycan q-value");
        int deltaMassCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Mass Shift");
        int fragmentStartCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Fragments:");

        // read all glycan info in
        String currentLine;
        while ((currentLine = in.readLine()) != null) {
            String[] splits = currentLine.split("\t", 0);       // limit 0 to discard extra empty cells if present
            // only read lines with glycan info (after column 5)
            if (splits.length > 5) {
                String glycanString = splits[glycanCol];
                boolean isDecoy = Double.parseDouble(splits[qValCol]) > finalGlycoFDR;
                String[] fragmentInfo = splits.length >= fragmentStartCol ? Arrays.copyOfRange(splits, fragmentStartCol, splits.length) : new String[]{};
                GlycanCandidate fragmentInfoContainer = new GlycanCandidate(glycanString, fragmentInfo);
                String glycanHash = fragmentInfoContainer.toString();
                // only include targets in fragment info
                if (!isDecoy) {
                    if (glycanInputMap.containsKey(glycanHash)) {
                        glycanInputMap.get(glycanHash).add(fragmentInfoContainer);
                    } else {
                        ArrayList<GlycanCandidate> newList = new ArrayList<>();
                        newList.add(fragmentInfoContainer);
                        glycanInputMap.put(glycanHash, newList);
                    }
                }
                // add to delta mass map for calculating glycan prevalence priors (targets and decoys)
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
                for (String fragmentHash : inputGlycan.Yfragments.keySet()) {
                    int count = YCounts.getOrDefault(fragmentHash, 0);
                    count++;
                    YCounts.put(fragmentHash, count);
                }
                for (String fragmentHash : inputGlycan.oxoniumFragments.keySet()) {
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
     * @param glycoFDR: desired FDR (typically 0.01 = 1%)
     */
    public void computeGlycanFDR(double glycoFDR) throws IOException {
        finalGlycoFDR = glycoFDR;
        BufferedReader in = new BufferedReader(new FileReader(glycoFile), 1 << 22);

        // read rawglyco file into map of spectrum index: full line (string)
        LinkedHashMap<String, String[]> glyLines = new LinkedHashMap<>();   // linkedHashMap to preserve spectrum order
        String cgline;

        // detect headers
        String[] headerSplits = in.readLine().split("\t");
        int gSpecCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Spectrum");
        int absScoreCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Glycan Score");
        int bestGlycanCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Best Glycan");
        int qValCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Glycan q-value");
        int bestNextScoreCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Best Target Score");

        if (absScoreCol <= 0 || bestGlycanCol <= 0 || qValCol <= 0) {
            PTMShepherd.print(String.format("Warning: rawglyco file headers not found! FDR calculation may fail for file %s\n", glycoFile));
        }

        // read file, accumulating scores
        HashMap<String, Double> scoreMap = new HashMap<>();
        ArrayList<GlycoScore> scoreDistribution = new ArrayList<>();
        int targets = 0;
        int decoys = 0;
        int scoreDistTargets = 0;
        int scoreDistDecoys = 0;
        while ((cgline = in.readLine()) != null) {
            if (cgline.equals("COMPLETE")) {
                break;
            }
            if (cgline.startsWith("ERROR"))
                continue;
            String[] splits = cgline.split("\t", -1);
            // skip non-glyco columns
            if (splits.length < bestGlycanCol + 1)
                continue;
            if (splits[bestGlycanCol].matches("ERROR"))
                continue;
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
                    // Best candidate was a decoy. Store both target and decoy scores, remembering that the decoy was top
                    scoreDistribution.add(new GlycoScore (Double.parseDouble(splits[absScoreCol]), true, spectrumID, true));
                    scoreDistribution.add(new GlycoScore (Double.parseDouble(splits[bestNextScoreCol]), false, spectrumID, false));
                } else {
                    targets++;
                    // Best candidate was a target. Store both target and decoy scores, remembering that the target was top
                    scoreDistribution.add(new GlycoScore (Double.parseDouble(splits[absScoreCol]), false, spectrumID, true));
                    scoreDistribution.add(new GlycoScore (Double.parseDouble(splits[bestNextScoreCol]), true, spectrumID, false));
                }
                scoreDistDecoys++;
                scoreDistTargets++;
                scoreMap.put(spectrumID, Double.parseDouble(splits[absScoreCol]));
            }
        }
        in.close();

        PTMShepherd.print("Calculating Glycan FDR");
        // sort scoreMap in order of ascending score
        scoreDistribution.sort(GlycoScore::compareTo);

        double targetDecoyRatio;
        double currentMinQ = 1;
        double scoreThreshold = -10000;
        HashMap<String, Double> qValMap = new HashMap<>();        // hashmap for saving q-values. NOTE: there are 2 score objects per spectrum, but only the one from the top candidate is used
        boolean foundScoreThresh = false;
        for (GlycoScore scoreObj : scoreDistribution) {
            if (!scoreObj.isDecoy) {
                scoreDistTargets--;
            } else {
                scoreDistDecoys--;
            }
            // compute TD ratio
            targetDecoyRatio = (2 * scoreDistDecoys) / (double) (scoreDistDecoys + scoreDistTargets);
            if (scoreDistDecoys > scoreDistTargets) {
                targetDecoyRatio = 1.0;     // cap FDR at 1
            } else if (scoreDistTargets == 0) {
                targetDecoyRatio = 0.0;     // min FDR = 0. Using else-if with the above block so that if decoys are nonzero with 0 targets, FDR = 1
            }

            // compute q-value and save for later
            double qval = Math.min(targetDecoyRatio, currentMinQ);
            if (qval < currentMinQ){
                currentMinQ = qval;
            }
            if (scoreObj.isFromTopCandidate) {
                // save q-val only for top candidates
                qValMap.put(scoreObj.spectrumID, qval);
            }

            // check for the score threshold that gives the requested FDR
            if (!foundScoreThresh) {
                if (targetDecoyRatio <= finalGlycoFDR) {
                    // stop here, found cutoff
                    scoreThreshold = scoreObj.score;
                    PTMShepherd.print(String.format("\tFound score threshold of %.2f for %.1f%% FDR with %d targets and %d decoys from non-competitive analysis (%d total inputs)", scoreThreshold, targetDecoyRatio * 100, scoreDistTargets, scoreDistDecoys, scoreDistribution.size()));
                    foundScoreThresh = true;
                }
            }
        }


        // sort scoreMap in order of ascending score
        List<Map.Entry<String, Double>> entries = new ArrayList<>(scoreMap.entrySet());
        entries.sort(Map.Entry.comparingByValue());
        Map<String, Double> sortedScoreMap = new LinkedHashMap<>();
        for (Map.Entry<String, Double> entry : entries) {
            sortedScoreMap.put(entry.getKey(), entry.getValue());
        }

        /*
            Find threshold at which target/decoy ratio hits desired value and update rawglyco lines
            NOTE: calculation is FDR <= desired ratio, so target/decoy counts updated AFTER the current PSM. This means the last
            decoy passes FDR, and the last target has the correct FDR instead of 0/0
            NOTE2: q-value is set to min(current FDR, FDR of all PSMs with lower score) to make it a step down rather than sawtooth shape
         */
        boolean foundThreshold = false;
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

            // Write q-value to output, and write q=1 for decoys
            if (rawGlycoLine[bestGlycanCol].toLowerCase(Locale.ROOT).contains("decoy")) {
                rawGlycoLine[qValCol] = "1";
            } else {
                rawGlycoLine[qValCol] = String.format("%s", qValMap.get(rawGlycoLine[gSpecCol]));
            }
            if (!foundThreshold) {
                // still below the threshold: continue checking decoys/targets and appending 'failfdr'
                if (scoreEntry.getValue() >= scoreThreshold) {
                    // stop here, found cutoff
                    foundThreshold = true;
                    PTMShepherd.print(String.format("\tUsed score threshold to obtain %.1f%% competitive FDR with %d targets and %d decoys (%d total inputs)", targetDecoyRatio * 100, targets, decoys, sortedScoreMap.size()));
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

    public void computeGlycanFDROld(double glycoFDR) throws IOException {
        finalGlycoFDR = glycoFDR;
        BufferedReader in = new BufferedReader(new FileReader(glycoFile), 1 << 22);

        // read rawglyco file into map of spectrum index: full line (string)
        LinkedHashMap<String, String[]> glyLines = new LinkedHashMap<>();   // linkedHashMap to preserve spectrum order
        String cgline;

        // detect headers
        String[] headerSplits = in.readLine().split("\t");
        int gSpecCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Spectrum");
        int absScoreCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Glycan Score");
        int bestGlycanCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Best Glycan");
        int qValCol = StaticGlycoUtilities.getHeaderColIndex(headerSplits, "Glycan q-value");

        if (absScoreCol <= 0 || bestGlycanCol <= 0 || qValCol <= 0) {
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
            // skip non-glyco columns
            if (splits.length < bestGlycanCol + 1)
                continue;
            if (splits[bestGlycanCol].matches("ERROR"))
                continue;
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
        if (targetDecoyRatio < finalGlycoFDR) {
            // not enough decoys to compute FDR - already above desired ratio. Do not update table
            PTMShepherd.print(String.format("\tNot enough decoys to compute FDR at %.1f%%, started at %.2f%%", finalGlycoFDR * 100, targetDecoyRatio * 100));
            // only missed by a little, try reducing desired FDR to accomodate
            finalGlycoFDR = targetDecoyRatio - (targetDecoyRatio * 0.1);
            PTMShepherd.print(String.format("\tFDR reduced to %.2f pct due to limited decoys", finalGlycoFDR * 100));
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
                if (targetDecoyRatio <= finalGlycoFDR) {
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

    /**
     * Run glycan assignment for a single PSM line from the PSM table
     * @param line String of a single PSM line
     * @return result container
     */
    public GlycanAssignmentResult processLine(String line) {
        // get basic info
        String[] sp = line.split("\\t");
        String seq = sp[pepCol];
        float dmass = Float.parseFloat(sp[deltaCol]);
        float pepMass = Float.parseFloat(sp[pmassCol]);
        String specName = sp[specCol];
        GlycanAssignmentResult glycoResult = new GlycanAssignmentResult(seq, dmass, pepMass, sp[modCol], specName);

        // read spectrum and condition
        Spectrum spec = mr.getSpectrum(reNormName(specName));
        if (spec == null) {
            this.lineWithoutSpectra.add(reNormName(specName));
            glycoResult.glycanAssignmentString = "ERROR";
            return glycoResult;
        }
        spec.conditionOptNorm(condPeaks, condRatio, false);

        // do glycan assignment
        glycoResult = assignGlycanToPSM(spec, glycoResult, glycanDatabase, massErrorWidth, meanMassError);
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
                for (GlycanFragment yFragment : candidate.Yfragments.values()) {
                    yFragment.foundIntensity = spec.findIonNeutral(yFragment.neutralMass + glycoResult.pepMass, ppmTol, spec.charge) / spec.basePeakInt;  // sum of charge state intensities if >1 found
                }
                for (GlycanFragment oxoniumFragment : candidate.oxoniumFragments.values()) {
                    // save oxonium ion intensity relative to base peak
                    oxoniumFragment.foundIntensity = spec.findIon(oxoniumFragment.neutralMass + AAMasses.protMass, ppmTol) / spec.basePeakInt;
                }
            }

            // score candidates and save results
            int bestCandidateIndex = 0;
            double[] scoresVsBestCandidate = new double[searchCandidates.size()];
            for (int i = 0; i < searchCandidates.size(); i++) {
                if (i == bestCandidateIndex) {
                    continue;
                }
                double comparisonScore;
                if (useFragmentSpecificProbs) {
                    comparisonScore = pairwiseCompareDynamic(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass);
                } else {
                    comparisonScore = pairwiseCompareStatic(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass, meanMassError);
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
                    scoresVsBestCandidate[i] = pairwiseCompareDynamic(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass);
                } else {
                    scoresVsBestCandidate[i] = pairwiseCompareStatic(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass, meanMassError);
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
            double absoluteScore;
            if (useFragmentSpecificProbs) {
                absoluteScore = computeAbsoluteScoreDynamic(searchCandidates.get(bestCandidateIndex), glycoResult.deltaMass, massErrorWidth, meanMassError);
            } else {
                absoluteScore = computeAbsoluteScore(searchCandidates.get(bestCandidateIndex), glycoResult.deltaMass, massErrorWidth, meanMassError);
            }
            glycoResult.bestCandidate = searchCandidates.get(bestCandidateIndex);
            glycoResult.glycanScore = absoluteScore;

            // output - best glycan, scores, etc back to PSM table
            // write top glycan info
            output = String.format("\t%s\t%.4f\t", searchCandidates.get(bestCandidateIndex).toString(), absoluteScore);

            // if top glycan is a decoy, also write best target and best target score to subsequent columns
            boolean bestWasTarget = !searchCandidates.get(bestCandidateIndex).isDecoy;
            glycoResult.isDecoyGlycan = !bestWasTarget;
            output = getNextGlycanScore(bestWasTarget, glycoResult, massErrorWidth, meanMassError, searchCandidates, output, sortedIndicesOfBestScores);
        } else {
            output = "\tNo Matches\t\t\t\t";
        }
        glycoResult.glycanAssignmentString = output;
        return glycoResult;
    }

    /**
     * Helper method to get output for next best glycan (i.e., best decoy if the top hit is a target and vice versa).
     * Returns formatted output string handling various cases
     * @param glycoResult
     * @param massErrorWidth
     * @param meanMassError
     * @param searchCandidates
     * @param output
     * @param sortedIndicesOfBestScores
     * @return
     */
    private String getNextGlycanScore(boolean bestWasTarget, GlycanAssignmentResult glycoResult, double massErrorWidth, double meanMassError, ArrayList<GlycanCandidate> searchCandidates, String output, int[] sortedIndicesOfBestScores) {
        boolean foundNext = false;
        for (int bestIndex : sortedIndicesOfBestScores) {
            GlycanCandidate nextCandidate = searchCandidates.get(bestIndex);
            // if Best hit was target, look for Next to be decoy and vice versa
            if (nextCandidate.isDecoy == bestWasTarget) {
                // add the next hit's information
                double bestNextScore;
                if (useFragmentSpecificProbs) {
                    bestNextScore = computeAbsoluteScoreDynamic(nextCandidate, glycoResult.deltaMass, massErrorWidth, meanMassError);
                } else {
                    bestNextScore = computeAbsoluteScore(nextCandidate, glycoResult.deltaMass, massErrorWidth, meanMassError);
                }
                output = String.format("%s\t%s\t%.4f", output, nextCandidate.toString(), bestNextScore);
                foundNext = true;
                if (bestWasTarget) {
                    glycoResult.bestDecoy = nextCandidate;
                    glycoResult.bestDecoyScore = bestNextScore;
                } else {
                    glycoResult.bestTarget = nextCandidate;
                    glycoResult.bestTargetScore = bestNextScore;
                }
                break;
            }
        }
        if (!foundNext) {
            // no target candidates found (only matched to a decoy)
            if (bestWasTarget) {
                output = output + "\tNo decoy matches\t";
            } else {
                output = output + "\tNo target matches\t";
            }
        }
        return output;
    }

    public double pairwiseCompareDynamic(GlycanCandidate glycan1, GlycanCandidate glycan2, double deltaMass) {
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
                glyc1Count = glyc1Count + glycanCountMap.getOrDefault(glycan1.toString(), 0);
                glyc2Count = glyc2Count + glycanCountMap.getOrDefault(glycan2.toString(), 0);
                for (int glycanCount : glycanCountMap.values()) {
                    totalGlycCount += glycanCount;
                }
            }
        }
        double propGlycan1 = glyc1Count / (double) totalGlycCount;
        double propGlycan2 = glyc2Count / (double) totalGlycCount;
        if (propGlycan1 == 0 || propGlycan2 == 0) {
            // ignore glycan proportion data if either is 0. todo: Could instead set some minimum value (empirically determined)?
            propGlycan1 = 1;
            propGlycan2 = 1;
        }

        // calculate fragment-specific prob estimates based on observed fragment ions
        double sumLogRatio = 0;
        for (GlycanFragment fragment1 : glycan1.Yfragments.values()) {
            double probRatio;
            if (fragment1.isAllowedFragment(glycan2)) {
                GlycanFragment fragment2 = glycan2.Yfragments.get(fragment1.toHashString());
                probRatio = computeFragmentPairwiseScore(fragment1, fragment2, propGlycan1, propGlycan2);
            } else {
                // fragment only possible for glycan 1, use glycan 1 only estimate
                probRatio = computeFragmentAbsoluteScore(fragment1, propGlycan1);
            }
            sumLogRatio += Math.log(probRatio);
        }
        // glycan 2 fragments - unique fragments get scored the same way as in the absolute method, and subtracted since they support glycan 2 not 1
        for (GlycanFragment fragment : glycan2.Yfragments.values()) {
            double probRatio = computeFragmentAbsoluteScore(fragment, propGlycan2);
            sumLogRatio -= Math.log(probRatio);
        }

        // oxonium ions
        sumLogRatio += pairwiseCompareOxoDynamic(glycan1, glycan2, propGlycan1, propGlycan2);

        // mass and isotope error
        sumLogRatio += determineIsotopeAndMassErrorProbs(glycan1, glycan2, deltaMass, meanMassError);

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
    public double pairwiseCompareOxoDynamic(GlycanCandidate glycan1, GlycanCandidate glycan2, double propGlycan1, double propGlycan2) {
        double sumLogRatio = 0;

        for (GlycanFragment fragment1 : glycan1.oxoniumFragments.values()) {
            double probRatio;
            if (fragment1.isAllowedFragment(glycan2)) {
                GlycanFragment fragment2 = glycan2.oxoniumFragments.get(fragment1.toHashString());
                probRatio = computeFragmentPairwiseScore(fragment1, fragment2, propGlycan1, propGlycan2);
            } else {
                // fragment only possible for glycan 1, use glycan 1 only (absolute) estimate
                probRatio = computeFragmentAbsoluteScore(fragment1, propGlycan1);
            }
            sumLogRatio += Math.log(probRatio);
        }
        // glycan 2 fragments - unique fragments get scored the same way as in the absolute method
        for (GlycanFragment fragment : glycan2.oxoniumFragments.values()) {
            double probRatio = computeFragmentAbsoluteScore(fragment, propGlycan2);
            sumLogRatio -= Math.log(probRatio);
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
    public double computeAbsoluteScoreDynamic(GlycanCandidate bestGlycan, double deltaMass, double massErrorWidth, double meanMassError) {
        double sumLogRatio = 0;
        // determine the overall likelihood priors of these glycans given the observed delta mass
        int glycCount = 0;
        int totalGlycCount = 0;
        HashMap<String, Integer> emptyMap = new HashMap<>();
        // count glycans in this delta mass bin and nearby allowed bins
        for (int isotope : glycoIsotopes) {
            int massBin = (int) Math.floor(deltaMass + isotope);
            HashMap<String, Integer> glycanCountMap = glycanMassBinMap.getOrDefault(massBin, emptyMap);
            if (glycanCountMap.size() > 0) {
                // count instances of glycan 1, glycan 2, and all glycans
                glycCount = glycCount + glycanCountMap.getOrDefault(bestGlycan.toString(), 0);
                for (int glycanCount : glycanCountMap.values()) {
                    totalGlycCount += glycanCount;
                }
            }
        }
        // ignore glycan proportion data if it's 0. todo: Could instead set some minimum value (empirically determined)?
        double propGlycan = glycCount == 0 ? 1 : glycCount / (double) totalGlycCount;

        // Y ions
        for (GlycanFragment fragment : bestGlycan.Yfragments.values()) {
            double probRatio = computeFragmentAbsoluteScore(fragment, propGlycan);
            sumLogRatio += Math.log(probRatio);
        }
        // oxonium ions
        for (GlycanFragment fragment : bestGlycan.oxoniumFragments.values()) {
            double probRatio = computeFragmentAbsoluteScore(fragment, propGlycan);
            sumLogRatio += Math.log(probRatio);
        }

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
     * Helper for computing absolute score of Y or oxonium fragment lists
     * @param fragment1 fragment from glycan 1
     * @param fragment2 fragment from glycan 2
     * @param propGlycan1 proportion of this glycan in the mass bin
     * @param propGlycan2 proportion of this glycan in the mass bin
     * @return ratio of fragment probs
     */
    public double computeFragmentPairwiseScore(GlycanFragment fragment1, GlycanFragment fragment2, double propGlycan1, double propGlycan2) {
        double probRatio;
        if (fragment1.foundIntensity > 0) {
            // "hit": fragment found in spectrum. Compute prob of glycans given the presence of this ion
            probRatio = fragment1.propensity * propGlycan1 / (fragment2.propensity * propGlycan2);
        } else {
            // "miss": fragment not found. Compute prob of glycans given absence of this ion. Miss propensity = 1 - hit propensity
            probRatio = (1 - fragment1.propensity) * propGlycan1 / ((1 - fragment2.propensity) * propGlycan2);
        }
        return probRatio;
    }

    /**
     * Helper for computing absolute score of Y or oxonium fragment lists
     * @param fragment fragment to consider
     * @param propGlycan proportion of this glycan in the mass bin
     * @return sum log ratio of fragment probs
     */
    public double computeFragmentAbsoluteScore(GlycanFragment fragment, double propGlycan) {
        double probRatio;
        if (fragment.foundIntensity > 0) {
            // "hit": fragment found in spectrum. Prob ratio is inverse of the miss probability in the absence of another glycan to compare against
            probRatio = 1 / ((1 - fragment.propensity) * propGlycan);
            // cap probRatio to prevent high-propensity fragments from going to infinity
            if (probRatio > MAX_PROB_RATIO_ABSOLUTE) {
                probRatio = MAX_PROB_RATIO_ABSOLUTE;
            }
        } else {
            // "miss": fragment not found. Compute prob of glycan given absence of this ion. Miss propensity = 1 - hit propensity
            probRatio = (1 - fragment.propensity) * propGlycan;
        }
        return probRatio;
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
    public double pairwiseCompareStatic(GlycanCandidate glycan1, GlycanCandidate glycan2, double deltaMass, double meanMassError) {
        double sumLogRatio = 0;
        // Y ions
        sumLogRatio += pairwiseCompareYstatic(glycan1, glycan2, normYions);

        // oxonium ions
        sumLogRatio += pairwiseCompareOxoStatic(glycan1, glycan2);

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
    public double pairwiseCompareYstatic(GlycanCandidate glycan1, GlycanCandidate glycan2, boolean normYions) {
        int cand1Misses = 0;
        int cand2Misses = 0;
        int cand1Hits = 0;
        int cand2Hits = 0;
        double sumLogRatio = 0;

        // Loop over each candidate's fragments, scoring unique (i.e., not in the other candidate) fragments as hit/miss if found/not in spectrum
        for (GlycanFragment fragment1 : glycan1.Yfragments.values()) {
            if (!fragment1.isAllowedFragment(glycan2)) {
                boolean foundInSpectrum = fragment1.foundIntensity > 0;
                if (foundInSpectrum) {
                    cand1Hits++;
                } else {
                    cand1Misses++;
                }
            }
        }
        for (GlycanFragment fragment2 : glycan2.Yfragments.values()) {
            if (!fragment2.isAllowedFragment(glycan1)) {
                boolean foundInSpectrum = fragment2.foundIntensity > 0;
                if (foundInSpectrum) {
                    cand2Hits++;
                } else {
                    cand2Misses++;
                }
            }
        }

        // get probability ratios from one fragment (same probs used for all fragments in this method)
        double cand1MissProb = glycan1.Yfragments.values().iterator().next().ruleProbabilities[1];
        double cand2MissProb = glycan2.Yfragments.values().iterator().next().ruleProbabilities[1];
        double cand1HitProb = glycan1.Yfragments.values().iterator().next().ruleProbabilities[0];
        double cand2HitProb = glycan2.Yfragments.values().iterator().next().ruleProbabilities[0];
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
    public double pairwiseCompareOxoStatic(GlycanCandidate glycan1, GlycanCandidate glycan2) {
        double sumLogRatio = 0;

        // Loop over each candidate's fragments, scoring unique (i.e., not in the other candidate) fragments as hit/miss if found/not in spectrum
        for (GlycanFragment fragment1 : glycan1.oxoniumFragments.values()) {
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
        for (GlycanFragment fragment2 : glycan2.oxoniumFragments.values()) {
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
        double hitProb = bestGlycan.Yfragments.values().iterator().next().ruleProbabilities[0];
        double missProb = bestGlycan.Yfragments.values().iterator().next().ruleProbabilities[1];

        for (GlycanFragment yFragment : bestGlycan.Yfragments.values()) {
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
        for (GlycanFragment yFragment : bestGlycan.Yfragments.values()) {
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
        for (GlycanFragment fragment : bestGlycan.oxoniumFragments.values()) {
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


    public boolean isGlycoComplete() throws Exception {
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


    public void completeGlyco() throws Exception {
        PrintWriter out = new PrintWriter(new FileWriter(glycoFile,true));
        out.println("COMPLETE");
        out.close();
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

// Container for holding results for FDR score threshold calculation, compares on score
class GlycoScore implements Comparable<GlycoScore> {
    public final double score;
    public final boolean isDecoy;
    public String spectrumID;
    public boolean isFromTopCandidate;

    public GlycoScore(double score, boolean isDecoy, String spectrumID, boolean isFromTopCandidate) {
        this.score = score;
        this.isDecoy = isDecoy;
        this.spectrumID = spectrumID;
        this.isFromTopCandidate = isFromTopCandidate;
    }

    @Override
    public int compareTo(GlycoScore o) {
        return Double.compare(score, o.score);
    }
}
