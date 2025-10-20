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

import edu.umich.andykong.ptmshepherd.PSM;
import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;
import ionquant.api.Entry;
import ionquant.api.IonQuantAPI;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.hipparchus.stat.descriptive.rank.Median;
import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanCandidate;
import umich.ms.glyco.GlycanFragment;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

public class GlycoAnalysis {
    String dsName;
    File glycoFile;                     // .rawglyco file
    MXMLReader mr;
    ArrayList<String> lineWithoutSpectra = new ArrayList<>();
    int totalLines;
    float ppmTol;
    int condPeaks;
    double condRatio;
    ArrayList<GlycanCandidate> glycanDatabase;
    Double meanMassError;
    double massErrorWidth;
    public static final double DEFAULT_GLYCO_PPM_TOL = 30;
    public static final double DEFAULT_GLYCO_FDR = 0.01;
    public static final int DEFAULT_GLYCO_DECOY_TYPE = 1;
    public static final double DEFAULT_GLYCO_ABS_SCORE_BASE = 5;
    public static final double DEFAULT_MASS_PROB_SCALING = 1;
    public static final String GLYCAN_COMP_COL_NAME = "Total Glycan Composition";
    public boolean useFragmentSpecificProbs;
    public HashMap<Integer, HashMap<String, Integer>> glycanMassBinMap;
    public static final int MIN_GLYCO_PSMS_FOR_BOOTSTRAP = 10;      // todo: param
    public double defaultPropensity;
    public static final double DEFAULT_GLYCO_PROPENSITY = 0.1;      // todo: param?
    private final GlycoParams glycoParams;
    public final ArrayList<GlycanAssignmentResult> allResults;
    private final String ldaHeader;
    private static IonQuantAPI api;
    private final boolean isFirstPass;

    // Default constructor
    public GlycoAnalysis(String dsName, ArrayList<GlycanCandidate> glycoDatabase, GlycoParams glycoParams, boolean isFirstPass) {
        this.dsName = dsName;
        this.isFirstPass = isFirstPass;
        String firstPassName = isFirstPass ? PTMShepherd.rawGlycoFirstPass : "";
        this.glycoFile = new File(PTMShepherd.normFName(dsName + firstPassName + PTMShepherd.rawGlycoName));
        this.glycanDatabase = glycoDatabase;
        this.glycoParams = glycoParams;
        this.useFragmentSpecificProbs = false;
        this.glycanMassBinMap = new HashMap<>();
        this.allResults = new ArrayList<>();
        ldaHeader = glycoParams.glycoLDA ? glycoParams.generateLDAheader() : "\t";
    }

    public void glycoPSMs(PSMFile psmFile,
                          HashMap<String, File> mzMappings,
                          HashMap<String, File> originalMzMappings,
                          ExecutorService executorService) {

        //open up output file
        HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
        ArrayList<String> linesWithoutSpectra = null;
        try {
            PrintWriter glycoOut = new PrintWriter(new FileWriter(glycoFile));
            linesWithoutSpectra = new ArrayList<>();

            //get necessary params
            ppmTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
            condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
            condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));

            //write header
            glycoOut.println(String.format("%s\t%s\t%s\t%s\t%s", "Spectrum", "Peptide", "Mods", "Pep Mass", "Mass Shift") + String.format("\t%s\tGlycan Score\tGlycan q-value\tBest Decoy Glycan\tBest Decoy Score", GLYCAN_COMP_COL_NAME) + ldaHeader + "\tFragments:");

            //map PSMs to file
            SiteLocalization.initSpectrumMappings(psmFile, mappings);

            /* Loop through spectral files -> indexed lines in PSM -> process each line */
            for (String mzFileName : mappings.keySet()) { //for file in relevant spectral files
                long t1 = System.currentTimeMillis();
                //System.out.println(cf);
                mr = new MXMLReader(mzMappings.get(mzFileName), glycoParams.numThreads);
                mr.readFully();
                api = indexBuilder(String.valueOf(originalMzMappings.get(mzFileName)), glycoParams);    // IonQuant API init for KL scoring

                long t2 = System.currentTimeMillis();
                ArrayList<Integer> clines = mappings.get(mzFileName); //lines corr to curr spec file

                if (meanMassError == null) {
                    getMassErrorsFirstPass(psmFile, clines);
                }

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
                    ArrayList<PSM> cBlock = new ArrayList<>();
                    for (int j = startInd; j < endInd; j++)
                        cBlock.add(psmFile.psms.get(clines.get(j)));
                    futureList.add(executorService.submit(() -> processLinesBlock(cBlock, glycoOut)));
                }
                /* Wait for all processes to finish */
                try {
                    for (Future future : futureList)
                        future.get();
                } catch (InterruptedException | ExecutionException e) {
                    e.printStackTrace();
                    PTMShepherd.die("Error in parallel processing glyco PSMs");
                }

                long t3 = System.currentTimeMillis();
                PTMShepherd.print(String.format("\t%s - %d (%d ms, %d ms)", mzFileName, clines.size(), t2 - t1, t3 - t2));
            }
            glycoOut.close();
        } catch (IOException e) {
            PTMShepherd.die("Error writing to glyco file " + glycoFile.getAbsolutePath() + "\n" + e.getMessage());
        }
        if (!linesWithoutSpectra.isEmpty()) {
            PTMShepherd.print(String.format("\tCould not find %d/%d (%.1f%%) spectra.\n", linesWithoutSpectra.size(), this.totalLines,
                    100.0 * ((double) linesWithoutSpectra.size() / this.totalLines)));
            int previewSize = Math.min(linesWithoutSpectra.size(), 5);
            PTMShepherd.print(String.format("\tShowing first %d of %d spectra IDs that could not be found: \n\t%s\n", previewSize, linesWithoutSpectra.size(),
                    String.join("\n\t\t", linesWithoutSpectra.subList(0, previewSize))));
        }

        // save results
        for (PSM psm: psmFile.psms) {
            allResults.add(psm.glycanAssignmentResult);
        }
    }

    /**
     * Initializes IonQuant and builds the feature index for a mass spectrometry file.
     * Only needed if using KL scoring.
     *
     * @param filePath Path to the MS data file
     * @param params Parameters controlling the feature detection
     */
    public static IonQuantAPI indexBuilder(String filePath, GlycoParams params) {
        if (params.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.kl)) {
            PTMShepherd.print("\tBuilding IonQuant index for " + filePath);
            api = new IonQuantAPI(
                    filePath,
                    params.numThreads,
                    (float) params.glycoPPMtol,
                    params.rtTol,
                    params.imTol,
                    params.minIsotopesIonQuant,
                    params.minScansIonQuant,
                    !params.isIMdata
            );
            api.buildIndex();
            return api;
        } else {
            return null;
        }
    }

    public void processLinesBlock(ArrayList<PSM> cBlock, PrintWriter fragmentOutWriter) {
        StringBuilder fragmentBlock = new StringBuilder();
        for (PSM psm : cBlock) {
            processPSM(psm);
            fragmentBlock.append(psm.glycanAssignmentResult.printGlycoFragmentInfo());
        }
        printLines(fragmentOutWriter, fragmentBlock.toString());
    }

    private static synchronized void printLines(PrintWriter out, String linesBlock) {
        out.print(linesBlock);
    }

    /**
     * Read the PSM level glycan assignment results to determine glycan fragment probabilities for each glycan in the database.
     * Option to save prevalence file for diagnostics/info to be added?
     *
     * @return Map of glycan string : fragment propensities container
     */
    public HashMap<String, GlycanCandidateFragments> computeGlycanFragmentProbs(GlycoParams glycoParams) {
        HashMap<String, GlycanCandidateFragments> glycanCandidateFragmentsMap = new HashMap<>();
        HashMap<String, ArrayList<GlycanCandidateResult>> glycanInputMap = new HashMap<>();    // container for glycan: glycan fragment info (read in from file)

        // read all glycan info in
        for (GlycanAssignmentResult result : allResults) {
            if (result.foundGlycan) {
                GlycanCandidateResult glycan = result.bestCandidate;
                GlycanCandidate fragmentInfoContainer = new GlycanCandidate(glycan.composition, 0, false, glycoParams.glycanResiduesMap, glycan.Yfragments, glycan.oxoniumFragments);

                String glycanHash = Glycan.toGlycanString(fragmentInfoContainer.composition);
                // only include good targets in fragment info
                if (result.glycanQval < glycoParams.glycoFDR) {
                    if (glycoParams.minYsForConsensus > 0) {
                        int foundYs = 0;
                        boolean notEnoughYs = true;
                        for (GlycanFragment yFragment : glycan.Yfragments.values()) {
                            if (yFragment.foundIntensity > 0) {
                                foundYs++;
                                if (foundYs >= glycoParams.minYsForConsensus) {
                                    notEnoughYs = false;
                                    break;
                                }
                            }
                        }
                        if (notEnoughYs) {
                            continue;
                        }
                    }
                    if (glycanInputMap.containsKey(glycanHash)) {
                        glycanInputMap.get(glycanHash).add(glycan);
                    } else {
                        ArrayList<GlycanCandidateResult> newList = new ArrayList<>();
                        newList.add(glycan);
                        glycanInputMap.put(glycanHash, newList);
                    }

                    // add to delta mass map for calculating glycan prevalence priors (targets and decoys)
                    double deltaMass = result.deltaMass;
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
        }


        // summarize results for each glycan to get final fragment propensities
        for (Map.Entry<String, ArrayList<GlycanCandidateResult>> glycanEntry : glycanInputMap.entrySet()) {
            // Determine the fragment likelihoods based on all PSMs for this entry
            HashMap<String, Integer> YCounts = new HashMap<>();
            HashMap<String, Integer> OxCounts = new HashMap<>();
            HashMap<String, ArrayList<Double>> YInts = new HashMap<>();
            HashMap<String, ArrayList<Double>> OxInts = new HashMap<>();
            ArrayList<GlycanCandidateResult> allPSMsWithThisGlycan = glycanEntry.getValue();
            // skip generating fragment information for glycans with too few PSMs to get reasonable values
            if (allPSMsWithThisGlycan.size() < MIN_GLYCO_PSMS_FOR_BOOTSTRAP) {
                continue;
            }
            if (glycoParams.topPctSpectraForConsensus < 1.0) {
                // sort PSMs by glycan score and keep only the top X%
                allPSMsWithThisGlycan.sort(Comparator.comparingDouble((GlycanCandidateResult result) -> result.glycanScore).reversed());
                int numToKeep = (int) Math.ceil(allPSMsWithThisGlycan.size() * glycoParams.topPctSpectraForConsensus);
                allPSMsWithThisGlycan = new ArrayList<>(allPSMsWithThisGlycan.subList(0, numToKeep));
            }

            for (GlycanCandidate inputGlycan : allPSMsWithThisGlycan) {
                // read all fragments from each input glycan into the count database
                for (String fragmentHash : inputGlycan.Yfragments.keySet()) {
                    int count = YCounts.getOrDefault(fragmentHash, 0);
                    count++;
                    YCounts.put(fragmentHash, count);
                    if (YInts.containsKey(fragmentHash)) {
                        YInts.get(fragmentHash).add(inputGlycan.Yfragments.get(fragmentHash).foundIntensity);
                    } else {
                        ArrayList<Double> newList = new ArrayList<>();
                        newList.add(inputGlycan.Yfragments.get(fragmentHash).foundIntensity);
                        YInts.put(fragmentHash, newList);
                    }
                }
                for (String fragmentHash : inputGlycan.oxoniumFragments.keySet()) {
                    int count = OxCounts.getOrDefault(fragmentHash, 0);
                    count++;
                    OxCounts.put(fragmentHash, count);
                    if (OxInts.containsKey(fragmentHash)) {
                        OxInts.get(fragmentHash).add(inputGlycan.oxoniumFragments.get(fragmentHash).foundIntensity);
                    } else {
                        ArrayList<Double> newList = new ArrayList<>();
                        newList.add(inputGlycan.oxoniumFragments.get(fragmentHash).foundIntensity);
                        OxInts.put(fragmentHash, newList);
                    }
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

            // save intensities
            // todo: test median vs average
            HashMap<String, Double> yFragmentIntensities = new HashMap<>();
            for (Map.Entry<String, ArrayList<Double>> fragmentEntry : YInts.entrySet()) {
                double[] intensities = new double[fragmentEntry.getValue().size()];
                for (int i = 0; i < fragmentEntry.getValue().size(); i++) {
                    intensities[i] = fragmentEntry.getValue().get(i);
                }
                Median median = new Median();
//                double testAvg = Arrays.stream(intensities).average().orElse(0);
//                double testMedian = median.evaluate(Arrays.stream(intensities).toArray());
                yFragmentIntensities.put(fragmentEntry.getKey(), median.evaluate(Arrays.stream(intensities).toArray()));
            }
            HashMap<String, Double> OxFragmentIntensities = new HashMap<>();
            for (Map.Entry<String, ArrayList<Double>> fragmentEntry : OxInts.entrySet()) {
                double[] intensities = new double[fragmentEntry.getValue().size()];
                for (int i = 0; i < fragmentEntry.getValue().size(); i++) {
                    intensities[i] = fragmentEntry.getValue().get(i);
                }
                Median median = new Median();
                OxFragmentIntensities.put(fragmentEntry.getKey(), median.evaluate(Arrays.stream(intensities).toArray()));
            }

            // save determined propensities to the output container
            GlycanCandidateFragments fragmentInfo = new GlycanCandidateFragments(yFragmentProps, OxFragmentProps, yFragmentIntensities, OxFragmentIntensities);
            glycanCandidateFragmentsMap.put(glycanEntry.getKey(), fragmentInfo);
        }
        return glycanCandidateFragmentsMap;
    }

    /**
     * Glycan FDR wrapper/main method. Called after the glycan assignment is done on PSMs to compute glycan FDR.
     * Handles various old and new methods and LDA.
     */
    public void runScoresAndFDR() {
        // LDA method
        if (glycoParams.glycoLDA) {
            ScoreLDA lda = new ScoreLDA();
            // add PSM results to LDA
            for (GlycanAssignmentResult result: allResults) {
                if (result.foundGlycan) {
                    if (result.bestTarget != null) {
                        lda.targetData.add(result.bestTarget.featureVec);
                    }
                    if (result.bestDecoy != null) {
                        lda.decoyData.add(result.bestDecoy.featureVec);
                    }
                }
            }
            lda.runLDA(allResults, ldaHeader, glycoParams.ldaTargetProp);
        }

        // Compute FDR
        PTMShepherd.print("Calculating Glycan FDR");
        boolean fdrSuccess = false;
        if (glycoParams.useNonCompFDR) {
            fdrSuccess = computeFDRNonCompetitive(allResults, glycoParams.glycoFDR);
        } else {
            boolean firstFDRsuccess = computeFDRcompetitive(allResults, glycoParams.glycoFDR);
            if (!firstFDRsuccess) {
                fdrSuccess = computeFDRNonCompetitive(allResults, glycoParams.glycoFDR);
            } else {
                fdrSuccess = true;
            }
        }

        try {
            PrintWriter glycoOut = new PrintWriter(new FileWriter(glycoFile));
            glycoOut.println(String.format("%s\t%s\t%s\t%s\t%s", "Spectrum", "Peptide", "Mods", "Pep Mass", "Mass Shift") + String.format("\t%s\tGlycan Score\tGlycan q-value\tBest Decoy Glycan\tBest Decoy Score", GLYCAN_COMP_COL_NAME) + ldaHeader + "Fragments:");

            printLines(glycoOut, allResults.stream().sorted(Comparator.comparingInt(result -> result.psmLineIndex))
                    .map(GlycanAssignmentResult::printGlycoFragmentInfo)
                    .collect(Collectors.joining("")));
            glycoOut.close();

            // for debugging
            PrintWriter glycoOut2 = new PrintWriter(new FileWriter(glycoFile + "2"));
            glycoOut2.println(String.format("%s\t%s\t%s\t%s\t%s", "Spectrum", "Peptide", "Mods", "Pep Mass", "Mass Shift") + String.format("\t%s\tGlycan Score\tGlycan q-value", GLYCAN_COMP_COL_NAME) + ldaHeader + "Fragments:");
            printLines(glycoOut2, allResults.stream().sorted(Comparator.comparingInt(result -> result.psmLineIndex))
                    .map(GlycanAssignmentResult::printAllCandidates)
                    .collect(Collectors.joining("")));
            glycoOut2.close();
        } catch (IOException e) {
            PTMShepherd.die("Error writing to glyco file " + glycoFile.getAbsolutePath() + "\n" + e.getMessage());
        }
        if (!fdrSuccess) {
            PTMShepherd.die("Stopping after failed glycan FDR estimation.");
        }
    }

    /**
     * Determines the score threshold based on the specified FDR.
     *
     * @param fdrCutOff Maximum acceptable FDR
     * @param results List of GlycanAssignmentResults containing target and decoy scores
     * @return The score threshold
     */
    public static boolean computeFDRcompetitive(List<GlycanAssignmentResult> results, double fdrCutOff) {
        // Sort target scores in descending order
        results.sort(Comparator.comparingDouble((GlycanAssignmentResult result) -> result.glycanScore).reversed());
        long glycoResultCount = results.stream().filter(r -> r.foundGlycan).count();

        // Get decoy indices in the combined sorted list
        List<Integer> decoyIndexes = getDecoyIndexes(results);

        int decoyCount = decoyIndexes.size();
        int targetCount = (int) glycoResultCount - decoyCount;
        // check if enough decoys were found (i.e., initial q-val is above the desired threshold)
        double initialFDR = calculateFDR(targetCount, decoyCount, false);
        if (initialFDR < fdrCutOff) {
            PTMShepherd.print(String.format("\tNot enough decoys to compute FDR at %.1f%% with competitive method, started at %.2f%%", fdrCutOff * 100, initialFDR * 100));
            return false;
        }

        double currentMinQ = Double.MAX_VALUE;
//        double scoreThreshold = Double.NaN;
        boolean foundThreshold = false;
        // Calculate FDR at each point
        for (int i = results.size() - 1; i >= 0; i--) {
            GlycanAssignmentResult result = results.get(i);
            if (!result.foundGlycan) {
                continue; // skip PSMs without glycan assignments
            }

            if (result.isDecoyGlycan) {
                decoyCount--;
            } else {
                targetCount--;
            }
            double fdr = Math.min(calculateFDR(targetCount, decoyCount, false), currentMinQ);        // q = (d+1)/t recommended per 10.1021/acs.jproteome.6b00144
            if (fdr < currentMinQ) {
                currentMinQ = fdr;
            }
            if (!foundThreshold) {
                if (fdr < fdrCutOff) {
//                    scoreThreshold = results.get(i).glycanScore;
                    PTMShepherd.print(String.format("Found glycan score threshold: %.2f with %d decoys, %d targets for %.2f%% estimated FDR (%d total inputs)",
                            result.glycanScore, decoyCount, targetCount, fdr * 100, glycoResultCount));
                    foundThreshold = true;
                }
            }
            result.glycanQval = result.isDecoyGlycan ? 1.0 : fdr;
        }
        if (!foundThreshold) {
            // could not reach threshold at all - stop the analysis
            PTMShepherd.print("Could not reach glycan FDR threshold of " + fdrCutOff * 100 + "% (insufficient target/decoy separation). " +
                    "Please check the search parameters and/or try a different glycan database.");
            return false;
        }
        return true;
    }

    /**
     * Determines the score threshold based on the specified FDR using the top target AND top decoy for each PSM.
     *
     * @param glycoFDR Maximum acceptable FDR
     * @param results List of GlycanAssignmentResults containing target and decoy scores
     */
    public static boolean computeFDRNonCompetitive(List<GlycanAssignmentResult> results, double glycoFDR) {
        HashMap<String, GlycanAssignmentResult> resultMap = new HashMap<>();

        ArrayList<GlycoScore> scoreDistribution = new ArrayList<>();
        int targets = 0;
        int decoys = 0;
        int totalGlycoResults = 0;
        for (GlycanAssignmentResult result: results) {
            resultMap.put(result.specName, result);

            if (result.foundGlycan) {
                totalGlycoResults++;
                // record top target and top decoy score
                if (result.isDecoyGlycan) {
                    decoys++;
                    scoreDistribution.add(new GlycoScore(result.glycanScore, true, result.specName, true));
                    if (result.bestTarget != null) {
                        targets++;
                        scoreDistribution.add(new GlycoScore(result.bestTarget.glycanScore, false, result.specName, false));
                    }
                } else {
                    targets++;
                    scoreDistribution.add(new GlycoScore(result.glycanScore, false, result.specName, true));
                    if (result.bestDecoy != null) {
                        decoys++;
                        scoreDistribution.add(new GlycoScore(result.bestDecoy.glycanScore, true, result.specName, false));
                    }
                }
            }
        }

        // sort scoreMap in order of ascending score
        scoreDistribution.sort(GlycoScore::compareTo);

        double targetDecoyRatio;
        double currentMinQ = 1;
        double scoreThreshold;
        boolean foundScoreThresh = false;
        for (GlycoScore scoreObj : scoreDistribution) {
            if (!scoreObj.isDecoy) {
                targets--;
            } else {
                decoys--;
            }
            // compute TD ratio
            targetDecoyRatio = calculateFDR(targets, decoys, true);
            if (decoys > targets) {
                targetDecoyRatio = 1.0;     // cap FDR at 1
            } else if (targets == 0) {
                targetDecoyRatio = 0.0;     // min FDR = 0. Using else-if with the above block so that if decoys are nonzero with 0 targets, FDR = 1
            }

            // compute q-value and save for later
            double qval = Math.min(targetDecoyRatio, currentMinQ);
            if (qval < currentMinQ) {
                currentMinQ = qval;
            }
            if (scoreObj.isFromTopCandidate) {
                resultMap.get(scoreObj.spectrumID).glycanQval = scoreObj.isDecoy ? 1.0 : qval;   // save q-val only for top candidates
            }

            // check for the score threshold that gives the requested FDR
            if (!foundScoreThresh) {
                if (targetDecoyRatio <= glycoFDR) {
                    // stop here, found cutoff
                    scoreThreshold = scoreObj.score;
                    PTMShepherd.print(String.format("\tFound score threshold of %.2f for %.1f%% FDR with %d targets and %d decoys from non-competitive method (%d total inputs)", scoreThreshold, targetDecoyRatio * 100, targets, decoys, totalGlycoResults));
                    foundScoreThresh = true;
                }
            }
        }
        if (!foundScoreThresh) {
            // could not reach threshold at all - stop the analysis
            PTMShepherd.print("Could not reach glycan FDR threshold of " + glycoFDR * 100 + "% (insufficient target/decoy separation). " +
                    "Please check the search parameters and/or try a different glycan database.");
            return false;
        }
        return true;
    }

    /**
     * Gets the indices of decoy scores in the combined sorted score list.
     */
    private static List<Integer> getDecoyIndexes(List<GlycanAssignmentResult> results) {
        List<Integer> decoyIndexes = new ArrayList<>();
        for (int i = 0; i < results.size(); i++) {
            if (results.get(i).isDecoyGlycan) {
                decoyIndexes.add(i);
            }
        }
        return decoyIndexes;
    }

    /**
     * Calculate FDR from target and decoy counts
     *
     * @param targets target count
     * @param decoys  decoy count
     * @return FDR
     */
    private static double calculateFDR(int targets, int decoys, boolean useNonCompFDR) {
        if (useNonCompFDR) {
            return (2 * decoys) / (double) (decoys + targets);
        } else {
            return (decoys + 1) / (double) targets;
        }
    }

    /**
     * Determine the width of mass errors in PSMs without delta mass to use for mass error probability
     * estimation. Returns the sigma of a Gaussian distribution fit to the mass errors of all PSMs without
     * delta masses (isotope corrected)
     *
     * @param psmFile PSM file to analyze
     * @param clines  line numbers in the PSM file?
     */
    public void getMassErrorsFirstPass(PSMFile psmFile, ArrayList<Integer> clines) {
        ArrayList<Double> massErrors = new ArrayList<>();

        // Get mass errors for PSMs with delta mass in exclusion range (-1.5 to 3.5)
        double minError = 10;
        double maxError = -10;
        for (Integer cline : clines) {//for relevant line in curr spec file
            float deltaMass = psmFile.psms.get(cline).getDMass();

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

        computeMassErrorsHelper(massErrors, maxError, minError);
    }


    /**
     * Calculate the mass and isotope error distributions for PSMs from the first pass of glycan assignment.
     * @param results results from 1st pass
     */
    public void getMassErrorsSecondPass(ArrayList<GlycanAssignmentResult> results) {
        ArrayList<Double> massErrors = new ArrayList<>();
        HashMap<Integer, Integer> isotopeCounts = new HashMap<>();
        double minError = 10;
        double maxError = -10;

        for (GlycanAssignmentResult result : results) {
            // only use target glyco PSMs that passed FDR
            if (result.foundGlycan && !result.isDecoyGlycan && result.glycanQval < glycoParams.glycoFDR) {
                // result.deltaMass is the PSM delta mass. Subtract the best candidate mass and isotope to get the final mass error
                double massError = result.deltaMass - result.bestCandidate.mass - result.bestCandidate.isotope * AAMasses.averagineIsotopeMass;
                massErrors.add(massError);
                if (massError > maxError) {
                    maxError = massError;
                }
                if (massError < minError) {
                    minError = massError;
                }
                isotopeCounts.put(result.bestCandidate.isotope, isotopeCounts.getOrDefault(result.bestCandidate.isotope, 0) + 1); // count the number of PSMs for each isotope
            }
        }
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.mass2nd)) {
            computeMassErrorsHelper(massErrors, maxError, minError);
        }
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.iso2nd)) {
            glycoParams.updateIsotopesProbsFromFirstPass(isotopeCounts);
        }
    }

    private void computeMassErrorsHelper(ArrayList<Double> massErrors, double maxError, double minError) {
        if (massErrors.size() < 200) {
            // not enough unmodified PSMs to compute mass error stats - use defaults
            PTMShepherd.print("\tNot enough unmodified PSMs to determine mass error distribution, using default values");
            massErrorWidth = 0.005;
            meanMassError = 0.0;
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
     *
     * @param psm String of a single PSM line
     * @return result container
     */
    public void processPSM(PSM psm) {
        // get basic info
        GlycanAssignmentResult glycoResult = new GlycanAssignmentResult(psm.lineNum, psm.getPeptide(), psm.getDMass(), psm.getCalcPepmass(), psm.printAssignedMods(), psm.getSpec());

        // read spectrum and condition
        Spectrum spec = mr.getSpectrum(psm.getSpec());
        if (spec == null) {
            this.lineWithoutSpectra.add(psm.getSpec());
            psm.glycanAssignmentResult = glycoResult;
            return;
        }
        spec.conditionOptNorm(condPeaks, condRatio, false);

        // do glycan assignment
        glycoResult = assignGlycanToPSM(spec, glycoResult, glycanDatabase, massErrorWidth, meanMassError);
        psm.glycanAssignmentResult = glycoResult;
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
            return glycoResult;
        }

        // Determine possible glycan candidates from mass
        ArrayList<GlycanCandidateResult> searchCandidates = getMatchingGlycansByMass(glycoResult.pepMass, glycoResult.deltaMass, glycanDatabase, glycoParams.glycoIsotopes, glycoParams.glycoPPMtol);
        if (!searchCandidates.isEmpty()) {
            // Search Y and oxonium ions in spectrum for each candidate
            float ppmTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
            double spectrumYIntensity = spec.getYintensity(glycoResult.pepMass);
            for (GlycanCandidateResult candidate : searchCandidates) {
                double foundYIntensity = 0;
                for (GlycanFragment yFragment : candidate.Yfragments.values()) {
                    yFragment.foundIntensity = spec.findIonNeutral(yFragment.neutralMass + glycoResult.pepMass, ppmTol, Integer.parseInt(PTMShepherd.getParam("spectra_maxPrecursorCharge"))) / spec.basePeakInt;  // sum of charge state intensities if >1 found
                    foundYIntensity += yFragment.foundIntensity;
                }
                for (GlycanFragment oxoniumFragment : candidate.oxoniumFragments.values()) {
                    // save oxonium ion intensity relative to base peak
                    oxoniumFragment.foundIntensity = spec.findIon(oxoniumFragment.neutralMass + AAMasses.protMass, ppmTol) / spec.basePeakInt;
                }
                candidate.YproportionScore = spectrumYIntensity == 0 ? 0 : (foundYIntensity * 100) / spectrumYIntensity;    // proportion of possible Y ions in the spectrum matched to the candidate
                if (glycoParams.normFragmentIntensities) {
                    GlycanCandidateResult.normalizeIntensities(candidate.Yfragments);
//                    GlycanCandidateResult.normalizeIntensities(candidate.oxoniumFragments);   // todo: enable once generalized oxos available
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
                    comparisonScore = pairwiseCompareDynamic(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass, glycoResult.pepMass, spec);
                } else {
                    comparisonScore = pairwiseCompareStatic(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass, meanMassError, glycoResult.pepMass, spec);
                }

                if (comparisonScore == 0) {
                    // exact same score (e.g., from target/decoy if no Y/oxo ions found and using decoy mass = target mass)
                    // Use a random tiebreaker to avoid bias from always picking the target
                    comparisonScore += (glycoParams.randomGenerator.nextDouble() - 0.5) * 1E-6; // random yields a value between 0-1, so subtracting 0.5 gives (approx) equal chance of being positive or negative
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
            for (int i = 0; i < bestCandidateIndex; i++) {
                if (useFragmentSpecificProbs) {
                    scoresVsBestCandidate[i] = pairwiseCompareDynamic(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass, glycoResult.pepMass, spec);
                } else {
                    scoresVsBestCandidate[i] = pairwiseCompareStatic(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult.deltaMass, meanMassError, glycoResult.pepMass, spec);
                }
                if (scoresVsBestCandidate[i] == 0) {
                    // this was from a tiebreak that was already chosen as not the best candidate, so set to a small positive value for sorting (since sort is done ascending)
                    scoresVsBestCandidate[i] = 1E-6;
                }
                if (scoresVsBestCandidate[i] < 0) {
                    // after direct comparison, this candidate is now better than the best candidate. But setting this as best and going back would cause an infinite loop, leave the original best candidate in place.
                    scoresVsBestCandidate[i] = 1E-5;
                }
            }
            scoresVsBestCandidate[bestCandidateIndex] = 0;

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
            if (useFragmentSpecificProbs) {
                computeAbsoluteScoreDynamic(spec, searchCandidates.get(bestCandidateIndex), glycoResult, massErrorWidth, meanMassError);
            } else {
                computeAbsoluteScore(spec, searchCandidates.get(bestCandidateIndex), glycoResult, massErrorWidth, meanMassError);
            }
            // save candidates to the result (in descending order of scores)
            glycoResult.bestCandidate = searchCandidates.get(bestCandidateIndex);
            for (int index : sortedIndicesOfBestScores) {
                glycoResult.allCandidates.add(searchCandidates.get(index));
            }
            glycoResult.summedScore = glycoResult.bestCandidate.summedScore;
            glycoResult.glycanScore = glycoResult.bestCandidate.summedScore;    // overwritten later if using LDA
            glycoResult.foundGlycan = true;

            // if top glycan is a decoy, also write best target and best target score to subsequent columns
            boolean bestWasTarget = !searchCandidates.get(bestCandidateIndex).isDecoy;
            glycoResult.isDecoyGlycan = !bestWasTarget;
            getNextGlycanScores(spec, bestWasTarget, glycoResult, massErrorWidth, meanMassError);
        }

        return glycoResult;
    }

    /**
     * Helper method to get output for next best glycan (i.e., best decoy if the top hit is a target and vice versa).
     * Returns formatted output string handling various cases
     *
     * @param glycoResult
     * @param massErrorWidth
     * @param meanMassError
     */
    private void getNextGlycanScores(Spectrum spec, boolean bestWasTarget, GlycanAssignmentResult glycoResult, double massErrorWidth, double meanMassError) {
        if (bestWasTarget) {
            glycoResult.bestTarget = glycoResult.bestCandidate;
        } else {
            glycoResult.bestDecoy = glycoResult.bestCandidate;
        }

        boolean foundNext = false;
        // compute scores for all candidates (except the best, since it was already computed)
        for (int i = 1; i < glycoResult.allCandidates.size(); i++) {
            GlycanCandidateResult nextCandidate = glycoResult.allCandidates.get(i);
            if (useFragmentSpecificProbs) {
                computeAbsoluteScoreDynamic(spec, nextCandidate, glycoResult, massErrorWidth, meanMassError);
            } else {
                computeAbsoluteScore(spec, nextCandidate, glycoResult, massErrorWidth, meanMassError);
            }

            // if Best hit was target, look for Next to be decoy and vice versa
            if (!foundNext) {
                if (nextCandidate.isDecoy == bestWasTarget) {
                    // add the next hit's information
                    foundNext = true;
                    if (bestWasTarget) {
                        // save best to target, next to decoy
                        glycoResult.bestDecoy = nextCandidate;
                        glycoResult.bestTarget = glycoResult.bestCandidate;
                    } else {
                        // save best to decoy, next to target
                        glycoResult.bestTarget = nextCandidate;
                        glycoResult.bestDecoy = glycoResult.bestCandidate;
                    }
                }
            }
        }
    }

    /**
     * Glycan frequency score calculator. Essentially a prior for how likely a given glycan is given delta mass bin.
     * Normalized to the most frequent glycan in the bin.
     * @param candidate   glycan candidate
     * @param deltaMass delta mass bin in question
     * @return frequency (between 0 and 1) of the glycan in the delta mass bin.
     */
    public double computeGlycanFrequencyScore(GlycanCandidate candidate, double deltaMass) {
        // determine the overall likelihood priors of these glycans given the observed delta mass
        int candidateCount = 0;
        int maxFrequency = 0;
        // Decoys are not included in the saved glycans. Use the frequency of the corresponding target
        String glycanHash = Glycan.toGlycanString(candidate.composition);
        HashMap<String, Integer> emptyMap = new HashMap<>();
        int massBin = (int) Math.floor(deltaMass);
        HashMap<String, Integer> glycanCountMap = glycanMassBinMap.getOrDefault(massBin, emptyMap);
        if (!glycanCountMap.isEmpty()) {
            // count instances of glycan 1, glycan 2, and all glycans
            candidateCount = candidateCount + glycanCountMap.getOrDefault(glycanHash, 0);
            for (int glycanCount : glycanCountMap.values()) {
                if (glycanCount > maxFrequency) {
                    maxFrequency = glycanCount; // save the max frequency of any glycan in this bin
                }
            }
        }
        if (maxFrequency == 0) {
            return 0.0;     // no glycans passed FDR in this bin (in 1st pass), return 0
        }
        return candidateCount / (double) maxFrequency;
    }

    public double pairwiseCompareDynamic(GlycanCandidateResult glycan1, GlycanCandidateResult glycan2, double deltaMass, double pepMass, Spectrum spec) {
        // calculate fragment-specific prob estimates based on observed fragment ions
        double sumLogRatio = 0;

        // Y ions
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.yscore))
        {
            if (glycoParams.glycoYnorm) {
                sumLogRatio += pairwiseCompareDynamicNormed(glycan1.Yfragments, glycan2.Yfragments, glycan1, glycan2);
            } else {
                sumLogRatio += pairwiseCompareDynamicNotNorm(glycan1.Yfragments, glycan2.Yfragments, glycan1, glycan2);
            }
        }

        // oxonium ions
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.oxo)) {
            sumLogRatio += pairwiseCompareDynamicNotNorm(glycan1.oxoniumFragments, glycan2.oxoniumFragments, glycan1, glycan2);
        }

        // mass and isotope error
        sumLogRatio += computeMassIsoScorePairwise(glycan1, glycan2, deltaMass, meanMassError);

        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.kl)) {
            double ms1score1 = calculateMS1score(glycan1, spec, pepMass);
            double ms1score2 = calculateMS1score(glycan2, spec, pepMass);
            sumLogRatio += (ms1score1 - ms1score2);
        }

        return sumLogRatio;
    }

    /**
     * Compute sum log probability ratios for the compared glycans for a particular fragment type.
     *
     * @param fragmentsMap1 Fragments from candidate 1
     * @param fragmentsMap2 Fragments from candidate 2
     * @param glycan2       candidate 2
     * @return sum log probability with normalization included
     */
    public double pairwiseCompareDynamicNormed(TreeMap<String, GlycanFragment> fragmentsMap1, TreeMap<String, GlycanFragment> fragmentsMap2, GlycanCandidate glycan1, GlycanCandidate glycan2) {
        double sumLogRatio = 0;
        double unique1score = 0;
        double unique2score = 0;
        int unique1count = 0;
        int unique2count = 0;
        for (GlycanFragment fragment1 : fragmentsMap1.values()) {
            double probRatio;
            if (fragment1.isAllowedFragment(glycan2, glycoParams.glycanResiduesMap)) {
                GlycanFragment fragment2 = fragmentsMap2.get(fragment1.hash);
                probRatio = computeFragmentPairwiseScore(fragment1, fragment2);
                sumLogRatio += Math.log(probRatio);
            } else {
                // fragment only possible for glycan 1, use glycan 1 only estimate
                probRatio = computeFragmentAbsoluteScore(fragment1);
                unique1score += Math.log(probRatio);
                unique1count++;
            }
        }
        // glycan 2 fragments - unique fragments get scored the same way as in the absolute method, and subtracted since they support glycan 2 not 1
        for (GlycanFragment fragment : fragmentsMap2.values()) {
            if (!fragment.isAllowedFragment(glycan1, glycoParams.glycanResiduesMap)) {
                double probRatio = computeFragmentAbsoluteScore(fragment);
                unique2score += Math.log(probRatio);
                unique2count++;
            }
        }
        double normScore1 = unique1count > 0 ? unique1score / Math.sqrt(unique1count) : 0;  // avoid divide by 0 if no unique fragments
        double normScore2 = unique2count > 0 ? unique2score / Math.sqrt(unique2count) : 0;
        sumLogRatio += (normScore1 - normScore2);
        return sumLogRatio;
    }

    /**
     * Compute sum log probability ratios for the compared glycans for a particular fragment type.
     *
     * @param fragmentsMap1 Fragments from candidate 1
     * @param fragmentsMap2 Fragments from candidate 2
     * @param glycan2       candidate 2
     * @return sum log probability with normalization included
     */
    public double pairwiseCompareDynamicNotNorm(TreeMap<String, GlycanFragment> fragmentsMap1, TreeMap<String, GlycanFragment> fragmentsMap2, GlycanCandidate glycan1, GlycanCandidate glycan2) {
        double sumLogRatio = 0;
        for (GlycanFragment fragment1 : fragmentsMap1.values()) {
            double probRatio;
            if (fragment1.isAllowedFragment(glycan2, glycoParams.glycanResiduesMap)) {
                GlycanFragment fragment2 = fragmentsMap2.get(fragment1.hash);
                probRatio = computeFragmentPairwiseScore(fragment1, fragment2);
            } else {
                // fragment only possible for glycan 1, use glycan 1 only estimate
                probRatio = computeFragmentAbsoluteScore(fragment1);
            }
            sumLogRatio += Math.log(probRatio);
        }
        // glycan 2 fragments - unique fragments get scored the same way as in the absolute method, and subtracted since they support glycan 2 not 1
        for (GlycanFragment fragment : fragmentsMap2.values()) {
            if (!fragment.isAllowedFragment(glycan1, glycoParams.glycanResiduesMap)) {
                double probRatio = computeFragmentAbsoluteScore(fragment);
                sumLogRatio -= Math.log(probRatio);
            }
        }
        return sumLogRatio;
    }

    /**
     * Compute the "absolute" score of this glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     *
     * @param candidate     glycan candidate to calculate score for
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     * @param meanMassError  mean mass error of non-delta mass peptides
     */
    public void computeAbsoluteScoreDynamic(Spectrum spec, GlycanCandidateResult candidate, GlycanAssignmentResult result, double massErrorWidth, double meanMassError) {
        // Y ions
        int index = 0;
        double yScore = 0;
        double[] foundYs = new double[candidate.Yfragments.size()];
        double[] expectedYs = new double[candidate.Yfragments.size()];
        for (GlycanFragment fragment : candidate.Yfragments.values()) {
            double probRatio = computeFragmentAbsoluteScore(fragment);
            yScore += Math.log(probRatio);
            foundYs[index] = fragment.foundIntensity;
            expectedYs[index] = fragment.expectedIntensity;
            index++;
        }
        if (glycoParams.glycoYnorm) {
            yScore = yScore / Math.sqrt(candidate.Yfragments.size());
        }
        candidate.YFragmentScore = yScore;
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ysim)) {
            candidate.ySpecSim = entropyScore(expectedYs, foundYs);
        }

        // oxonium ions
        index = 0;
        double oxoScore = 0;
        double[] foundOxos = new double[candidate.oxoniumFragments.size()];
        double[] expectedOxos = new double[candidate.oxoniumFragments.size()];
        for (GlycanFragment fragment : candidate.oxoniumFragments.values()) {
            double probRatio = computeFragmentAbsoluteScore(fragment);
            oxoScore += Math.log(probRatio);
            foundOxos[index] = fragment.foundIntensity;
            expectedOxos[index] = fragment.expectedIntensity;
            index++;
        }
        candidate.OxFragmentScore = oxoScore;
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.oxsim)) {
            candidate.oxSpecSim = entropyScore(expectedOxos, foundOxos);
        }

        // isotope and mass errors. Isotope is ratio relative to no isotope error (0)
        candidate.isotopeScore = computeIsoScoreAbs(candidate, result);
        candidate.massErrorScore = computeMassScoreAbs(candidate, result, massErrorWidth, meanMassError);

        // only calculate MS1 score if requested because it requires slow index building
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.kl)) {
            candidate.ms1Score = calculateMS1score(candidate, spec, result.pepMass);
        }
        if (!isFirstPass) {
            candidate.frequencyPrior = computeGlycanFrequencyScore(candidate, result.deltaMass);
        }
        generateScores(candidate);
    }

    /**
     * Compute propensity-specific score for fragment ion that is NOT unique (i.e., shared between two candidates).
     * Score is the ratio of the propensities for the two candidates, positive towards candidate 1 if found in the
     * spectrum or towards candidate 2 if not.
     *
     * @param fragment1 fragment from glycan 1
     * @param fragment2 fragment from glycan 2
     * @return ratio of fragment probs
     */
    public double computeFragmentPairwiseScore(GlycanFragment fragment1, GlycanFragment fragment2) {
        double probRatio;
        if (fragment1.propensity > 0 && fragment2.propensity > 0) {
            if (fragment1.foundIntensity > 0) {
                // "hit": fragment found in spectrum. Compute prob of glycans given the presence of this ion
                probRatio = computePropensityRatio(fragment1, fragment2);
                // todo: optional, test multiplying by intensity ratio (obs/exp) for fragment1 only
            } else {
                // "miss": fragment not found. Compute prob of glycans given absence of this ion. Miss propensity = 1 - hit propensity
                probRatio = computePropensityRatio(fragment2, fragment1);
            }
        } else {
            probRatio = 1;
        }
        return probRatio;
    }

    /**
     * Helper for computing absolute score of Y or oxonium fragments. Uses empirical probability for this fragment
     * type and weights it by intensity vs expected
     *
     * @param fragment fragment to consider
     * @return sum log ratio of fragment probs
     */
    public double computeFragmentAbsoluteScore(GlycanFragment fragment) {
        double probRatio;
        if (fragment.foundIntensity > 0) {
            // only compute fragment intensity ratio for oxonium ions, not Y
            double intensityRatio = fragment.fragType == GlycanFragment.FragType.Ox ? computeIntensityRatio(fragment) : 1.0;
            probRatio = fragment.ruleProbabilities[0] * intensityRatio;     // found in spectrum - ion supports this glycan
        } else {
            if (fragment.fragType == GlycanFragment.FragType.Y) {
                if (fragment.propensity > 0) {
                    // Weight misses by propensity (if known), so that unlikely Y ions don't over-penalize a reasonable glycan
                    // formula sets range from prob = 1 as prop -> 0 to prob = ruleprobs[1] as prop -> 1
                    probRatio = 1 - fragment.propensity * (1 - fragment.ruleProbabilities[1]);
                } else {
                    probRatio = fragment.ruleProbabilities[1];     // not found in spectrum - ion does not support this glycan
                }
            } else {
                probRatio = fragment.ruleProbabilities[1];     // not found in spectrum - ion does not support this glycan
            }
        }
        return probRatio;
    }

    /**
     * Perform pairwise comparison of two glycans. Uses sum of log probability ratios between candidates for
     * each category (mass/iso error and fragment ion) being considered. Returns a single score of combined
     * probability of first glycan candidate over second.
     *
     * @param glycan1   candidate 1
     * @param glycan2   candidate 2
     * @param deltaMass observed delta mass
     * @return output probability score (sum of log ratios)
     */
    public double pairwiseCompareStatic(GlycanCandidateResult glycan1, GlycanCandidateResult glycan2, double deltaMass, double meanMassError, double pepMass, Spectrum spec) {
        double sumLogRatio = 0;
        // Y ions
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.yscore)) {
            sumLogRatio += pairwiseCompareYstatic(glycan1, glycan2, glycoParams.glycoYnorm);
        }

        // oxonium ions
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.oxo)) {
            sumLogRatio += pairwiseCompareOxoStatic(glycan1, glycan2);
        }

        // isotope and mass errors
        sumLogRatio += computeMassIsoScorePairwise(glycan1, glycan2, deltaMass, meanMassError);

        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.kl)) {
            double ms1score1 = calculateMS1score(glycan1, spec, pepMass);
            double ms1score2 = calculateMS1score(glycan2, spec, pepMass);
            sumLogRatio += (ms1score1 - ms1score2);
        }

        return sumLogRatio;
    }

    /**
     * Compute sum log probability ratios for the compared glycans for a particular fragment type with normalization
     * of misses.
     * NOTE: does NOT allow fragment specific probabilities (for missed ions)
     *
     * @param glycan1 candidate 1
     * @param glycan2 candidate 2
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
            if (!fragment1.isAllowedFragment(glycan2, glycoParams.glycanResiduesMap)) {
                boolean foundInSpectrum = fragment1.foundIntensity > 0;
                if (foundInSpectrum) {
                    cand1Hits++;
                } else {
                    cand1Misses++;
                }
            }
        }
        for (GlycanFragment fragment2 : glycan2.Yfragments.values()) {
            if (!fragment2.isAllowedFragment(glycan1, glycoParams.glycanResiduesMap)) {
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
     * @param glycan1 candidate 1
     * @param glycan2 candidate 2
     * @return sum log probability with normalization included
     */
    public double pairwiseCompareOxoStatic(GlycanCandidate glycan1, GlycanCandidate glycan2) {
        double sumLogRatio = 0;

        // Loop over each candidate's fragments, scoring unique (i.e., not in the other candidate) fragments as hit/miss if found/not in spectrum
        for (GlycanFragment fragment1 : glycan1.oxoniumFragments.values()) {
            if (!fragment1.isAllowedFragment(glycan2, glycoParams.glycanResiduesMap)) {
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
            if (!fragment2.isAllowedFragment(glycan1, glycoParams.glycanResiduesMap)) {
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
     * @param glycan1       glycan 1
     * @param glycan2       glycan 2
     * @param deltaMass     observed delta mass
     * @param meanMassError mean mass error of non-delta mass peptides
     * @return probability ratio (glycan 1 over 2)
     */
    public double computeMassIsoScorePairwise(GlycanCandidate glycan1, GlycanCandidate glycan2, double deltaMass, double meanMassError) {
        // Determine isotopes
        float iso1 = (float) (deltaMass - glycan1.mass);
        int roundedIso1 = Math.round(iso1);
        float iso2 = (float) (deltaMass - glycan2.mass);
        int roundedIso2 = Math.round(iso2);

        double isotopeProbRatio;
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.iso)) {
            isotopeProbRatio = glycoParams.isotopeProbTable.get(roundedIso1) / glycoParams.isotopeProbTable.get(roundedIso2);
        } else {
            isotopeProbRatio = 1;
        }

        // mass error calc
        double massProbRatio = 1.0;
        if (glycoParams.massProbScaling != 0) {
            double minMassError = deltaMass * (glycoParams.glycoPPMtol * 0.01) * 1e-6;  // min mass error is ppmTol / 100
            double massError1 = deltaMass - glycan1.mass - (roundedIso1 * AAMasses.averagineIsotopeMass);
            double massStDevs1 = massError1 - meanMassError;
            if (Math.abs(massStDevs1) < minMassError)
                massStDevs1 = minMassError;
            double massError2 = deltaMass - glycan2.mass - (roundedIso2 * AAMasses.averagineIsotopeMass);
            double massStDevs2 = massError2 - meanMassError;
            if (Math.abs(massStDevs2) < minMassError)
                massStDevs2 = minMassError;
            if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.mass)) {
                massProbRatio = Math.abs(massStDevs2 / massStDevs1);     // divide #2 by #1 to get ratio for likelihood of #1 vs #2, adjust by scaling factor
            }
        }
        return Math.log(isotopeProbRatio) + Math.log(massProbRatio) * glycoParams.massProbScaling;
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
     *
     * @param bestGlycan glycan candidate to calculate score for
     * @return absolute score
     */
    public double computeYAbsoluteScore(GlycanCandidate bestGlycan) {
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
     *
     * @param bestGlycan glycan candidate to calculate score for
     * @return absolute score
     */
    public double computeOxoAbsoluteScore(GlycanCandidate bestGlycan) {
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
     *
     * @param candidate     glycan candidate to calculate score for
     * @param result         result container
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     * @param meanMassError  mean mass error of non-delta mass peptides
     */
    public void computeAbsoluteScore(Spectrum spec, GlycanCandidateResult candidate, GlycanAssignmentResult result, double massErrorWidth, double meanMassError) {
        if (glycoParams.glycoYnorm) {
            candidate.YFragmentScore = computeYAbsoluteScoreNormed(candidate);
        } else {
            candidate.YFragmentScore = computeYAbsoluteScore(candidate);
        }
        candidate.OxFragmentScore = computeOxoAbsoluteScore(candidate);
        candidate.isotopeScore = computeIsoScoreAbs(candidate, result);
        candidate.massErrorScore = computeMassScoreAbs(candidate, result, massErrorWidth, meanMassError);

        // only calculate MS1 score if requested because it requires slow index building
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.kl)) {
            candidate.ms1Score = calculateMS1score(candidate, spec, result.pepMass);
        }

        generateScores(candidate);
    }

    private double computeMassScoreAbs(GlycanCandidateResult candidate, GlycanAssignmentResult result, double massErrorWidth, double meanMassError) {
        float iso1 = (float) (result.deltaMass - candidate.mass);
        int roundedIso1 = Math.round(iso1);
        double massError1 = result.deltaMass - candidate.mass - (roundedIso1 * AAMasses.averagineIsotopeMass);
        candidate.massError = massError1;
        if (!(glycoParams.massProbScaling == 0)) {
            // mass error is computed in the absolute sense - the number of std devs from mean is used instead of the ratio of two such numbers
            double massStDevs1 = (massError1 - meanMassError) / massErrorWidth;
            double massDist = Math.abs(massStDevs1);
            // Compare to "default" mass error, set to 5 std devs since glycopeps tend to have larger error than regular peps, AND we're more concerned about penalizing large misses
            return Math.log(glycoParams.absScoreErrorParam / massDist) * glycoParams.massProbScaling;
        } else {
            return 1;
        }
    }

    private double computeIsoScoreAbs(GlycanCandidateResult candidate, GlycanAssignmentResult result) {
        float iso1 = (float) (result.deltaMass - candidate.mass);
        int roundedIso1 = Math.round(iso1);
        candidate.isotope = roundedIso1;
        double isotopeProbRatio = glycoParams.isotopeProbTable.get(roundedIso1) / glycoParams.isotopeProbTable.get(0);
        return Math.log(isotopeProbRatio);
    }

    /**
     * Compute propensity adjustment for fragments in common between candidates. Caps it so that low-propensity fragments being found
     * can't reduce the score (same as is done for intensity ratio). Assumes both propensities are non-zero!
     *
     * @param fragment1
     * @param fragment2
     * @return
     */
    public double computePropensityRatio(GlycanFragment fragment1, GlycanFragment fragment2) {
        double propensityRatio;
        if (fragment1.propensity > 0) {
            if (fragment2.propensity > 0) {
                propensityRatio = fragment1.propensity / fragment2.propensity;
                propensityRatio = Math.sqrt(propensityRatio);   // todo: param, test
            } else {
                // only propensity for fragment 1 - score against default min prop (if prop > min prop)
                propensityRatio = fragment1.propensity > defaultPropensity ? fragment1.propensity / defaultPropensity : 1;
            }
        } else {
            if (fragment2.propensity > 0) {
                // only propensity for fragment 2 - score against default min prop as a negative for candidate 1
                propensityRatio = fragment2.propensity > defaultPropensity ? defaultPropensity / fragment2.propensity : 1;
            } else {
                // ignore if no propensity present for this fragment
                // todo: use all-glycan fragment lookup in this case?
                propensityRatio = 1;
            }
        }
        return propensityRatio;
    }

    /**
     * Helper method to compute the ratio of observed to expected intensity for a fragment.
     * NOTE: If the observed intensity is less than the "critical point" value (i.e., the product
     * of hit probability * intensity ratio < 1), the intensity ratio is set to the critical point value.
     * This makes it so that finding a low-intensity peak will have no effect on score rather than causing a
     * reduction in score.
     *
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
     *
     * @param pepMass          peptide mass - needed for correct PPM error calculation
     * @param deltaMass        delta mass being searched
     * @param glycanDatabase   list of glycan candidates
     * @param isotopesToSearch list of isotope errors to consider
     * @param ms1TolerancePPM  MS1 tolerance to consider around delta mass and isotope errors
     * @return list of glycan candidates with masses within the delta mass + iso errors and tolerance
     */
    public ArrayList<GlycanCandidateResult> getMatchingGlycansByMass(double pepMass, double deltaMass, ArrayList<GlycanCandidate> glycanDatabase, Integer[] isotopesToSearch, double ms1TolerancePPM) {
        ArrayList<GlycanCandidateResult> matchingGlycans = new ArrayList<>();
        for (int isotope : isotopesToSearch) {
            // add isotope error, which is recorded as an increase relative to delta mass
            double isotopeCorrMass = deltaMass - (isotope * AAMasses.averagineIsotopeMass) + pepMass;   // add peptide mass to get correct PPM calculation
            double massRangeDa = isotopeCorrMass * 0.000001 * ms1TolerancePPM;
            double massLo = isotopeCorrMass - massRangeDa - pepMass;    // remove pep mass after PPM calc for final calculation
            double massHi = isotopeCorrMass + massRangeDa - pepMass;
            for (GlycanCandidate glycan : glycanDatabase) {
                // see if mass within specified ranges
                if (glycan.mass >= massLo && glycan.mass <= massHi) {
                    // match. todo: check duplicates (could be if user inputs them)
                    // add copy of candidate to allow multi-threading without competing access
                    matchingGlycans.add(new GlycanCandidateResult(glycan, glycoParams.glycanResiduesMap));
                }
            }
        }
        return matchingGlycans;
    }

    /**
     * Calculate KL divergence score for MS1 spectrum match to glycan candidate using the IonQuant API.
     *
     * KL divergence near 0 is a good match, so taking the log of the 1/absolute value gives a positive score
     * that increases with better matches.
     * @param candidate glycan candidate to score
     * @param spec spectrum
     * @return score
     */
    private double calculateMS1score(GlycanCandidate candidate, Spectrum spec, double pepmass) {
        double candidateMass = candidate.mass + pepmass;    // pep mass + glycan mass
        double mz = Spectrum.calcMZ(candidateMass, spec.charge);
        Entry quantifiedEntry = api.quantXIC123((float) mz, (float) spec.rt, (float) spec.im, spec.charge, spec.cv);

        double score;
        if (quantifiedEntry == null) {
            score = 0;      // no peak found at this m/z
        } else {
            score = Math.log(1.0 / Math.abs(quantifiedEntry.kl));
        }
        return score;
    }

    /**
     * Compute unweighted spectral entropy between the two spectrum vectors, as in
     * www.nature.com/articles/s41592-021-01331-z. Assumes spectrum vectors are of equal
     * length.
     */
    private double entropyScore(double[] theoreticalPks, double[] exptPks) {
        double[] SabVector = new double[theoreticalPks.length];
        int numFrags = 0;
        for (double j : exptPks) {
            if (j != 0) {
                numFrags += 1;
            }
        }
        double[] normThyPks = normalize(theoreticalPks);
        double[] normExptPks = normalize(exptPks);

        if (numFrags < 2) {
            return 0;
        } else {
            for (int i = 0; i < SabVector.length; i++) {
                SabVector[i] = (normThyPks[i] + normExptPks[i]) / 2;
            }
        }
        return 1 - ( ((2 * spectralEntropy(SabVector)) - spectralEntropy(normExptPks) - spectralEntropy(normThyPks)) / Math.log(4));
    }

    private double spectralEntropy(double[] vector) {
        double entropy = 0;
        for (double f : vector) {
            if (f != 0) {
                entropy += (f * Math.log(f));
            }
        }
        return -1 * entropy;
    }

    /**
     * Normalize input vector so that the sum of all entries is 1
     */
    private double[] normalize(double[] vector) {
        double total = 0;
        for (double i : vector) {
            total += i;
        }
        if (total == 0) {
            return vector;
        }
        double[] output = new double[vector.length];
        for (int i=0; i < vector.length; i++) {
            output[i] = vector[i] / total;
        }
        return output;
    }


    /**
     * Generate the feature vector for LDA and the summed score for the glycan candidate result using the
     * list of features to use.
     */
    private void generateScores(GlycanCandidateResult candidate) {
        ArrayList<Double> features = new ArrayList<>();
        double sumLogRatio = 0;
        for (GlycoParams.LDAFeature feature : glycoParams.ldaFeaturesToUse) {
            switch (feature) {
                case kl: // KL score
                    features.add(candidate.ms1Score);
                    sumLogRatio += candidate.ms1Score;
                    break;
                case yprop: // Y proportion score
                    features.add(candidate.YproportionScore);
                    sumLogRatio += candidate.YproportionScore;
                    break;
                case yscore: // Y fragment score
                    features.add(candidate.YFragmentScore);
                    sumLogRatio += candidate.YFragmentScore;
                    break;
                case oxo: // Oxonium ion score
                    features.add(candidate.OxFragmentScore);
                    sumLogRatio += candidate.OxFragmentScore;
                    break;
                case mass: // Mass error score
                    features.add(candidate.massErrorScore);
                    sumLogRatio += candidate.massErrorScore;
                    break;
                case iso: // Isotope score
                    features.add(candidate.isotopeScore);
                    sumLogRatio += candidate.isotopeScore;
                    break;
                case glycanfreq:
                    if (!isFirstPass) {     // frequency can only be computed in the second pass
                        features.add(candidate.frequencyPrior);
                        sumLogRatio += candidate.frequencyPrior;
                    }
                    break;
                case ysim:
                    if (!isFirstPass) {
                        features.add(candidate.ySpecSim);
                        sumLogRatio += candidate.ySpecSim;
                    }
                    break;
                case oxsim:
                    if (!isFirstPass) {
                        features.add(candidate.oxSpecSim);
                        sumLogRatio += candidate.oxSpecSim;
                    }
                    break;
            }
        }
        candidate.summedScore = sumLogRatio;
        if (Double.isNaN(sumLogRatio)) {
            int x=0;
        }
        if (!glycoParams.glycoLDA) {
            candidate.glycanScore = candidate.summedScore;
        }
        candidate.featureVec = features.stream().mapToDouble(Double::doubleValue).toArray();
    }


    public boolean isGlycoComplete() {
        try {
            if (glycoFile.exists()) {
                RandomAccessFile raf = new RandomAccessFile(glycoFile, "r");
                raf.seek(Math.max(0, glycoFile.length() - 20));
                String cline;
                while ((cline = raf.readLine()) != null)
                    if (cline.equals("COMPLETE")) {
                        raf.close();
                        return true;
                    }
                raf.close();
                glycoFile.delete();
            }
        } catch (IOException e) {
            PTMShepherd.die("Error checking glyco file for completion: " + e.getMessage());
        }
        return false;
    }


    public void completeGlyco() {
        try {
            PrintWriter out = new PrintWriter(new FileWriter(glycoFile, true));
            out.println("COMPLETE");
            out.close();
        } catch (IOException e) {
            PTMShepherd.die("Error completing glyco file: " + e.getMessage());
        }
    }


    /**
     * Find the first location of N-glycan sequon (N-X-S/T, X is not P) in the provided peptide sequence.
     *
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
        if (pepseq_glyco[pepSeq.length() - 2] == 'N') {
            return pepSeq.length() - 2;
        } else {
            return -1;
        }
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
