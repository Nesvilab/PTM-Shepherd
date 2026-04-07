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

import edu.umich.andykong.ptmshepherd.Mod;
import edu.umich.andykong.ptmshepherd.PSM;
import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;
import ionquant.api.Entry;
import ionquant.api.IonQuantAPI;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.hipparchus.stat.descriptive.rank.Median;
import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanCandidate;
import umich.ms.glyco.GlycanFragment;
import umich.ms.util.AminoAcids;
import umich.ms.util.ElementalComposition;
import precise.IsotopeDistributionApi;

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
    public ArrayList<GlycanCandidate> glycanDatabase;
    LinkedHashMap<String, GlycanCandidate> glycanDBmap;
    Double meanMassError;
    double massErrorWidth;
    public static final double DEFAULT_GLYCO_PPM_TOL = 30;
    public static final double DEFAULT_GLYCO_FDR = 0.01;
    public static final int DEFAULT_GLYCO_DECOY_TYPE = 1;
    public static final String GLYCAN_COMP_COL_NAME = "Total Glycan Composition";
    public LinkedHashMap<Integer, LinkedHashMap<String, Integer>> glycanMassBinMap;
    public static final double DEFAULT_GLYCO_PROPENSITY = 0.1;
    public static final double MIN_SIMILARITY = 0.01;
    private final GlycoParams glycoParams;
    public final ArrayList<GlycanAssignmentResult> allResults;
    public LinkedHashMap<String, ArrayList<GlycanCandidateResult>> highConfidenceResultMap;
    private ArrayList<GlycanCandidateResult> allHighConfidenceResults;
    public LinkedHashMap<String, GlycanCandidateFragments> targetGlycanFragmentProps;
    public LinkedHashMap<String, GlycanCandidateFragments> decoyGlycanFragmentProps;
    private final String ldaHeader;
    private static IonQuantAPI api;
    public final boolean isFirstPass;
    private static final float ISOTOPE_MASS_DIFF = 1.00235f;
    private static final int MAX_ISOTOPE_PEAKS = 6;
    public static final IsotopeDistributionApi isotopeDistributionApi = new IsotopeDistributionApi();

    // Default constructor
    public GlycoAnalysis(String dsName, ArrayList<GlycanCandidate> glycoDatabase, GlycoParams glycoParams, boolean isFirstPass) {
        this.dsName = dsName;
        this.isFirstPass = isFirstPass;
        String firstPassName = isFirstPass ? PTMShepherd.rawGlycoFirstPass : "";
        this.glycoFile = new File(PTMShepherd.normFName(dsName + firstPassName + PTMShepherd.rawGlycoName));
        this.glycanDatabase = glycoDatabase;
        glycanDBmap = new LinkedHashMap<>();
        for (GlycanCandidate glycan : glycanDatabase) {
            if (!glycan.isDecoy) {
                String glycanHash = Glycan.toGlycanString(glycan.composition);
                glycanDBmap.put(glycanHash, glycan);
            }
        }
        this.glycoParams = glycoParams;
        this.glycanMassBinMap = new LinkedHashMap<>();
        this.allResults = new ArrayList<>();
        this.highConfidenceResultMap = new LinkedHashMap<>();
        this.targetGlycanFragmentProps = new LinkedHashMap<>();
        this.decoyGlycanFragmentProps = new LinkedHashMap<>();
        ldaHeader = glycoParams.glycoLDA ? glycoParams.generateLDAheader() : "\t";
    }

    public void glycoPSMs(PSMFile psmFile,
                          HashMap<String, File> mzMappings,
                          HashMap<String, File> originalMzMappings,
                          ExecutorService executorService) {

        //open up output file
        LinkedHashMap<String, ArrayList<Integer>> mappings = new LinkedHashMap<>();
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
        if (params.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.kl) || params.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ms1) || params.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ms1delta)) {
            PTMShepherd.print("\tBuilding IonQuant index for " + filePath);
            api = new IonQuantAPI(
                    filePath,
                    params.numThreads,
                    (float) params.glycoPPMtol,
                    params.rtTol,
                    params.imTol,
                    params.minIsotopesIonQuant,
                    params.minScansIonQuant,
                    !filePath.toLowerCase().endsWith(".d")      // auto-detect IM data for the "noPASEF" parameter
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
     * Converged when glycan database stops changing size between passes.
     * @param passNum current pass number
     * @param newGlycanDB new glycan database generated from prev pass
     * @return true if converged
     */
    public boolean checkConvergence(int passNum, ArrayList<GlycanCandidate> newGlycanDB) {
        if (passNum == 1) {
            return false;   // always run at least 2 passes
        }
        return newGlycanDB.size() == glycanDatabase.size();
    }

    /**
     * Generate decoy glycan candidates for 2nd pass using glycan fragment information from 1st pass.
     * Decoys are generated using spectra/masses from alternate possible composition matches for a target.
     * Called after target glycan spectra are initialized (computeGlycanFragmentProbs).
     */
    public void generateDecoys2ndPass() {
        for (String glycanKey : targetGlycanFragmentProps.keySet()) {
            GlycanCandidate target = glycanDBmap.get(glycanKey);
            GlycanCandidateFragments decoyFragmentInfo = getGlycanCandidateFragmentsRandom(target, allHighConfidenceResults);
            decoyGlycanFragmentProps.put(glycanKey, decoyFragmentInfo);
        }
    }

    /**
     * Process summarized glycan fragment information from the first pass to compute fragment intensities
     * for each glycan as a target and decoy. Target fragments are from filtered (FDR passed, sufficient PSMs, at least
     * min Y fragments) glycans. Decoy fragments are from all PSMs not passing those filters.
     * If insufficient target PSMs, the glycan is not considered in the second pass.
     * If insufficient decoy PSMs, the lowest scoring target PSMs for that glycan are used to fill up to the min PSMs.
     */
    public void computeGlycanFragmentProbs() {
        // filter target glycans to only those with sufficient PSMs
        LinkedHashMap<String, ArrayList<GlycanCandidateResult>> filteredTargetInputs = highConfidenceResultMap.entrySet().stream()
                .filter(e -> e.getValue().size() > glycoParams.minPSMsForConsensus)
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        Map.Entry::getValue,
                        (v1, v2) -> v1,  // merge function (not needed here, but required)
                        LinkedHashMap::new  // supplier to create LinkedHashMap
                ));

        // Generate target fragment intensity profiles
        for (Map.Entry<String, ArrayList<GlycanCandidateResult>> glycanEntry : filteredTargetInputs.entrySet()) {
            // Determine the fragment likelihoods based on all PSMs for this entry
            ArrayList<GlycanCandidateResult> targetGlycoPSMs = glycanEntry.getValue();
            // skip generating fragment information for glycans with too few PSMs to get reasonable values
            if (targetGlycoPSMs.size() < glycoParams.minPSMsForConsensus) {
                continue;
            }
            if (glycoParams.topPctSpectraForConsensus < 1.0) {
                // sort PSMs by glycan score and keep only the top X%
                targetGlycoPSMs.sort(Comparator.comparingDouble((GlycanCandidateResult result) -> result.glycanScore).reversed());
                int numToKeep = (int) Math.ceil(targetGlycoPSMs.size() * glycoParams.topPctSpectraForConsensus);
                if (numToKeep > glycoParams.minPSMsForConsensus) {
                    targetGlycoPSMs = new ArrayList<>(targetGlycoPSMs.subList(0, numToKeep));
                }
            }
            GlycanCandidateFragments fragmentInfo = getGlycanCandidateFragments(targetGlycoPSMs);
            targetGlycanFragmentProps.put(glycanEntry.getKey(), fragmentInfo);
        }

        // generate decoy fragment intensity profiles
        generateDecoys2ndPass();
    }

    /**
     * Compute glycan candidate fragment intensities from a list of PSMs for that glycan.
     * @param glycanPSMlist list of PSMs for the glycan
     * @return container with computed fragment intensities
     */
    private GlycanCandidateFragments getGlycanCandidateFragments(ArrayList<GlycanCandidateResult> glycanPSMlist) {
        LinkedHashMap<String, Double> yFragmentIntensities;
        LinkedHashMap<String, Double> OxFragmentIntensities;
        LinkedHashMap<String, Double> generalOxFragmentIntensities;

        // use fragment intensities from all PSMs and compute avg or median. Assumes all PSMs are of the same glycan.
        LinkedHashMap<String, ArrayList<Double>> YInts = new LinkedHashMap<>();
        LinkedHashMap<String, ArrayList<Double>> OxInts = new LinkedHashMap<>();
        LinkedHashMap<String, ArrayList<Double>> generalOxInts = new LinkedHashMap<>();
        for (GlycanCandidate inputGlycan : glycanPSMlist) {
            // read all fragments from each input glycan into the intensity lists
            for (String fragmentHash : inputGlycan.Yfragments.keySet()) {
                YInts.computeIfAbsent(fragmentHash, k -> new ArrayList<>()).add(inputGlycan.Yfragments.get(fragmentHash).foundIntensity);
            }
            for (String fragmentHash : inputGlycan.oxoniumFragments.keySet()) {
                OxInts.computeIfAbsent(fragmentHash, k -> new ArrayList<>()).add(inputGlycan.oxoniumFragments.get(fragmentHash).foundIntensity);
            }
            for (String fragmentHash : inputGlycan.generalOxoniumFragments.keySet()) {
                generalOxInts.computeIfAbsent(fragmentHash, k -> new ArrayList<>()).add(inputGlycan.generalOxoniumFragments.get(fragmentHash).foundIntensity);
            }
        }
        // combine intensities
        if (glycoParams.glycoAvgInts) {
            yFragmentIntensities = calculateFragmentAvgInts(YInts);
            OxFragmentIntensities = calculateFragmentAvgInts(OxInts);
            generalOxFragmentIntensities = calculateFragmentAvgInts(generalOxInts);
        } else {
            yFragmentIntensities = calculateFragmentMedianInts(YInts);
            OxFragmentIntensities = calculateFragmentMedianInts(OxInts);
            generalOxFragmentIntensities = calculateFragmentMedianInts(generalOxInts);
        }

        // save determined propensities to the output container
        return new GlycanCandidateFragments(yFragmentIntensities, OxFragmentIntensities, generalOxFragmentIntensities);
    }

    /**
     * Compute glycan candidate fragment intensities from a list of PSMs for that glycan by randomly selecting
     * fragment intensities from the input PSMs.
     * @param targetGlycan the glycan candidate to generate fragments for
     * @param PSMlist list of PSMs to pull fragment intensities from
     * @return
     */
    private GlycanCandidateFragments getGlycanCandidateFragmentsRandom(GlycanCandidate targetGlycan, ArrayList<GlycanCandidateResult> PSMlist) {
        LinkedHashMap<String, Double> yFragmentIntensities = new LinkedHashMap<>();
        LinkedHashMap<String, Double> OxFragmentIntensities = new LinkedHashMap<>();
        LinkedHashMap<String, Double> generalOxFragmentIntensities = new LinkedHashMap<>();

        // Pre-filter PSMs to exclude target glycan - this is much faster than checking in every iteration
        ArrayList<GlycanCandidateResult> eligiblePSMs = new ArrayList<>();
        for (GlycanCandidateResult psm : PSMlist) {
            if (!psm.composition.equals(targetGlycan.composition)) {
                eligiblePSMs.add(psm);
            }
        }

        // If no eligible PSMs, return empty fragment info
        if (eligiblePSMs.isEmpty()) {
            return new GlycanCandidateFragments(yFragmentIntensities, OxFragmentIntensities, generalOxFragmentIntensities);
        }

        final int MAX_RETRIES = 100;  // Prevent infinite loops

        // randomly select fragment intensities from input PSMs. Input PSMs do not need to be of the same glycan.
        for (String fragmentHash : targetGlycan.Yfragments.keySet()) {
            double intensity = -1;
            int retries = 0;
            while (intensity == -1 && retries < MAX_RETRIES) {
                GlycanCandidateResult randomPSM = eligiblePSMs.get(glycoParams.randomGenerator.nextInt(eligiblePSMs.size()));
                if (randomPSM.Yfragments.containsKey(fragmentHash)) {
                    intensity = randomPSM.Yfragments.get(fragmentHash).foundIntensity;
                    yFragmentIntensities.put(fragmentHash, intensity);
                }
                retries++;
            }
            // If no match found after max retries, use intensity of 0
            if (intensity == -1) {
                yFragmentIntensities.put(fragmentHash, 0.0);
            }
        }
        for (String fragmentHash : targetGlycan.oxoniumFragments.keySet()) {
            double intensity = -1;
            int retries = 0;
            while (intensity == -1 && retries < MAX_RETRIES) {
                GlycanCandidateResult randomPSM = eligiblePSMs.get(glycoParams.randomGenerator.nextInt(eligiblePSMs.size()));
                if (randomPSM.oxoniumFragments.containsKey(fragmentHash)) {
                    intensity = randomPSM.oxoniumFragments.get(fragmentHash).foundIntensity;
                    OxFragmentIntensities.put(fragmentHash, intensity);
                }
                retries++;
            }
            if (intensity == -1) {
                OxFragmentIntensities.put(fragmentHash, 0.0);
            }
        }
        for (String fragmentHash : targetGlycan.generalOxoniumFragments.keySet()) {
            double intensity = -1;
            int retries = 0;
            while (intensity == -1 && retries < MAX_RETRIES) {
                GlycanCandidateResult randomPSM = eligiblePSMs.get(glycoParams.randomGenerator.nextInt(eligiblePSMs.size()));
                if (randomPSM.generalOxoniumFragments.containsKey(fragmentHash)) {
                    intensity = randomPSM.generalOxoniumFragments.get(fragmentHash).foundIntensity;
                    generalOxFragmentIntensities.put(fragmentHash, intensity);
                }
                retries++;
            }
            if (intensity == -1) {
                generalOxFragmentIntensities.put(fragmentHash, 0.0);
            }
        }

        // save determined propensities to the output container
        return new GlycanCandidateFragments(yFragmentIntensities, OxFragmentIntensities, generalOxFragmentIntensities);
    }

    private void shuffleValues(HashMap<String, Double> map) {
        List<String> keys = new ArrayList<>(map.keySet());
        Collections.sort(keys);     // Ensure consistent key order so that shuffling is reproducible

        List<Double> values = new ArrayList<>(map.values());
        Collections.shuffle(values, glycoParams.randomGenerator);

        for (int i = 0; i < keys.size(); i++) {
            map.put(keys.get(i), values.get(i));
        }
    }

    // compute the avg intensity for each fragment using the associated intensity list
    private static LinkedHashMap<String, Double> calculateFragmentAvgInts(HashMap<String, ArrayList<Double>> intensityList) {
        LinkedHashMap<String, Double> fragmentAvgInts = new LinkedHashMap<>();
        for (Map.Entry<String, ArrayList<Double>> fragmentEntry : intensityList.entrySet()) {
            double sum = 0;
            for (int i = 0; i < fragmentEntry.getValue().size(); i++) {
                sum += fragmentEntry.getValue().get(i);
            }
            fragmentAvgInts.put(fragmentEntry.getKey(), sum / (double) fragmentEntry.getValue().size());
        }
        return fragmentAvgInts;
    }

    // compute the median intensity for each fragment using the associated intensity list
    private static LinkedHashMap<String, Double> calculateFragmentMedianInts(HashMap<String, ArrayList<Double>> intensityList) {
        LinkedHashMap<String, Double> fragmentMedianInts = new LinkedHashMap<>();
        for (Map.Entry<String, ArrayList<Double>> fragmentEntry : intensityList.entrySet()) {
            double[] intensities = new double[fragmentEntry.getValue().size()];
            for (int i = 0; i < fragmentEntry.getValue().size(); i++) {
                intensities[i] = fragmentEntry.getValue().get(i);
            }
            Median median = new Median();
            fragmentMedianInts.put(fragmentEntry.getKey(), median.evaluate(Arrays.stream(intensities).toArray()));
        }
        return fragmentMedianInts;
    }

    public void summarizeGlycanResults() {
        // read all glycan info in
        for (GlycanAssignmentResult result : allResults) {
            if (result.foundGlycan) {
                GlycanCandidateResult glycan = result.bestCandidate;
                GlycanCandidate fragmentInfoContainer = new GlycanCandidate(glycan.composition, 0, false, glycoParams.glycanResiduesMap, glycan.Yfragments, glycan.oxoniumFragments, glycan.generalOxoniumFragments);

                String glycanHash = Glycan.toGlycanString(fragmentInfoContainer.composition);
                // only include good targets in fragment info
                if (glycoParams.noFDR || result.glycanQval < glycoParams.glycoFDR) {
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
                    addGlycanToMap(highConfidenceResultMap, glycanHash, glycan);

                    // add to delta mass map for calculating glycan prevalence priors (targets and decoys)
                    double deltaMass = result.deltaMass;
                    int massBin = (int) Math.floor(deltaMass);
                    if (glycanMassBinMap.containsKey(massBin)) {
                        // seen this mass bin before. Get the count-by-glycan dict and increment the count for this glycan
                        LinkedHashMap<String, Integer> massBinGlycanCounts = glycanMassBinMap.get(massBin);
                        int glycanCount = massBinGlycanCounts.getOrDefault(glycanHash, 0);
                        glycanCount++;
                        massBinGlycanCounts.put(glycanHash, glycanCount);
                    } else {
                        // New mass bin. Create a new count-by-glycan dict
                        LinkedHashMap<String, Integer> massBinGlycanCounts = new LinkedHashMap<>();
                        massBinGlycanCounts.put(glycanHash, 1);
                        glycanMassBinMap.put(massBin, massBinGlycanCounts);
                    }
                }
            }
        }
        // filter high confidence results to only those with sufficient PSMs
        highConfidenceResultMap = highConfidenceResultMap.entrySet().stream()
                .filter(e -> e.getValue().size() >= glycoParams.minPSMsForConsensus)
                .collect(Collectors.toMap(
                        Map.Entry::getKey,
                        Map.Entry::getValue,
                        (v1, v2) -> v1,
                        LinkedHashMap::new
                ));
        allHighConfidenceResults = new ArrayList<>();
        for (ArrayList<GlycanCandidateResult> glycanList : highConfidenceResultMap.values()) {
            allHighConfidenceResults.addAll(glycanList);
        }
    }

    private static void addGlycanToMap(LinkedHashMap<String, ArrayList<GlycanCandidateResult>> targetInputGlycans, String glycanHash, GlycanCandidateResult glycan) {
        if (targetInputGlycans.containsKey(glycanHash)) {
            targetInputGlycans.get(glycanHash).add(glycan);
        } else {
            ArrayList<GlycanCandidateResult> newList = new ArrayList<>();
            newList.add(glycan);
            targetInputGlycans.put(glycanHash, newList);
        }
    }

    /**
     * Glycan FDR wrapper/main method. Called after the glycan assignment is done on PSMs to compute glycan FDR.
     * Handles various old and new methods and LDA.
     */
    public void runScoresAndFDR() {
        // LDA method
        if (glycoParams.glycoLDA && !isFirstPass) {
            ScoreLDA lda = new ScoreLDA();
            if (glycoParams.noFDR) {
                // In no-FDR mode there are no decoy glycans; use low-scoring target results as decoy training data
                for (GlycanAssignmentResult result : allResults) {
                    if (result.foundGlycan && result.bestTarget != null) {
                        lda.targetData.add(result.bestTarget.featureVec);
                    }
                }
                // Sort by sum of feature scores ascending so the lowest-scoring results come first
                List<double[]> sortedByScore = new ArrayList<>(lda.targetData);
                sortedByScore.sort(Comparator.comparingDouble(fv -> Arrays.stream(fv).sum()));
                // Bottom (1 - ldaTargetProp) fraction becomes the decoy training data
                int decoyCount = (int) Math.floor(sortedByScore.size() * (1.0 - glycoParams.ldaTargetProp));
                for (int i = 0; i < decoyCount; i++) {
                    lda.decoyData.add(sortedByScore.get(i));
                }
            } else {
                // Standard mode: use decoy glycan results as decoy training data
                for (GlycanAssignmentResult result : allResults) {
                    if (result.foundGlycan) {
                        if (result.bestTarget != null) {
                            lda.targetData.add(result.bestTarget.featureVec);
                        }
                        if (result.bestDecoy != null) {
                            lda.decoyData.add(result.bestDecoy.featureVec);
                        }
                    }
                }
            }
            lda.runLDA(allResults, ldaHeader, glycoParams.ldaTargetProp);
        }

        // Compute FDR (skipped in no-FDR mode)
        boolean fdrSuccess = true;
        if (!glycoParams.noFDR) {
            PTMShepherd.print("\tCalculating Glycan FDR");
            fdrSuccess = computeFDRcompetitive(allResults, glycoParams.glycoFDR, glycoParams.numDecoysPerTarget);
            if (!fdrSuccess) {
                fdrSuccess = computeFDRNonCompetitive(allResults, glycoParams.glycoFDR, glycoParams.numDecoysPerTarget);
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
    public static boolean computeFDRcompetitive(List<GlycanAssignmentResult> results, double fdrCutOff, int numDecoysPerTarget) {
        // Sort target scores in descending order
        results.sort(Comparator.comparingDouble((GlycanAssignmentResult result) -> result.glycanScore).reversed());
        long glycoResultCount = results.stream().filter(r -> r.foundGlycan).count();

        // Get decoy indices in the combined sorted list
        List<Integer> decoyIndexes = getDecoyIndexes(results);

        int decoyCount = decoyIndexes.size();
        int targetCount = (int) glycoResultCount - decoyCount;
        // check if enough decoys were found (i.e., initial q-val is above the desired threshold)
        double initialFDR = calculateFDR(targetCount, decoyCount, false, numDecoysPerTarget);
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
            double fdr = Math.min(calculateFDR(targetCount, decoyCount, false, numDecoysPerTarget), currentMinQ);        // q = (d+1)/t recommended per 10.1021/acs.jproteome.6b00144
            if (fdr < currentMinQ) {
                currentMinQ = fdr;
            }
            if (!foundThreshold) {
                if (fdr < fdrCutOff) {
//                    scoreThreshold = results.get(i).glycanScore;
                    PTMShepherd.print(String.format("\tFound glycan score threshold: %.2f with %d decoys, %d targets for %.2f%% estimated FDR (%d total inputs)",
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
    public static boolean computeFDRNonCompetitive(List<GlycanAssignmentResult> results, double glycoFDR, int numDecoysPerTarget) {
        LinkedHashMap<String, GlycanAssignmentResult> resultMap = new LinkedHashMap<>();

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
            targetDecoyRatio = calculateFDR(targets, decoys, true, numDecoysPerTarget);
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
     * Calculate FDR from target and decoy counts, corrected for the number of decoys per target in the database.
     * The dbRatio (database targets / database decoys = 1 / numDecoysPerTarget) scales the observed decoy count
     * to account for inflated decoy numbers when multiple decoys per target are used.
     *
     * @param targets target count
     * @param decoys  decoy count
     * @param useNonCompFDR true for non-competitive FDR formula, false for competitive
     * @param numDecoysPerTarget number of decoys generated per target glycan
     * @return FDR
     */
    private static double calculateFDR(int targets, int decoys, boolean useNonCompFDR, int numDecoysPerTarget) {
        double dbRatio = 1.0 / numDecoysPerTarget;
        if (useNonCompFDR) {
            return dbRatio * (2.0 * decoys) / (decoys + targets);
        } else {
            return dbRatio * (decoys + 1) / (double) targets;
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
            float deltaMass = (float) psmFile.psms.get(cline).getDMass();

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
        LinkedHashMap<Integer, Integer> isotopeCounts = new LinkedHashMap<>();
        double minError = 10;
        double maxError = -10;

        for (GlycanAssignmentResult result : results) {
            // only use target glyco PSMs that passed FDR (or all targets in no-FDR mode)
            if (result.foundGlycan && !result.isDecoyGlycan && (glycoParams.noFDR || result.glycanQval < glycoParams.glycoFDR)) {
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
        computeMassErrorsHelper(massErrors, maxError, minError);
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
        GlycanAssignmentResult glycoResult = new GlycanAssignmentResult(psm.lineNum, psm.getPeptide(), (float) psm.getDMass(), (float) (double) psm.getCalcPepmass(), psm.printAssignedMods(), psm.getAssignedMods(), psm.getSpec());

        // read spectrum and condition
        Spectrum spec = mr.getSpectrum(psm.getSpec());
        if (spec == null) {
            this.lineWithoutSpectra.add(psm.getSpec());
            psm.glycanAssignmentResult = glycoResult;
            return;
        }
        spec.conditionOptNorm(condPeaks, condRatio, false);

        // do glycan assignment with original delta mass
        glycoResult = assignGlycanToPSM(spec, glycoResult, glycanDatabase, massErrorWidth, meanMassError);

        // optionally test removing each assigned variable mod to see if a better glycan assignment is possible
        if (!isFirstPass && glycoParams.checkVariableMods && psm.getAssignedMods() != null && psm.getDMass() > 3.5) {
            for (Mod mod : psm.getAssignedMods()) {
                // compute alternative masses: removing mod adds its mass to delta, subtracts from pep mass
                float altDeltaMass = (float) (psm.getDMass() + mod.mass);
                float altPepMass = (float) (psm.getCalcPepmass() - mod.mass);

                GlycanAssignmentResult altResult = new GlycanAssignmentResult(psm.lineNum, psm.getPeptide(), altDeltaMass, altPepMass, psm.printAssignedMods(), psm.getAssignedMods(), psm.getSpec());
                altResult = assignGlycanToPSM(spec, altResult, glycanDatabase, massErrorWidth, meanMassError);

                // compare by glycanScore (summed score before LDA) - keep the better result
                if (altResult.foundGlycan && altResult.glycanScore > glycoResult.glycanScore) {
                    // format mod description
                    String posStr;
                    if (mod.position == 0) {
                        posStr = "N-term";
                    } else if (mod.position == psm.getPeptide().length() + 1) {
                        posStr = "C-term";
                    } else {
                        posStr = mod.position + "" + psm.getPeptide().charAt(mod.position - 1);
                    }
                    altResult.modChangeDescription = String.format("removed %s(%+.4f)", posStr, mod.mass);
                    altResult.removedMod = mod;
                    glycoResult = altResult;
                }
            }
        }

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
            for (GlycanCandidateResult candidate : searchCandidates) {
                matchFragmentsToSpectra(spec, glycoResult, candidate, ppmTol);
            }

            // score candidates and save results
            int bestCandidateIndex = 0;
            double[] scoresVsBestCandidate = new double[searchCandidates.size()];
            for (int i = 0; i < searchCandidates.size(); i++) {
                if (i == bestCandidateIndex) {
                    continue;
                }
                double comparisonScore;
                if (!isFirstPass) {
                    comparisonScore = pairwiseCompare2ndPass(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult, spec);
                } else {
                    comparisonScore = pairwiseCompare1stPass(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult, spec);
                }

                if (comparisonScore == 0) {
                    // exact same score (e.g., from target/decoy if no Y/oxo ions found and using decoy mass = target mass)
                    // Use a deterministic tiebreaker based on candidate composition strings.
                    // NOTE: do NOT use glycoParams.randomGenerator here - this method runs in parallel threads,
                    // and concurrent access to Random causes non-deterministic results even with a fixed seed.
                    String compositionA = Glycan.toGlycanString(searchCandidates.get(bestCandidateIndex).composition);
                    String compositionB = Glycan.toGlycanString(searchCandidates.get(i).composition);
                    int cmp = compositionA.compareTo(compositionB);
                    comparisonScore = (cmp != 0 ? Math.signum(cmp) : 1.0) * 1E-6;
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
                if (!isFirstPass) {
                    scoresVsBestCandidate[i] = pairwiseCompare2ndPass(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult, spec);
                } else {
                    scoresVsBestCandidate[i] = pairwiseCompare1stPass(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), glycoResult, spec);
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
            if (!isFirstPass) {
                computeAbsoluteScore2ndPass(spec, searchCandidates.get(bestCandidateIndex), glycoResult);
            } else {
                computeAbsoluteScore1stPass(spec, searchCandidates.get(bestCandidateIndex), glycoResult);
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
        } else {
            // no glycan candidates found for this delta mass - filter out of results
            glycoResult.foundGlycan = false;
        }

        return glycoResult;
    }

    private static void matchFragmentsToSpectra(Spectrum spec, GlycanAssignmentResult glycoResult, GlycanCandidateResult candidate, float ppmTol) {
        // Search all Y fragment ions in a single spectrum pass
        List<GlycanFragment> yFragList = new ArrayList<>(candidate.Yfragments.values());
        double[] yMasses = new double[yFragList.size()];
        for (int i = 0; i < yFragList.size(); i++) {
            yMasses[i] = yFragList.get(i).neutralMass + glycoResult.pepMass;
        }
        double[] yIntensities = spec.findIonsNeutral(yMasses, ppmTol, Integer.parseInt(PTMShepherd.getParam("spectra_maxPrecursorCharge")));
        for (int i = 0; i < yFragList.size(); i++) {
            yFragList.get(i).foundIntensity = yIntensities[i] / spec.basePeakInt;
        }

        // Search all oxonium and general oxonium ions in a single spectrum pass
        List<GlycanFragment> oxoFragList = new ArrayList<>(candidate.oxoniumFragments.values());
        List<GlycanFragment> genOxoFragList = new ArrayList<>(candidate.generalOxoniumFragments.values());
        double[] oxoMasses = new double[oxoFragList.size() + genOxoFragList.size()];
        for (int i = 0; i < oxoFragList.size(); i++) {
            oxoMasses[i] = oxoFragList.get(i).neutralMass + AAMasses.protMass;
        }
        for (int j = 0; j < genOxoFragList.size(); j++) {
            oxoMasses[oxoFragList.size() + j] = genOxoFragList.get(j).neutralMass + AAMasses.protMass;
        }
        double[] oxoIntensities = spec.findIons(oxoMasses, ppmTol);
        for (int i = 0; i < oxoFragList.size(); i++) {
            oxoFragList.get(i).foundIntensity = oxoIntensities[i] / spec.basePeakInt;
        }
        for (int j = 0; j < genOxoFragList.size(); j++) {
            genOxoFragList.get(j).foundIntensity = oxoIntensities[oxoFragList.size() + j] / spec.basePeakInt;
        }

        // normalize intensities for Y and oxonium ions separately, so that the ratio of Y to oxo (or unfragmented precursor/etc) does not impact scores
        GlycanCandidateResult.normalizeIntensities(candidate.Yfragments);
        GlycanCandidateResult.normalizeIntensities(candidate.generalOxoniumFragments);
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
            if (!isFirstPass) {
                computeAbsoluteScore2ndPass(spec, nextCandidate, glycoResult);
            } else {
                computeAbsoluteScore1stPass(spec, nextCandidate, glycoResult);
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
                // ensure same targets and decoys match by removing the decoy mass shift (mass shift is for scoring only). Target mass shift = 0
                double correctedMass = glycan.mass - glycan.decoyMassShift;
                if (correctedMass >= massLo && correctedMass <= massHi) {
                    // add copy of candidate to allow multi-threading without competing access
                    matchingGlycans.add(new GlycanCandidateResult(glycan, glycoParams.glycanResiduesMap));
                }
            }
        }
        return matchingGlycans;
    }

    // ---------------------------------------------------------------------------------------------
    // 1st pass scoring methods
    // ---------------------------------------------------------------------------------------------

    /**
     * Perform pairwise comparison of two glycans. Uses sum of log probability ratios between candidates for
     * each category (mass/iso error and fragment ion) being considered. Returns a single score of combined
     * probability of first glycan candidate over second.
     *
     * @param glycan1   candidate 1
     * @param glycan2   candidate 2
     * @return output probability score (sum of log ratios)
     */
    public double pairwiseCompare1stPass(GlycanCandidateResult glycan1, GlycanCandidateResult glycan2, GlycanAssignmentResult glycoResult, Spectrum spec) {
        double sumLogRatio = 0;
        // Y ions
        sumLogRatio += pairwiseScoreY1stPass(glycan1, glycan2);

        // oxonium ions
        sumLogRatio += pairwiseScoreOxo1stPass(glycan1, glycan2);

        // isotope and mass errors
        double massScore1 = Math.max(normedMassScore(glycan1, glycoResult), MIN_SIMILARITY);  // prevent extreme values from divide by zero/small value
        double massScore2 = Math.max(normedMassScore(glycan2, glycoResult), MIN_SIMILARITY);
        sumLogRatio += Math.log(massScore1 / massScore2);

        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.kl)) {
            double ms1score1 = klMS1score(glycan1, spec, glycoResult.pepMass);
            double ms1score2 = klMS1score(glycan2, spec, glycoResult.pepMass);
            sumLogRatio += (ms1score1 - ms1score2);
        }
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ms1) || glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ms1delta)) {
            double ms1score1 = ms1CentroidScore(glycan1, spec, glycoResult);
            double ms1score2 = ms1CentroidScore(glycan2, spec, glycoResult);
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
    public double pairwiseScoreY1stPass(GlycanCandidate glycan1, GlycanCandidate glycan2) {
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
        sumLogRatio += Math.sqrt(cand1Hits) * Math.log(cand1HitProb);       // candidate 1 hits - added
        sumLogRatio -= Math.sqrt(cand2Hits) * Math.log(cand2HitProb);       // candidate 2 hits - subtracted
        sumLogRatio += Math.sqrt(cand1Misses) * Math.log(cand1MissProb);    // candidate 1 misses - negative value added
        sumLogRatio -= Math.sqrt(cand2Misses) * Math.log(cand2MissProb);    // candidate 2 misses - negative value subtracted
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
    public double pairwiseScoreOxo1stPass(GlycanCandidate glycan1, GlycanCandidate glycan2) {
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
     * Compute the "absolute" score of this glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     *
     * @param candidate     glycan candidate to calculate score for
     * @param result         result container
     */
    public void computeAbsoluteScore1stPass(Spectrum spec, GlycanCandidateResult candidate, GlycanAssignmentResult result) {
        candidate.YFragmentScore = absoluteYScore1stPass(candidate);
        candidate.OxFragmentScore = absoluteOxoScore(candidate);
        candidate.massErrorScore = normedMassScore(candidate, result);

        // only calculate MS1 score if requested because it requires slow index building
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.kl)) {
            candidate.klScore = klMS1score(candidate, spec, result.pepMass);
        }
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ms1) || glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ms1delta)) {
            candidate.ms1Score = ms1CentroidScore(candidate, spec, result);
        }

        generateScores(candidate);
    }

    /**
     * Compute the "absolute" score of the Y ions of the provided glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     * Normalized to account for different glycan sizes. Not compatible with fragment-specific probabilities
     *
     * @param bestGlycan glycan candidate to calculate score for
     * @return absolute score
     */
    public double absoluteYScore1stPass(GlycanCandidate bestGlycan) {
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
    public double absoluteOxoScore(GlycanCandidate bestGlycan) {
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

    public double pairwiseCompare2ndPass(GlycanCandidateResult glycan1, GlycanCandidateResult glycan2, GlycanAssignmentResult glycoResult, Spectrum spec) {
        // calculate fragment-specific prob estimates based on observed fragment ions
        double sumLogRatio = 0;
        // Y ions (similarity score)
        double ySim1 = similarityScore(new ArrayList<>(glycan1.Yfragments.values()));
        double ySim2 = similarityScore(new ArrayList<>(glycan2.Yfragments.values()));
        sumLogRatio += Math.log(ySim1 / ySim2);

        // general oxonium ions (similarity score)
        double oxSim1 = similarityScore(new ArrayList<>(glycan1.generalOxoniumFragments.values()));
        double oxSim2 = similarityScore(new ArrayList<>(glycan2.generalOxoniumFragments.values()));
        sumLogRatio += Math.log(oxSim1 / oxSim2);

        // mass error score
        double massScore1 = Math.max(normedMassScore(glycan1, glycoResult), MIN_SIMILARITY);  // prevent extreme values from divide by zero/small value
        double massScore2 = Math.max(normedMassScore(glycan2, glycoResult), MIN_SIMILARITY);
        sumLogRatio += Math.log(massScore1 / massScore2);

        // MS1 score
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.kl)) {
            double ms1score1 = klMS1score(glycan1, spec, glycoResult.pepMass);
            double ms1score2 = klMS1score(glycan2, spec, glycoResult.pepMass);
            sumLogRatio += (ms1score1 - ms1score2);
        }
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ms1) || glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ms1delta)) {
            double ms1score1 = ms1CentroidScore(glycan1, spec, glycoResult);
            double ms1score2 = ms1CentroidScore(glycan2, spec, glycoResult);
            sumLogRatio += (ms1score1 - ms1score2);
        }

        // todo: glycan freq score for pairwise? I think it's good not to bias the comparison, but this is to remember that it is NOT used in pairwise

        return sumLogRatio;
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
     * Compute the "absolute" score of this glycan for the given spectrum, meaning the score if all ions are distinguishing
     * (i.e. the sum total evidence for/against this glycan, not relative to another glycan).
     */
    public void computeAbsoluteScore2ndPass(Spectrum spec, GlycanCandidateResult candidate, GlycanAssignmentResult result) {
        // Y ions
        candidate.ySpecSim = similarityScore(new ArrayList<>(candidate.Yfragments.values()));

        // oxonium ions
        candidate.OxFragmentScore = similarityScore(new ArrayList<>(candidate.oxoniumFragments.values()));
        candidate.oxSpecSim = similarityScore(new ArrayList<>(candidate.generalOxoniumFragments.values()));

        // mass error score
        candidate.massErrorScore = normedMassScore(candidate, result);

        // only calculate MS1 score if requested because it requires slow index building
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.kl)) {
            candidate.klScore = klMS1score(candidate, spec, result.pepMass);
        }
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ms1) || glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.ms1delta)) {
            candidate.ms1Score = ms1CentroidScore(candidate, spec, result);
        }
        if (glycoParams.ldaFeaturesToUse.contains(GlycoParams.LDAFeature.glycanfreq)) {
            candidate.frequencyPrior = glycanFrequencyScore(candidate, result.deltaMass);
        }
        generateScores(candidate);
    }

    /**
     * Compute similarity score between found and expected intensities for a set of fragments. Requires that
     * the fragment found and expected intensities are already populated.
     * @param fragments list of fragments to compute similarity
     * @return (double) similarity score
     */
    private double similarityScore(ArrayList<GlycanFragment> fragments) {
        double[] foundYs = new double[fragments.size()];
        double[] expectedYs = new double[fragments.size()];
        boolean foundNonZero = false;
        for (int i = 0; i < fragments.size(); i++) {
            foundYs[i] = fragments.get(i).foundIntensity;
            if (foundYs[i] > 0) {
                foundNonZero = true;
            }
            expectedYs[i] = fragments.get(i).expectedIntensity;
        }
        if (!foundNonZero) {
            return MIN_SIMILARITY; // no matching ions found, return minimum similarity
        }

        double score;
        if (glycoParams.cosineSimilarityScoring) {
            score = cosineSimilarity(expectedYs, foundYs);
        } else {
            score = entropyScore(expectedYs, foundYs);
        }
        if (score < MIN_SIMILARITY) {
            score = MIN_SIMILARITY;     // cap at minimum similarity to prevent extreme values and log(0) issues
        }
        return score;
    }


    /**
     * Glycan frequency score calculator. Essentially a prior for how likely a given glycan is given delta mass bin.
     * Normalized to the most frequent glycan in the bin.
     * @param candidate   glycan candidate
     * @param deltaMass delta mass bin in question
     * @return frequency (between 0 and 1) of the glycan in the delta mass bin.
     */
    public double glycanFrequencyScore(GlycanCandidate candidate, double deltaMass) {
        // determine the overall likelihood priors of these glycans given the observed delta mass
        int candidateCount = 0;
        int maxFrequency = 0;
        // Decoys are not included in the saved glycans. Use the frequency of the corresponding target
        String glycanHash = Glycan.toGlycanString(candidate.composition);
        LinkedHashMap<String, Integer> emptyMap = new LinkedHashMap<>();
        int massBin = (int) Math.floor(deltaMass);
        LinkedHashMap<String, Integer> glycanCountMap = glycanMassBinMap.getOrDefault(massBin, emptyMap);
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

    /**
     * Compute normalized mass error score for a glycan candidate.
     * Score = probability from normal distribution centered at mean mass error with width 2 * mass error width.
     */
    private double normedMassScore(GlycanCandidateResult candidate, GlycanAssignmentResult result) {
        int isotopeErr = Math.toIntExact(Math.round(result.deltaMass - candidate.mass));
        double massError = result.deltaMass - candidate.mass - (isotopeErr * AAMasses.averagineIsotopeMass);
        candidate.massError = massError;
        return Math.exp(-0.5 * Math.pow((massError - meanMassError) / (2 * massErrorWidth), 2));
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
    private double klMS1score(GlycanCandidate candidate, Spectrum spec, double pepmass) {
        double candidateMass = candidate.mass + pepmass;    // pep mass + glycan mass
        double mz = Spectrum.calcMZ(candidateMass, spec.charge);
        // todo: change to use exact comp (plus compare) (needs IonQuant API to expose that option)
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
     * Calculate MS1 isotope envelope centroid score for glycan candidate. Uses IonQuant API to get the experimental
     * isotope envelope and compares the centroid (weighted average) m/z to the theoretical centroid m/z for
     * the glycan candidate. Requires glycan elemental compositions to be provided.
     * Note: modification elemental compositions are NOT used, but would improve the results. Instead, mod masses
     * are assumed to be composed of averagine, which is not ideal. Peptide sequence and glycan elemental compositions
     * are used to compute the theoretical isotope envelope.
     *
     * Score: log(1 / |theoretical centroid - experimental centroid|)
     * A score of 0 is given if no experimental envelope is found.
     *
     * @param candidate current glycan candidate to score
     * @param spec spectrum
     * @param glycoResult glycan assignment result container with peptide sequence and mods info
     * @return double: MS1 score
     */
    private double ms1CentroidScore(GlycanCandidateResult candidate, Spectrum spec, GlycanAssignmentResult glycoResult) {
        // Calculate theoretical isotope envelope for peptide + glycan + mods
        ElementalComposition glycanComp = candidate.getElementalCompositionOfIon();
        ElementalComposition peptideComp = ElementalComposition.getElementalCompositionOfPeptide(glycoResult.peptide);
        peptideComp.addComposition(new ElementalComposition("H2O"), false, 1);   // add terminal H2O
        glycanComp.addComposition(peptideComp, false, 1);
        double modsMass = 0;
        for (Mod mod : glycoResult.assignedModsList) {
            modsMass += mod.mass;
        }
        double monoMassThy = glycanComp.getMass() + modsMass + candidate.decoyMassShift;
        double[] thyIsoIntensities = isotopeDistributionApi.getTheoIsotopeDistribution(glycanComp.getCompositionMap());
        // clean up theoretical intensities (normalize, remove near-zeroes)
        thyIsoIntensities = normalize(thyIsoIntensities, 0.1, 0);  // todo: make minValue a param

        // Find and normalize experimental isotope envelope (assumes 1st theoretical peak is mono, but might not be true if there are unusual elements in the glycan)
        LinkedHashMap<Integer, Double> isotopeEnvelope = new LinkedHashMap<>();
        double theoMonoMz = Spectrum.calcMZ(monoMassThy, spec.charge);
        Entry monoisotopicPeak = api.quantXIC((float) theoMonoMz, (float) spec.rt, (float) spec.im, spec.charge, spec.cv);
        if (monoisotopicPeak != null) {
            isotopeEnvelope.put(0, (double) monoisotopicPeak.intensity);
        }
        isotopeEnvelope.putAll(searchForIsotopePeaks(spec, monoMassThy, 1, 1));
        isotopeEnvelope.putAll(searchForIsotopePeaks(spec, monoMassThy, -1, -1));
        if (isotopeEnvelope.isEmpty()) {
            return 0;
        }
        double maxIntensity = Collections.max(isotopeEnvelope.values());

        LinkedHashMap<Integer, Double> entropyScores =  new LinkedHashMap<>();
        for (int isoOffset: glycoParams.glycoIsotopes) {
            double[] obsIntensities = new double[thyIsoIntensities.length];
            for (int i = 0; i < thyIsoIntensities.length; i++) {
                obsIntensities[i] = isotopeEnvelope.getOrDefault(i + isoOffset, 0.0);
            }
            obsIntensities = normalize(obsIntensities, 0, maxIntensity);    // normalize observed intensities using the peak intensity of the cluster
            entropyScores.put(isoOffset, entropyScore(thyIsoIntensities, obsIntensities));
        }

        double maxScore = Collections.max(entropyScores.values());
        double deltaScore = 1 + entropyScores.get(0) - maxScore;     // 1 if this is the best score, <1 if not (ent - maxScore always negative)
        candidate.ms1Score = entropyScores.get(0);
        candidate.ms1DeltaScore = deltaScore;
        return entropyScores.get(0);
    }

    /**
     * Check for isotope peaks in the spectrum starting at the given m/z and isotope and continuing in the provided
     * step (direction) until a peak is not found.
     * @return LinkedHashMap of isotope index to intensity for all found peaks
     */
    /**
     * Search for isotope peaks in the envelope. Peaks may initially increase to a maximum and
     * then decrease, but once the sequence has started decreasing any increase terminates the
     * envelope (no second bump is allowed).
     */
    private static LinkedHashMap<Integer, Double> searchForIsotopePeaks(Spectrum spec, double neutralMass, int startIso, int step) {
        LinkedHashMap<Integer, Double> isotopeEnvelope = new LinkedHashMap<>();
        boolean foundPeak = true;
        int isoIndex = startIso;
        int count = 0;
        Double prevPeakIntensity = null;    // intensity of the most recently added peak
        boolean wasDecreasing = false;      // true once the sequence has started decreasing
        while (foundPeak) {
            if (count > MAX_ISOTOPE_PEAKS) {
                break;
            }
            double mz = Spectrum.calcMZ(neutralMass, spec.charge);
            double newMz = (mz + ((isoIndex + step) * ISOTOPE_MASS_DIFF ) / (double) spec.charge);
            double newMz2 = Spectrum.calcMZ(mz + ((isoIndex + step) * ISOTOPE_MASS_DIFF), spec.charge);
            double newMzGood = Spectrum.calcMZ(neutralMass + ((isoIndex + step) * ISOTOPE_MASS_DIFF), spec.charge);

            Entry isotopePeak = api.quantXIC((float) newMz, (float) spec.rt, (float) spec.im, spec.charge, spec.cv);
            if (isotopePeak != null) {
                double intensity = isotopePeak.intensity;
                if (prevPeakIntensity != null) {
                    // Once decreasing has started, any increase is invalid (outside a 5% tolerance to account for noise)
                    if (wasDecreasing && intensity > (prevPeakIntensity * 1.05)) {
                        break;
                    }
                    if (intensity <= prevPeakIntensity) {
                        wasDecreasing = true;
                    }
                }
                prevPeakIntensity = intensity;
                isotopeEnvelope.put(isoIndex, intensity);
                isoIndex += step;
            } else {
                foundPeak = false;
            }
            count++;
        }
        return isotopeEnvelope;
    }

    /**
     * Compute unweighted spectral entropy between the two spectrum vectors, as in
     * www.nature.com/articles/s41592-021-01331-z. Assumes spectrum vectors are of equal
     * length.
     */
    private double entropyScore(double[] theoreticalPks, double[] exptPks) {
        if (Arrays.stream(theoreticalPks).sum() == 0) {
            return 0;
        }
        double[] SabVector = new double[theoreticalPks.length];
        int numFrags = 0;
        for (double j : exptPks) {
            if (j != 0) {
                numFrags += 1;
            }
        }
        double[] normThyPks = normalize(theoreticalPks, 0, 0);
        double[] normExptPks = normalize(exptPks, 0, 0);

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
     * Compute cosine similarity between two spectrum vectors. Returns a score between 0 and 1,
     * where 1 indicates identical spectra and 0 indicates completely dissimilar spectra.
     * Assumes spectrum vectors are of equal length.
     *
     * @param theoreticalPks theoretical peak intensities
     * @param exptPks experimental peak intensities
     * @return cosine similarity score (0 to 1)
     */
    private double cosineSimilarity(double[] theoreticalPks, double[] exptPks) {
        if (Arrays.stream(theoreticalPks).sum() == 0 || Arrays.stream(exptPks).sum() == 0) {
            return 0;
        }
        double dotProduct = 0;
        double normThySq = 0;
        double normExptSq = 0;

        for (int i = 0; i < theoreticalPks.length; i++) {
            dotProduct += theoreticalPks[i] * exptPks[i];
            normThySq += theoreticalPks[i] * theoreticalPks[i];
            normExptSq += exptPks[i] * exptPks[i];
        }

        double magnitude = Math.sqrt(normThySq) * Math.sqrt(normExptSq);
        if (magnitude == 0) {
            return 0;
        }
        return dotProduct / magnitude;
    }

    /**
     * Normalize input vector so that the sum of all entries is 1. If normalizationFactor is provided (nonzero), use it
     * instead of the max of the input array
     */
    private double[] normalize(double[] vector, double minNormalizedValue, double normalizationFactor) {
        double max = 0;
        if (normalizationFactor > 0) {
            max = normalizationFactor;
        } else {
            for (double i : vector) {
                if (i > max) {
                    max = i;
                }
            }
            if (max == 0) {
                return vector;

            }
        }

        ArrayList<Double> normalized = new ArrayList<>();
        for (double v : vector) {
            double normVal = v / max;
            if (normVal >= minNormalizedValue) {
                normalized.add(v / max);
            }
        }
        return normalized.stream().mapToDouble(Double::doubleValue).toArray();
    }


    /**
     * Generate the feature vector for LDA and the summed score for the glycan candidate result using the
     * list of features to use.
     */
    private void generateScores(GlycanCandidateResult candidate) {
        ArrayList<Double> features = new ArrayList<>();
        double summedScore = 0;

        // some scores are always included
        // Y ions
        if (isFirstPass) {
            features.add(candidate.YFragmentScore);
            summedScore += candidate.YFragmentScore;
        } else {
            features.add(candidate.ySpecSim);
            summedScore += candidate.ySpecSim;
        }
        // oxonium ions
        features.add(candidate.OxFragmentScore);
        summedScore += candidate.OxFragmentScore;
        if (!isFirstPass) {
            features.add(candidate.oxSpecSim);
            summedScore += candidate.oxSpecSim;
        }
        // mass error
        features.add(candidate.massErrorScore);
        summedScore += candidate.massErrorScore;

        // some scores are only included if requested
        for (GlycoParams.LDAFeature feature : glycoParams.ldaFeaturesToUse) {
            switch (feature) {
                case kl: // KL score
                    features.add(candidate.klScore);
                    summedScore += candidate.klScore;
                    break;
                case ms1: // MS1 centroid score
                    features.add(candidate.ms1Score);
                    summedScore += candidate.ms1Score;
                    break;
                case ms1delta:
                    features.add(candidate.ms1DeltaScore);
                    summedScore += candidate.ms1DeltaScore;
                    break;
                case glycanfreq:
                    if (!isFirstPass) {     // frequency can only be computed in the second pass
                        features.add(candidate.frequencyPrior);
                        summedScore += candidate.frequencyPrior;
                    }
                    break;
            }
        }
        candidate.summedScore = summedScore;
        if (!glycoParams.glycoLDA || isFirstPass) {
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
