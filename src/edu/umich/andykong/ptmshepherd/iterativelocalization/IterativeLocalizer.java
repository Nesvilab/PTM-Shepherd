package edu.umich.andykong.ptmshepherd.iterativelocalization;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.utils.Peptide;
import edu.umich.andykong.ptmshepherd.utils.StringParsingUtils;

import java.io.File;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

import static edu.umich.andykong.ptmshepherd.PTMShepherd.reNormName;
import static edu.umich.andykong.ptmshepherd.utils.StringParsingUtils.subString;


public class IterativeLocalizer {
    MXMLReader mr;
    int nThreads;
    TreeMap<String, ArrayList<String[]>> datasets;
    HashMap<String, HashMap<String, File>> mzMap;

    FastLocator locate;
    double[][] peaks;
    int zeroBin;

    BinPriorProbabilities[] priorProbs;
    MatchedIonDistribution matchedIonDist;

    double convCriterion;
    int maxEpoch;

    String allowedAAs;
    String decoyAAs;
    String ionTypes;
    float fragTol;

    static String pep;
    static int scanNum;
    static boolean debugFlag;
    boolean printIonDistribution = true; // TODO make this a parameter

    //1. Learn distribution of intensities from unmodified peptides //TODO change this description...
    //   Loop through files
    //2. Loop through bins
    //3.    While bin has not converged
    //4.        Loop through spectra
    //5.            Evaluate probability
    //6.            Cache change in gross probability function
    //7.        Update

    public IterativeLocalizer(double[][] peakBounds, double peakTol, int precursorMassUnits,
                              TreeMap<String, ArrayList<String[]>> ds, HashMap<String, HashMap<String, File>> mzs,
                              int nThreads, String allowedAAs, float fragTol, String ionTypes, double convCriterion,
                              int maxEpoch) {
        initPrecursorPeakBounds(peakBounds, peakTol, precursorMassUnits);
        this.priorProbs = new BinPriorProbabilities[this.peaks[0].length];
        for (int i = 0; i < this.priorProbs.length; i++)
            this.priorProbs[i] = new BinPriorProbabilities();
        this.datasets = ds;
        this.mzMap = mzs;
        this.nThreads = nThreads;
        this.allowedAAs = allowedAAs;
        this.ionTypes = ionTypes;
        this.fragTol = fragTol;
        this.convCriterion = convCriterion;
        this.maxEpoch = maxEpoch;
        this.decoyAAs = "ALIV"; // Well, you can tell by the way I use my walks //TODO make this a parameter if its worth it
    }

    private void initPrecursorPeakBounds(double[][] peakBounds, double peakTol, int precursorMassUnits) {
        this.peaks = peakBounds;
        this.locate = new FastLocator(peakBounds, peakTol, precursorMassUnits);
        this.zeroBin = this.locate.getIndex(0.0);
    }

    public void localize() throws Exception {
        // Step 1: use matched intensities of unmodified peptides to fit nonparametric distribution
        fitMatchedIonDistribution();
        // Step 2: draw the rest of the fucking owl
        calculateLocalizationProbabilities();
        // Step 3: calculate FLRs
        calculateFalseLocalizationRates();
    }

    private void fitMatchedIonDistribution() throws Exception {
        System.out.println("\tFitting distribution to matched zero-bin fragments");
        // Set up distribution
        this.matchedIonDist = new MatchedIonDistribution(0.01f);
        // Set up zero bin boundaries for faster comp
        double zbL = this.peaks[1][this.zeroBin];
        double zbR = this.peaks[2][this.zeroBin];
        // Set up missing spectra error handling
        ArrayList<String> linesWithoutSpectra = new ArrayList<>(); //todo

        // Loop through datasets
        long t1 = System.currentTimeMillis();
        for (String ds : this.datasets.keySet()) {
            ArrayList<String[]> dsData = this.datasets.get(ds);
            // Loop through PSM files
            for (int i = 0; i < dsData.size(); i++) {
                PSMFile psmf = new PSMFile(new File(dsData.get(i)[0]));
                HashMap<String, ArrayList<Integer>> runToLine = psmf.getRunMappings();
                // Loop through runs
                for (String cf : runToLine.keySet()) {
                    // Load current run
                    mr = new MXMLReader(mzMap.get(ds).get(cf), this.nThreads);
                    mr.readFully();
                    // Get matched ion intensities for unmodified peptides
                    for (int j : runToLine.get(cf)) {
                        PSMFile.PSM psm = psmf.getLine(j);
                        float dMass = psm.getDMass();

                        // Limit to unmodified peptides
                        if ((dMass <= zbL) || (dMass >= zbR))
                            continue;
                        String specName = psm.getSpec();
                        String pep = psm.getPep();
                        float[] mods = psm.getModsAsArray();

                        Spectrum spec = mr.getSpectrum(specName);
                        if(spec == null) {
                            linesWithoutSpectra.add(specName); //TODO handle this
                            continue;
                        }

                        float [] matchedInts = spec.getMatchedFrags(pep, mods, this.fragTol, this.ionTypes, 0.0F);
                        this.matchedIonDist.addIons(matchedInts);
                    }
                }
            }
        }
        this.matchedIonDist.calculateCdf();
        if (this.printIonDistribution)
            this.matchedIonDist.printHisto("matched_ion_distribution.tsv");

        long t2 = System.currentTimeMillis();
        System.out.printf("\tDone fitting distribution to matched zero-bin fragments (%d ms processing)\n", t2-t1);
    }

    private void calculateLocalizationProbabilities() throws Exception {
        System.out.println("\tCalculating PSM-level localization probabilities");

        // Set up missing spectra error handling
        ArrayList<String> linesWithoutSpectra = new ArrayList<>(); //todo

        // Faster access to zero bin spectra to be ignored
        double zbL = this.peaks[1][this.zeroBin]; // TODO: set up custom bounds?
        double zbR = this.peaks[2][this.zeroBin]; // TODO: set up custom bounds?

        long t1 = System.currentTimeMillis();
        int epoch = 1;
        int totalBins = this.peaks[0].length - 1;
        int convergedBins = 0;
        boolean finalPass = false;
        boolean complete = false;

        // Check whether there is still data to be processed
        while (!complete) {
            long t2 = System.currentTimeMillis();
            for (String ds : this.datasets.keySet()) {
                ArrayList<String[]> dsData = this.datasets.get(ds);

                // Loop through PSM files
                for (int i = 0; i < dsData.size(); i++) {
                    PSMFile psmf = new PSMFile(new File(dsData.get(i)[0]));
                    HashMap<String, ArrayList<Integer>> runToLine = psmf.getRunMappings();
                    ArrayList<String> specNames = new ArrayList<>();
                    ArrayList<String> strOutputProbs = new ArrayList<>();
                    ArrayList<String> strMaxProbs = new ArrayList<>();
                    ArrayList<String> strEntropies = new ArrayList<>();
                    ArrayList<String> strMaxProbs2 = new ArrayList<>();

                    // Loop through runs
                    for (String cf : runToLine.keySet()) {
                        // Load current run
                        mr = new MXMLReader(mzMap.get(ds).get(cf), this.nThreads);
                        mr.readFully();

                        // Calculate PSM-level localization probabilities
                        for (int j : runToLine.get(cf)) {
                            this.pep = null;
                            this.scanNum = -1;
                            debugFlag = false;

                            PSMFile.PSM psm = psmf.getLine(j);
                            float dMass = psm.getDMass();
                            String pep = psm.getPep();
                            String specName = psm.getSpec();
                            int cBin = this.locate.getIndex(dMass);

                            // Ignore cases that are not in any bin, unless on the last pass and writing results
                            float dMassApex;
                            if (cBin == -1) {
                                if (finalPass)
                                    dMassApex = dMass;
                                else
                                    continue;
                            }
                            else
                                dMassApex = (float) this.peaks[0][cBin]; // TODO which ones to use?

                            // Ignore zero bin, unless on the last pass and writing results TODO: set up custom bounds?
                            if ((zbL <= dMass) && (dMass <= zbR)) {
                                if (finalPass) {
                                    specNames.add(specName);
                                    strOutputProbs.add(""); // Add empty string if zero bin
                                    strMaxProbs.add("");
                                    strEntropies.add("");
                                    strMaxProbs2.add("");
                                    continue;
                                }
                                else
                                    continue;
                            }

                            // If converged skip unless final round
                            if (!finalPass && this.priorProbs[cBin].getIsConverged()) // Safe because left is first
                                continue; // TODO XXX this is breaking

                            Spectrum spec = mr.getSpectrum(specName);
                            if (spec == null) {
                                linesWithoutSpectra.add(specName);
                                continue;
                            }

                            //if (pep.equals("LYGYQAQHYVCMK") && spec.scanNum == 24865) {
                            //    this.debugFlag = true;
                            //}

                            // Calculate site-specific localization probabilities
                            float[] mods = psm.getModsAsArray();
                            boolean[] allowedPoses = parseAllowedPositions(pep, this.allowedAAs, mods);
                            double[] siteProbs = localizePsm(spec, pep, mods, dMassApex, cBin, allowedPoses); // TODO check whether raw, theoretical, or peakapex is better

                            // Update prior probabilities
                            if (!finalPass)
                                this.priorProbs[cBin].update(pep, siteProbs, allowedPoses);
                            else {
                                specNames.add(specName);
                                strOutputProbs.add(probabilitiesToPepString(pep, dMass, siteProbs, allowedPoses));
                                double maxProb = findMaxLocalizationProbability(siteProbs);
                                strMaxProbs2.add(Double.toString(maxProb));
                                String maxProbAA = findMaxLocalizationProbabilitySite(siteProbs, pep);
                                strMaxProbs.add(maxProbToString(maxProb, maxProbAA));
                                double locEntropy = calculateLocalizationEntropy(siteProbs, allowedPoses);
                                strEntropies.add(entropyToString(locEntropy));
                            }

                        }
                    }
                    // If models have converged, update PSM tables
                    if (finalPass) {
                        // Update PSM table with new columns
                        psmf.addColumn(psmf.getColumn("Observed Modifications") + 1,
                                "delta_mass_loc", specNames, strOutputProbs);
                        psmf.addColumn(psmf.getColumn("delta_mass_loc") + 1,
                                "delta_mass_maxloc", specNames, strMaxProbs);
                        psmf.addColumn(psmf.getColumn("delta_mass_maxloc") + 1,
                                "delta_mass_entropy", specNames, strEntropies);
                        psmf.addColumn(psmf.getColumn("delta_mass_entropy") + 1,
                                "delta_mass_maxprob", specNames, strMaxProbs2);
                        psmf.save(true); // Do not overwrite
                        complete = true;
                    }
                }
            }
            long t3 = System.currentTimeMillis();

            // If analysis is complete and results have been written, exit loop
            if (finalPass) {
                complete = true;
                break;
            }

            // Check bins for convergence
            advanceEpochBins(this.convCriterion);
            convergedBins = checkBinConvergences();
            double avgNorm = calcAvgDistToConvergence();

            // Check if analysis is finished
            if (convergedBins == totalBins || epoch == this.maxEpoch)
                finalPass = true;

            // Output status
            if (convergedBins < totalBins && epoch < this.maxEpoch) {
                System.out.printf("\t\tEpoch %d of %d: %d of %d bins have converged, average norm in remaining bins %f (%d" +
                        " ms processing)\n", epoch, this.maxEpoch, convergedBins, totalBins, avgNorm, t3 - t2);
            } else if (convergedBins == totalBins) {
                System.out.printf("\t\tEpoch %d of %d: %d of %d bins have converged, average norm in remaining bins 0.0 (%d" +
                        " ms processing)\n", epoch, this.maxEpoch, convergedBins, totalBins, t3-t2);
                System.out.printf("\tConverged in %d epochs (%d ms processing)\n", epoch, t3-t1);
            } else if (convergedBins < totalBins && epoch == this.maxEpoch) {
                System.out.printf("\t\tEpoch %d of %d: %d of %d bins have converged, average norm in remaining bins %f (%d" +
                        " ms processing)\n", epoch, this.maxEpoch, convergedBins, totalBins, avgNorm, t3-t2);
                System.out.printf("\t%d of %d bins did not converge in %d epochs (%d ms processing)\n",
                        totalBins-convergedBins, totalBins, epoch, t3-t2);
            }

            if (finalPass)
                System.out.println("\tUpdating PSM tables with localization information");

            epoch++;
        }
    }

    //todo mods should be parsed here if we don't want to localized on top of var mods
    //todo test parsing with all different AAs
    private boolean[] parseAllowedPositions(String seq, String allowedAAs, float[] mods) {
        boolean[] allowedPoses = new boolean[seq.length()+2];
        if (allowedAAs.equals("all") || allowedAAs.equals(""))
            Arrays.fill(allowedPoses, true);
        else {
            Arrays.fill(allowedPoses, false);
            for (int i = 0; i < allowedAAs.length(); i++) {
                if (allowedAAs.charAt(i) == 'n') {
                    i++;
                    if (allowedAAs.charAt(i) == '*')
                        allowedPoses[0] = true;
                    else if (allowedAAs.charAt(i) == seq.charAt(0))
                        allowedPoses[1] = true;
                } else if (allowedAAs.charAt(i) == 'c') {
                    i++;
                    if (allowedAAs.charAt(i) == '*')
                        allowedPoses[seq.length()+1] = true;
                    else if (allowedAAs.charAt(i) == seq.charAt(seq.length()-1))
                        allowedPoses[seq.length()+1] = true;
                } else {
                    for (int j = 0; j < seq.length()-1; j++) {
                        if (seq.charAt(j) == allowedAAs.charAt(i))
                            allowedPoses[j+1] = true;
                    }
                }
            }
        }

        /** TODO this will be necessary to debug the allowed positions if we decide not allow loc on top of mods
        for (int i = 0; i < allowedPoses.length; i++) {
            System.out.print((allowedPoses[i] == true ? 1 : 0) + " ");
        }
        System.out.println();
         **/
        return allowedPoses;
    }

    /**
     * Estimates false localization rates in two ways
     * The first assumes that the probabilities computed by the model are well-calibrated and computed the BH
     * adjusted q-values \hat{FLR} (sum_0^i 1-(P(loc_i)) / i for each PSM i.
     * The second does not assume that the probabilities computed by the model are correct, and instead estimates the
     * FLR using decoy amino acids. This automatically removes any modifications on these residues from consideration,
     * which may not be desirable under all circumstances. This is calculated using the unbiased estimator (d+1)/t
     * from Levitsky J Proteome Res. (2017), but adjusted by the ratio of decoy AAs/target AAs.
     *
     * @return
     */
    private void calculateFalseLocalizationRates() throws Exception { //TODO this needs to be modularized so it can be unit tested
        System.out.println("\tEstimating false localization rates");

        long t1 = System.currentTimeMillis();

        int nTargetAAs = 0;
        int nDecoyAAs = 0;
        // Build histos of max localization probs, 4 digit accuracy plus one bin for 1.0000
        int[] targetProbs = new int[1000+1];
        int[] decoyProbs = new int[1000+1];
        int totalPsms = 0;
        // Build histos of localization entropies, 4 digit accuracy plus one bin for 1.0000
        int[] targetEntropies = new int[1000+1]; //4 digit accuracy, plus one bin for 1.0000
        int[] decoyEntropies = new int[1000+1];

        // Loop through datasets
        for (String ds : this.datasets.keySet()) {
            ArrayList<String[]> dsData = this.datasets.get(ds);
            // Loop through PSM files
            for (int i = 0; i < dsData.size(); i++) {

                // Get values we're working with on first pass
                PSMFile psmf = new PSMFile(new File(dsData.get(i)[0]));
                ArrayList<String> peps = psmf.getColumnValues("Peptide");
                ArrayList<String> maxProbs = psmf.getColumnValues("delta_mass_maxloc");
                ArrayList<String> entropies = psmf.getColumnValues("delta_mass_entropy");

                // Add probs to target and decoy histos and calculate nTarget and nDecoyAAs
                for (int j = 0; j < peps.size(); j++) {
                    String pep = peps.get(j);
                    // Skip if PSM is not localized, otherwise count it
                    String maxProb = maxProbs.get(j);
                    if (maxProb.equals(""))
                        continue;
                    else
                        totalPsms++;

                    // Count decoy and target AA counts for this line
                    for (int k = 0; k < pep.length(); k++) {
                        if(isDecoyAA(pep.charAt(k)))
                            nDecoyAAs++;
                        else
                            nTargetAAs++;
                    }

                    // Add max probability and entropy
                    float p = Float.parseFloat(subString(maxProb, "(", ")"));
                    float entropy = Float.parseFloat(entropies.get(j));
                    if (isDecoyAA(maxProb.charAt(0))) {
                        decoyProbs[(int) (p * 1000.0)]++;
                        decoyEntropies[(int) (entropy * 1000.0)]++;
                    } else {
                        targetProbs[(int) (p * 1000.0)]++;
                        targetEntropies[(int) (entropy * 1000.0)]++;
                    }
                }
            }
        }

        // Calculate the FLRs of each type
        // Three different FLRs to compute
        double[] flrProb = new double[1000+1]; // Assumes model probabilities are valid, does not use decoys
        double[] flrProbDecoy = new double[1000+1]; // Uses decoys instead
        double[] flrEntropyDecoy = new double[1000+1]; // Uses decoys and entropy
        double tarDecRatio = (double) nTargetAAs / (double) nDecoyAAs;

        // Compute FLR based on assumption that model is true using BH correction in blocks of 0.0001 probabilities
        double runningFlr = 0.0;
        int runningHits = 0;
        for (int i = flrProb.length-1; i >= 0; i--) {
            // Note: we are assuming decoys are masked, and since they are real AAs rather than reversed sequences,
            // they are included in the number of hits
            int cHits = targetProbs[i] + decoyProbs[i];
            double proportionOfPSMs = safeDivide(cHits, totalPsms);
            double cLocP = (double) i / 1000; // Current probability
            double binPVal = 1.0 - cLocP;
            double binFlr = binPVal * proportionOfPSMs;

            // Add new IDs, but weight it by their contribution to total number of IDs above them in the ranking
            runningFlr += binFlr * safeDivide(cHits, cHits + runningHits);

            flrProb[i] = runningFlr;
            runningHits += cHits;
        }

        // Calculate FLR based on decoy AAs and maxProb, don't forget to multipy by tar/dec ratio
        int cDecoys = 1; // unbiased estimator has d+1 in numerator
        int cTargets = 0;
        for (int i = flrProbDecoy.length-1; i >= 0; i--) {
            cDecoys += decoyProbs[i];
            cTargets += targetProbs[i];
            double qVal = safeDivide(cDecoys, cTargets) * tarDecRatio;
            flrProbDecoy[i] = qVal;
        }

        // Calculate FLR based on decoy AAs and entropy, don't forget to multipy by tar/dec ratio
        cDecoys = 1; // unbiased estimator has d+1 in numerator
        cTargets = 0;
        for (int i = 0; i < flrEntropyDecoy.length; i++) {
            cDecoys += decoyEntropies[i];
            cTargets += targetEntropies[i];
            double qVal = safeDivide(cDecoys, cTargets) * tarDecRatio;
            flrEntropyDecoy[i] = qVal;
        }

        // Reestimate using the minimum q val of items > i to make q-values monotonic
        double[] qProbModel = new double[1000+1]; // Assumes model probabilities are valid, does not use decoys
        double[] qProbDecoyModel = new double[1000+1]; // Uses decoys instead TODO returning all 0
        double[] qEntropyDecoyModel = new double[1000+1]; // Uses decoys and entropy TODO returning all 0.0001
        double min = 1000.0;
        for (int i = 0; i < flrProb.length; i++) {
            if (flrProb[i] < min)
                min = flrProb[i];
            qProbModel[i] = min;
        }

        min = 1000.0;
        for (int i = 0; i < flrProbDecoy.length; i++) {
            if (flrProbDecoy[i] < min)
                min = flrProbDecoy[i];
            qProbDecoyModel[i] = min;
        }

        min = 1000.0;
        for (int i = flrEntropyDecoy.length - 1; i >= 0; i--) {
            if (flrEntropyDecoy[i] < min)
                min = flrEntropyDecoy[i];
            qEntropyDecoyModel[i] = min;
        }

        /**
        // Print to test
        for (int i = flrProb.length-1; i >= 0; i--) {
            double cProb = (double) i / 1000.0;
            System.out.printf("%.4f\t%.4f\t%d\t%d%n", cProb, qProbModel[i], targetProbs[i], decoyProbs[i]);
        }
        System.out.println("*****");

        for (int i = flrProb.length-1; i >= 0; i--) {
            double cProb = (double) i / 1000.0;
            System.out.printf("%.4f\t%.4f\t%d\t%d%n", cProb, qProbDecoyModel[i], targetProbs[i], decoyProbs[i]);
        }

        System.out.println("*****");
        // Print to test
        for (int i = flrProb.length-1; i >= 0; i--) {
            double cProb = (double) i / 1000;
            System.out.printf("%.4f\t%.4f\t%d\t%d%n", cProb, flrEntropyDecoy[i], targetEntropies[i], decoyEntropies[i]);
        }
        System.out.println("*****");
        for (int i = flrProb.length-1; i >= 0; i--) {
            double cProb = (double) i / 1000;
            System.out.printf("%.4f\t%.4f\t%d\t%d%n", cProb, qEntropyDecoyModel[i], targetEntropies[i], decoyEntropies[i]);
        }
         **/

        // Loop through datasets
        for (String ds : this.datasets.keySet()) {
            ArrayList<String[]> dsData = this.datasets.get(ds);
            // Loop through PSM files
            for (int i = 0; i < dsData.size(); i++) {

                // Get values to map to q-vals
                PSMFile psmf = new PSMFile(new File(dsData.get(i)[0]));
                ArrayList<String> specNames = psmf.getColumnValues("Spectrum");
                ArrayList<String> maxProbs = psmf.getColumnValues("delta_mass_maxloc");
                ArrayList<String> entropies = psmf.getColumnValues("delta_mass_entropy");

                // Assign each q-Val
                ArrayList<String> probModelVals = new ArrayList<>(specNames.size());
                ArrayList<String> probDecoyModelVals = new ArrayList<>(specNames.size());
                ArrayList<String> entropyDecoyModelVals = new ArrayList<>(specNames.size());
                for (int j = 0; j < specNames.size(); j++) {

                    // Check if peptide is unmod and deal with it if it is
                    boolean unmodFlag = maxProbs.get(i).equals("") ? true : false;
                    if (unmodFlag) {
                        // prob model
                        probModelVals.add("");
                        // prob model with decoys
                        probDecoyModelVals.add("");
                        // entropy model
                        entropyDecoyModelVals.add("");
                    } else {
                        // prob model
                        double maxProb = Double.parseDouble(subString(maxProbs.get(i), "(", ")"));
                        probModelVals.add(new DecimalFormat("0.0000").format(qProbModel[(int) (maxProb * 1000)]));
                        // prob model with decoys
                        probDecoyModelVals.add(new DecimalFormat("0.0000").format(qProbDecoyModel[(int) (maxProb * 1000)]));
                        // entropy model
                        double entropy = Double.parseDouble(entropies.get(i));
                        entropyDecoyModelVals.add(new DecimalFormat("0.0000").format(qEntropyDecoyModel[(int) (entropy * 1000)]));
                    }
                }

                // Send to PSM file
                psmf.addColumn(psmf.getColumn("delta_mass_maxloc") + 1, "delta_mass_BH_loc_q",
                        specNames, probModelVals);
                /** //TODO figure out what's going on with these before implementing them, assuming they're even worth doing
                psmf.addColumn(psmf.getColumn("delta_mass_BH_loc_q") + 1, "delta_mass_prob_decoyAA_q",
                        specNames, probDecoyModelVals);
                psmf.addColumn(psmf.getColumn("delta_mass_entropy") + 1, "delta_mass_entropy_decoyAA_q",
                        specNames, entropyDecoyModelVals);
                psmf.save(true);
                 **/
            }
        }

        long t2 = System.currentTimeMillis();
        System.out.printf("\tDone estimating false localization rates (%d ms processing)\n", t2-t1);

    }

    boolean isDecoyAA(char aa) {
        for (int i = 0; i < this.decoyAAs.length(); i++) {
            if (aa == this.decoyAAs.charAt(i))
                return true;
        }
        return false;
    }

    /**
     * This function computes the P(Pep_{ij}|Spec_i) =
     * P(Pep_{ij}*P(Spec_i|Pep_{ij}) / Sum_{k=0}^{{L_i}+1} P(Pep_{ik})*P(Spec_i|Pep_{ik})
     *
     * P(Pep_{ij}|Spec_i)                                       ->  Posterior probability
     * P(Pep_{ij}                                               ->  Prior probability
     * P(Spec_i|Pep_{ij})                                       ->  Likelihood
     * Sum_{k=0}^{{L_i}+1} P(Pep_{ik})*P(Spec_i|Pep_{ik})       ->  Marginal probability
     *
     * @param spec          Spectrum class opject containing pre-process mass spectrum
     * @param pep           pep sequence
     * @param mods          array containing masses to be added on to pep sequence at mods[i] position
     * @param dMass         delta mass of PSM
     * @param cBin          current MS1 mass shift bin index
     * @param allowedPoses  array of allowed positions based on peptide sequence localization restrictions TODO add mods
     * @return double[] of localization probabilities
     */
    private double[] localizePsm (Spectrum spec, String pep, float[] mods, float dMass, int cBin, boolean[] allowedPoses) {
        double[] sitePriorProbs;
        double[] siteLikelihoods = new double[pep.length()+2];
        double marginalProb = 0.0;
        double[] sitePosteriorProbs = new double[pep.length()+2];

        // Compute prior probability or use uniform prior if PSM does not belong to mass shift bin
        if (cBin == -1)
            sitePriorProbs = BinPriorProbabilities.computeUniformPriorProbs(pep, allowedPoses); // TODO unclustered mass shifts should take the conditional prob as their final prob, not computed with a uniform prior
        else
            sitePriorProbs = this.priorProbs[cBin].computePriorProbs(pep, allowedPoses);

        // Iterate through sites to compute likelihood for each site P(Spec_i|Pep_{ij})
        // There are no ions that can differentiate termini and terminal AAs, so the likelihood for each terminus
        // is equal to the proximal AA
        // System.out.println(spec.toString());
        /**
        for(int i = 0; i < pep.length(); i++) {
            if (allowedPoses[i+1] == false) {
                siteLikelihoods[i+1] = 0.0;
                continue;
            }
            mods[i] += dMass;
            ArrayList pepFrags = Peptide.calculatePeptideFragments(pep, mods, this.ionTypes, 1);
            if (debugFlag) {
                System.out.println(pep);
                System.out.println("Position " + (i + 1));
                System.out.println(spec.toString());
            }
            float[] peakMzs = spec.getPeakMZ();
            float[] peakInts = spec.getPeakInt();
            siteLikelihoods[i+1] = computeLikelihood(pepFrags, peakMzs, peakInts);
            mods[i] -= dMass;
        }
         **/
        //NEW FUNCTION STARTS HERE
        // First calculate the set of shifted and unshifted ions
        ArrayList<Float> pepFrags = Peptide.calculatePeptideFragments(pep, mods, this.ionTypes, 1);
        ArrayList<Float> shiftedPepFrags = new ArrayList<Float>(pepFrags.size());
        for (Float frag : pepFrags)
            shiftedPepFrags.add(frag + dMass);
        pepFrags.addAll(shiftedPepFrags);
        if (debugFlag)
            System.out.println(pepFrags.stream().map(Object::toString)
                .collect(Collectors.joining(", ")));

        // Filter peakMzs and peakInts to only those that match at least one ion
        float[] peakMzs = spec.getPeakMZ();
        float[] peakInts = spec.getPeakInt();
        float[] matchedIons = findMatchedIons(pepFrags, peakMzs, peakInts); // Returns -1 if unmatched, intensity otherwise
        int matchedCount = 0;
        for (int i = 0; i < matchedIons.length; i++) {
            if (matchedIons[i] > 0.0)
                matchedCount++;
        }
        float[] reducedMzs = new float[matchedCount];
        float[] reducedInts = new float[matchedCount];
        int j = 0;
        for (int i = 0; i < matchedIons.length; i++) {
            if (matchedIons[i] > 0.0) {
                reducedMzs[j] = peakMzs[i];
                reducedInts[j] = peakInts[i];
                j++;
            }
        }

        // Iterate through sites to compute likelihood for each site P(Spec_i|Pep_{ij})
        // There are no ions that can differentiate termini and terminal AAs, so the likelihood for each terminus
        // is equal to the proximal AA
        for(int i = 0; i < pep.length(); i++) {
            if (allowedPoses[i+1] == false) {
                siteLikelihoods[i+1] = 0.0;
                continue;
            }
            mods[i] += dMass;
            ArrayList<Float> sitePepFrags = Peptide.calculatePeptideFragments(pep, mods, this.ionTypes, 1);
            if (debugFlag) {
                System.out.println(pep);
                System.out.println("Position " + (i + 1));
                System.out.println(spec.toString());
            }
            siteLikelihoods[i+1] = computeLikelihood(sitePepFrags, reducedMzs, reducedInts);
            mods[i] -= dMass;
        }

        // Propagate terminal AA likelihoods to each terminus //TODO this isn't going to handle cases when termini are allowed but the first residue isnt
        if (allowedPoses[0])
            siteLikelihoods[0] = siteLikelihoods[1]; // N-term
        if (allowedPoses[allowedPoses.length-1])
            siteLikelihoods[siteLikelihoods.length-1] = siteLikelihoods[siteLikelihoods.length-2]; // C-term

        // Compute marginal probability for peptide Sum_{k=0}^{{L_i}+1} P(Pep_{ik})*P(Spec_i|Pep_{ik})
        // Sum (site prior * site likelihoods)
        for (int i = 0; i < siteLikelihoods.length; i++)
            marginalProb += (sitePriorProbs[i] * siteLikelihoods[i]);

        // Compute the site-level posterior probabilities
        for(int i = 0; i < sitePosteriorProbs.length; i++)
            sitePosteriorProbs[i] = (sitePriorProbs[i] * siteLikelihoods[i]) / marginalProb;

        // TODO remove once unit tests are in place
        if (debugFlag) {
            System.out.print("LH\t");
            for (int i = 0; i < siteLikelihoods.length; i++) {
                System.out.print(siteLikelihoods[i] + "\t");
            }
            System.out.println();
            System.out.print("Prior\t");
            for (int i = 0; i < sitePriorProbs.length; i++) {
                System.out.print(sitePriorProbs[i] + "\t");
            }
            System.out.println();
            System.out.print("Post\t");
            for (int i = 0; i < sitePosteriorProbs.length; i++) {
                System.out.print(sitePosteriorProbs[i] + "\t");
            }
            System.out.println("\n");
        }

        return sitePosteriorProbs;
    }

    /**
     * Computes the likelihood P(Spec_i|Pep_{ij}) of a particular localization site. Matched and unmatched ion
     * intensities are retrieved from the spectrum. Unmatched ions have their intensities multiplied by -1. Matched ions
     * are mapped to the MatchedIonDistribution with negative (unmatched) intensities ignored, then the array is flipped
     * and the unmatched ions in the spectrum are mapped to the MatchedIonDistribution.
     *
     * P(Spec_i|Pep_{ij}) is given by prod(Peak_k|Pep_{ij}),
     * where P(Peak_k|Pep_{ij}) = CDF(intensity_i) if matched and 1-CDF(intensity_i) if unmatched.
     * @param pepFrags  theoretical peptide fragments
     * @param peakMzs   spectrum peak M/Zs
     * @param peakInts  spectrum peak intensities
     * @return
     */
    private double computeLikelihood(ArrayList<Float> pepFrags, float[] peakMzs, float[] peakInts) {
        // Find matched ion intensities. Matched ions will have positive intensities, unmatched ions will have negative
        float[] matchedIonIntensities = findMatchedIons(pepFrags, peakMzs, peakInts);
        // Map matched ion intensities to MatchedIonDistribution, negative intensities will be returned as -1
        double[] matchedIonProbabilities = this.matchedIonDist.calcIonProbabilities(matchedIonIntensities);

        if (debugFlag) {
            System.out.println("*"+Arrays.toString(matchedIonIntensities));
            System.out.println("**"+Arrays.toString(matchedIonProbabilities));
        }

        // Map unmatched ion intensities to MatchedIonDistribution, negative values will be returned as -1
        for (int i = 0; i < matchedIonIntensities.length; i++) // Flip sign to compute unmatched ions
            matchedIonIntensities[i] *= -1;

        double[] unmatchedIonProbabilities = this.matchedIonDist.calcIonProbabilities(matchedIonIntensities);

        if (debugFlag) {
            System.out.println("***"+Arrays.toString(matchedIonIntensities));
            System.out.println("****"+Arrays.toString(unmatchedIonProbabilities));
        }

        //Compute likelihood, product of matched ions and 1-unmatched ions
        double likelihood = 1.0;
        for (int i = 0; i < matchedIonProbabilities.length; i++) {
            if (matchedIonProbabilities[i] >= 0.0) {
                likelihood *= matchedIonProbabilities[i];
            }
        }

        if (debugFlag)
            System.out.println( "*****"+likelihood);

        double unlikelihood = 1.0;
        for (int i = 0; i < unmatchedIonProbabilities.length; i++) {
            if (unmatchedIonProbabilities[i] >= 0.0) {
                unlikelihood *= (1.0 - unmatchedIonProbabilities[i]);
            }
        }

        if (debugFlag) {
            System.out.println("******" + unlikelihood);
        }

        likelihood *= unlikelihood;

        if (debugFlag) {
            System.out.println("*******" + likelihood);
        }

        return likelihood;
    }

    /**
     * To speed up search, sort the pep ions and the spec ions in parallel. Until the next pep ion is reached, multiply
     * by -1 (to signify unmatched intensity) and add to matched ion intensities. Once an ion is matched, find the max
     * intensity matching ion, add to matched ion intensities, and multiply others by -1 to signify unmatched.
     *
     * @param pepFrags  theoretical peptide fragments
     * @param peakMzs   spectrum peak M/Zs
     * @param peakInts  spectrum peak intensities
     * @return float[] of size spectrum_ions, matched ions show intensity, unmatched ions show -1 * intensity
     */
    private float[] findMatchedIons(ArrayList<Float> pepFrags, float[] peakMzs, float[] peakInts) {
        float[] matchedIonIntensities = new float[peakMzs.length];
        Arrays.sort(matchedIonIntensities);
        Collections.sort(pepFrags);

        if (debugFlag)
            System.out.println(pepFrags.stream().map(Object::toString)
                    .collect(Collectors.joining(", ")));

        int specIdx = 0;
        int maxSpecIdx = peakMzs.length;

        // Reduce the pep_frags * spec_size array search by sorting arrays and searching in parallel
        for (Float frag : pepFrags) {
            double minVal = frag - (frag * this.fragTol / 1000000.0f);
            double maxVal = frag + (frag * this.fragTol / 1000000.0f);
            // Until we get to a matched ion, multiply intensities by -1 to signify that they are unmatched
            while (specIdx < maxSpecIdx) {
                if (peakMzs[specIdx] < minVal) {
                    matchedIonIntensities[specIdx] = -1 * peakInts[specIdx];
                    specIdx++;
                } else
                    break;
            }
            // Find the maximum intensity of ions between fragment's minVal and maxVal
            int wIndex = specIdx; // This pointer holds the index within the window to search for max intensity fragment
            int startIdx = specIdx; // This point is the index at which the window begins
            float maxInt = 0.0f;
            int matchedIons = 0;
            while (wIndex < maxSpecIdx && peakMzs[wIndex] <= maxVal) { // This is only safe because left side is evald first
                if (peakInts[wIndex] > maxInt) {
                    maxInt = peakInts[wIndex];
                    matchedIons++;
                }
                wIndex++;
            }
            wIndex--;
            if (matchedIons == 1) { // If only matched one ion in minVal->maxVal window
                if (debugFlag)
                    System.out.println(peakMzs[wIndex]);
                matchedIonIntensities[wIndex] = maxInt;
                specIdx = wIndex + 1;
            } else { // If more than one ion in minVal->maxVal window, max is matched and rest are unmatched
                for (int i = startIdx; i <= wIndex; i++) {
                    if (peakInts[i] == maxInt) {
                        if (debugFlag)
                            System.out.println(peakMzs[wIndex]);
                        matchedIonIntensities[i] = maxInt;
                        specIdx = i + 1;
                        break;
                    }
                }
            }
        }

        if (debugFlag)
            System.out.println(Arrays.toString(matchedIonIntensities));

        return matchedIonIntensities;
    }

    /**
     * Calculate the Shannon entropy of posterior localization probabilities for a peptide.
     * H_B = -sum_{i=0}^{i=N} p_i*log(p_i)
     * The formula, taken from PTMProphet Shteynberg et al. J Proteome Res. (2019) is
     * H_t^{norm} = - (sum_{i=1}^S p_i*log_{s/m}p_i) / m
     *      t is the modification type
     *      s is the number of potential sites of modification type t
     *      m is the number of modifications of type t
     *      p_i is the site probability
     *
     * @param locProbs      peptide length+2 probabilities
     * @param allowedPoses   peptide length+2 length bool array of allowed localizations
     * @return localization entropy
     */
    private double calculateLocalizationEntropy(double[] locProbs, boolean[] allowedPoses) {
        double locEntropy = 0.0;
        int m = 1;

        // Calculate potential sites
        int s = 0;
        for (int i  = 0; i < allowedPoses.length; i++) {
            if (allowedPoses[i] == true) {
                s++;
            }
        }
        double base = (double) s / (double) m;

        double max = 0;
        for (int i  = 0; i < locProbs.length; i++) {
            if (locProbs[i] > max)
                max = locProbs[i];
        }

        // Calculate entropy
        for (int i = 0; i < locProbs.length; i++)
            locEntropy += locProbs[i] * (Math.log(locProbs[i]) / Math.log(base));
        locEntropy *= -1;

        return locEntropy;
    }

    private double findMaxLocalizationProbability(double[] locProbs) {
        double max = 0.0;
        for (int i  = 0; i < locProbs.length; i++) {
            if (locProbs[i] > max)
                max = locProbs[i];
        }
        return max;
    }

    private String findMaxLocalizationProbabilitySite(double[] locProbs, String pep) {
        double max = 0.0;
        int maxI = -1;
        for (int i  = 0; i < locProbs.length; i++) {
            if (locProbs[i] > max) {
                max = locProbs[i];
                maxI = i;
            }
        }

        String site;
        if (maxI == 0)
            site = "N-term";
        else if (maxI == locProbs.length-1)
            site = "C-term";
        else
            site = pep.charAt(maxI-1) + Integer.toString(maxI);

        return site;
    }

    private void advanceEpochBins(double convCriterion) {
        for (int i = 0; i < this.peaks[0].length; i++) {
            if (i == this.zeroBin)
                continue;
            this.priorProbs[i].calcConvergence(convCriterion);
        }
    }

    private int checkBinConvergences () {
        int nConverged = 0;
        for (int i = 0; i < this.peaks[0].length; i++) {
            if (i == this.zeroBin)
                continue;
            if (this.priorProbs[i].getIsConverged() == true)
                nConverged++;
        }
        return nConverged;
    }

    private double calcAvgDistToConvergence () {
        int nUnconverged = 0;
        double avgNorm = 0.0;
        for (int i = 0; i < this.peaks[0].length; i++) {
            if (i == this.zeroBin)
                continue;
            if (this.priorProbs[i].getIsConverged() == false) {
                nUnconverged++;
                avgNorm += this.priorProbs[i].getPriorVectorNorm();
            }
        }
        return avgNorm / nUnconverged;
    }

    private String probabilitiesToPepString(String pep, float dMass, double[] probs, boolean [] allowedPoses) {
        StringBuffer sb = new StringBuffer();

        sb.append(new DecimalFormat("0.0000").format(dMass));
        sb.append("@");
        if (allowedPoses[0] == true)
            sb.append("("+new DecimalFormat("0.0000").format(probs[0])+")");
        for (int i = 0; i < pep.length(); i++) {
            sb.append(pep.charAt(i));
            if (allowedPoses[i+1] == true)
                sb.append("("+new DecimalFormat("0.0000").format(probs[i+1])+")");
        }
        if (allowedPoses[allowedPoses.length-1] == true)
            sb.append("("+new DecimalFormat("0.0000").format(probs[0])+")");

        return sb.toString();
    }

    private String entropyToString(double entropy) {
        return new DecimalFormat("0.0000").format(entropy);
    }

    private String maxProbToString(double maxProb, String maxProbSite) {
        String strProb = new DecimalFormat("0.0000").format(maxProb);
        return maxProbSite + "(" + strProb + ")";
    }

    private double safeDivide(int x, int y) {
        if (y == 0)
            return 0;
        else
            return (double) x / (double) y;
    }

}