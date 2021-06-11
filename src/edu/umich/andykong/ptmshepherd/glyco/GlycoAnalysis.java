package edu.umich.andykong.ptmshepherd.glyco;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;
import edu.umich.andykong.ptmshepherd.specsimilarity.SimRTProfile;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

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
    public static final int NUM_ADDED_GLYCO_PSM_COLUMNS = 5;
    boolean normYions;
    double defaultMassErrorAbsScore;

    public GlycoAnalysis(String dsName, ArrayList<GlycanCandidate> glycoDatabase, ProbabilityTables inputProbabilityTable, boolean normYs, double absMassErrorDefault) {
        this.dsName = dsName;
        this.glycoFile = new File(PTMShepherd.normFName(dsName + ".rawglyco"));
        this.glycanDatabase = glycoDatabase;

        // init with default values, can be changed by params
        this.probabilityTable = inputProbabilityTable;
        this.normYions = normYs;
        this.defaultMassErrorAbsScore = absMassErrorDefault;
    }

    public void glycoPSMs(PSMFile pf, HashMap<String, File> mzMappings) throws Exception {
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
        StringBuffer headbuff = new StringBuffer(String.format("%s\t%s\t%s\t%s\t%s", "Spectrum", "Peptide", "Mods", "Pep Mass", "Mass Shift"));
        headbuff.append("\tBest Glycan\tLog Delta Score\t2nd Best Glycan\tAbs Score\tGlycan q-value");
        for (int i = 0; i < capYShifts.length; i++)
            headbuff.append(String.format("\tY_%.4f_intensity", capYShifts[i]));
        for (int i = 0; i < oxoniumIons.length; i++)
            headbuff.append(String.format("\tox_%.4f_intensity", oxoniumIons[i]));
        for (int i = 0; i < remainderMasses.length; i++)
            headbuff.append(String.format("\tdeltascore_%.4f\tlocalization_%.4f", remainderMasses[i], remainderMasses[i]));
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
            for (int i = 0; i < clines.size(); i++) {//for relevant line in curr spec file
                String newline = processLine(pf.data.get(clines.get(i)));
                this.totalLines++;
                if (newline.equals("ERROR"))
                    continue;
                else
                    out.println(newline);
            }
            long t3 = System.currentTimeMillis();
            PTMShepherd.print(String.format("\t%s - %d (%d ms, %d ms)", cf, clines.size(), t2 - t1, t3 - t2));
        }
        out.close();

        if (!linesWithoutSpectra.isEmpty()) {
            System.out.printf("Could not find %d/%d (%.1f%%) spectra.\n", linesWithoutSpectra.size(), this.totalLines,
                    100.0*((double)linesWithoutSpectra.size()/this.totalLines));
            int previewSize = Math.min(linesWithoutSpectra.size(), 5);
            System.out.printf("Showing first %d of %d spectra IDs that could not be found: \n\t%s\n", previewSize, linesWithoutSpectra.size(),
                    String.join("\n\t", linesWithoutSpectra.subList(0, previewSize)));
        }
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
                case "Abs Score":
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
            System.out.printf("Warning: rawglyco file headers not found! FDR calculation may fail for file %s\n", glycoFile);
        }

        // read file, accumulating scores
        HashMap<String, Double> scoreMap = new HashMap<>();
        int targets = 0;
        int decoys = 0;
        int numLines = 0;
        while ((cgline = in.readLine()) != null) {
            numLines++;
            if (cgline.equals("COMPLETE")) {
                break;
            }
            String[] splits = cgline.split("\t", -1);
            String spectrumID = splits[gSpecCol];
            glyLines.put(spectrumID, splits);     // save full line for later editing/writing
            // only consider columns with actual glycan info
            if (!splits[bestGlycanCol].matches("") && !splits[bestGlycanCol].toLowerCase(Locale.ROOT).matches("no matches")) {
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
            System.out.printf("Not enough decoys to compute FDR at %.1f pct, initial pct %.2f\n", desiredRatio * 100, targetDecoyRatio * 100);
            if (desiredRatio < targetDecoyRatio * 0.1) {
                System.out.print("Not enough decoys to compute FDR at 0.1 * initial ratio. Check data and parameters. No FDR calculation performed!\n");
                return;
            } else {
                // only missed by a little, try reducing desired FDR to accomodate
                desiredRatio = targetDecoyRatio - (targetDecoyRatio * 0.1);
                System.out.printf("FDR reduced to %.1f pct due to limited decoys\n", desiredRatio * 100);
            }
        }

        // find threshold at which target/decoy ratio hits desired value and update rawglyco lines
        boolean foundThreshold = false;
        for (Map.Entry<String, Double> scoreEntry : sortedScoreMap.entrySet()) {
            String[] rawGlycoLine = glyLines.get(scoreEntry.getKey());
            if (!foundThreshold) {
                // still below the threshold: continue checking decoys/targets and updating counts
                if (rawGlycoLine[bestGlycanCol].toLowerCase(Locale.ROOT).contains("decoy")) {
                    decoys--;
                } else {
                    targets--;
                }
                targetDecoyRatio = decoys / (double) targets;
                if (targetDecoyRatio <= desiredRatio) {
                    // stop here, found cutoff
                    foundThreshold = true;
                    System.out.printf("Converged to %.1f pct FDR with %d targets and %d decoys\n", targetDecoyRatio * 100, targets, decoys);
                }
                rawGlycoLine[qValCol] = String.format("%s", targetDecoyRatio);
                rawGlycoLine[bestGlycanCol] = "FailFDR_" + rawGlycoLine[bestGlycanCol];
            } else {
                // passed FDR threshold - all entries above this point pass. Continue updating q-value
                if (rawGlycoLine[bestGlycanCol].toLowerCase(Locale.ROOT).contains("decoy")) {
                    decoys--;
                } else {
                    targets--;
                }
                targetDecoyRatio = decoys / (double) targets;
                rawGlycoLine[qValCol] = String.format("%s", targetDecoyRatio);
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
        for (int i = 0; i < clines.size(); i++) {//for relevant line in curr spec file
            String line = psmFile.data.get(clines.get(i));
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


    public StringBuffer processLine(String line) {
    public String processLine(String line) {
        StringBuffer sb = new StringBuffer();
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

        sb.append(assignGlycanToPSM(spec, pepMass, dmass, glycanDatabase, massErrorWidth, meanMassError));

        //System.out.println("got spec");
        double[] capYIonIntensities;
        double[] oxoniumIonIntensities;
        capYIonIntensities = findCapitalYIonMasses(spec, pepMass);
        oxoniumIonIntensities = findOxoniumIonMasses(spec, pepMass);

        for (int i = 0; i < capYIonIntensities.length; i++)
            sb.append(String.format("\t%.2f", capYIonIntensities[i]));
        for (int i = 0; i < oxoniumIonIntensities.length; i++)
            sb.append(String.format("\t%.2f", oxoniumIonIntensities[i]));
        float[] deltaScores = new float[remainderMasses.length];
        boolean[][] isMaxScores = localizeRemainderFragments(spec, sp[pepCol], smods, deltaScores);

        for (int i = 0; i < remainderMasses.length; i++) {
            sb.append(String.format("\t%.1f", deltaScores[i]));
            StringBuffer locSb = new StringBuffer("\t");
            for (int j = 0; j < seq.length(); j++) {
                if (isMaxScores[i][j] == true) {
                    locSb.append(String.format("%d%c", j + 1, seq.charAt(j))); //position (1 indexed), character
                }
            }
            sb.append(locSb.toString());
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
        // skip non-delta mass PSMs
        if (deltaMass < 3.5 && deltaMass > -1.5) {
            return "\t\t\t\t\t";
        }

        // Determine possible glycan candidates from mass
        int[] isotopesToSearch = {-1, 0, 1, 2, 3};
        double ms1TolPPM = 50;  // todo: connect to existing param?
        ArrayList<GlycanCandidate> searchCandidates = getMatchingGlycansByMass(deltaMass, glycanDatabase, isotopesToSearch, ms1TolPPM);

        // Get Y/oxo ions possible for these candidates (and decoy ions)
        GlycanFragment[] possibleYIons = initializeYFragments(searchCandidates);
        GlycanFragment[] possibleOxoniums = initializeOxoniumFragments(searchCandidates);
        double[] yMasses = new double[possibleYIons.length];
        for (int i = 0; i < possibleYIons.length; i++) {
            yMasses[i] = possibleYIons[i].neutralMass + pepMass;
        }
        double[] oxoMasses = new double[possibleOxoniums.length];
        for (int i = 0; i < possibleOxoniums.length; i++) {
            oxoMasses[i] = possibleOxoniums[i].neutralMass;
        }

        // Search Y and oxonium ions in spectrum
        float ppmTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
        for (int i = 0; i < possibleYIons.length; i++) {
            possibleYIons[i].foundIntensity = spec.findIonNeutral(yMasses[i], ppmTol);  // sum of charge state intensities if >1 found
        }
        for (int i = 0; i < possibleOxoniums.length; i++) {
            // only search 1+ oxonium ions, not all possible charge states
            possibleOxoniums[i].foundIntensity = spec.findIon(oxoMasses[i] + AAMasses.protMass, ppmTol);
        }

        // score candidates and save results
        int bestCandidateIndex = 0;
        int nextBestCandidateIndex = 1;
        double[] scoresVsBestCandidate = new double[searchCandidates.size()];

        for (int i = 1; i < searchCandidates.size(); i++) {
            if (i == bestCandidateIndex) {
                continue;
            }
            double comparisonScore = pairwiseCompareGlycans(searchCandidates.get(bestCandidateIndex), searchCandidates.get(i), possibleYIons, possibleOxoniums, pepMass, deltaMass, massErrorWidth, meanMassError);
            if (comparisonScore > 0) {
                // best candidate obtained better score and remains unchanged. Put comparison score at position i to indicate the score of this candidate relative to current best candidate
                scoresVsBestCandidate[i] = comparisonScore;
            } else {
                // new best candidate - reset best candidate position and update scores at all other positions
                for (int j = 0; j < i; j++) {
                    scoresVsBestCandidate[j] -= comparisonScore;    // subtract score vs previous best candidate from existing scores to update them
                }
                nextBestCandidateIndex = bestCandidateIndex;
                bestCandidateIndex = i;
                scoresVsBestCandidate[i] = 0;
            }
        }
        // it's possible for the 2nd best score to be worse than those for candidates searched AFTER the last change to best candidate. Compare those against 2nd best score to get correct result
        for (int i = bestCandidateIndex + 1; i < searchCandidates.size(); i++) {
            double comparisonScore2 = pairwiseCompareGlycans(searchCandidates.get(nextBestCandidateIndex), searchCandidates.get(i), possibleYIons, possibleOxoniums, pepMass, deltaMass, massErrorWidth, meanMassError);
            if (comparisonScore2 < 0) {
                // currently listed 2nd best glycan had worse score than this one: change index to this glycan
                nextBestCandidateIndex = i;
                // don't need to recalculate comparison score because we still have it relative to the best glycan from before
            }
        }

        // compute absolute score for best glycan
        double absoluteScore = 0;
        if (searchCandidates.size() > 0) {
            absoluteScore = computeAbsoluteScore(searchCandidates.get(bestCandidateIndex), possibleYIons, possibleOxoniums, deltaMass, massErrorWidth, meanMassError);
        }
        // output - best glycan, scores, etc back to PSM table
        String output;
        if (searchCandidates.size() == 0) {
            output = "\tNo Matches\t\t\t\t";
        } else if (searchCandidates.size() == 1) {
            output = String.format("\t%s\t\t\t%.1f\t", searchCandidates.get(bestCandidateIndex).toString(), absoluteScore);
        } else {
            output = String.format("\t%s\t%.2f\t%s\t%.1f\t", searchCandidates.get(bestCandidateIndex).toString(), scoresVsBestCandidate[nextBestCandidateIndex], searchCandidates.get(nextBestCandidateIndex).toString(), absoluteScore);
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
     * @param yFragments     array of Y fragments with spectrum intensities already matched
     * @param oxoFragments   array of oxonium fragments with spectrum intensities already matched
     * @param pepMass        peptide neutral mass
     * @param deltaMass      observed delta mass
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     * @return output probability score (sum of log ratios)
     */
    public double pairwiseCompareGlycans(GlycanCandidate glycan1, GlycanCandidate glycan2, GlycanFragment[] yFragments, GlycanFragment[] oxoFragments, double pepMass, double deltaMass, double massErrorWidth, double meanMassError) {
        double sumLogRatio = 0;
        // Y ions
        if (normYions) {
            sumLogRatio += pairwiseCompareFragmentsNorm(glycan1, glycan2, yFragments);
        } else {
            sumLogRatio += pairwiseCompareFragments(glycan1, glycan2, yFragments);
        }

        // oxonium ions
        sumLogRatio += pairwiseCompareFragments(glycan1, glycan2, oxoFragments);

        // isotope and mass errors
        sumLogRatio += determineIsotopeAndMassErrorProbs(glycan1, glycan2, deltaMass, massErrorWidth, meanMassError);

        return sumLogRatio;
    }

    /**
     * Compute sum log probability ratios for the compared glycans for a particular fragment type with normalization
     * of misses.
     * NOTE: does NOT allow fragment specific probabilities (for missed ions)
     *
     * @param glycan1       candidate 1
     * @param glycan2       candidate 2
     * @param fragments     array of fragments with spectrum intensities already matched
     * @return sum log probability with normalization included
     */
    public double pairwiseCompareFragmentsNorm(GlycanCandidate glycan1, GlycanCandidate glycan2, GlycanFragment[] fragments) {
        // probabilities calculated as in non-normed method, except for misses, which are normalized
        int cand1Misses = 0;
        int cand2Misses = 0;
        int cand1Hits = 0;
        int cand2Hits = 0;
        double sumLogRatio = 0;
        if (fragments.length == 0){
            return sumLogRatio;
        }
        for (GlycanFragment matchedFragment : fragments) {
            boolean foundInSpectrum = matchedFragment.foundIntensity > 0;
            if (foundInSpectrum) {
                if (matchedFragment.isAllowedFragment(glycan1)) {
                    if (!matchedFragment.isAllowedFragment(glycan2)) {
                        cand1Hits++;
                    }
                } else {
                    if (matchedFragment.isAllowedFragment(glycan2)) {
                        // allowed in 2 but not 1. Found in spectrum. Prob is 1/found probability
                        cand2Hits++;
                    }
                }
            } else {
                if (matchedFragment.isAllowedFragment(glycan1)) {
                    if (!matchedFragment.isAllowedFragment(glycan2)) {
                        // allowed in glycan 1, but NOT glycan 2. Unique miss for glycan 1
                        cand1Misses++;
                    }
                } else {
                    if (matchedFragment.isAllowedFragment(glycan2)) {
                        // allowed in 2 but not 1. Not found in spectrum. Unique miss for glycan 2
                        cand2Misses++;
                    }
                }
            }
        }

        // add prob normalized for number of misses for each candidate
        double cand1MissProb = fragments[0].ruleProbabilities[1];
        double cand2MissProb = 1 / cand1MissProb;
        double cand1HitProb = fragments[0].ruleProbabilities[0];
        double cand2HitProb = 1 / cand1HitProb;
        sumLogRatio += Math.sqrt(cand1Misses) * Math.log(cand1MissProb);    // candidate 1 misses
        sumLogRatio += Math.sqrt(cand2Misses) * Math.log(cand2MissProb);    // candidate 2 misses
        sumLogRatio += Math.sqrt(cand1Hits) * Math.log(cand1HitProb);    // candidate 1 hits
        sumLogRatio += Math.sqrt(cand2Hits) * Math.log(cand2HitProb);    // candidate 2 hits

        return sumLogRatio;
    }

    /**
     * Compute sum log probability ratios for the compared glycans for a particular fragment type. Does NOT
     * normalize fragment miss rate. Allows fragment specific probabilities.
     *
     * @param glycan1       candidate 1
     * @param glycan2       candidate 2
     * @param fragments     array of fragments with spectrum intensities already matched
     * @return sum log probability with normalization included
     */
    public double pairwiseCompareFragments(GlycanCandidate glycan1, GlycanCandidate glycan2, GlycanFragment[] fragments) {
        double sumLogRatio = 0;
        for (GlycanFragment matchedFragment : fragments) {
            double probRatio;
            if (matchedFragment.foundIntensity > 0) {
                if (matchedFragment.isAllowedFragment(glycan1)) {
                    if (matchedFragment.isAllowedFragment(glycan2)) {
                        // allowed in both - not distinguishing
                        probRatio = 1.0;
                    } else {
                        // allowed in glycan 1, but NOT glycan 2. Found in spectrum. Position 0 in rules array
                        probRatio = matchedFragment.ruleProbabilities[0];
                    }
                } else {
                    if (matchedFragment.isAllowedFragment(glycan2)) {
                        // allowed in 2 but not 1. Found in spectrum. Prob is 1/found probability
                        probRatio = 1.0 / matchedFragment.ruleProbabilities[0];
                    } else {
                        // allowed in neither candidate - not distinguishing
                        probRatio = 1.0;
                    }
                }
            } else {
                if (matchedFragment.isAllowedFragment(glycan1)) {
                    if (matchedFragment.isAllowedFragment(glycan2)) {
                        // allowed in both, but not found in spectrum - not distinguishing
                        probRatio = 1.0;
                    } else {
                        // allowed in glycan 1, but NOT glycan 2. Not found in spectrum. Position 1 in rules array
                        probRatio = matchedFragment.ruleProbabilities[1];
                    }
                } else {
                    if (matchedFragment.isAllowedFragment(glycan2)) {
                        // allowed in 2 but not 1. Not found in spectrum. Prob is 1/not-found probability
                        probRatio = 1.0 / matchedFragment.ruleProbabilities[1];
                    } else {
                        // allowed in neither candidate. Not found in spectrum. Not distinguishing
                        probRatio = 1.0;
                    }
                }
            }
            sumLogRatio += Math.log(probRatio);
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
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     * @param meanMassError  mean mass error of non-delta mass peptides
     * @return probability ratio (glycan 1 over 2)
     */
    public double determineIsotopeAndMassErrorProbs(GlycanCandidate glycan1, GlycanCandidate glycan2, double deltaMass, double massErrorWidth, double meanMassError) {
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
            double massStDevs1 = (massError1 - meanMassError) / massErrorWidth;
            double massError2 = deltaMass - glycan2.monoisotopicMass - (roundedIso2 * AAMasses.averagineIsotopeMass);
            double massStDevs2 = (massError2 - meanMassError) / massErrorWidth;
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
        if (yFragments.length == 0){
            return sumLogRatio;
        }
        // Y ions - check if allowed for this composition and score if so (ignore if not)
        int hitCount = 0;
        int missCount = 0;
        int disallowedHitCount = 0;
        // probabilities MUST be the same for all Y ions for this method to work
        double hitProb = yFragments[0].ruleProbabilities[0];
        double missProb = yFragments[0].ruleProbabilities[1];

        for (GlycanFragment yFragment : yFragments) {
            boolean foundInSpectrum = yFragment.foundIntensity > 0;
            if (yFragment.isAllowedFragment(bestGlycan)) {
                if (foundInSpectrum) {
                    hitCount++;
                } else {
                    missCount++;
                }
            } else {
                // fragment not allowed for this glycan, but if found in spectrum is still evidence against this composition
                if (foundInSpectrum) {
                    // only allow target-target and decoy-decoy matches to affect absolute score
                    if ((bestGlycan.isDecoy && yFragment.isDecoy) || (!bestGlycan.isDecoy && !yFragment.isDecoy)) {
                        disallowedHitCount++;
                    }
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
    public double computeFragmentAbsoluteScore(GlycanCandidate bestGlycan, GlycanFragment[] fragments){
        double sumLogRatio = 0;
        for (GlycanFragment fragment : fragments) {
            boolean foundInSpectrum = fragment.foundIntensity > 0;
            double probRatio;
            if (fragment.isAllowedFragment(bestGlycan)) {
                if (foundInSpectrum) {
                    probRatio = fragment.ruleProbabilities[0];     // found in spectrum - ion supports this glycan
                } else {
                    probRatio = fragment.ruleProbabilities[1];     // not found in spectrum - ion does not support this glycan
                }
            } else {
                // fragment not allowed for this glycan, but if found in spectrum is still evidence against this composition
                if (foundInSpectrum) {
                    // only allow target-target and decoy-decoy matches to affect absolute score
                    if ((bestGlycan.isDecoy && fragment.isDecoy) || (!bestGlycan.isDecoy && !fragment.isDecoy)) {
                        probRatio = 1 / fragment.ruleProbabilities[0];
                    } else {
                        probRatio = 1.0;
                    }
                } else {
                    // not allowed, and not found in spectrum - ignore
                    probRatio = 1.0;
                }
            }
            sumLogRatio += Math.log(probRatio);
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
            sumLogRatio = computeFragmentAbsoluteScore(bestGlycan, yFragments);
        }
        sumLogRatio += computeFragmentAbsoluteScore(bestGlycan, oxoFragments);

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
     * @param deltaMass delta mass being searched
     * @param glycanDatabase list of glycan candidates
     * @param isotopesToSearch list of isotope errors to consider
     * @param ms1TolerancePPM MS1 tolerance to consider around delta mass and isotope errors
     * @return list of glycan candidates with masses within the delta mass + iso errors and tolerance
     */
    public ArrayList<GlycanCandidate> getMatchingGlycansByMass(double deltaMass, ArrayList<GlycanCandidate> glycanDatabase, int[] isotopesToSearch, double ms1TolerancePPM) {
        ArrayList<GlycanCandidate> matchingGlycans = new ArrayList<>();
        for (int isotope : isotopesToSearch) {
            // add isotope error, which is recorded as an increase relative to delta mass
            double isotopeCorrMass = deltaMass - (isotope * AAMasses.averagineIsotopeMass);
            double massRangeDa = isotopeCorrMass * 0.000001 * ms1TolerancePPM;
            for (GlycanCandidate glycan : glycanDatabase) {
                // see if mass within specified ranges
                if (glycan.monoisotopicMass >= isotopeCorrMass - massRangeDa && glycan.monoisotopicMass <= isotopeCorrMass + massRangeDa) {
                    // match. todo: check duplicates (could be if user inputs them)
                    matchingGlycans.add(glycan);
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
            capYIonIntensities[i] = spec.findIonNeutral(capYIons[i], Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"))); //todo simplify parameter calling
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
        boolean [] allowedPoses = parseAllowedPositions(seq, PTMShepherd.getParam("localization_allowed_res"));
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
                if (allowedPoses[j] == true)
                    mods[j] += dmass;
                scores[j] = spec.getHyper(seq, mods, ppmTol);
                //System.out.println(scores[j] + "score");
                frags[j] = spec.getFrags(seq, mods, ppmTol);
                if(frags[j] > maxFrags[i])
                    maxFrags[i] = frags[j];
                if(scores[j] > maxScores[i])
                    maxScores[i] = scores[j];
                if (allowedPoses[j] == true)
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

    private boolean[] parseAllowedPositions(String seq, String allowedReses) {
        boolean [] allowedPoses = new boolean[seq.length()];
        if (allowedReses.equals("all") || allowedReses.equals(""))
            Arrays.fill(allowedPoses, true);
        else {
            Arrays.fill(allowedPoses, false);
            for (int i = 0; i < seq.length(); i++) {
                for (int j = 0; j < allowedReses.length(); j++) {
                    if (seq.charAt(i) == allowedReses.charAt(j)) {
                        allowedPoses[i] = true;
                        break;
                    }
                }
            }
        }
        return allowedPoses;
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
                    profiles[i].records[cind].updateWithLine(sp);
                }
            }
        }
        in.close();
    }

    /**
     * Initialize array of all fragment ions to search from the list of possible candidates. Each candidate has
     * fragment Ys up to the max HexNAc, Hex, and dHex present in the candidate. Duplicates are ignored (only 1
     * instance of a fragment in the list). Decoy fragments are generated for decoy candidates.
     * @param glycanCandidates input glycan database
     * @return array of GlycanFragments to search - all Y and oxonium ion
     */
    public GlycanFragment[] initializeYFragments(ArrayList<GlycanCandidate> glycanCandidates) {
        // Initialize list of all Y fragments to consider. Currently using only HexNAc, Hex, and dHex in Y ions
        ArrayList<GlycanFragment> yFragments = new ArrayList<>();
        HashMap<String, Boolean> fragmentsInList = new HashMap<>();     // record generated fragments to prevent adding duplicates
        for (GlycanCandidate candidate : glycanCandidates) {
            for (int hexnac = 0; hexnac <= candidate.glycanComposition.get(GlycanResidue.HexNAc); hexnac++) {
                for (int hex = 0; hex <= candidate.glycanComposition.get(GlycanResidue.Hex); hex++) {
                    if (hexnac == 0 && hex == 0) {
                        continue;
                    }
                    // add "regular" (no dHex) Y fragment for this HexNAc/Hex combination
                    Map<GlycanResidue, Integer> composition = new HashMap<>();
                    composition.put(GlycanResidue.HexNAc, hexnac);
                    composition.put(GlycanResidue.Hex, hex);
                    GlycanFragment fragment = new GlycanFragment(composition, probabilityTable.regularYrules, candidate.isDecoy);
                    // only add if not already in list (prevent duplicates)
                    if (!fragmentsInList.containsKey(fragment.toHashString())) {
                        yFragments.add(fragment);
                        fragmentsInList.put(fragment.toHashString(), true);
                    }
                    for (int dHex = 1; dHex <= candidate.glycanComposition.get(GlycanResidue.dHex); dHex++) {
                        // add dHex fragments (if allowed)
                        Map<GlycanResidue, Integer> dHexcomposition = new HashMap<>();
                        dHexcomposition.put(GlycanResidue.HexNAc, hexnac);
                        dHexcomposition.put(GlycanResidue.Hex, hex);
                        dHexcomposition.put(GlycanResidue.dHex, dHex);
                        GlycanFragment dHexfragment = new GlycanFragment(dHexcomposition, probabilityTable.dHexYrules, candidate.isDecoy);
                        if (!fragmentsInList.containsKey(dHexfragment.toHashString())) {
                            yFragments.add(dHexfragment);
                            fragmentsInList.put(dHexfragment.toHashString(), true);
                        }
                    }
                }
            }
        }
        return yFragments.toArray(new GlycanFragment[0]);
    }

    /**
     * Helper method to initialize hard-coded oxonium fragment rules. Only initializes fragments for a
     * residue type if at least one candidate contains that residue type (no need to consider if not).
     * Decoys generated for all residue types that have at least one decoy candidate containing that type.
     * @return list of GlycanFragments for oxonium ions
     */
    public GlycanFragment[] initializeOxoniumFragments(ArrayList<GlycanCandidate> searchCandidates){
        boolean targetNeuAc = false;
        boolean decoyNeuAc = false;
        boolean targetNeuGc = false;
        boolean decoyNeuGc = false;
        boolean targetPhospho = false;
        boolean decoyPhospho = false;
        boolean targetSulfo = false;
        boolean decoySulfo = false;

        for (GlycanCandidate candidate : searchCandidates) {
            if (candidate.isDecoy) {
                if (candidate.glycanComposition.get(GlycanResidue.NeuAc) > 0) {
                    decoyNeuAc = true;
                }
                if (candidate.glycanComposition.get(GlycanResidue.NeuGc) > 0) {
                    decoyNeuGc = true;
                }
                if (candidate.glycanComposition.get(GlycanResidue.Phospho) > 0) {
                    decoyPhospho = true;
                }
                if (candidate.glycanComposition.get(GlycanResidue.Sulfo) > 0) {
                    decoySulfo = true;
                }
            } else {
                if (candidate.glycanComposition.get(GlycanResidue.NeuAc) > 0) {
                    targetNeuAc = true;
                }
                if (candidate.glycanComposition.get(GlycanResidue.NeuGc) > 0) {
                    targetNeuGc = true;
                }
                if (candidate.glycanComposition.get(GlycanResidue.Phospho) > 0) {
                    targetPhospho = true;
                }
                if (candidate.glycanComposition.get(GlycanResidue.Sulfo) > 0) {
                    targetSulfo = true;
                }
            }
        }

        ArrayList<GlycanFragment> oxoniumList = new ArrayList<>();
        // HexNAc, Hex oxoniums

        // NeuAc
        if (targetNeuAc) {
            Map<GlycanResidue, Integer> neuacComposition = new HashMap<>();
            neuacComposition.put(GlycanResidue.NeuAc, 1);
            oxoniumList.add(new GlycanFragment(neuacComposition, probabilityTable.neuacRules, 273.0848565, false));     // NeuAc - H20
            oxoniumList.add(new GlycanFragment(neuacComposition, probabilityTable.neuacRules, 291.0954165, false));     // NeuAc
            Map<GlycanResidue, Integer> neuacHexComposition = new HashMap<>();
            neuacHexComposition.put(GlycanResidue.NeuAc, 1);
            neuacHexComposition.put(GlycanResidue.Hex, 1);
            neuacHexComposition.put(GlycanResidue.HexNAc, 1);
            oxoniumList.add(new GlycanFragment(neuacHexComposition, probabilityTable.neuacRules, 656.227624, false));     // NeuAc + HexNAc + Hex
        }
        if (decoyNeuAc) {
            Map<GlycanResidue, Integer> neuacComposition = new HashMap<>();
            neuacComposition.put(GlycanResidue.NeuAc, 1);
            oxoniumList.add(new GlycanFragment(neuacComposition, probabilityTable.neuacRules, 273.0848565, true));     // NeuAc - H20
            oxoniumList.add(new GlycanFragment(neuacComposition, probabilityTable.neuacRules, 291.0954165, true));     // NeuAc
            Map<GlycanResidue, Integer> neuacHexComposition = new HashMap<>();
            neuacHexComposition.put(GlycanResidue.NeuAc, 1);
            neuacHexComposition.put(GlycanResidue.Hex, 1);
            neuacHexComposition.put(GlycanResidue.HexNAc, 1);
            oxoniumList.add(new GlycanFragment(neuacHexComposition, probabilityTable.neuacRules, 656.227624, true));     // NeuAc + HexNAc + Hex
        }
        // NeuGc
        if (targetNeuGc) {
            Map<GlycanResidue, Integer> neugcComposition = new HashMap<>();
            neugcComposition.put(GlycanResidue.NeuGc, 1);
            oxoniumList.add(new GlycanFragment(neugcComposition, probabilityTable.neugcRules, 291.0954165, false));     // NeuGc - H20
            oxoniumList.add(new GlycanFragment(neugcComposition, probabilityTable.neugcRules, 307.090334, false));     // NeuGc
            Map<GlycanResidue, Integer> neugcHexComposition = new HashMap<>();
            neugcHexComposition.put(GlycanResidue.NeuGc, 1);
            neugcHexComposition.put(GlycanResidue.Hex, 1);
            neugcHexComposition.put(GlycanResidue.HexNAc, 1);
            oxoniumList.add(new GlycanFragment(neugcHexComposition, probabilityTable.neugcRules, 672.222524, false));     // NeuGc + HexNAc + Hex
        }
        if (decoyNeuGc) {
            Map<GlycanResidue, Integer> neugcComposition = new HashMap<>();
            neugcComposition.put(GlycanResidue.NeuGc, 1);
            oxoniumList.add(new GlycanFragment(neugcComposition, probabilityTable.neugcRules, 291.0954165, true));     // NeuGc - H20
            oxoniumList.add(new GlycanFragment(neugcComposition, probabilityTable.neugcRules, 307.090334, true));     // NeuGc
            Map<GlycanResidue, Integer> neugcHexComposition = new HashMap<>();
            neugcHexComposition.put(GlycanResidue.NeuGc, 1);
            neugcHexComposition.put(GlycanResidue.Hex, 1);
            neugcHexComposition.put(GlycanResidue.HexNAc, 1);
            oxoniumList.add(new GlycanFragment(neugcHexComposition, probabilityTable.neugcRules, 672.222524, true));     // NeuGc + HexNAc + Hex
        }
        // Phospho-Hex
        if (targetPhospho) {
            Map<GlycanResidue, Integer> phosphoHexComposition = new HashMap<>();
            phosphoHexComposition.put(GlycanResidue.Hex, 1);
            phosphoHexComposition.put(GlycanResidue.Phospho, 1);
            oxoniumList.add(new GlycanFragment(phosphoHexComposition, probabilityTable.phosphoRules, 242.01915, false));
            Map<GlycanResidue, Integer> phospho2HexComposition = new HashMap<>();
            phospho2HexComposition.put(GlycanResidue.Hex, 2);
            phospho2HexComposition.put(GlycanResidue.Phospho, 1);
            oxoniumList.add(new GlycanFragment(phospho2HexComposition, probabilityTable.phosphoRules, 404.07197, false));
            Map<GlycanResidue, Integer> twoPhosphoHexComposition = new HashMap<>();
            twoPhosphoHexComposition.put(GlycanResidue.Hex, 2);
            twoPhosphoHexComposition.put(GlycanResidue.Phospho, 2);
            oxoniumList.add(new GlycanFragment(twoPhosphoHexComposition, probabilityTable.phosphoRules, 484.0383, false));
        }
        if (decoyPhospho) {
            Map<GlycanResidue, Integer> phosphoHexComposition = new HashMap<>();
            phosphoHexComposition.put(GlycanResidue.Hex, 1);
            phosphoHexComposition.put(GlycanResidue.Phospho, 1);
            oxoniumList.add(new GlycanFragment(phosphoHexComposition, probabilityTable.phosphoRules, 242.01915, true));
            Map<GlycanResidue, Integer> phospho2HexComposition = new HashMap<>();
            phospho2HexComposition.put(GlycanResidue.Hex, 2);
            phospho2HexComposition.put(GlycanResidue.Phospho, 1);
            oxoniumList.add(new GlycanFragment(phospho2HexComposition, probabilityTable.phosphoRules, 404.07197, true));
            Map<GlycanResidue, Integer> twoPhosphoHexComposition = new HashMap<>();
            twoPhosphoHexComposition.put(GlycanResidue.Hex, 2);
            twoPhosphoHexComposition.put(GlycanResidue.Phospho, 2);
            oxoniumList.add(new GlycanFragment(twoPhosphoHexComposition, probabilityTable.phosphoRules, 484.0383, true));
        }
        // Sulfo
        if (targetSulfo) {
            Map<GlycanResidue, Integer> sulfoComposition = new HashMap<>();
            sulfoComposition.put(GlycanResidue.HexNAc, 1);
            sulfoComposition.put(GlycanResidue.Sulfo, 1);
            oxoniumList.add(new GlycanFragment(sulfoComposition, probabilityTable.sulfoRules, 283.036724, false));
            Map<GlycanResidue, Integer> sulfoHexNAcHexComposition = new HashMap<>();
            sulfoHexNAcHexComposition.put(GlycanResidue.HexNAc, 1);
            sulfoHexNAcHexComposition.put(GlycanResidue.Hex, 1);
            sulfoHexNAcHexComposition.put(GlycanResidue.Sulfo, 1);
            oxoniumList.add(new GlycanFragment(sulfoHexNAcHexComposition, probabilityTable.sulfoRules, 445.090724, false));
            Map<GlycanResidue, Integer> sulfoTwoHexNAcHexComposition = new HashMap<>();
            sulfoTwoHexNAcHexComposition.put(GlycanResidue.HexNAc, 2);
            sulfoHexNAcHexComposition.put(GlycanResidue.Hex, 2);
            sulfoTwoHexNAcHexComposition.put(GlycanResidue.Sulfo, 1);
            oxoniumList.add(new GlycanFragment(sulfoTwoHexNAcHexComposition, probabilityTable.sulfoRules, 810.204724, false));
        }
        if (decoySulfo) {
            Map<GlycanResidue, Integer> sulfoComposition = new HashMap<>();
            sulfoComposition.put(GlycanResidue.HexNAc, 1);
            sulfoComposition.put(GlycanResidue.Sulfo, 1);
            oxoniumList.add(new GlycanFragment(sulfoComposition, probabilityTable.sulfoRules, 283.036724, true));
            Map<GlycanResidue, Integer> sulfoHexNAcHexComposition = new HashMap<>();
            sulfoHexNAcHexComposition.put(GlycanResidue.HexNAc, 1);
            sulfoHexNAcHexComposition.put(GlycanResidue.Hex, 1);
            sulfoHexNAcHexComposition.put(GlycanResidue.Sulfo, 1);
            oxoniumList.add(new GlycanFragment(sulfoHexNAcHexComposition, probabilityTable.sulfoRules, 445.090724, true));
            Map<GlycanResidue, Integer> sulfoTwoHexNAcHexComposition = new HashMap<>();
            sulfoTwoHexNAcHexComposition.put(GlycanResidue.HexNAc, 2);
            sulfoHexNAcHexComposition.put(GlycanResidue.Hex, 2);
            sulfoTwoHexNAcHexComposition.put(GlycanResidue.Sulfo, 1);
            oxoniumList.add(new GlycanFragment(sulfoTwoHexNAcHexComposition, probabilityTable.sulfoRules, 810.204724, true));
        }
        // dHex

        return oxoniumList.toArray(new GlycanFragment[0]);
    }

}
