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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

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

    public GlycoAnalysis(String dsName, ArrayList<GlycanCandidate> glycoDatabase, ProbabilityTables inputProbabilityTable) {
        this.dsName = dsName;
        this.glycoFile = new File(PTMShepherd.normFName(dsName+".rawglyco"));
        this.glycanDatabase = glycoDatabase;
        this.probabilityTable = inputProbabilityTable;    // init with default values, can be changed by params
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
        for(int i = 0; i < capYstrs.length; i++)
            capYShifts[i] = Double.parseDouble(capYstrs[i]);
        //oxonium ions
        String[] oxStrs;
        if (PTMShepherd.getParam("diag_ions").length() > 0)
            oxStrs = PTMShepherd.getParam("diag_ions").split(",| |/");
        else
            oxStrs = new String[0];
        oxoniumIons = new double[oxStrs.length];
        for(int i = 0; i < oxStrs.length; i++)
            oxoniumIons[i] = Double.parseDouble(oxStrs[i]);
        //remainder masses
        String[] remainderStrs;
        if (PTMShepherd.getParam("remainder_masses").length() > 0)
            remainderStrs = PTMShepherd.getParam("remainder_masses").split(",| |/");
        else
            remainderStrs = new String[0];
        remainderMasses = new double[remainderStrs.length];
        for(int i = 0; i < remainderStrs.length; i++)
            remainderMasses[i] = Double.parseDouble(remainderStrs[i]);

        //write header
        StringBuffer headbuff = new StringBuffer(String.format("%s\t%s\t%s\t%s\t%s", "Spectrum", "Peptide", "Mods", "Pep Mass", "Mass Shift"));
        headbuff.append("\tBest Glycan\tLog Delta Score\t2nd Best Glycan");
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
            PTMShepherd.print(String.format("\t%s - %d (%d ms, %d ms)", cf, clines.size(), t2-t1,t3-t2));
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
            String [] sp = line.split("\\t");
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
        for (int i=0; i < binCounts.length; i++) {
            // x-value is minError + binSize*i + binSize/2 (middle of the bin), y-value is counts
            double xval = minError + i*binSize + binSize / 2.0;
            massErrorObservations.add(xval, binCounts[i]);
        }
        double[] fitParameters = fitter.fit(massErrorObservations.toList());
        meanMassError = fitParameters[1];   // output is amplitude, mean, sigma of fitted curve
        massErrorWidth = fitParameters[2];
    }


    public StringBuffer processLine(String line) {
    public String processLine(String line) {
        StringBuffer sb = new StringBuffer();
        String [] sp = line.split("\\t");
        String seq = sp[pepCol];
        float dmass = Float.parseFloat(sp[deltaCol]);
        float pepMass = Float.parseFloat(sp[pmassCol]);
        String[] smods = sp[modCol].split(",");
        String specName = sp[specCol];

        sb.append(String.format("%s\t%s\t%s\t%.4f\t%.4f", specName,seq,sp[modCol],pepMass, dmass));

        Spectrum spec = mr.getSpectrum(reNormName(specName));
        if (spec == null) {
            this.lineWithoutSpectra.add(reNormName(specName));
            return "ERROR";
        }
        spec.conditionOptNorm(condPeaks, condRatio, false);

        sb.append(assignGlycanToPSM(spec, pepMass, dmass, glycanDatabase, massErrorWidth, meanMassError));

        //System.out.println("got spec");
        double [] capYIonIntensities;
        double [] oxoniumIonIntensities;

        capYIonIntensities = findCapitalYIonMasses(spec, pepMass);
        oxoniumIonIntensities = findOxoniumIonMasses(spec, pepMass);

        for (int i = 0; i < capYIonIntensities.length; i++)
            sb.append(String.format("\t%.2f", capYIonIntensities[i]));
        for (int i = 0; i < oxoniumIonIntensities.length; i++)
            sb.append(String.format("\t%.2f", oxoniumIonIntensities[i]));

        float [] deltaScores = new float[remainderMasses.length];
        boolean [][] isMaxScores = localizeRemainderFragments(spec, sp[pepCol], smods, deltaScores);

        for (int i = 0; i < remainderMasses.length; i++) {
            sb.append(String.format("\t%.1f", deltaScores[i]));
            StringBuffer locSb = new StringBuffer("\t");
            for (int j = 0; j < seq.length(); j++) {
                if (isMaxScores[i][j] == true) {
                    locSb.append(String.format("%d%c",j+1, seq.charAt(j))); //position (1 indexed), character
                }
            }
            sb.append(locSb.toString());
        }
        return sb.toString();
    }

    /**
     * Main glycan assignment method at PSM level. Searches Y/Oxonium ions (and eventually exact mass/isotope) to compare
     * to possible glycan candidates. Goal is to return best glycan candidate and score.
     * @param spec spectrum being searched
     * @param pepMass peptide mass (without glycan)
     * @param glycanDatabase possible glycan candidates
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     * @param deltaMass observed delta mass from PSM
     */
    public String assignGlycanToPSM(Spectrum spec, double pepMass, double deltaMass, ArrayList<GlycanCandidate> glycanDatabase, double massErrorWidth, double meanMassError) {
        // skip non-delta mass PSMs
        if (deltaMass < 3.5 && deltaMass > -1.5) {
            return "\t\t\t";
        }

        // Determine possible glycan candidates from mass
        int[] isotopesToSearch = {-1, 0, 1, 2, 3};
        double ms1TolPPM = 50;  // todo: connect to existing param?
        ArrayList<GlycanCandidate> searchCandidates = getMatchingGlycansByMass(deltaMass, glycanDatabase, isotopesToSearch, ms1TolPPM);

        // Get Y ions possible for these candidates
        // todo: move these out to use same ions for all? Or keep in here to search specific peaks for each candidate set?
        GlycanFragment[] possibleYIons = initializeYFragments(searchCandidates);
        GlycanFragment[] possibleOxoniums = initializeOxoniumFragments();
        double[] yMasses = new double[possibleYIons.length];
        for (int i=0; i < possibleYIons.length; i++) {
            yMasses[i] = possibleYIons[i].neutralMass + pepMass;
        }
        double[] oxoMasses = new double[possibleOxoniums.length];
        for (int i=0; i < possibleOxoniums.length; i++){
            oxoMasses[i] = possibleOxoniums[i].neutralMass;
        }

        // Search Y and oxonium ions in spectrum
        float ppmTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
        for (int i=0; i < possibleYIons.length; i++) {
            possibleYIons[i].foundIntensity = spec.findIonNeutral(yMasses[i], ppmTol);  // sum of charge state intensities if >1 found
        }
        for (int i=0; i < possibleOxoniums.length; i++) {
            // only search 1+ oxonium ions, not all possible charge states
            possibleOxoniums[i].foundIntensity = spec.findIon(oxoMasses[i] + AAMasses.protMass, ppmTol);
        }

        // score candidates and save results
        int bestCandidateIndex = 0;
        int nextBestCandidateIndex = 1;
        double[] scoresVsBestCandidate = new double[searchCandidates.size()];
        // todo: need to do something if only 1 candidate

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
                for (int j=0; j < i; j++) {
                    scoresVsBestCandidate[j] -= comparisonScore;    // subtract score vs previous best candidate from existing scores to update them
                }
                nextBestCandidateIndex = bestCandidateIndex;
                bestCandidateIndex = i;
                scoresVsBestCandidate[i] = 0;
            }
        }

        // output - best glycan, scores, etc back to PSM table
        String output;
        if (searchCandidates.size() == 0) {
            output = "\tNo Matches\t\t";
        } else if (searchCandidates.size() == 1) {
            output = String.format("\t%s\t\t", searchCandidates.get(bestCandidateIndex).toString());
        } else {
            output = String.format("\t%s\t%.2f\t%s", searchCandidates.get(bestCandidateIndex).toString(), scoresVsBestCandidate[nextBestCandidateIndex], searchCandidates.get(nextBestCandidateIndex).toString());
        }
        return output;
    }

    /**
     * Perform pairwise comparison of two glycans. Uses sum of log probability ratios between candidates for
     * each category (mass/iso error and fragment ion) being considered. Returns a single score of combined
     * probability of first glycan candidate over second.
     * @param glycan1 candidate 1
     * @param glycan2 candidate 2
     * @param yFragments array of Y fragments with spectrum intensities already matched
     * @param oxoFragments array of oxonium fragments with spectrum intensities already matched
     * @param pepMass peptide neutral mass
     * @param deltaMass observed delta mass
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     * @return output probability score (sum of log ratios)
     */
    public double pairwiseCompareGlycans(GlycanCandidate glycan1, GlycanCandidate glycan2, GlycanFragment[] yFragments, GlycanFragment[] oxoFragments, double pepMass, double deltaMass, double massErrorWidth, double meanMassError) {
        double sumLogRatio = 0;
        // Y ions
        for (GlycanFragment yFragment : yFragments) {
            boolean foundInSpectrum = yFragment.foundIntensity > 0;
            double probRatio = determineProbRatio(yFragment, glycan1, glycan2, foundInSpectrum);
            sumLogRatio += Math.log(probRatio);
        }
        // oxonium ions
        for (GlycanFragment oxoFragment : oxoFragments) {
            boolean foundInSpectrum = oxoFragment.foundIntensity > 0;
            double probRatio = determineProbRatio(oxoFragment, glycan1, glycan2, foundInSpectrum);
            sumLogRatio += Math.log(probRatio);
        }

        // isotope and mass errors
        sumLogRatio += determineIsotopeAndMassErrorProbs(glycan1, glycan2, deltaMass, massErrorWidth, meanMassError);

        return sumLogRatio;
    }

    /**
     * Determine the probability ratio for this pairwise comparison based on isotope error. Currently
     * uses hard-coded isotope probabilities, but could be updated to get rate from dataset
     * @param glycan1 glycan 1
     * @param glycan2 glycan 2
     * @param deltaMass observed delta mass
     * @param massErrorWidth Width of the mass error distribution for non-delta mass peptides to use for determining probability of glycan candidates
     * @return probability ratio (glycan 1 over 2)
     */
    public double determineIsotopeAndMassErrorProbs(GlycanCandidate glycan1, GlycanCandidate glycan2, double deltaMass, double massErrorWidth, double meanMassError) {
        // Determine isotopes
        float iso1 = (float) (glycan1.monoisotopicMass - deltaMass);
        int roundedIso1 = Math.round(iso1);
        float iso2 = (float) (glycan2.monoisotopicMass - deltaMass);
        int roundedIso2 = Math.round(iso2);

        double isotopeProbRatio = probabilityTable.isotopeProbTable.get(roundedIso1) / probabilityTable.isotopeProbTable.get(roundedIso2);

        // mass error calc
        double massError1 = deltaMass - glycan1.monoisotopicMass - (roundedIso1 * AAMasses.averagineIsotopeMass);
        double massStDevs1 = (massError1 - meanMassError) / massErrorWidth;
        double massError2 = deltaMass - glycan2.monoisotopicMass - (roundedIso2 * AAMasses.averagineIsotopeMass);
        double massStDevs2 = (massError2 - meanMassError) / massErrorWidth;
        double massProbRatio = Math.abs(massStDevs2 / massStDevs1) * probabilityTable.massProbScaling;     // divide #2 by #1 to get ratio for likelihood of #1 vs #2, adjust by scaling factor

        return Math.log(isotopeProbRatio) + Math.log(massProbRatio);
    }

    /**
     * Given a fragment matched in the spectrum, determine the ratio of probabilities for a pair of glycan candidates
     * by whether that fragment is allowed in both, 1 but not 2, 2 but not 1, or neither.
     * @param matchedFragment fragment of interest
     * @param glycan1 candidate 1
     * @param glycan2 candidate 2
     * @return probability ratio
     */
    public double determineProbRatio(GlycanFragment matchedFragment, GlycanCandidate glycan1, GlycanCandidate glycan2, boolean foundInSpectrum) {
        double probRatio;
        if (foundInSpectrum){
            if (matchedFragment.isAllowedFragment(glycan1.glycanComposition)) {
                if (matchedFragment.isAllowedFragment(glycan2.glycanComposition)) {
                    // allowed in both - not distinguishing
                    probRatio = 1.0;
                } else {
                    // allowed in glycan 1, but NOT glycan 2. Found in spectrum. Position 0 in rules array
                    probRatio = matchedFragment.ruleProbabilities[0];
                }
            } else {
                if (matchedFragment.isAllowedFragment(glycan2.glycanComposition)) {
                    // allowed in 2 but not 1. Found in spectrum. Prob is 1/found probability
                    probRatio = 1.0 / matchedFragment.ruleProbabilities[0];
                } else {
                    // allowed in neither candidate - not distinguishing
                    probRatio = 1.0;
                }
            }
        } else {
            if (matchedFragment.isAllowedFragment(glycan1.glycanComposition)) {
                if (matchedFragment.isAllowedFragment(glycan2.glycanComposition)) {
                    // allowed in both, but not found in spectrum - not distinguishing
                    probRatio = 1.0;
                } else {
                    // allowed in glycan 1, but NOT glycan 2. Not found in spectrum. Position 1 in rules array
                    probRatio = matchedFragment.ruleProbabilities[1];
                }
            } else {
                if (matchedFragment.isAllowedFragment(glycan2.glycanComposition)) {
                    // allowed in 2 but not 1. Not found in spectrum. Prob is 1/not-found probability
                    probRatio = 1.0 / matchedFragment.ruleProbabilities[1];
                } else {
                    // allowed in neither candidate. Not found in spectrum. Not distinguishing
                    probRatio = 1.0;
                }
            }
        }
        return probRatio;
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
            double isotopeCorrMass = deltaMass + (isotope * AAMasses.averagineIsotopeMass);
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
     * Generate a list of all possible Y ion shifts to look for, given the list of glycan candidates
     * being considered. Output list is the set of all Y ion shifts in all candidates without duplicates.
     * @param glycanCandidates list of glycan candidates
     * @return array of Y shifts
     */
    public double[] getAllPossibleYShifts(ArrayList<GlycanCandidate> glycanCandidates) {
        ArrayList<Double> yShifts = new ArrayList<>();
        for (GlycanCandidate glycan : glycanCandidates) {
            for (double yMass : glycan.expectedYIons) {
                if (!yShifts.contains(yMass)) {
                    yShifts.add(yMass);
                }
            }
            // need to look for disallowed Y ions in spectra as well
            for (double yMass : glycan.disallowedYIons) {
                if (!yShifts.contains(yMass)) {
                    yShifts.add(yMass);
                }
            }
        }

        // convert to array
        final double[] allYShifts = new double[yShifts.size()];
        for (int i=0; i < allYShifts.length; i++) {
            allYShifts[i] = yShifts.get(i);
        }
        return allYShifts;
    }

    /**
     * Initialize array of all fragment ions to search from the list of possible candidates.
     * Determines all Y ions (all combinations up to max number of a given residue in all candidates,
     * plus extra 'disallowed' possible number)
     * @param glycanCandidates input glycan database
     * @return array of GlycanFragments to search - all Y and oxonium ion
     */
    public GlycanFragment[] initializeYFragments(ArrayList<GlycanCandidate> glycanCandidates) {
        // init hashmap with all types of glycans
        HashMap<GlycanResidue, Integer> maxGlycanResidues = new HashMap<>();
        for (GlycanResidue residue : GlycanResidue.values()){
            maxGlycanResidues.put(residue, 0);
        }

        // determine maximum number of each residue type in any candidate being considered
        for (GlycanCandidate glycanCandidate : glycanCandidates) {
            for (Map.Entry<GlycanResidue, Integer> residue : glycanCandidate.glycanComposition.entrySet()) {
                if (maxGlycanResidues.get(residue.getKey()) < residue.getValue()) {
                    // new max found, update
                    maxGlycanResidues.put(residue.getKey(), residue.getValue());
                }
            }
        }

        // Initialize list of all Y fragments to consider. Currently using only HexNAc, Hex, and dHex in Y ions
        ArrayList<GlycanFragment> yFragments = new ArrayList<>();
        for (int hexnac=0; hexnac <= maxGlycanResidues.get(GlycanResidue.HexNAc); hexnac++) {
            for (int hex=0; hex <= maxGlycanResidues.get(GlycanResidue.Hex); hex++) {
                if (hexnac == 0 && hex == 0) {
                    continue;
                }
                // add "regular" (no dHex) Y fragment for this HexNAc/Hex combination
                Map<GlycanResidue, Integer> composition = new HashMap<>();
                composition.put(GlycanResidue.HexNAc, hexnac);
                composition.put(GlycanResidue.Hex, hex);
                GlycanFragment fragment = new GlycanFragment(composition, probabilityTable.regularYrules);
                yFragments.add(fragment);
                for (int dHex=1; dHex <= maxGlycanResidues.get(GlycanResidue.dHex); dHex++) {
                    // add dHex fragments (if allowed)
                    Map<GlycanResidue, Integer> dHexcomposition = new HashMap<>();
                    dHexcomposition.put(GlycanResidue.HexNAc, hexnac);
                    dHexcomposition.put(GlycanResidue.Hex, hex);
                    dHexcomposition.put(GlycanResidue.dHex, dHex);
                    GlycanFragment dHexfragment = new GlycanFragment(dHexcomposition, probabilityTable.dHexYrules);
                    yFragments.add(dHexfragment);
                }
            }
        }
        return yFragments.toArray(new GlycanFragment[0]);
    }

    /**
     * Helper method to initialize hard-coded oxonium fragment rules.
     * @return list of GlycanFragments for oxonium ions
     */
    public GlycanFragment[] initializeOxoniumFragments(){
        ArrayList<GlycanFragment> oxoniumList = new ArrayList<>();
        // HexNAc, Hex oxoniums

        // NeuAc
        Map<GlycanResidue, Integer> neuacComposition = new HashMap<>();
        neuacComposition.put(GlycanResidue.NeuAc, 1);
        oxoniumList.add(new GlycanFragment(neuacComposition, probabilityTable.neuacRules, 273.0848565));     // NeuAc - H20
        oxoniumList.add(new GlycanFragment(neuacComposition, probabilityTable.neuacRules, 291.0954165));     // NeuAc
        Map<GlycanResidue, Integer> neuacHexComposition = new HashMap<>();
        neuacHexComposition.put(GlycanResidue.NeuAc, 1);
        neuacHexComposition.put(GlycanResidue.Hex, 1);
        neuacHexComposition.put(GlycanResidue.HexNAc, 1);
        oxoniumList.add(new GlycanFragment(neuacHexComposition, probabilityTable.neuacRules, 656.227624));     // NeuAc + HexNAc + Hex

        // NeuGc
        Map<GlycanResidue, Integer> neugcComposition = new HashMap<>();
        neugcComposition.put(GlycanResidue.NeuGc, 1);
        oxoniumList.add(new GlycanFragment(neugcComposition, probabilityTable.neugcRules, 291.0954165));     // NeuGc - H20
        oxoniumList.add(new GlycanFragment(neugcComposition, probabilityTable.neugcRules, 307.090334));     // NeuGc
        Map<GlycanResidue, Integer> neugcHexComposition = new HashMap<>();
        neugcHexComposition.put(GlycanResidue.NeuGc, 1);
        neugcHexComposition.put(GlycanResidue.Hex, 1);
        neugcHexComposition.put(GlycanResidue.HexNAc, 1);
        oxoniumList.add(new GlycanFragment(neugcHexComposition, probabilityTable.neugcRules, 672.222524));     // NeuGc + HexNAc + Hex

        // Phospho-Hex
        Map<GlycanResidue, Integer> phosphoHexComposition = new HashMap<>();
        phosphoHexComposition.put(GlycanResidue.Hex, 1);
        phosphoHexComposition.put(GlycanResidue.Phospho, 1);
        oxoniumList.add(new GlycanFragment(phosphoHexComposition, probabilityTable.phosphoRules, 242.01915));
        Map<GlycanResidue, Integer> phospho2HexComposition = new HashMap<>();
        phospho2HexComposition.put(GlycanResidue.Hex, 2);
        phospho2HexComposition.put(GlycanResidue.Phospho, 1);
        oxoniumList.add(new GlycanFragment(phospho2HexComposition, probabilityTable.phosphoRules, 404.07197));
        Map<GlycanResidue, Integer> twoPhosphoHexComposition = new HashMap<>();
        twoPhosphoHexComposition.put(GlycanResidue.Hex, 2);
        twoPhosphoHexComposition.put(GlycanResidue.Phospho, 2);
        oxoniumList.add(new GlycanFragment(twoPhosphoHexComposition, probabilityTable.phosphoRules, 484.0383));

        // Sulfo
        Map<GlycanResidue, Integer> sulfoComposition = new HashMap<>();
        sulfoComposition.put(GlycanResidue.HexNAc, 1);
        sulfoComposition.put(GlycanResidue.Sulfo, 1);
        oxoniumList.add(new GlycanFragment(sulfoComposition, probabilityTable.sulfoRules, 283.036724));
        Map<GlycanResidue, Integer> sulfoHexNAcHexComposition = new HashMap<>();
        sulfoHexNAcHexComposition.put(GlycanResidue.HexNAc, 1);
        sulfoHexNAcHexComposition.put(GlycanResidue.Hex, 1);
        sulfoHexNAcHexComposition.put(GlycanResidue.Sulfo, 1);
        oxoniumList.add(new GlycanFragment(sulfoHexNAcHexComposition, probabilityTable.sulfoRules, 445.090724));
        Map<GlycanResidue, Integer> sulfoTwoHexNAcHexComposition = new HashMap<>();
        sulfoTwoHexNAcHexComposition.put(GlycanResidue.HexNAc, 2);
        sulfoHexNAcHexComposition.put(GlycanResidue.Hex, 2);
        sulfoTwoHexNAcHexComposition.put(GlycanResidue.Sulfo, 1);
        oxoniumList.add(new GlycanFragment(sulfoTwoHexNAcHexComposition, probabilityTable.sulfoRules, 810.204724));

        // dHex

        return oxoniumList.toArray(new GlycanFragment[0]);
    }

}
