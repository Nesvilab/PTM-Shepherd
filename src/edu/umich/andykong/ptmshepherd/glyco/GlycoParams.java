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

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
import umich.ms.glyco.*;

import java.io.*;
import java.util.*;

/**
 * Organization - putting a bunch of glyco-specific utilities here rather than cluttering the main PTM-S.java
 */
public class GlycoParams {

    public Random randomGenerator;
    public HashMap<String, GlycanResidue> glycanResiduesMap;
    public ArrayList<GlycanResidue> glycanResidues;
    public ArrayList<GlycanCandidate> glycoDatabase;
    public int decoyType;
    public double glycoPPMtol;
    public Integer[] glycoIsotopes;
    public boolean nGlycan;
    public boolean glycoYnorm;
    public double absScoreErrorParam;
    public double glycoFDR;
    public boolean printFullParams;
    public boolean writeGlycansToAssignedMods;
    public boolean removeGlycanDeltaMass;
    public boolean printGlycoDecoys;
    public int numThreads;
    public boolean useGlycanFragmentProbs;
    public boolean useNonCompFDR ;
    public double defaultProp;
    public String allowedLocalizationResidues;
    public HashMap<GlycanResidue, ArrayList<GlycanFragment>> glycoOxoniumDatabase;
    public HashMap<Integer, Double> isotopeProbTable;
    public double massProbScaling;
    public boolean glycoLDA;
    public double topPctSpectraForConsensus;
    public int minYsForConsensus;
    public int minPSMsForConsensus;
    // IonQuant MS1 params
    public float rtTol = 0.4f;
    public float imTol = 0.05f;
    public int minIsotopesIonQuant = 2;
    public int minScansIonQuant = 1;
    public boolean isIMdata = false;
    public ArrayList<LDAFeature> ldaFeaturesToUse;
    public double ldaTargetProp;
    public boolean useShuffledIntensities;
    public boolean simOnly2ndPass;
    public boolean useNormMassScore;

    private static final String defaultResiduePath = "glycan_residues.txt";
    private static final String defaultModsPath = "glycan_mods.txt";
    public static final String defaultOxoPath = "oxonium_ion_list.txt";

    public static final double MAX_CANDIDATE_DECOY_SHIFT_DA = 3;
    public static final double DEFAULT_PEPTIDE_MASS = 1500;

    public GlycoParams(String glycanResiduesPath, String glycanModsPath, String oxoniumListPath) {
        // parse the glycan residues and mods tables, using internal defaults if no paths provided from FragPipe or user
        if (glycanResiduesPath.isEmpty()) {
            glycanResiduesMap = GlycanParser.parseGlycoResiduesDBStream(GlycoParams.class.getResourceAsStream(defaultResiduePath));
        } else {
            glycanResiduesMap = GlycanParser.parseGlycoResiduesDB(glycanResiduesPath);
        }
        HashMap<String, GlycanMod> modsMap;
        if (glycanModsPath.isEmpty()) {
            modsMap = GlycanParser.parseGlycoModsDBStream(GlycoParams.class.getResourceAsStream(defaultModsPath), glycanResiduesMap.size(), glycanResiduesMap);
        } else {
            modsMap = GlycanParser.parseGlycoModsDB(glycanModsPath, glycanResiduesMap.size(), glycanResiduesMap);
        }
        glycanResiduesMap.putAll(modsMap);
        glycanResidues = new ArrayList<>(glycanResiduesMap.values());
        glycanResidues.sort(GlycanResidue::compareTo);

        if (oxoniumListPath.isEmpty()) {
            glycoOxoniumDatabase = GlycanParser.parseOxoDBStream(GlycoParams.class.getResourceAsStream(defaultOxoPath), glycanResiduesMap, randomGenerator);
        } else {
            glycoOxoniumDatabase = GlycanParser.parseOxoDB(oxoniumListPath, glycanResiduesMap, randomGenerator);
        }
        ldaFeaturesToUse = new ArrayList<>();
    }


    /**
     * Write a list of masses for all glycan candidates to pass to IonQuant. Format is one mass per line.
     * @param glycanDB list of glycan candidates (whole glycan database)
     * @param outputPath where to save the file
     */
    public static void writeGlycanMassList(ArrayList<GlycanCandidate> glycanDB, String outputPath) {
        try {
            PrintWriter out = new PrintWriter(new FileWriter(outputPath));
            HashSet<Integer> writtenMasses = new HashSet<>();
            for (GlycanCandidate candidate : glycanDB) {
                if (candidate.isDecoy) {
                    continue;       // do not write decoy masses to list - only target masses are reported in PSM table, even if decoy is assigned
                }
                int roundedMass = (int) Math.round(candidate.mass * 100);
                if (!writtenMasses.contains(roundedMass)) {
                    out.write(String.format("%.4f\n",candidate.mass));
                    writtenMasses.add(roundedMass);
                }
            }
            out.flush();
            out.close();
        } catch (IOException e) {
            PTMShepherd.die("Could not write glycan mass list to file " + outputPath + "\ndue to error: " + e.getMessage());
        }
    }

    /**
     * Parse FragPipe glycan database string to list of glycan candidates
     * @param glycanDBString
     * @return
     */
    public ArrayList<GlycanCandidate> parseGlycanDatabaseString(String glycanDBString) {
        ArrayList<Glycan> glycans = GlycanParser.parseGlycanDatabaseString(glycanDBString, glycanResiduesMap);
        return convertGlycansToCandidates(glycans, glycanResiduesMap, nGlycan, glycoOxoniumDatabase, decoyType, glycoPPMtol, glycoIsotopes, randomGenerator);
    }

    /**
     * Deprecated method, retained for backwards compatibility and command line access.
     * Parse input glycan database file. Formatting: 1 glycan per line, "Residue1-count_Residue2-count_...\n"
     * @param inputPath path to input file
     * @return list of glycans to consider
     */
    public ArrayList<GlycanCandidate> parseGlycanDatabaseFile(String inputPath) {
        ArrayList<Glycan> glycans = GlycanParser.loadGlycansFromText(inputPath, GlycanParser.detectDBtype(inputPath), glycanResiduesMap);
        return convertGlycansToCandidates(glycans, glycanResiduesMap, nGlycan, glycoOxoniumDatabase, decoyType, glycoPPMtol, glycoIsotopes, randomGenerator);
    }

    /**
     * Generate GlycanCandidates for composition assignment from Glycans, preventing duplicates and adding decoy Candidates.
     * @param glycans
     * @return
     */
    public static ArrayList<GlycanCandidate> convertGlycansToCandidates(ArrayList<Glycan> glycans,
                                                                        HashMap<String, GlycanResidue> glycanResiduesMap,
                                                                        boolean nGlycan,
                                                                        HashMap<GlycanResidue, ArrayList<GlycanFragment>> glycoOxoniumDatabase,
                                                                        int decoyType,
                                                                        double glycoPPMtol,
                                                                        Integer[] glycoIsotopes,
                                                                        Random randomGenerator) {
        HashMap<String, Boolean> glycansInDB = new HashMap<>();
        ArrayList<GlycanCandidate> glycanDB = new ArrayList<>();
        for (Glycan glycan: glycans) {
            if (glycan.composition.isEmpty()) {
                continue;
            }
            // generate a new candidate from this composition and add to DB
            GlycanCandidate candidate = GlycanCandidate.initGlycanCandidate(glycan.composition, 0,false, glycanResiduesMap, nGlycan, randomGenerator, glycoOxoniumDatabase);
            String compositionHash = candidate.toString();
            // prevent addition of duplicates if user has them in database
            if (!glycansInDB.containsKey(compositionHash)) {
                glycanDB.add(candidate);
                glycansInDB.put(compositionHash, Boolean.TRUE);
                // also add a decoy for this composition
                double decoyMassShift = setDecoyShift(candidate.mass, decoyType, glycoPPMtol, glycoIsotopes, randomGenerator);
                GlycanCandidate decoy = GlycanCandidate.initGlycanCandidate(glycan.composition, decoyMassShift, true, glycanResiduesMap, nGlycan, randomGenerator, glycoOxoniumDatabase);
                glycanDB.add(decoy);
            }
        }
        return glycanDB;
    }

    // Helper method for determining decoy masses for various decoy mass generation settings
    private static double setDecoyShift(double baseMass, int decoyType, double glycoPPMtol, Integer[] glycoIsotopes, Random randomGenerator) {
        double randomShift = 0;
        switch (decoyType) {
            case 0:
                // simple mass window
                randomShift = GlycanFragment.randomMassShift(MAX_CANDIDATE_DECOY_SHIFT_DA, randomGenerator);
                break;
            case 1:
                // random isotope and mass error
                randomShift = getRandomShiftIsotopes(baseMass, glycoIsotopes, glycoPPMtol, randomGenerator);
                break;
            case 2:
                // random mass error, no isotope error
                Integer[] noIsotopes = {0};
                randomShift = getRandomShiftIsotopes(baseMass, noIsotopes, glycoPPMtol, randomGenerator);
            case 3:
                // exact target mass - random shift left at 0
                break;
        }
        return randomShift;
    }

    /**
     * Generate a random mass shift within tolerancePPM about a randomly selected isotope peak in the
     * provided isotopes list.
     * @param isotopes list of isotopes
     * @param tolerancePPM Match tolerance (ppm) for glycan matching (from input parameter)
     * @param randomGenerator single random generator instance for whole glycan analysis
     * @return random mass shift within specified ranges
     */
    public static double getRandomShiftIsotopes(double glycanMass, Integer[] isotopes, double tolerancePPM, Random randomGenerator) {
        // randomly select isotope (must be sorted in ascending order)
        int minIso = isotopes[0];
        int maxIso = isotopes[isotopes.length - 1];
        // randomInt(0, max - min) + min yields correct range of min : max (including if min < 0)
        int isotope = randomGenerator.nextInt(maxIso + 1 - minIso) + minIso;  // upper bound is not inclusive, need to add 1 to get to max isotope

        // randomly generate mass shift within tolerance and add to chosen isotope
        double random = randomGenerator.nextDouble();       // between 0 and 1
        double baseMassEstimate = glycanMass + DEFAULT_PEPTIDE_MASS + isotope;
        double toleranceDa = baseMassEstimate * 1e-6 * tolerancePPM;
        double randomShift = -toleranceDa + random * (2 * toleranceDa);     // shift to range (min - random * (max - min)), where min = -toleranceDa and max = +toleranceDa
        return isotope * AAMasses.averagineIsotopeMass + randomShift;
    }

    /**
     * Generate a new glycan candidate database from the provided database of glycan fragment-specific propensities
     * and the input database for this glycoAnalysis
     * @param fragmentDB map of glycan string: fragment container
     * @param oldGlycoDB original glycan DB used for bootstrap analysis (or just initial DB provided by user)
     * @return glycan candidate arraylist
     */
    public ArrayList<GlycanCandidate> updateGlycanDatabase(HashMap<String, GlycanCandidateFragments> fragmentDB, ArrayList<GlycanCandidate> oldGlycoDB) {
        ArrayList<GlycanCandidate> newGlycoDB = new ArrayList<>();
        for (GlycanCandidate oldCandidate : oldGlycoDB) {
            GlycanCandidate newCandidate;
            String currentGlycanHash = Glycan.toGlycanString(oldCandidate.composition);
            if (oldCandidate.isDecoy) {
                // use target glycan propensity information for decoys as well
                currentGlycanHash = currentGlycanHash.replace("Decoy_", "");
            }

            // Do not include candidates without any fragment info (i.e., not passing filters after 1st pass)
            if (!fragmentDB.containsKey(currentGlycanHash)) {
                continue;
            }

            // initialize new candidate
            if (useShuffledIntensities) {
                // new method: shuffle fragment intensities rather than giving random masses. Init as target, then change
                newCandidate = GlycanCandidate.initGlycanCandidate(oldCandidate.composition,
                        0.0,
                        false,
                        this.glycanResiduesMap,
                        this.nGlycan,
                        this.randomGenerator,
                        this.glycoOxoniumDatabase);
            } else {
                // original method: shift decoy masses with same propensities/intensities as target
                newCandidate = GlycanCandidate.copyCandidate(oldCandidate, this.glycanResiduesMap);
            }

            GlycanCandidateFragments fragmtInfo = fragmentDB.getOrDefault(currentGlycanHash, new GlycanCandidateFragments());
            newCandidate.Yfragments = initFragmentsFromConsensus(newCandidate.Yfragments, fragmtInfo.yFragmentProps, fragmtInfo.yFragmentIntensities);
            newCandidate.oxoniumFragments = initFragmentsFromConsensus(newCandidate.oxoniumFragments, fragmtInfo.OxFragmentProps, fragmtInfo.OxFragmentIntensities);
            newCandidate.generalOxoniumFragments = initFragmentsFromConsensus(oldCandidate.generalOxoniumFragments, fragmtInfo.generalOxFragmentIntensities, fragmtInfo.generalOxFragmentIntensities);

            if (useShuffledIntensities) {
                if (oldCandidate.isDecoy) {
                    newCandidate.isDecoy = true;
                    newCandidate.mass = oldCandidate.mass;
                    newCandidate.decoyMassShift = oldCandidate.decoyMassShift;
                    shuffleFragmentIntensities(newCandidate.Yfragments);
                    shuffleFragmentIntensities(newCandidate.oxoniumFragments);
                    shuffleFragmentIntensities(newCandidate.generalOxoniumFragments);
                }
            }
            newGlycoDB.add(newCandidate);
        }
        return newGlycoDB;
    }

    private TreeMap<String, GlycanFragment> initFragmentsFromConsensus(TreeMap<String, GlycanFragment> originalFragments,
                                                                      HashMap<String, Double> fragmentPropensities,
                                                                      HashMap<String, Double> fragmentIntensities) {
        TreeMap<String, GlycanFragment> fragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> originalFragEntry : originalFragments.entrySet()) {
            String fragmentKey = originalFragEntry.getKey().replace("Decoy_", "");  // give decoys same fragment info as targets
            double expectedIntensity;
            double propensity;
            GlycanFragment origFrag = originalFragEntry.getValue();
            if (fragmentPropensities.containsKey(fragmentKey)) {
                // have propensity/intensity info for this fragment - read from input fragmentInfo
                expectedIntensity = fragmentIntensities.get(fragmentKey);
                propensity = fragmentPropensities.get(fragmentKey);
            } else {
                // no added info - copy the original
                expectedIntensity = origFrag.expectedIntensity;
                propensity = origFrag.propensity;
            }
            GlycanFragment newFragment = GlycanFragment.copyFragmentWithPropensity(origFrag, expectedIntensity, propensity);
            fragments.put(originalFragEntry.getKey(), newFragment);
        }
        return fragments;
    }

    private void shuffleFragmentIntensities(TreeMap<String, GlycanFragment> fragments) {
        // extract intensities
        ArrayList<Double> intensities = new ArrayList<>();
        for (GlycanFragment fragment : fragments.values()) {
            intensities.add(fragment.expectedIntensity);
        }
        // shuffle
        Collections.shuffle(intensities, randomGenerator);
        // reassign
        int index = 0;
        for (GlycanFragment fragment : fragments.values()) {
            fragment.expectedIntensity = intensities.get(index);
            index++;
        }
    }

    // Print glycan database (including decoys and associated mass shifts) to file
    public void printGlycanDatabase(String outputPath) {
        try {
            PrintWriter out = new PrintWriter(new FileWriter(outputPath));
            out.write("Composition\tMass\tIs Decoy\tDecoy Shift (Da)\tFragment Ions\n");
            for (GlycanCandidate candidate : glycoDatabase) {
                StringBuilder sb = new StringBuilder();
                sb.append(Glycan.toGlycanString(candidate.composition)).append("\t");
                sb.append(candidate.mass).append("\t");
                sb.append(candidate.isDecoy).append("\t");
                sb.append(candidate.decoyMassShift).append("\t");
                for (GlycanFragment ion : candidate.Yfragments.values()) {
                    sb.append(String.format("\tY~%s", ion));      // format is [ion type] [ion comp] [found intensity]
                }
                // oxonium ions
                for (GlycanFragment ion : candidate.oxoniumFragments.values()) {
                    sb.append(String.format("\tOx~%s", ion));      // format is [ion type] [ion comp] [found intensity]
                }
                sb.append("\n");
                out.write(sb.toString());
            }
            out.flush();
            out.close();
        } catch (IOException e) {
            PTMShepherd.die("Could not write glycan database to file " + outputPath + "\ndue to error: " + e.getMessage());
        }
    }


    /**
     * Parse isotopes parameter of format 'min,max' into Integer[] to pass to glyco analysis
     * @return Integer[] of all isotopes to consider (from min to max)
     */
    public static Integer[] parseGlycoIsotopesParam() {
        String isoLowStr = PTMShepherd.getParam("glyco_isotope_min");
        String isoHighStr = PTMShepherd.getParam("glyco_isotope_max");
        if (isoLowStr.length() > 0 && isoHighStr.length() > 0) {
            int minIso = Integer.parseInt(isoLowStr);
            int maxIso = Integer.parseInt(isoHighStr);
            // check for min/max swap from user input
            if (maxIso < minIso) {
                int saveMin = minIso;
                minIso = maxIso;
                maxIso = saveMin;
            }
            ArrayList<Integer> isotopes = new ArrayList<>();
            for (int i = minIso; i <= maxIso; i++) {
                isotopes.add(i);
            }
            return isotopes.toArray(new Integer[0]);

        } else {
            // todo: deprecated. Keeping for now for compatibility with old param files, but will remove at some point
            String paramValue = PTMShepherd.getParam("glyco_isotope_range");
            String[] paramStrs;
            ArrayList<Integer> isotopes = new ArrayList<Integer>();
            if (paramValue.length() > 0) {
                paramStrs = paramValue.split(",| |/");
                if (paramStrs.length == 2) {
                    int minIso = Integer.parseInt(paramStrs[0]);
                    int maxIso = Integer.parseInt(paramStrs[1]);
                    if (maxIso < minIso) {    // reverse if input flipped
                        int oldMax = maxIso;
                        maxIso = minIso;
                        minIso = oldMax;
                    }
                    for (int i = minIso; i <= maxIso; i++) {
                        isotopes.add(i);
                    }
                } else {
                    // invalid input: warn user
                    PTMShepherd.die(String.format("Invalid isotopes string %s input to glyco mode: must be in format 'min,max'", paramValue));
                }
            } else {
                return new Integer[]{-1, 0, 1, 2, 3};    // return default value if not specified
            }
            return isotopes.toArray(new Integer[0]);
        }
    }

    /**
     * Initialize a default glyco probability table and update it with probabilities from
     * parameters, if any are present.
     * Param formats:
     * -Fragments (Y/Oxo): 8 values, comma separated
     * -mass/isotope: key1:value1,key2:value2,etc
     * @return ProbabilityTable to use
     */
    public void initIsotopeProbs(String paramStr) {
        isotopeProbTable = new HashMap<>();
        if (!paramStr.matches("")) {
            String[] isoSplits = paramStr.split(",");
            for (String split : isoSplits) {
                String[] keyValue = split.trim().split(":");
                try {
                    isotopeProbTable.put(Integer.parseInt(keyValue[0]), Double.parseDouble(keyValue[1]));
                } catch (NumberFormatException ex) {
                    System.out.printf("Illegal character %s in parameter prob_isotopes, must pairs of numbers like '0:1.5,1:0.75'", split);
                }
            }
        }
        if (isotopeProbTable.isEmpty()) {
            // init defaults
            isotopeProbTable.put(-2, 0.125);
            isotopeProbTable.put(-1, 0.25);
            isotopeProbTable.put(0, 1.0);
            isotopeProbTable.put(1, 0.95);
            isotopeProbTable.put(2, 0.5);
            isotopeProbTable.put(3, 0.25);
            isotopeProbTable.put(4, 0.125);
        }
    }

    /**
     * Update isotope probabilities from the first pass of glycan assignment.
     * New probability is the proportion of the Results with this isotope, scaled so that the most common
     * isotope has a probability of 1.0.
     * @param isotopeCounts map of isotope: count
     */
    public void updateIsotopesProbsFromFirstPass(HashMap<Integer, Integer> isotopeCounts) {
        int maxCount = isotopeCounts.values().stream().max(Integer::compare).orElse(0);
        // set min prob to the lowest observed prob divided by 2
        double minProb = ((double) isotopeCounts.values().stream().min(Integer::compare).orElse(0) / maxCount) / 2.0;

        for (Map.Entry<Integer, Double> entry : isotopeProbTable.entrySet()) {
            int isotope = entry.getKey();
            if (isotopeCounts.containsKey(isotope)) {
                double newProb = (double) isotopeCounts.get(isotope) / maxCount;  // normalize to max count
                isotopeProbTable.put(isotope, newProb);
            } else {
                isotopeProbTable.put(isotope, minProb);
            }
        }
    }

    /**
     * Helper method to parse probability arrays from provided parameters with error checking. Returns
     * double[] of length 2 if successful or empty array if failed.
     * @param paramStr param to parse
     * @return double[] of length 2 if successful or empty array if failed.
     */
    public static double[] parseProbParam(String paramStr, String paramName) {
        String[] splits = paramStr.split(",");
        if (splits.length == 2 || splits.length == 3) {
            double[] values = new double[splits.length];
            for (int i = 0; i < splits.length; i++) {
                try {
                    values[i] = Double.parseDouble(splits[i]);
                } catch (NumberFormatException ex) {
                    System.out.printf("Invalid character in value %s, must be a number", splits[i]);
                    return new double[0];
                }
            }
            return values;
        } else {
            System.out.printf("Invalid format for parameter %s, must have 2 or 3 comma-separated values. Param was: %s\n", paramName, paramStr);
            return new double[0];
        }
    }

    /**
     * Helper method for determining column numbers in the .rawglyco file
     * @param headerSplits file header string, split on the appropriate delimiter
     * @param columnName name of the column to find
     * @return column index
     */
    public static int getHeaderColIndex(String[] headerSplits, String columnName) {
        for (int i = 0; i < headerSplits.length; i++) {
            if (headerSplits[i].trim().matches(columnName)) {
                return i;
            }
        }
        // column not found, return -1
        return -1;
    }

    /**
     * Print glyco params used
     */
    public void printGlycoParams() {
        PTMShepherd.print("Glycan Assignment params:");
        PTMShepherd.print(String.format("\tGlycan FDR: %.1f%%", glycoFDR * 100));
        PTMShepherd.print(String.format("\tMass error (ppm): %.1f", glycoPPMtol));
        PTMShepherd.print(String.format("\tIsotope errors: %s", Arrays.toString(glycoIsotopes)));
        PTMShepherd.print(String.format("\tGlycan Database size (including adducts): %d", glycoDatabase.size() / 2));
        if (nGlycan) {
            PTMShepherd.print("\tmode: N-glycan");
            PTMShepherd.print("\tAllowed Sites: N in N-X-S/T sequons only");
        } else {
            PTMShepherd.print(String.format("\tAllowed Sites: %s", allowedLocalizationResidues));
        }
        // print residue and mod lists
        if (printFullParams) {
            PTMShepherd.print("\tGlycan residue definitions:");
            for (GlycanResidue residue : glycanResidues) {
                PTMShepherd.print("\t\t" + residue.printToDatabaseFile());
            }
        }

        if (printFullParams) {
            PTMShepherd.print(String.format("\tMass prob (score scaling): %.1f", massProbScaling));
            StringBuilder isoString = new StringBuilder("\tIsotope probs:");
            for (Map.Entry<Integer, Double> isoEntry: isotopeProbTable.entrySet()) {
                isoString.append(String.format(" %d:%.1f", isoEntry.getKey(), isoEntry.getValue()));
            }
            PTMShepherd.print(isoString.toString());
            PTMShepherd.print(String.format("\tNormalize Y ion counts: %s", glycoYnorm));
            PTMShepherd.print(String.format("\tTypical mass error std devs (for absolute score): %.1f", absScoreErrorParam));
            PTMShepherd.print(String.format("\tDecoy type: %d", decoyType));
            if (printGlycoDecoys) {
                PTMShepherd.print("\tPrinting decoy glycans");
            }
            if (removeGlycanDeltaMass) {
                PTMShepherd.print("\tRemoving glycan delta mass from PSM table");
            }
        }
    }

    public String generateLDAheader() {
        StringBuilder sb = new StringBuilder();
        sb.append("\t");
        for (LDAFeature feature : ldaFeaturesToUse) {
            sb.append(feature.name()).append("\t");
        }
        return sb.toString();
    }

    public static ArrayList<LDAFeature> parseLDAfeatures(String featuresStr) {
        ArrayList<LDAFeature> features = new ArrayList<>();
        if (featuresStr.isEmpty()) {
            // use default (original) set of features
            PTMShepherd.print("No score features specified, using default set.");
            features.add(LDAFeature.yscore);
            features.add(LDAFeature.oxo);
            features.add(LDAFeature.mass);
            features.add(LDAFeature.iso);
            return features;
        }
        String[] featureSplits = featuresStr.split("[,\\s;/]+");
        for (String str : featureSplits) {
            try {
                LDAFeature f = LDAFeature.valueOf(str.trim());
                features.add(f);
            } catch (IllegalArgumentException e) {
                PTMShepherd.print("Invalid LDA feature: " + str + ", skipping.");
            }
        }
        return features;
    }

    public enum LDAFeature {
        kl,
        ms1,
        yprop,
        yscore,
        oxo,
        mass,
        iso,
        iso2nd,
        glycanfreq,
        ysim,
        oxsim,
    }
}
