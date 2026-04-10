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

import static umich.ms.glyco.Glycan.toGlycanString;

/**
 * Organization - putting a bunch of glyco-specific utilities here rather than cluttering the main PTM-S.java
 */
public class GlycoParams {

    public Random randomGenerator;
    public HashMap<String, GlycanResidue> glycanResiduesMap;
    public ArrayList<GlycanResidue> glycanResidues;
    public ArrayList<GlycanCandidate> glycoDatabase;
    public ArrayList<GlycanCandidate> sortedGlycansByMass;
    public int decoyType;
    public double glycoPPMtol;
    public static final double DECOY_PPM_SHIFT = 30;
    public Integer[] glycoIsotopes;
    public static final Integer[] DECOY_ISOTOPES = {-1, 0, 1, 2, 3};
    public boolean nGlycan;
    public double glycoFDR;
    public boolean printFullParams;
    public boolean writeGlycansToAssignedMods;
    public boolean removeGlycanDeltaMass;
    public boolean printGlycoDecoys;
    public int numThreads;
    public String allowedLocalizationResidues;
    public HashMap<GlycanResidue, ArrayList<GlycanFragment>> glycoOxoniumDatabase;
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
    public boolean twoPassMode;
    public boolean removeGlycans2ndPass;
    public boolean cosineSimilarityScoring;
    public boolean glycoAvgInts;
    public int numDecoysPerTarget;
    public boolean checkVariableMods;
    public boolean noFDR;
    public double minDecoyFragmentDiff;
    public String glycoLibPath;
    public HashMap<String, GlycanCandidateFragments> glycoLibFragments;
    public HashMap<String, Integer> glycoLibCounts;
    public boolean useGlycoLibFirstPass;
    public boolean glycoSkipPairwise;
    public boolean updateGlycoLib;
    public boolean incr1stFDR;

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
     * Parse glyco library (.glycolib) file. Glycolib file format:
     *  ion type~composition [tab] expected intensity [tab] count in library
     *  Lines for each glycan are grouped together, starting with a line with the glycan composition alone (starting with "GLYCAN"), followed by lines for each fragment ion and associated info. Example:
     *   GLYCAN        HexNAc(4)Hex(5)NeuAc(1)
     *   Y~HexNAc(1)   0.999900        0.000141        2
     *   Y~HexNAc(2)   0.529000        0.008556        2
     *   ...
     *   Ox~NeuAc(1)   0.034600        0.008202        2
     *   END
     * Populates glycoLibFragments (keyed by glycan composition string) and glycoLibCounts (count of PSMs per glycan).
     */
    public void parseGlycoLib() {
        glycoLibFragments = new HashMap<>();
        glycoLibCounts = new HashMap<>();
        if (glycoLibPath == null || glycoLibPath.isEmpty()) {
            return;
        }
        File libFile = new File(glycoLibPath);
        if (!libFile.exists()) {
            PTMShepherd.print("\tGlyco library file not found at " + glycoLibPath + "; will create new library if glyco_update_lib is enabled");
            return;
        }
        HashSet<String> diagnosticIonHashes = buildDiagnosticIonSet();
        try (BufferedReader reader = new BufferedReader(new FileReader(glycoLibPath))) {
            String line;
            String currentGlycan = null;
            LinkedHashMap<String, Double> yFragments = null;
            LinkedHashMap<String, Double> oxFragments = null;
            LinkedHashMap<String, Double> generalOxFragments = null;
            int currentCount = 0;

            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                String[] parts = line.split("\t");
                if (parts[0].equals("GLYCAN")) {
                    currentGlycan = parts.length > 1 ? parts[1].trim() : null;
                    yFragments = new LinkedHashMap<>();
                    oxFragments = new LinkedHashMap<>();
                    generalOxFragments = new LinkedHashMap<>();
                    currentCount = 0;
                } else if (parts[0].equals("END")) {
                    if (currentGlycan != null) {
                        glycoLibFragments.put(currentGlycan, new GlycanCandidateFragments(yFragments, oxFragments, generalOxFragments));
                        glycoLibCounts.put(currentGlycan, currentCount);
                    }
                    currentGlycan = null;
                } else if (currentGlycan != null && parts.length >= 4) {
                    String ionKey = parts[0].trim();
                    double intensity = Double.parseDouble(parts[1].trim());
                    int count = Integer.parseInt(parts[3].trim());
                    if (count > currentCount) {
                        currentCount = count;
                    }
                    if (ionKey.startsWith("Y~")) {
                        yFragments.put(ionKey.substring(2), intensity);
                    } else if (ionKey.startsWith("Ox~")) {
                        String fragKey = ionKey.substring(3);
                        // All Ox~ ions are general oxonium ions
                        generalOxFragments.put(fragKey, intensity);
                        // Also add to diagnostic oxonium ions if this key is in the diagnostic set
                        if (diagnosticIonHashes.contains(fragKey)) {
                            oxFragments.put(fragKey, intensity);
                        }
                    }
                }
            }
        } catch (IOException e) {
            PTMShepherd.die("Error reading glyco library file " + glycoLibPath + ": " + e.getMessage());
        }
        PTMShepherd.print(String.format("\tLoaded glyco library with %d glycan entries", glycoLibFragments.size()));
    }

    /**
     * Build a set of fragment hashes for all diagnostic oxonium ions in the glycoOxoniumDatabase.
     * Uses the fragment's hash field (which encodes composition and any composition comment).
     * Only non-decoy fragments are included since glycolib keys are target-oriented.
     * @return set of hash strings for diagnostic oxonium ions
     */
    public HashSet<String> buildDiagnosticIonSet() {
        HashSet<String> diagnosticHashes = new HashSet<>();
        for (ArrayList<GlycanFragment> fragments : glycoOxoniumDatabase.values()) {
            for (GlycanFragment fragment : fragments) {
                if (fragment.isDiagnostic && !fragment.isDecoy) {
                    diagnosticHashes.add(fragment.hash);
                }
            }
        }
        return diagnosticHashes;
    }

    /**
     * Find the GlycanCandidateFragments for the given glycan composition from the glycolib.
     * If an exact match is not found, returns the library entry with the smallest number of total
     * residue differences to the query glycan. On ties, returns the entry with the largest PSM count.
     * @param glycanString glycan composition string (e.g., "HexNAc(4)Hex(5)NeuAc(1)")
     * @return GlycanCandidateFragments for the nearest library entry, or empty if library is empty
     */
    public GlycanCandidateFragments findNearestGlycanInLib(String glycanString) {
        if (glycoLibFragments == null || glycoLibFragments.isEmpty()) {
            return new GlycanCandidateFragments();
        }
        if (glycoLibFragments.containsKey(glycanString)) {
            return glycoLibFragments.get(glycanString);
        }
        HashMap<String, Integer> queryComp = parseGlycanStringToResidueMap(glycanString);
        String nearestKey = null;
        int minDistance = Integer.MAX_VALUE;
        int maxCount = -1;
        for (String libKey : glycoLibFragments.keySet()) {
            HashMap<String, Integer> libComp = parseGlycanStringToResidueMap(libKey);
            int distance = computeCompositionDistance(queryComp, libComp);
            int count = glycoLibCounts.getOrDefault(libKey, 0);
            if (distance < minDistance || (distance == minDistance && count > maxCount)) {
                minDistance = distance;
                maxCount = count;
                nearestKey = libKey;
            }
        }
        return nearestKey != null ? glycoLibFragments.get(nearestKey) : new GlycanCandidateFragments();
    }

    /**
     * Parse a glycan composition string (e.g., "HexNAc(4)Hex(5)NeuAc(1)") into a map of residue name to count.
     */
    private static HashMap<String, Integer> parseGlycanStringToResidueMap(String glycanString) {
        HashMap<String, Integer> compositionMap = new HashMap<>();
        java.util.regex.Pattern pattern = java.util.regex.Pattern.compile("(\\w+)\\((\\d+)\\)");
        java.util.regex.Matcher matcher = pattern.matcher(glycanString);
        while (matcher.find()) {
            compositionMap.put(matcher.group(1), Integer.parseInt(matcher.group(2)));
        }
        return compositionMap;
    }

    /**
     * Compute the total number of residue differences between two glycan compositions.
     * Counts each residue type's count difference, including residues present in one but not the other.
     */
    private static int computeCompositionDistance(HashMap<String, Integer> comp1, HashMap<String, Integer> comp2) {
        int distance = 0;
        HashSet<String> allKeys = new HashSet<>(comp1.keySet());
        allKeys.addAll(comp2.keySet());
        for (String key : allKeys) {
            distance += Math.abs(comp1.getOrDefault(key, 0) - comp2.getOrDefault(key, 0));
        }
        return distance;
    }

    /**
     * Parse FragPipe glycan database string to list of glycan candidates
     * @param glycanDBString
     * @return
     */
    public ArrayList<GlycanCandidate> parseGlycanDatabaseString(String glycanDBString) {
        ArrayList<Glycan> glycans = GlycanParser.parseGlycanDatabaseString(glycanDBString, glycanResiduesMap);
        return convertGlycansToCandidates(glycans, glycanResiduesMap, nGlycan, glycoOxoniumDatabase, decoyType, DECOY_PPM_SHIFT, DECOY_ISOTOPES, randomGenerator, numDecoysPerTarget, useGlycoLibFirstPass);
    }

    /**
     * Deprecated method, retained for backwards compatibility and command line access.
     * Parse input glycan database file. Formatting: 1 glycan per line, "Residue1-count_Residue2-count_...\n"
     * @param inputPath path to input file
     * @return list of glycans to consider
     */
    public ArrayList<GlycanCandidate> parseGlycanDatabaseFile(String inputPath) {
        ArrayList<Glycan> glycans = GlycanParser.loadGlycansFromText(inputPath, GlycanParser.detectDBtype(inputPath), glycanResiduesMap);
        return convertGlycansToCandidates(glycans, glycanResiduesMap, nGlycan, glycoOxoniumDatabase, decoyType, DECOY_PPM_SHIFT, DECOY_ISOTOPES, randomGenerator, numDecoysPerTarget, useGlycoLibFirstPass);
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
                                                                        Random randomGenerator,
                                                                        int numDecoysPerTarget,
                                                                        boolean useGlycoLib) {
        LinkedHashMap<String, Boolean> glycansInDB = new LinkedHashMap<>();
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
                // add numDecoysPerTarget decoys for this composition
                for (int d = 0; d < numDecoysPerTarget; d++) {
                    double decoyMassShift = setDecoyShift(candidate.mass, decoyType, glycoPPMtol, glycoIsotopes, randomGenerator);
                    GlycanCandidate decoy;
                    if (!useGlycoLib) {
                        // original method: shift fragment masses (up to 20)
                        decoy = GlycanCandidate.initGlycanCandidate(glycan.composition, decoyMassShift, true, glycanResiduesMap, nGlycan, randomGenerator, glycoOxoniumDatabase);
                    } else {
                        // glycolib method: alter fragment intensities rather than giving random masses. Init as target, then change to decoy later so avoid fragment mass shifting
                        decoy = GlycanCandidate.initGlycanCandidate(candidate.composition,
                                0.0,
                                false,
                                glycanResiduesMap,
                                nGlycan,
                                randomGenerator,
                                glycoOxoniumDatabase);
                        decoy.isDecoy = true;
                        decoy.mass = candidate.mass + decoyMassShift;
                        decoy.decoyMassShift = decoyMassShift;
                    }
                    glycanDB.add(decoy);
                }
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
     * Initialize sorted glycan database by mass for efficient nearest mass lookups.
     * Call this after glycoDatabase is populated.
     */
    public void initSortedGlycanDatabase(ArrayList<GlycanCandidate> glycoDatabase) {
        sortedGlycansByMass = new ArrayList<>(glycoDatabase);
        sortedGlycansByMass.sort(Comparator.comparingDouble(g -> g.mass));
    }

    /**
     * Generate a new glycan candidate database from the provided database of glycan fragment-specific propensities
     * and the input database for this glycoAnalysis
     * @param fragmentDB map of glycan string: fragment container
     * @param oldGlycoDB original glycan DB used for bootstrap analysis (or just initial DB provided by user)
     * @return glycan candidate arraylist
     */
    public ArrayList<GlycanCandidate> updateGlycanDatabase(HashMap<String, GlycanCandidateFragments> fragmentDB, LinkedHashMap<String, GlycanCandidateFragments> decoyFragmentDB, ArrayList<GlycanCandidate> oldGlycoDB) {
        ArrayList<GlycanCandidate> newGlycoDB = new ArrayList<>();
        if (removeGlycans2ndPass) {
            ArrayList<GlycanCandidate> reducedDB = new ArrayList<>();
            for (GlycanCandidate candidate : oldGlycoDB) {
                String glycanHash = toGlycanString(candidate.composition);
                if (candidate.isDecoy) {
                    glycanHash = glycanHash.replace("Decoy_", "");
                }
                if (fragmentDB.containsKey(glycanHash)) {
                    reducedDB.add(candidate);
                }
            }
            oldGlycoDB = reducedDB;
        }

        // update the glycan candidates to generate the new database
        for (GlycanCandidate oldCandidate : oldGlycoDB) {
            GlycanCandidate newCandidate;
            String currentGlycanHash = toGlycanString(oldCandidate.composition);
            if (oldCandidate.isDecoy) {
                // use target glycan propensity information for decoys as well
                currentGlycanHash = currentGlycanHash.replace("Decoy_", "");
            }

            // init new candidate: alter fragment intensities rather than giving random masses. Init as target, then change to decoy later so avoid fragment mass shifting
            newCandidate = GlycanCandidate.initGlycanCandidate(oldCandidate.composition,
                    0.0,
                    false,
                    this.glycanResiduesMap,
                    this.nGlycan,
                    this.randomGenerator,
                    this.glycoOxoniumDatabase);


            GlycanCandidateFragments fragmtInfo;
            if (oldCandidate.isDecoy) {
                fragmtInfo = decoyFragmentDB.getOrDefault(currentGlycanHash, new GlycanCandidateFragments());
                newCandidate.isDecoy = true;
                // todo: change to mass of the alternate matched glycan?
                newCandidate.mass = oldCandidate.mass;
                newCandidate.decoyMassShift = oldCandidate.decoyMassShift;
            } else {
                fragmtInfo = fragmentDB.getOrDefault(currentGlycanHash, new GlycanCandidateFragments());
            }
            newCandidate.Yfragments = initFragmentsFromConsensus(newCandidate.Yfragments, fragmtInfo.yFragmentIntensities);
            newCandidate.oxoniumFragments = initFragmentsFromConsensus(newCandidate.oxoniumFragments, fragmtInfo.OxFragmentIntensities);
            newCandidate.generalOxoniumFragments = initFragmentsFromConsensus(oldCandidate.generalOxoniumFragments, fragmtInfo.generalOxFragmentIntensities);

            newGlycoDB.add(newCandidate);
        }
        return newGlycoDB;
    }

    private TreeMap<String, GlycanFragment> initFragmentsFromConsensus(TreeMap<String, GlycanFragment> originalFragments,
                                                                      LinkedHashMap<String, Double> fragmentIntensities) {
        TreeMap<String, GlycanFragment> fragments = new TreeMap<>();
        for (Map.Entry<String, GlycanFragment> originalFragEntry : originalFragments.entrySet()) {
            String fragmentKey = originalFragEntry.getKey().replace("Decoy_", "");  // give decoys same fragment info as targets
            GlycanFragment origFrag = originalFragEntry.getValue();
            double expectedIntensity = fragmentIntensities.getOrDefault(fragmentKey, 0.0);
            GlycanFragment newFragment = GlycanFragment.copyFragmentWithPropensity(origFrag, expectedIntensity, 0);
            fragments.put(originalFragEntry.getKey(), newFragment);
        }
        return fragments;
    }

    // Print glycan database (including decoys and associated mass shifts) to file
    public void printGlycanDatabase(String outputPath) {
        try {
            PrintWriter out = new PrintWriter(new FileWriter(outputPath));
            out.write("Composition\tMass\tIs Decoy\tDecoy Shift (Da)\tFragment Ions\n");
            for (GlycanCandidate candidate : glycoDatabase) {
                StringBuilder sb = new StringBuilder();
                sb.append(toGlycanString(candidate.composition)).append("\t");
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
     * Print glyco params used
     */
    public void printGlycoParams() {
        PTMShepherd.print("Glycan Assignment params:");
        PTMShepherd.print(String.format("\tGlycan FDR: %.1f%%", glycoFDR * 100));
        PTMShepherd.print(String.format("\tMass error (ppm): %.1f", glycoPPMtol));
        PTMShepherd.print(String.format("\tIsotope errors: %s", Arrays.toString(glycoIsotopes)));
        PTMShepherd.print(String.format("\tGlycan Database size (including adducts): %d", glycoDatabase.size() / (1 + numDecoysPerTarget)));
        PTMShepherd.print(String.format("\tDecoys per target glycan: %d", numDecoysPerTarget));
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
            // todo: add new params
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
        sb.append("Y score\tOx diagnostic score\tOx similarity score\tMass score\t");      // always used
        for (LDAFeature feature : ldaFeaturesToUse) {
            sb.append(feature.name()).append("\t");
        }
        return sb.toString();
    }

    public static ArrayList<LDAFeature> parseLDAfeatures(String featuresStr) {
        ArrayList<LDAFeature> features = new ArrayList<>();
        if (featuresStr.isEmpty()) {
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
        ms1delta,
        glycanfreq,
    }

    /**
     * Update the in-memory glycolib with new fragment intensity results from a completed glyco analysis.
     * New glycans not already in the library are added directly. Existing glycans receive a weighted merge:
     * for each fragment, if the new median intensity is > 0, the updated intensity is the weighted average
     * of the existing value (weighted by existing library count) and the new value (weighted by new count).
     * The library PSM count is incremented by the new count.
     * @param newFragments map of glycan composition string → fragment intensities from new results
     * @param newCounts map of glycan composition string → number of PSMs in new results
     */
    public void updateGlycoLibInMemory(LinkedHashMap<String, GlycanCandidateFragments> newFragments, HashMap<String, Integer> newCounts) {
        if (glycoLibFragments == null) glycoLibFragments = new HashMap<>();
        if (glycoLibCounts == null) glycoLibCounts = new HashMap<>();
        for (Map.Entry<String, GlycanCandidateFragments> entry : newFragments.entrySet()) {
            String glycanKey = entry.getKey();
            GlycanCandidateFragments newFragData = entry.getValue();
            int newCount = newCounts.getOrDefault(glycanKey, 0);
            if (newCount == 0) {
                continue;
            }

            if (!glycoLibFragments.containsKey(glycanKey)) {
                glycoLibFragments.put(glycanKey, newFragData);
                glycoLibCounts.put(glycanKey, newCount);
            } else {
                GlycanCandidateFragments existing = glycoLibFragments.get(glycanKey);
                int existingCount = glycoLibCounts.getOrDefault(glycanKey, 0);
                glycoLibFragments.put(glycanKey, weightedMergeFragments(existing, existingCount, newFragData, newCount));
                glycoLibCounts.put(glycanKey, existingCount + newCount);
            }
        }
    }

    private GlycanCandidateFragments weightedMergeFragments(GlycanCandidateFragments existing, int existingCount,
                                                             GlycanCandidateFragments newFrags, int newCount) {
        return new GlycanCandidateFragments(
                weightedMergeFragmentMap(existing.yFragmentIntensities, existingCount, newFrags.yFragmentIntensities, newCount),
                weightedMergeFragmentMap(existing.OxFragmentIntensities, existingCount, newFrags.OxFragmentIntensities, newCount),
                weightedMergeFragmentMap(existing.generalOxFragmentIntensities, existingCount, newFrags.generalOxFragmentIntensities, newCount)
        );
    }

    private LinkedHashMap<String, Double> weightedMergeFragmentMap(LinkedHashMap<String, Double> existing, int existingCount,
                                                                     LinkedHashMap<String, Double> newMap, int newCount) {
        LinkedHashMap<String, Double> result = new LinkedHashMap<>(existing);
        for (Map.Entry<String, Double> entry : newMap.entrySet()) {
            String key = entry.getKey();
            double newIntensity = entry.getValue();
            if (newIntensity <= 0.0) {
                continue;  // only update if new median intensity > 0
            }

            if (result.containsKey(key)) {
                double existingIntensity = result.get(key);
                result.put(key, (existingIntensity * existingCount + newIntensity * newCount) / (double) (existingCount + newCount));
            } else {
                result.put(key, newIntensity);
            }
        }
        return result;
    }

    /**
     * Write the current in-memory glycolib (glycoLibFragments and glycoLibCounts) to the glycolib file at glycoLibPath.
     * Format: tab-delimited, one glycan entry per GLYCAN...END block.
     */
    public void writeGlycoLib() {
        if (glycoLibPath == null || glycoLibPath.isEmpty()) {
            PTMShepherd.print("Warning: glyco_update_lib is true but no glyco_lib_path is set; skipping glycolib update");
            return;
        }
        try (PrintWriter out = new PrintWriter(new FileWriter(glycoLibPath))) {
            for (String glycanKey : glycoLibFragments.keySet()) {
                GlycanCandidateFragments frags = glycoLibFragments.get(glycanKey);
                int count = glycoLibCounts.getOrDefault(glycanKey, 0);
                out.write("GLYCAN\t" + glycanKey + "\n");
                for (Map.Entry<String, Double> fragEntry : frags.yFragmentIntensities.entrySet()) {
                    out.write(String.format("Y~%s\t%.6f\t0.000000\t%d\n", fragEntry.getKey(), fragEntry.getValue(), count));
                }
                for (Map.Entry<String, Double> fragEntry : frags.generalOxFragmentIntensities.entrySet()) {
                    out.write(String.format("Ox~%s\t%.6f\t0.000000\t%d\n", fragEntry.getKey(), fragEntry.getValue(), count));
                }
                out.write("END\n");
            }
        } catch (IOException e) {
            PTMShepherd.die("Error writing updated glyco library to " + glycoLibPath + ": " + e.getMessage());
        }
        PTMShepherd.print(String.format("\tUpdated glyco library written to %s with %d glycan entries", glycoLibPath, glycoLibFragments.size()));
    }
}
