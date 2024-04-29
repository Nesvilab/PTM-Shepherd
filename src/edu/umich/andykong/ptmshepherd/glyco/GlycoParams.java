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
import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanMod;
import umich.ms.glyco.GlycanParser;
import umich.ms.glyco.GlycanResidue;

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
    public boolean useNewFDR;
    public boolean useNonCompFDR ;
    public double defaultProp;
    public String allowedLocalizationResidues;
    public HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase;
    public HashMap<Integer, Double> isotopeProbTable;
    public double massProbScaling;

    private static final String defaultResiduePath = "glycan_residues.txt";
    private static final String defaultModsPath = "glycan_mods.txt";
    public static final String defaultOxoPath = "oxonium_ion_list.txt";

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

        glycoOxoniumDatabase = GlycoAnalysis.parseOxoniumDatabase(oxoniumListPath, this);
    }

//    private ArrayList<GlycanResidue> parseGlycoResiduesDB(String glycanResiduesPath, String defaultDefinitionsPath) {
//        ArrayList<GlycanResidue> residues = new ArrayList<>();
//        BufferedReader in;
//        try {
//            // no glycan database provided - fall back to default glycan list in PeakAnnotator
//            if (glycanResiduesPath.matches("")) {
//                in = new BufferedReader(new InputStreamReader(GlycoParams.class.getResourceAsStream(defaultDefinitionsPath)));
//            } else {
//                in = new BufferedReader(new FileReader(glycanResiduesPath));
//            }
//            String line;
//            while ((line = in.readLine()) != null) {
//                if (line.startsWith("#")) {
//                    continue;
//                }
//                GlycanResidue residue = new GlycanResidue(line, glycanResidueCounter);
//                glycanResidueCounter++;
//                residues.add(residue);
//            }
//        } catch (IOException ex) {
//            ex.printStackTrace();
//            PTMShepherd.die(String.format("Error parsing input monosaccharide table %s", glycanResiduesPath));
//        }
//        return residues;
//    }


    /**
     * Write a list of masses for all glycan candidates to pass to IonQuant. Format is one mass per line.
     * @param glycanDB list of glycan candidates (whole glycan database)
     * @param outputPath where to save the file
     */
    public static void writeGlycanMassList(ArrayList<GlycanCandidate> glycanDB, String outputPath) throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(outputPath));
        HashSet<Integer> writtenMasses = new HashSet<>();
        for (GlycanCandidate candidate : glycanDB) {
            if (candidate.isDecoy) {
                continue;       // do not write decoy masses to list - only target masses are reported in PSM table, even if decoy is assigned
            }
            int roundedMass = (int) Math.round(candidate.monoisotopicMass * 100);
            if (!writtenMasses.contains(roundedMass)) {
                out.write(String.format("%.4f\n",candidate.monoisotopicMass));
                writtenMasses.add(roundedMass);
            }
        }
        out.flush();
        out.close();
    }

    /**
     * Parse FragPipe glycan database string to list of glycan candidates
     * @param glycanDBString
     * @return
     */
    public ArrayList<GlycanCandidate> parseGlycanDatabaseString(String glycanDBString) {
        ArrayList<Glycan> glycans = GlycanParser.parseGlycanDatabaseString(glycanDBString, glycanResiduesMap);
        return convertGlycansToCandidates(glycans);
    }

    /**
     * Deprecated method, retained for backwards compatibility and command line access.
     * Parse input glycan database file. Formatting: 1 glycan per line, "Residue1-count_Residue2-count_...\n"
     * @param inputPath path to input file
     * @return list of glycans to consider
     */
    public ArrayList<GlycanCandidate> parseGlycanDatabaseFile(String inputPath) {
        ArrayList<Glycan> glycans = GlycanParser.loadGlycansFromText(inputPath, GlycanParser.detectDBtype(inputPath), glycanResiduesMap);
        return convertGlycansToCandidates(glycans);
    }

    /**
     * Generate PTM-S GlycanCandidates from Glycans, preventing duplicates and adding decoy Candidates.
     * @param glycans
     * @return
     */
    private ArrayList<GlycanCandidate> convertGlycansToCandidates(ArrayList<Glycan> glycans) {
        HashMap<String, Boolean> glycansInDB = new HashMap<>();
        ArrayList<GlycanCandidate> glycanDB = new ArrayList<>();
        for (Glycan glycan: glycans) {
            if (glycan.composition.isEmpty()) {
                continue;
            }
            // generate a new candidate from this composition and add to DB
            GlycanCandidate candidate = new GlycanCandidate(glycan.composition, false, this);
            String compositionHash = candidate.toString();
            // prevent addition of duplicates if user has them in database
            if (!glycansInDB.containsKey(compositionHash)) {
                glycanDB.add(candidate);
                glycansInDB.put(compositionHash, Boolean.TRUE);
                // also add a decoy for this composition
                GlycanCandidate decoy = new GlycanCandidate(glycan.composition, true, this);
                glycanDB.add(decoy);
            }
        }
        return glycanDB;
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
            String currentGlycanHash = oldCandidate.toString();
            if (oldCandidate.isDecoy) {
                // use target glycan propensity information for decoys as well
                currentGlycanHash = currentGlycanHash.replace("Decoy_", "");
            }

            // Get fragment info if present and initialize new candidate based on the old and fragment info (if present)
            GlycanCandidateFragments fragmtInfo = fragmentDB.getOrDefault(currentGlycanHash, new GlycanCandidateFragments());
            newCandidate = new GlycanCandidate(oldCandidate, fragmtInfo, this);
            newGlycoDB.add(newCandidate);
        }
        return newGlycoDB;
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
        if (isotopeProbTable.size() == 0) {
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
        PTMShepherd.print("Assigning glycans:");
    }
}
