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
import edu.umich.andykong.ptmshepherd.peakpicker.PeakAnnotator;
import org.hipparchus.random.RandomGenerator;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Organization - putting a bunch of glyco-specific utilities here rather than cluttering the main PTM-S.java
 */
public class GlycoParams {

    public Random randomGenerator;
    public ArrayList<GlycanMod> glycanMods;
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
    public ProbabilityTables glycoProbabilityTable;     // todo: remove?

    private static final String defaultResiduePath = "glycan_residues.tsv";
    private static final String defaultModsPath = "glycan_mods.tsv";

    private static final Pattern numberPattern = Pattern.compile("[0-9]");
    private static final Pattern letterPattern = Pattern.compile("[a-zA-Z]+");
    private static final Pattern pGlycoPattern = Pattern.compile("[()]?[AGFHNP]+[()]");         // e.g., (N(H))
    private static final Pattern byonicPattern = Pattern.compile("[A-Za-z]+\\([0-9]+\\)");      // e.g., HexNAc(1)Hex(1)
    private static final Pattern metamorpheusPattern = Pattern.compile("[A-Z]+[0-9]+");         // e.g., N1H1
    private static final Pattern oldPTMSPattern = Pattern.compile("[A-Za-z]+-[0-9]+");          // e.g., HexNAc-1_Hex-1

    public GlycoParams(String glycanResiduesPath, String glycanModsPath) {
        // parse the glycan residues and mods tables, using internal defaults if no paths provided from FragPipe or user
        glycanResidues = parseGlycoResiduesDB(glycanResiduesPath);
        glycanMods = parseGlycoModsDB(glycanModsPath);
    }

    private ArrayList<GlycanResidue> parseGlycoResiduesDB(String glycanResiduesPath) {
        ArrayList<GlycanResidue> residues = new ArrayList<>();
        BufferedReader in;
        try {
            // no glycan database provided - fall back to default glycan list in PeakAnnotator
            if (glycanResiduesPath.matches("")) {
                in = new BufferedReader(new InputStreamReader(GlycoParams.class.getResourceAsStream(glycanResiduesPath)));
            } else {
                in = new BufferedReader(new FileReader(glycanResiduesPath));
            }
            String line;
            while ((line = in.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue;
                }
                GlycanResidue residue = GlycanResidue.parseResidue(line);
                residues.add(residue);
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            PTMShepherd.die(String.format("Error parsing input monosaccharide table %s", glycanResiduesPath));
        }
        return residues;
    }

    private ArrayList<GlycanMod> parseGlycoModsDB(String glycanModsPath) {
        ArrayList<GlycanMod> mods = new ArrayList<>();
        BufferedReader in;
        try {
            // no glycan database provided - fall back to default glycan list in PeakAnnotator
            if (glycanModsPath.matches("")) {
                in = new BufferedReader(new InputStreamReader(GlycoParams.class.getResourceAsStream(glycanModsPath)));
            } else {
                in = new BufferedReader(new FileReader(glycanModsPath));
            }
            String line;
            while ((line = in.readLine()) != null) {
                if (line.startsWith("#")) {
                    continue;
                }
//                GlycanMod mod = GlycanMod.parseMod(line);
//                mods.add(mod);
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            PTMShepherd.die(String.format("Error parsing input monosaccharide table %s", glycanModsPath));
        }
        return mods;
    }


    /**
     * Write a list of masses for all glycan candidates to pass to IonQuant. Format is one mass per line.
     * @param glycanDB list of glycan candidates (whole glycan database)
     * @param outputPath where to save the file
     */
    public static void writeGlycanMassList(ArrayList<GlycanCandidate> glycanDB, String outputPath) throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(outputPath));
        for (GlycanCandidate candidate : glycanDB) {
            if (candidate.isDecoy) {
                continue;       // do not write decoy masses to list - only target masses are reported in PSM table, even if decoy is assigned
            }
            out.write(String.format("%.4f\n",candidate.monoisotopicMass));
        }
        out.flush();
        out.close();
    }

    /**
     * Parse input glycan database file. Formatting: 1 glycan per line, "Residue1-count_Residue2-count_...\n"
     * @param inputPath path to input file
     * @return list of glycans to consider
     */
    public ArrayList<GlycanCandidate> parseGlycanDatabase(String inputPath) {
        // read input glycan database or default database if none provided
        ArrayList<GlycanCandidate> glycanDB = new ArrayList<>();
        DatabaseType dbType;
        try {
            BufferedReader in;
            if (inputPath.equals("")) {
                // no glycan database provided - fall back to default glycan list
                String defaultDB;
                if (nGlycan) {
                    defaultDB = "glyco_mods_20210127.txt";
                } else {
                    defaultDB = "glyco_mods_N-O_20211025.txt";
                }
                in = new BufferedReader(new InputStreamReader(PeakAnnotator.class.getResourceAsStream(defaultDB)));
                dbType = DatabaseType.byonic;   // internal DBs are in Byonic format
            } else {
                Path path = null;
                try {
                    path = Paths.get(inputPath.replaceAll("['\"]", ""));
                } catch (Exception e) {
                    System.out.println(e);
                    PTMShepherd.die(String.format("Malformed glycan path string: [%s]", inputPath));
                }
                if (path == null || !Files.exists(path)) {
                    PTMShepherd.die(String.format("Glycan database file does not exist: [%s]", inputPath));
                }
                in = new BufferedReader(new FileReader(path.toFile()));
                dbType = detectDBtype(path);
                if (dbType == DatabaseType.unknown) {
                    PTMShepherd.die(String.format("Error: could not detect glycan database type for file %s. Please check the formatting and file path", path));
                }
            }

            HashMap<String, Boolean> glycansInDB = new HashMap<>();
            // parse decoy param

            String readline;
            while ((readline = in.readLine()) != null) {
                String[] lineSplits;
                if (!(dbType == DatabaseType.pGlyco)) {
                    // handle csv and tsv inputs that may contain multiple entries per line
                    if (readline.contains(",")) {
                        lineSplits = readline.split(",");
                    } else if (readline.contains("\t") && readline.contains("%")) {
                        lineSplits = readline.split("\t");
                    } else {
                        lineSplits = new String[]{readline};
                    }
                } else {
                    // single line only for pGlyco
                    lineSplits = new String[] {readline};
                }

                // parse glycans
                for (String line : lineSplits) {
                    TreeMap<GlycanResidue, Integer> glycanComp = new TreeMap<>();
                    switch (dbType){
                        case byonic:
                            glycanComp = parseGlycan(line, byonicPattern.matcher(line));
                            break;
                        case pGlyco:
                            glycanComp = parsePGlycoGlycan(line);
                            break;
                        case metamorpheus:
                            glycanComp = parseGlycan(line, metamorpheusPattern.matcher(line));
                            break;
                        case oldPTMShepherd:
                            glycanComp = parseGlycan(line, oldPTMSPattern.matcher(line));
                            break;
                    }

                    if (glycanComp.size() == 0) {
                        continue;
                    }

                    // generate a new candidate from this composition and add to DB
                    GlycanCandidate candidate = new GlycanCandidate(glycanComp, false, this);
                    String compositionHash = candidate.toString();
                    // prevent addition of duplicates if user has them in database
                    if (!glycansInDB.containsKey(compositionHash)) {
                        glycanDB.add(candidate);
                        glycansInDB.put(compositionHash, Boolean.TRUE);
                        // also add a decoy for this composition
                        GlycanCandidate decoy = new GlycanCandidate(glycanComp, true, this);
                        glycanDB.add(decoy);

                        // todo: add mods/adducts handling
                    }
                }
            }

        } catch (FileNotFoundException e) {
            PTMShepherd.die(String.format("Glycan database file not found: [%s]", inputPath));
        } catch (IOException e) {
            e.printStackTrace();
            PTMShepherd.die("IO Exception while reading database file");
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
     * Parse glycans of varying formats given an input Regex Pattern.Matcher
     * @param glycanFormatRegex from Pattern.Matcher(string)
     * @return
     */
    private TreeMap<GlycanResidue, Integer> parseGlycan(String line, Matcher glycanFormatRegex) {
        TreeMap<GlycanResidue, Integer> glycanComp = new TreeMap<>();
        while(glycanFormatRegex.find()) {
            String glycanToken = glycanFormatRegex.group();
            Matcher glycanName = letterPattern.matcher(glycanToken);
            GlycanResidue residue;
            if (glycanName.find()) {
                String glycan = glycanName.group();
                residue = findResidueName(glycan);
                if (residue == null) {
                    // todo: try to match to a modification
                    PTMShepherd.print(String.format("Warning: glycan type %s in glycan %s not recognized. This glycan will not be included in search! Please add %s to the glycan_residues.tsv table.", glycanName, line, glycanName));
                    return new TreeMap<>();
                }
            } else {
                PTMShepherd.print(String.format("Warning: glycan type %s in glycan %s not recognized. This glycan will not be included in search! Please add %s to the glycan_residues.tsv table.", glycanName, line, glycanName));
                return new TreeMap<>();
            }
            Matcher glycanCount = numberPattern.matcher(glycanToken);
            if (glycanCount.find()) {
                int count = Integer.parseInt(glycanCount.group());
                glycanComp.put(residue, count);
            } else {
                PTMShepherd.print(String.format("Warning: glycan count not found in glycan %s and was ignored", line));
            }
        }
        return glycanComp;
    }

    /**
     * Search all GlycanResidue names and alternate names to try to match this input to one of them.
     * @param inputName string
     * @return
     */
    public GlycanResidue findResidueName(String inputName) {
        for (GlycanResidue residue: glycanResidues) {
            if (residue.name.equalsIgnoreCase(inputName)) {
                return residue;
            }
        }
        // no exact match: check alternate names
        for (GlycanResidue residue: glycanResidues) {
            for (String altName : residue.alternateNames) {
                if (altName.equalsIgnoreCase(inputName)) {
                    return residue;
                }
            }
        }
        // no match - return null and warn user
        return null;
    }

    /**
     * Parse pGlyco glycan from String
     * @param line
     * @return
     */
    private TreeMap<GlycanResidue, Integer> parsePGlycoGlycan(String line) {
        TreeMap<GlycanResidue, Integer> glycanComp = new TreeMap<>();

        // skip empty lines and lines without a glycan structure
        if (!line.contains("(")) {
            return glycanComp;
        }

        // parse glycans from tokens
        Matcher matcher = pGlycoPattern.matcher(line);
        while(matcher.find()) {
            String glycanToken = matcher.group();
            GlycanResidue residue = findResidueName(glycanToken);
            if (residue == null) {
                // todo: try to match to a modification
                PTMShepherd.print(String.format("Warning: glycan type %s in glycan %s not recognized. This glycan will not be included in search! Please add %s to the glycan_residues.tsv table.", glycanToken, line, glycanToken));
                return new TreeMap<>();
            }
            if (glycanComp.containsKey(residue)) {
                glycanComp.put(residue, glycanComp.get(residue) + 1);
            } else {
                glycanComp.put(residue, 1);
            }
        }
        return glycanComp;
    }

    /**
     * Peek at first line(s) of the database file and determine what type it is for parsing
     * @param databasePath reader
     * @return db type
     */
    public static DatabaseType detectDBtype(Path databasePath) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(databasePath.toFile()));
        String line;
        int i = 0;
        while ((line = in.readLine()) != null) {
            if (i == 0) {
                // headers are not required, but may be present
                if (line.contains("#") || line.contains("\\\\") || line.contains(",") || line.startsWith("%")) {
                    i++;
                    continue;
                }
            }
            if (line.contains("%")) {
                line = line.split("%")[0];
            }

            // check if line matches any of the known patterns
            Matcher byonicMatcher = byonicPattern.matcher(line);
            Matcher metamMatcher = metamorpheusPattern.matcher(line);
            Matcher pGlycoMatcher = pGlycoPattern.matcher(line);
            Matcher oldPTMSMatcher = oldPTMSPattern.matcher(line);
            if (byonicMatcher.find()) {
                return DatabaseType.byonic;
            } else if (oldPTMSMatcher.find()) {
                return DatabaseType.oldPTMShepherd;
            } else if (metamMatcher.find()) {
                return DatabaseType.metamorpheus;
            } else if (pGlycoMatcher.find()) {
                return DatabaseType.pGlyco;
            }
            i++;
        }
        return DatabaseType.unknown;
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
    public static ProbabilityTables initGlycoProbTable(GlycoParams glycoParams) {
        ProbabilityTables probabilityTable = new ProbabilityTables();

        // Read params for probabilities if present, otherwise use default values (already init'd in the constructor)
        for (String paramName : ProbabilityTables.probabilityParams) {
            String paramStr = PTMShepherd.getParam(paramName);
            if (paramStr.length() > 0) {
                // parameter provided, read probability
                double[] values;
                switch (paramName) {
                    case "prob_neuacOx":
                        values = parseProbParam(paramStr, paramName);
                        if (values.length > 0)
                            probabilityTable.neuacRules = values;
                        break;
                    case "prob_neugcOx":
                        values = parseProbParam(paramStr, paramName);
                        if (values.length > 0)
                            probabilityTable.neugcRules = values;
                        break;
                    case "prob_phosphoOx":
                        values = parseProbParam(paramStr, paramName);
                        if (values.length > 0)
                            probabilityTable.phosphoRules = values;
                        break;
                    case "prob_sulfoOx":
                        values = parseProbParam(paramStr, paramName);
                        if (values.length > 0)
                            probabilityTable.sulfoRules = values;
                        break;
                    case "prob_dhexOx":
                        values = parseProbParam(paramStr, paramName);
                        if (values.length > 0)
                            probabilityTable.dhexOxoRules = values;
                        break;
                    case "prob_regY":
                        values = parseProbParam(paramStr, paramName);
                        if (values.length > 0)
                            probabilityTable.regularYrules = values;
                        break;
                    case "prob_dhexY":
                        values = parseProbParam(paramStr, paramName);
                        if (values.length > 0)
                            probabilityTable.dHexYrules = values;
                        break;
                    case "prob_mass":
                        probabilityTable.massProbScaling = Double.parseDouble(paramStr);
                        break;
                    case "prob_isotope":
                        String[] isoSplits = paramStr.split(",");
                        HashMap<Integer, Double> isoProbRules = new HashMap<>();
                        for (String split: isoSplits) {
                            String[] keyValue = split.trim().split(":");
                            try {
                                isoProbRules.put(Integer.parseInt(keyValue[0]), Double.parseDouble(keyValue[1]));
                            } catch (NumberFormatException ex) {
                                System.out.printf("Illegal character %s in parameter %s, must pairs of numbers like '0:1.5,1:0.75'", split, paramName);
                            }
                        }
                        if (isoProbRules.size() > 0) {
                            probabilityTable.isotopeProbTable = isoProbRules;
                        }
                        break;
                }
            }
        }
        probabilityTable.updateRulesByResidue(glycoParams);
        return probabilityTable;
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

    public TreeMap<GlycanResidue, Integer> parseGlycanStringOld(String glycanString) {
        Matcher matcher = oldPTMSPattern.matcher(glycanString);
        return parseGlycan(glycanString, matcher);
    }

    /**
     * Updated for Byonic base format (Glyc1(#)Glyc2(#) ... % Mass)
     * @param glycanString Byonic format glycan string
     * @return composition map
     */
    public TreeMap<GlycanResidue, Integer> parseGlycanString(String glycanString) {
        glycanString = glycanString.replace("FailFDR_", "");
        glycanString = glycanString.replace("Decoy_", "");
        Matcher matcher = byonicPattern.matcher(glycanString);
        return parseGlycan(glycanString, matcher);
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
        PTMShepherd.print(String.format("\tY ion probability ratio: %.1f,%.2f; dHex-containing: %.1f,%.2f", glycoProbabilityTable.regularYrules[0], glycoProbabilityTable.regularYrules[1], glycoProbabilityTable.dHexYrules[0], glycoProbabilityTable.dHexYrules[1]));
        String neuac = Arrays.toString(glycoProbabilityTable.neuacRules);
        String neugc = Arrays.toString(glycoProbabilityTable.neugcRules);
        String dhex = Arrays.toString(glycoProbabilityTable.dhexOxoRules);
        String phospho = Arrays.toString(glycoProbabilityTable.phosphoRules);
        String sulfo = Arrays.toString(glycoProbabilityTable.sulfoRules);
        PTMShepherd.print(String.format("\tOxonium probability ratios: NeuAc %s; NeuGc %s; dHex %s; Phospho %s; Sulfo %s", neuac, neugc, dhex,phospho, sulfo));
        if (printFullParams) {
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
