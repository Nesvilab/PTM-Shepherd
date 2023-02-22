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
public class StaticGlycoUtilities {

    private static final Pattern numberPattern = Pattern.compile("[0-9]");
    private static final Pattern pGlycoPattern = Pattern.compile("[AGFHNP]");
    private static final HashMap<String, GlycanResidue> pGlycoTokenMap;    // map pGlyco tokens to our internal Glycan strings
    static
    {
        pGlycoTokenMap = new HashMap<>();
        pGlycoTokenMap.put("A", GlycanResidue.NeuAc);
        pGlycoTokenMap.put("G", GlycanResidue.NeuGc);
        pGlycoTokenMap.put("F", GlycanResidue.dHex);
        pGlycoTokenMap.put("H", GlycanResidue.Hex);
        pGlycoTokenMap.put("N", GlycanResidue.HexNAc);
        pGlycoTokenMap.put("P", GlycanResidue.Phospho);
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
     * @param glycoIsotopes
     * @return list of glycans to consider
     */
    public static ArrayList<GlycanCandidate> parseGlycanDatabase(String inputPath, ArrayList<GlycanResidue> adductList, int maxAdducts, Random randomGenerator, int decoyType, double glycoTolPPM, Integer[] glycoIsotopes, ProbabilityTables probabilityTable, HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase, boolean nGlycan) {
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
                    TreeMap<GlycanResidue, Integer> glycanComp;
                    if (dbType == DatabaseType.pGlyco) {
                        glycanComp = parsePGlycoGlycan(line);
                    } else {
                        glycanComp = parseDefaultGlycan(line);
                    }

                    if (glycanComp.size() == 0) {
                        continue;
                    }

                    // generate a new candidate from this composition and add to DB
                    GlycanCandidate candidate = new GlycanCandidate(glycanComp, false, decoyType, glycoTolPPM, glycoIsotopes, probabilityTable, glycoOxoniumDatabase, randomGenerator);
                    String compositionHash = candidate.toString();
                    // prevent addition of duplicates if user has them in database
                    if (!glycansInDB.containsKey(compositionHash)) {
                        glycanDB.add(candidate);
                        glycansInDB.put(compositionHash, Boolean.TRUE);
                        // also add a decoy for this composition
                        GlycanCandidate decoy = new GlycanCandidate(glycanComp, true, decoyType, glycoTolPPM, glycoIsotopes, probabilityTable, glycoOxoniumDatabase, randomGenerator);
                        glycanDB.add(decoy);

                        // add adducts from adduct list to each composition
                        for (GlycanResidue adduct : adductList) {
                            for (int numAdducts = 1; numAdducts <= maxAdducts; numAdducts++) {
                                // deep copy the original composition and add the adduct to it
                                TreeMap<GlycanResidue, Integer> adductComp = new TreeMap<>();

                                adductComp.putAll(glycanComp);
                                adductComp.put(adduct, numAdducts);

                                GlycanCandidate adductCandidate = new GlycanCandidate(adductComp, false, decoyType, glycoTolPPM, glycoIsotopes, probabilityTable, glycoOxoniumDatabase, randomGenerator);
                                String adductCompositionHash = adductCandidate.toString();
                                if (!glycansInDB.containsKey(adductCompositionHash)) {
                                    glycanDB.add(adductCandidate);
                                    glycansInDB.put(adductCompositionHash, Boolean.TRUE);
                                    // also add a decoy for this composition
                                    GlycanCandidate adductDecoy = new GlycanCandidate(adductComp, true, decoyType, glycoTolPPM, glycoIsotopes, probabilityTable, glycoOxoniumDatabase, randomGenerator);
                                    glycanDB.add(adductDecoy);
                                }
                            }
                        }
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
    public static ArrayList<GlycanCandidate> updateGlycanDatabase(HashMap<String, GlycanCandidateFragments> fragmentDB, ArrayList<GlycanCandidate> oldGlycoDB, Random randomGenerator, int decoyType, double glycoTolPPM, Integer[] glycoIsotopes, ProbabilityTables probabilityTable, HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase) {
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
            newCandidate = new GlycanCandidate(oldCandidate, fragmtInfo, decoyType, glycoTolPPM, glycoIsotopes, randomGenerator, glycoOxoniumDatabase);
            newGlycoDB.add(newCandidate);
        }
        return newGlycoDB;
    }

    /**
     * Parse new (Byonic) or old PTM-Shepherd format glycan from String
     * @param line
     * @return
     */
    private static TreeMap<GlycanResidue, Integer> parseDefaultGlycan(String line) {
        TreeMap<GlycanResidue, Integer> glycanComp;
        String glycanName;
        // handle new and old formats
        if (line.contains("%") || line.contains("(")) {
            if (line.contains("%")) {
                // default database includes mass in addition to glycan name - ignore mass and only take name
                glycanName = line.split("%")[0];
            } else {
                glycanName = line;
            }
            glycanComp = parseGlycanString(glycanName);
        } else {
            if (line.contains("\t")) {
                // default database includes mass in addition to glycan name - ignore mass and only take name
                glycanName = line.split("\t")[0];
            } else {
                glycanName = line;
            }
            glycanComp = parseGlycanStringOld(glycanName);
        }
        return glycanComp;
    }

    /**
     * Parse pGlyco glycan from String
     * @param line
     * @return
     */
    private static TreeMap<GlycanResidue, Integer> parsePGlycoGlycan(String line) {
        TreeMap<GlycanResidue, Integer> glycanComp = new TreeMap<>();

        // skip empty lines and lines without a glycan structure
        if (!line.contains("(")) {
            return glycanComp;
        }

        // parse glycans from tokens
        line = line.replace("Hp", "H(P");      // remove double characters, if present
        Matcher matcher = pGlycoPattern.matcher(line);
        boolean valid = true;
        while(matcher.find()) {
            String glycanToken = matcher.group();
            if (pGlycoTokenMap.containsKey(glycanToken)) {
                GlycanResidue residue = pGlycoTokenMap.get(glycanToken);
                if (glycanComp.containsKey(residue)) {
                    glycanComp.put(residue, glycanComp.get(residue) + 1);
                } else {
                    glycanComp.put(residue, 1);
                }

            } else {
                PTMShepherd.print(String.format("Invalid token %s in line %s during glycan database parsing. This glycan will be skipped", glycanToken, line));
                valid = false;
                break;
            }
        }
        if (!valid) {
            glycanComp = new TreeMap<>();
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
                if (line.contains("#") || line.contains("\\\\")) {
                    i++;
                    continue;
                }
                // pGlyco header is of format "N,H,F,...". This won't always detect it if only N and H types are present
                if (line.contains("N,") && line.contains("H,")) {
                    return DatabaseType.pGlyco;
                }
            } else {
                // non-header line
                if (line.contains("%")) {
                    return DatabaseType.byonic;
                } else {
                    if (line.contains("-")) {
                        return DatabaseType.oldPTMShepherd;
                    } else {
                        // pGlyco format has no numbers, byonic format has
                        Matcher matcher = numberPattern.matcher(line);
                        if (matcher.find()) {
                            return DatabaseType.byonic;
                        } else {
                            return DatabaseType.pGlyco;
                        }
                    }
                }
            }
            i++;
        }
        return DatabaseType.unknown;
    }

    /**
     * Read the desired adducts from the parameters file.
     * Parameter key: glyco_adducts
     * Parameter values: must match GlycanResidue adduct types (case insensitive), comma separated
     * @return list of adduct GlycanResidues to add to compositions in glycan database
     */
    public static ArrayList<GlycanResidue> parseGlycoAdductParam() {
        ArrayList<GlycanResidue> adducts = new ArrayList<>();
        String adductParamValue = PTMShepherd.getParam("glyco_adducts");
        String[] adductStrs;
        if (adductParamValue.length() > 0)
            adductStrs = adductParamValue.split(",| |/");
        else
            adductStrs = new String[0];

        for (String adductStr : adductStrs) {
            if (adductStr.equals("")) {
                continue;
            }
            if (GlycanMasses.glycoNames.containsKey(adductStr.trim())) {
                GlycanResidue adduct = GlycanMasses.glycoNames.get(adductStr.trim());
                adducts.add(adduct);
            } else {
                System.out.printf("Invalid glyco adduct %s, ignored\n", adductStr);
            }
        }
        return adducts;
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
    public static ProbabilityTables initGlycoProbTable() {
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
        probabilityTable.updateRulesByResidue();
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

    public static TreeMap<GlycanResidue, Integer> parseGlycanStringOld(String glycanString) {
        String[] splits = glycanString.split("_");
        TreeMap<GlycanResidue, Integer> glycanComp = new TreeMap<>();
        // Read all residue counts into the composition container
        for (String split: splits) {
            String[] glycanSplits = split.split("-");
            parserHelper(glycanComp, glycanSplits[0], glycanSplits[1], glycanString);
        }
        return glycanComp;
    }

    /**
     * Updated for Byonic base format (Glyc1(#)Glyc2(#) ... % Mass)
     * @param glycanString Byonic format glycan string
     * @return composition map
     */
    public static TreeMap<GlycanResidue, Integer> parseGlycanString(String glycanString) {
        glycanString = glycanString.replace("FailFDR_", "");
        glycanString = glycanString.replace("Decoy_", "");
        String[] massSplits = glycanString.split(" % ");
        String[] compositionSplits = massSplits[0].trim().split("\\)");
        TreeMap<GlycanResidue, Integer> glycanComp = new TreeMap<>();
        // Read all residue counts into the composition container
        for (String split: compositionSplits) {
            String[] glycanSplits = split.split("\\(");
            parserHelper(glycanComp, glycanSplits[0], glycanSplits[1], glycanString);
        }
        return glycanComp;
    }

    /**
     * Helper for parsing glycan compositions with error checking. NOTE: ignores unknown glycan residues, but continues
     * assembling the rest of the composition, leading to altered compositions when unknowns are present
     * @param glycanComp
     * @param residueStr
     * @param countStr
     * @param glycanString
     */
    private static void parserHelper(TreeMap<GlycanResidue, Integer> glycanComp, String residueStr, String countStr, String glycanString) {
        GlycanResidue residue = GlycanMasses.glycoNames.get(residueStr.trim().toLowerCase(Locale.ROOT));
        if (residue == null) {
            PTMShepherd.print(String.format("Warning: glycan type %s in glycan %s not recognized and was ignored", residueStr, glycanString));
            return;
        }
        int count;
        try {
            count = Integer.parseInt(countStr.trim());
        } catch (ArrayIndexOutOfBoundsException ex) {
            return;   // empty string or incorrect format - ignore
        }
        glycanComp.put(residue, count);
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
    public static void printGlycoParams(ArrayList<GlycanResidue> adductList, int maxAdducts, ArrayList<GlycanCandidate> glycoDatabase,
                                        ProbabilityTables glycoProbabilityTable, boolean glycoYnorm, double absScoreErrorParam,
                                        Integer[] glycoIsotopes, double glycoPPMtol, double glycoFDR, boolean printFullParams,
                                        boolean nGlycanMode, String allowedLocalizationResidues, int decoyType,
                                        boolean printDecoys, boolean removeGlycanDelta) {
        PTMShepherd.print("Glycan Assignment params:");
        PTMShepherd.print(String.format("\tGlycan FDR: %.1f%%", glycoFDR * 100));
        PTMShepherd.print(String.format("\tMass error (ppm): %.1f", glycoPPMtol));
        PTMShepherd.print(String.format("\tIsotope errors: %s", Arrays.toString(glycoIsotopes)));
        // adducts
        if (maxAdducts > 0 && adductList.size() > 0) {
            StringBuilder adductStr = new StringBuilder();
            int count = 0;
            for (GlycanResidue adduct: adductList) {
                adductStr.append(adduct.name());
                count++;
                if (count < adductList.size()) {
                    adductStr.append(", ");
                }
            }
            PTMShepherd.print(String.format("\tAdducts: %s, max %d", adductStr, maxAdducts));
        } else {
            PTMShepherd.print("\tAdducts: none");
        }
        PTMShepherd.print(String.format("\tGlycan Database size (including adducts): %d", glycoDatabase.size() / 2));
        if (nGlycanMode) {
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
            if (printDecoys) {
                PTMShepherd.print("\tPrinting decoy glycans");
            }
            if (removeGlycanDelta) {
                PTMShepherd.print("\tRemoving glycan delta mass from PSM table");
            }
        }
        PTMShepherd.print("Assigning glycans:");
    }
}
