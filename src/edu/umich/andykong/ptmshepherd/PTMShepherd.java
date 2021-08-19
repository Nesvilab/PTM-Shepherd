package edu.umich.andykong.ptmshepherd;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.io.*;
import java.util.concurrent.*;
import java.util.Random;

import edu.umich.andykong.ptmshepherd.cleaner.CombinedExperimentsSummary;
import edu.umich.andykong.ptmshepherd.cleaner.CombinedTable;
import edu.umich.andykong.ptmshepherd.diagnosticmining.DiagnosticAnalysis;
import edu.umich.andykong.ptmshepherd.diagnosticmining.DiagnosticPeakPicker;
import edu.umich.andykong.ptmshepherd.glyco.*;
import edu.umich.andykong.ptmshepherd.localization.*;
import edu.umich.andykong.ptmshepherd.peakpicker.*;
import edu.umich.andykong.ptmshepherd.specsimilarity.*;

import static java.util.concurrent.Executors.newFixedThreadPool;

public class PTMShepherd {

	public static final String name = "PTM-Shepherd";
 	public static final String version = "1.2.0";

	static HashMap<String,String> params;
	static TreeMap<String,ArrayList<String []>> datasets;
	static HashMap<String,HashMap<String,File>> mzMap;
	static HashMap<String,Integer> datasetMS2;
	static ArrayList<GlycanCandidate> glycoDatabase;
	private static String outputPath;
	public static ExecutorService executorService;
	private static final long glycoRandomSeed = 1364955171;

	public static String getParam(String key) {
		if(params.containsKey(key))
			return params.get(key);
		else 
			return "";
	}
	public static void die(String s) {
		System.err.println("Fatal error: " + s);
		System.exit(1);
	}
	
	public static synchronized void print(String s) {
		System.out.println(s);
	}

	public static void parseParamFile(String fn) throws Exception {
		Path path = null;
		//String fn2 = fn.replaceAll("\"", "");
		//for (int i = 0; i < fn2.length(); i++) {
		//	System.out.println(fn.charAt(i) + "*");
		//}
		try {
			path = Paths.get(fn.replaceAll("['\"]", ""));
		} catch (Exception e) {
			System.out.println(e);
			die(String.format("Malformed parameter path string: [%s]", fn));
		}
		if (path == null || !Files.exists(path)) {
			die(String.format("Parameter file does not exist: [%s]", fn));
		}

		BufferedReader in = new BufferedReader(new FileReader(path.toFile()));
		String cline;
		while((cline = in.readLine())!= null) {
			int comments = cline.indexOf("//");
			if(comments >= 0)
				cline = cline.substring(0, comments);
			cline = cline.trim();
			if(cline.length() == 0 || cline.indexOf("=") < 0)
				continue;
			String key = cline.substring(0,cline.indexOf("=")).trim();
			String value = cline.substring(cline.indexOf("=")+1).trim();
			if(key.equals("dataset")) {
				StringTokenizer st = new StringTokenizer(value);
				String dsName = st.nextToken();
				String tsvTxt = st.nextToken();
				String mzPath = st.nextToken();
				if(!datasets.containsKey(dsName))
					datasets.put(dsName, new ArrayList<>());
				datasets.get(dsName).add(new String[] {tsvTxt,mzPath});
			} else {
				params.put(key, value.trim());
			}
		}
		in.close();

		if (Integer.parseInt(params.get("threads")) <= 0) {
			params.put("threads", String.valueOf(Runtime.getRuntime().availableProcessors()));
		}
	}

	/**
	 * Parse input glycan database file. Formatting: 1 glycan per line, "Residue1-count_Residue2-count_...\n"
	 * @param inputPath path to input file
	 * @return list of glycans to consider
	 */
	public static ArrayList<GlycanCandidate> parseGlycanDatabase(String inputPath, ArrayList<GlycanResidue> adductList, int maxAdducts, Random randomGenerator) {
		// read input glycan database or default database if none provided
		ArrayList<GlycanCandidate> glycanDB = new ArrayList<>();
		try {
			BufferedReader in;
			if (inputPath.equals("")) {
				// no glycan database provided - fall back to default glycan list in PeakAnnotator
				String defaultDB = "glyco_mods_20210127.txt";
				in = new BufferedReader(new InputStreamReader(PeakAnnotator.class.getResourceAsStream(defaultDB)));
			} else {
				Path path = null;
				try {
					path = Paths.get(inputPath.replaceAll("['\"]", ""));
				} catch (Exception e) {
					System.out.println(e);
					die(String.format("Malformed glycan path string: [%s]", inputPath));
				}
				if (path == null || !Files.exists(path)) {
					die(String.format("Glycan database file does not exist: [%s]", inputPath));
				}
				in = new BufferedReader(new FileReader(path.toFile()));
			}

			HashMap<String, Boolean> glycansInDB = new HashMap<>();
			// parse decoy param
			String decoyParam = getParam("decoy_type");
			int decoyType = decoyParam.length() > 0 ? Integer.parseInt(decoyParam): 0;

			String line;
			while ((line = in.readLine()) != null) {
				String glycanName;
				if (line.contains("\t")) {
					// default database includes mass in addition to glycan name - ignore mass and only take name
					glycanName = line.split("\t")[0];
				} else {
					glycanName = line;
				}
				TreeMap<GlycanResidue, Integer> glycanComp = parseGlycanString(glycanName);
				// generate a new candidate from this composition and add to DB
				GlycanCandidate candidate = new GlycanCandidate(glycanComp, false, decoyType, randomGenerator);
				String compositionHash = candidate.toHashString();
				// prevent addition of duplicates if user has them in database
				if (!glycansInDB.containsKey(compositionHash)) {
					glycanDB.add(candidate);
					glycansInDB.put(compositionHash, Boolean.TRUE);
					// also add a decoy for this composition
					GlycanCandidate decoy = new GlycanCandidate(glycanComp, true, decoyType, randomGenerator);
					glycanDB.add(decoy);

					// add adducts from adduct list to each composition
					for (GlycanResidue adduct : adductList) {
						for (int numAdducts = 1; numAdducts <= maxAdducts; numAdducts++) {
							// deep copy the original composition and add the adduct to it
							TreeMap<GlycanResidue, Integer> adductComp = new TreeMap<>();

							for (Map.Entry<GlycanResidue, Integer> previousResidue : glycanComp.entrySet()) {
								adductComp.put(previousResidue.getKey(), previousResidue.getValue());
							}
							adductComp.put(adduct, numAdducts);

							GlycanCandidate adductCandidate = new GlycanCandidate(adductComp, false, decoyType, randomGenerator);
							String adductCompositionHash = adductCandidate.toHashString();
							if (!glycansInDB.containsKey(adductCompositionHash)) {
								glycanDB.add(adductCandidate);
								glycansInDB.put(adductCompositionHash, Boolean.TRUE);
								// also add a decoy for this composition
								GlycanCandidate adductDecoy = new GlycanCandidate(adductComp, true, decoyType, randomGenerator);
								glycanDB.add(adductDecoy);
							}
						}
					}
				}
			}

		} catch (FileNotFoundException e) {
			die(String.format("Glycan database file not found: [%s]", inputPath));
		} catch (IOException e) {
			e.printStackTrace();
			die("IO Exception while reading database file");
		}
		return glycanDB;
	}

	/**
	 * Read the desired adducts from the parameters file.
	 * Parameter key: glyco_adducts
	 * Parameter values: must match GlycanResidue adduct types (case insensitive), comma separated
	 * @return list of adduct GlycanResidues to add to compositions in glycan database
	 */
	public static ArrayList<GlycanResidue> parseGlycoAdductParam() {
		ArrayList<GlycanResidue> adducts = new ArrayList<>();
		String adductParamValue = getParam("glyco_adducts");
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
		String isoLowStr = getParam("glyco_isotope_min");
		String isoHighStr = getParam("glyco_isotope_max");
		if (isoLowStr.length() > 0 && isoHighStr.length() > 0) {
			int minIso = Integer.parseInt(isoLowStr);
			int maxIso = Integer.parseInt(isoHighStr);
			ArrayList<Integer> isotopes = new ArrayList<>();
			for (int i = minIso; i <= maxIso; i++) {
				isotopes.add(i);
			}
			return isotopes.toArray(new Integer[0]);

		} else {
			// todo: deprecated. Keeping for now for compatibility with old param files, but will remove at some point
			String paramValue = getParam("glyco_isotope_range");
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
					die(String.format("Invalid isotopes string %s input to glyco mode: must be in format 'min,max'", paramValue));
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
			String paramStr = getParam(paramName);
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
		if (splits.length == 2) {
			double[] values = new double[2];
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
			System.out.printf("Invalid format for parameter %s, must have 2 comma-separated values. Param was: %s\n", paramName, paramStr);
			return new double[0];
		}
	}

	public static TreeMap<GlycanResidue, Integer> parseGlycanString(String glycanString) {
		String[] splits = glycanString.split("_");
		TreeMap<GlycanResidue, Integer> glycanComp = new TreeMap<>();
		// Read all residue counts into the composition container
		for (String split: splits) {
			String[] glycanSplits = split.split("-");
			// todo: add error catching
			GlycanResidue residue = GlycanMasses.glycoNames.get(glycanSplits[0].trim().toLowerCase(Locale.ROOT));
			int count = Integer.parseInt(glycanSplits[1].trim());
			glycanComp.put(residue, count);
		}
		return glycanComp;
	}
	
	public static void init(String [] args) throws Exception {
		if (args.length == 1) {
			if (args[0].equals("--config")) {
				printConfigFiles();
				System.exit(1);
			}
			if (args[0].equals("--annotate")) {
				printAnnotationFiles();
				System.exit(1);
			}
		}

		HashMap<String,String> overrides = new HashMap<>();
		params = new HashMap<>();
		datasets = new TreeMap<>();
		mzMap = new HashMap<>();
		datasetMS2 = new HashMap<>();
		
		//default values
		params.put("threads", String.valueOf(Runtime.getRuntime().availableProcessors()));
		params.put("histo_bindivs", "5000"); //number of divisions in histogram
		params.put("histo_smoothbins", "2"); //smoothing factor
		params.put("histo_normalizeTo", "psms"); //changing default normalization to "psms" instead of "scans"
		
		params.put("peakpicking_promRatio", "0.3"); //prominence ratio for peakpicking
		params.put("peakpicking_mass_units", "0");
		params.put("peakpicking_width", "0.002"); //width for peakpicking
		//params.put("peakpicking_background", "0.005");
		params.put("peakpicking_topN", "500"); //num peaks
        params.put("peakpicking_minPsm", "10");
        params.put("localization_background", "4");
        params.put("localization_allowed_res", "all"); //all or ABCDEF

        params.put("varmod_masses", "");
        params.put("precursor_mass_units", "0"); // 0 = Da, 1 = ppm
		params.put("precursor_tol", "0.01"); //unimod peakpicking width collapse these two parameters to one (requires redefining precursor tol in peakannotation module)
		//params.put("precursor_tol_ppm", "20.0"); //for use in mass offset and glyco modes
		params.put("annotation_tol", "0.01"); //annotation tolerance (in daltons) for unimod matching
        params.put("mass_offsets", "");
        params.put("isotope_error", "0");

		params.put("spectra_tol", "20.0"); //unused //todo
		params.put("spectra_mass_units", "1"); //unused //todo
		params.put("spectra_ppmtol", "20.0"); //obvious, used in localization and simrt
		params.put("spectra_condPeaks", "100"); //
		params.put("spectra_condRatio", "0.01");
		params.put("spectra_maxfragcharge", "2");//todo

		params.put("compare_betweenRuns", "false");
		params.put("annotation_file", "");

		params.put("glyco_mode", "false");
		params.put("cap_y_ions", "0,203.07937,406.15874,568.21156,730.26438,892.3172,349.137279");
		params.put("glyco_cap_y_ions_normalize", "1"); //0 = off, 1 = base peak
		params.put("max_cap_y_charge", "0");
		params.put("diag_ions", "204.086646,186.076086,168.065526,366.139466,144.0656,138.055,512.197375,292.1026925,274.0921325,657.2349,243.026426,405.079246,485.045576,308.09761");
		params.put("glyco_diag_ions_normalize", "1"); //0 = off, 1 = base peak
		params.put("remainder_masses", "203.07937");//,406.15874,568.21156,730.26438,892.3172,349.137279");
		params.put("remainder_mass_allowed_res", "all"); //unused

		params.put("iontype_a", "0");
		params.put("iontype_b", "1");
		params.put("iontype_c", "0");
		params.put("iontype_x", "0");
		params.put("iontype_y", "1");
		params.put("iontype_z", "0");

		params.put("output_extended", "false");
		params.put("output_path", "");
		params.put("run_from_old", "false");
		params.put("max_adducts", "1");
		
		//load parameters
		for(int i = 0; i < args.length; i++) {
			if(args[i].contains("--")) {
				overrides.put(args[i].substring(2),args[i+1]);
				i++;
			} else
				parseParamFile(args[i].trim());
		}
		params.put("peakpicking_background", Double.toString(2.5*Double.parseDouble(params.get("peakpicking_width"))));
		//replace overrides
		for(Iterator<String> it = overrides.keySet().iterator(); it.hasNext();) {
			String ckey = it.next();
			params.put(ckey, overrides.get(ckey));
		}
		
		//assertions
//		if(!params.containsKey("database"))
//			die("no database specified!");
		if(datasets.size() == 0)
			die("no datasets specified!");
		if(datasets.containsKey("combined"))
			die("combined is a reserved keyword and cannot be a dataset name!");

		if (params.get("threads").equals("0"))
			params.put("threads", String.valueOf(Runtime.getRuntime().availableProcessors()));

		if (!params.get("output_path").endsWith("/"))
			if (!params.get("output_path").equals(""))
				outputPath = params.get("output_path") + "/";
			else outputPath = "./";
		else
			outputPath = params.get("output_path");

		makeOutputDir(outputPath, Boolean.parseBoolean(params.get("output_extended")));

		executorService = Executors.newFixedThreadPool(Integer.parseInt(params.get("threads")));
	}

	private static void printAnnotationFiles() throws Exception {
		extractFile("peakpicker/glyco_mods_20210127.txt", "glyco_annotation.txt");
		extractFile("peakpicker/common_mods_20200813.txt", "common_mods_annotation.txt");
		extractFile("peakpicker/unimod_20191002.txt", "unimod_annotation.txt");
	}

	private static void printConfigFiles() throws Exception {
		extractFile("utils/open_default_params.txt", "shepherd_open_params.txt");
	}

	/**
	 * Print glyco params used
	 */
	private static void printGlycoParams(ArrayList<GlycanResidue> adductList, int maxAdducts, ArrayList<GlycanCandidate> glycoDatabase,
										 ProbabilityTables glycoProbabilityTable, boolean glycoYnorm, double absScoreErrorParam,
										 Integer[] glycoIsotopes, double glycoPPMtol, double glycoFDR, boolean printFullParams,
										 boolean nGlycanMode, String allowedLocalizationResidues) {
		print("Glycan Assignment params:");
		print(String.format("\tGlycan FDR: %.1f%%", glycoFDR * 100));
		print(String.format("\tMass error (ppm): %.1f", glycoPPMtol));
		print(String.format("\tIsotope errors: %s", Arrays.toString(glycoIsotopes)));
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
			print(String.format("\tAdducts: %s, max %d", adductStr, maxAdducts));
		} else {
			print("\tAdducts: none");
		}
		print(String.format("\tGlycan Database size (including adducts): %d", glycoDatabase.size() / 2));
		if (nGlycanMode) {
			print("\tmode: N-glycan");
			print("\tAllowed Sites: N in N-X-S/T sequons only");
		} else {
			print("\tmode: O-glycan");
			print(String.format("\tAllowed Sites: %s", allowedLocalizationResidues));
		}
		if (printFullParams) {
			print(String.format("\tNormalize Y ion counts: %s", glycoYnorm));
			print(String.format("\tTypical mass error std devs (for absolute score): %.1f", absScoreErrorParam));
			print(String.format("\tY ion probability ratio: %.1f,%.2f; dHex-containing: %.1f,%.2f", glycoProbabilityTable.regularYrules[0], glycoProbabilityTable.regularYrules[1], glycoProbabilityTable.dHexYrules[0], glycoProbabilityTable.dHexYrules[1]));
			print(String.format("\tOxonium probability ratios: NeuAc %.1f,%.2f; NeuGc %.1f,%.2f; Phospho %.1f,%.2f; Sulfo %.1f,%.2f", glycoProbabilityTable.neuacRules[0], glycoProbabilityTable.neuacRules[1], glycoProbabilityTable.neugcRules[0], glycoProbabilityTable.neugcRules[1], glycoProbabilityTable.phosphoRules[0], glycoProbabilityTable.phosphoRules[1], glycoProbabilityTable.sulfoRules[0], glycoProbabilityTable.sulfoRules[1]));
		}
		print("Assigning glycans:");
	}

	public static void main(String [] args) throws Exception {
		Locale.setDefault(new Locale("en","US"));
		System.out.println();
		System.out.printf("%s version %s",name,version);
		System.out.println("(c) University of Michigan\n");
		System.out.printf("Using Java %s on %dMB memory\n\n", System.getProperty("java.version"),(int)(Runtime.getRuntime().maxMemory()/Math.pow(2, 20)));
		
		if(args.length == 0) {
			System.out.printf("%s %s\n", name, version);
			System.out.println();
			System.out.printf("Usage:\n");
			System.out.printf("\tTo print the parameter files:\n" +
					"\t\tjava -jar ptmshepherd-%s-.jar --config\n", version);
			System.out.printf("\tTo print the annotation files:\n" +
					"\t\tjava -jar ptmshepherd-%s-.jar --annotate\n", version);
			System.out.printf("\tTo run PTM-Shepherd:\n" +
					"\t\tjava -jar ptmshepherd-%s-.jar config_file.txt\n", version);
			System.out.println();
			System.exit(0);
		}

		init(args);

		if(!Boolean.parseBoolean(params.get("run_from_old"))) {
			List<String> filesToDelete = Arrays.asList("peaks.tsv",
				"peaksummary.annotated.tsv", "peaksummary.tsv", "combined.tsv");

			for (String f : filesToDelete) {
				Path p = Paths.get(normFName(f)).toAbsolutePath().normalize();
				deleteFile(p, true);
			}

			// delete dataset files with specific extensions
			List<String> extsToDelete = Arrays
					.asList(".histo", ".locprofile.txt", ".glycoprofile.txt", ".ms2counts", ".simrtprofile.txt", ".rawlocalize", ".rawsimrt", ".rawglyco", ".modsummary.tsv");
			for (String ds : datasets.keySet()) {
				//System.out.println("Writing combined table for dataset " + ds);
				//CombinedTable.writeCombinedTable(ds);
				for (String ext : extsToDelete) {
					Path p = Paths.get(normFName(ds + ext)).toAbsolutePath().normalize();
					deleteFile(p, true);
				}
			}
			String dsTmp = "combined";
			for (String ext : extsToDelete) {
				Path p = Paths.get(outputPath + dsTmp + ext).toAbsolutePath().normalize();
				deleteFile(p, true);
			}
			dsTmp = "global";
			for (String ext : extsToDelete) {
				Path p = Paths.get(outputPath + dsTmp + ext).toAbsolutePath().normalize();
				deleteFile(p, true);
			}

		}

		//Get mzData mapping
		ArrayList<String> cacheFiles = new ArrayList<>();

		for(String ds : datasets.keySet()) {
			ArrayList<String []> dsData = datasets.get(ds);
			mzMap.put(ds, new HashMap<>());
			for(int i = 0; i < dsData.size(); i++) {
				File tpf = new File(dsData.get(i)[0]);
				String crc = PSMFile.getCRC32(tpf);
				File cacheFile = new File(normFName("cache-"+crc+".txt"));
				cacheFiles.add(crc);
				HashSet<String> fNames;
				if(!cacheFile.exists()) {
					PSMFile pf = new PSMFile(tpf);
					fNames = pf.getRunNames();
					PrintWriter out = new PrintWriter(new FileWriter(cacheFile));
					for(String cn : fNames)
						out.println(cn);
					out.close();
				} else {
					String cline;
					BufferedReader in = new BufferedReader(new FileReader(cacheFile));
					fNames = new HashSet<>();
					while((cline = in.readLine())!= null)
						fNames.add(cline);
					in.close();
				}
				for(String cname: fNames) {
					mzMap.get(ds).put(cname, null);
				}
				PSMFile.getMappings(new File(dsData.get(i)[1]), mzMap.get(ds));
				//System.out.println(mzMap.get(ds));
			}
			for(String crun : mzMap.get(ds).keySet()) {
				if(mzMap.get(ds).get(crun) == null) {
					die("In dataset \""+ds+"\" could not find mzData for run " +  crun);
				}
			}
		}
		
		//Count MS2 scans
		for(String ds : datasets.keySet()) {
			File countsFile = new File(normFName(ds+".ms2counts"));
			int sumMS2 = 0;
			TreeMap<String,Integer> counts = new TreeMap<>();
			if(!countsFile.exists()) {
				print("Counting MS2 scans for dataset " + ds);
				if (params.get("histo_normalizeTo").equals("psms")) {
					ArrayList<String []> dsData = datasets.get(ds);
					for (int i  = 0; i < dsData.size(); i++) {
						File tpf = new File(dsData.get(i)[0]);
						PSMFile pf = new PSMFile(tpf);
						counts = pf.getMS2Counts();
					}
				}
				else if (params.get("histo_normalizeTo").equals("scans")) {
					for (String crun : mzMap.get(ds).keySet()) {
						File tf = mzMap.get(ds).get(crun);
						int cnt = MS2Counts.countMS2Scans(tf, Integer.parseInt(params.get("threads")));
						print(String.format("\t%s - %d scans", crun, cnt));
						counts.put(crun, cnt);
					}
				}
				PrintWriter out = new PrintWriter(new FileWriter(countsFile));
				for (String cf : counts.keySet()) {
					sumMS2 += counts.get(cf);
					out.printf("%s\t%d\n", cf, counts.get(cf));
				}
				out.close();
			} else {
				BufferedReader in = new BufferedReader(new FileReader(countsFile));
				String cline;
				while((cline = in.readLine())!= null) {
					String [] sp = cline.split("\t");
					int v = Integer.parseInt(sp[1]);
					counts.put(sp[0], v);
					sumMS2 += v;
				}
				in.close();
			}
			for(String crun : mzMap.get(ds).keySet()) {
				if(!counts.containsKey(crun) || counts.get(crun) <= 0)
					die("Invalid MS2 counts for run " + crun + " in dataset " + ds);
			}
			datasetMS2.put(ds, sumMS2);
			print("\t" + datasetMS2.get(ds) +" MS2 scans present in dataset " + ds + "\n");
		}

		//Generate histograms
		File combinedHisto = new File(normFName("combined.histo"));
		if(!combinedHisto.exists()) {
			print("\nCreating combined histogram");
			int min = 1 << 30;
			int max = -1*(1<<30);

			for(String ds : datasets.keySet()) {
				File histoFile = new File(normFName(ds+".histo"));
				if(!histoFile.exists()) {
					ArrayList<String []> dsData = datasets.get(ds);
					ArrayList<Float> vals = new ArrayList<>();
					for(int i = 0; i < dsData.size(); i++) {
						PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
						vals.addAll(pf.getMassDiffs());
					}
					Histogram chisto = new Histogram(vals, datasetMS2.get(ds), Integer.parseInt(params.get("histo_bindivs")),Integer.parseInt(params.get("histo_smoothbins"))*2+1);
					min = Math.min(min, chisto.start);
					max = Math.max(max, chisto.end);
					chisto.writeHistogram(histoFile);
					print(String.format("\tGenerated histogram file for dataset %s [%d - %d]",ds,chisto.start,chisto.end));
				} else {
					Histogram h = Histogram.readHistogramHeader(histoFile);
					min = Math.min(min, h.start);
					max = Math.max(max, h.end);
					print(String.format("\tFound histogram file for dataset %s [%d - %d]",ds,h.start,h.end));
				}
			}

			Histogram combined = new Histogram(min,max,Integer.parseInt(params.get("histo_bindivs")));
			for(String ds : datasets.keySet()) {
				File histoFile = new File(normFName(ds+".histo"));
				Histogram h = Histogram.readHistogram(histoFile);
				combined.mergeHistogram(h, ds);
			}
			combined.writeHistogram(combinedHisto);
			combined.writeCombinedTSV(new File(normFName("combined.tsv")));
			print("Created combined histogram!\n");
		} else
			print("Combined histogram found\n");

		//Perform peak detection
		File peaks = new File(normFName("peaks.tsv"));
		if (peaks.exists()) {
			print("Deleting old peaks.tsv file at: " + peaks.toString() + "\n");
		}
		print("Running peak picking\n");
		Histogram combined = Histogram.readHistogram(combinedHisto);
		PeakPicker pp = new PeakPicker();
		pp.pickPeaks(combined.getOffsets(), combined.histo, Double.parseDouble(params.get("peakpicking_promRatio")),
				Double.parseDouble(params.get("peakpicking_background")),
				Integer.parseInt(params.get("peakpicking_topN")), params.get("mass_offsets"), params.get("isotope_error"),
				Integer.parseInt(params.get("peakpicking_mass_units")), Double.parseDouble(params.get("peakpicking_width")),
				Integer.parseInt(params.get("precursor_mass_units")), Double.parseDouble(params.get("precursor_tol")));
		pp.writeTSV(new File(normFName("peaks.tsv")));
		print("Picked top " + Integer.parseInt(params.get("peakpicking_topN")) + " peaks\n");

		//PSM assignment
		File peaksummary = new File(normFName("peaksummary.tsv"));
		if(!peaksummary.exists()) {
			PeakSummary ps = new PeakSummary(peaks,Integer.parseInt(params.get("precursor_mass_units")),
					Double.parseDouble(params.get("precursor_tol")), params.get("mass_offsets"),
					Integer.parseInt(params.get("peakpicking_minPsm")));
			for(String ds : datasets.keySet()) {
				ps.reset();
				ArrayList<String []> dsData = datasets.get(ds);
				for(int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					ps.appendPSMs(pf);
				}
				ps.commit(ds,datasetMS2.get(ds));
			}
			ps.writeTSVSummary(peaksummary);
			print("Created summary table\n");
		}

		//Assign peak IDs
		File peakannotated = new File(normFName("peaksummary.annotated.tsv"));
		if(!peakannotated.exists()) {
			PeakAnnotator pa = new PeakAnnotator();
			pa.init(params.get("varmod_masses"), params.get("annotation_file").trim());
			pa.annotateTSV(peaksummary, peakannotated, params.get("mass_offsets"), params.get("isotope_error"), Double.parseDouble(params.get("annotation_tol")));
			print("Annotated summary table\n");
			print("Mapping modifications back into PSM lists");
			for(String ds : datasets.keySet()) {
				ArrayList<String []> dsData = datasets.get(ds);
				for(int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					pa.loadAnnotatedFile(peakannotated, Double.parseDouble(params.get("precursor_tol")), Integer.parseInt(params.get("precursor_mass_units")));
					ArrayList<Float> dmasses = pf.getMassDiffs();
					ArrayList<Float> precs = pf.getPrecursorMasses();
					pf.annotateMassDiffs(pa.getDeltaMassMappings(dmasses, precs, Double.parseDouble(params.get("precursor_tol")), Integer.parseInt(params.get("precursor_mass_units"))));
				}
			}
		}

		//Mod-centric quantification
		File modSummary = new File(normFName("global.modsummary.tsv"));
		if(!modSummary.exists()) {
			ModSummary ms = new ModSummary(peakannotated, datasets.keySet());
			ms.toFile(modSummary);
			print("Created modification summary");
		}

		//Localization analysis

		//Perform initial annotation
		print("Begin localization annotation");
		for(String ds : datasets.keySet()) {
			SiteLocalization sl = new SiteLocalization(ds);
			if(sl.isComplete())
				continue;
			ArrayList<String []> dsData = datasets.get(ds);
			for(int i = 0; i < dsData.size(); i++) {
				PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
				sl.localizePSMs(pf, mzMap.get(ds));
			}
			sl.complete();
		}
		print("Done\n");

		//Localization summaries
		double [][] peakBounds = PeakSummary.readPeakBounds(peaksummary);
		LocalizationProfile loc_global = new LocalizationProfile(peakBounds, Double.parseDouble(params.get("precursor_tol")), Integer.parseInt(params.get("precursor_mass_units"))); //TODO
		for(String ds : datasets.keySet()) {
			System.out.println();
			LocalizationProfile loc_current = new LocalizationProfile(peakBounds, Double.parseDouble(params.get("precursor_tol")), Integer.parseInt(params.get("precursor_mass_units"))); //TODO
			SiteLocalization sl = new SiteLocalization(ds);
			LocalizationProfile [] loc_targets = {loc_global, loc_current};
			sl.updateLocalizationProfiles(loc_targets); //this is where the localization is actually happening
			loc_current.writeProfile(normFName(ds+".locprofile.txt"));
		}
		loc_global.writeProfile(normFName("global.locprofile.txt"));
		print("Created localization reports\n");

		//Spectra similarity analysis with retention time analysis

		//Perform similarity and RT annotation
		print("Begin similarity and retention time annotation");
		for(String ds : datasets.keySet()) {
			SimRTAnalysis sra = new SimRTAnalysis(ds);
			if(sra.isComplete())
				continue;
			ArrayList<String []> dsData = datasets.get(ds);
			for(int i = 0; i < dsData.size(); i++) {
				PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
				sra.simrtPSMs(pf, mzMap.get(ds),Boolean.parseBoolean(params.get("compare_betweenRuns")));
			}
			sra.complete();
		}
		print("Done\n");

		//SimRT summaries
		SimRTProfile simrt_global = new SimRTProfile(peakBounds, Double.parseDouble(params.get("precursor_tol")), Integer.parseInt(params.get("precursor_mass_units"))); //TODO add units
		for(String ds : datasets.keySet()) {
			SimRTProfile simrt_current = new SimRTProfile(peakBounds, Double.parseDouble(params.get("precursor_tol")), Integer.parseInt(params.get("precursor_mass_units"))); //TODO add units
			SimRTAnalysis sra = new SimRTAnalysis(ds);
			SimRTProfile [] simrt_targets = {simrt_global, simrt_current};
			sra.updateSimRTProfiles(simrt_targets);
			simrt_current.writeProfile(normFName(ds+".simrtprofile.txt"));
		}
		simrt_global.writeProfile(normFName("global.simrtprofile.txt"));
		print("Created similarity/RT reports\n");

		//Diagnostic mining
		//boolean diagMineMode = Boolean.parseBoolean(params.get("mine_diag_ions"));
		boolean diagMineMode = false;
		if (diagMineMode) {
			System.out.println("Beginning mining diagnostic ions");
			long t1 = System.currentTimeMillis();
			double peakBoundaries[][] = PeakSummary.readPeakBounds(peaksummary);
			for (String ds : datasets.keySet()) {
				System.out.println("\tPreprocessing dataset " + ds);
				DiagnosticAnalysis da = new DiagnosticAnalysis(ds);
				if (da.isComplete()) {
					System.out.println("\tFound existing data for dataset " + ds);
					continue;
				}
				da.initializeBinBoundaries(peakBoundaries);
				ArrayList<String[]> dsData = datasets.get(ds);
				for (int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					da.diagIonsPSMs(pf, mzMap.get(ds), executorService, Integer.parseInt(params.get("threads"))); //todo this is where multithreading should be done
				}
				da.complete();
				long t2 = System.currentTimeMillis();
				System.out.printf("\tFinished preprocessing dataset %s - %d ms total\n", ds, t2-t1);
			}
			float mineNoise = 0.01f; //todo cast to param
			System.out.println("\tBuilding ion histograms");
			DiagnosticPeakPicker dpp = new DiagnosticPeakPicker(mineNoise, peakBoundaries, Double.parseDouble(params.get("precursor_tol")),
					Integer.parseInt(params.get("precursor_mass_units")), "by"); //make ion types parameter
			for (String ds : datasets.keySet())
				dpp.addFileToIndex(ds);
			dpp.process();
			System.out.println("Done mining diagnostic ions");
		}

		//Combine tables
		System.out.println("Combining and cleaning reports");
		CombinedTable gct = new CombinedTable("global");
		gct.writeCombinedTable();
		for (String ds : datasets.keySet()){
			System.out.println("Writing combined table for dataset " + ds);
			CombinedTable cct = new CombinedTable(ds);
			cct.writeCombinedTable();
		}

		//Glyco analyses
		boolean glycoMode = Boolean.parseBoolean(params.get("glyco_mode"));
		if (glycoMode) {
			System.out.println("Beginning glyco/labile analysis");
			// parse glyco parameters and initialize database and ratio tables
			Random randomGenerator = new Random(glycoRandomSeed);
			ArrayList<GlycanResidue> adductList = parseGlycoAdductParam();
			int maxAdducts = Integer.parseInt(params.get("max_adducts"));
			glycoDatabase = parseGlycanDatabase(getParam("glycodatabase"), adductList, maxAdducts, randomGenerator);
			ProbabilityTables glycoProbabilityTable = initGlycoProbTable();
			boolean glycoYnorm = getParam("norm_Ys").equals("") || Boolean.parseBoolean(getParam("norm_Ys"));		// default to True if not specified
			double absScoreErrorParam = getParam("glyco_abs_score_base").equals("") ? 5.0 : Double.parseDouble(getParam("glyco_abs_score_base"));
			double glycoPPMtol = getParam("glyco_ppm_tol").equals("") ? 50.0 : Double.parseDouble(getParam("glyco_ppm_tol"));
			Integer[] glycoIsotopes = parseGlycoIsotopesParam();
			String glycoFDRParam = getParam("glyco_fdr");
			double glycoFDR = glycoFDRParam.equals("") ? 0.01 : Double.parseDouble(glycoFDRParam) / 100.0; 	// default 0.01 if param not provided, otherwise read provided value as % and convert to ratio
			boolean alreadyPrintedParams = false;
			boolean runGlycanAssignment = getParam("assign_glycans").equals("") || Boolean.parseBoolean(getParam("assign_glycans"));		// default true
			boolean printFullParams = !getParam("print_full_glyco_params").equals("") && Boolean.parseBoolean(getParam("print_full_glyco_params"));		// default false - for diagnostics
			boolean writeGlycansToAssignedMods = !getParam("put_glycans_to_assigned_mods").equals("") && Boolean.parseBoolean(getParam("put_glycans_to_assigned_mods"));	// default false
			boolean nGlycan = getParam("n_glyco").equals("") || Boolean.parseBoolean(getParam("n_glyco"));		// default true
			String allowedLocRes = PTMShepherd.getParam("localization_allowed_res");

			for (String ds : datasets.keySet()) {
				GlycoAnalysis ga = new GlycoAnalysis(ds, runGlycanAssignment, glycoDatabase, glycoProbabilityTable, glycoYnorm, absScoreErrorParam, glycoIsotopes, glycoPPMtol, randomGenerator);
				if (ga.isComplete()) {
					print(String.format("\tGlyco/labile analysis already done for dataset %s, skipping", ds));
					continue;
				}

				// print params here to avoid printing if the analysis is already done/not being run
				if (!alreadyPrintedParams && runGlycanAssignment) {
					printGlycoParams(adductList, maxAdducts, glycoDatabase, glycoProbabilityTable, glycoYnorm, absScoreErrorParam, glycoIsotopes, glycoPPMtol, glycoFDR, printFullParams, nGlycan, allowedLocRes);
					alreadyPrintedParams = true;
				}
				ArrayList<String[]> dsData = datasets.get(ds);
				for (int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					ga.glycoPSMs(pf, mzMap.get(ds));
				}
				ga.complete();
			}
			GlycoProfile glyProGLobal = new GlycoProfile(peakBounds, Integer.parseInt(params.get("precursor_mass_units")), Double.parseDouble(params.get("precursor_tol")));
			for (String ds : datasets.keySet()) {
				GlycoProfile glyProCurr = new GlycoProfile(peakBounds, Integer.parseInt(params.get("precursor_mass_units")), Double.parseDouble(params.get("precursor_tol")));
				GlycoAnalysis ga = new GlycoAnalysis(ds, runGlycanAssignment, glycoDatabase, glycoProbabilityTable, glycoYnorm, absScoreErrorParam, glycoIsotopes, glycoPPMtol, randomGenerator);
				GlycoProfile[] gaTargets = {glyProGLobal, glyProCurr};
				ga.updateGlycoProfiles(gaTargets);
				glyProCurr.writeProfile(PTMShepherd.normFName(ds + ".glycoprofile.txt"));
			}
			glyProGLobal.writeProfile(PTMShepherd.normFName("global.glycoprofile.txt"));

			// second pass: calculate glycan FDR and update results
			if (runGlycanAssignment) {
				for (String ds : datasets.keySet()) {
					GlycoAnalysis ga = new GlycoAnalysis(ds, true, glycoDatabase, glycoProbabilityTable, glycoYnorm, absScoreErrorParam, glycoIsotopes, glycoPPMtol, randomGenerator);
					ga.computeGlycanFDR(glycoFDR);
					ga.complete();
				}

				/* Save best glycan information from glyco report to psm tables */
				for (String ds : datasets.keySet()) {
					ArrayList<String[]> dsData = datasets.get(ds);
					for (int i = 0; i < dsData.size(); i++) {
						PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
						pf.mergeGlycoTable(new File(normFName(ds + ".rawglyco")), GlycoAnalysis.NUM_ADDED_GLYCO_PSM_COLUMNS, writeGlycansToAssignedMods, nGlycan, allowedLocRes);
					}
				}
			}
			print("Created glyco reports");
			print("Done with glyco/labile analysis\n");
		}

		/* Make psm table IonQuant compatible */
		if (Boolean.parseBoolean(params.get("prep_for_ionquant"))) {
			System.out.println("Prepping PSM tables for IonQuant");
			for (String ds : datasets.keySet()) {
				ArrayList<String[]> dsData = datasets.get(ds);
				for (int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					pf.preparePsmTableForIonQuant(peakBounds, Integer.parseInt(params.get("precursor_mass_units")), Double.parseDouble(params.get("precursor_tol")));
				}
			}
			System.out.println("Done");
		}

		/* Make experiment-level table */
		if (Boolean.parseBoolean(params.get("output_extended"))) {
			System.out.println("Creating experiment-level profile report");
			CombinedExperimentsSummary cs = new CombinedExperimentsSummary(normFName("combined_experiment_profile.tsv"));
			cs.initializeExperimentSummary(normFName("global.profile.tsv"));
			cs.addLocalizationSummary(normFName("global.locprofile.txt"), "combined");
			cs.addSimilarityRTSummary(normFName("global.simrtprofile.txt"), "combined");
			for (String ds : datasets.keySet()) {
				cs.addExperimentSummary(normFName(ds + ".profile.tsv"), ds);
				cs.addLocalizationSummary(normFName(ds + ".locprofile.txt"), ds);
				cs.addSimilarityRTSummary(normFName(ds+".simrtprofile.txt"), ds);
			}
			cs.printFile();
		}

		/* Move main tables out of extended subfolder */
		if (Boolean.parseBoolean(params.get("output_extended"))) {
			moveFileFromExtendedDir(normFName("global.profile.tsv"));
			moveFileFromExtendedDir(normFName("global.modsummary.tsv"));
		}

		List<String> filesToDelete = Arrays.asList(normFName("peaks.tsv"), normFName("peaksummary.annotated.tsv"),
				normFName("peaksummary.tsv"), normFName("combined.tsv"), normFName("combined.histo"));
		//delete redundant files
		for (String f : filesToDelete) {
			Path p = Paths.get(f).toAbsolutePath().normalize();
			deleteFile(p, true);
		}

		List<String> extsToDelete = Arrays.asList(".locprofile.txt", ".simrtprofile.txt", ".profile.tsv",
				".histo", ".ms2counts");
		for (String ext : extsToDelete) {
			Path p = Paths.get(normFName("global" + ext)).toAbsolutePath().normalize();
			if (!ext.equals(".profile.tsv"))
				deleteFile(p, true);
		}
		for (String ds : datasets.keySet()) {
			for (String ext : extsToDelete) {
				Path p = Paths.get(normFName(ds + ext)).toAbsolutePath().normalize();
				deleteFile(p, true);
			}
		}

		if(!Boolean.parseBoolean(params.get("output_extended"))) {
			// delete dataset files with specific extensions
			extsToDelete = Arrays
					.asList(".rawlocalize", ".rawsimrt", ".rawglyco", ".histo");
			for (String ds : datasets.keySet()) {
				//System.out.println("Writing combined table for dataset " + ds);
				//CombinedTable.writeCombinedTable(ds);
				for (String ext : extsToDelete) {
					Path p = Paths.get(normFName(ds + ext)).toAbsolutePath().normalize();
					deleteFile(p, true);
				}
			}
			String dsTmp = "combined";
			for (String ext : extsToDelete) {
				if (ext.equals(".glycoprofile.txt"))
					continue;
				Path p = Paths.get(normFName(dsTmp + ext)).toAbsolutePath().normalize();
				deleteFile(p, true);
			}
			dsTmp = "global";
			for (String ext : extsToDelete) {
				if (ext.equals(".glycoprofile.txt"))
					continue;
				Path p = Paths.get(normFName(dsTmp + ext)).toAbsolutePath().normalize();
				deleteFile(p, true);
			}

		}

		for (String crc : cacheFiles) {
			File crcf = new File(normFName("cache-"+crc +".txt"));
			Path crcpath = crcf.toPath().toAbsolutePath();
			deleteFile(crcpath, true);
		}

		executorService.shutdown();

	}

	private static void deleteFile(Path p, boolean printOnDeletion) throws IOException {
		if (Files.deleteIfExists(p) && printOnDeletion) {
			print("Deleted file: " + p.toAbsolutePath().normalize().toString());
		}
	}

	/* This method extracts compiled resources from the jar */
	private static void extractFile(String jarFilePath, String fout) throws Exception {
		BufferedReader in;
		PrintWriter out;
		System.out.printf("Copying internal resource to %s.\n", fout);
		try {
			in = new BufferedReader(new InputStreamReader(PTMShepherd.class.getResourceAsStream(jarFilePath)));
			out = new PrintWriter(new FileWriter(fout));
			String cline;
			while ((cline = in.readLine()) != null)
				out.println(cline);
			in.close();
			out.close();
		} catch (Exception ex) {
			System.out.println(ex);
			System.exit(1);
		}
	}

	/* This method adds the output directory path to file strings */
	public static String normFName(String fpath) {
		String newFPath;
		if (Boolean.parseBoolean(params.get("output_extended")))
			newFPath = new String(outputPath + "ptm-shepherd_extended_output/" + fpath);
		else
			newFPath = new String(outputPath + fpath);
		return newFPath;
	}

	/* This method will move main files out of subdirectory and into main directory */
	public static void moveFileFromExtendedDir(String oldFpath) {
		String newFpath = oldFpath.replace("ptm-shepherd_extended_output/", "");
		File oldF = new File(oldFpath);
		File newF = new File(newFpath);

		if (newF.exists()) {
			try {
				deleteFile(newF.toPath(), true);
				System.out.printf("Moving %s\n", oldF);
			} catch (Exception e) {
				System.out.printf("Error moving %s\n", oldF);
			}
		}
		try {
			oldF.renameTo(newF);
		} catch (Exception e) {
			System.out.printf("Error deleting old %s\n", oldF);
		}
	}

	/* Makes output directory */
	public static void makeOutputDir(String dpath, boolean extendedOutput) throws Exception {
		try {
			if (!dpath.equals("")) {
				File dir = new File(dpath);
				if (!dir.exists())
					dir.mkdirs();
				if (extendedOutput) {
					File eoDir = new File(dpath + "ptm-shepherd_extended_output/");
					if (!eoDir.exists())
						eoDir.mkdir();
				}
			}
		} catch (Exception e) {
			System.out.println("Error creating output directory. Terminating.");
			System.exit(1);
		}
	}
}
