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

package edu.umich.andykong.ptmshepherd;

import static java.lang.System.out;

import edu.umich.andykong.ptmshepherd.cleaner.CombinedExperimentsSummary;
import edu.umich.andykong.ptmshepherd.cleaner.CombinedTable;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.diagnosticanalysis.DiagnosticExtractor;
import edu.umich.andykong.ptmshepherd.diagnosticmining.DiagnosticAnalysis;
import edu.umich.andykong.ptmshepherd.diagnosticmining.DiagnosticPeakPicker;
import edu.umich.andykong.ptmshepherd.glyco.GlycanCandidate;
import edu.umich.andykong.ptmshepherd.glyco.GlycanCandidateFragments;
import edu.umich.andykong.ptmshepherd.glyco.GlycanFragmentDescriptor;
import edu.umich.andykong.ptmshepherd.glyco.GlycanResidue;
import edu.umich.andykong.ptmshepherd.glyco.GlycoAnalysis;
import edu.umich.andykong.ptmshepherd.glyco.GlycoProfile;
import edu.umich.andykong.ptmshepherd.glyco.ProbabilityTables;
import edu.umich.andykong.ptmshepherd.glyco.StaticGlycoUtilities;
import edu.umich.andykong.ptmshepherd.iterativelocalization.IterativeLocalizer;
import edu.umich.andykong.ptmshepherd.localization.LocalizationProfile;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;
import edu.umich.andykong.ptmshepherd.peakpicker.Histogram;
import edu.umich.andykong.ptmshepherd.peakpicker.MS2Counts;
import edu.umich.andykong.ptmshepherd.peakpicker.ModSummary;
import edu.umich.andykong.ptmshepherd.peakpicker.PeakAnnotator;
import edu.umich.andykong.ptmshepherd.peakpicker.PeakPicker;
import edu.umich.andykong.ptmshepherd.peakpicker.PeakSummary;
import edu.umich.andykong.ptmshepherd.specsimilarity.SimRTAnalysis;
import edu.umich.andykong.ptmshepherd.specsimilarity.SimRTProfile;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import umich.ms.fileio.filetypes.mzbin.MZBINFile;
import umich.ms.fileio.filetypes.mzbin.MZBINFile.MZBINSpectrum;

public class PTMShepherd {

	public static final String name = "PTM-Shepherd";
 	public static final String version = "2.1.0";

	static HashMap<String,String> params;
	static TreeMap<String,ArrayList<String []>> datasets;
	static HashMap<String,HashMap<String,File>> mzMap;
	static HashMap<String,Integer> datasetMS2;
	static ArrayList<String> cacheFiles;
	static ArrayList<GlycanCandidate> glycoDatabase;
	private static String outputPath;
	public static ExecutorService executorService;
	private static final long glycoRandomSeed = 1364955171;

	// filenames for output files
	public static final String outputDirName = "ptm-shepherd-output/";
	public static final String globalName = "global";
	public static final String combinedName = "combined";
	public static final String mzBinFilename = ".mzBIN_cache";
	public static final String ms2countsName = ".ms2counts";
	public static final String locProfileName =  ".locprofile.txt";
	public static final String simRTProfileName =  ".simrtprofile.txt";
	public static final String profileName =  ".profile.tsv";
	public static final String histoName = ".histo";
	public static final String peaksName = "peaks.tsv";
	public static final String peakSummaryAnnotatedName = "peaksummary.annotated.tsv";
	public static final String peakSummaryName = "peaksummary.tsv";
	public static final String combinedTSVName = "combined.tsv";
	public static final String combinedHistoName = "combined.histo";
	public static final String glycoProfileName = ".diagnosticProfile.txt";
	public static final String rawLocalizeName = ".rawlocalize";
	public static final String rawSimRTName = ".rawsimrt";
	public static final String rawGlycoName = ".rawglyco";
	public static final String modSummaryName = ".modsummary.tsv";
	public static final String diagBinFilename = ".diagBIN";
	public static final String diagMineName = ".diagmine.tsv";
	public static final String diagIonsExtractName = ".diagnosticIons.tsv";
	public static final String glycoMassListName = "glyco_masses_list.txt";


	public static String getParam(String key) {
		return params.getOrDefault(key, "");
	}
	public static void die(String s) {
		System.err.println("Fatal error: " + s);
		System.exit(1);
	}
	
	public static synchronized void print(String s) {
		out.println(s);
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
			out.println(e);
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
		params.put("histo_intensity", "0"); //use MS1 intensities instead of
		params.put("histo_Tmt", "0"); //uses TMT values to create multiple experiments
		
		params.put("peakpicking_promRatio", "0.3"); //prominence ratio for peakpicking
		params.put("peakpicking_mass_units", "0");
		params.put("peakpicking_width", "0.002"); //width for peakpicking
		//params.put("peakpicking_background", "0.005");
		params.put("peakpicking_topN", "500"); //num peaks
        params.put("peakpicking_minPsm", "10");
        params.put("localization_background", "4");
        params.put("localization_allowed_res", ""); //all or ABCDEF

        params.put("varmod_masses", "");
        params.put("precursor_mass_units", "0"); // 0 = Da, 1 = ppm
		params.put("precursor_tol", "0.01"); //unimod peakpicking width collapse these two parameters to one (requires redefining precursor tol in peakannotation module)
		params.put("precursor_maxCharge", "1");
		//params.put("precursor_tol_ppm", "20.0"); //for use in mass offset and glyco modes

        params.put("mass_offsets", "");
        params.put("isotope_error", "0");

		params.put("spectra_tol", "20.0"); //unugsed //todo
		params.put("spectra_mass_units", "1"); //unused //todo
		params.put("spectra_ppmtol", "20.0"); //obvious, used in localization and simrt
		params.put("spectra_condPeaks", "150"); //
		params.put("spectra_condRatio", "0.00001");
		params.put("spectra_maxfragcharge", "2");//todo
		params.put("spectra_maxPrecursorCharge", "4");
		//TODO add precursor removal

		params.put("compare_betweenRuns", "true");
		params.put("annotation_file", "");
		params.put("annotation_tol", "0.01"); //annotation tolerance (in daltons) for unimod matching

		params.put("glyco_mode", "false");
		params.put("cap_y_ions", "0,203.07937,349.137279,406.15874,568.21156,730.26438,892.3172");
		params.put("glyco_cap_y_ions_normalize", "1"); //0 = off, 1 = base peak
		params.put("max_cap_y_charge", "0");
		params.put("diag_ions", "144.0656,138.055,168.065526,186.076086,204.086646,243.026426,274.0921325,292.1026925,308.09761,366.139466,405.079246,485.045576,512.197375,657.2349");
		params.put("glyco_diag_ions_normalize", "1"); //0 = off, 1 = base peak
		params.put("remainder_masses", "203.07937");//,406.15874,568.21156,730.26438,892.3172,349.137279");
		params.put("remainder_mass_allowed_res", "all"); //unused

		params.put("iontype_a", "0");
		params.put("iontype_b", "1");
		params.put("iontype_c", "0");
		params.put("iontype_x", "0");
		params.put("iontype_y", "1");
		params.put("iontype_z", "0");
		
		params.put("iterloc_mode", "true");
		params.put("iterloc_convergeCriterion", "0.01");
		params.put("iterloc_maxEpoch", "100"); // Todo hangs when set to 0

		params.put("diagmine_mode", "false");
		params.put("diagmine_minSignal", "0.001");
		params.put("diagmine_filterIonTypes", "aby");
		params.put("diagmine_ionTypes", "by");
		params.put("diagmine_maxP", "0.05");
		params.put("diagmine_minAuc", "0.01"); //todo deprecated
		params.put("diagmine_minSpecDiff", "0.25"); //todo deprecated
		params.put("diagmine_minFoldChange", "2.0"); //todo deprecated
		params.put("diagmine_pepMinSpecDiff", "25.0");
		params.put("diagmine_pepMinFoldChange", "3.0");
		params.put("diagmine_diagMinFoldChange", "3.0");
		params.put("diagmine_diagMinSpecDiff", "25.0");
		params.put("diagmine_fragMinSpecDiff", "15.0");
		params.put("diagmine_fragMinPropensity", "12.5");
		params.put("diagmine_fragMinFoldChange", "3.0");
		params.put("diagmine_minPeps", "25");
		params.put("diagmine_twoTailedTests", "0");
		params.put("diagmine_maxPsms", "1000");
		params.put("diagmine_minIonsPerSpec", "2");
		params.put("diagmine_printHistos", "false");
		params.put("diagmine_printDebug", "false");
		params.put("diagmine_printDebugFile", "false");
		params.put("diagmine_printIsotopes", "false");
		params.put("diagmine_printRedundantTests", "0");
		params.put("diagmine_printAuc", "true");

		params.put("output_extended", "false");
		params.put("output_extended_diagmine", "false");
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

		if(datasets.size() == 0)
			die("no datasets specified!");
		if(datasets.containsKey(combinedName))
			die(String.format("%s is a reserved keyword and cannot be a dataset name!", combinedName));

		if (params.get("threads").equals("0"))
			params.put("threads", String.valueOf(Runtime.getRuntime().availableProcessors()));

		if (!params.get("output_path").endsWith("/"))
			if (!params.get("output_path").equals(""))
				outputPath = params.get("output_path") + "/";
			else outputPath = "./";
		else
			outputPath = params.get("output_path");

		makeOutputDir(outputPath);

		executorService = Executors.newFixedThreadPool(Integer.parseInt(params.get("threads")));
	}

	private static void printAnnotationFiles() throws Exception {
		extractFile("peakpicker/glyco_mods_20210127.txt", "glyco_annotation.txt");
		extractFile("peakpicker/common_mods_20200813.txt", "common_mods_annotation.txt");
		extractFile("peakpicker/unimod_20221028.txt", "unimod_annotation.txt");
	}

	private static void printConfigFiles() throws Exception {
		extractFile("utils/open_default_params.txt", "shepherd_open_params.txt");
	}

	public static void main(String [] args) throws Exception {
		Locale.setDefault(new Locale("en","US"));
		out.println();
		out.printf("%s version %s\n",name,version);
		out.println("(c) 2022 University of Michigan\n");
		out.printf("Using Java %s on %dMB memory\n\n", System.getProperty("java.version"),(int)(Runtime.getRuntime().maxMemory()/Math.pow(2, 20)));
		
		if(args.length == 0) {
			out.printf("%s %s\n", name, version);
			out.println();
			out.printf("Usage:\n");
			out.printf("\tTo print the parameter files:\n" +
					"\t\tjava -jar ptmshepherd-%s-.jar --config\n", version);
			out.printf("\tTo print the annotation files:\n" +
					"\t\tjava -jar ptmshepherd-%s-.jar --annotate\n", version);
			out.printf("\tTo run PTM-Shepherd:\n" +
					"\t\tjava -jar ptmshepherd-%s-.jar config_file.txt\n", version);
			out.println();
			System.exit(0);
		}

		init(args);

		//TODO initialize program blocks here so that they can be accessed outside their modules and stored internally
		PeakAnnotator pa = new PeakAnnotator();

		//Get mzData mapping
		print("Finding spectral data");
		getMzDataMapping();
		print("Done finding spectral data\n");

		//After knowing where all files should be, remove ones from old runs
		if(!Boolean.parseBoolean(params.get("run_from_old"))) {
			List<String> filesToDelete = Arrays.asList(peaksName,
				peakSummaryAnnotatedName, peakSummaryName, combinedTSVName);
			for (String f : filesToDelete) {
				Path p = Paths.get(normFName(f)).toAbsolutePath().normalize();
				deleteFile(p, true);
			}

			// delete dataset files with specific extensions
			List<String> extsToDelete = Arrays
					.asList(histoName, locProfileName, glycoProfileName, ms2countsName, simRTProfileName, rawLocalizeName, rawSimRTName, rawGlycoName, modSummaryName);
			for (String ds : datasets.keySet()) {
				for (String ext : extsToDelete) {
					Path p = Paths.get(normFName(ds + ext)).toAbsolutePath().normalize();
					deleteFile(p, true);
				}
				for (String mzFile : mzMap.get(ds).keySet()) {
					Path p = Paths.get(normFName(mzFile + diagBinFilename)).toAbsolutePath().normalize();
					deleteFile(p, true);
				}
			}
			String dsTmp = combinedName;
			for (String ext : extsToDelete) {
				Path p = Paths.get(outputPath + dsTmp + ext).toAbsolutePath().normalize();
				deleteFile(p, true);
			}
			dsTmp = globalName;
			for (String ext : extsToDelete) {
				Path p = Paths.get(outputPath + dsTmp + ext).toAbsolutePath().normalize();
				deleteFile(p, true);
			}
			print("");
		}

		// Cache spectral data
		print("Caching spectral data");
		rewriteDataToMzBin();
		print("Done caching spectral data\n");

		
		//Count MS2 scans
		for(String ds : datasets.keySet()) {
			File countsFile = new File(normFName(ds+ms2countsName));
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
				} else if (params.get("histo_normalizeTo").equals("scans")) {
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
			print("\t" + datasetMS2.get(ds) +" MS2 scans present in dataset " + ds);
			print("Done counting MS2 scans for dataset " + ds + "\n");
		}

		//Generate histograms
		File combinedHisto = new File(normFName(combinedHistoName));
		if(!combinedHisto.exists()) {
			print("Creating combined histogram");
			int min = 1 << 30;
			int max = -1*(1<<30);

			for(String ds : datasets.keySet()) {
				File histoFile = new File(normFName(ds+histoName));
				if(!histoFile.exists()) {
					ArrayList<String []> dsData = datasets.get(ds);
					ArrayList<Float> vals = new ArrayList<>();
					ArrayList<Double> ints =  new ArrayList<>();
					for(int i = 0; i < dsData.size(); i++) {
						PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
						vals.addAll(pf.getMassDiffs());
						ints.addAll(pf.getIntensities());
					}
					Histogram chisto = new Histogram(vals, ints, datasetMS2.get(ds), Integer.parseInt(params.get("histo_bindivs")),Integer.parseInt(params.get("histo_smoothbins"))*2+1);
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
				File histoFile = new File(normFName(ds+histoName));
				Histogram h = Histogram.readHistogram(histoFile);
				combined.mergeHistogram(h, ds);
			}
			combined.writeHistogram(combinedHisto);
			combined.writeCombinedTSV(new File(normFName(combinedTSVName)));
			print("Done creating combined histogram\n");
		} else
			print("Combined histogram found\n");

		//Perform peak detection
		File peaks = new File(normFName(peaksName));
		if (peaks.exists()) {
			print("Deleting old peaks.tsv file at: " + peaks.toString() + "\n");
		}
		print("Running peak picking");
		Histogram combined = Histogram.readHistogram(combinedHisto);
		PeakPicker pp = new PeakPicker();
		pp.pickPeaks(combined.getOffsets(), combined.histo, Double.parseDouble(params.get("peakpicking_promRatio")),
				Double.parseDouble(params.get("peakpicking_background")),
				Integer.parseInt(params.get("peakpicking_topN")), params.get("mass_offsets"), params.get("isotope_error"),
				Integer.parseInt(params.get("peakpicking_mass_units")), Double.parseDouble(params.get("peakpicking_width")),
				Integer.parseInt(params.get("precursor_mass_units")), Double.parseDouble(params.get("precursor_tol")));
		pp.writeTSV(new File(normFName(peaksName)));
		print("\tPicked top " + pp.getPeaks().length + " peaks");

		//PSM assignment
		File peaksummary = new File(normFName(peakSummaryName));
		if(!peaksummary.exists()) {
			PeakSummary ps = new PeakSummary(peaks,Integer.parseInt(params.get("precursor_mass_units")),
					Double.parseDouble(params.get("precursor_tol")), params.get("mass_offsets"),
					Integer.parseInt(params.get("peakpicking_minPsm")), Integer.parseInt(params.get("histo_intensity")));
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
			print("\tCreated summary table");
		}
		print("Done running peak picking\n");

		//Assign peak IDs
		File peakannotated = new File(normFName(peakSummaryAnnotatedName));
		if(!peakannotated.exists()) {
			pa.init(params.get("varmod_masses"), params.get("annotation_file").trim());
			pa.annotateTSV(peaksummary, peakannotated, params.get("mass_offsets"), params.get("isotope_error"), Double.parseDouble(params.get("annotation_tol")));
			print("Annotated summary table\n");
			print("Mapping modifications back into PSM lists\n");
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
			File modSummary = new File(normFName(globalName + modSummaryName));
			ModSummary ms = new ModSummary(peakannotated, datasets.keySet());
			ms.toFile(modSummary);
			print("Created modification summary\n");
		}

		//PTMiner-style iterative localization
		Boolean iterLocMode = Boolean.parseBoolean(params.get("iterloc_mode"));
		if (iterLocMode) {
			out.println("Beginning iterative localization");
			long t1 = System.currentTimeMillis();
			double peakBoundaries[][] = PeakSummary.readPeakBounds(peaksummary);
			IterativeLocalizer IterLoc = new IterativeLocalizer(peakBoundaries,
					Double.parseDouble(params.get("precursor_tol")),
					Integer.parseInt(params.get("precursor_mass_units")),
					datasets, mzMap, Integer.parseInt(params.get("threads")),
					params.get("localization_allowed_res"), Float.parseFloat(params.get("spectra_tol")),
					concatIonTypes(), Double.parseDouble(params.get("iterloc_convergeCriterion")),
					Integer.parseInt(params.get("iterloc_maxEpoch"))
			);
			IterLoc.localize();
			out.println("Done\n");
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
			LocalizationProfile loc_current = new LocalizationProfile(peakBounds, Double.parseDouble(params.get("precursor_tol")), Integer.parseInt(params.get("precursor_mass_units"))); //TODO
			SiteLocalization sl = new SiteLocalization(ds);
			LocalizationProfile [] loc_targets = {loc_global, loc_current};
			sl.updateLocalizationProfiles(loc_targets); //this is where the localization is actually happening
			loc_current.writeProfile(normFName(ds+locProfileName));
		}
		loc_global.writeProfile(normFName(globalName + locProfileName));
		print("Created localization reports\n");

		//Spectra similarity analysis with retention time analysis

		//Perform similarity and RT annotation
		print("Begin similarity and retention time annotation");
		boolean calcIntensity = false;
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
			calcIntensity = sra.getCalcIntensity();
		}
		print("Done\n");

		//SimRT summaries
		SimRTProfile simrt_global = new SimRTProfile(peakBounds, Double.parseDouble(params.get("precursor_tol")), Integer.parseInt(params.get("precursor_mass_units")), calcIntensity); //TODO add units
		for(String ds : datasets.keySet()) {
			SimRTProfile simrt_current = new SimRTProfile(peakBounds, Double.parseDouble(params.get("precursor_tol")), Integer.parseInt(params.get("precursor_mass_units")), calcIntensity); //TODO add units
			SimRTAnalysis sra = new SimRTAnalysis(ds);
			SimRTProfile [] simrt_targets = {simrt_global, simrt_current};
			sra.updateSimRTProfiles(simrt_targets);
			simrt_current.writeProfile(normFName(ds+simRTProfileName));
		}
		simrt_global.writeProfile(normFName(globalName + simRTProfileName));
		print("Created similarity/RT reports\n");

		//Diagnostic mining
		boolean diagMineMode = Boolean.parseBoolean(params.get("diagmine_mode"));
		if (diagMineMode) {
			out.println("Beginning diagnostic ion mining");
			long t1 = System.currentTimeMillis();
			double peakBoundaries[][] = PeakSummary.readPeakBounds(peaksummary);
			for (String ds : datasets.keySet()) {
				out.println("\tPreprocessing dataset " + ds);
				DiagnosticAnalysis da = new DiagnosticAnalysis(ds);
				//if (da.isComplete()) {
				//	System.out.println("\tFound existing data for dataset " + ds);
				//	continue;
				//}
				da.initializeBinBoundaries(peakBoundaries);
				ArrayList<String[]> dsData = datasets.get(ds);
				for (int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					da.diagIonsPSMs(pf, mzMap.get(ds), executorService, Integer.parseInt(params.get("threads")));
				}
				//da.complete();
				long t2 = System.currentTimeMillis();
				out.printf("\tDone preprocessing dataset %s - %d ms total\n", ds, t2-t1);
			}
			DiagnosticPeakPicker dpp = new DiagnosticPeakPicker(Double.parseDouble(getParam("diagmine_minSignal")), peakBoundaries, Double.parseDouble(params.get("precursor_tol")),
					Integer.parseInt(params.get("precursor_mass_units")), params.get("diagmine_ionTypes"),Float.parseFloat(params.get("spectra_tol")), Integer.parseInt(params.get("precursor_maxCharge")),
					Double.parseDouble(params.get("diagmine_maxP")), Double.parseDouble(params.get("diagmine_minAuc")), Double.parseDouble(params.get("diagmine_minSpecDiff")), Double.parseDouble(params.get("diagmine_minFoldChange")), Integer.parseInt(params.get("diagmine_minIonsPerSpec")),
					Integer.parseInt(params.get("diagmine_twoTailedTests")), Integer.parseInt(params.get("spectra_condPeaks")), Double.parseDouble(params.get("spectra_condRatio")), pa.getPeakApexMappings());
			for (String ds : datasets.keySet()) {
				ArrayList<String[]> dsData = datasets.get(ds);
				for (int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					dpp.addPepkeysToIndex(pf);
				}
			}
			out.println("\tIdentifying candidate ions");
			dpp.filterPepkeys();
			dpp.process(executorService, Integer.parseInt(getParam("threads")));
			out.println("\tDone identifying candidate ions");
			out.println("\tExtracting ions from spectra");
			dpp.initDiagProfRecs();
			for (String ds : datasets.keySet()) {
				ArrayList<String[]> dsData = datasets.get(ds);
				for (int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					dpp.diagIonsPSMs(pf, mzMap.get(ds), executorService);
				}
			}
			dpp.print(normFName(globalName + diagMineName));
			out.println("Done mining diagnostic ions\n");
		}

		// diagnostic ion extraction (original glyco/labile mode)
		boolean extractDiagnosticIons = Boolean.parseBoolean(params.get("diagextract_mode"));
		if (extractDiagnosticIons) {
			System.out.println("Beginning diagnostic ion extraction");
			int numThreads = Integer.parseInt(params.get("threads"));
			for (String ds : datasets.keySet()) {
				DiagnosticExtractor da = new DiagnosticExtractor(ds);
				if (da.isDiagnosticComplete()) {
					print(String.format("\tDiagnostic extraction already done for dataset %s, skipping", ds));
					continue;
				}
				ArrayList<String[]> dsData = datasets.get(ds);
				for (int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					da.extractDiagPSMs(pf, mzMap.get(ds), executorService, numThreads);
				}
				da.completeDiagnostic();
			}
			// calculate glycoprofile after all other diagnostic extraction analysis is done
			GlycoProfile glyProGLobal = new GlycoProfile(peakBounds, Integer.parseInt(params.get("precursor_mass_units")), Double.parseDouble(params.get("precursor_tol")));
			for (String ds : datasets.keySet()) {
				GlycoProfile glyProCurr = new GlycoProfile(peakBounds, Integer.parseInt(params.get("precursor_mass_units")), Double.parseDouble(params.get("precursor_tol")));
				DiagnosticExtractor da = new DiagnosticExtractor(ds);
				GlycoProfile[] gaTargets = {glyProGLobal, glyProCurr};
				da.updateGlycoProfiles(gaTargets);
				glyProCurr.writeProfile(PTMShepherd.normFName(ds + glycoProfileName));
			}
			glyProGLobal.writeProfile(PTMShepherd.normFName(globalName + glycoProfileName));
			System.out.println("Done with diagnostic ion extraction\n");
		}

		//Combine tables
		out.println("Combining and cleaning reports");
		CombinedTable gct = new CombinedTable(globalName);
		gct.writeCombinedTable(Integer.parseInt(params.get("histo_intensity")), calcIntensity);
		for (String ds : datasets.keySet()){
			out.println("Writing combined table for dataset " + ds);
			CombinedTable cct = new CombinedTable(ds);
			cct.writeCombinedTable(Integer.parseInt(params.get("histo_intensity")), calcIntensity);
		}

		//Glycan assignment
		boolean glycoMode = Boolean.parseBoolean(params.get("glyco_mode"));
		if (glycoMode) {
			System.out.println("Beginning glycan assignment");
			// parse glyco parameters and initialize database and ratio tables
			Random randomGenerator = new Random(glycoRandomSeed);
			ArrayList<GlycanResidue> adductList = StaticGlycoUtilities.parseGlycoAdductParam();
			int maxAdducts = Integer.parseInt(params.get("max_adducts"));
			String decoyParam = getParam("decoy_type");
			int decoyType = decoyParam.length() > 0 ? Integer.parseInt(decoyParam): GlycoAnalysis.DEFAULT_GLYCO_DECOY_TYPE;
			ProbabilityTables glycoProbabilityTable = StaticGlycoUtilities.initGlycoProbTable();
			HashMap<GlycanResidue, ArrayList<GlycanFragmentDescriptor>> glycoOxoniumDatabase = GlycoAnalysis.parseOxoniumDatabase(glycoProbabilityTable);
			double glycoPPMtol = getParam("glyco_ppm_tol").equals("") ? GlycoAnalysis.DEFAULT_GLYCO_PPM_TOL : Double.parseDouble(getParam("glyco_ppm_tol"));
			Integer[] glycoIsotopes = StaticGlycoUtilities.parseGlycoIsotopesParam();
			boolean nGlycan = getParam("n_glyco").equals("") || Boolean.parseBoolean(getParam("n_glyco"));		// default true
			glycoDatabase = StaticGlycoUtilities.parseGlycanDatabase(getParam("glycodatabase"), adductList, maxAdducts, randomGenerator, decoyType, glycoPPMtol, glycoIsotopes, glycoProbabilityTable, glycoOxoniumDatabase, nGlycan);
			String glycoMassFilePath = normFName(glycoMassListName);
			StaticGlycoUtilities.writeGlycanMassList(glycoDatabase, glycoMassFilePath);
			boolean glycoYnorm = getParam("norm_Ys").equals("") || Boolean.parseBoolean(getParam("norm_Ys"));		// default to True if not specified
			double absScoreErrorParam = getParam("glyco_abs_score_base").equals("") ? GlycoAnalysis.DEFAULT_GLYCO_ABS_SCORE_BASE : Double.parseDouble(getParam("glyco_abs_score_base"));
			String glycoFDRParam = getParam("glyco_fdr");
			double glycoFDR = glycoFDRParam.equals("") ? GlycoAnalysis.DEFAULT_GLYCO_FDR : Double.parseDouble(glycoFDRParam); 	// default 0.01 if param not provided, otherwise read provided value
			boolean alreadyPrintedParams = false;
			boolean printFullParams = !getParam("print_full_glyco_params").equals("") && Boolean.parseBoolean(getParam("print_full_glyco_params"));		// default false - for diagnostics
			boolean writeGlycansToAssignedMods = getParam("put_glycans_to_assigned_mods").equals("") || Boolean.parseBoolean(getParam("put_glycans_to_assigned_mods"));	// default true
			boolean removeGlycanDeltaMass = getParam("remove_glycan_delta_mass").equals("") || Boolean.parseBoolean(getParam("remove_glycan_delta_mass"));	// default true
			boolean printGlycoDecoys = !getParam("print_decoys").equals("") && Boolean.parseBoolean(getParam("print_decoys"));	// default false
			String allowedLocRes = getParam("localization_allowed_res");
			int numThreads = Integer.parseInt(params.get("threads"));
			boolean useGlycanFragmentProbs = !getParam("use_glycan_fragment_probs").equals("") && Boolean.parseBoolean(getParam("use_glycan_fragment_probs"));	// default false
			boolean useNewFDR = getParam("use_new_glycan_fdr").equals("") || Boolean.parseBoolean(getParam("use_new_glycan_fdr"));	// default true
			boolean useNonCompFDR = !getParam("use_noncomp_glycan_fdr").equals("") && Boolean.parseBoolean(getParam("use_noncomp_glycan_fdr"));	// default false
			double defaultProp = getParam("glyco_default_propensity").equals("") ? GlycoAnalysis.DEFAULT_GLYCO_PROPENSITY : Double.parseDouble(getParam("glyco_default_propensity"));

			// Glyco: first pass
			for (String ds : datasets.keySet()) {
				GlycoAnalysis ga = new GlycoAnalysis(ds, glycoDatabase, glycoProbabilityTable, glycoYnorm, absScoreErrorParam, glycoIsotopes, glycoPPMtol);
				if (ga.isGlycoComplete()) {
					print(String.format("\tGlyco analysis already done for dataset %s, skipping", ds));
					continue;
				}

				// print params here to avoid printing if the analysis is already done/not being run
				if (!alreadyPrintedParams) {
					StaticGlycoUtilities.printGlycoParams(adductList, maxAdducts, glycoDatabase, glycoProbabilityTable, glycoYnorm, absScoreErrorParam, glycoIsotopes, glycoPPMtol, glycoFDR, printFullParams, nGlycan, allowedLocRes, decoyType, printGlycoDecoys, removeGlycanDeltaMass);
					alreadyPrintedParams = true;
				}
				ArrayList<String[]> dsData = datasets.get(ds);
				for (int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					ga.glycoPSMs(pf, mzMap.get(ds), executorService, numThreads);
				}
				ga.completeGlyco();
			}

			// second pass: calculate glycan FDR and update results
			for (String ds : datasets.keySet()) {
				GlycoAnalysis ga = new GlycoAnalysis(ds, glycoDatabase, glycoProbabilityTable, glycoYnorm, absScoreErrorParam, glycoIsotopes, glycoPPMtol);
				if (!useNewFDR) {
					ga.computeGlycanFDROld(glycoFDR,true);
				} else {
					if (useNonCompFDR) {
						ga.useNonCompFDR = true;
						ga.computeGlycanFDR(glycoFDR);
					} else {
						boolean firstFDRsuccess = ga.computeGlycanFDROld(glycoFDR, false);
						if (!firstFDRsuccess) {
							ga.useNonCompFDR = true;
							ga.computeGlycanFDR(glycoFDR);
						}
					}
				}

				if (useGlycanFragmentProbs) {
					// second pass - calculate fragment propensities, regenerate database, and re-run
					HashMap<String, GlycanCandidateFragments> fragmentDB = ga.computeGlycanFragmentProbs();
					ArrayList<GlycanCandidate> propensityGlycanDB = StaticGlycoUtilities.updateGlycanDatabase(fragmentDB, glycoDatabase, randomGenerator, decoyType, glycoPPMtol, glycoIsotopes, glycoProbabilityTable, glycoOxoniumDatabase);

					// run glyco PSM-level analysis with the new database
					GlycoAnalysis ga2 = new GlycoAnalysis(ds, propensityGlycanDB, glycoProbabilityTable, glycoYnorm, absScoreErrorParam, glycoIsotopes, glycoPPMtol);
					ga2.glycanMassBinMap = ga.glycanMassBinMap;
					ga2.useFragmentSpecificProbs = true;
					ga2.defaultPropensity = defaultProp;
					ArrayList<String[]> dsData = datasets.get(ds);
					for (String[] dsDatum : dsData) {
						PSMFile pf = new PSMFile(new File(dsDatum[0]));
						ga2.glycoPSMs(pf, mzMap.get(ds), executorService, numThreads);
					}
					if (!useNewFDR) {
						ga2.computeGlycanFDROld(glycoFDR, true);
					} else {
						if (useNonCompFDR) {
							ga2.useNonCompFDR = true;
							ga2.computeGlycanFDR(glycoFDR);
						} else {
							boolean firstFDRsuccess2 = ga2.computeGlycanFDROld(glycoFDR, false);
							if (!firstFDRsuccess2) {
								ga2.useNonCompFDR = true;
								ga2.computeGlycanFDR(glycoFDR);
							}
						}
					}
					ga2.completeGlyco();
				}
			}

			/* Save best glycan information from glyco report to psm tables */
			for (String ds : datasets.keySet()) {
				ArrayList<String[]> dsData = datasets.get(ds);
				for (int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					pf.mergeGlycoTable(new File(normFName(ds + rawGlycoName)), GlycoAnalysis.NUM_ADDED_GLYCO_PSM_COLUMNS, writeGlycansToAssignedMods, nGlycan, allowedLocRes, removeGlycanDeltaMass, printGlycoDecoys);
				}
			}

			print("Created glyco reports");
			print("Done with glycan assignment\n");
		}

		/* Make psm table IonQuant compatible */
		if (Boolean.parseBoolean(params.get("prep_for_ionquant"))) {
			out.println("Prepping PSM tables for IonQuant");
			for (String ds : datasets.keySet()) {
				ArrayList<String[]> dsData = datasets.get(ds);
				for (int i = 0; i < dsData.size(); i++) {
					PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
					pf.preparePsmTableForIonQuant(peakBounds, Integer.parseInt(params.get("precursor_mass_units")), Double.parseDouble(params.get("precursor_tol")));
				}
			}
			out.println("Done");
		}

		/* Make experiment-level table */
		if (Boolean.parseBoolean(params.get("output_extended"))) {
			out.println("Creating experiment-level profile report");
			CombinedExperimentsSummary cs = new CombinedExperimentsSummary(normFName("combined_experiment_profile.tsv"));
			cs.initializeExperimentSummary(normFName(globalName + profileName), Integer.parseInt(params.get("histo_intensity")));
			cs.addLocalizationSummary(normFName(globalName + locProfileName), combinedName);
			cs.addSimilarityRTSummary(normFName(globalName + simRTProfileName), combinedName, calcIntensity);
			for (String ds : datasets.keySet()) {
				cs.addExperimentSummary(normFName(ds + profileName), ds);
				cs.addLocalizationSummary(normFName(ds + locProfileName), ds);
				cs.addSimilarityRTSummary(normFName(ds + simRTProfileName), ds, calcIntensity);
			}
			cs.printFile();
		}

		List<String> filesToDelete = Arrays.asList(normFName(peaksName), normFName(peakSummaryAnnotatedName),
				normFName(peakSummaryName), normFName(combinedTSVName), normFName(combinedHistoName));
		//delete redundant files
		for (String f : filesToDelete) {
			Path p = Paths.get(f).toAbsolutePath().normalize();
			deleteFile(p, true);
		}

		List<String> extsToDelete = Arrays.asList(locProfileName, simRTProfileName, profileName,
				histoName, ms2countsName);
		for (String ext : extsToDelete) {
			Path p = Paths.get(normFName(globalName + ext)).toAbsolutePath().normalize();
			if (!ext.equals(profileName))
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
					.asList(rawLocalizeName, rawSimRTName, rawGlycoName, histoName, diagBinFilename, mzBinFilename);
			for (String ds : datasets.keySet()) {
				//System.out.println("Writing combined table for dataset " + ds);
				//CombinedTable.writeCombinedTable(ds);
				for (String ext : extsToDelete) {
					Path p = Paths.get(normFName(ds + ext)).toAbsolutePath().normalize();
					deleteFile(p, true);
				}
				for (String mzFile : mzMap.get(ds).keySet()) {
					Path p = Paths.get(normFName(mzFile + mzBinFilename)).toAbsolutePath().normalize();
					deleteFile(p, true);
				}
			}
			String dsTmp = combinedName;
			for (String ext : extsToDelete) {
				if (ext.equals(glycoProfileName))
					continue;
				Path p = Paths.get(normFName(dsTmp + ext)).toAbsolutePath().normalize();
				deleteFile(p, true);
			}
			dsTmp = globalName;
			for (String ext : extsToDelete) {
				if (ext.equals(glycoProfileName))
					continue;
				Path p = Paths.get(normFName(dsTmp + ext)).toAbsolutePath().normalize();
				deleteFile(p, true);
			}

		}
		// delete diagbin files even if output extended is on (hidden parameter to keep them due to processing time)
		if (!Boolean.parseBoolean(params.get("output_extended_diagmine"))) {
			for (String ds : datasets.keySet()) {
				for (String mzFile : mzMap.get(ds).keySet()) {
					Path p = Paths.get(normFName(mzFile + diagBinFilename)).toAbsolutePath().normalize();
					deleteFile(p, true);
				}
			}
		}

		for (String crc : cacheFiles) {
			File crcf = new File(normFName("cache-"+crc +".txt"));
			Path crcpath = crcf.toPath().toAbsolutePath();
			deleteFile(crcpath, true);
		}

		for (String ds : datasets.keySet()) {
			for (String crun : mzMap.get(ds).keySet()) {
				File mzbf = mzMap.get(ds).get(crun);
				Path mzbpath = mzbf.toPath().toAbsolutePath();
				deleteFile(mzbpath, true);
			}
		}

		executorService.shutdown();
	}

	private static void deleteFile(Path p, boolean printOnDeletion) throws IOException {
		if (Files.deleteIfExists(p) && printOnDeletion) {
			print("Deleted file: " + p.toAbsolutePath().normalize().toString());
		}
	}

	private static void getMzDataMapping() throws Exception {
		cacheFiles = new ArrayList<>();

		// Get true paths to mzData
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
				PTMShepherd.print("\tIndexing data from " + ds);
				PSMFile pf = new PSMFile(dsData.get(i)[0]);
				PSMFile.getMappings(new File(dsData.get(i)[1]), mzMap.get(ds), pf.getRunNames());
			}
			// Assure that mzData was found
			for(String crun : mzMap.get(ds).keySet()) {
				if(mzMap.get(ds).get(crun) == null) {
					die("In dataset \""+ds+"\" could not find mzData for run " +  crun);
				}
			}
		}
	}

	// Rewrite mzData to MZBIN files
	private static void rewriteDataToMzBin() throws Exception {
		for(String ds : datasets.keySet()) {
			ArrayList<String []> dsData = datasets.get(ds);
			for(int i = 0; i < dsData.size(); i++) {
				PTMShepherd.print("\tCaching data from " + ds);
				PSMFile pf = new PSMFile(dsData.get(i)[0]);
				rewriteMzDataToMzBin(pf, mzMap.get(ds), Integer.parseInt(params.get("spectra_condPeaks")), Float.parseFloat(params.get("spectra_condRatio")));
				PTMShepherd.print("\tDone caching data from " + ds);
			}
		}
	}

	private static void rewriteMzDataToMzBin(PSMFile pf, HashMap<String, File> mzMappings, int topNPeaks, float minPeakRatio) throws Exception {
		// Get PSM scan num -> spectral file mapping
		HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
		int specCol = pf.getColumn("Spectrum");
		for (int i = 0; i < pf.data.size(); i++) {
			String[] sp = pf.data.get(i).split("\t");
			String bn = sp[specCol].substring(0, sp[specCol].indexOf(".")); //fraction
			if (!mappings.containsKey(bn))
				mappings.put(bn, new ArrayList<>());
			mappings.get(bn).add(i);
		}

		// Loop through spectral files -> indexed lines in PSM -> process each line
		HashMap<Integer, String> linesWithoutSpectra = new HashMap<>();
		for (String cf : mappings.keySet()) { //for file in relevant spectral files
			// Check to see if an mzBIN_cache file was discovered
			if (mzMappings.get(cf).toString().endsWith(mzBinFilename)) {
				System.out.println("\t\tFound existing cached spectral data for " + cf);
				continue;
			}
			// Check to see if mzBIN_cache file exists in output directory
			else if (new File(normFName(cf + mzBinFilename)).exists()) {
				System.out.println("\t\tFound existing cached spectral data for " + cf);
				mzMappings.put(cf, new File(normFName(cf + mzBinFilename)));
				continue;
			}
			long t1 = System.currentTimeMillis();
			MXMLReader mr = new MXMLReader(mzMappings.get(cf), Integer.parseInt(PTMShepherd.getParam("threads")));
			ArrayList<MZBINSpectrum> specs = new ArrayList<>(); // Holds parsed spectra
			mr.readFully();
			long t2 = System.currentTimeMillis();
			ArrayList<Integer> clines = mappings.get(cf); //lines corr to curr spec file
			for (int i = 0; i < clines.size(); i++) {//for relevant line in curr spec file
				String line = pf.data.get(clines.get(i));
				String [] sp = line.split("\\t");
				String specName = sp[specCol];
				Spectrum spec =  mr.getSpectrum(reNormName(specName));
				if (spec == null)
					linesWithoutSpectra.put(i, line);
				else {
					spec.condition(topNPeaks, minPeakRatio); // TODO Why aren't these being saved as conditioned spectra?
					specs.add(spec.toMZBINSpectrum());
				}
			}
			long t3 = System.currentTimeMillis();
			PTMShepherd.print(String.format("\t\t%s - %d (%d ms, %d ms)", mzMappings.get(cf), clines.size(), t2-t1,t3-t2));
			MZBINFile mzbinFile = new MZBINFile(normFName(cf + mzBinFilename), specs, "", "");
			mzbinFile.writeMZBIN();
			mzMappings.put(cf, new File(normFName(cf + mzBinFilename)));
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
		return outputPath + outputDirName + fpath;
	}

	/* This method will move main files out of subdirectory and into main directory */
	public static void moveFileFromExtendedDir(String oldFpath) {
		String newFpath = oldFpath.replace(outputDirName, "");
		File oldF = new File(oldFpath);
		File newF = new File(newFpath);

		if (newF.exists()) {
			try {
				deleteFile(newF.toPath(), true);
				out.printf("Moving %s\n", oldF);
			} catch (Exception e) {
				out.printf("Error moving %s\n", oldF);
			}
		}
		try {
			oldF.renameTo(newF);
		} catch (Exception e) {
			out.printf("Error deleting old %s\n", oldF);
		}
	}

	/* Makes output directory */
	public static void makeOutputDir(String dpath) throws Exception {
		try {
			if (!dpath.equals("")) {
				File dir = new File(dpath);
				if (!dir.exists())
					dir.mkdirs();

				File eoDir = new File(dpath + outputDirName);
				if (!eoDir.exists())
					eoDir.mkdirs();

			}
		} catch (Exception e) {
			out.println("Error creating output directory. Terminating.");
			System.exit(1);
		}
	}

	/**
	 * Removes charge state from spec name to match spectra that have had charge state reannotated.
	 * If the spec name is already renormed, returns input string
	 * @param s
	 * @return specName without charge state
	 */
	public static String reNormName(String s) {
		String[] sp = s.split("\\.");
		if (sp[sp.length-2].equals(sp[sp.length-1])) // If the format is XXX.Spec.Spec already, return input string
			return s;
		int sn = Integer.parseInt(sp[1]);
		//without charge state
		return String.format("%s.%d.%d", sp[0], sn, sn);
	}

	public static String reNormNameWithCharge(String s) {
		String[] sp = s.split("\\.");
		int sn = Integer.parseInt(sp[1]);
		//with charge state
		return String.format("%s.%d.%d.%s",sp[0],sn,sn,sp[3]);
	}

	public static String concatIonTypes() {
		StringBuffer sb = new StringBuffer();
		if (params.get("iontype_a").equals("1"))
			sb.append("a");
		if (params.get("iontype_b").equals("1"))
			sb.append("b");
		if (params.get("iontype_c").equals("1"))
			sb.append("c");
		if (params.get("iontype_x").equals("1"))
			sb.append("x");
		if (params.get("iontype_y").equals("1"))
			sb.append("y");
		if (params.get("iontype_z").equals("1"))
			sb.append("z");
		return sb.toString();
	}
}
