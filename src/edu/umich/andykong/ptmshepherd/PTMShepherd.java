package edu.umich.andykong.ptmshepherd;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.io.*;

import edu.umich.andykong.ptmshepherd.cleaner.CombinedTable;
import edu.umich.andykong.ptmshepherd.localization.*;
import edu.umich.andykong.ptmshepherd.peakpicker.*;
import edu.umich.andykong.ptmshepherd.specsimilarity.*;

public class PTMShepherd {

	public static final String name = "PTM-Shepherd";
 	public static final String version = "0.2.5";

	static HashMap<String,String> params;
	static TreeMap<String,ArrayList<String []>> datasets;
	static HashMap<String,HashMap<String,File>> mzMap;
	static HashMap<String,Integer> datasetMS2;
	
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
		File f = new File(fn);
		if(!f.exists()) 
			die(String.format("Parameter file %s does not exist",fn));
		BufferedReader in = new BufferedReader(new FileReader(f));
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
				params.put(key, value);
			}
		}
		in.close();
	}
	
	public static void init(String [] args) throws Exception {
		HashMap<String,String> overrides = new HashMap<>();
		params = new HashMap<>();
		datasets = new TreeMap<>();
		mzMap = new HashMap<>();
		datasetMS2 = new HashMap<>();
		
		//default values
		params.put("threads", ""+Math.min(8, Runtime.getRuntime().availableProcessors()));
		params.put("histo_bindivs", "5000"); //number of divisions in histogram
		params.put("histo_smoothbins", "5"); //smoothing factor
		
		params.put("peakpicking_promRatio", "0.3"); //prominence ratio for peakpicking
		params.put("peakpicking_width", "0.005"); //width for peakpicking
		//params.put("peakpicking_background", "0.005");
		params.put("peakpicking_topN", "500"); //num peaks
        params.put("peakpicking_minPsm", "10");
        params.put("localization_background", "4");
        params.put("varmod_masses", ":0");
		params.put("precursor_tol", "0.01"); //unimod peakpicking width
		//params.put("precursor_tol_ppm", "20.0"); //unused
		
		params.put("spectra_ppmtol", "20.0"); //obvious, used in localization and simrt
		params.put("spectra_condPeaks", "100"); //
		params.put("spectra_condRatio", "0.01"); //

		params.put("output_extended", "false");
		//params.put("output_extended", "true");
		
		//load parameters
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals("--")) {
				overrides.put(args[i].substring(2),args[i+1]);
				i++;
			} else
				parseParamFile(args[i]);
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
	}
	
	public static void main(String [] args) throws Exception {
		Locale.setDefault(new Locale("en","US"));
		System.out.printf("%s version %s",name,version);
		System.out.println("(c) University of Michigan\n");
		System.out.printf("Using Java %s on %dMB memory\n\n", System.getProperty("java.version"),(int)(Runtime.getRuntime().maxMemory()/Math.pow(2, 20)));
		
		if(args.length == 0) {
			System.out.printf("Usage: %s [config.txt]\n",name);
			System.exit(1);
		}
		
		init(args);
		
		//Get mzData mapping
		ArrayList<String> cacheFiles = new ArrayList<>();
		for(String ds : datasets.keySet()) {
			ArrayList<String []> dsData = datasets.get(ds);
			mzMap.put(ds, new HashMap<>());
			for(int i = 0; i < dsData.size(); i++) {
				File tpf = new File(dsData.get(i)[0]);
				String crc = PSMFile.getCRC32(tpf);
				File cacheFile = new File("cache-"+crc+".txt");
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
			}
			for(String crun : mzMap.get(ds).keySet()) {
				if(mzMap.get(ds).get(crun) == null) {
					die("In dataset \""+ds+"\" could not find mzData for run " +  crun);
				}
			}
		}
		
		//Count MS2 scans
		for(String ds : datasets.keySet()) {
			File countsFile = new File(ds+".ms2counts");
			int sumMS2 = 0;
			TreeMap<String,Integer> counts = new TreeMap<>();
			if(!countsFile.exists()) {
				print("Counting MS2 scans for dataset " + ds);
				for(String crun : mzMap.get(ds).keySet()) {
					File tf = mzMap.get(ds).get(crun);
					int cnt = MS2Counts.countMS2Scans(tf);
					print(String.format("\t%s - %d scans", crun,cnt));
					counts.put(crun, cnt);
				}
				PrintWriter out = new PrintWriter(new FileWriter(countsFile));
				for(String cf : counts.keySet()) {
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
			print(datasetMS2.get(ds) +" MS2 scans present in dataset " + ds + "\n");
		}
		
		//Generate histograms
		File combinedHisto = new File("combined.histo");
		if(!combinedHisto.exists()) {
			print("\nCreating combined histogram");
			int min = 1 << 30;
			int max = -1*(1<<30);
			
			for(String ds : datasets.keySet()) {
				File histoFile = new File(ds+".histo");
				if(!histoFile.exists()) {
					ArrayList<String []> dsData = datasets.get(ds);
					ArrayList<Float> vals = new ArrayList<>();
					for(int i = 0; i < dsData.size(); i++) {
						PSMFile pf = new PSMFile(new File(dsData.get(i)[0]));
						vals.addAll(pf.getMassDiffs());
					}
					Histogram chisto = new Histogram(vals, datasetMS2.get(ds), Integer.parseInt(params.get("histo_bindivs")),Integer.parseInt(params.get("histo_smoothbins")));
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
				File histoFile = new File(ds+".histo");
				Histogram h = Histogram.readHistogram(histoFile);
				combined.mergeHistogram(h, ds);
			}
			combined.writeHistogram(combinedHisto);
			combined.writeCombinedTSV(new File("combined.tsv"));
			print("Created combined histogram!\n");
		} else
			print("Combined histogram found\n");
		
		//Perform peak detection
		File peaks = new File("peaks.tsv");
		if (peaks.exists()) {
			print("Deleting old peaks.tsv file at: " + peaks.toString() + "\n");
		}
		print("Running peak picking\n");
		Histogram combined = Histogram.readHistogram(combinedHisto);
		PeakPicker pp = new PeakPicker();
		pp.pickPeaks(combined.getOffsets(), combined.histo, Double.parseDouble(params.get("peakpicking_promRatio")),
				Double.parseDouble(params.get("peakpicking_width")),Double.parseDouble(params.get("peakpicking_background")),Integer.parseInt(params.get("peakpicking_topN")));
		pp.writeTSV(new File("peaks.tsv"));
		print("Picked top " + Integer.parseInt(params.get("peakpicking_topN")) + " peaks\n");
		
		//PSM assignment
		File peaksummary = new File("peaksummary.tsv");
		if(!peaksummary.exists()) {
			PeakSummary ps = new PeakSummary(peaks,Double.parseDouble(params.get("precursor_tol")));
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
			print("created summary table\n");
		}

		//Assign peak IDs
		File peakannotated = new File("peaksummary.annotated.tsv");
		if(!peakannotated.exists()) {
			PeakAnnotator pa = new PeakAnnotator();
			pa.init(params.get("varmod_masses"));
			pa.annotateTSV(peaksummary, peakannotated);
			print("annotated summary table\n");
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
		print("done\n");
		
		//Localization summaries
		double [] peakCenters = PeakSummary.readPeakCenters(peaksummary);
		LocalizationProfile loc_global = new LocalizationProfile(peakCenters, Double.parseDouble(params.get("precursor_tol")));
		for(String ds : datasets.keySet()) {
			LocalizationProfile loc_current = new LocalizationProfile(peakCenters, Double.parseDouble(params.get("precursor_tol")));
			SiteLocalization sl = new SiteLocalization(ds);
			LocalizationProfile [] loc_targets = {loc_global, loc_current};
			sl.updateLocalizationProfiles(loc_targets); //this is where the localization is actually happening
			loc_current.writeProfile(ds+".locprofile.txt");
		} 
		loc_global.writeProfile("global.locprofile.txt");
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
				sra.simrtPSMs(pf, mzMap.get(ds));
			}
			sra.complete();
		}		
		print("Done\n");
		
		//SimRT summaries
		SimRTProfile simrt_global = new SimRTProfile(peakCenters, Double.parseDouble(params.get("precursor_tol")));
		for(String ds : datasets.keySet()) {
			SimRTProfile simrt_current = new SimRTProfile(peakCenters, Double.parseDouble(params.get("precursor_tol")));
			SimRTAnalysis sra = new SimRTAnalysis(ds);
			SimRTProfile [] simrt_targets = {simrt_global, simrt_current};
			sra.updateSimRTProfiles(simrt_targets);
			simrt_current.writeProfile(ds+".simrtprofile.txt");
		} 
		simrt_global.writeProfile("global.simrtprofile.txt");
		print("Created similarity/RT reports\n");

		//Combine tables
		System.out.println("Combining and cleaning reports");
		CombinedTable.writeCombinedTable("global");
		for (String ds : datasets.keySet()){
			System.out.println("Writing combined table for dataset " + ds);
			CombinedTable.writeCombinedTable(ds);
		}

		List<String> filesToDelete = Arrays.asList("peaks.tsv", "peaksummary.annotated.tsv", "peaksummary.tsv", "combined.tsv");

		for (String f : filesToDelete) {
			Path p = Paths.get(f).toAbsolutePath().normalize();
			deleteFile(p, true);
		}

		if(!Boolean.parseBoolean(params.get("output_extended"))) {
			// delete dataset files with specific extensions
			List<String> extsToDelete = Arrays
					.asList(".histo", ".locprofile.txt", ".ms2counts", ".simrtprofile.txt", ".rawlocalize", ".rawsimrt");
			for (String ds : datasets.keySet()) {
				//System.out.println("Writing combined table for dataset " + ds);
				//CombinedTable.writeCombinedTable(ds);
				for (String ext : extsToDelete) {
					Path p = Paths.get(ds + ext).toAbsolutePath().normalize();
					deleteFile(p, true);
				}
			}
			String dsTmp = "combined";
			for (String ext : extsToDelete) {
				Path p = Paths.get(dsTmp + ext).toAbsolutePath().normalize();
				deleteFile(p, true);
			}
			dsTmp = "global";
			for (String ext : extsToDelete) {
				Path p = Paths.get(dsTmp + ext).toAbsolutePath().normalize();
				deleteFile(p, true);
			}

		}

		for (String crc : cacheFiles) {
			File crcf = new File("cache-"+crc +".txt");
			Path crcpath = crcf.toPath().toAbsolutePath();
			deleteFile(crcpath, true);
		}

	}

	private static void deleteFile(Path p, boolean printOnDeletion) throws IOException {
		if (Files.deleteIfExists(p) && printOnDeletion) {
			print("Deleted file: " + p.toAbsolutePath().normalize().toString());
		}
	}

}
