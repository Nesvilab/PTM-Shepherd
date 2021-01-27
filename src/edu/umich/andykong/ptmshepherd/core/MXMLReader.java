package edu.umich.andykong.ptmshepherd.core;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.filetypes.LCMSDataSource;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;
import umich.ms.fileio.filetypes.thermo.ThermoRawFile;


public class MXMLReader {

	File f;
	final int threads;
	
	int nSpecs;
	public Spectrum [] specs;
	HashMap<String,Spectrum> specsByName;
	HashMap<String,Spectrum> specsByStrippedName;

	public Spectrum getSpectrum(String specName) {
		return specsByStrippedName.get(specName);
//		return specsByName.get(specName);
	}
	
	public String stripName(String specName) {
		String [] sp = specName.split(".");
		if(sp.length == 4) 
			return sp[0]+"."+sp[1];
		else
			return specName;
	}
	
	public MXMLReader(File f, int threads) {
		this.f = f;
		this.threads = threads;
	}
	
	public void readFully() throws Exception {
		String fn = f.toPath().getFileName().toString().toLowerCase();
		LCMSDataSource<?> source = null;
		MZBINFile mzbinSource = null;
		MSFMGFFile mgfSource = null;

		if (fn.endsWith("_calibrated.mgf")) {
			mgfSource = new MSFMGFFile(PTMShepherd.executorService, Integer.parseInt(PTMShepherd.getParam("threads")), f, true);
		} else if (fn.endsWith(".mzxml")) {
			source = new MZXMLFile(f.getAbsolutePath());
		} else if (fn.endsWith(".mzml")) {
			source = new MZMLFile(f.getAbsolutePath());
		} else if (fn.endsWith(".raw")) {
			source = new ThermoRawFile(f.getAbsolutePath());
		} else if (fn.endsWith(".mzbin")) {
			mzbinSource = new MZBINFile(PTMShepherd.executorService, Integer.parseInt(PTMShepherd.getParam("threads")), f, true);
		}
		if ((mzbinSource == null) && (mgfSource == null) && (source == null)) {
			System.out.println("Cannot read mzFile with unrecognized extension: " + f.getName());
			System.exit(1);
		}

		specsByName = new HashMap<>();
		specsByStrippedName = new HashMap<>();
		if (mzbinSource == null && mgfSource == null) { //if filetype is not mzBin
			readFully(source);
			for(int i = 0; i < specs.length; i++) {
				specsByName.put(specs[i].scanName, specs[i]);
				specsByStrippedName.put(specs[i].scanName, specs[i]);
			}
		} else if (mgfSource == null) {
			readAsMzBIN(mzbinSource);
			for(int i = 0; i < specs.length; i++){
				specsByName.put(specs[i].scanName, specs[i]);
				specsByStrippedName.put(specs[i].scanName, specs[i]);
			}
		} else {
			readAsFraggerMGF(mgfSource);
			for(int i = 0; i < specs.length; i++) {
				specsByName.put(specs[i].scanName, specs[i]);
				specsByStrippedName.put(specs[i].scanName, specs[i]);
			}
		}
	}

	//TN_CSF_062617_32.46509.46509
	public void readFully(LCMSDataSource<?> source) throws Exception {
		String baseName = f.getName();
		baseName = baseName.substring(0,baseName.lastIndexOf("."));

		source.setExcludeEmptyScans(false);
		source.setNumThreadsForParsing(threads);

		ScanCollectionDefault scans = new ScanCollectionDefault();
		scans.setDataSource(source);
		scans.loadData(LCMSDataSubset.MS2_WITH_SPECTRA);
		TreeMap<Integer,IScan> num2scan = scans.getMapNum2scan();
		Set<Map.Entry<Integer, IScan>> scanEntries = num2scan.entrySet();

		ArrayList<Spectrum> cspecs = new ArrayList<Spectrum>();
		for (Map.Entry<Integer, IScan> scanEntry : scanEntries) {
			Integer scanNum = scanEntry.getKey();
			IScan scan = scanEntry.getValue();
			if(scan.getMsLevel() != 2)
				continue;
			int scanNumRaw = scanNum;
			ISpectrum spectrum = scan.fetchSpectrum();
			int clen = (spectrum == null)?0:spectrum.getMZs().length;
			Spectrum ns = new Spectrum(clen);
			try {
				ns.charge = (byte)scan.getPrecursor().getCharge().intValue();
			} catch(Exception e) {
				ns.charge = 0;
			}
			ns.scanNum = scanNumRaw;
			ns.rt = scan.getRt();
			//includes charge state in scan num
			//ns.scanName = baseName+"."+scanNumRaw+"."+scanNumRaw+"."+ns.charge;
			//does not include charge state in scan num
			ns.scanName = baseName+"."+scanNumRaw+"."+scanNumRaw;
			if (scan.getPrecursor().getMzTargetMono() != null) {
				ns.precursorMass = scan.getPrecursor().getMzTargetMono();
			} else if (scan.getPrecursor().getMzTarget() != null) {
				ns.precursorMass = scan.getPrecursor().getMzTarget();
			} else {
				throw new IllegalStateException("No precursor mz information found");
			}
			if (scan.getPrecursor().getMzTargetMono() != null)
				ns.monoMass = scan.getPrecursor().getMzTargetMono();
			if (scan.getPrecursor().getMzTarget() != null)
				ns.targetMass = scan.getPrecursor().getMzTarget();
			for(int i = 0; i < clen; i++) {
				ns.peakMZ[i] = (float)spectrum.getMZs()[i];
				ns.peakInt[i] = (float)spectrum.getIntensities()[i];
			}
			cspecs.add(ns);
		}
		nSpecs = cspecs.size();
		specs = new Spectrum[nSpecs];
		cspecs.toArray(specs);
		scans.reset();
		if (source instanceof ThermoRawFile) {
			((ThermoRawFile) source).close();
		};
	}

	//400ngHeLaosmoothCE20-52lowguessSRIG450easy4_30tbl1_0NOexp12scansi_A1_01_3366.109793.109793.2
	private void readAsMzBIN(MZBINFile mf) throws Exception {
		ArrayList cspecs = new ArrayList<Spectrum>();
		for (Spectrum spectrum : mf.specs) {
			if (spectrum.msLevel == 2)
				cspecs.add(spectrum);
		}
		nSpecs =  cspecs.size();
		specs = new Spectrum[nSpecs];
		cspecs.toArray(specs);
	}

	private void readAsFraggerMGF(MSFMGFFile mf) throws Exception {
		ArrayList cspecs = new ArrayList<Spectrum>();
		for (Spectrum spectrum : mf.specs) {
			if (spectrum.msLevel == 2)
				cspecs.add(spectrum);
		}
		nSpecs = cspecs.size();
		specs = new Spectrum[nSpecs];
		cspecs.toArray(specs);
	}

	/*
	private void readAsMzBIN(ExecutorService executorService) throws Exception {
		ArrayList scans = new ArrayList<>();
		MZBINFile mf = new MZBINFile(executorService, threads, f, true);
		for (Spectrum spectrum : mf.specs) {
			if (spectrum.msLevel == 2 && (sp.excludedScanSet.isEmpty() || !sp.excludedScanSet.contains(spectrum.scanName)))
				scans.add(spectrum);
		}
		results = new SpecResults[scans.size()];
		nScans = results.length;
	}
	*/

	public static void main(String [] args) throws Exception {
		MXMLReader mr = new MXMLReader(new File("E:\\q01507.mzXML"), 7);
		mr.readFully();
		//for(int i = 0; i < mr.specs.length; i++)
		//	System.out.println(mr.specs[i].scanName);
	}
	
}
