package edu.umich.andykong.ptmshepherd.core;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
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
	
	int nSpecs;
	Spectrum [] specs;
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
	
	public MXMLReader(File f) {
		this.f = f;
	}
	
	public void readFully() throws Exception {

		String fn = f.toPath().getFileName().toString().toLowerCase();
		LCMSDataSource<?> source = null;
		if(fn.endsWith(".mzxml")) {
			source = new MZXMLFile(f.getAbsolutePath());
		} else if(fn.endsWith(".mzml")) {
			source = new MZMLFile(f.getAbsolutePath());
		} else if (fn.endsWith(".raw")) {
			source = new ThermoRawFile(f.getAbsolutePath());
			/*
		} else if (fn.endsWith(".d")) {
			String tempFn = fn.replaceFirst(".d$", ".mzBIN");
			File tempF = new File(tempFn);
			if (tempF.exists()) {
				String f = tempFn;
				source = new MZBINFile(tempF.getAbsolutePath());
			} else {
				System.out.println("Cannot read .d files without associated .mzBIN");
				System.exit(1);
			}
			*/

		}
		if (source == null) {
			System.out.println("Cannot read mzFile with unrecognized extension: " + f.getName());
			System.exit(1);
		}

		readFully(source);
		specsByName = new HashMap<>();
		specsByStrippedName = new HashMap<>();
		for(int i = 0; i < specs.length; i++) {
			specsByName.put(specs[i].scanName, specs[i]);
			specsByStrippedName.put(specs[i].scanName, specs[i]);
		}
	}

	public void readFully(LCMSDataSource<?> source) throws Exception {
		String baseName = f.getName();
		baseName = baseName.substring(0,baseName.lastIndexOf("."));

		source.setExcludeEmptyScans(false);
		source.setNumThreadsForParsing(Math.min(8,Runtime.getRuntime().availableProcessors()));

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
			ns.scanName = baseName+"."+scanNumRaw+"."+scanNumRaw+"."+ns.charge;
			ns.precursorMass = scan.getPrecursor().getMzTarget();
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
	}
	
	public static void main(String [] args) throws Exception {
		MXMLReader mr = new MXMLReader(new File("E:\\q01507.mzXML"));
		mr.readFully();
		for(int i = 0; i < mr.specs.length; i++)
			System.out.println(mr.specs[i].scanName);
	}
	
}
