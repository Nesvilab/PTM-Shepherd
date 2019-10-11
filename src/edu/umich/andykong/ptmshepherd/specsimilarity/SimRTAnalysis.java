package edu.umich.andykong.ptmshepherd.specsimilarity;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.localization.LocalizationProfile;
import java.util.List;

public class SimRTAnalysis {

	String dsName;
	File simRTFile;
	MXMLReader mr;
	double ppmTol, condRatio, peakTol;
	int condPeaks;
	int specCol, pepCol, modpepCol, chargeCol, deltaCol, rtCol, intCol;
	
	static final int MAX_ZERO_COMPARE = 20;
	
	public SimRTAnalysis(String dsName) {
		this.dsName = dsName;
		simRTFile = new File(dsName+".rawsimrt");
	}
	
	public String reNormName(String s) {
		String [] sp = s.split("\\.");
		int sn = Integer.parseInt(sp[1]);
		return String.format("%s.%d.%d.%s",sp[0],sn,sn,sp[3]);
	}
	
	public boolean isComplete() throws Exception {
		if(simRTFile.exists()) {
			RandomAccessFile raf = new RandomAccessFile(simRTFile, "r");
			raf.seek(Math.max(0, simRTFile.length() - 20));
			String cline;
			while((cline = raf.readLine())!=null)
				if(cline.equals("COMPLETE")) {
					raf.close();
					return true;
				}
			raf.close();
			simRTFile.delete();
		}
		return false;
	}
	
	public void complete() throws Exception {
		PrintWriter out = new PrintWriter(new FileWriter(simRTFile,true));
		out.println("COMPLETE");
		out.close();
	}
	
	public void simrtPSMs(PSMFile pf, HashMap<String,File> mzMappings) throws Exception {
		//assemble PSMs into per file groupings
		HashMap<String,ArrayList<Integer>> mappings = new HashMap<>();
		PrintWriter out = new PrintWriter(new FileWriter(simRTFile,true));
		
		specCol = pf.getColumn("Spectrum");
		pepCol = pf.getColumn("Peptide");
		modpepCol = pf.getColumn("Modified Peptide");
		chargeCol = pf.getColumn("Charge");
		deltaCol = pf.getColumn("Adjusted Delta Mass");
		rtCol = pf.getColumn("Retention");
		intCol = pf.getColumn("Intensity");
		
		ppmTol = Double.parseDouble(PTMShepherd.getParam("spectra_ppmtol"));
		condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
		condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));
		peakTol = Double.parseDouble(PTMShepherd.getParam("precursor_tol"));
		
		for(int i = 0; i < pf.data.size(); i++) {
			String [] sp = pf.data.get(i).split("\t");
			String bn = sp[specCol].substring(0,sp[specCol].indexOf("."));
			if(!mappings.containsKey(bn))
				mappings.put(bn, new ArrayList<>());
			mappings.get(bn).add(i);
		}

		//iterate and calculate similarity/retention time deltas for each file
		for(String cf : mappings.keySet()) {
			long t1 = System.currentTimeMillis();
			mr = new MXMLReader(mzMappings.get(cf));
			mr.readFully();
			long t2 = System.currentTimeMillis();
			ArrayList<Integer> clines = mappings.get(cf);
			HashMap<String,ArrayList<Integer>> zTolLines = new HashMap<>();
			HashMap<String,ArrayList<Spectrum>> zTolSpecs = new HashMap<>();
			HashMap<String,ArrayList<Double>> zTolRT = new HashMap<>();
			HashMap<String,Double> avgzSim = new HashMap<>();
			HashMap<String,Double> avgzRT = new HashMap<>();
			
			//get zero bin data
			for(int i = 0; i < clines.size(); i++) {
				String [] crow = pf.data.get(clines.get(i)).split("\t");
				boolean isZero = (Math.abs(Double.parseDouble(crow[deltaCol])) <= peakTol);
				if(!isZero)
					continue;
				
				String key = crow[pepCol].trim();
				if(crow[modpepCol].trim().length() != 0)
					key = crow[modpepCol].trim();
				
				if(!zTolRT.containsKey(key)) 
					zTolRT.put(key, new ArrayList<>());
				zTolRT.get(key).add(Double.parseDouble(crow[rtCol]));
				
				key += "." + crow[chargeCol];
				if(!zTolLines.containsKey(key))
					zTolLines.put(key, new ArrayList<>());
				zTolLines.get(key).add(clines.get(i));
			}
			
			//calculate zeroSim
			int totalLines = 0;
			List<String> linesWithoutSpectra = new ArrayList<>();
			for(String pepZ : zTolLines.keySet()) {
				ArrayList<Integer> relLines = zTolLines.get(pepZ);
				Collections.shuffle(zTolLines.get(pepZ));
				int nComp = Math.min(relLines.size(), MAX_ZERO_COMPARE);
				zTolSpecs.put(pepZ, new ArrayList<>());
				for(int i = 0; i < nComp; i++) {
					String [] crow = pf.data.get(relLines.get(i)).split("\t");
					zTolSpecs.get(pepZ).add(mr.getSpectrum(reNormName(crow[specCol])));
				}
				
				double zSimSum = 0;
				totalLines += relLines.size();
				for(int i = 0; i < relLines.size(); i++) {
					String [] crow = pf.data.get(relLines.get(i)).split("\t");
					String specNormName = reNormName(crow[specCol]);
					Spectrum cspec = mr.getSpectrum(specNormName);
					if (cspec == null) {
						linesWithoutSpectra.add(specNormName);
						continue;
					}
					zSimSum += cspec.averageSimilarity(zTolSpecs.get(pepZ), ppmTol);
				}
				avgzSim.put(pepZ, zSimSum / relLines.size());
			}

			if (!linesWithoutSpectra.isEmpty()) {
				System.out.printf("Could not find %d/%d (%.1f%%) spectra.\n", linesWithoutSpectra.size(), totalLines,
						100.0*((double)linesWithoutSpectra.size()/totalLines));
				int previewSize = Math.min(linesWithoutSpectra.size(), 5);
				System.out.printf("Showing first %d of %d spectra IDs that could not be found: \n\t%s\n", previewSize, linesWithoutSpectra.size(),
						String.join("\n\t", linesWithoutSpectra.subList(0, previewSize)));
			}
			
			//calculate zeroRT
			for(String pep : zTolRT.keySet()) {
				ArrayList<Double> rts = zTolRT.get(pep);
				double rtsum = 0;
				for(double v : rts)
					rtsum += v;
				avgzRT.put(pep, rtsum / rts.size());
			}
			
			//calculate metrics
			for(int i = 0; i < clines.size(); i++) {
				String [] crow = pf.data.get(clines.get(i)).split("\t");
				boolean isZero = (Math.abs(Double.parseDouble(crow[deltaCol])) <= peakTol);
				
				String key = crow[pepCol].trim();
				if(crow[modpepCol].trim().length() != 0)
					key = crow[modpepCol].trim();

				int rtSize = 0, specSimSize = 0;
				double rtDelta = -1e10;
				double avgSim = -1e10, avgZeroSim = -1e10;
				
				if(zTolRT.containsKey(key)) {
					rtDelta = Double.parseDouble(crow[rtCol]) - avgzRT.get(key);
					rtSize = zTolRT.get(key).size();
				}
				
				key += "." + crow[chargeCol];
				if(zTolSpecs.containsKey(key)) {
					Spectrum cspec = mr.getSpectrum(reNormName(crow[specCol]));
					if(cspec != null) {
						avgSim = cspec.averageSimilarity(zTolSpecs.get(key), ppmTol);
						avgZeroSim = avgzSim.get(key);
						specSimSize = zTolSpecs.get(key).size();
					}
				}
				out.printf("%s\t%s\t%s\t%s\t%d\t%.5f\t%d\t%.5f\t%.5f\t%d\n",crow[specCol],crow[pepCol],crow[modpepCol],crow[deltaCol],isZero?1:0,
						rtDelta, rtSize, avgSim, avgZeroSim, specSimSize);
			}
			
			out.flush();
			long t3 = System.currentTimeMillis();
			
			PTMShepherd.print(String.format("\t%s - %d (%d ms, %d ms)", cf, clines.size(), t2-t1,t3-t2));
		}
		
		out.close();
	}
	
	public void updateSimRTProfiles(SimRTProfile [] profiles) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(simRTFile));
		String cline;
		while((cline = in.readLine())!= null) {
			if(cline.equals("COMPLETE"))
				break;
			String [] sp = cline.split("\\t");
			double md = Double.parseDouble(sp[3]);
			for(int i = 0; i < profiles.length; i++) {
				int cind = profiles[i].locate.getIndex(md);
				if(cind != -1) 
					profiles[i].records[cind].updateWithLine(sp);
			}
		}
		in.close();
	}
}
