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

package edu.umich.andykong.ptmshepherd.specsimilarity;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import edu.umich.andykong.ptmshepherd.PSM;
import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;

import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;
import edu.umich.andykong.ptmshepherd.utils.Variance;

import java.util.List;

public class SimRTAnalysis {

	String dsName;
	File simRTFile;
	MXMLReader mr;
	HashMap<String, MXMLReader> multiMr;
	double ppmTol, condRatio, peakTol;
	int condPeaks, precursorUnits;

	boolean calcIntensity;
	
	static final int MAX_ZERO_COMPARE = 20;
	
	public SimRTAnalysis(String dsName) {
		this.dsName = dsName;
		this.simRTFile = new File(PTMShepherd.normFName(dsName+".rawsimrt"));
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
	
	public void simrtPSMs(PSMFile pf, HashMap<String,File> mzMappings, boolean interRunComparisons) throws Exception {
		//assemble PSMs into per file groupings
		HashMap<String,ArrayList<Integer>> mappings = new HashMap<>();
		PrintWriter out = new PrintWriter(new FileWriter(simRTFile,true));
        calcIntensity = pf.intensityCol != -1;

		//Write header
		if (calcIntensity) {
			out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Spectrum", "Peptide", "Mod_Peptide", "Shift", "Is_Zero_Pep",
					"rt_shift", "nZeroSpecs_RT_shift", "Avg_Sim", "Avg_ZeroSim", "nZeroSpecs_Sim", "Avg_IntChange", "nZeroSpecs_Int");
		} else {
			out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Spectrum", "Peptide", "Mod_Peptide", "Shift", "Is_Zero_Pep",
					"rt_shift", "nZeroSpecs_RT_shift", "Avg_Sim", "Avg_ZeroSim", "nZeroSpecs_Sim");
		}
		
		ppmTol = Double.parseDouble(PTMShepherd.getParam("spectra_ppmtol"));
		condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
		condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));
		peakTol = Double.parseDouble(PTMShepherd.getParam("precursor_tol")); //determines zero bin
		double cPeakTol = peakTol;
		precursorUnits = Integer.parseInt(PTMShepherd.getParam("precursor_mass_units"));

		if (!interRunComparisons) {
			SiteLocalization.initSpectrumMappings(pf, mappings);
		} else {
			for (int i = 0; i < pf.psms.size(); i++) {
				String bn = pf.psms.get(i).getFileName(); //fraction
				if (!mappings.containsKey(bn))
					mappings.put(bn, new ArrayList<>());
			}
			for (int i = 0; i < pf.psms.size(); i++) {
				for (String fraction : mappings.keySet())
					mappings.get(fraction).add(i);
			}
		}

		//read all files into memory
		multiMr = new HashMap<>();
		long t1 = System.currentTimeMillis();
		for (String cf : mappings.keySet()) {
			mr = new MXMLReader(mzMappings.get(cf), Integer.parseInt(PTMShepherd.getParam("threads")));
			mr.readFully();
			multiMr.put(cf, mr);
		}
		long t2 = System.currentTimeMillis();
		PTMShepherd.print(String.format("\tSpectral data read into memory (%d ms reading)", t2-t1));

		//initialize these outside of run interactions for inter run comparisons
		HashMap<String,ArrayList<Integer>> zTolLines = new HashMap<>();
		HashMap<String,ArrayList<Spectrum>> zTolSpecs = new HashMap<>();
		HashMap<String,ArrayList<Double>> zTolRT = new HashMap<>();
		HashMap<String,ArrayList<Double>> zTolInt = new HashMap<>();
		HashMap<String,Double> avgzSim = new HashMap<>(); //{modPep.charge:avg zero sim in bin}
		HashMap<String,Double> avgzRT = new HashMap<>();
		HashMap<String,Double> avgzInt = new HashMap<>();

		//iterate and calculate similarity/retention time deltas for each file
		//if inter run comparisons is turned on, this look will be broken after first iteration
		for(String cf : mappings.keySet()) { //cf = fraction
			long t3 = System.currentTimeMillis();
			//leaving this here in case we need a less memory intensive application at some point
			//mr = new MXMLReader(mzMappings.get(cf), Integer.parseInt(PTMShepherd.getParam("threads")));
			//mr.readFully();
			ArrayList<Integer> clines = mappings.get(cf);
			//if not performing inter run comparisons, need to reset every time
			if (!interRunComparisons) {
				zTolLines = new HashMap<>();
				zTolSpecs = new HashMap<>();
				zTolRT = new HashMap<>();
				zTolInt = new HashMap<>();
				avgzSim = new HashMap<>(); //{modPep.charge:avg zero sim in bin}
				avgzRT = new HashMap<>();
				avgzInt = new HashMap<>();
			}

			//get zero bin data and calculate baselines
			for(int i = 0; i < clines.size(); i++) {
				PSM psm = pf.psms.get(clines.get(i));
				if(precursorUnits == 1)//ppm
					cPeakTol = calculatePeakTol(1500, peakTol, 0.0);
				boolean isZero = (psm.getDMass() <= cPeakTol);
				if(!isZero)
					continue;
				
				String key = generateKey(psm); //using pep seq as key

				if(!zTolRT.containsKey(key)) //structure {modpep:<rt>}
					zTolRT.put(key, new ArrayList<>());
				zTolRT.get(key).add(Double.parseDouble(psm.spLine.get(pf.retentionCol)));

				if (calcIntensity) {
					if (!zTolInt.containsKey(key)) //structure {modpep:<rt>}
						zTolInt.put(key, new ArrayList<>());
					zTolInt.get(key).add(Double.parseDouble(psm.spLine.get(pf.intensityCol)));
				}
				
				key += "." + psm.getCharge(); //structure {modpep.charge:<spec line>}
				if(!zTolLines.containsKey(key))
					zTolLines.put(key, new ArrayList<>());
				zTolLines.get(key).add(clines.get(i));
			}
			
			//calculate zeroSim
			int totalLines = 0;
			List<String> linesWithoutSpectra = new ArrayList<>();
			for(String pepZ : zTolLines.keySet()) {
				ArrayList<Integer> relLines = zTolLines.get(pepZ);
				Collections.shuffle(zTolLines.get(pepZ)); //shuffled to remove bias
				int nComp = Math.min(relLines.size(), MAX_ZERO_COMPARE);
				zTolSpecs.put(pepZ, new ArrayList<>());
				for(int i = 0; i < nComp; i++) {
					PSM psm = pf.psms.get(relLines.get(i));
					String targetFrac = psm.getFileName();
					zTolSpecs.get(pepZ).add(multiMr.get(targetFrac).getSpectrum(psm.getSpec()));
				}
				
				double zSimSum = 0;
				totalLines += relLines.size();
				for(int i = 0; i < relLines.size(); i++) {
					PSM psm = pf.psms.get(relLines.get(i));
					String specNormName = psm.getSpec();
					String targetFrac = psm.getFileName();
					Spectrum cspec = multiMr.get(targetFrac).getSpectrum(specNormName);
					if (cspec == null) {
						linesWithoutSpectra.add(specNormName);
						continue;
					}
					zSimSum += cspec.averageSimilarity(zTolSpecs.get(pepZ), ppmTol);
				}
				avgzSim.put(pepZ, zSimSum / relLines.size());
			}

            SiteLocalization.warnLinesWithoutSpectra(totalLines, linesWithoutSpectra);

            //calculate zeroRT
			for(String pep : zTolRT.keySet()) {
				ArrayList<Double> rts = zTolRT.get(pep);
				double rtsum = 0;
				for(double v : rts)
					rtsum += v;
				avgzRT.put(pep, rtsum / rts.size());
			}

			//calculate zeroQuant
			if (calcIntensity) {
				for(String pep : zTolInt.keySet()) {
					Variance ints = new Variance();
					ints.update(zTolInt.get(pep));
					//ints.logTransform();
					avgzInt.put(pep, ints.getMedian());
				}
			}
			
			//calculate metrics
			for(int i = 0; i < clines.size(); i++) {
				PSM psm = pf.psms.get(clines.get(i));
				if(precursorUnits == 1)//ppm
					cPeakTol = calculatePeakTol(1500, peakTol, 0.0);
				boolean isZero = (psm.getDMass() <= cPeakTol);

				String key = generateKey(psm); //using pep seq as key

				int rtSize = 0, specSimSize = 0, intSize = 0;
				double rtDelta = -1e20;
				double intDelta = -1e20;
				double avgSim = -1e20, avgZeroSim = -1e20;
				
				if(zTolRT.containsKey(key)) { //calculated against average RT time
					rtDelta = Double.parseDouble(psm.spLine.get(pf.retentionCol)) - avgzRT.get(key);
					rtSize = zTolRT.get(key).size();
				}

				if (calcIntensity) {
					if(zTolInt.containsKey(key)) { //calculated against average RT time
						intDelta = (Double.parseDouble(psm.spLine.get(pf.intensityCol))) / (avgzInt.get(key));
						intSize = zTolInt.get(key).size();
					}
				}
				
				key += "." + psm.getCharge(); //based on charge state
				if(zTolSpecs.containsKey(key)) {
					String targetFrac = psm.getFileName();
					Spectrum cspec = multiMr.get(targetFrac).getSpectrum(psm.getSpec());
					if(cspec != null) {
						avgSim = cspec.averageSimilarity(zTolSpecs.get(key), ppmTol); //all v all comparison
						avgZeroSim = avgzSim.get(key);
						specSimSize = zTolSpecs.get(key).size();
					}
				}
				if (calcIntensity) {
					out.printf("%s\t%s\t%s\t%s\t%d\t%.5f\t%d\t%.5f\t%.5f\t%d\t%.5f\t%d\n", psm.getSpec(), psm.getPeptide(), psm.getModifiedPeptide(), psm.getDMass(), isZero ? 1 : 0,
							rtDelta, rtSize, avgSim, avgZeroSim, specSimSize, intDelta, intSize);
				} else {
					out.printf("%s\t%s\t%s\t%s\t%d\t%.5f\t%d\t%.5f\t%.5f\t%d\n", psm.getSpec(), psm.getPeptide(), psm.getModifiedPeptide(), psm.getDMass(),isZero?1:0,
							rtDelta, rtSize, avgSim, avgZeroSim, specSimSize);
				}
			}
			
			out.flush();
			long t4 = System.currentTimeMillis();

			if (interRunComparisons) {
				PTMShepherd.print(String.format("\tProcessed - %d lines (%d ms processing)", clines.size(), t4-t3));
				break;
			}
			PTMShepherd.print(String.format("\t%s - %d lines (%d ms processing)", cf, clines.size(), t4-t3));
		}
		
		out.close();
	}
	
	public void updateSimRTProfiles(SimRTProfile [] profiles) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(simRTFile));
		String cline;
		in.readLine();
		while((cline = in.readLine())!= null) {
			if(cline.equals("COMPLETE"))
				break;
			if(cline.startsWith("Spectrum"))
				continue;
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

	public double calculatePeakTol(double pepmass, double ppmtol, double modmass){
		double peakTol = ((pepmass + modmass) / 1000000.0) * ppmtol;
		return peakTol;
	}

	public boolean getCalcIntensity() {
		return calcIntensity;
	}

	// key is peptide_assignedMods. Replaced modified peptide to support massdiff-to-varmod and reruns
	private String generateKey(PSM psm) {
		String key = psm.getPeptide(); //using pep seq as key
		if(!psm.getAssignedMods().isEmpty()) {
			key += "_" + psm.printAssignedMods();
		}
		return key;
	}
}
