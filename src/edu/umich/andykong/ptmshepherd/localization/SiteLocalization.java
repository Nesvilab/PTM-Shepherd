package edu.umich.andykong.ptmshepherd.localization;

import java.io.*;
import java.util.*;

import edu.umich.andykong.ptmshepherd.*;
import edu.umich.andykong.ptmshepherd.core.*;


public class SiteLocalization {

	String dsName;
	File localizationFile;
	MXMLReader mr;
	double ppmTol, condRatio;
	int condPeaks;
	int specCol, pepCol, modCol, deltaCol;
	
	public SiteLocalization(String dsName) {
		this.dsName = dsName;
		localizationFile = new File(dsName+".rawlocalize");
	}
	
	public String reNormName(String s) {
		String [] sp = s.split("\\.");
		int sn = Integer.parseInt(sp[1]);
		//with charge state
		//return String.format("%s.%d.%d.%s",sp[0],sn,sn,sp[3]);
		//without charge state
		return String.format("%s.%d.%d",sp[0],sn,sn);
	}
	
	public boolean isComplete() throws Exception {
		if(localizationFile.exists()) {
			try (RandomAccessFile raf = new RandomAccessFile(localizationFile, "r")) {
				raf.seek(Math.max(0, localizationFile.length() - 20));
				String cline;
				while ((cline = raf.readLine()) != null)
					if (cline.equals("COMPLETE")) {
						raf.close();
						return true;
					}
			}
			localizationFile.delete();
		}
		return false;
	}
	
	public void complete() throws Exception {
		PrintWriter out = new PrintWriter(new FileWriter(localizationFile,true));
		out.println("COMPLETE");
		out.close();
	}
	
	
	public void localizePSMs(PSMFile pf, HashMap<String,File> mzMappings) throws Exception {
		//assemble PSMs into per file groupings
		HashMap<String,ArrayList<Integer>> mappings = new HashMap<>();
		PrintWriter out = new PrintWriter(new FileWriter(localizationFile,true));
		
		specCol = pf.getColumn("Spectrum");
		pepCol = pf.getColumn("Peptide");
		modCol = pf.getColumn("Assigned Modifications");
		deltaCol = pf.dMassCol;
		ppmTol = Double.parseDouble(PTMShepherd.getParam("spectra_ppmtol"));
		condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
		condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));

		for(int i = 0; i < pf.data.size(); i++) {
			String [] sp = pf.data.get(i).split("\t");
			String bn = sp[specCol].substring(0,sp[specCol].indexOf("."));
			if(!mappings.containsKey(bn))
				mappings.put(bn, new ArrayList<>());
			mappings.get(bn).add(i);
		}
		
		//iterate and localize each file
		for(String cf : mappings.keySet()) {
			long t1 = System.currentTimeMillis();
			mr = new MXMLReader(mzMappings.get(cf));
			mr.readFully();
			long t2 = System.currentTimeMillis();
			ArrayList<Integer> clines = mappings.get(cf);
			for(int i = 0; i < clines.size(); i++) {
				try {
					out.println(annotateLine(pf.data.get(clines.get(i))));
				} catch(Exception e) {
					e.printStackTrace();
					System.out.println("Error in: " +pf.data.get(clines.get(i)));
				}
			}
			out.flush();
			long t3 = System.currentTimeMillis();
			
			PTMShepherd.print(String.format("\t%s - %d (%d ms, %d ms)", cf, clines.size(), t2-t1,t3-t2));
		}
		out.close();
	}
	
	public String annotateLine(String line) throws Exception {
		StringBuffer sb = new StringBuffer();
		String [] sp = line.split("\\t");
		String seq = sp[pepCol];
		float dmass = Float.parseFloat(sp[deltaCol]);
		float [] scores = new float[seq.length()];
		int [] frags = new int[seq.length()];
		String specName = sp[specCol];
		String [] smods = sp[modCol].split(",");
		
		sb.append(String.format("%s\t%s\t%s\t%.4f", specName,seq,sp[modCol],dmass));
		
		Spectrum spec = mr.getSpectrum(reNormName(specName));
		
		if(spec == null) {
			sb.append("\tMISSINGSPECTRA");
			PTMShepherd.print("Cannot get spec: " + specName);
			return sb.toString();
		}
		
		spec.condition(condPeaks, condRatio);
		
		float [] mods = new float[seq.length()];
		
		for(int i = 0; i < smods.length; i++) {
			smods[i] = smods[i].trim();
			if(smods[i].length() == 0)
				continue;
			int p = smods[i].indexOf("(");
			int q = smods[i].indexOf(")");
			String spos = smods[i].substring(0, p).trim();
			double mass = Double.parseDouble(smods[i].substring(p+1, q).trim());
			int pos = -1;
			if(spos.equals("N-term")) {
				pos = 0;
				
//				This subtraction is necessary when the over mass is reported instead of the mass difference
//				mass -= AAMasses.monoisotopic_nterm_mass;
			}
			else if(spos.equals("c")) {
				pos = mods.length - 1;
//				This subtraction is necessary when the over mass is reported instead of the mass difference
//				mass -= (AAMasses.monoisotopic_cterm_mass + AAMasses.protMass);
			}
			else
				pos = Integer.parseInt(spos.substring(0,spos.length()-1)) - 1;
			mods[pos] += mass;
		}
		
		
		float baseScore = spec.getHyper(seq, mods, ppmTol);
		int baseFrags = spec.getFrags(seq,mods, ppmTol);
		float maxScore = baseScore;
		int maxFrags = baseFrags;
		for(int i = 0; i < seq.length(); i++) {
			mods[i] += dmass;
			scores[i] = spec.getHyper(seq, mods, ppmTol);
			frags[i] = spec.getFrags(seq,mods, ppmTol);
			if(frags[i] > maxFrags) 
				maxFrags = frags[i];
			if(scores[i] > maxScore)
				maxScore = scores[i];
			mods[i] -= dmass;
		}
		
		StringBuffer annoSeq = new StringBuffer();
		for(int i  = 0; i < seq.length(); i++) {
			if(frags[i] == maxFrags)
				annoSeq.append(seq.charAt(i));
			else
				annoSeq.append((char)(seq.charAt(i)+('a'-'A')));
		}
		
		sb.append(String.format("\t%s\t%.2f\t%.2f\t%d\t%d",annoSeq.toString(),baseScore,maxScore,baseFrags,maxFrags));
		
		for(int i = 0; i < scores.length; i++)
			sb.append(String.format("\t%.2f\t%d", scores[i],frags[i]));
		
		return sb.toString();
	}
	
	public void updateLocalizationProfiles(LocalizationProfile [] profiles) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(localizationFile));
		String cline;
		while((cline = in.readLine())!= null) {
			if(cline.equals("COMPLETE"))
				break;
			if(cline.endsWith("MISSINGSPECTRA"))
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
}
