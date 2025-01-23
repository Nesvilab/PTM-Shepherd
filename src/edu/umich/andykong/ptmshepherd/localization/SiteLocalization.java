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

package edu.umich.andykong.ptmshepherd.localization;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.umich.andykong.ptmshepherd.*;
import edu.umich.andykong.ptmshepherd.core.*;
import static edu.umich.andykong.ptmshepherd.PTMShepherd.reNormName;


public class SiteLocalization {

	String dsName;
	File localizationFile;
	MXMLReader mr;
	double ppmTol, condRatio;
	int condPeaks;
	int specCol, pepCol, assignedModCol, deltaCol;
	List<String> linesWithoutSpectra;
	
	public SiteLocalization(String dsName) {
		this.dsName = dsName;
		this.localizationFile = new File(PTMShepherd.normFName(dsName+".rawlocalize"));
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
	
	
	public void localizePSMs(PSMFile pf, HashMap<String,File> mzMappings, boolean useMSFraggerLoc) throws Exception {
		//assemble PSMs into per file groupings
		HashMap<String,ArrayList<Integer>> mappings = new HashMap<>();
		PrintWriter out = new PrintWriter(new FileWriter(localizationFile,true));

		//write headers
		out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n","Spectrum","Peptide","Mods","Shift","Localized_Pep",
				"MaxHyper_Unloc", "MaxHyper_Loc", "MaxPeaks_Unloc", "MaxPeaks_Loc");

		specCol = pf.specCol;
		pepCol = pf.peptideCol;
		assignedModCol = pf.assignedModCol;
		deltaCol = pf.dMassCol;
		ppmTol = Double.parseDouble(PTMShepherd.getParam("spectra_ppmtol"));
		condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
		condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));
		linesWithoutSpectra = new ArrayList<>();
		int totalLines;

		if (useMSFraggerLoc && pf.msfraggerLocalizationCol == -1) {
			throw new Exception(String.format("MSFragger localization requested, but localization columns not found in PSM file %s.", pf.fname.toString()));
		}

		initSpectrumMappings(pf, mappings, specCol);

		if (useMSFraggerLoc) {
			for (String cf : mappings.keySet()) { //cf = fraction
				long t1 = System.currentTimeMillis();
				ArrayList<Integer> clines = mappings.get(cf);
				totalLines = 0;
				for (Integer cline : clines) {
					out.println(annotateLineUsingMSFragger(pf, pf.data.get(cline)));
					totalLines++;
				}
				totalLines--;
				out.flush();

				long t2 = System.currentTimeMillis();
				warnLinesWithoutSpectra(totalLines, linesWithoutSpectra);
				PTMShepherd.print(String.format("\t%s - %d lines (%d ms processing)", cf, clines.size(), t2 - t1));
			}
		} else {
			//iterate and localize each file
			for (String cf : mappings.keySet()) { //cf = fraction
				long t1 = System.currentTimeMillis();
				mr = new MXMLReader(mzMappings.get(cf), Integer.parseInt(PTMShepherd.getParam("threads")));
				mr.readFully();
				long t2 = System.currentTimeMillis();
				ArrayList<Integer> clines = mappings.get(cf);
				totalLines = 0;
				for (Integer cline : clines) {
					out.println(annotateLine(pf.data.get(cline)));
					totalLines++;
				}
				totalLines--;
				out.flush();
				long t3 = System.currentTimeMillis();

				warnLinesWithoutSpectra(totalLines, linesWithoutSpectra);

				PTMShepherd.print(String.format("\t%s - %d lines (%d ms reading, %d ms processing)", cf, clines.size(), t2 - t1, t3 - t2));
			}
		}
		out.close();
	}

	public static void initSpectrumMappings(PSMFile pf, HashMap<String, ArrayList<Integer>> mappings, int specCol) {
		for(int i = 0; i < pf.data.size(); i++) {
			String [] sp = pf.data.get(i).split("\t");
			String bn = sp[specCol].substring(0,sp[specCol].indexOf("."));
			if(!mappings.containsKey(bn))
				mappings.put(bn, new ArrayList<>());
			mappings.get(bn).add(i);
		}
	}

	public static void warnLinesWithoutSpectra(int totalLines, List<String> linesWithoutSpectra) {
		if (!linesWithoutSpectra.isEmpty()) {
			System.out.printf("\tCould not find %d/%d (%.1f%%) spectra.\n", linesWithoutSpectra.size(), totalLines,
					100.0*((double) linesWithoutSpectra.size()/totalLines));
			int previewSize = Math.min(linesWithoutSpectra.size(), 5);
			System.out.printf("\tShowing first %d of %d spectra IDs that could not be found: \n\t%s\n", previewSize, linesWithoutSpectra.size(),
					String.join("\n\t\t", linesWithoutSpectra.subList(0, previewSize)));
		}
	}

	public String annotateLine(String line) {
		StringBuilder sb = new StringBuilder();
		String [] sp = line.split("\\t");
		String seq = sp[pepCol];
		float dmass = Float.parseFloat(sp[deltaCol]);
		float [] scores = new float[seq.length()];
		int [] frags = new int[seq.length()];
		String specName = sp[specCol];
		String [] assignedMods = sp[assignedModCol].split(",");

		sb.append(String.format("%s\t%s\t%s\t%.4f", specName, seq, sp[assignedModCol], dmass));
		Spectrum spec = mr.getSpectrum(reNormName(specName));

		boolean [] allowedPoses = parseAllowedPositions(seq, PTMShepherd.getParam("localization_allowed_res"));

		if(spec == null) {
			sb.append("\tMISSINGSPECTRA");
			linesWithoutSpectra.add(specName);
			//PTMShepherd.print("Cannot get spec: " + specName);
			return sb.toString();
		}
		
		spec.condition(condPeaks, condRatio);
		
		float [] mods = new float[seq.length()];

		localizeMods(assignedMods, mods);

		float baseScore = spec.getHyper(seq, mods, ppmTol);
		int baseFrags = spec.getFrags(seq, mods, ppmTol);
		float maxScore = baseScore;
		int maxFrags = baseFrags;
		for(int i = 0; i < seq.length(); i++) {
			if (allowedPoses[i])
				mods[i] += dmass;
			scores[i] = spec.getHyper(seq, mods, ppmTol);
			frags[i] = spec.getFrags(seq,mods, ppmTol);
			if(frags[i] > maxFrags) 
				maxFrags = frags[i];
			if(scores[i] > maxScore)
				maxScore = scores[i];
			if (allowedPoses[i])
				mods[i] -= dmass;
		}
		
		StringBuilder annoSeq = new StringBuilder();
		for(int i  = 0; i < seq.length(); i++) {
			if(frags[i] == maxFrags)
				annoSeq.append(seq.charAt(i));
			else
				annoSeq.append((char)(seq.charAt(i)+('a'-'A')));
		}
		
		sb.append(String.format("\t%s\t%.2f\t%.2f\t%d\t%d",annoSeq,baseScore,maxScore,baseFrags,maxFrags));
		
		for(int i = 0; i < scores.length; i++)
			sb.append(String.format("\t%.2f\t%d", scores[i],frags[i]));
		
		return sb.toString();
	}

	/**
	 * Use the MSFragger localization result instead of recalculating the localization, but format
	 * so that it can be used for the downstream localization summary.
	 * @param line PSM line
	 * @return string: [Spectrum, Peptide, Assigned Mods, Delta Mass, Localized_Pep, MaxHyper_Unloc, MaxHyper_Loc, MaxPeaks_Unloc, MaxPeaks_Loc, scores, frags]
	 */
	public String annotateLineUsingMSFragger(PSMFile psmFile, String line) {
		StringBuilder sb = new StringBuilder();
		String [] sp = line.split("\\t");
		String seq = sp[pepCol];
		float dmass = Float.parseFloat(sp[deltaCol]);
		float [] scores = new float[seq.length()];
		int [] frags = new int[seq.length()];
		String specName = sp[specCol];
		String [] assignedMods = sp[assignedModCol].split(",");

		sb.append(String.format("%s\t%s\t%s\t%.4f", specName, seq, sp[assignedModCol], dmass));

		double baseScore, maxScore;
		int baseFrags, maxFrags;
		String annoSeq;
		if (sp[psmFile.positionScoresCol].isEmpty()) {
			// no localization result from MSFragger
			baseScore = 0;
			maxScore = 0;
			baseFrags = 0;
			maxFrags = 0;
			annoSeq = seq;
		} else {
			baseScore = Float.parseFloat(sp[psmFile.scoreAllUnshiftedCol]);
			baseFrags = Integer.parseInt(sp[psmFile.ionsAllUnshiftedCol]);
			maxScore = Float.parseFloat(sp[psmFile.scoreBestPositionCol]);
			maxFrags = Integer.parseInt(sp[psmFile.ionsBestPosCol]);
			annoSeq = swapCase(sp[psmFile.msfraggerLocalizationCol]);
			scores = extractMSFraggerScores(sp[psmFile.positionScoresCol], seq.length());
		}

		sb.append(String.format("\t%s\t%.2f\t%.2f\t%d\t%d",annoSeq,baseScore,maxScore,baseFrags,maxFrags));

		for(int i = 0; i < scores.length; i++)
			sb.append(String.format("\t%.2f\t%d", scores[i],frags[i]));		// note: frags (number of ions at each site) is not recorded by MSFragger, and is left blank

		return sb.toString();
	}

	public String swapCase(String input) {
		StringBuilder swapped = new StringBuilder(input.length());
		for (int i = 0; i < input.length(); i++) {
			char c = input.charAt(i);
			if (Character.isUpperCase(c)) {
				swapped.append(Character.toLowerCase(c));
			} else {
				swapped.append(Character.toUpperCase(c));
			}
		}
		return swapped.toString();
	}

	/**
	 * Extract the position scores from the MSFragger scores string. Format is "A(0.1)B(0.2)C(0.3)..."
	 */
	private float[] extractMSFraggerScores(String positionScores, int length) {
		Pattern pattern = Pattern.compile("\\(([\\d.]+)\\)");
		Matcher matcher = pattern.matcher(positionScores);
		float[] scores = new float[length];

		int i = 0;
		while (matcher.find()) {
			scores[i] = Float.parseFloat(matcher.group(1));
			i++;
		}
		return scores;
	}

	public static void localizeMods(String[] smods, float[] mods) {
		for(int i = 0; i < smods.length; i++) {
			smods[i] = smods[i].trim();
			if(smods[i].isEmpty())
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
	}

	public void updateLocalizationProfiles(LocalizationProfile [] profiles) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(localizationFile));
		String cline;
		in.readLine();
		while((cline = in.readLine())!= null) {
			if(cline.equals("COMPLETE"))
				break;
			if(cline.endsWith("MISSINGSPECTRA"))
				continue;
			if(cline.startsWith("Spectrum"))
				continue;
			String [] sp = cline.split("\\t");
			double md = Double.parseDouble(sp[3]);
			for(int i = 0; i < profiles.length; i++) {
				int cind = profiles[i].locate.getIndex(md);
				if(cind != -1)
					profiles[i].records[cind].updateWithLine(sp);
				else /* ensure we count the specs that aren't matched to a bin for global calcs */
					profiles[i].records[profiles[i].records.length - 1].updateWithLine(sp);
			}
		}
		in.close();
	}

	public static boolean [] parseAllowedPositions(String seq, String allowedReses){
		boolean [] allowedPoses = new boolean[seq.length()];
		if (allowedReses.equals("all") || allowedReses.equals(""))
			Arrays.fill(allowedPoses, true);
		else {
			Arrays.fill(allowedPoses, false);
			for (int i = 0; i < seq.length(); i++) {
				for (int j = 0; j < allowedReses.length(); j++) {
					if (seq.charAt(i) == allowedReses.charAt(j)) {
						allowedPoses[i] = true;
						break;
					}
				}
			}
		}
		return allowedPoses;
	}
}
