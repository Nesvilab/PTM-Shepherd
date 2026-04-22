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
	List<String> linesWithoutSpectra;
	
	public SiteLocalization(String dsName) {
		this.dsName = dsName;
		this.localizationFile = new File(PTMShepherd.normFName(dsName+".rawlocalize"));
	}

	
	public boolean isComplete() {
		try {
			if (localizationFile.exists()) {
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
		} catch (IOException e) {
			PTMShepherd.die("Error writing localization file: " + localizationFile.getAbsolutePath() + "\n" + e.getMessage());
		}
		return false;
	}
	
	public void complete() {
		try {
			PrintWriter out = new PrintWriter(new FileWriter(localizationFile,true));
			out.println("COMPLETE");
			out.close();
		} catch (IOException e) {
			PTMShepherd.die("Error writing localization file: " + localizationFile.getAbsolutePath() + "\n" + e.getMessage());
		}
	}
	
	
	public void localizePSMs(PSMFile pf, HashMap<String,File> mzMappings, boolean useMSFraggerLoc) {
		//assemble PSMs into per file groupings
		HashMap<String,ArrayList<Integer>> mappings = new HashMap<>();
		try {
			PrintWriter out = new PrintWriter(new FileWriter(localizationFile, true));

			//write headers
			out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Spectrum", "Peptide", "Mods", "Shift", "Localized_Pep",
					"MaxHyper_Unloc", "MaxHyper_Loc", "MaxPeaks_Unloc", "MaxPeaks_Loc");

			ppmTol = Double.parseDouble(PTMShepherd.getParam("spectra_ppmtol"));
			condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
			condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));
			linesWithoutSpectra = new ArrayList<>();
			int totalLines;

			if (useMSFraggerLoc && pf.msfraggerLocalizationCol == -1) {
				PTMShepherd.print(String.format("Warning! MSFragger localization requested, but localization columns not found in PSM file %s. MSFragger localization will not be used", pf.fname.toString()));
				useMSFraggerLoc = false;
			}

			initSpectrumMappings(pf, mappings);

			if (useMSFraggerLoc) {
				for (String cf : mappings.keySet()) { //cf = fraction
					long t1 = System.currentTimeMillis();
					ArrayList<Integer> clines = mappings.get(cf);
					totalLines = 0;
					for (Integer cline : clines) {
						out.println(annotateLineUsingMSFragger(pf, cline));
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
						out.println(annotateLine(pf, cline));
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
		} catch (IOException e) {
			PTMShepherd.die("IOError writing localized PSMs to localization file: " + localizationFile.getAbsolutePath() + "\n" + e.getMessage());
		}
	}

	public static void initSpectrumMappings(PSMFile pf, HashMap<String, ArrayList<Integer>> mappings) {
		for(int i = 0; i < pf.psms.size(); i++) {
			String bn = pf.psms.get(i).getFileName();
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

	public String annotateLine(PSMFile psmFile, int lineIndex) {
		StringBuilder sb = new StringBuilder();
		PSM psm = psmFile.psms.get(lineIndex);
		String seq = psm.getPeptide();
		float dmass = (float) psm.getDMass();
		float [] scores = new float[seq.length()];
		int [] frags = new int[seq.length()];
		String specName = psm.getSpec();

		sb.append(String.format("%s\t%s\t%s\t%.4f", specName, seq, psm.printAssignedMods(), dmass));
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

		localizeMods(psm.getAssignedMods(), mods);

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
	 * @return string: [Spectrum, Peptide, Assigned Mods, Delta Mass, Localized_Pep, MaxHyper_Unloc, MaxHyper_Loc, MaxPeaks_Unloc, MaxPeaks_Loc, scores, frags]
	 */
	public String annotateLineUsingMSFragger(PSMFile psmFile, int lineIndex) {
		StringBuilder sb = new StringBuilder();
		PSM psm = psmFile.psms.get(lineIndex);
		String seq = psm.getPeptide();
		float dmass = (float) psm.getDMass();
		float [] scores = new float[seq.length()];
		int [] frags = new int[seq.length()];
		String specName = psm.getSpec();

		sb.append(String.format("%s\t%s\t%s\t%.4f", specName, seq, psm.printAssignedMods(), dmass));

		double baseScore, maxScore;
		int baseFrags, maxFrags;
		String annoSeq;
		if (psm.spLine.get(psmFile.positionScoresCol).isEmpty()) {
			// no localization result from MSFragger
			baseScore = 0;
			maxScore = 0;
			baseFrags = 0;
			maxFrags = 0;
			annoSeq = seq;
		} else {
			baseScore = Float.parseFloat(psm.spLine.get(psmFile.scoreAllUnshiftedCol));
			baseFrags = Integer.parseInt(psm.spLine.get(psmFile.ionsAllUnshiftedCol));
			maxScore = Float.parseFloat(psm.spLine.get(psmFile.scoreBestPositionCol));
			maxFrags = Integer.parseInt(psm.spLine.get(psmFile.ionsBestPosCol));
			annoSeq = swapCase(psm.spLine.get(psmFile.msfraggerLocalizationCol));
			scores = extractMSFraggerScores(psm.spLine.get(psmFile.positionScoresCol), seq.length());
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
			try {
				scores[i] = Float.parseFloat(matcher.group(1));
			} catch (ArrayIndexOutOfBoundsException ex) {
				PTMShepherd.die("It appears that Crystal-C was run and using MSFragger localization in PTM-Shepherd was requested. This is not currently supported due to an issue with Crystal-C. Please disable the Use MSFragger Localization option (or re-run without Crystal-C) and try again.");
			}
			i++;
		}
		return scores;
	}

	public static void localizeMods(ArrayList<Mod> assignedMods, float[] mods) {
		for (Mod mod : assignedMods) {
			// Mods are 1-indexed in PSMs, except for terminal mods which are 0 (N-term) or length+1 (C-term)
			int pos;
            if (mod.position == 0) {
                pos = 0;                    // N-term
            } else if (mod.position == mods.length + 1) {
                pos = mods.length - 1;      // C-term
            } else {
                pos = mod.position - 1;     // all other mods
            }
			float mass = (float) mod.mass;
            mods[pos] += mass;
		}
	}

	public void updateLocalizationProfiles(LocalizationProfile [] profiles) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(localizationFile));
			String cline;
			in.readLine();
			while ((cline = in.readLine()) != null) {
				if (cline.equals("COMPLETE"))
					break;
				if (cline.endsWith("MISSINGSPECTRA"))
					continue;
				if (cline.startsWith("Spectrum"))
					continue;
				String[] sp = cline.split("\\t");
				double md = Double.parseDouble(sp[3]);
				for (int i = 0; i < profiles.length; i++) {
					int cind = profiles[i].locate.getIndex(md);
					if (cind != -1)
						profiles[i].records[cind].updateWithLine(sp);
					else /* ensure we count the specs that aren't matched to a bin for global calcs */
						profiles[i].records[profiles[i].records.length - 1].updateWithLine(sp);
				}
			}
			in.close();
		} catch (IOException e) {
			PTMShepherd.die("Error updating localization profile for localization file: " + localizationFile.getAbsolutePath() + "\n" + e.getMessage());
		}
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
