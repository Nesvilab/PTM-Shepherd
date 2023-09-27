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

import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.glyco.GlycanCandidate;
import edu.umich.andykong.ptmshepherd.glyco.GlycanResidue;
import edu.umich.andykong.ptmshepherd.glyco.GlycoAnalysis;
import edu.umich.andykong.ptmshepherd.glyco.StaticGlycoUtilities;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;
import org.apache.commons.lang3.tuple.ImmutablePair;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.zip.CRC32;

import static edu.umich.andykong.ptmshepherd.PTMShepherd.reNormName;

public class PSMFile {

	String [] headers;
	public ArrayList<String> data;
	public ArrayList<String> mappedRuns;
	public int dMassCol, precursorCol, assignedModCol, observedModCol, fraggerLocCol, peptideCol, modPeptideCol, deltaMassCol, calcMZcol, peptideCalcMassCol, chargeCol;
	public String prefType;

	private HashMap<String, Integer> scanToLine;
	public File fname;

	public PSMFile(String fn) throws Exception {
		this(new File(fn));
	}

	/**
	 * PSM class to hold parsed line info.
	 */
	public class PSM {
		private int lineNum; // 0 indexed starting from header, 1 indexed starting from data
		private ArrayList<String> spLine;
		private String spec;
		int specNum;
		private String pep;
		private ArrayList<ImmutablePair<Integer, Float>> mods;
		private float [] modArr;
		private Float dMass;

		PSM(int lineNum, String line) {
			this.lineNum = lineNum;
			this.spLine = new ArrayList<>(Arrays.asList(line.replace("\n","").split("\t", -1)));
			this.spec = null;
			this.pep = null;
			this.mods = null;
			this.modArr = null;
			this.dMass = null;
		}

		public void addValAtColumn(int colIdx, String val){
			this.spLine.add(colIdx, val);
		}

		public void addValAtColumn(String colName, String val){
			this.spLine.add(getColumn(colName), val);
		}

		public void replaceValAtColumn(int colIdx, String val){
			this.spLine.set(colIdx, val);
		}

		public void replaceValAtColumn(String colName, String val){
			this.spLine.set(getColumn(colName), val);
		}

		public String getSpec() {
			if (this.spec == null)
				this.spec = reNormName(spLine.get(getColumn("Spectrum")));
			return this.spec;
		}
		public String getPep() {
			if (this.pep == null)
				this.pep = spLine.get(getColumn("Peptide"));
			return this.pep;
		}

		public ArrayList<ImmutablePair<Integer,Float>> getMods() {
			if (this.mods == null) {
				this.mods = new ArrayList<>();
				String strMods = spLine.get(getColumn("Assigned Modifications"));
				if (!strMods.isEmpty()) {
					String[] spMods = strMods.split(",", -1);
					for (int i = 0; i < spMods.length; i++) {
						int p = spMods[i].indexOf("(");
						int q = spMods[i].indexOf(")");
						String spos = spMods[i].substring(0, p).trim();
						float mass = Float.parseFloat(spMods[i].substring(p + 1, q).trim());
						int pos;
						if (spos.equals("N-term"))
							pos = 0;
						else if (spos.equals("c"))
							pos = this.getPep().length();
						else
							pos = Integer.parseInt(spos.substring(0, spos.length() - 1));
						this.mods.add(new ImmutablePair<>(pos, mass));
					}
				}
			}
			return this.mods;
		}

		public float [] getModsAsArray() {
			if (this.modArr == null) {
				this.modArr = new float[this.getPep().length()];
				Arrays.fill(this.modArr, 0.0f);
				for (ImmutablePair<Integer, Float> mod : this.getMods()) {
					if (mod.getLeft() == 0)
						this.modArr[0] = mod.getRight();
					else
						this.modArr[mod.getLeft()-1] = mod.getRight();
				}
			}
			return this.modArr;
		}

		public float getDMass() {
			if (this.dMass == null)
				this.dMass = Float.parseFloat(this.spLine.get(getColumn("Delta Mass")));
			return this.dMass;
		}

		public String getColumnValue(String colName) {
			return this.spLine.get(getColumn(colName));
		}

		public String toString() {
			String newStr = String.join("\t", this.spLine);
			return newStr;
		}

		public void updateLine() {
			PSMFile.this.data.set(this.lineNum, this.toString());
		}

	}
	public int getColumn(String head) {
		for(int i = 0; i < headers.length; i++)
			if(headers[i].equals(head))
				return i;
		return -1;
	}

	public static String [] splitName(String fn) {
		String [] res = new String[2];
		if(fn.indexOf(".") < 0) {
			res[0] = fn;
			res[1] = "";
		} else {
			res[0] = fn.substring(0, fn.lastIndexOf("."));
			res[1] = fn.substring(fn.lastIndexOf(".")+1);
		}
		return res;
	}
	
	public static String getCRC32(File f) throws Exception {
		CRC32 crc = new CRC32();
		byte [] buf = new byte[1024*1024];
		int nread = 0;
		DataInputStream dis = new DataInputStream(new FileInputStream(f));
		while(dis.available() > 0) {
			nread = dis.read(buf);
			crc.update(buf,0,nread);
		}
		dis.close();
		return Long.toHexString(crc.getValue()) + Long.toHexString(f.length());
	}

	public PSM getLine(int i) {
		PSM psm = new PSM(i, data.get(i));
		return psm;
	}

	public HashSet<String> getRunNames() {
		HashSet<String> res = new HashSet<>();
		int col = getColumn("Spectrum");
		for(int i = 0; i < data.size(); i++) {
			String [] sp = data.get(i).split("\t");
			String d = sp[col];
			res.add(d.substring(0, d.indexOf(".")));
		}
		return res;
	}

	/**
	 * Makes a map of run names to line indices within the PSM file
	 * @return HashMap<String, ArrayList<Integer>> ran name -> list of line indices
	 * TODO make this automatically happen upon construction
	 */
	public HashMap<String,ArrayList<Integer>> getRunMappings() {
		HashMap<String,ArrayList<Integer>> mappings = new HashMap<>();
		for(int i = 0; i < this.data.size(); i++) {
			String [] sp = this.data.get(i).split("\t");
			String bn = sp[getColumn("Spectrum")].substring(0,sp[getColumn("Spectrum")].indexOf("."));
			if(!mappings.containsKey(bn))
				mappings.put(bn, new ArrayList<>());
			mappings.get(bn).add(i);
		}
		return mappings;
	}
	
	public ArrayList<Float> getMassDiffs() {
		ArrayList<Float> res = new ArrayList<>();
		int col = getColumn("Delta Mass");
		if (col == -1) {
			col = getColumn("Adjusted Delta Mass");
		}
		if (col == -1) {
			col = getColumn("Original Delta Mass");
		}
		this.dMassCol = col;
		for(int i = 0; i < data.size(); i++) {
			String [] sp = data.get(i).split("\t");
			res.add(Float.parseFloat(sp[col]));
		}
		return res;
	}

	public ArrayList<Double> getIntensities() {
		/* set up initial properties */
		int intPeaks = Integer.parseInt(PTMShepherd.getParam("histo_intensity"));
		int intCol = getColumn("Intensity");
		ArrayList<Double> ints = new ArrayList<>();

		/* if no intensity column or intensity no wanted, just get spectral counts */
		if (intCol == -1 || intPeaks == 0) {
			if (intPeaks == 1)
				System.out.printf("\tCould not identify 'Intensity' column in %s. Defaulting to spectral counts.", this.fname);
			for (int i = 0; i < this.data.size(); i++) {
				String[] sp = this.data.get(i).split("\t");
				ints.add(1.0);
			}
		}
		/* if intensity column found, collect counts */
		else {
			double total = 0;
			for (int i = 0; i < this.data.size(); i++) {
				String[] sp = this.data.get(i).split("\t");
				double cInt = Double.parseDouble(sp[intCol]);
				ints.add(cInt);
				total += cInt;
			}
			/* if column is invalid, redo calculation with spectral counts */
			if (total < 1) {
				ints = new ArrayList<>();
				System.out.printf("\tEmpty 'Intensity' column in %s. Defaulting to spectral counts.", this.fname);
				for (int i = 0; i < this.data.size(); i++) {
					String[] sp = this.data.get(i).split("\t");
					ints.add(1.0);
				}
			}
		}

		return ints;
	}

	public ArrayList<Float> getPrecursorMasses() {
		ArrayList<Float> precs = new ArrayList<>();
		int col = getPrecursorCol();
		for (int i = 0; i < data.size(); i++) {
			String [] sp = data.get(i).split("\t");
			precs.add(Float.parseFloat(sp[col]));
		}
		return precs;
	}

	public int getPrecursorCol() {
		int col = getColumn("Calibrated Observed Mass");
		if (col == -1)
			col = getColumn("Observed Mass");
		return col;
	}

	public static void getMappings(File path, HashMap<String,File> mappings, HashSet<String> runNames) {
		HashMap<String, Integer> datTypes = new HashMap<>();
		datTypes.put("mgf", 4);
		datTypes.put("mzBIN_cache", 7);
		datTypes.put("mzBIN", 3);
		datTypes.put("mzML", 2);
		datTypes.put("mzXML", 2);
		datTypes.put("raw", 1);

		int fileScore = 0;
		if(path.isDirectory()) {		
			File [] ls = path.listFiles();
			//get mapping for each file
			for(int i = 0; i < ls.length; i++) {
				getMappings(ls[i],mappings, runNames);
			}
		} else {
			String [] ns = splitName(path.getName());
			if (ns[0].contains("_calibrated")) {
				fileScore += 6;
				ns[0] = ns[0].substring(0, ns[0].indexOf("_calibrated"));
			}
			else if (ns[0].contains("_uncalibrated")) {
				fileScore += 5;
				ns[0] = ns[0].substring(0, ns[0].indexOf("_uncalibrated"));
			}
			if (!runNames.contains(ns[0]))
				return;
			if (mappings.containsKey(ns[0]) && (ns[1].equals("mzXML") || ns[1].equals("mzML") || ns[1].equals("raw") || ns[1].equals("mzBIN") || ns[1].equals("mgf")) || ns[1].equals("mzBIN_cache")) {
				if (mappings.get(ns[0]) == null)
					mappings.put(ns[0], path);
				else {
					File storedPath = mappings.get(ns[0]);
					String[] storedNs = splitName(storedPath.getName());
					fileScore += datTypes.get(ns[1]);
					if (fileScore > datTypes.get(storedNs[1]))
						mappings.put(ns[0], path);
				}
			}
		}
	}

	public TreeMap<String, Integer> getMS2Counts() {
		TreeMap<String, Integer> cnts = new TreeMap<>();
		int col = getColumn("Spectrum");
		for(int i = 0; i < data.size(); i++) {
			String [] sp = data.get(i).split("\t");
			String d = sp[col];
			String crun = d.substring(0, d.indexOf("."));
			if (!cnts.containsKey(crun))
				cnts.put(crun, 0);
			cnts.put(crun, cnts.get(crun) + 1);
		}
		return cnts;
	}

	/* Merges the rawglyco table onto the existing psm.tsv
	*  Update: only writes some columns rather than full rawglyco table
	*  numColsToUse gives the number of columns to take
	*/
	public void mergeGlycoTable(File glyf, int numColsToUse, boolean writeGlycansToAssignedMods, boolean nGlycan, String allowedResidues, boolean removeGlycanDeltaMass, boolean printDecoys) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(glyf), 1 << 22);
		String tempFoutName = this.fname + ".glyco.tmp";
		PrintWriter out = new PrintWriter(new FileWriter(tempFoutName));
		String[] glyHeaders = in.readLine().split("\t");

		/* Get glyco data */
		HashMap<String, String[]> glyLines = new HashMap<>();
		String cgline;
		int gSpecCol = 0; //todo: should be dynamically calculated
		while ((cgline = in.readLine()) != null) {
			String[] sp = cgline.split("\t", -1);
			glyLines.put(sp[gSpecCol], sp);
		}
		in.close();

		if (glyLines.size() < 2) {
			// no glycan information found (empty file or only line is "COMPLETE") - do not edit PSM table
			PTMShepherd.print("Warning: no modified spectra found, no glycans written to PSM table. Check input data and parameters");
			return;
		}

		/* Find headers, dynamically detect columns */
		observedModCol = getColumn("Observed Modifications");
		assignedModCol = getColumn("Assigned Modifications");
		fraggerLocCol = getColumn("MSFragger Localization");
		peptideCol = getColumn("Peptide");
		modPeptideCol = getColumn("Modified Peptide");
		deltaMassCol = getColumn("Delta Mass");
		chargeCol = getColumn("Charge");
		peptideCalcMassCol = getColumn("Calculated Peptide Mass");
		calcMZcol = getColumn("Calculated M/Z");
		if (warnPSMcolNotFound(observedModCol, "Observed Modifications")) {
			observedModCol = 27;		// default is 27
		}
		if (warnPSMcolNotFound(assignedModCol, "Assigned Modifications") || warnPSMcolNotFound(fraggerLocCol, "MSFragger Localization") || warnPSMcolNotFound(peptideCol, "Peptide") || warnPSMcolNotFound(modPeptideCol, "Modified Peptide") || warnPSMcolNotFound(deltaMassCol, "Delta Mass")) {
			writeGlycansToAssignedMods = false;		// can't write to assigned mods without these columns, disable
		}

		// get rawglyco file headers
		int mergeFromCol = -1;
		for (int i=0; i < glyHeaders.length; i++) {
			if (glyHeaders[i].matches("Best Glycan")) {
				mergeFromCol = i;
			}
		}
		if (mergeFromCol == -1) {
			System.out.println("Warning: Could not find Mass Shift column in rawglyco table, using default insert point");
			mergeFromCol = 5;
		}
		int bestTargetGlycanCol = mergeFromCol + 3;
		int bestTargetScoreCol = bestTargetGlycanCol + 1;
		int glycanScoreCol = mergeFromCol + 1;
		int glycanQvalCol = mergeFromCol + 2;

		ArrayList<String> newHeaders = new ArrayList<>();
		for (int i = 0; i < this.headers.length; i++)
			newHeaders.add(this.headers[i]);
		// add after observed mod column (observedModCol + 1)
		boolean hasPreviousGlycoInfo = hasGlycanAssignmentsWritten();
		if (!hasPreviousGlycoInfo) {
			// do not add glyco headers if the PSM table already has them
			newHeaders.addAll(observedModCol + 1, Arrays.asList(glyHeaders).subList(mergeFromCol + 1, mergeFromCol + numColsToUse));
			fraggerLocCol = fraggerLocCol + numColsToUse - 1;
		}
		out.println(String.join("\t", newHeaders));

		/* Match glycolines on PSM spectrum keys */
		int pSpecCol = getColumn("Spectrum");
		for (String cpline : this.data) {
			ArrayList<String> newLine = new ArrayList<>(Arrays.asList(cpline.split("\t")));
			String pSpec = newLine.get(pSpecCol);
			if (!hasPreviousGlycoInfo) {
				// add columns for glycan score and q-value before proceeding (for all lines, whether glycan-containing or not)
				for (int i=0; i < numColsToUse - 1; i++) {
					newLine.add(observedModCol + 1 + i, "");
				}
			}

			// check if a glycan was found
			if (glyLines.containsKey(pSpec)) {
				ArrayList<String> glyLine = new ArrayList<>(Arrays.asList(glyLines.get(pSpec)));
				if (glyLine.get(mergeFromCol).length() > 0) {
					String rawGlycan = glyLine.get(mergeFromCol);
					String observedGlycan;
					String glycanScore = glyLine.get(glycanScoreCol);
					if (rawGlycan.contains("Decoy")) {
						if (!printDecoys) {
							// report best target glycan instead of decoy (q-value will be reported as 1)
							observedGlycan = glyLine.get(bestTargetGlycanCol);
							glycanScore = glyLine.get(bestTargetScoreCol);
						} else {
							// user requested printing decoys, save decoy glycan (removing FailFDR if present)
							if (rawGlycan.contains("FailFDR")) {
								observedGlycan = rawGlycan.replace("FailFDR_", "");
							} else {
								observedGlycan = rawGlycan;
							}
						}
					} else if (rawGlycan.contains("FailFDR")) {
						observedGlycan = rawGlycan.replace("FailFDR_", "");
					} else {
						observedGlycan = rawGlycan;
					}

					// add the final glycan info to the line
					newLine.set(observedModCol, observedGlycan);
					newLine.set(observedModCol + 1, glycanScore);
					newLine.set(observedModCol + 2, glyLine.get(glycanQvalCol));
					int charge = Integer.parseInt(newLine.get(chargeCol));
					if (writeGlycansToAssignedMods) {
						newLine = writeGlycanToAssignedMod(newLine, rawGlycan, nGlycan, allowedResidues, removeGlycanDeltaMass, charge, writeGlycansToAssignedMods);
					}
				}
			}
			out.println(String.join("\t", newLine));
		}
		out.close();

		/* Move old file onto new file */
		Files.move(Paths.get(tempFoutName), Paths.get(String.valueOf(this.fname)), StandardCopyOption.REPLACE_EXISTING);
	}

	/**
	 * Determine if this PSM file has already had glycan assignment info written to it.
	 * Checks for the glycan header column name
	 * @return true if glycan info present
	 */
	public boolean hasGlycanAssignmentsWritten() {
		for (String header : headers) {
			if (header.contains("Glycan Score")) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Helper method for warning that given column in the PSM table was not found. If column is < 0, assumes the
	 * column was not found in the PSM table
	 * @param column index of the column in the PSM table
	 * @param columnName name of this column
	 * @return true if column NOT found (index < 0), false if found
	 */
	public boolean warnPSMcolNotFound(int column, String columnName) {
		if (column == -1) {
			// did not find correct column! Use default and warn user
			PTMShepherd.print(String.format("Warning: Could not find %s column in PSM table, PSM table editing may fail", columnName));
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Converts glycan ID to mass and position for Assigned Mods (may change to string for quant later). Uses
	 * MSFragger localization string from Philosopher (4.0.0+) if that contains localization info OR places
	 * glycan on the first allowed position if localization is ambiguous.
	 * Also writes to modified peptide and delta mass columns.
	 * Handles cases where information was previously written to the PSM table by removing/replacing the previous ID if present
	 */
	public ArrayList<String> writeGlycanToAssignedMod(ArrayList<String> newLine, String rawGlycan, boolean nGlycan, String allowedResidues, boolean removeGlycanDeltaMass, int charge, boolean writeAllGlycans) {
		// parse glycan composition from recently edited observed mods col
		TreeMap<GlycanResidue, Integer> glycanComp;
		String glycanOnly = rawGlycan.replace("FailFDR_", "").replace("Decoy_", "");
		boolean failOrDecoy = rawGlycan.contains("FailFDR") || rawGlycan.contains("Decoy");

		// Default is to always write glycan mass to assigned mod, unless option to only write glycans passing FDR is specified
		boolean editPSMGlycoEntry = true;
		if (!writeAllGlycans) {
			editPSMGlycoEntry = !failOrDecoy;
			if (removeGlycanDeltaMass) {
				// if removing delta mass for quant, ALWAYS edit PSM entry, even if it did not pass FDR (needed for quant)
				PTMShepherd.print("Note: remove_glycan_delta_mass requires that put_glycans_to_assigned_mods be set to True. All glycans will be written to assigned mods, regardless of FDR");
				editPSMGlycoEntry = true;
			}
		}

		/* Get glycan mass */
		try {
			glycanComp = StaticGlycoUtilities.parseGlycanString(glycanOnly);
		} catch (Exception ex) {
			try {
				// try old format in case of old PSM file
				glycanComp = StaticGlycoUtilities.parseGlycanStringOld(glycanOnly);
			} catch (Exception ex2) {
				// Not a glycan (PTM-S or Philosopher may put other string formats here) - ignore and continue
				return newLine;
			}
		}
		double glycanMass = GlycanCandidate.computeMonoisotopicMass(glycanComp);
		if (glycanMass == 0) {
			// do not edit entry if no matches found to the delta mass
			editPSMGlycoEntry = false;
			removeGlycanDeltaMass = false;
		}

		/* Get glycan location */
		int glycanLocation = readMSFraggerGlycanLocation(newLine, fraggerLocCol, peptideCol, nGlycan, allowedResidues);
		String glycanAA;
		String fraggerPepLocStr = newLine.get(fraggerLocCol);
		// skip missing loc column for now (quant will fail, but not needed for basic ID)
		try {
			glycanAA = fraggerPepLocStr.substring(glycanLocation, glycanLocation + 1).toUpperCase();
		} catch (StringIndexOutOfBoundsException ex) {
			PTMShepherd.print(String.format("ERROR: MSFragger localization not reported for spectrum %s. Spectrum will NOT have glycan put to assigned mods", newLine.get(0)));
			return newLine;
		}

		/* write mass and location to Assigned Mods */
		String currentAssignedMods = newLine.get(assignedModCol);
		String glycanMod = String.format("%d%s(%.4f)", glycanLocation + 1, glycanAA, glycanMass);	// site is 1-indexed in PSM table, not 0-indexed
		String[] modSplits = currentAssignedMods.split(", ");
		ArrayList<String> newModSplits = new ArrayList<>();

		// check if glycanAssignedMod already present and avoid double adding if so (in case of re-runs on same file)
		boolean prevGlycanWritten = false;
		double previousGlycanMass = 0;
		for (String mod: modSplits) {
			if (mod.length() > 0) {
				int modLocation = parseModLocation(mod);
				if (modLocation == glycanLocation + 1) {
					// a mod is already present at this site in the PSM table. Probably an existing glycan mod, but check the mass too to confirm
					double existingModMass = parseModMass(mod);
					if (glycanMass - 10 <= existingModMass && glycanMass + 10 >= existingModMass) {
						// found existing glycan annotation - do not add it to the updated mod list
						prevGlycanWritten = true;
						previousGlycanMass = existingModMass;
					} else {
						newModSplits.add(mod);
					}
				} else {
					// other mod
					newModSplits.add(mod);
				}
			}
		}
		// add the assigned glycan to the updated mod list (from which we removed any old glycan mods) if not failed FDR or is decoy
		if (editPSMGlycoEntry) {
			newModSplits.add(glycanMod);
		}
		String newAssignedMods = String.join(", ", newModSplits);
		newLine.set(assignedModCol, newAssignedMods);

		/* Write mass and location to Modified Pep */
		String modifiedPep = newLine.get(modPeptideCol);
		if (modifiedPep.length() == 0) {
			// no previous mods here, copy from peptide
			modifiedPep = newLine.get(peptideCol);
		}
		String newModPep = modifiedPep;
		int residueIndex = -1;
		for (int i=0; i < modifiedPep.length(); i++) {
			if (modifiedPep.charAt(i) >= 'A' && modifiedPep.charAt(i) <= 'Z') {
				// actual peptide residue, not modification info
				residueIndex++;
			}
			if (residueIndex == glycanLocation) {
				// found the glycan location. Check for existing glycan
				int AAindex = modifiedPep.charAt(i) - 65;		// capital alphabet starts at 65
				int roundedGlycanMass = (int) Math.round(glycanMass + AAMasses.monoisotopic_masses[AAindex]);
				int secondSubstringStart = i + 1;
				// Skip checking for modification if we've already reached the end of the modified peptide string
				if (!(secondSubstringStart == modifiedPep.length())) {
					if (modifiedPep.charAt(i + 1) > 'Z' || modifiedPep.charAt(i + 1) < 'A') {
						// non-residue character is next, so a modification has previously been written here. Remove it and replace
						int j = i + 1;
						// find length of the modification at this site (could vary if glycan has mass < 1000 or not)
						while (modifiedPep.charAt(j) > 'Z' || modifiedPep.charAt(j) < 'A') {
							j++;
							if (j == modifiedPep.length())
								break;
						}
						secondSubstringStart = j;
					}
				}
				if (editPSMGlycoEntry) {
					newModPep = modifiedPep.substring(0, i + 1) + String.format("[%d]", roundedGlycanMass) + modifiedPep.substring(secondSubstringStart);
				} else {
					// if failed or decoy, write no modification at this position (but do remove previous glycan if one was present from an older analysis)
					newModPep = modifiedPep.substring(0, i + 1) + modifiedPep.substring(secondSubstringStart);
				}
				break;
			}
		}
		newLine.set(modPeptideCol, newModPep);

		/* Update delta mass AND calc m/z columns */
		if (removeGlycanDeltaMass) {
			double prevDeltaMass = Double.parseDouble(newLine.get(deltaMassCol));
			double prevCalcMZ = Double.parseDouble(newLine.get(calcMZcol));
			double prevNeutralMass = Spectrum.mzToNeutralMass((float) prevCalcMZ, charge);
			double prevCalcPeptideMass = Double.parseDouble(newLine.get(peptideCalcMassCol));
			if (prevCalcPeptideMass - prevNeutralMass > 0.1 || prevCalcPeptideMass - prevNeutralMass < -0.1) {
				PTMShepherd.print(String.format("Warning: corrupted PSM table, calc peptide neutral mass and m/z do not match! PSM: %s", newLine));
			}

			if (prevGlycanWritten) {
				// A glycan was previously written to this line, and thus may have been subtracted from delta mass and added to peptide calc mass and MZ.
				// If DeltaMass was removed, add it back before subtracting new glycan's mass
				double correctedDelta, correctedMass;
				if (prevDeltaMass - previousGlycanMass < 10 && prevDeltaMass - previousGlycanMass > -10) {
					// a previous glycan was written, but the delta mass was NOT changed (i.e. full glycan mass still present in delta mass in psm table). No need to correct. +/- 10 to allow for isotope errors in delta mass
					correctedDelta = prevDeltaMass;
					correctedMass = prevCalcPeptideMass;
				} else {
					// delta mass WAS changed before. Add back previous glycan mass prior to subtracting the new glycan (in case the assigned glycan has changed)
					correctedDelta = prevDeltaMass + previousGlycanMass;
					correctedMass = prevCalcPeptideMass - previousGlycanMass;
				}
				newLine.set(deltaMassCol, String.format("%.4f", correctedDelta - glycanMass));
				newLine.set(peptideCalcMassCol, String.format("%.4f", correctedMass + glycanMass));
				newLine.set(calcMZcol, String.format("%.4f", Spectrum.neutralMassToMZ((float) (correctedMass + glycanMass), charge)));

			} else {
				// subtract glycan mass from delta mass, add to calc peptide mass and MZ
				newLine.set(deltaMassCol, String.format("%.4f", prevDeltaMass - glycanMass));
				newLine.set(peptideCalcMassCol, String.format("%.4f", prevCalcPeptideMass + glycanMass));
				newLine.set(calcMZcol, String.format("%.4f", Spectrum.neutralMassToMZ((float) (prevCalcPeptideMass + glycanMass), charge)));
			}
		}
		return newLine;
	}

	/**
	 * Read MSFragger localization string for a glycan and return the location (0-indexed)
	 * @param psmLineList ArrayList of splits from the line of the PSM table being read
	 * @param fraggerLocCol index of MSFragger localization column
	 * @param peptideCol index of the peptide column
	 * @param nGlycan if using N-glycan mode
	 * @param allowedResidues if not using N-glycan mode, the allowed residues (as a single string with no delimiters)
	 * @return 0-indexed location of the glycan (or of first allowed site if not localized by MSFragger)
	 */
	public int readMSFraggerGlycanLocation(ArrayList<String> psmLineList, int fraggerLocCol, int peptideCol, boolean nGlycan, String allowedResidues) {
		int glycanLocation = -1;
		String fraggerPepLocStr = psmLineList.get(fraggerLocCol);
		ArrayList<Integer> allowedPositions = new ArrayList<>();
		for (int i = 0; i < fraggerPepLocStr.length(); i++) {
			if (Character.isLowerCase(fraggerPepLocStr.charAt(i))) {
				allowedPositions.add(i);
			}
		}
		if (allowedPositions.size() >= 1) {
			// 1 or more positions - take first position if ambiguous
			glycanLocation = allowedPositions.get(0);
		} else {
			// no localization info provided by MSFragger. Take first allowed position
			if (nGlycan) {
				// find first sequon if Nglycan mode
				glycanLocation = GlycoAnalysis.findNGlycSequon(psmLineList.get(peptideCol));
			} else {
				// find first allowed residue if not NGlycan mode
				boolean[] allowedPos = SiteLocalization.parseAllowedPositions(psmLineList.get(peptideCol), allowedResidues);
				for (int i = 0; i < allowedPos.length; i++) {
					if (allowedPos[i]) {
						glycanLocation = i;
						break;
					}
				}
			}
		}
		return glycanLocation;
	}

	/* Add new column to PSM table in place to make it IonQuant compatible */
	public void preparePsmTableForIonQuant(double[][] peakBounds, int precUnits, double precTol) throws Exception {
		/* Check to make sure IonQuant column doesnt already exist */
		if (!(getColumn("Theoretical Modification Mass") == -1)) { /* Already exists */
			System.out.printf("\tPSM table at %s already IonQuant compatible\n",this.fname);
			return;
		}

		FastLocator locator = new FastLocator(peakBounds, precTol, precUnits);
		String tempFoutName = this.fname + ".iq.tmp";
		PrintWriter out = new PrintWriter(new FileWriter(tempFoutName));

		/* Write the new header */
		out.println(String.join("\t", this.headers) + "\tTheoretical Modification Mass");

		/* For each line in the file, find the peak apex from the delta mass*/
		for (int i = 0; i < this.data.size(); i++) {
			String[] sp = this.data.get(i).split("\t");
			double dmass = Double.parseDouble(sp[this.dMassCol]);
			double theoreticalDmass;
			if (locator.getIndex(dmass) == -1)
				theoreticalDmass = dmass;
			else
				theoreticalDmass = peakBounds[0][locator.getIndex(dmass)];
			out.println(this.data.get(i) + "\t" + String.format("%.4f", theoreticalDmass));
		}

		out.close();

		this.fname.delete();
		File newFileName = new File(tempFoutName);
		newFileName.renameTo(this.fname);
	}

	public PSMFile(File f) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(f), 1 << 22);
		this.fname = f;
		this.headers = in.readLine().split("\t");
		//find delta mass column for different philosopher versions
		int col = getColumn("Delta Mass");
		if (col == -1)
			col = getColumn("Adjusted Delta Mass");
		if (col == -1)
			col = getColumn("Original Delta Mass");
		this.dMassCol = col;
		this.data = new ArrayList<>();
		String cline;
		while((cline = in.readLine()) != null) {
			//if(cline.trim().length() > 0)
			if (cline.length() > 0)
				this.data.add(cline);
		}
		in.close();

		// Build index of PSM file lines to spectrum names
		this.scanToLine = new HashMap<>();
		for (int i = 0; i < this.data.size(); i++) {
			PSM tPSM = this.getLine(i);
			this.scanToLine.put(tPSM.getSpec(), i);
		}
	}

    public void annotateMassDiffs(String [] annotations) throws IOException {
		/* find column to modify, overwrite Observed Modifications col if exists */
		int annoCol = getColumn("Observed Modifications");
		boolean overwrite = true;
		if (annoCol == -1) {
			annoCol = getColumn("Assigned Modifications");
			overwrite = false;
		}
		if (annoCol == -1)
			annoCol = this.headers.length - 1;

		/* write annotations to lines */
		ArrayList<String> newLines = new ArrayList<>();
		for (int i = 0; i < this.data.size(); i++) {
			ArrayList<String> sp = new ArrayList<String>(Arrays.asList(this.data.get(i).split("\t")));
			if (overwrite == true)
				sp.set(annoCol, annotations[i]);
			else
				sp.add(annoCol, annotations[i]);
			newLines.add(String.join("\t", sp));
		}

		/* fix up headers */
		if (overwrite == false) {
			ArrayList<String> heads = new ArrayList<>(Arrays.asList(this.headers));
			heads.add(annoCol, "Observed Modifications");
			heads.toArray(this.headers);
		}

		/* write output */
		String tempFoutName = this.fname + ".anno.tmp";
		PrintWriter out = new PrintWriter(new FileWriter(tempFoutName));

		/* write the new header */
		out.println(String.join("\t", this.headers));
		/* write file lines */
		for (int i = 0; i < newLines.size(); i ++)
			out.println(newLines.get(i));

		/* close and rename temp file */
		out.close();
		this.fname.delete();
		File newFileName = new File(tempFoutName);
		newFileName.renameTo(this.fname);
    }

	/**
	 * Parse the location index from an Assigned Modification in the PSM table. If location is C-term,
	 * return -1 and if location is N-term, return -2
	 * @param mod string to parse
	 * @return mod location
	 */
    public static int parseModLocation(String mod) {
		// get all chars of the mod location and parse it
		String trimmedMod = mod.trim();
		int stopChar = 0;
		for (int i=0; i < trimmedMod.length(); i++) {
			if (!Character.isDigit(trimmedMod.charAt(i))) {
				// stop here
				stopChar = i;
				break;
			}
		}

		if (stopChar == 0) {
			// special case: N-term or C-term mod
			if (trimmedMod.startsWith("N-term")) {
				return -2;
			} else if (trimmedMod.startsWith("C-term")) {
				return -1;
			} else {
				PTMShepherd.print(String.format("Warning: invalid Assigned Modification format for mod %s. Not removing existing glycans", mod));
				return -3;
			}
		} else {
			return Integer.parseInt(trimmedMod.substring(0, stopChar));
		}
	}

	/**
	 * Parse the mass of an assigned modification in the PSM table from between the parentheses [e.g., 4N(1000.0000)]
	 * @param mod string to parse
	 * @return mod mass
	 */
	public static double parseModMass(String mod) {
    	String[] splits = mod.split("\\(");
    	String[] massSplits = splits[1].split("\\)");
    	return Double.parseDouble(massSplits[0]);
	}

	/**
	 * @param indx		insertion index
	 * @param newHead	new header to be inserted
	 * @return			-1 if header not already present and successfully inserted, header index if already present
	 */
	private int addHeader(int indx, String newHead) {
		// Check if column already exists
		int oldHeaderIndx = getColumn(newHead);

		// If header doesn't already exist, edit headers
		if (oldHeaderIndx == -1) {
			ArrayList<String> newHeaders = new ArrayList<String>(Arrays.asList(this.headers));
			newHeaders.add(indx, newHead);
			this.headers = newHeaders.toArray(new String[newHeaders.size()]);
		} else {
			System.out.printf("\t%s found in headers, overwriting existing column%n", newHead);
		}

		return oldHeaderIndx;
	}

	public void addColumn(int colIndx, String newHead, ArrayList<String> keys, ArrayList<String> vals) {
		// Check that PSM table editing will not fail
		if (vals.size() != this.data.size() || keys.size() != data.size()) {
			throw new ArrayIndexOutOfBoundsException("\tInput arrays and PSM table are not the same length, " +
					"editing PSM table will fail\n");
		}

		//TODO if column not found and inserting to the right, it will insert at the beginning of the table
		// This fixes it, but the error message should be different. Not sure how to check colIndx
		if (colIndx == 0)
			throw new ArrayIndexOutOfBoundsException("\tSpectrum is a protected column, refusing to overwrite\n");

		// Check that header doesn't already exist
		int existingHeaderIndx = addHeader(colIndx, newHead);
		if (existingHeaderIndx == -1) {
			for (int i = 0; i < keys.size(); i++) {
				String specKey = reNormName(keys.get(i)); // Automatically check for renormed name and apply
				int rowIndx = this.scanToLine.get(specKey);
				PSM tPSM = this.getLine(rowIndx);
				tPSM.addValAtColumn(newHead, vals.get(i));
				tPSM.updateLine();
			}
		} else {
			for (int i = 0; i < keys.size(); i++) {
				String specKey = reNormName(keys.get(i)); // Automatically check for renormed name and apply
				int rowIndx = this.scanToLine.get(specKey);
				PSM tPSM = this.getLine(rowIndx);
				tPSM.replaceValAtColumn(newHead, vals.get(i));
				tPSM.updateLine();
			}
		}
	}

	/**
	 * Returns the values of the column at colIndx.
	 * @param header column name
	 * @return	generic ArrayList of type Object
	 */
	public ArrayList<String> getColumnValues(String header) { //TODO I think you can make this a generic parser by making one param a callable of DataDtype
		int colIndx = getColumn(header);

		if (colIndx < 0 || colIndx >= (this.headers.length)) {
			throw new ArrayIndexOutOfBoundsException(String.format("Column index %d is out of bounds for a %d column" +
					"wide table", colIndx, this.headers.length));
		}

		ArrayList<String> values = new ArrayList<>(this.data.size());

		for (int i = 0; i < data.size(); i++) {
			PSM tPSM = this.getLine(i);
			values.add(tPSM.getColumnValue(header));
		}

		return values;
	}

	/**
	 * Returns the values of the column at colIndx.
	 * @param header column name
	 * @return	generic ArrayList of type Object
	 */
	public HashMap<String,String> getColumnValuesAndSpecs(String header) {
		int colIndx = getColumn(header);

		if (colIndx < 0 || colIndx >= (this.headers.length)) {
			throw new ArrayIndexOutOfBoundsException(String.format("Column index %d is out of bounds for a %d column" +
					"wide table", colIndx, this.headers.length));
		}

		HashMap<String,String> values = new HashMap<>(this.data.size());

		for (int i = 0; i < data.size(); i++) {
			PSM tPSM = this.getLine(i);
			values.put(tPSM.getSpec(), tPSM.getColumnValue(header));
		}

		return values;
	}

	public void save(boolean overwrite) throws IOException {
		String tempFoutName = this.fname + ".tmp";
		PrintWriter out = new PrintWriter(new FileWriter(tempFoutName));

		// Write lines to .tmp file file
		out.println(headersToString() + '\n');
		for (int i = 0; i < this.data.size(); i++) {
			out.println(this.getLine(i).toString());
		}

		out.close();

		if (overwrite)
			Files.move(Paths.get(tempFoutName), Paths.get(String.valueOf(this.fname)), StandardCopyOption.REPLACE_EXISTING);
	}

	private String headersToString() {
		String [] headers = Arrays.copyOf(this.headers, this.headers.length);
		//TODO Check with Dan, make sure we can remove trailing newlines character in this.headers
		headers[headers.length-1].replace("\n","");
		return String.join("\t", headers);
	}

	public void writeToPsmCache(File cacheFname, ArrayList<Double> varModMasses) throws FileNotFoundException {
		PrintWriter out = new PrintWriter(cacheFname);

		/* Write the new header */
		out.println(String.join("\t", this.headers));

		/* For each line in the file, modify accordingly */
		for (int i = 0; i < this.data.size(); i++) {
			ArrayList<String> sp = new ArrayList<>(Arrays.asList(this.data.get(i).split("\t")));
			//Shift variable mod info into delta mass
			if (varModMasses.size() > 0) {
				double oldDmass = Double.parseDouble(sp.get(this.dMassCol));
				double newDmass = oldDmass;
				double oldPepmass = Double.parseDouble(sp.get(getColumn("Calculated Peptide Mass")));
				double newPepmass = oldPepmass;
				String[] oldAssignedMods = sp.get(getColumn("Assigned Modifications")).split(",");
				ArrayList<String> newAssignedMods = new ArrayList<>(Arrays.asList(oldAssignedMods));
				ArrayList<Integer> dropMods = new ArrayList<>();
				for (int j = 0; j < oldAssignedMods.length; j++) {
					System.out.println(oldAssignedMods[j]);
					String modStr = oldAssignedMods[j];
					int p = modStr.indexOf("(");
					int q = modStr.indexOf(")");
					double modMass = Double.parseDouble(modStr.substring(p + 1, q).trim());
					for (int k = 0; k < varModMasses.size(); k++) {
						if (Math.abs(modMass - varModMasses.get(k)) < 0.001) {
							dropMods.add(j);
							newDmass += modMass;
							newPepmass -= modMass;
							break;
						}
					}
				}
				for (int k = dropMods.size() - 1; k > 0; k--)
					newAssignedMods.remove(k);
				// reassign adjusted values
				sp.set(this.deltaMassCol, Double.toString(newDmass));
				sp.set(getColumn("Calculated Peptide Mass"), Double.toString(newPepmass));
				sp.set(getColumn("Assigned Modifications"), String.join(", ", newAssignedMods));
				sp.add(Double.toString(oldDmass));
				sp.add(Double.toString(oldPepmass));
				sp.add(String.join(", ", newAssignedMods));
			}
			out.println(this.data.get(i));
		}

		out.close();
	}
}

