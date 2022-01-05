package edu.umich.andykong.ptmshepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.glyco.GlycanCandidate;
import edu.umich.andykong.ptmshepherd.glyco.GlycanResidue;
import edu.umich.andykong.ptmshepherd.glyco.GlycoAnalysis;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.zip.CRC32;

public class PSMFile {

	String [] headers;
	public ArrayList<String> data;
	public ArrayList<String> mappedRuns;
	public int dMassCol, precursorCol, assignedModCol, observedModCol, fraggerLocCol, peptideCol, modPeptideCol, deltaMassCol, calcMZcol, peptideCalcMassCol, chargeCol;
	public String prefType;
	public File fname;

	public PSMFile(String fn) throws Exception {
		this(new File(fn));
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
			datTypes.put("mzBIN_cache", 5);
			datTypes.put("mzBIN", 3);
			datTypes.put("mzML", 2);
			datTypes.put("mzXML", 2);
			datTypes.put("raw", 1);

		if(path.isDirectory()) {		
			File [] ls = path.listFiles();
			//get mapping for each file
			for(int i = 0; i < ls.length; i++) {
				getMappings(ls[i],mappings, runNames);
			}
		} else {
			String [] ns = splitName(path.getName());
			if (ns[0].contains("_calibrated"))
				ns[0] = ns[0].substring(0, ns[0].indexOf("_calibrated"));
			else if (ns[0].contains("_uncalibrated"))
				ns[0] = ns[0].substring(0, ns[0].indexOf("_uncalibrated"));
			if (!runNames.contains(ns[0]))
				return;
			if (mappings.containsKey(ns[0]) && (ns[1].equals("mzXML") || ns[1].equals("mzML") || ns[1].equals("raw") || ns[1].equals("mzBIN") || ns[1].equals("mgf")) || ns[1].equals("mzBIN_cache")) {
				if (mappings.get(ns[0]) == null)
					mappings.put(ns[0], path);
				else {
					File storedPath = mappings.get(ns[0]);
					String[] storedNs = splitName(storedPath.getName());
					if (datTypes.get(ns[1]) > datTypes.get(storedNs[1]))
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
			ArrayList<String> glyLine = new ArrayList<>(Arrays.asList(glyLines.get(pSpec)));

			if (!hasPreviousGlycoInfo) {
				// add columns for glycan score and q-value before proceeding
				newLine.addAll(observedModCol + 1, glyLine.subList(mergeFromCol + 1, mergeFromCol + numColsToUse));
			}

			// check if a glycan was found
			if (glyLine.get(mergeFromCol).length() > 0){
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
				newLine = writeGlycanToAssignedMod(newLine, rawGlycan, nGlycan, allowedResidues, removeGlycanDeltaMass, charge, writeGlycansToAssignedMods);

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
			glycanComp = PTMShepherd.parseGlycanString(glycanOnly);
		} catch (Exception ex) {
			// Not a glycan (PTM-S or Philosopher may put other string formats here) - ignore and continue
			return newLine;
		}
		double glycanMass = GlycanCandidate.computeMonoisotopicMass(glycanComp);

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
					if (glycanMass - 5 <= existingModMass && glycanMass + 5 >= existingModMass) {
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
			if (prevGlycanWritten) {
				// A glycan was previously written to this line, and thus was subtracted from delta mass. Add it back, then subtract new glycan's mass
				double correctedDelta = Double.parseDouble(newLine.get(deltaMassCol)) + previousGlycanMass;
				if (correctedDelta - 2*previousGlycanMass < 10 && correctedDelta - 2*previousGlycanMass > -10) {
					// a previous glycan was written, but the delta mass was not changed (i.e. full glycan mass still present in delta mass in psm table). No need to correct. +/- 10 to allow for isotope errors in delta mass
					correctedDelta = Double.parseDouble(newLine.get(deltaMassCol));
				}
				newLine.set(deltaMassCol, String.format("%.4f", correctedDelta - glycanMass));

				// Also correct calc m/z. A glycan was previously written, and thus its mass may have been added to calcMZ. Remove the old mass and add the new
				double prevCalcMZ = Double.parseDouble(newLine.get(calcMZcol));
				double prevNeutralMass = Spectrum.mzToNeutralMass((float) prevCalcMZ, charge);
				double calcPeptideMass = Double.parseDouble(newLine.get(peptideCalcMassCol));
				double pepMassMinusPrevGlyc = calcPeptideMass - previousGlycanMass;
				if (prevNeutralMass - pepMassMinusPrevGlyc > -5 && prevNeutralMass - calcPeptideMass - previousGlycanMass < 5) {
					// calc m/z had previous glycan mass added to it, so subtracting previous glycan mass from it equals 0 (+/- isotope errors)
					// remove the previous glycan mass from the calcMZ before adding the new one
					prevNeutralMass = prevNeutralMass - previousGlycanMass;
				}
				// previous mass either did not have previous glycan mass added OR it has already been removed above. Add glycan mass and correct the calc mz column
				double newMZ = Spectrum.neutralMassToMZ((float) (prevNeutralMass + glycanMass), charge);
				newLine.set(calcMZcol, String.format("%.4f", newMZ));
			} else {
				// subtract glycan mass from delta mass
				double observedDelta = Double.parseDouble(newLine.get(deltaMassCol));
				newLine.set(deltaMassCol, String.format("%.4f", observedDelta - glycanMass));

				// add glycan mass to the calculated m/z column. Convert to neutral masses to avoid double-adding protons
				double prevCalcMZ = Double.parseDouble(newLine.get(calcMZcol));
				double peptideMass = Spectrum.mzToNeutralMass((float) prevCalcMZ, charge);
				double newMZ = Spectrum.neutralMassToMZ((float) (peptideMass + glycanMass), charge);
				newLine.set(calcMZcol, String.format("%.4f", newMZ));
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
}
