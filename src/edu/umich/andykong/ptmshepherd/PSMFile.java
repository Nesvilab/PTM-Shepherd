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
import edu.umich.andykong.ptmshepherd.glyco.GlycoAnalysis;
import edu.umich.andykong.ptmshepherd.glyco.GlycoParams;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;
import org.apache.commons.lang3.tuple.ImmutablePair;
import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanParser;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.CRC32;

import static edu.umich.andykong.ptmshepherd.PTMShepherd.reNormName;

public class PSMFile {

	String [] headers;
	public ArrayList<PSM> psms;
	public int dMassCol, precursorCol, assignedModCol, observedModCol, fraggerLocCol, peptideCol, modPeptideCol,
			calcMZcol, peptideCalcMassCol, chargeCol, intensityCol, specCol, msfraggerLocalizationCol, positionScoresCol,
			bestPositionsCol, ionsBestPosCol, scoreBestPositionCol, scoreAllUnshiftedCol, ionsAllUnshiftedCol,
			eValCol, retentionCol, glycanCompCol, glycanScoreCol, glycanQvalCol;

	private final HashMap<String, Integer> scanToLineMap;
	public File fname;
	boolean alreadyWarned;
	public static final Pattern massPattern = Pattern.compile("\\(([-.\\d]+)\\)");

	/**
	 * PSM class to hold parsed line info.
	 */
	public class PSM {
		private int lineNum; // 0 indexed starting from header, 1 indexed starting from data
		public ArrayList<String> spLine;
		private final String spec;
		private String fileName;
		private int specNum;
		private String pep;
		private ArrayList<ImmutablePair<Integer, Float>> mods;
		private float [] modArr;
		private Float dMass;

		PSM(int lineNum, String line) {
			this.lineNum = lineNum;
			this.spLine = new ArrayList<>(Arrays.asList(line.replace("\n","").split("\t")));
			this.fileName = null;
			this.spec = reNormName(spLine.get(getColumn("Spectrum")));
			this.specNum = -1;
			this.pep = null;
			this.mods = null;
			this.modArr = null;
			this.dMass = null;
		}

		public String printLine() {
			return String.join("\t", spLine);
		}

		public void addValAtColumn(int colIdx, String val){
			spLine.add(colIdx, val);
		}

		public void replaceValAtColumn(int colIdx, String val){
			spLine.set(colIdx, val);
		}

		public String getSpec() {
			return spec;
		}

		public int getSpecNum() {
			if (specNum == -1) {
				String[] spSpec = spec.split("\\.", -1);
				specNum = Integer.parseInt(spSpec[spSpec.length-2]);
			}
			return specNum;
		}

		public String getFileName() {
			if (fileName == null)
				fileName = spec.substring(0, spec.indexOf("."));
			return fileName;
		}

		public String getPep() {
			if (pep == null)
				pep = spLine.get(getColumn("Peptide"));
			return pep;
		}

		public int getCharge() {
			return Integer.parseInt(spLine.get(getColumn("Charge")));
		}

		public String[] getAssignedMods() {
			String mods = spLine.get(getColumn("Assigned Modifications"));
			if (mods.isEmpty())
				return new String[0];
			return mods.split(",");
		public Float getCalcPepmass() {
			if (calcPepMass == null) {
				calcPepMass = Float.parseFloat(spLine.get(getColumn("Calculated Peptide Mass")));
			}
			return calcPepMass;
		}

		public String getAssignedModsStr() {
			return spLine.get(getColumn("Assigned Modifications"));
		}

		public ArrayList<ImmutablePair<Integer,Float>> getMods() {
			if (mods == null) {
				mods = new ArrayList<>();
				String strMods = spLine.get(getColumn("Assigned Modifications"));
				if (!strMods.isEmpty()) {
					String[] spMods = strMods.split(",", -1);
                    for (String spMod : spMods) {
                        int p = spMod.indexOf("(");
                        int q = spMod.indexOf(")");
                        String spos = spMod.substring(0, p).trim();
                        float mass = Float.parseFloat(spMod.substring(p + 1, q).trim());
                        int pos;
                        if (spos.equals("N-term"))
                            pos = 0;
                        else if (spos.equals("c"))
                            pos = this.getPep().length();
                        else
                            pos = Integer.parseInt(spos.substring(0, spos.length() - 1));
                        mods.add(new ImmutablePair<>(pos, mass));
                    }
				}
			}
			return mods;
		}

		public float [] getModsAsArray() {
			if (modArr == null) {
				modArr = new float[getPep().length()];
				Arrays.fill(modArr, 0.0f);
				for (ImmutablePair<Integer, Float> mod : getMods()) {
					if (mod.getLeft() == 0)
						modArr[0] = mod.getRight();
					else
						modArr[mod.getLeft()-1] = mod.getRight();
				}
			}
			return modArr;
		}

		// todo: handle adjusted/original and PTV
		public float getDMass() {
			if (dMass == null)
				dMass = Float.parseFloat(spLine.get(getColumn("Delta Mass")));
			return dMass;
		}

		public String getColumnValue(String colName) {
			return spLine.get(getColumn(colName));
		}

		public ArrayList<String> getSpLine() {
			return spLine;
		}

		public String toString() {
            return String.join("\t", this.spLine);
		}

	}

	public int getColumn(String head) {
		for(int i = 0; i < headers.length; i++)
			if(headers[i].equals(head))
				return i;
		return -1;
	}

	private void initColumns() {
		/* dynamically detect columns */
		observedModCol = getColumn("Observed Modifications");
		assignedModCol = getColumn("Assigned Modifications");
		fraggerLocCol = getColumn("MSFragger Localization");
		peptideCol = getColumn("Peptide");
		modPeptideCol = getColumn("Modified Peptide");
		chargeCol = getColumn("Charge");
		peptideCalcMassCol = getColumn("Calculated Peptide Mass");
		calcMZcol = getColumn("Calculated M/Z");
		intensityCol = getColumn("Intensity");
		specCol = getColumn("Spectrum");
		retentionCol = getColumn("Retention");
		eValCol = getColumn("Expectation");
		msfraggerLocalizationCol = getColumn("MSFragger Localization");
		positionScoresCol = getColumn("Position Scores");
		bestPositionsCol = getColumn("Best Positions");
		scoreBestPositionCol = getColumn("Score Best Position");
		ionsBestPosCol = getColumn("Ions Best Position");
		scoreAllUnshiftedCol = getColumn("Score All Unshifted");
		ionsAllUnshiftedCol = getColumn("Ions All Unshifted");
		glycanCompCol = getColumn("Total Glycan Composition");
		glycanScoreCol = getColumn("Glycan Score");
		glycanQvalCol = getColumn("Glycan q-value");
		precursorCol = getPrecursorCol();

		//find delta mass column for different philosopher versions
		int col = getColumn("Delta Mass");
		if (col == -1)
			col = getColumn("Adjusted Delta Mass");
		if (col == -1)
			col = getColumn("Original Delta Mass");
		dMassCol = col;
	}

	public static String [] splitName(String fn) {
		String [] res = new String[2];
		if(!fn.contains(".")) {
			res[0] = fn;
			res[1] = "";
		} else {
			res[0] = fn.substring(0, fn.lastIndexOf("."));
			res[1] = fn.substring(fn.lastIndexOf(".")+1);
		}
		return res;
	}

	/**
	 * Uses name with file extension to support handling user-named files like "_calibrated.raw", which became
	 * "_calibrated_calibrated.mzML" and crashed before. Returns filename without extension after removing
	 * _(un)calibrated.mzML/MGF.
	 */
	private static String removeCalTag(String fileBaseName) {
		String nameWithExt;
		if (fileBaseName.contains("_calibrated.mzML"))
			nameWithExt = fileBaseName.replace("_calibrated.mzML", ".mzML");
		else if (fileBaseName.contains("_uncalibrated.mzML"))
			nameWithExt = fileBaseName.replace("_uncalibrated.mzML", ".mzML");
		else if (fileBaseName.contains("_calibrated.MGF"))
			nameWithExt = fileBaseName.replace("_calibrated.MGF", ".MGF");
		else if (fileBaseName.contains("_uncalibrated.MGF"))
			nameWithExt = fileBaseName.replace("_uncalibrated.MGF", ".MGF");
		else
			nameWithExt = fileBaseName;
		return splitName(nameWithExt)[0];
	}
	
	public static String getCRC32(File f) throws Exception {
		CRC32 crc = new CRC32();
		byte [] buf = new byte[1024*1024];
		int nread = 0;
		DataInputStream dis = new DataInputStream(Files.newInputStream(f.toPath()));
		while(dis.available() > 0) {
			nread = dis.read(buf);
			crc.update(buf,0,nread);
		}
		dis.close();
		return Long.toHexString(crc.getValue()) + Long.toHexString(f.length());
	}

	/**
	 * Returns a PSM based on the file line
	 * @param i file line, doesn't include headers
	 * @return PSM
	 */
	public PSM getLine(int i) {
        return psms.get(i);
	}

	public HashSet<String> getRunNames() {
		HashSet<String> res = new HashSet<>();
        for (PSM psm : psms) {
			res.add(psm.getFileName());
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
		for(int i = 0; i < psms.size(); i++) {
			String fileName = psms.get(i).getFileName();
			if(!mappings.containsKey(fileName))
				mappings.put(fileName, new ArrayList<>());
			mappings.get(fileName).add(i);
		}
		return mappings;
	}
	
	public ArrayList<Float> getMassDiffs() {
		ArrayList<Float> res = new ArrayList<>();
		for (PSM psm : psms) {
            res.add(Float.parseFloat(psm.spLine.get(dMassCol)));
        }
		return res;
	}

	public ArrayList<ArrayList<Float>> getMassDiffsWithVarmods(boolean useAssignedMods) {
		ArrayList<ArrayList<Float>> res = new ArrayList<>();
		for (PSM psm : psms) {
			ArrayList<Float> psmMods = new ArrayList<>();
			psmMods.add(Float.parseFloat(psm.spLine.get(dMassCol)));
			if (useAssignedMods) {
				String mods = psm.spLine.get(assignedModCol);
				if (!mods.isEmpty()) {
					String[] modArr = mods.split(",");
					for (String mod : modArr) {
						Matcher m = massPattern.matcher(mod);
						if (m.find()) {
							psmMods.add(Float.parseFloat(m.group(1)));
						}
					}
				}
			}
			res.add(psmMods);
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
			for (PSM psm : psms) {
                ints.add(1.0);
            }
		}
		/* if intensity column found, collect counts */
		else {
			double total = 0;
			for (PSM psm : psms) {
                double cInt = Double.parseDouble(psm.spLine.get(intCol));
                ints.add(cInt);
                total += cInt;
            }
			/* if column is invalid, redo calculation with spectral counts */
			if (total < 1) {
				ints = new ArrayList<>();
				System.out.printf("\tEmpty 'Intensity' column in %s. Defaulting to spectral counts.", fname);
				for (PSM psm : psms) {
                    ints.add(1.0);
                }
			}
		}

		return ints;
	}

	public ArrayList<Float> getPrecursorMasses() {
		ArrayList<Float> precs = new ArrayList<>();
		for (PSM psm : psms) {
            precs.add(Float.parseFloat(psm.spLine.get(precursorCol)));
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
		// File priority list, this will return the first one matched so insertion order must be consistent
		LinkedHashMap<String, Integer> priorities = new LinkedHashMap<>();
		priorities.put(".mzBIN_cache", 20);
		priorities.put("_calibrated.mzML", 19);
		priorities.put("_uncalibrated.mzML", 18);
		priorities.put("_calibrated.mgf", 17);
		priorities.put("_uncalibrated.mgf", 16);
		priorities.put(".mzML", 15);
		priorities.put(".mzXML", 14);
		priorities.put(".mzBIN", 13);
		priorities.put(".mgf", 4);
		priorities.put(".raw", 1);
		priorities.put("None", 0);

		// Recursively search all directories
		if(path.isDirectory()) {
			File [] ls = path.listFiles();
			// get mapping for each file
            for (File l : ls) {
                getMappings(l, mappings, runNames);
            }
		} else { // see if valid file ext
			String matchedKey = getMatchingExtension(path, priorities);
			if (matchedKey != null) { // end of name exists in priorities map
				String rawFileName = removeCalTag(path.getName());
				// If raw file not part of this analysis, continue
				if (!runNames.contains(rawFileName))
					return;
				// New raw file info
				int newRawFilePriority = priorities.get(matchedKey);
				// Existing raw file info
				File previousRawFile = mappings.getOrDefault(rawFileName, null);
				if (previousRawFile == null) // no file mapped yet
					mappings.put(rawFileName, path);
				else { // Compare and replace if appropriate
					int previousRawFilePriority = priorities.get(getMatchingExtension(previousRawFile, priorities));
					if (newRawFilePriority > previousRawFilePriority) {
						mappings.put(rawFileName, path);
					}
				}
			}
		}
	}

	private static String getMatchingExtension(File file, HashMap<String, Integer> priority) {
		for (String key : priority.keySet()) {
			if (file.getName().endsWith(key)) {
				return key;  // Returns the extension key if the file name ends with it
			}
		}
		return null;  // No matching extension found
	}


	public TreeMap<String, Integer> getMS2Counts() {
		TreeMap<String, Integer> cnts = new TreeMap<>();
		for (PSM psm : psms) {
            String crun = psm.getFileName();
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
	public void mergeGlycoTable(File glyf, int numColsToUse, GlycoParams glycoParams) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(glyf), 1 << 22);
		String tempFoutName = this.fname + ".glyco.tmp";
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


		if (warnPSMcolNotFound(observedModCol, "Observed Modifications")) {
			observedModCol = 27;		// default is 27
		}
		if (warnPSMcolNotFound(assignedModCol, "Assigned Modifications") || warnPSMcolNotFound(fraggerLocCol, "MSFragger Localization") || warnPSMcolNotFound(peptideCol, "Peptide") || warnPSMcolNotFound(modPeptideCol, "Modified Peptide") || warnPSMcolNotFound(dMassCol, "Delta Mass")) {
			glycoParams.writeGlycansToAssignedMods = false;		// can't write to assigned mods without these columns, disable
		}

		// get rawglyco file headers
		int mergeFromCol = -1;
		for (int i=0; i < glyHeaders.length; i++) {
			if (glyHeaders[i].matches(GlycoAnalysis.GLYCAN_COMP_COL_NAME)) {
				mergeFromCol = i;
			}
		}
		if (mergeFromCol == -1) {
			System.out.println("Warning: Could not find Mass Shift column in rawglyco table, using default insert point");
			mergeFromCol = 5;
		}
		int bestTargetGlycanCol = mergeFromCol + 3;
		int bestTargetScoreCol = bestTargetGlycanCol + 1;
		int rawGlycoScoreCol = mergeFromCol + 1;
		int rawGlycoQvalCol = mergeFromCol + 2;

		boolean hasPreviousGlycoInfo = hasGlycanAssignmentsWritten();

		/* Match glycolines on PSM spectrum keys */
		ArrayList<String> psmKeys = new ArrayList<>();
		ArrayList<String> glycanComps = new ArrayList<>();
		ArrayList<String> glycanScores = new ArrayList<>();
		ArrayList<String> glycanQvals = new ArrayList<>();
		for (PSM psm : psms) {
			psmKeys.add(psm.getSpec());

			// check if a glycan was found
			if (glyLines.containsKey(psm.getSpec())) {
				ArrayList<String> glyLine = new ArrayList<>(Arrays.asList(glyLines.get(psm.getSpec())));
				if (!glyLine.get(mergeFromCol).isEmpty()) {
					String rawGlycan = glyLine.get(mergeFromCol);
					String observedGlycan;
					String glycanScore = glyLine.get(rawGlycoScoreCol);
					if (rawGlycan.contains("Decoy")) {
						if (!glycoParams.printGlycoDecoys) {
							// report best target glycan instead of decoy (q-value will be reported as 1)
							rawGlycan = glyLine.get(bestTargetGlycanCol);
							glycanScore = glyLine.get(bestTargetScoreCol);
						}
					} 
					if (rawGlycan.contains("FailFDR")) {
						observedGlycan = rawGlycan.replace("FailFDR_", "");
					} else {
						observedGlycan = rawGlycan;
					}

					// save glycan info directly or to the lists to add columns to the PSM table later
					if (!hasPreviousGlycoInfo) {
						glycanComps.add(observedGlycan);
						glycanScores.add(glycanScore);
						glycanQvals.add(glyLine.get(rawGlycoQvalCol));
					} else {
						psm.spLine.set(glycanCompCol, observedGlycan);
						psm.spLine.set(glycanScoreCol, glycanScore);
						psm.spLine.set(glycanQvalCol, glyLine.get(rawGlycoQvalCol));
					}
					// update assigned mods column
					if (glycoParams.writeGlycansToAssignedMods) {
						writeGlycanToAssignedMod(psm, rawGlycan, glycoParams);
					}
				}
			} else {
				// no glycan found for this spectrum
				if (!hasPreviousGlycoInfo) {
					glycanComps.add("");
					glycanScores.add("");
					glycanQvals.add("");
				} else {
					psm.spLine.set(glycanCompCol, "");
					psm.spLine.set(glycanScoreCol, "");
					psm.spLine.set(glycanQvalCol, "");
				}
			}
		}

		if (!hasPreviousGlycoInfo) {
			// add new glycan columns to the PSM table if not previously written
			addColumn(observedModCol + 1, "Glycan q-value", psmKeys, glycanQvals);
			addColumn(observedModCol + 1, "Glycan Score", psmKeys, glycanScores);
			addColumn(observedModCol + 1, "Total Glycan Composition", psmKeys, glycanComps);
		}

		save(true);
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
	public void writeGlycanToAssignedMod(PSM psm, String rawGlycan, GlycoParams glycoParams) {
		// parse glycan composition from recently edited observed mods col
		String glycanOnly = rawGlycan.replace("FailFDR_", "").replace("Decoy_", "");
		boolean failOrDecoy = rawGlycan.contains("FailFDR") || rawGlycan.contains("Decoy");

		// Default is to always write glycan mass to assigned mod, unless option to only write glycans passing FDR is specified
		boolean editPSMGlycoEntry = true;
		if (!glycoParams.writeGlycansToAssignedMods) {
			editPSMGlycoEntry = !failOrDecoy;
			if (glycoParams.removeGlycanDeltaMass) {
				// if removing delta mass for quant, ALWAYS edit PSM entry, even if it did not pass FDR (needed for quant)
				PTMShepherd.print("Note: remove_glycan_delta_mass requires that put_glycans_to_assigned_mods be set to True. All glycans will be written to assigned mods, regardless of FDR");
				editPSMGlycoEntry = true;
			}
		}

		if (glycanOnly.contains("no target matches") || glycanOnly.contains("No Glycan Matched")) {
			return;		// skip, no target glycan info to propagate
		}

		/* Get glycan mass */
		Glycan glyc;
		glyc = GlycanParser.parseGlycanString(glycanOnly, glycoParams.glycanResiduesMap);
		if (glyc == null) {
			// try old format in case of old PSM file
			glyc = GlycanParser.parseOldPTMSGlycanString(glycanOnly, glycoParams.glycanResiduesMap);
			if (glyc == null) {
				// Not a glycan (PTM-S may put other string formats here) - ignore and continue
				return;
			}
		}
		double glycanMass = glyc.mass;
		if (glycanMass == 0) {
			// do not edit entry if no matches found to the delta mass
			editPSMGlycoEntry = false;
			glycoParams.removeGlycanDeltaMass = false;
		}

		/* Get glycan location */
		int glycanLocation = readMSFraggerGlycanLocation(psm.spLine, glycoParams.nGlycan, glycoParams.allowedLocalizationResidues);
		String glycanAA;
		// skip missing loc column for now (quant will fail, but not needed for basic ID)
		try {
			glycanAA = psm.getPep().substring(glycanLocation, glycanLocation + 1).toUpperCase();
		} catch (StringIndexOutOfBoundsException ex) {
			if (!alreadyWarned) {
				PTMShepherd.print(String.format("WARNING: invalid MSFragger localization reported for spectrum %s. Spectrum will NOT have glycan put to assigned mods. Please check the Allowed Residues parameter if not in N-glyco mode.", psm.getSpec()));
				alreadyWarned = true;
			}
			return;
		}

		/* write mass and location to Assigned Mods */
		String currentAssignedMods = psm.getAssignedModsStr();
		String glycanMod = String.format("%d%s(%.4f)", glycanLocation + 1, glycanAA, glycanMass);	// site is 1-indexed in PSM table, not 0-indexed
		String[] modSplits = currentAssignedMods.split(", ");
		ArrayList<String> newModSplits = new ArrayList<>();

		// check if glycanAssignedMod already present and avoid double adding if so (in case of re-runs on same file)
		boolean prevGlycanWritten = false;
		double previousGlycanMass = 0;
		for (String mod: modSplits) {
			if (!mod.isEmpty()) {
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
		psm.spLine.set(assignedModCol, newAssignedMods);
		// todo: update PSM's Assigned mods!

		/* Write mass and location to Modified Pep */
		String modifiedPep = psm.spLine.get(modPeptideCol);
		if (modifiedPep.isEmpty()) {
			// no previous mods here, copy from peptide
			modifiedPep = psm.spLine.get(peptideCol);
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
		psm.spLine.set(modPeptideCol, newModPep);

		/* Update delta mass AND calc m/z columns */
		if (glycoParams.removeGlycanDeltaMass) {
			// todo: update with standardized delta mass
			double prevDeltaMass = psm.getDMass();
			double prevCalcMZ = Double.parseDouble(psm.spLine.get(calcMZcol));
			double prevNeutralMass = Spectrum.mzToNeutralMass((float) prevCalcMZ, psm.getCharge());
			double prevCalcPeptideMass = Double.parseDouble(psm.spLine.get(peptideCalcMassCol));
			if (prevCalcPeptideMass - prevNeutralMass > 0.1 || prevCalcPeptideMass - prevNeutralMass < -0.1) {
				PTMShepherd.print(String.format("Warning: corrupted PSM table, calc peptide neutral mass and m/z do not match! PSM: %s", psm.printLine()));
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
				// todo: handle PTV
				psm.spLine.set(dMassCol, String.format("%.4f", correctedDelta - glycanMass));
				psm.spLine.set(peptideCalcMassCol, String.format("%.4f", correctedMass + glycanMass));
				psm.spLine.set(calcMZcol, String.format("%.4f", Spectrum.neutralMassToMZ((float) (correctedMass + glycanMass), psm.getCharge())));

			} else {
				// todo: handle PTV
				// subtract glycan mass from delta mass, add to calc peptide mass and MZ
				psm.spLine.set(dMassCol, String.format("%.4f", prevDeltaMass - glycanMass));
				psm.spLine.set(peptideCalcMassCol, String.format("%.4f", prevCalcPeptideMass + glycanMass));
				psm.spLine.set(calcMZcol, String.format("%.4f", Spectrum.neutralMassToMZ((float) (prevCalcPeptideMass + glycanMass), psm.getCharge())));
			}
		}
	}

	/**
	 * Read MSFragger localization string for a glycan and return the location (0-indexed)
	 * @param psmLineList ArrayList of splits from the line of the PSM table being read
	 * @param nGlycan if using N-glycan mode
	 * @param allowedResidues if not using N-glycan mode, the allowed residues (as a single string with no delimiters)
	 * @return 0-indexed location of the glycan (or of first allowed site if not localized by MSFragger)
	 */
	public int readMSFraggerGlycanLocation(ArrayList<String> psmLineList, boolean nGlycan, String allowedResidues) {
		int glycanLocation = -1;
		String fraggerPepLocStr = psmLineList.get(fraggerLocCol);
		ArrayList<Integer> allowedPositions = new ArrayList<>();
		for (int i = 0; i < fraggerPepLocStr.length(); i++) {
			if (Character.isLowerCase(fraggerPepLocStr.charAt(i))) {
				allowedPositions.add(i);
			}
		}
		if (!allowedPositions.isEmpty()) {
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
		String tempFoutName = fname + ".iq.tmp";
		PrintWriter out = new PrintWriter(new FileWriter(tempFoutName));

		/* Write the new header */
		out.println(String.join("\t", headers) + "\tTheoretical Modification Mass");

		/* For each line in the file, find the peak apex from the delta mass*/
        for (PSM psm: psms) {
            double dmass = psm.getDMass();
            double theoreticalDmass;
            if (locator.getIndex(dmass) == -1)
                theoreticalDmass = dmass;
            else
                theoreticalDmass = peakBounds[0][locator.getIndex(dmass)];
            out.println(String.join("\t", psm.spLine) + "\t" + String.format("%.4f", theoreticalDmass));
        }

		out.close();

		fname.delete();
		File newFileName = new File(tempFoutName);
		newFileName.renameTo(fname);
	}

	/**
	 * Reads all lines into PSMs in the psms array and stores the scan to line index mapping in scanToLineMap.
	 */
	public PSMFile(File f) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(f), 1 << 22);
		fname = f;
		headers = in.readLine().split("\t");
		initColumns();

		psms = new ArrayList<>();
		scanToLineMap = new HashMap<>();
		int i = 0;
		String cline;
		while((cline = in.readLine()) != null) {
			if (!cline.isEmpty()) {
				PSM thisPSM = new PSM(i, cline);
				psms.add(thisPSM);
				scanToLineMap.put(thisPSM.spec, i);
				i++;
			}
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
		for (PSM psm: psms) {
			if (overwrite)
				psm.replaceValAtColumn(annoCol, annotations[psm.lineNum]);
			else
				psm.addValAtColumn(annoCol, annotations[psm.lineNum]);
			newLines.add(String.join("\t", psm.spLine));
		}

		/* fix up headers */
		if (!overwrite) {
			ArrayList<String> heads = new ArrayList<>(Arrays.asList(this.headers));
			heads.add(annoCol, "Observed Modifications");
			heads.toArray(this.headers);
		}

		/* write output */
		// todo: skip this and only write out at the end??
		String tempFoutName = this.fname + ".anno.tmp";
		PrintWriter out = new PrintWriter(new FileWriter(tempFoutName));

		/* write the new header */
		out.println(String.join("\t", this.headers));
		/* write file lines */
        for (String newLine : newLines) out.println(newLine);

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
			ArrayList<String> newHeaders = new ArrayList<>(Arrays.asList(headers));
			newHeaders.add(indx, newHead);
			headers = newHeaders.toArray(new String[newHeaders.size()]);
		} else {
			System.out.printf("\t%s found in headers, overwriting existing column%n", newHead);
		}

		return oldHeaderIndx;
	}

	public void addColumn(int colIndx, String newHead, ArrayList<String> keys, ArrayList<String> vals) {
		// Check that PSM table editing will not fail
		if (vals.size() != psms.size() || keys.size() != psms.size()) {
			throw new ArrayIndexOutOfBoundsException("Input arrays and PSM table are not the same length, " +
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
				int rowIndx = scanToLineMap.get(specKey);
				psms.get(rowIndx).addValAtColumn(colIndx, vals.get(i));
			}
		} else {
			for (int i = 0; i < keys.size(); i++) {
				String specKey = reNormName(keys.get(i)); // Automatically check for renormed name and apply
				int rowIndx = scanToLineMap.get(specKey);
				psms.get(rowIndx).replaceValAtColumn(colIndx, vals.get(i));
			}
		}

		// re-initialize column indices in case any have changed
		initColumns();
	}

	/**
	 * Returns the values of the column at colIndx.
	 * @param header column name
	 * @return	generic ArrayList of type Object
	 */
	public ArrayList<String> getColumnValues(String header) { //TODO I think you can make this a generic parser by making one param a callable of DataDtype
		int colIndx = getColumn(header);

		if (colIndx < 0 || colIndx >= (headers.length)) {
			throw new ArrayIndexOutOfBoundsException(String.format("Cannot fetch %s column. Column index %d is out of bounds for a %d column" +
					"wide table", header, colIndx, headers.length));
		}

		ArrayList<String> values = new ArrayList<>(psms.size());

		for (PSM psm : psms) {
			values.add(psm.getColumnValue(header));
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

		if (colIndx < 0 || colIndx >= (headers.length)) {
			throw new ArrayIndexOutOfBoundsException(String.format("Column index %d is out of bounds for a %d column" +
					"wide table", colIndx, headers.length));
		}

		HashMap<String,String> values = new HashMap<>(psms.size());

		for (PSM psm : psms) {
			values.put(psm.spec, psm.getColumnValue(header));
		}

		return values;
	}

	public void save(boolean overwrite) throws IOException {
		String tempFoutName = fname + ".tmp";
		PrintWriter out = new PrintWriter(new FileWriter(tempFoutName));

		// Write lines to .tmp file
		out.println(String.join("\t", headers));
		for (PSM psm : psms) {
			out.println(psm.printLine());
		}

		out.close();

		if (overwrite)
			Files.move(Paths.get(tempFoutName), Paths.get(String.valueOf(fname)), StandardCopyOption.REPLACE_EXISTING);
	}

//	public void writeToPsmCache(File cacheFname, ArrayList<Double> varModMasses) throws FileNotFoundException {
//		PrintWriter out = new PrintWriter(cacheFname);
//
//		/* Write the new header */
//		out.println(String.join("\t", this.headers));
//
//		/* For each line in the file, modify accordingly */
//        for (String datum : this.data) {
//            ArrayList<String> sp = new ArrayList<>(Arrays.asList(datum.split("\t")));
//            //Shift variable mod info into delta mass
//            if (!varModMasses.isEmpty()) {
//                double oldDmass = Double.parseDouble(sp.get(this.dMassCol));
//                double newDmass = oldDmass;
//                double oldPepmass = Double.parseDouble(sp.get(getColumn("Calculated Peptide Mass")));
//                double newPepmass = oldPepmass;
//                String[] oldAssignedMods = sp.get(getColumn("Assigned Modifications")).split(",");
//                ArrayList<String> newAssignedMods = new ArrayList<>(Arrays.asList(oldAssignedMods));
//                ArrayList<Integer> dropMods = new ArrayList<>();
//                for (int j = 0; j < oldAssignedMods.length; j++) {
//                    System.out.println(oldAssignedMods[j]);
//                    String modStr = oldAssignedMods[j];
//                    int p = modStr.indexOf("(");
//                    int q = modStr.indexOf(")");
//                    double modMass = Double.parseDouble(modStr.substring(p + 1, q).trim());
//                    for (Double varModMass : varModMasses) {
//                        if (Math.abs(modMass - varModMass) < 0.001) {
//                            dropMods.add(j);
//                            newDmass += modMass;
//                            newPepmass -= modMass;
//                            break;
//                        }
//                    }
//                }
//                for (int k = dropMods.size() - 1; k > 0; k--)
//                    newAssignedMods.remove(k);
//                // reassign adjusted values
//                sp.set(dMassCol, Double.toString(newDmass));
//                sp.set(getColumn("Calculated Peptide Mass"), Double.toString(newPepmass));
//                sp.set(getColumn("Assigned Modifications"), String.join(", ", newAssignedMods));
//                sp.add(Double.toString(oldDmass));
//                sp.add(Double.toString(oldPepmass));
//                sp.add(String.join(", ", newAssignedMods));
//            }
//            out.println(datum);
//        }
//
//		out.close();
//	}
}

