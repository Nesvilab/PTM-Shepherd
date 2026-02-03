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

import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.glyco.GlycanAssignmentResult;
import edu.umich.andykong.ptmshepherd.glyco.GlycoAnalysis;
import edu.umich.andykong.ptmshepherd.glyco.GlycoParams;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;
import umich.ms.glyco.Glycan;
import umich.ms.glyco.GlycanParser;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.regex.Pattern;
import java.util.zip.CRC32;

import static edu.umich.andykong.ptmshepherd.PTMShepherd.normFName;

public class PSMFile {

	public String [] headers;
	public ArrayList<PSM> psms;
    public ArrayList<PSM> glycoRemovedPsms;
    public int dMassCol, precursorCol, assignedModCol, observedModCol, fraggerLocCol, peptideCol, modPeptideCol,
			calcMZcol, peptideCalcMassCol, chargeCol, intensityCol, specCol, msfraggerLocalizationCol, positionScoresCol,
			bestPositionsCol, ionsBestPosCol, scoreBestPositionCol, scoreAllUnshiftedCol, ionsAllUnshiftedCol,
			eValCol, retentionCol, glycanCompCol, glycanScoreCol, glycanQvalCol;

	public int massdiffToVarmod;
	public File fname;
    public static final Pattern massPattern = Pattern.compile("\\(([-.\\d]+)\\)");



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
            res.add(psm.getDMass());
        }
		return res;
	}

	public ArrayList<ArrayList<Float>> getMassDiffsWithVarmods(boolean useAssignedMods) {
		ArrayList<ArrayList<Float>> res = new ArrayList<>();
		for (PSM psm : psms) {
			ArrayList<Float> psmMods = new ArrayList<>();
			psmMods.add(psm.getDMass());
			if (useAssignedMods) {
				for (Mod mod : psm.getAssignedMods()) {
					psmMods.add((float) mod.mass);
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

	public static void getMappings(File path, HashMap<String,File> mappings, HashMap<String,File> originalMappings, HashSet<String> runNames) {
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
		LinkedHashMap<String, Integer> originalPriorities = new LinkedHashMap<>();		// same as priorities but without mzBIN_cache
		for (String key : priorities.keySet()) {
			if (key.equals(".mzBIN_cache")) {
				continue;
			}
			originalPriorities.put(key, priorities.get(key));
		}

		// Recursively search all directories
		if(path.isDirectory()) {
			File [] ls = path.listFiles();
			if (ls == null) {
				return;		// empty directory or IO error. Ignore directory
			}
			// get mapping for each file
            for (File l : ls) {
                getMappings(l, mappings, originalMappings, runNames);
            }
		} else { // see if valid file ext
			tryMapFile(path, mappings, runNames, priorities);
			tryMapFile(path, originalMappings, runNames, originalPriorities);
		}
	}

	private static void tryMapFile(File path, HashMap<String, File> mappings, HashSet<String> runNames, LinkedHashMap<String, Integer> priorities) {
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
	 */
	public void mergeGlycoTable(String datasetName, GlycoParams glycoParams) {
		boolean hasPreviousGlycoInfo = hasGlycanAssignmentsWritten();
        ArrayList<PSM> filteredPSMs = new ArrayList<>();

		/* Match glycolines on PSM spectrum keys */
		ArrayList<Integer> psmKeys = new ArrayList<>();
		ArrayList<String> glycanComps = new ArrayList<>();
		ArrayList<String> glycanScores = new ArrayList<>();
		ArrayList<String> glycanQvals = new ArrayList<>();
		for (PSM psm : psms) {
			psmKeys.add(psm.lineNum);

			// check if a glycan was found
			if (psm.glycanAssignmentResult != null && psm.glycanAssignmentResult.foundGlycan) {
				GlycanAssignmentResult result = psm.glycanAssignmentResult;
				String assignedGlycan;
				String glycanScore = String.format("%.4f", result.glycanScore);
				String glycanQval = String.format("%.6f", result.glycanQval);
				if (result.isDecoyGlycan && !glycoParams.printGlycoDecoys) {
					if (result.bestTarget == null) {
						assignedGlycan = "No target matches";
						glycanScore = "";
						glycanQval = "";
					} else {
						// report best target glycan instead of decoy (q-value will be reported as 1)
						assignedGlycan = result.bestTarget.toString();
						glycanScore = String.format("%.4f", result.bestTarget.glycanScore);
						glycanQval = "1";
					}
				} else {
					assignedGlycan = result.bestCandidate.toString();
				}
				// save glycan info directly or to the lists to add columns to the PSM table later
				if (!hasPreviousGlycoInfo) {
					glycanComps.add(assignedGlycan);
					glycanScores.add(glycanScore);
					glycanQvals.add(glycanQval);
				} else {
					psm.spLine.set(glycanCompCol, assignedGlycan);
					psm.spLine.set(glycanScoreCol, glycanScore);
					psm.spLine.set(glycanQvalCol, glycanQval);
				}
                // update assigned mods column
                boolean failOrDecoy = result.isDecoyGlycan || result.glycanQval >= glycoParams.glycoFDR;
                if (glycoParams.writeGlycansToAssignedMods) {
                    writeGlycanToAssignedMod(psm, assignedGlycan, failOrDecoy, glycoParams);
                }

                // keep glyco PSMs that pass FDR only
                if (!failOrDecoy) {
                    filteredPSMs.add(psm);
                } else {
                    glycoRemovedPsms.add(psm);
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
                filteredPSMs.add(psm);  // keep non-glyco PSMs
			}
		}

		if (!hasPreviousGlycoInfo) {
			// add new glycan columns to the PSM table if not previously written
			addColumn(observedModCol + 1, "Glycan q-value", psmKeys, glycanQvals);
			addColumn(observedModCol + 1, "Glycan Score", psmKeys, glycanScores);
			addColumn(observedModCol + 1, "Total Glycan Composition", psmKeys, glycanComps);
		}

        // save unfiltered PSMs to separate file if glyco FDR filtering is being applied
        if (filteredPSMs.size() != psms.size()) {
            File unfilteredPsmFile = new File(normFName(datasetName + "_unfiltered_psm.tsv"));
            save(unfilteredPsmFile, headers, psms, true);
        }
		save(fname, headers, filteredPSMs, true);

        try {
            assert psms.size() == filteredPSMs.size() + glycoRemovedPsms.size();
        } catch (AssertionError e) {
            PTMShepherd.die(String.format("Error: Mismatch in PSM counts after glyco FDR filtering: %d removed, %d filtered, %d total. Please report this issue to the developers.", glycoRemovedPsms.size(), filteredPSMs.size(), psms.size()));
        }
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
	 * Converts glycan ID to mass and position for Assigned Mods (may change to string for quant later). Uses
	 * MSFragger localization string from Philosopher (4.0.0+) if that contains localization info OR places
	 * glycan on the first allowed position if localization is ambiguous.
	 * Also writes to modified peptide and delta mass columns.
	 * Handles cases where information was previously written to the PSM table by removing/replacing the previous ID if present
	 */
	public void writeGlycanToAssignedMod(PSM psm, String glycanStr, boolean failOrDecoy, GlycoParams glycoParams) {
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

		if (glycanStr.equalsIgnoreCase("no target matches") || glycanStr.equalsIgnoreCase("No Glycan Matched")) {
			return;		// skip, no target glycan info to propagate
		}

		/* Get glycan mass */
		Glycan glyc;
		glyc = GlycanParser.parseGlycanString(glycanStr, glycoParams.glycanResiduesMap);
		if (glyc == null) {
			// try old format in case of old PSM file
			glyc = GlycanParser.parseOldPTMSGlycanString(glycanStr, glycoParams.glycanResiduesMap);
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

		/* write mass and location to Assigned Mods */
        Mod prevGlycanMod = null;
        for (Mod mod: psm.getOriginalAssignedMods()) {
            if (mod.position == glycanLocation + 1) {
                prevGlycanMod = mod;
                break;
            }
        }
        double prevGlycanMass = prevGlycanMod == null? 0 : prevGlycanMod.mass;

		// add the assigned glycan to the updated mod list (from which we removed any old glycan mods) if not failed FDR or is decoy
		if (editPSMGlycoEntry) {
            int originalIndex = psm.getOriginalAssignedMods().indexOf(prevGlycanMod);
            if (originalIndex == -1) {
                // glycan not previously placed, add at appropriate position
                int addIndex = 0;
                for (Mod mod : psm.getAssignedMods()) {
                    if (mod.position < glycanLocation + 1) {
                        addIndex++;
                    } else {
                        break;
                    }
                }
                psm.getAssignedMods().add(addIndex, new Mod(glycanLocation + 1, glycanMass));
            } else {
                psm.getAssignedMods().add(originalIndex, new Mod(glycanLocation + 1, glycanMass));
            }
        }
		psm.spLine.set(assignedModCol, psm.printAssignedMods());

		psm.editModifiedPeptide(glycanLocation + 1, glycanMass, massdiffToVarmod, modPeptideCol);

		/* Update delta mass AND calc m/z columns */
		if (glycoParams.removeGlycanDeltaMass) {
			psm.updateDeltaMass(glycanMass, massdiffToVarmod, prevGlycanMass, peptideCalcMassCol, calcMZcol, dMassCol, assignedModCol);
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
	public void preparePsmTableForIonQuant(double[][] peakBounds, int precUnits, double precTol) {
		/* Check to make sure IonQuant column doesnt already exist */
		if (!(getColumn("Theoretical Modification Mass") == -1)) { /* Already exists */
			System.out.printf("\tPSM table at %s already IonQuant compatible\n",this.fname);
			return;
		}

		FastLocator locator = new FastLocator(peakBounds, precTol, precUnits);
		String tempFoutName = fname + ".iq.tmp";
        try {
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
        } catch (IOException e) {
            PTMShepherd.die("error writing IonQuant compatible PSM file: " + fname + "\n" + e.getMessage());
        }

        fname.delete();
		File newFileName = new File(tempFoutName);
		newFileName.renameTo(fname);
	}

	/**
	 * Reads all lines into PSMs in the psms array and stores the scan to line index mapping in scanToLineMap.
	 */
	public PSMFile(File f, int massdiffToVarmod) {
		psms = new ArrayList<>();
        glycoRemovedPsms = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new FileReader(f), 1 << 22);
			fname = f;
			headers = in.readLine().split("\t");
			initColumns();

			if (massdiffToVarmod > 0 && msfraggerLocalizationCol == -1) {
				PTMShepherd.print("Error: Delta masses were removed by MSFragger but localization was not performed. This mode is not supported. Please either enable localize_delta_mass or disable mass_diff_to_variable_mod in MSFragger and try again.");
				PTMShepherd.die("Invalid MSFragger mass-diff-to-variable-mod configuration");
			}
			if (glycanScoreCol != -1) {
				// handle previous glyco run removing mass diffs
				massdiffToVarmod = 1;
			}
			this.massdiffToVarmod = massdiffToVarmod;

			int i = 0;
			String cline;
			while ((cline = in.readLine()) != null) {
				if (!cline.isEmpty()) {
					PSM thisPSM = new PSM(i, cline, this.massdiffToVarmod, specCol, peptideCol, modPeptideCol, chargeCol, peptideCalcMassCol, dMassCol, assignedModCol, msfraggerLocalizationCol);
					psms.add(thisPSM);
					i++;
				}
			}
			in.close();
		} catch (IOException e) {
			PTMShepherd.die("error reading PSM file: " + f.getName() + "\n" + e.getMessage());
		}
	}

    public void annotateMassDiffs(String [] annotations) {
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
		try {
			PrintWriter out = new PrintWriter(new FileWriter(tempFoutName));

			/* write the new header */
			out.println(String.join("\t", this.headers));
			/* write file lines */
			for (String newLine : newLines) out.println(newLine);

			/* close and rename temp file */
			out.close();
		} catch (IOException e) {
			PTMShepherd.die("error writing annotated mass diffs file: " + this.fname + "\n" + e.getMessage());
		}
		this.fname.delete();
		File newFileName = new File(tempFoutName);
		newFileName.renameTo(this.fname);
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

	public void addColumn(int colIndx, String newHead, ArrayList<Integer> psmLineNumbers, ArrayList<String> vals) {
		// Check that PSM table editing will not fail
		if (vals.size() != psms.size() || psmLineNumbers.size() != psms.size()) {
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
			for (int i = 0; i < psmLineNumbers.size(); i++) {
				psms.get(psmLineNumbers.get(i)).addValAtColumn(colIndx, vals.get(i));
			}
		} else {
			for (int i = 0; i < psmLineNumbers.size(); i++) {
				psms.get(psmLineNumbers.get(i)).replaceValAtColumn(colIndx, vals.get(i));
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
			values.add(psm.getColumnValue(colIndx));
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
			values.put(psm.getSpec(), psm.getColumnValue(colIndx));
		}

		return values;
	}

	public ArrayList<Integer> getPSMlineNums() {
		ArrayList<Integer> indices = new ArrayList<>(psms.size());
		for (PSM psm : psms) {
			indices.add(psm.lineNum);
		}
		return indices;
	}

	public static void save(File fname, String[] headers, ArrayList<PSM> psms, boolean overwrite) {
		String tempFoutName = fname + ".tmp";
		try {
			PrintWriter out = new PrintWriter(new FileWriter(tempFoutName));

			// Write lines to .tmp file
			out.println(String.join("\t", headers));
			for (PSM psm : psms) {
				out.println(psm.printLine());
			}
			out.close();

			if (overwrite) {
				Files.move(Paths.get(tempFoutName), Paths.get(String.valueOf(fname)), StandardCopyOption.REPLACE_EXISTING);
			}

		} catch (IOException e) {
			PTMShepherd.die("error writing PSM file: " + fname + "\n" + e.getMessage());
		}
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

