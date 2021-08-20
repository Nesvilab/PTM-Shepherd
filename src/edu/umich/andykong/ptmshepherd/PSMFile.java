package edu.umich.andykong.ptmshepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.glyco.GlycanCandidate;
import edu.umich.andykong.ptmshepherd.glyco.GlycanResidue;
import edu.umich.andykong.ptmshepherd.glyco.GlycoAnalysis;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;
import edu.umich.andykong.ptmshepherd.peakpicker.PeakAnnotator;

import java.io.*;
import java.lang.reflect.Array;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.zip.CRC32;

public class PSMFile {

	String [] headers;
	public ArrayList<String> data;
	public ArrayList<String> mappedRuns;
	public int dMassCol, precursorCol;
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



	public static void getMappings(File path, HashMap<String,File> mappings) {
		HashMap<String, Integer> datTypes = new HashMap<>();
			datTypes.put("mgf", 4);
			datTypes.put("mzBIN", 3);
			datTypes.put("mzML", 2);
			datTypes.put("mzXML", 2);
			datTypes.put("raw", 1);

		if(path.isDirectory()) {		
			File [] ls = path.listFiles();
			//get mapping for each file
			for(int i = 0; i < ls.length; i++) {
				getMappings(ls[i],mappings);
			}
		} else {
			String [] ns = splitName(path.getName());
			if (ns[0].contains("_calibrated"))
				ns[0] = ns[0].substring(0, ns[0].indexOf("_calibrated"));
			else if (ns[0].contains("_uncalibrated"))
				ns[0] = ns[0].substring(0, ns[0].indexOf("_uncalibrated"));
			if (mappings.containsKey(ns[0]) && (ns[1].equals("mzXML") || ns[1].equals("mzML") || ns[1].equals("raw") || ns[1].equals("mzBIN") || ns[1].equals("mgf"))) {
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
	public void mergeGlycoTable(File glyf, int numColsToUse, boolean writeGlycansToAssignedMods, boolean nGlycan, String allowedResidues) throws Exception {
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

		/* Find headers, dynamically detect columns */
		int observedModCol = getColumn("Observed Modifications");
		int assignedModCol = getColumn("Assigned Modifications");
		int fraggerLocCol = getColumn("MSFragger Localization");
		int peptideCol = getColumn("Peptide");
		if (observedModCol == -1) {
			// did not find correct column! Use default and warn user
			System.out.println("Warning: Could not find Observed Modifications column in PSM table, using default insert point");
			observedModCol = 27;
		}
		if (assignedModCol == -1 && writeGlycansToAssignedMods) {
			PTMShepherd.print("ERROR: Could not find Assigned Modifications column(s) in PSM table. Glycans NOT written to Assigned Modifications");
			writeGlycansToAssignedMods = false;
		}
		if (fraggerLocCol == -1 && writeGlycansToAssignedMods) {
			PTMShepherd.print("ERROR: Could not find MSFragger localization column in PSM table. Make sure you are using Philosopher 4.0.0+ and have delta mass localization enabled in MSFragger. Glycans NOT written to Assigned Modifications");
			writeGlycansToAssignedMods = false;
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

		ArrayList<String> newHeaders = new ArrayList<>();
		for (int i = 0; i < this.headers.length; i++)
			newHeaders.add(this.headers[i]);
		// add after observed mod column (observedModCol + 1)
		boolean hasPreviousGlycoInfo = hasGlycanAssignmentsWritten();
		if (!hasPreviousGlycoInfo) {
			// do not add glyco headers if the PSM table already has them
			newHeaders.addAll(observedModCol + 1, Arrays.asList(glyHeaders).subList(mergeFromCol + 1, mergeFromCol + numColsToUse));
		}
		out.println(String.join("\t", newHeaders));

		/* Match glycolines on PSM spectrum keys */
		int pSpecCol = getColumn("Spectrum");
		for (String cpline : this.data) {
			ArrayList<String> newLine = new ArrayList<>(Arrays.asList(cpline.split("\t")));
			String pSpec = newLine.get(pSpecCol);
			ArrayList<String> glyLine = new ArrayList<>(Arrays.asList(glyLines.get(pSpec)));

			// insert into and after Observed Modifications column
			if (hasPreviousGlycoInfo) {
				// only replace observed mods if a glycan was found in a given line
				if (glyLine.get(mergeFromCol).length() > 0){
					// replace old glyco info with the new, including the observed mods column if a glycan is present
					for (int i=0; i < numColsToUse; i++) {
						newLine.set(observedModCol + i, glyLine.get(mergeFromCol + i));		// existing info is at obs mod column, new info is in glyline at mergeFromCol
					}
					if (writeGlycansToAssignedMods) {
						newLine = getGlycanAssignedMod(newLine, observedModCol, nGlycan, allowedResidues, assignedModCol, fraggerLocCol, peptideCol);
					}
				}
			} else {
				// overwrite the observed mods column if a glycan is present, then add the two new columns for scores
				if (glyLine.get(mergeFromCol).length() > 0){
					newLine.set(observedModCol, glyLine.get(mergeFromCol));
					if (writeGlycansToAssignedMods) {
						newLine = getGlycanAssignedMod(newLine, observedModCol, nGlycan, allowedResidues, assignedModCol, fraggerLocCol, peptideCol);
					}
				}
				newLine.addAll(observedModCol + 1, glyLine.subList(mergeFromCol + 1, mergeFromCol + numColsToUse));
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
	 * Converts glycan ID to mass and position for Assigned Mods (may change to string for quant later). Uses
	 * MSFragger localization string from Philosopher (4.0.0+) if that contains localization info OR places
	 * glycan on the first allowed position if localization is ambiguous.
	 * NOTE: must come AFTER reading glycan composition from the rawglyco file
	 */
	public ArrayList<String> getGlycanAssignedMod(ArrayList<String> newLine, int observedModCol, boolean nGlycan, String allowedResidues, int assignedModCol, int fraggerLocCol, int peptideCol) {
		// parse glycan composition from recently edited observed mods col
		TreeMap<GlycanResidue, Integer> glycanComp;
		String observedMods = newLine.get(observedModCol);
		try {
			glycanComp = PTMShepherd.parseGlycanString(observedMods);
		} catch (Exception ex) {
			// Not a glycan OR failedFDR/decoy (PTM-S or Philosopher may put other string formats here). In all cases, ignore
			// todo: not great design - should have better checking to make sure legit glycan parsing didn't fail
			return newLine;
		}
		double glycanMass = GlycanCandidate.computeMonoisotopicMass(glycanComp);

		// Get glycan location from lower case position
		int glycanLocation = -1;
		String fraggerPepLocStr = newLine.get(fraggerLocCol);
		String glycanAA;
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
				glycanLocation = GlycoAnalysis.findNGlycSequon(newLine.get(peptideCol));
			} else {
				// find first allowed residue if not NGlycan mode
				boolean[] allowedPos = SiteLocalization.parseAllowedPositions(newLine.get(peptideCol), allowedResidues);
				for (int i = 0; i < allowedPos.length; i++) {
					if (allowedPos[i]) {
						glycanLocation = i;
						break;
					}
				}
			}
		}
		glycanAA = fraggerPepLocStr.substring(glycanLocation, glycanLocation + 1).toUpperCase();

		// write mass and location to Assigned Mods
		String currentAssignedMods = newLine.get(assignedModCol);
		String glycanMod = String.format("%d%s(%.4f)", glycanLocation + 1, glycanAA, glycanMass);	// site is 1-indexed in PSM table, not 0-indexed
		// check if glycanAssignedMod already present and avoid double adding if so (in case of re-runs on same file)
		boolean alreadyAddedGlycan = false;
		String[] modSplits = currentAssignedMods.split(",");
		for (String mod: modSplits) {
			if (mod.equals(glycanMod)){
				alreadyAddedGlycan = true;
				break;
			}
		}
		String newAssignedMods = alreadyAddedGlycan ? currentAssignedMods : currentAssignedMods.length() == 0 ? glycanMod : currentAssignedMods + ", " + glycanMod;
		newLine.set(assignedModCol, newAssignedMods);
		return newLine;
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
}
