package edu.umich.andykong.ptmshepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
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
	*  Update: only writes 3 columns (best glycan and score, 2nd best glycan) rather than full table
	*  numColsToUse gives the number of columns to take (intended to be 3)
	*  todo: NOTE - if glyco information has been written to this table already, it will get doubled
	*/
	public void mergeGlycoTable(File glyf, int numColsToUse) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(glyf), 1 << 22);
		String tempFoutName = this.fname + ".glyco.tmp";
		PrintWriter out = new PrintWriter(new FileWriter(tempFoutName));
		String[] glyHeaders = in.readLine().split("\t");

		/* Get glyco data */
		HashMap<String, String[]> glyLines = new HashMap<>();
		String cgline;
		int gSpecCol = 0; //should be dynamically calculated
		while ((cgline = in.readLine()) != null) {
			String[] sp = cgline.split("\t", -1);
			glyLines.put(sp[gSpecCol], sp);
		}
		in.close();

		/* Merge headers */
		int mergeFromCol = 5; // todo this should be dynamically calculated
		int observedModCol = 23; 	// todo this should be dynamically detected?
		ArrayList<String> newHeaders = new ArrayList<>();
		for (int i = 0; i < this.headers.length; i++)
			newHeaders.add(this.headers[i]);
		// add after observed mod column
		newHeaders.addAll(observedModCol, Arrays.asList(glyHeaders).subList(mergeFromCol, mergeFromCol + numColsToUse));
		out.println(String.join("\t", newHeaders));

		/* Match glycolines on PSM spectrum keys */
		int pSpecCol = getColumn("Spectrum");
		for (String cpline : this.data) {
			ArrayList<String> newLine = new ArrayList<>(Arrays.asList(cpline.split("\t")));
			String pSpec = newLine.get(pSpecCol);
			ArrayList<String> glyLine = new ArrayList<>(Arrays.asList(glyLines.get(pSpec)));
			// insert after Observed Modifications column
			newLine.addAll(observedModCol, glyLine.subList(mergeFromCol, mergeFromCol + numColsToUse));
			out.println(String.join("\t", newLine));
		}
		out.close();

		/* Move old file onto new file */
		Files.move(Paths.get(tempFoutName), Paths.get(String.valueOf(this.fname)), StandardCopyOption.REPLACE_EXISTING);
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
