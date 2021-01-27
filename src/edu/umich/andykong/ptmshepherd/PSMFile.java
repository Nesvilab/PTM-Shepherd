package edu.umich.andykong.ptmshepherd;
import java.io.*;
import java.util.*;
import java.util.zip.CRC32;

public class PSMFile {

	String [] headers;
	public ArrayList<String> data;
	public ArrayList<String> mappedRuns;
	public static int dMassCol;
	public String prefType;
	
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
			//System.out.println(data.get(i));
			String [] sp = data.get(i).split("\t");
			res.add(Float.parseFloat(sp[col]));
		}
		return res;
	}
	
	public static void getMappings(File path, HashMap<String,File> mappings) {
		HashMap<String, Integer> datTypes = new HashMap<>();
			datTypes.put("mgf", 0);
			datTypes.put("raw", 1);
			datTypes.put("mzBIN", 2);
			datTypes.put("mzXML", 3);
			datTypes.put("mzML", 4);

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

	public PSMFile(File f) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(f), 1 << 22);
		headers = in.readLine().split("\t");
		//find delta mass column for different philosopher versions
		int col = getColumn("Delta Mass");
		if (col == -1) {
			col = getColumn("Adjusted Delta Mass");
		}
		if (col == -1) {
			col = getColumn("Original Delta Mass");
		}
		this.dMassCol = col;
		data = new ArrayList<>();
		String cline;
		while((cline = in.readLine()) != null) {
			if(cline.trim().length() > 0)
				data.add(cline);
		}
		in.close();
	}
}
