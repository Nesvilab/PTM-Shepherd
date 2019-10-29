package edu.umich.andykong.ptmshepherd.specsimilarity;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;

import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.localization.LocalizationRecord;

public class SimRTProfile {
	public SimRTRecord [] records;
	
	public FastLocator locate;
	double [] masses;
	double peakTol;
	
	public SimRTProfile(double [] massRanges, double peakTol) {
		masses = Arrays.copyOf(massRanges, massRanges.length);
		this.peakTol = peakTol;
		
		locate = new FastLocator(masses, peakTol);
		records = new SimRTRecord[masses.length];
		for(int i = 0; i < masses.length; i++)
			records[i] = new SimRTRecord(masses[i], i);
	}
	
	public void writeProfile(String path) throws Exception {
		PrintWriter out = new PrintWriter(new FileWriter(path));
		out.printf("%s\t%s\t%s\t%s\t%s\n",
				"Peak","Matched PSMs","Similarity (mean)","Similarity (variance)","DeltaRT_Stat (mean/variance)");
		for(int i = 0; i < records.length; i++) {
			out.println(records[i].toString());
		}
		out.close();
	}

}
