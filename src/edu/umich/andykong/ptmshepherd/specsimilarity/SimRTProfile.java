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
	double [][] peaks;
	double peakTol;
	int precursorUnits;
	
	public SimRTProfile(double [][] peakVals, double peakTol, int precursorUnits) {
		this.masses = new double[peakVals[0].length];
		for (int i = 0; i < peakVals[0].length; i++) {
			this.masses[i] = peakVals[0][i];
			//System.out.println(this.masses[i]);
		}

		peaks = peakVals;

		this.peakTol = peakTol;
		this.precursorUnits = precursorUnits;
		locate = new FastLocator(peaks, peakTol, precursorUnits);
		records = new SimRTRecord[masses.length];
		for(int i = 0; i < masses.length; i++)
			records[i] = new SimRTRecord(masses[i], i);
	}
	
	public void writeProfile(String path) throws Exception {
		PrintWriter out = new PrintWriter(new FileWriter(path));
		out.printf("%s\t%s\t%s\t%s\t%s\t%s\n",
				"peak","PSMs","similarity","similarity_(variance)","rt_shift", "rt_shift_(variance)");
		for(int i = 0; i < records.length; i++) {
			out.println(records[i].toString());
		}
		out.close();
	}

}
