package edu.umich.andykong.ptmshepherd.localization;

import java.io.*;
import java.util.*;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;

public class LocalizationProfile {

	public LocalizationRecord [] records;
	
	public FastLocator locate;
	double [] masses;
	double peakTol;

	static final int [] AAcnts = {3637222,0,1163038,2477586,3690290,1854622,
			  3426107,1356881,2222643,0,2959209,5141499,
			  1134389,1840802,0,3301694,2489112,2946921,
			  4383423,2856550,0,3117149,647263,0,1358211,0}; //this is unnecessary now but an easy check for letters that aren't real

	public LocalizationProfile(double [] massRanges, double peakTol) {
		masses = Arrays.copyOf(massRanges, massRanges.length);
		this.peakTol = peakTol;

		locate = new FastLocator(masses, peakTol);
		records = new LocalizationRecord[masses.length];

		for(int i = 0; i < masses.length; i++)
			records[i] = new LocalizationRecord(masses[i], i);
	}
	
	public void writeProfile(String path) throws Exception {
		PrintWriter out = new PrintWriter(new FileWriter(path));
		out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s",
				"Peak","Localized PSMs","Total PSMs","N-term rate",
				"EnrichedAA_1","EnrichedAA_2","EnrichedAA_3");
		for(int i = 0; i < 26; i++)
			if(LocalizationProfile.AAcnts[i] != 0)
				out.printf("\t%c_enrichment", 'A' + i);
		out.println();
		for(int i = 0; i < records.length; i++) {
			out.println(records[i].toString());
		}
		out.close();
	}

	public static void printNormsFromCounts() {
		double [] AAnorm = new double[26];
		double AAsum = 0;
		for(int i = 0; i < 26; i++)
			AAsum += AAcnts[i];
		for(int i = 0; i < 26; i++) {
			AAnorm[i] = AAcnts[i];
			AAnorm[i] /= AAsum;
			if(AAnorm[i] == 0)
				AAnorm[i] = 1;
		}
		for(int i = 0; i < 26; i++)
			System.out.printf("%.5f\n", AAnorm[i]);
	}
	
}
