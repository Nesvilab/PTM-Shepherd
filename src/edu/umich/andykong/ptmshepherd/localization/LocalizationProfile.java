package edu.umich.andykong.ptmshepherd.localization;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;

public class LocalizationProfile {

	public LocalizationRecord [] records;
	
	public FastLocator locate;
	double [] masses;
	double [][] peaks; //apex, left, right
	double peakTol;
	int precursorUnits;
	ArrayList<String> globalPepSeqs;
	ArrayList<String> uniqueGlobalPepSeqs;

	static final int [] AAcnts = {3637222,0,1163038,2477586,3690290,1854622,
			  3426107,1356881,2222643,0,2959209,5141499,
			  1134389,1840802,0,3301694,2489112,2946921,
			  4383423,2856550,0,3117149,647263,0,1358211,0}; //this is unnecessary now but an easy check for letters that aren't real

	public LocalizationProfile(double [][] peakVals, double peakTol, int precursorUnits) {
		this.masses = new double[peakVals[0].length];
		for (int i = 0; i < peakVals[0].length; i++)
			this.masses[i] = peakVals[0][i];
		peaks = peakVals;
		this.peakTol = peakTol;
		this.precursorUnits = precursorUnits;

		locate = new FastLocator(peaks, peakTol, precursorUnits);
		records = new LocalizationRecord[masses.length + 1]; // + 1 as dumping ground for unmatched psms

		for(int i = 0; i < masses.length; i++)
			records[i] = new LocalizationRecord(masses[i], i);
		records[masses.length] = new LocalizationRecord(masses[0], masses.length);
	}

	public void collectPeptideSequences() {
		ArrayList<String> globalPepSeqs = new ArrayList<>();
		for (int i = 0; i < records.length; i++)
			globalPepSeqs.addAll(records[i].localPepSeqs);

		ArrayList<String> uniqueGlobalPepSeqs = new ArrayList<>(globalPepSeqs.stream().distinct().collect(Collectors.toList()));

		this.globalPepSeqs = globalPepSeqs;
		this.uniqueGlobalPepSeqs = uniqueGlobalPepSeqs;
	}

	public void writeProfile(String path) throws Exception {
		PrintWriter out = new PrintWriter(new FileWriter(path));
		out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
				"peak","localized_PSMs","PSMs","n-term_localization_rate",
				"AA1","AA1_enrichment_score", "AA1_psm_count",
				"AA2","AA2_enrichment_score", "AA2_psm_count",
				"AA3","AA3_enrichment_score", "AA3_psm_count");
		for(int i = 0; i < 26; i++)
			if(LocalizationProfile.AAcnts[i] != 0)
				out.printf("\t%c_enrichment", 'A' + i);
		out.println();
		collectPeptideSequences();
		for(int i = 0; i < records.length - 1; i++) {
			out.println(records[i].toString(this.globalPepSeqs, this.uniqueGlobalPepSeqs));
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
