package edu.umich.andykong.ptmshepherd.peakpicker;

import java.io.*;
import java.util.*;

import edu.umich.andykong.ptmshepherd.PSMFile;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class PeakSummary {
	private static final Logger log = LoggerFactory.getLogger(PeakSummary.class);

	ArrayList<PeakFeature> features;
	PeakFeature topFeature;
	TreeMap<String,int [][]> counts;
	TreeMap<String,Integer> dsSize;
	int minPsms;
	boolean prepForIonQuant;
	
	public PeakSummary(File peakTSV, int precursorUnits, double pt, String massOffsets, int psmFilter) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(peakTSV));
		features = new ArrayList<>();
		counts = new TreeMap<>();
		dsSize = new TreeMap<>();
		boolean offsetMode = massOffsets.contains("/");
		double precursorTol = pt;
		this.minPsms = psmFilter;

		String cline;
		int cnt = 0;
		while((cline = in.readLine())!= null) {
			String [] sp = cline.split("\t");
			features.add(new PeakFeature(Double.parseDouble(sp[0]), Double.parseDouble(sp[3]), cnt++));
		}
		in.close();

		if (!features.isEmpty()) {
			topFeature = features.get(0);
		} else {
			log.warn("Empty peak features encountered, error in peakpicking. Exiting program.");
			System.exit(1);
		}
		
		Collections.sort(features);
		for (int i = 0; i < features.size(); i++) {
			double left = -1e10;
			double right = 1e10;
			if(i != 0)
				left = features.get(i-1).peakCenter;
			if(i != (features.size()-1))
				right = features.get(i+1).peakCenter;
			if(precursorUnits == 1) { //ppm units
				precursorTol = ((1500.0 + features.get(i).peakCenter) / 1000000.0) * pt;
				features.get(i).peakLower = Math.max((left + features.get(i).peakCenter) / 2, features.get(i).peakCenter - precursorTol);
				features.get(i).peakUpper = Math.min((right + features.get(i).peakCenter) / 2, features.get(i).peakCenter + precursorTol);
			} else {
				features.get(i).peakLower = Math.max((left + features.get(i).peakCenter) / 2, features.get(i).peakCenter - precursorTol);
				features.get(i).peakUpper = Math.min((right + features.get(i).peakCenter) / 2, features.get(i).peakCenter + precursorTol);
			}
		}
		for(int i = 0; i < features.size() - 1; i++) {
			if(features.get(i).peakUpper == features.get(i+1).peakLower)
				features.get(i).peakUpper -= 1e-6;
		}
		
//		Collections.sort(features,new Comparator<PeakFeature>() {
//			public int compare(PeakFeature o1, PeakFeature o2) {
//				return Integer.valueOf(o1.order).compareTo(o2.order);
//			}
//		});
	}
	
	static double nzr(double a, double b) {
		if(b == 0)
			return 0;
		return a / b;
	}
	
	public static double [] readPeakCenters(File f) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(f));
		String cline;
		ArrayList<String> vals = new ArrayList<>();
		//Skip header line
		in.readLine();
		while((cline = in.readLine())!= null) {
			String [] sp = cline.split("\t");
			vals.add(sp[0]);
		}
		in.close();
		
		double [] res = new double[vals.size()];
		for(int i = 0; i < res.length; i++)
			res[i] = Double.parseDouble(vals.get(i));
		return res;
	}

	public static double [][] readPeakBounds(File f) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(f));
		String cline;
		ArrayList<String[]> vals = new ArrayList<>();
		//Skip header line
		in.readLine();
		while((cline = in.readLine())!= null) {
			String [] sp = cline.split("\t");
			vals.add(new String[]{sp[0], sp[1], sp[2]});
		}
		in.close();

		double [][] res = new double[3][vals.size()];
		for(int i = 0; i < 3; i++) {
			for (int j = 0; j < res[0].length; j++)
				res[i][j] = Double.parseDouble(vals.get(j)[i]);
		}

		return res;
	}

	public void writeTSVSummary(File f) throws Exception {
		int nInsuffPeaks = 0;

		String [] exps = new String[counts.size()];
		int cnt = 0;
		for(String ds : counts.keySet())
			exps[cnt++] = ds;
		
		PrintWriter out = new PrintWriter(new FileWriter(f));
		
		out.print("peak_apex\tpeak_lower\tpeak_upper");
		out.print("\t"+"peak_signal");

		for(int i = 0; i < exps.length; i++)
			out.print("\t"+exps[i] + "_PSMs");
		for(int i = 0; i < exps.length; i++)
			out.print("\t"+exps[i] + "_percent_PSMs");
		for(int i = 0; i < exps.length; i++)
			out.print("\t"+exps[i] + "_peptides\t" + exps[i] + "_percent_also_in_unmodified");
		out.print("\t"+"percent_also_in_unmodified");
		out.print("\t"+"PSMs");
		out.println();

		for(int i = 0; i < features.size(); i++) {
			int totPsm = 0;
			int pt = -1;
			for(int j = 0; j < features.size() && pt == -1; j++)
				if(features.get(j).order == i)
					pt = j;
			//purge peaks with too low PSM count
			for(int j = 0; j < exps.length; j++) //count PSMs
				totPsm += counts.get(exps[j])[pt][1];
			if (totPsm < minPsms) {
				nInsuffPeaks += 1;
				continue;
			}
			/** Get and print counts to file */
			out.printf("%.5f\t%.5f\t%.5f\t%.2f", features.get(pt).peakCenter,features.get(pt).peakLower,
					features.get(pt).peakUpper, (features.get(pt).snr / 1000000) * 100);
			for(int j = 0; j < exps.length; j++) //count PSMs
				out.printf("\t%d", counts.get(exps[j])[pt][1]);
			for(int j = 0; j < exps.length; j++) //count percent psms
				out.printf("\t%.3f", (100.0*counts.get(exps[j])[pt][1])/dsSize.get(exps[j]));
			for(int j = 0; j < exps.length; j++) //count peps
				out.printf("\t%d\t%.2f", counts.get(exps[j])[pt][0],100*nzr(counts.get(exps[j])[pt][2],counts.get(exps[j])[pt][0]));
			double weightedPeps = 0;
			double totPep = 0;
			for(int j = 0; j < exps.length; j++){
				weightedPeps += counts.get(exps[j])[pt][0]*nzr(counts.get(exps[j])[pt][2],counts.get(exps[j])[pt][0]);
				totPep += counts.get(exps[j])[pt][0];
			}
			out.printf("\t%.2f", 100*(weightedPeps/totPep));
			out.printf("\t%d", totPsm);
			out.println();
		}
		out.close();

		if (minPsms > 0)
			System.out.printf("\tRemoved %d peaks with insufficient PSMs\n", nInsuffPeaks);
	}
	
	public void reset() {
		for(int i = 0; i < features.size(); i++)
			features.get(i).reset();
	}
	
	public void commit(String dsName, int sz) {
		int [][] cnts = new int[features.size()][3];
		for(int i = 0; i < features.size(); i++) {
			int countZeroPeps = 0;
			if (topFeature != null) {
				for (String cp : features.get(i).peps)
					if (topFeature.peps.contains(cp))
						countZeroPeps++;
			}
			cnts[i][0] = features.get(i).peps.size();
			cnts[i][1] = features.get(i).psms;
			cnts[i][2] = countZeroPeps;
		}
		dsSize.put(dsName, sz);
		counts.put(dsName, cnts);
	}
	
	public void appendPSMs(PSMFile pf) {
		int seqcol = pf.getColumn("Peptide");
		int mdcol = pf.dMassCol;
//		long stime = System.currentTimeMillis();
		for(int i = 0; i < pf.data.size(); i++) {
			String [] sp = pf.data.get(i).split("\t");
			double md = Double.parseDouble(sp[mdcol]);

//			for(int j = 0; j < features.size(); j++) {
//				if(features.get(j).peakLower <= md && features.get(j).peakUpper >= md) {
//					features.get(j).peps.add(sp[seqcol]);
//					features.get(j).psms++;
//					break;
//				}
//			}
			if (topFeature != null) {
				if (md >= topFeature.peakLower && md <= topFeature.peakUpper) {
					topFeature.peps.add(sp[seqcol]);
					topFeature.psms++;
				} else {
					PeakFeature fast = PeakFeature.getMatchedFeature(features, md);
					if (fast != null) {
						fast.peps.add(sp[seqcol]);
						fast.psms++;
					}
				}
			}
		}
//		System.out.println(pf + " " + (System.currentTimeMillis()-stime));
	}

}
