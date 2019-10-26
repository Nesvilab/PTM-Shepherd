package edu.umich.andykong.ptmshepherd.peakpicker;

import java.io.*;
import java.util.*;

public class PeakPicker {

	public double [][] peaks;
	public double [] peakCenters;
	
	public static double nzr(double a, double b) {
		if(b == 0)
			return 0;
		return (a/b);
	}
	
	public void writeTSV(File f) throws Exception {
		try (PrintWriter out = new PrintWriter(new FileWriter(f))) {
			for (double[] peak : peaks) {
				for (int j = 0; j < 3; j++) {
					if (j != 0)
						out.print("\t");
					out.printf("%.8f", peak[j]);
				}
				out.println();
			}
		}
	}

	//inputs: bin divisions, bin weights, input prominence ratio, peakwidth, peakbackground, number of top peaks to report
	public void pickPeaks(double [] offsets, double [] sum, double promRatio, double peakWidth, double peakBackground, int nBins) throws Exception {
		//inner double pair class
		class DPair implements Comparable<DPair> {
			double key, value;
			int pos;
			double sum;
			public int compareTo(DPair arg0) {
				return -1* Double.compare(key, arg0.key);
			}
		}

		//calculates prominence of every bin
		double [] prom = Prominence.computeProminence(sum);
		ArrayList<DPair> dps = new ArrayList<DPair>();
		for(int i = 0; i < prom.length; i++) {
			int p;
			double cpr = nzr(prom[i],sum[i]);
			if(cpr > promRatio) {
				double ncpr = 0;
				double mh = 0;
				double isum = sum[i], osum = 0;
				int icount = 1, ocount = 0; //inner sum/outer sumer
				p = i-1;
				while(p >= 0 && Math.abs(offsets[p]-offsets[i]) <= peakBackground) {
					mh = Math.max(mh,sum[p]);
					ncpr = Math.max(ncpr,nzr(prom[p],sum[p]));
					if(Math.abs(offsets[p]-offsets[i]) <= peakWidth) {
						isum += sum[p];
						icount++;
					}
					else {
						osum += sum[p];
						ocount++;
					}
					p--;
				}
				p = i+1;
				while(p < offsets.length && Math.abs(offsets[p]-offsets[i]) <= peakBackground) {
					mh = Math.max(mh,sum[p]);
					ncpr = Math.max(ncpr,nzr(prom[p],sum[p]));
					if(Math.abs(offsets[p]-offsets[i]) <= peakWidth) {
						isum += sum[p];
						icount++;
					}
					else {
						osum += sum[p];
						ocount++;
					}
					p++;
				}

				if(cpr >= (ncpr-1e3) && sum[i] >= (mh-1e-3)) {
					DPair dp = new DPair();
					dp.pos = i;
					dp.key = isum/icount - osum/ocount;
					dp.value = offsets[i];
					dp.sum = isum;
					dps.add(dp);
				}
			}
		}
		Collections.sort(dps);
		
		nBins = Math.min(dps.size(),nBins);
		
		peaks = new double[nBins][3];
		for(int i = 0; i < nBins; i++) {
			peaks[i][0] = dps.get(i).value;
			peaks[i][1] = sum[dps.get(i).pos];
			peaks[i][2] = prom[dps.get(i).pos];
			//System.out.println(dps.get(i).sum);
		}

	}
	

}
