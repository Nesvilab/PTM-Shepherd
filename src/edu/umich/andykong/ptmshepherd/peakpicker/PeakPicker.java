package edu.umich.andykong.ptmshepherd.peakpicker;

import java.io.*;
import java.util.*;

public class PeakPicker {

	public double [][] peaks;
	public double [] peakCenters;
	public final static double c13Mass = 1.003355; //1.00235
	
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
	public void pickPeaks(double [] offsets, double [] sum, double promRatio, double peakWidth, double peakBackground, int nBins, String massOffsets, String isotopes, double precursorTol) throws Exception {

		boolean offsetMode = false;
		double[] mos = new double[0];
		if (!massOffsets.equals("") & (!massOffsets.equals("None"))) {
			offsetMode = true;
		}

		if (offsetMode) {
			String[] tmpIsos; //temp str to hold isotopes
			String[] tmpMos = massOffsets.split("/"); //temp string to hold mass offsets
			if (isotopes.equals("")) {
				tmpIsos = new String[]{"0"};
			} else {
				tmpIsos = isotopes.split("/");
			}
			mos = new double[tmpMos.length*tmpIsos.length+1]; //list of mass offsets
			mos[0] = 0.0; //include unmodified peps
			for (int i = 0; i < tmpMos.length; i++) {
				for (int j = 0; j < tmpIsos.length; j++) {
					double mo = Double.parseDouble(tmpMos[i]) + Integer.parseInt(tmpIsos[j])*c13Mass;
					mos[i*tmpIsos.length+j+1] = mo;
				}
			}
		}


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
					dp.key = isum / icount - osum / ocount;
					dp.value = offsets[i];
					dp.sum = isum;
					dps.add(dp);
				}
			}
		}
		Collections.sort(dps);

		//filter out peaks that don't match offset
		if (offsetMode) {
			//if a peak isnt within precursor_tol range of a mass offset, prep it to be removed
			ArrayList<Integer> removeIs = new ArrayList<>();
			for (int i = 0; i < dps.size(); i++) {
				boolean keepDp = false;
				for (double mo : mos) {
					if ((mo - precursorTol <= dps.get(i).value) && (dps.get(i).value <= mo + precursorTol)) {
						keepDp = true;
					}
				}
				if (!keepDp) {
					removeIs.add(i);
				}
			}
			//iterate backwards over removeIs list to remove indices that don't correspond to a mass offset
			for (int i = removeIs.size()-1; i >= 0; i--) {
				int iToRemove = removeIs.get(i);
				dps.remove(iToRemove);
			}
		}

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
