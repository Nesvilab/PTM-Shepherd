/*
 *    Copyright 2022 University of Michigan
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

package edu.umich.andykong.ptmshepherd.peakpicker;

//import com.sun.deploy.util.ArrayUtil;

import java.io.*;
import java.sql.Array;
import java.util.*;

public class PeakPicker {

	public double [][] peaks;
	public double [] peakCenters;
	public final static double c13Mass = 1.00235;
	
	public static double nzr(double a, double b) {
		if(b == 0)
			return 0;
		return (a/b);
	}
	
	public void writeTSV(File f) throws Exception {
		try (PrintWriter out = new PrintWriter(new FileWriter(f))) {
			for (double[] peak : peaks) {
				for (int j = 0; j < 4; j++) {
					if (j != 0)
						out.print("\t");
					out.printf("%.8f", peak[j]);
				}
				out.println();
			}
		}
	}

	public double calculatePeakTol(double pepmass, double ppmtol, double modmass){
		double peakTol = ((pepmass + modmass) / 1000000.0) * ppmtol;
		return peakTol;
	}

	//inputs: bin divisions, bin weights, input prominence ratio, peakwidth, peakbackground, number of top peaks to report
	public void pickPeaks(double [] offsets, double [] sum, double promRatio, double peakBackground, int nBins, String massOffsets, String isotopes, int peakUnits, double pw, int precursorUnits, double pt) throws Exception {
		double pepmass = 1500.0;
		double precursorTol = pt;
		double peakWidth = pw;

		boolean offsetMode = false;
		double[] mos = new double[0];
		if(massOffsets.contains("/")){ //todo make this an explicit parameter
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
			double snr;
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
				if(peakUnits == 1) { //ppm
					peakWidth = calculatePeakTol(pepmass, pw, offsets[i]);
					peakBackground = 2.5 * peakWidth;
				}
				//filter out peaks that don't match offset
				if (offsetMode) {
					if(precursorUnits == 1) { //ppm
						precursorTol = calculatePeakTol(pepmass, pt, offsets[i]);
					}
					//if a peak isnt within precursor_tol range of a mass offset, skip it
					for (double mo : mos) {
						if (!((mo - precursorTol <= offsets[i]) && (offsets[i] <= mo + precursorTol))) {
							continue;
						}
					}
				}
				double ncpr = 0;
				double mh = 0;
				double isum = sum[i], osum = 0;
				int icount = 1, ocount = 0; //inner sum/outer sum
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
					dp.key = (isum / icount - osum / ocount) * icount; //old mode excluded icount
					dp.value = offsets[i];
					dp.sum = isum;
					dp.snr = dp.key;// * icount; multiple by icount for original algorithm
					dps.add(dp);
				}
			}
		}

		//this will merge peaks that should not be split if operating in massoffset mode
		if (offsetMode == true) {
			ArrayList<Integer> popDPs = new ArrayList<>();
			Arrays.sort(mos);
			int[] nmos = new int[mos.length]; //holds the max number of peaks that should be present near an mo
			int[] npeaks = new int[mos.length]; //holds the number of peaks that are present near an mo
			//cleaning of peaks is not performed in this loop in case we need to remove peaks in crowded areas at a later point
			for (int i = 0; i < mos.length; i++) {
				double mo = mos[i];
				if (peakUnits == 1) //ppm
					peakWidth = calculatePeakTol(pepmass, pw, mo);
				else
					peakWidth = pw;
				for (int j = 0; j < mos.length; j++) {
					if (Math.abs(mo - mos[j]) <= peakWidth)
						nmos[i]++;
				}
				for (DPair dp : dps) {
					if (Math.abs(mo - dp.value) <= peakWidth)
						npeaks[i]++;
				}
			}

			//remove peaks
			for (int i = 0; i < mos.length; i++) {
				if (nmos[i] == 1) { //if only one mo in the region
					if (npeaks[i] > 1) { //if more than one peak detected
						double mo = mos[i];
						double maxSig = 0.0;
						double maxVal = -10000;
						int maxIndex = -1;
						if (peakUnits == 1) //ppm
							peakWidth = calculatePeakTol(pepmass, pw, mo);
						else
							peakWidth = pw;
						ArrayList<Integer> popDps = new ArrayList<>();
						for (int j = 0; j < dps.size(); j++) {
							if (Math.abs(mo - dps.get(j).value) <= peakWidth)
								popDPs.add(j);
						}
						for (int j = 0; j < popDPs.size(); j++) {
							if (dps.get(j).sum >= maxSig) {
								maxVal = dps.get(j).value;
								maxSig = dps.get(j).sum;
								maxIndex = j;
							}
						}
						popDPs.remove(maxIndex);
						Collections.sort(popDPs);
						for (int j = popDPs.size() - 1; i >= 0; i--) {
							dps.remove(popDPs.get(j));
						}
					}
				}
			}
		}

		Collections.sort(dps);

		nBins = Math.min(dps.size(),nBins);

		peaks = new double[nBins][5];
		for(int i = 0; i < nBins; i++) {
			peaks[i][0] = dps.get(i).value;
			peaks[i][1] = sum[dps.get(i).pos];
			peaks[i][2] = prom[dps.get(i).pos];
			peaks[i][3] = dps.get(i).snr;
			//System.out.println(dps.get(i).sum);
		}

	}

	public double[][] getPeaks() {
		return this.peaks;
	}
}
