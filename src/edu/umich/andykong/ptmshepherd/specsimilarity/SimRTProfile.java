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
	boolean calcIntensity;
	
	public SimRTProfile(double [][] peakVals, double peakTol, int precursorUnits, boolean calcIntensity) {
		this.masses = new double[peakVals[0].length];
		for (int i = 0; i < peakVals[0].length; i++) {
			this.masses[i] = peakVals[0][i];
		}

		this.calcIntensity = calcIntensity;
		peaks = peakVals;

		this.peakTol = peakTol;
		this.precursorUnits = precursorUnits;
		locate = new FastLocator(peaks, peakTol, precursorUnits);
		records = new SimRTRecord[masses.length];
		for(int i = 0; i < masses.length; i++)
			records[i] = new SimRTRecord(masses[i], i, calcIntensity);
	}
	
	public void writeProfile(String path) throws Exception {
		PrintWriter out = new PrintWriter(new FileWriter(path));
		if (calcIntensity) {
			out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					"peak", "PSMs", "similarity", "similarity_(variance)", "rt_shift", "rt_shift_(variance)", "int_log2fc", "int_log2fc_(variance)");
		} else {
			out.printf("%s\t%s\t%s\t%s\t%s\t%s\n",
					"peak", "PSMs", "similarity", "similarity_(variance)", "rt_shift", "rt_shift_(variance)");
		}
		for(int i = 0; i < records.length; i++) {
			out.println(records[i].toString());
		}
		out.close();
	}

}
