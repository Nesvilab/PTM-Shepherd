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


public class RMQLinearithmic {

	double [][] mins;
	
	public RMQLinearithmic(double [] v) {
		int levels = (int)(Math.log(v.length)/Math.log(2) + 1e-9) + 1;
		mins = new double[levels][v.length];
		for(int i = 0; i < levels; i++)
			for(int j = 0; j < v.length; j++) {
				if(i == 0)
					mins[i][j] = v[j];
				else {
					if((j + (1 << (i-1))) < v.length)
						mins[i][j] = Math.min(mins[i-1][j],mins[i-1][j + (1 << (i-1))]);
				}
			}
	}
	
	public double query(int l, int r) {
		int sz = r - l + 1;
		int level = (int)(Math.log(sz)/Math.log(2) + 1e-9);
		return Math.min(mins[level][l],mins[level][r+1 - (1 << level)]);
	}
	
}
