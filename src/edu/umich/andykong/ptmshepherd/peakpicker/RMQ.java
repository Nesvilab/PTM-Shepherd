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

import java.util.Arrays;


public class RMQ {

	double [] vals;
	int BLOCKSZ;
	RMQLinearithmic lmql;
	
	public RMQ(double [] v) {
		vals = v;
		BLOCKSZ = (int)((Math.log(v.length)/Math.log(2))/4)+1;
		double [] blockMins = new double[vals.length / BLOCKSZ + 1];
		Arrays.fill(blockMins, 1e100);
		for(int i = 0; i < v.length; i++) //
			blockMins[i / BLOCKSZ] = Math.min(blockMins[i / BLOCKSZ], v[i]); //the smallest value from adjacent BLOCKSZ bins
		lmql = new RMQLinearithmic(blockMins);
	}
	
	public double query(int l, int r) {
		double res = 1e100;
		int a = l / BLOCKSZ + 1;
		int b = r / BLOCKSZ - 1;
		for(int i = l; i <= r && i < a*BLOCKSZ; i++) //find minimum value in block for fuzzy vals
			res = Math.min(res,vals[i]);
		for(int i = Math.max(l,(b+1)*BLOCKSZ); i <= r; i++) //
			res = Math.min(res,vals[i]);
		if(b >= a)
			res = Math.min(res,lmql.query(a, b));
		return res;
	}
	
}
