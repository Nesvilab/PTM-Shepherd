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

public class Prominence {

	//takes histo heights as argument
	public static double [] computeProminence(double [] v) {
		double [] pv = new double[v.length + 2];
		int [] qPos = new int[v.length + 2];
		double [] res = new double[v.length];
		Arrays.fill(res,1e100);
		int pos = 0;
		pv[0] = pv[pv.length-1] = 1e100; //sets terminal values of pv
		for(int i = 0; i < v.length; i++) 
			pv[i+1] = v[i] + (Math.random() / 100000.0); //sets pv to be noisy v
		RMQ rmq = new RMQ(pv); //takes pv (noisy v)
		
		qPos[pos++] = 0;
		for(int i = 1; i < pv.length - 1; i++) { //for bin in noisy heights length
			while(pv[qPos[pos-1]] <= pv[i])  //while noisy histo heights[pos] noisy histo height[i], move cursor left
				pos--;
			res[i-1] = Math.min(res[i-1],pv[i] - rmq.query(qPos[pos-1], i));
			qPos[pos++] = i;
		}
		pos = 0;
		qPos[pos++] = pv.length - 1;
		for(int i = pv.length - 2; i > 0; i--) {
			while(pv[qPos[pos-1]] <= pv[i])
				pos--;
			res[i-1] = Math.min(res[i-1],pv[i] - rmq.query(i, qPos[pos-1]));
			qPos[pos++] = i;
		}
		return res;
	}
}
