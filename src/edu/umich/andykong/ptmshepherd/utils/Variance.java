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

package edu.umich.andykong.ptmshepherd.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.stream.Collectors;

public class Variance {

	int N;
	double mean, M2, median;
	ArrayList<Double> vals;
	
	public Variance() {
		N = 0;
		mean = M2 = median = 0;
		vals = new ArrayList<>();
	}
	
	public void update(double v) {
		N++;
		vals.add(v);
		//double delta = v - mean;
		//mean += (delta / N);
		//double delta2 = v - mean;
		//M2 += delta*delta2;
	}

	public void update(ArrayList<Double> vs) {
		this.N += vs.size();
		vals.addAll(vs);
	}
	
	public double getMean() {
		if (N == 0)
			return 0;
		for (double v : vals){
			mean += v;
		}
		mean /= N;
		return mean;
	}

	public void logTransform() {
		this.vals = this.vals.stream()
				.map(d -> (Math.log(d) / Math.log(2)))
				.filter(d -> (Double.isFinite(d)))
				.collect(Collectors.toCollection(ArrayList::new));
	}

	public double getMedian() {
		Collections.sort(this.vals);
		int ind = this.vals.size() / 2;
		if (this.vals.size() % 2 == 1)
			median = this.vals.get(ind);
		else if (this.vals.size() == 0)
			return Double.NaN;
		else
			median = (this.vals.get(ind - 1) + this.vals.get(ind) / 2.0);
		return median;
	}
	
	public double getVariance() {
		if(N < 2)
			return 0;
		getMean();
		for (double v : vals) {
			M2 += ((v - mean)*(v-mean));
		}
		return (M2 / (N-1));
	}
	
}
