package edu.umich.andykong.ptmshepherd.specsimilarity;

import java.util.ArrayList;

public class Variance {

	int N;
	double mean, M2;
	ArrayList<Double> vals;
	
	public Variance() {
		N = 0;
		mean = M2 = 0;
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
	
	public double getMean() {
		if (N == 0)
			return 0;
		for (double v : vals){
			mean += v;
		}
		mean /= N;
		return mean;
	}
	
	public double getVariance() {
		if(N < 2)
			return 0;
		for (double v : vals) {
			M2 += ((v - mean)*(v-mean));
		}
		return (M2 / (N-1));
	}
	
}
