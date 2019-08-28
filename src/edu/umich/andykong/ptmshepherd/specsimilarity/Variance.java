package edu.umich.andykong.ptmshepherd.specsimilarity;

public class Variance {

	int N;
	double mean, M2;
	
	public Variance() {
		N = 0;
		mean = M2 = 0;
	}
	
	public void update(double v) {
		N++;
		double delta = v - mean;
		mean += (delta / N);
		double delta2 = v - mean;
		M2 += delta*delta2;
	}
	
	public double getMean() {
		return mean;
	}
	
	public double getVariance() {
		if(N < 2)
			return 0;
		return (M2 / (N-1));
	}
	
}
