package edu.umich.andykong.ptmshepherd.peakpicker;

import java.util.*;

public class PeakFeature implements Comparable<PeakFeature> {
	
	double peakCenter, peakLower, peakUpper, snr;
	
	HashSet<String> peps;
	int psms, order;

	public PeakFeature(double peakCenter, double snr, int order) {
		this.peakCenter = peakCenter;
		this.order = order;
		this.snr = snr;
		peps = new HashSet<>();
		psms = 0;
	}
	
	public static PeakFeature getMatchedFeature(ArrayList<PeakFeature> features, double v) {
		int lb = 0, ub = features.size() - 1;
		while(ub > lb) {
			int mid = (lb + ub) / 2;
			if(v > features.get(mid).peakUpper)
				lb = mid+1;
			else if(v < features.get(mid).peakLower)
				ub = mid-1;
			else
				return features.get(mid);
		}
		if(v >= features.get(lb).peakLower && v <= features.get(lb).peakUpper)
			return features.get(lb);
		else
			return null;
	}
	
	public void reset() {
		peps = new HashSet<>();
		psms = 0;
	}
	
	public int compareTo(PeakFeature o) {
		return Double.valueOf(peakCenter).compareTo(o.peakCenter);
	}

}
