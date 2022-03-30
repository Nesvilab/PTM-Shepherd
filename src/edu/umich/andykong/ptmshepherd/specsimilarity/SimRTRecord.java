package edu.umich.andykong.ptmshepherd.specsimilarity;

public class SimRTRecord {
	double mass;
	int originalOrder;
	
	int count;
	Variance sim, deltart, intensity;
	boolean calcIntensity;
	
	public SimRTRecord(double mass, int originalOrder, boolean calcIntensity) {
		this.mass = mass;
		this.originalOrder = originalOrder;
		count = 0;
		sim = new Variance();
		deltart = new Variance();
		intensity = new Variance();
		this.calcIntensity = calcIntensity;
	}
	
	public void updateWithLine(String [] sp) {
		count++;
		if (Integer.parseInt(sp[6]) > 0) //num zero bin spectra > 0
			deltart.update(Double.parseDouble(sp[5]));
		if (Integer.parseInt(sp[9]) > 0) //num zero bin spectra > 0
			sim.update(Double.parseDouble(sp[7]));
		if (calcIntensity) {
			if (Integer.parseInt(sp[11]) > 0) // num zero bin spectra > 0
				intensity.update(Double.parseDouble(sp[10]));
		}
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		
		//General stats
		sb.append(String.format("%.4f\t%d", mass, count));
		double tmean = sim.getMean();
		double tvar = sim.getVariance();
		if (sim.getMean() == 0.0)
			tmean = Double.NaN;
		if(sim.getVariance() == 0.0)
			tvar = Double.NaN;
		sb.append(String.format("\t%.3f\t%.3f", tmean, tvar));
		tmean = deltart.getMean();
		tvar = deltart.getVariance();
		if (deltart.getMean() == 0.0)
			tmean = Double.NaN;
		if(deltart.getVariance() == 0.0)
			tvar = Double.NaN;
		sb.append(String.format("\t%.1f\t%.0f", tmean, tvar));

		if (calcIntensity) {
			intensity.logTransform();
			double tmed = intensity.getMedian();
			tvar = intensity.getVariance();
			sb.append(String.format("\t%2.4e\t%2.4e", tmed, tvar));
		}

		return sb.toString();
	}
	
}
