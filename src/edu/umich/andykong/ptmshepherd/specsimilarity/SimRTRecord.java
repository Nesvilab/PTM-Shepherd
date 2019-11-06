package edu.umich.andykong.ptmshepherd.specsimilarity;

public class SimRTRecord {
	double mass;
	int originalOrder;
	
	int count;
	Variance sim, deltart;
	
	public SimRTRecord(double mass, int originalOrder) {
		this.mass = mass;
		this.originalOrder = originalOrder;
		count = 0;
		sim = new Variance();
		deltart = new Variance();
	}
	
	public void updateWithLine(String [] sp) {
		count++;
		if(Integer.parseInt(sp[4]) == 0) { //not zero pep
			if (Integer.parseInt(sp[6]) > 0) //num spectra > 0
				deltart.update(Double.parseDouble(sp[5]));
			if (Integer.parseInt(sp[9]) > 0) //num spectra > 0
				sim.update(Double.parseDouble(sp[7]));
		}
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		
		//General stats
		sb.append(String.format("%.4f\t%d", mass,count));
		double tmean = sim.getMean();
		double tvar = sim.getVariance();
		if (sim.getMean() == 0.0)
			tmean = Double.NaN;
		if(sim.getVariance() == 0.0)
			tvar = Double.NaN;
		sb.append(String.format("\t%.4f\t%.4f",tmean,tvar));
		tmean = deltart.getMean();
		tvar = deltart.getVariance();
		if (deltart.getMean() == 0.0)
			tmean = Double.NaN;
		if(deltart.getVariance() == 0.0)
			tvar = Double.NaN;
		sb.append(String.format("\t%.4f\t%.4f", tmean, tvar));
		
		return sb.toString();
	}
	
}
