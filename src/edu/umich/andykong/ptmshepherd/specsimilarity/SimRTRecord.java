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
		if(Integer.parseInt(sp[6]) > 0)
			deltart.update(Double.parseDouble(sp[5]));
		if(Integer.parseInt(sp[9]) > 0)
			sim.update(Double.parseDouble(sp[7]));
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		
		//General stats
		sb.append(String.format("%.4f\t%d", mass,count));
		sb.append(String.format("\t%.4f\t%.4f", sim.getMean(),sim.getVariance()));
		sb.append(String.format("\t%.4f", deltart.getMean()/deltart.getVariance()));
		
		return sb.toString();
	}
	
}
