package edu.umich.andykong.ptmshepherd.localization;

public class LocalizationRecord {

	double mass;
	int originalOrder;
	
	double [] aaScores;
	
	String mostImprovedSpectrum;
	double mostImprovedScore;
	
	int nTerm, improved, total;
	
	public LocalizationRecord(double mass, int originalOrder) {
		this.mass = mass;
		this.originalOrder = originalOrder;
		
		aaScores = new double[26];
		nTerm = improved = total = 0;
	}
	
	public void updateWithLine(String [] sp) {
		String pep = sp[4];
		
		double origScore = Double.parseDouble(sp[5]);
		double bestScore = Double.parseDouble(sp[6]);
		
		int origFrag = Integer.parseInt(sp[7]);
		int bestFrag = Integer.parseInt(sp[8]);

		if(bestFrag > origFrag) {
			double improv = bestScore-origScore;
			if(improv > mostImprovedScore) {
				mostImprovedScore = improv;
				mostImprovedSpectrum = sp[0];
			}
			improved++;
			int cnt = 0;
			int cntN = 0;
			boolean noUpper = false;
			for(int i = 0; i < pep.length(); i++) {
				if(Character.isUpperCase(pep.charAt(i))) {
					cnt++;
					if(!noUpper)
						cntN++;
				} else
					noUpper = true;
			}
			for(int i = 0; i < pep.length(); i++)
				if(Character.isUpperCase(pep.charAt(i))) {
					aaScores[pep.charAt(i)-'A'] += 1.0/cnt;
				}
			if(cnt == cntN)
				nTerm++;
		}
		total++;
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		
		//General stats
		sb.append(String.format("%.4f\t%d\t%d\t%.2f", mass,improved,total,((double)100*nTerm)/((total==0)?1:total)));
		
		//Top N
		int outputTopN = 3;
		
		//Normalize
		double sum = 0;
		for(int i = 0; i < 26; i++)
			sum += aaScores[i];
		if(sum != 0)
			for(int i = 0; i < 26; i++)
				aaScores[i] /= (LocalizationProfile.AAnorm[i] * sum);
		
		double [] taScores = new double[aaScores.length];
		for(int j = 0; j < aaScores.length; j++)
			taScores[j] = aaScores[j];
		double bscore = 0;
		for(int j = 0; j < outputTopN; j++) {
			int best = 0;
			for(int i = 0; i < 26; i++)
				if(aaScores[i] > aaScores[best])
					best = i;
			char b = (char)('A'+best);
			bscore = Math.max(bscore,aaScores[best]);
			if(Double.isNaN(aaScores[best]))
				sb.append("\tERROR");
			else
				sb.append(String.format("\t%c - %.2f",b,(sum == 0)?0:(aaScores[best])));
			aaScores[best] = 0;
		}
		
		//Output rest of localized scores
		for(int i = 0; i < 26; i++)
			if(LocalizationProfile.AAcnts[i] != 0)
				sb.append(String.format("\t%.4f", taScores[i]));
		
		return sb.toString();
	}
	
}
