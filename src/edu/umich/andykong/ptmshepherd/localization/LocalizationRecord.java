package edu.umich.andykong.ptmshepherd.localization;

import edu.umich.andykong.ptmshepherd.PTMShepherd;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class LocalizationRecord {

	double mass;
	int originalOrder;
	static int backgroundEnrich;
	static boolean globalEnrich = false;

	double [] aaScores;
	double [] aaScoresSafe;

	String mostImprovedSpectrum;
	double mostImprovedScore;
	
	int nTerm, improved, total;

	ArrayList<String> localPepSeqs = new ArrayList<>();
	static ArrayList<String> globalPepSeqs = new ArrayList<>();
	ArrayList<String> uniqueLocalPepSeqs = new ArrayList<>();
	static ArrayList<String> uniqueGlobalPepSeqs = new ArrayList<>();
	double [] localAANorm;
	static double [] globalAANorm;

	double [] AAnorms;

	public LocalizationRecord(double mass, int originalOrder) {
		this.mass = mass;
		this.originalOrder = originalOrder;
		backgroundEnrich = Integer.parseInt(PTMShepherd.getParam("localization_background"));
		
		aaScores = new double[26];
		aaScoresSafe = new double[26];
		nTerm = improved = total = 0;
	}
	
	public void updateWithLine(String [] sp) {

		String pep = sp[4];
		String pepSeq = (sp[1]);
		//System.out.println(pepSeq);
		
		double origScore = Double.parseDouble(sp[5]);
		double bestScore = Double.parseDouble(sp[6]);
		
		int origFrag = Integer.parseInt(sp[7]);
		int bestFrag = Integer.parseInt(sp[8]);

		//check for local or global enrichment calculation, default is global (4)
		if(backgroundEnrich > 2)
			globalEnrich = true;

		//todo implement nonlocalized psm background

		//if "localized"
		if(bestFrag > origFrag) {
			if(globalEnrich)
				globalPepSeqs.add(pepSeq);
			else
				localPepSeqs.add(pepSeq);

			//calc score diff
			double improv = bestScore-origScore;
			if(improv > mostImprovedScore) {
				mostImprovedScore = improv;
				mostImprovedSpectrum = sp[0];
			}
			improved++;
			int cnt = 0;
			int cntN = 0;
			boolean noUpper = false;
			//iterate left to right through sequences
			for(int i = 0; i < pep.length(); i++) {
				if(Character.isUpperCase(pep.charAt(i))) {
					//if localized
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
		
		/* General stats */
		sb.append(String.format("%.4f\t%d\t%d\t%.2f", mass , improved, total,
				((double)100*nTerm)/((total==0)?1:total)));

		//Top N
		int outputTopN = 3;
		
		//Normalize
		if(globalEnrich) { //if norm at global level
			if (globalAANorm == null) { //if global norms uncomputed
				if (backgroundEnrich == 3) { //if norm at peptide level
					HashSet<String> t = new HashSet<>(globalPepSeqs);
					uniqueGlobalPepSeqs = new ArrayList<>();
					uniqueGlobalPepSeqs.addAll(t);
					//calc norm from unique global peptide seqs
					globalAANorm = calcNorms(uniqueGlobalPepSeqs);
				} else { //if norm at psm level
					//calc norms from global peptide seuqneces
					globalAANorm = calcNorms(globalPepSeqs);
				}
			}
			AAnorms = globalAANorm;
		} else { //if norm at local level
			if (localAANorm == null) { //if local norms uncomputed
				if (backgroundEnrich == 1) { //if norm at peptide level
					HashSet<String> t = new HashSet<>(localPepSeqs);
					uniqueLocalPepSeqs = new ArrayList<>();
					uniqueLocalPepSeqs.addAll(t);
					//calc norm from uniqe local peptide sequences
					localAANorm = calcNorms(uniqueLocalPepSeqs);
				} else { //if norm at psm level
					//calc norm from local peptide sequences
					localAANorm = calcNorms(localPepSeqs);
				}
			}
			AAnorms = localAANorm;
		}

		double sum = 0;
		for(int i = 0; i < 26; i++) {
			sum += aaScores[i]; // == number of successfully localized peptides
			aaScoresSafe[i] = aaScores[i];
		}

		if(sum != 0) //if at least 1 peptide was successfully localized
			for(int i = 0; i < 26; i++)
				aaScores[i] /= (AAnorms[i] * sum);
		
		double [] taScores = new double[aaScores.length];
		for(int j = 0; j < aaScores.length; j++)
			taScores[j] = aaScores[j];
		double bscore = 0;
		for(int j = 0; j < outputTopN; j++) {
			int best = 25;
			for (int i = 0; i < 26; i++) {
				if ((aaScores[i] > aaScores[best]) &&		/* 	If better */
						(aaScoresSafe[i] / total > 0.01) &&	/* and bin > 5% localized */
						(aaScores[i] > 1.5) && 				/* and enrichment score > 1.5 */
						(aaScoresSafe[i] > 5)) 				/* and weighted PSMs > 5 */
					best = i;
			}
			char b = (char) ('A' + best);
			bscore = Math.max(bscore, aaScores[best]);

			if (aaScores[best] == 0.0 || (double) improved / (double) total < 0.05)
				sb.append(String.format("\t\t\t"));//%c - NaN - NaN",b));//,b,(sum == 0)?0:(aaScores[best]),(sum == 0)?0:(aaScoresSafe[best])));
			else
				sb.append(String.format("\t%c\t%.1f\t%.0f", b, (sum == 0) ? 0 : (aaScores[best]), (sum == 0) ? 0 : (aaScoresSafe[best])));
			aaScores[best] = 0;
		}
		
		//Output rest of localized scores
		for(int i = 0; i < 26; i++)
			if(LocalizationProfile.AAcnts[i] != 0)
				sb.append(String.format("\t%.4f - %.4f", taScores[i], aaScoresSafe[i]));
		
		return sb.toString();
	}

	//This function takes in a list of relevant background peptides, then calculates its background amino acid content
	public static double[] calcNorms(ArrayList<String> peps) {
		int [] AAcnts = new int[26];

		for(int i = 0; i < peps.size(); i++){
			String pepSeq = peps.get(i);
			//System.out.println(i + pepSeq);
			for(int j = 0; j < pepSeq.length(); j++){
				//System.out.println("************");
				//System.out.println(pepSeq + pepSeq.charAt(j));
				//System.out.println(pepSeq.charAt(j) - 'A');
				//System.out.println(AAcnts[pepSeq.charAt(j)-'A']);
				AAcnts[pepSeq.charAt(j)-'A'] += 1;
				//System.out.println(AAcnts[pepSeq.charAt(j)-'A']);
			}
		}

		double[] AAnorm = new double[26];
		double AAsum = 0;

		for (int i = 0; i < 26; i++) {
			AAsum += AAcnts[i];
		}
		for (int i = 0; i < 26; i++) {
			AAnorm[i] = AAcnts[i];
			//System.out.println(AAnorm[i]);
			AAnorm[i] /= AAsum;
			//System.out.println(AAnorm[i]);
			if (AAnorm[i] == 0)
				AAnorm[i] = 1;
		}

		//System.out.println("**********");
		//for (int i = 0; i < 26; i++)
		//	System.out.printf("%.5f\n", AAnorm[i]);

		return AAnorm;
		//for(int i = 0; i < 26; i++)
			//System.out.println(AAnorm[i]);
		//for (int i = 0; i < 26; i++)
			//System.out.printf("%.5f\n", AAnorm[i]);
	}


}

