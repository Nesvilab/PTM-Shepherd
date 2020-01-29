package edu.umich.andykong.ptmshepherd.core;

import java.util.*;

import edu.umich.andykong.ptmshepherd.localization.LocalizationProfile;

public class FastLocator {

	double [] masses;
	int [] massOrder;
	int [] hintIndex;
	
	static final int scale = 10000;
	int offset;
	
	double peakTol;
	
	public FastLocator(double [] masses, double peakTol) {
		this.masses = masses;
		this.peakTol = peakTol;
		massOrder = new int[masses.length];
		
		ArrayList<Integer> alMassOrder = new ArrayList<>();
		for(int i = 0; i < masses.length; i++) {
			massOrder[i] = i;
			alMassOrder.add(i);
		}
		
		Collections.sort(alMassOrder,new Comparator<Integer>() {
			public int compare(Integer o1, Integer o2) {
				return Double.valueOf(masses[o1]).compareTo(masses[o2]);
			}
		});
		for(int i = 0; i < masses.length; i++)
			massOrder[i] = alMassOrder.get(i);
		
		hintIndex = new int[scale];
		offset = scale / 2;
		
		Arrays.fill(hintIndex, -1);
		int curPos = 0;
		for(int i = 0; i < masses.length; i++) {
			int tpos = (int)Math.floor(masses[massOrder[i]]) + offset;
			for(int j = curPos; j <= tpos && tpos < hintIndex.length; j++)
				if(hintIndex[j] == -1)
					hintIndex[j] = i;
			curPos = tpos;
		}
		for(int i = curPos; i < hintIndex.length; i++)
			if(hintIndex[i] == -1)
				hintIndex[i] = masses.length - 1;
	}
	
	public int getIndex(double mass) {
		int fV = (int)Math.floor(mass) + offset;
		int best = -1;
		double bV = peakTol, cV;
		int lb = hintIndex[Math.max(0, Math.min(hintIndex.length-1, fV-1))];
		int ub = hintIndex[Math.max(0, Math.min(hintIndex.length-1, fV+1))];
		for(int i = lb; i <= ub; i++) {
			cV = Math.abs(masses[massOrder[i]]-mass);
			if(cV < bV) {
				best = massOrder[i];
				bV = cV; 
			}
		}
		return best;
	}
	
//	public static void main(String [] args) throws Exception {
//		int scale = 10;
//		double [] randNums = new double[100];
//		for(int i = 0; i < randNums.length; i++) {
//			randNums[i] = Math.random()*scale*2 - scale;
//			System.out.println(randNums[i]);
//		}
//		System.out.println();
//		LocalizationProfile np = new LocalizationProfile(randNums, 0.1);
//		for(int i = 0; i < np.massOrder.length; i++)
//			System.out.println(i + " " + np.masses[np.massOrder[i]]);
//
//		System.out.println();
//		
//		for(int i = 0; i < np.hintIndex.length; i++) {
//			System.out.printf("%d %d %f\n",i - np.offset, np.hintIndex[i], np.masses[np.massOrder[np.hintIndex[i]]]);
//		}
//		
//		int samps = 1000;
//		for(int i = 0; i < 1000; i++) {
//			double v = ((i + Math.random())*2*scale)/samps - scale;
//			System.out.printf("%.6f %d %.6f\n", v,np.getIndex(v),(np.getIndex(v)==-1?0:np.masses[np.getIndex(v)]));
//		}
//		
//	}
	
}
