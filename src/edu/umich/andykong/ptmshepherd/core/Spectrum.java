package edu.umich.andykong.ptmshepherd.core;

import edu.umich.andykong.ptmshepherd.PTMShepherd;

import java.lang.reflect.Array;
import java.util.*;

public class Spectrum implements Comparable<Spectrum> {

	int scanNum, charge;
	double precursorMass, rt;
	double monoMass, targetMass;
	float [] peakMZ;
	float [] peakInt;
	double norm;
	public String scanName;
	public int msLevel;
	
	static float [] fact;

	//this constructor parses batmass compatible filetypes
	public Spectrum(int nfrags) {
		peakMZ = new float[nfrags];
		peakInt = new float[nfrags];
		norm = -1;
		if(fact == null) {
			fact = new float[128];
			fact[0] = 0;
			for(int i = 1; i < 128; i++)
				fact[i] = (float)(fact[i-1] + Math.log(i));
		}
	}

	//this constructor parses MZBin files and MGF files
	public Spectrum(String scanname, int scannum, int z, int mslevel, double precursormass, double rettime, float[] peakmz, float[] peakint) {
		scanName = scanname;
		scanNum = scannum;
		charge = z;
		precursorMass = precursormass;
		rt = rettime;
		peakMZ = peakmz;
		peakInt = peakint;
		msLevel = mslevel;
		norm = -1;
		if(fact == null) {
			fact = new float[128];
			fact[0] = 0;
			for(int i = 1; i < 128; i++)
				fact[i] = (float)(fact[i-1] + Math.log(i));
		}
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(scanNum+" -");
		for(int i = 0; i < peakMZ.length; i++)
			sb.append(String.format(" [%.2f %.2f]",peakMZ[i],peakInt[i]));
		return sb.toString();
	}
	
	public double getSimilarity(Spectrum s, double ppmTol) {
		double sum = 0;
		double ppm = ppmTol*(1.0/1000000);
		int pos = 0;
		for(int i = 0; i < peakMZ.length; i++) {
			for(int j = pos; j < s.peakMZ.length; j++) {
				if(s.peakMZ[j] < peakMZ[i]*(1-ppm))
					pos = j;
				else if(s.peakMZ[j] > peakMZ[i]*(1+ppm))
					break;
				else {
					sum += peakInt[i]*s.peakInt[j];
				}
			}
		}
		return sum;
	}
	
	public double cosineSimilarity(Spectrum s, double ppmTol) {
		double n1 = norm(ppmTol);
		double n2 = s.norm(ppmTol);
		return (getSimilarity(s, ppmTol) / (n1 * n2));
	}
	
	public double averageSimilarity(Collection<Spectrum> specs, double ppmTol) {
		double sum = 0;
		for(Spectrum s : specs) {
			if (s == null) {
				continue;
			}
			sum += cosineSimilarity(s, ppmTol);
		}
		return sum / specs.size();
	}
	
	public double norm(double ppmTol) {
		if(norm >= 0)
			return norm;
		norm = Math.sqrt(getSimilarity(this, ppmTol));
		return norm;
	}
	
	public void condition(int topN, double ratio) {
		ArrayList<Peak> peaks = new ArrayList<Peak>();
		float mv = 0;
		for(int i = 0 ; i < peakInt.length; i++) {
			peaks.add(new Peak(peakMZ[i],peakInt[i]));
			if(peakInt[i] > mv)
				mv = peakInt[i];
		}
				
		Collections.sort(peaks, new Comparator<Peak>() {
			public int compare(Peak o1, Peak o2) {
				return -1*Float.valueOf(o1.Int).compareTo(o2.Int);
			}
		});
		
		while(peaks.size() > topN)
			peaks.remove(peaks.size()-1);
		while(peaks.size() > 0 && (peaks.get(peaks.size()-1).Int < peaks.get(0).Int*ratio))
			peaks.remove(peaks.size()-1);
		
		Collections.sort(peaks, new Comparator<Peak>() {
			public int compare(Peak o1, Peak o2) {
				return Float.valueOf(o1.MZ).compareTo(o2.MZ);
			}
		});
		
		peakMZ = new float[peaks.size()];
		peakInt = new float[peaks.size()];
		for(int i = 0; i < peaks.size(); i++) {
			peakMZ[i] = peaks.get(i).MZ;
			peakInt[i] = (100*peaks.get(i).Int) / mv;
		}
	}

	public void conditionOptNorm(int topN, double ratio, boolean normalize) {
		ArrayList<Peak> peaks = new ArrayList<Peak>();
		float mv = 0;
		for(int i = 0 ;i < peakInt.length; i++) {
			peaks.add(new Peak(peakMZ[i],peakInt[i]));
			if(peakInt[i] > mv)
				mv = peakInt[i];
		}

		Collections.sort(peaks, new Comparator<Peak>() {
			public int compare(Peak o1, Peak o2) {
				return -1*Float.valueOf(o1.Int).compareTo(o2.Int);
			}
		});

		while(peaks.size() > topN)
			peaks.remove(peaks.size()-1);
		while(peaks.size() > 0 && (peaks.get(peaks.size()-1).Int < peaks.get(0).Int*ratio))
			peaks.remove(peaks.size()-1);

		Collections.sort(peaks, new Comparator<Peak>() {
			public int compare(Peak o1, Peak o2) {
				return Float.valueOf(o1.MZ).compareTo(o2.MZ);
			}
		});

		peakMZ = new float[peaks.size()];
		peakInt = new float[peaks.size()];
		for(int i = 0; i < peaks.size(); i++) {
			peakMZ[i] = peaks.get(i).MZ;
			if (normalize)
				peakInt[i] = (100*peaks.get(i).Int) / mv;
			else
				peakInt[i] = peaks.get(i).Int;
		}
	}
	
	public float getHyper(String seq, float [] mods, double ppmTol) {
		float score = 0.0f;
		float iB = 0.0f, iY = 0.0f;
		int nB = 0, nY = 0;
		int maxCharge = Math.min(Integer.parseInt(PTMShepherd.getParam("spectra_maxfragcharge")), charge);

		float [] aaMasses = AAMasses.monoisotopic_masses;
		float [] fragTypeShifts = AAMasses.ionTypeShifts;

		float cM = 0;
		double tol;
		int fP = 0;
		int cLen = seq.length();

		ArrayList<Character> nIonTypes = new ArrayList<>();
		ArrayList<Character> cIonTypes = new ArrayList<>();

		if (PTMShepherd.getParam("iontype_a").trim().equals("1"))
			nIonTypes.add('a');
		if (PTMShepherd.getParam("iontype_b").trim().equals("1"))
			nIonTypes.add('b');
		if (PTMShepherd.getParam("iontype_c").trim().equals("1"))
			nIonTypes.add('c');
		if (PTMShepherd.getParam("iontype_x").trim().equals("1"))
			cIonTypes.add('x');
		if (PTMShepherd.getParam("iontype_y").trim().equals("1"))
			cIonTypes.add('y');
		if (PTMShepherd.getParam("iontype_z").trim().equals("1"))
			cIonTypes.add('z');

		float nTermMass;
		for (Character iType : nIonTypes) {
			nTermMass = fragTypeShifts[iType - 'a'];
			for (int ccharge = 1; ccharge <= maxCharge; ccharge++) {
				double cmass = AAMasses.monoisotopic_nterm_mass + nTermMass;
				fP = 0;
				for (int i = 0; i < cLen - 1; i++) {
					cmass += (aaMasses[seq.charAt(i) - 'A'] + mods[i]) / ccharge;
					tol = cmass * (ppmTol / 1000000.0);
					while (fP < peakMZ.length && peakMZ[fP] < (cmass - tol))
						fP++;
					if (fP < peakMZ.length && Math.abs(peakMZ[fP] - cmass) < tol) {
						cM = peakInt[fP];
					} else
						cM = 0;
					if (cM > 0) {
						nB++;
						iB += cM;
					}
				}
			}
		}
		float cTermMass;
		for (Character iType : cIonTypes) {
			cTermMass = fragTypeShifts[iType - 'x' + 3];
			for (int ccharge = 1; ccharge <= maxCharge; ccharge++) {
				//double cmass = (AAMasses.monoisotopic_cterm_mass + (ccharge+1)*AAMasses.monoisotopic_nterm_mass) / ccharge;
				double cmass = (cTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
				fP = 0;
				for (int i = 0; i < cLen - 1; i++) {
					cmass += (aaMasses[seq.charAt(cLen - 1 - i) - 'A'] + mods[cLen - 1 - i]) / ccharge;
					tol = cmass * (ppmTol / 1000000.0);
					while (fP < peakMZ.length && peakMZ[fP] < (cmass - tol))
						fP++;
					if (fP < peakMZ.length && Math.abs(peakMZ[fP] - cmass) < tol) {
						cM = peakInt[fP];
					} else
						cM = 0;
					if (cM > 0) {
						nY++;
						iY += cM;
					}
				}
			}
		}
		
		score = fact[nB] + fact[nY];
		if(iB > 1)
			score += Math.log(iB);
		if(iY > 1)
			score += Math.log(iY);
		
		return score;
	}
	
	public int getFrags(String seq, float [] mods, double ppmTol) {
		float iB = 0.0f, iY = 0.0f;
		int nB = 0, nY = 0;
		int maxCharge = (charge==2)?1:2;

		float [] aaMasses = AAMasses.monoisotopic_masses;
		float [] fragTypeShifts = AAMasses.ionTypeShifts;

		float cM = 0;
		double tol;
		int fP = 0;
		int cLen = seq.length();

		ArrayList<Character> nIonTypes = new ArrayList<>();
		ArrayList<Character> cIonTypes = new ArrayList<>();

		if (PTMShepherd.getParam("iontype_a").trim().equals("1"))
			nIonTypes.add('a');
		if (PTMShepherd.getParam("iontype_b").trim().equals("1"))
			nIonTypes.add('b');
		if (PTMShepherd.getParam("iontype_c").trim().equals("1"))
			nIonTypes.add('c');
		if (PTMShepherd.getParam("iontype_x").trim().equals("1"))
			cIonTypes.add('x');
		if (PTMShepherd.getParam("iontype_y").trim().equals("1"))
			cIonTypes.add('y');
		if (PTMShepherd.getParam("iontype_z").trim().equals("1"))
			cIonTypes.add('z');

		float nTermMass;
		for (Character iType : nIonTypes) {
			nTermMass = fragTypeShifts[iType - 'a'];
			for (int ccharge = 1; ccharge <= maxCharge; ccharge++) {
				//double cmass = AAMasses.monoisotopic_nterm_mass;
				double cmass = AAMasses.monoisotopic_nterm_mass + nTermMass;
				fP = 0;
				for (int i = 0; i < cLen - 1; i++) {
					cmass += (aaMasses[seq.charAt(i) - 'A'] + mods[i]) / ccharge;
					tol = cmass * (ppmTol / 1000000.0);
					while (fP < peakMZ.length && peakMZ[fP] < (cmass - tol))
						fP++;
					if (fP < peakMZ.length && Math.abs(peakMZ[fP] - cmass) < tol) {
						cM = peakInt[fP];
					} else
						cM = 0;
					if (cM > 0) {
						nB++;
						iB += cM;
					}
				}
			}
		}
		float cTermMass;
		for (Character iType : cIonTypes) {
			cTermMass = fragTypeShifts[iType - 'x' + 3];
			for (int ccharge = 1; ccharge <= maxCharge; ccharge++) {
				//double cmass = (AAMasses.monoisotopic_cterm_mass + (ccharge + 1) * AAMasses.monoisotopic_nterm_mass) / ccharge;
				double cmass = (cTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
				fP = 0;
				for (int i = 0; i < cLen - 1; i++) {
					cmass += (aaMasses[seq.charAt(cLen - 1 - i) - 'A'] + mods[cLen - 1 - i]) / ccharge;
					tol = cmass * (ppmTol / 1000000.0);
					while (fP < peakMZ.length && peakMZ[fP] < (cmass - tol))
						fP++;
					if (fP < peakMZ.length && Math.abs(peakMZ[fP] - cmass) < tol) {
						cM = peakInt[fP];
					} else
						cM = 0;
					if (cM > 0) {
						nY++;
						iY += cM;
					}
				}
			}
		}
		return nB+nY;
	}

	public double findIon(double ion, double ppmTol) {
		double ppmRange = ppmTol * (1.0 / 1000000) * ion;
		double max = ion + ppmRange;
		double ionIntensity = 0;
		for (int i = 0; i < peakMZ.length; i++) {
			if (peakMZ[i] > max)
				break;
			if (Math.abs(peakMZ[i] - ion) < ppmRange)
				ionIntensity += peakInt[i];
		}
		return ionIntensity;
	}

	public double findIonNeutral(double neutralIonMass, double ppmTol) {
		double ionIntensity = 0;
		for (int z = 1; z <= charge; z++) {
			double ion = (neutralIonMass + (1.00727 * (double) z)) / (double) z;
			double ppmRange = ppmTol * (1.0 / 1000000) * ion;
			double max = ion + ppmRange;
			for (int i = 0; i < peakMZ.length; i++) {
				if (peakMZ[i] > max)
					break;
				if (Math.abs(peakMZ[i] - ion) < ppmRange)
					ionIntensity += peakInt[i];
			}
		}
		return ionIntensity;
	}

	public double findBasePeakInt() {
		double bpInt = 0;
		for (int i = 0; i < peakInt.length; i++) {
			if (peakInt[i] > bpInt)
				bpInt = peakInt[i];
		}
		return bpInt;
	}

	/**
	 * Convert neutral mass to m/z [M+H+]x+, where x is the provided charge
	 * @param neutralMass neutral mass
	 * @param charge charge
	 * @return
	 */
	public static float neutralMassToMZ(float neutralMass, int charge) {
		return (neutralMass + AAMasses.protMass) / (float) charge;
	}

	/**
	 * Convert m/z [M+H+]x+ to neutral mass M
	 * @param mz m/z
	 * @param charge z
	 * @return
	 */
	public static float mzToNeutralMass(float mz, int charge) {
		return (mz - AAMasses.protMass) * charge;
	}

	public float[][] calcImmoniumIons(int min, int max) { //todo I don't think these mins and maxes are very well informed
		ArrayList<Peak> ps = new ArrayList<>();
		for (int i = 0; i < peakMZ.length; i++) {
			if (peakMZ[i] > max)
				break;
			if (peakMZ[i] > min)
				ps.add(new Peak(peakMZ[i], peakInt[i]));
		}
		float[][] peaks =  new float[ps.size()][2];
		for (int i = 0; i < ps.size(); i++) {
			peaks[i][0] = ps.get(i).MZ;
			peaks[i][1] = ps.get(i).Int;
		}
		return peaks;
	}

	public float[][] calcCapYIons(float pepMass) { //todo include tolerance
		ArrayList<Peak> ps = new ArrayList<>();
		float shift;
		for (int i = 0; i < peakMZ.length; i++) {
			shift = peakMZ[i] - pepMass;
			ps.add(new Peak(shift, peakInt[i]));
		}
		float[][] peaks =  new float[ps.size()][2];
		for (int i = 0; i < ps.size(); i++) {
			peaks[i][0] = ps.get(i).MZ;
			peaks[i][1] = ps.get(i).Int;
		}
		return peaks;
	}

	public HashMap<Character, float[][]> calcSquigglePeaks(float precursorMass, float ppmTol, String seq, float[] mods, String ionTypes, int maxCharge) { //todo max charge? //todo min max boundaries?
		HashMap<Character, float[][]> squigglePeaks = new HashMap<>();

		ArrayList<Character> nIonTypes = new ArrayList<>();
		ArrayList<Character> cIonTypes = new ArrayList<>();

		for (int i = 0; i < ionTypes.length(); i++) {
			char curIonType = ionTypes.charAt(i);
			if (curIonType == 'a' || curIonType == 'b' || curIonType == 'c')
				nIonTypes.add(curIonType);
			else if (curIonType == 'x' || curIonType == 'y' || curIonType == 'z')
				cIonTypes.add(curIonType);
		}

		ArrayList<Peak> ps = new ArrayList<>();

		//float iB = 0.0f, iY = 0.0f;
		//int nB = 0, nY = 0;
		//int maxCharge = (charge==2)?1:2;
		float [] aaMasses = AAMasses.monoisotopic_masses;
		float [] fragTypeShifts = AAMasses.ionTypeShifts;
        //float[][] tempPeaks;
		int cLen = seq.length();

		float nTermMass;

		for (Character iType : nIonTypes) {
		    ps = new ArrayList<>();
			nTermMass = fragTypeShifts[iType - 'a'];
			for (int ccharge = 1; ccharge <= maxCharge; ccharge++) { //loop through charge states
			    float cmass = AAMasses.monoisotopic_nterm_mass + nTermMass;
				for (int i = 0; i < cLen - 1; i++) { //loop through positions on the peptide
					cmass += (aaMasses[seq.charAt(i) - 'A'] + mods[i]) / ccharge;
					for (int j = 0; j < peakMZ.length; j++) { //loop through peaks in spectrum
						ps.add(new Peak(peakMZ[j] - cmass, peakInt[j]));
					}
                }
			}
			squigglePeaks.put(iType, peaksToArray(ps));
		}

		float cTermMass;
		for (Character iType : cIonTypes) {
            ps = new ArrayList<>();
			cTermMass = fragTypeShifts[iType - 'x' + 3];
			for (int ccharge = 1; ccharge <= maxCharge; ccharge++) {
                float cmass = (cTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
                for (int i = 0; i < cLen - 1; i++) {
                    cmass += (aaMasses[seq.charAt(cLen - 1 - i) - 'A'] + mods[cLen - 1 - i]) / ccharge;
                    for (int j = 0; j < peakMZ.length; j++) //loop through peaks in spectrum
                        ps.add(new Peak(peakMZ[j] - cmass, peakInt[j]));
                }
            }
            squigglePeaks.put(iType, peaksToArray(ps));
		}
		return squigglePeaks;
	}

	float[][] peaksToArray (ArrayList<Peak> peaks) {
	    float[][] peakArr = new float[peaks.size()][2];
	    for (int i = 0; i < peaks.size(); i ++) {
            peakArr[i][0] = peaks.get(i).MZ;
            peakArr[i][1] = peaks.get(i).Int;

        }
	    return  peakArr;
    }

	public int compareTo(Spectrum o) {
		return scanNum - o.scanNum;
	}
	public double getPrecursorMass() {
		return precursorMass;
	}

}

class Peak {
	float MZ, Int;
	public Peak(float MZ, float Int) {
		this.MZ = MZ;
		this.Int = Int;
	}
}
