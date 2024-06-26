/*
 *    Copyright 2022 University of Michigan
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

package edu.umich.andykong.ptmshepherd.core;

import static edu.umich.andykong.ptmshepherd.core.AAMasses.protMass;
import static edu.umich.andykong.ptmshepherd.core.MXMLReader.stripChargeState;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Floats;
import edu.umich.andykong.ptmshepherd.PTMShepherd;

import java.util.*;

import umich.ms.fileio.filetypes.mzbin.MZBINFile;

public class Spectrum implements Comparable<Spectrum> {

	public int scanNum;
	public int charge;
	double precursorMass, rt;
	double monoMass, targetMass;
	public float [] peakMZ;
	public float [] peakInt;
	double norm;
	public String scanName;
	public int msLevel;
	public double[] averageFragMass;
	public double averageIonMass;
	public double basePeakInt;
	
	static final float [] fact;
	static {
		fact = new float[128];
		fact[0] = 0;
		for(int i = 1; i < 128; i++)
			fact[i] = (float)(fact[i-1] + Math.log(i));
	}

	//this constructor parses batmass compatible filetypes
	public Spectrum(int nfrags) {
		peakMZ = new float[nfrags];
		peakInt = new float[nfrags];
		norm = -1;
		basePeakInt = findBasePeakInt();
	}

	//this constructor parses MZBin files and MGF files
	public Spectrum(String scanname, int scannum, int z, int mslevel, double precursormass, double rettime, float[] peakmz, float[] peakint) {
		scanName = stripChargeState(scanname);
		scanNum = scannum;
		charge = z;
		precursorMass = precursormass;
		rt = rettime;
		peakMZ = peakmz;
		peakInt = peakint;
		msLevel = mslevel;
		norm = -1;
		basePeakInt = findBasePeakInt();
	}

	public Spectrum(MZBINFile.MZBINSpectrum s, String runName) {
		scanNum = s.scanNum;
		scanName = runName + "." + this.scanNum + "." + this.scanNum;
		charge = s.charge;
		msLevel = s.msLevel;
		rt = s.retentionTime;
		precursorMass = s.selectedMz;
		monoMass = s.selectedMz;
		targetMass = s.isolatedMz;
		peakMZ = Floats.toArray(Doubles.asList(s.peakMZ));
		peakInt = s.peakInt;
		basePeakInt = findBasePeakInt();
		norm = -1;
	}

	public MZBINFile.MZBINSpectrum toMZBINSpectrum() {
		return new MZBINFile.MZBINSpectrum(scanNum, scanName, peakMZ.length, (float) rt, (float) precursorMass, msLevel, charge, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, precursorMass, 0, 0, Doubles.toArray(Floats.asList(peakMZ)), peakInt);
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

	public void logTransform() {
		for (int i = 0; i < this.peakInt.length; i++)
			peakInt[i] = (float) Math.log10(peakInt[i]);
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

	public int getFrags(String seq, float [] mods, double ppmTol, String it, float dMass) {
		int nB = 0, nY = 0;
		int maxCharge = 1;

		char ionType = it.charAt(0);

		float [] aaMasses = AAMasses.monoisotopic_masses;
		float [] fragTypeShifts = AAMasses.ionTypeShifts;

		float cM = 0;
		double tol;
		int fP = 0;
		int cLen = seq.length();

		ArrayList<Character> nIonTypes = new ArrayList<>();
		ArrayList<Character> cIonTypes = new ArrayList<>();

		if (ionType == 'a')
			nIonTypes.add('a');
		if (ionType == 'b')
			nIonTypes.add('b');
		if (ionType == 'c')
			nIonTypes.add('c');
		if (ionType == 'x')
			cIonTypes.add('x');
		if (ionType == 'y')
			cIonTypes.add('y');
		if (ionType == 'z')
			cIonTypes.add('z');

		float nTermMass;
		for (Character iType : nIonTypes) {
			nTermMass = fragTypeShifts[iType - 'a'] + dMass;
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
					}
				}
			}
		}
		float cTermMass;
		for (Character iType : cIonTypes) {
			cTermMass = fragTypeShifts[iType - 'x' + 3] + dMass;
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
					}
				}
			}
		}
		return nB+nY;
	}

	public float[] getMatchedFrags(String seq, float [] mods, double ppmTol, String it, float dMass) {
		int maxCharge = 1;

		float [] aaMasses = AAMasses.monoisotopic_masses;
		float [] fragTypeShifts = AAMasses.ionTypeShifts;

		float cM = 0;
		double tol;
		int fP = 0;
		int cLen = seq.length();

		float [] matchedFragInts = new float[cLen * it.length()];
		Arrays.fill(matchedFragInts, -1.0F);

		ArrayList<Character> nIonTypes = new ArrayList<>();
		ArrayList<Character> cIonTypes = new ArrayList<>();

		for (char ionType : it.toCharArray()) {
			if (ionType == 'a')
				nIonTypes.add('a');
			if (ionType == 'b')
				nIonTypes.add('b');
			if (ionType == 'c')
				nIonTypes.add('c');
			if (ionType == 'x')
				cIonTypes.add('x');
			if (ionType == 'y')
				cIonTypes.add('y');
			if (ionType == 'z')
				cIonTypes.add('z');
		}

		float nTermMass;
		int cIType = 0;
		for (Character iType : nIonTypes) {
			nTermMass = fragTypeShifts[iType - 'a'] + dMass;
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
						matchedFragInts[cIType * cLen + i] = peakInt[fP];
					}
				}
			}
			cIType++;
		}
		float cTermMass;
		for (Character iType : cIonTypes) {
			cTermMass = fragTypeShifts[iType - 'x' + 3] + dMass;
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
						matchedFragInts[cIType * cLen + i] = peakInt[fP];
					}
				}
			}
			cIType++;
		}
		return matchedFragInts;
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

	public double findIonNeutral(double neutralIonMass, double ppmTol, int maxCharge) {
		double ionIntensity = 0;
		int charge = Math.min(this.charge, maxCharge);
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
		return (neutralMass + charge * protMass) / (float) charge;
	}

	/**
	 * Convert m/z [M+H+]x+ to neutral mass M
	 * @param mz m/z
	 * @param charge z
	 * @return
	 */
	public static float mzToNeutralMass(float mz, int charge) {
		return (mz - protMass) * charge;
	}

	public float[][] calcImmoniumPeaks(int min, int max, String seq, float[] mods, String filterIonTypes, int maxCharge, float dmass, float tol) { //todo I don't think these mins and maxes are very well informed
		ArrayList<Float> knownPeaks = calculatePeptideFragments(seq, mods, filterIonTypes, maxCharge, dmass);
		//ArrayList<Float> shiftedPeaks = new ArrayList<>();
		//for (Float peak : knownPeaks) //todo test?
		//	shiftedPeaks.add(peak + dmass);
		//knownPeaks.addAll(shiftedPeaks);

		this.averageIonMass = 0;
		ArrayList<Peak> ps = new ArrayList<>();
		for (int i = 0; i < peakMZ.length; i++) {
			if (peakMZ[i] > max) //todo remove min max
				break;
			if (peakMZ[i] > min) { //todo remove min max
				double absTol = peakMZ[i] * tol / 1000000;
				boolean skipFlag = false;
				for (Float peak : knownPeaks) {
					if (Math.abs(peak - peakMZ[i]) < absTol) { //todo ppm tol
						skipFlag = true;
						break;
					}
				}
				if (skipFlag == false) {
					ps.add(new Peak(peakMZ[i], peakInt[i], (float)absTol));
					this.averageIonMass += peakMZ[i];
				}
			}
		}
		float[][] peaks =  new float[ps.size()][3];
		for (int i = 0; i < ps.size(); i++) {
			peaks[i][0] = ps.get(i).MZ;
			peaks[i][1] = ps.get(i).Int;
			peaks[i][2] = ps.get(i).Tol;
		}
		return peaks;
	}

	public float[][] calcCapYPeaks(String seq, float[] mods, String filterIonTypes, int maxCharge, float pepMass, float dmass, float tol) {
		double minPeptideIon = -250.0;
		ArrayList<Float> knownPeaks = calculatePeptideFragments(seq, mods, filterIonTypes, maxCharge, 0.0f);
		knownPeaks.addAll(calculatePeptideFragments(seq, mods, filterIonTypes, maxCharge, dmass));
		ArrayList<Peak> ps = new ArrayList<>();
		int localMaxCharge = Math.min(this.charge, maxCharge);

		for (int z = 1; z <= localMaxCharge; z++) { //todo charge
		//for (int z = 1; z <= 1; z++) {
			double adjPepMass = (pepMass + 1.00727 * z) / z;
			for (int i = 0; i < peakMZ.length; i++) {
				double absTol = peakMZ[i] * tol / 1000000;
				boolean skipFlag = false;
				for (Float peak : knownPeaks) {
					if (Math.abs(peak - peakMZ[i]) < absTol) {
						skipFlag = true;
						break;
					}
				}
				//if (Math.abs(peakMZ[i] - adjPepMass < 0.01))
				//	skipFlag = true;
				if (skipFlag == false)
					if ((peakMZ[i] - adjPepMass) < minPeptideIon)
						ps.add(new Peak((float)((peakMZ[i] - adjPepMass) * z), 0, (float)absTol));
					else
						ps.add(new Peak((float)((peakMZ[i] - adjPepMass) * z), peakInt[i], (float)absTol));
			}
		}
		float[][] peaks =  new float[ps.size()][3];
		for (int i = 0; i < ps.size(); i++) {
			peaks[i][0] = ps.get(i).MZ;
			peaks[i][1] = ps.get(i).Int;
			peaks[i][2] = ps.get(i).Tol;
		}
		return peaks;
	}

	private ArrayList<Float> calculatePeptideFragments(String seq, float[] mods, String ionTypes, int maxCharge, float dmass) {
		ArrayList<Float> knownFrags = new ArrayList<>();

		ArrayList<Character> nIonTypes = new ArrayList<>();
		ArrayList<Character> cIonTypes = new ArrayList<>();
		for (int i = 0; i < ionTypes.length(); i++) {
			char curIonType = ionTypes.charAt(i);
			if (curIonType == 'a' || curIonType == 'b' || curIonType == 'c')
				nIonTypes.add(curIonType);
			else if (curIonType == 'x' || curIonType == 'y' || curIonType == 'z')
				cIonTypes.add(curIonType);
		}

		float [] aaMasses = AAMasses.monoisotopic_masses;
		float [] fragTypeShifts = AAMasses.ionTypeShifts;
		int cLen = seq.length();

		float nTermMass;
		for (Character iType : nIonTypes) {
			nTermMass = fragTypeShifts[iType - 'a'];
			for (int ccharge = 1; ccharge <= maxCharge; ccharge++) { //loop through charge states
				float cmass = AAMasses.monoisotopic_nterm_mass + nTermMass;
				for (int i = 0; i < cLen - 1; i++) { //loop through positions on the peptide
					cmass += (aaMasses[seq.charAt(i) - 'A'] + mods[i]) / ccharge;
					knownFrags.add(cmass);
					if (Math.abs(dmass) > 0.001) { //add fragments with mass shift to known fragments
						cmass += dmass / ccharge;
						knownFrags.add(cmass);
						cmass -= dmass / ccharge;
					}
				}
			}
		}
		float cTermMass;
		for (Character iType : cIonTypes) {
			cTermMass = fragTypeShifts[iType - 'x' + 3];
			for (int ccharge = 1; ccharge <= maxCharge; ccharge++) {
				float cmass = (cTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
				for (int i = 0; i < cLen - 1; i++) {
					cmass += (aaMasses[seq.charAt(cLen - 1 - i) - 'A'] + mods[cLen - 1 - i]) / ccharge;
					knownFrags.add(cmass);
					if (Math.abs(dmass) > 0.1) { //add fragments with mass shift to known fragments
						cmass += dmass / ccharge;
						knownFrags.add(cmass);
						cmass -= dmass / ccharge;
					}
				}
			}
		}

		return knownFrags;
	}

	public HashMap<Character, float[][]> calcSquigglePeaks(float ppmTol, String seq, float[] mods, String ionTypes, String filterIonTypes, int maxCharge) {
		HashMap<Character, float[][]> squigglePeaks = new HashMap<>();

		ArrayList<Character> nIonTypes = new ArrayList<>();
		ArrayList<Character> cIonTypes = new ArrayList<>();

		this.averageFragMass = new double[ionTypes.length()];
		double minSquiggleIon = -250.0;
		int iTypeIndx = 0;

		for (int i = 0; i < ionTypes.length(); i++) {
			char curIonType = ionTypes.charAt(i);
			if (curIonType == 'a' || curIonType == 'b' || curIonType == 'c')
				nIonTypes.add(curIonType);
			else if (curIonType == 'x' || curIonType == 'y' || curIonType == 'z')
				cIonTypes.add(curIonType);
		}

		ArrayList<Peak> ps;

		float [] aaMasses = AAMasses.monoisotopic_masses;
		float [] fragTypeShifts = AAMasses.ionTypeShifts;
		int cLen = seq.length();

		ArrayList<Float> knownFrags = calculatePeptideFragments(seq, mods, filterIonTypes, maxCharge, 0.0f);

		float nTermMass;
		for (Character iType : nIonTypes) {
			ps = new ArrayList<>();
			nTermMass = fragTypeShifts[iType - 'a'];
			for (int j = 0; j < peakMZ.length; j++) {
				boolean skipFlag = false;
				float trueTol = ppmTol * peakMZ[j] / (float)1000000;
				for (Float ion : knownFrags) {
					if (Math.abs(peakMZ[j] - ion) < trueTol) {
						skipFlag = true;
						break;
					}
				}
				ArrayList<Peak> cPeaksNaked = new ArrayList<>();
				for (int ccharge = 1; ccharge <= maxCharge; ccharge++) { //loop through charge states
					float cmass = AAMasses.monoisotopic_nterm_mass + nTermMass;
					for (int i = 0; i < cLen - 1; i++) { //loop through positions on the peptide
						cmass += (aaMasses[seq.charAt(i) - 'A'] + mods[i]) / ccharge;
						if (skipFlag)
							cPeaksNaked.add(new Peak(peakMZ[j] - cmass, 0, trueTol));
						else {
							if (peakMZ[j] - cmass < minSquiggleIon)
								cPeaksNaked.add(new Peak(peakMZ[j] - cmass, 0, trueTol));
							else {
								cPeaksNaked.add(new Peak(peakMZ[j] - cmass, peakInt[j], trueTol));
								this.averageFragMass[iTypeIndx] += peakMZ[j];
							}
						}
					}
				}
				for(int i = 0; i < cPeaksNaked.size(); i++)
					ps.add(cPeaksNaked.get(i));
			}
			squigglePeaks.put(iType, peaksWTolToArray(ps));
			iTypeIndx++;
		}

		float cTermMass;
		for (Character iType : cIonTypes) {
			ps = new ArrayList<>();
			cTermMass = fragTypeShifts[iType - 'x' + 3];
			for (int j = 0; j < peakMZ.length; j++) {
				boolean skipFlag = false;
				float trueTol = ppmTol * peakMZ[j] / (float) 1000000;
				for (Float ion : knownFrags) {
					if (Math.abs(peakMZ[j] - ion) < trueTol) {
						skipFlag = true;
						break;
					}
				}
				ArrayList<Peak> cPeaksNaked = new ArrayList<>();
				for (int ccharge = 1; ccharge <= maxCharge; ccharge++) {
					float cmass = (cTermMass + ccharge * AAMasses.monoisotopic_nterm_mass) / ccharge;
					for (int i = 0; i < cLen - 1; i++) {
						cmass += (aaMasses[seq.charAt(cLen - 1 - i) - 'A'] + mods[cLen - 1 - i]) / ccharge;
						if (skipFlag)
							cPeaksNaked.add(new Peak(peakMZ[j] - cmass, 0, trueTol));
						else {
							if (peakMZ[j] - cmass < minSquiggleIon)
								cPeaksNaked.add(new Peak(peakMZ[j] - cmass, 0, trueTol));
							else {
								cPeaksNaked.add(new Peak(peakMZ[j] - cmass, peakInt[j], trueTol));
								this.averageFragMass[iTypeIndx] += peakMZ[j];
							}
						}
					}
				}
				for(int i = 0; i < cPeaksNaked.size(); i++)
					ps.add(cPeaksNaked.get(i));
			}
			squigglePeaks.put(iType, peaksWTolToArray(ps));
			iTypeIndx++;
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

	float[][] peaksWTolToArray (ArrayList<Peak> peaks) {
		float[][] peakArr = new float[peaks.size()][3];
		for (int i = 0; i < peaks.size(); i ++) {
			peakArr[i][0] = peaks.get(i).MZ;
			peakArr[i][1] = peaks.get(i).Int;
			peakArr[i][2] = peaks.get(i).Tol;
		}
		return peakArr;
	}

	public double getAverageIonMass() {
		return this.averageIonMass;
	}

	public double[] getAverageFragMass() {
		return this.averageFragMass;
	}

	public int compareTo(Spectrum o) {
		return scanNum - o.scanNum;
	}
	public double getPrecursorMass() {
		return precursorMass;
	}

	public float[] getPeakMZ() {
		return this.peakMZ;
	}

	public float[] getPeakInt() {
		return this.peakInt;
	}
}

class Peak implements Comparable<Peak> {
	float MZ, Int, Tol;
	public Peak(float MZ, float Int) {
		this.MZ = MZ;
		this.Int = Int;
	}
	public Peak(float MZ, float Int, float Tol) {
		this.MZ = MZ;
		this.Int = Int;
		this.Tol = Tol;
	}
	public int compareTo(Peak p) {
		return Double.compare(Math.abs(this.MZ), Math.abs(p.MZ));
	}
}
