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

package edu.umich.andykong.ptmshepherd.peakpicker;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.FastLocator;


public class PeakAnnotator {

	public static ArrayList<String> mods;
	public static ArrayList<Double> mod_diffs;
	ArrayList<Integer> allowed_list;
	String userMods = "";
	String modSource;
	ArrayList<String> vModNames; //privileged mods that can be matched with any mass in unimod
	ArrayList<Double> vModMasses;

	static final double C13delta = 1.00235;
	static final double modEqual_tol = 0.001;
	//static final double mod_tol = 0.01;
	static final double mod_tol = Double.parseDouble(PTMShepherd.getParam(("annotation_tol")));
	
	int [] indices;
	static final int maxDepth = 2;
	static int cmaxDepth;
	static boolean found, debug;

	/* Loaded file variables */
	public String [] headers;
	public double[][] peaks;
	public String[][] modMappings;
	public FastLocator fastLocator;
	
	public void annotateTSV(File inTSV, File outTSV, String mos, String isos, Double mod_tol) throws Exception {
        //ArrayList<String> vModNames = new ArrayList<>();
        //ArrayList<Double> vModMasses = new ArrayList<>();

		//if(mods == null)
		//	init(varMods);

		//buildOffsetList(mos, isos);

		ArrayList<String> inFile = new ArrayList<>();
		String cline;
		BufferedReader in = new BufferedReader(new FileReader(inTSV));
		while((cline = in.readLine())!= null) {
			inFile.add(cline);
		}
		in.close();
		double [] masses = new double[inFile.size()-1];
		for(int i = 1; i < inFile.size(); i++) {
			String [] sp = inFile.get(i).split("\t");
			masses[i-1] = Double.parseDouble(sp[0]);
		}
		String [][] annotations = annotate(masses);
		
		PrintWriter out = new PrintWriter(new FileWriter(outTSV));
		out.print(inFile.get(0));
		for(int i = 1; i <= maxDepth; i++)
			out.print("\tmapped_mass_"+i);
		out.println();
		for(int i = 0; i < masses.length; i++) {
			out.print(inFile.get(i+1));
			for(int j = 0; j < maxDepth; j++)
				out.print("\t"+annotations[i][j]);
			out.println();
		}
		out.close();
	}
	
	public String [][] annotate(double [] masses) throws Exception {
		String [][] res = new String[masses.length][maxDepth];
		//init(userMods);
		for(int i = 0; i < masses.length; i++) {
			if(Math.abs(masses[i]) < mod_tol){
				for(int j = 0; j < maxDepth; j++){
					res[i][j] = "";
				}
			} else {
				int [] cres = checkMod(masses[i]);
				for(int j = 0; j < cres.length; j++) {
					if (cres[j] != -1) {
						res[i][j] = mods.get(cres[j]);
					} else
						res[i][j] = "";
				}
			}
		}
		return res;
	}
	
	public void dfs(int depth, double rem) {
		if(Math.abs(rem) < mod_tol) {
			found = true;
			return;
		}
//		if(debug) {
//			System.out.print("depth: " + depth + " rem: " + rem);
//			for(int i = 0; i < depth; i++) 
//				System.out.printf(" [%.4f]",mod_diffs.get(indices[i]));
//			System.out.println();
//		}
		if(depth == cmaxDepth)
			return;
		for(int i = 0; i < allowed_list.size(); i++) {
			indices[depth] = allowed_list.get(i);
			dfs(depth+1,rem - mod_diffs.get(allowed_list.get(i)));
			if(found)
				return;
			indices[depth] = -1;
		}
	}

	public static void buildOffsetList(String massOffsets, String isotopes) {
		boolean offsetMode = false;
		double[] mos = new double[0];

		if (!massOffsets.equals("") & (!massOffsets.equals("None"))) {
			offsetMode = true;
		}

		if (offsetMode) {
			String[] tmpIsos; //temp str to hold isotopes
			String[] tmpMos = massOffsets.split("/"); //temp string to hold mass offsets
			if (isotopes.equals("")) {
				tmpIsos = new String[]{"0"};
			} else {
				tmpIsos = isotopes.split("/");
			}
			mos = new double[tmpMos.length * tmpIsos.length]; //list of mass offsets
			//mos[0] = 0.0; //include unmodified peps
			for (int i = 0; i < tmpMos.length; i++) {
				for (int j = 0; j < tmpIsos.length; j++) {
					double mo = Double.parseDouble(tmpMos[i]) + Integer.parseInt(tmpIsos[j]) * C13delta;
					mos[i * tmpIsos.length + j] = mo;
				}
			}

			Arrays.sort(mos);

			for(int i = 0; i < mos.length; i++){
				System.out.println(mos[i]);
			}

			for (int i = 0; i < mos.length; i++) {
				double massToCheck = mos[i];
				for (int j = 1; j < mos.length; j++) {
					for (int k = j; k < mos.length; k++) {
						if (Math.abs(massToCheck - (mos[j] + mos[k])) < modEqual_tol){
							//System.out.println(massToCheck + "\t" + mos[j] + "\t" + mos[k] + "***");
						}
					}
				}
			}
		}
	}

	public int [] checkMod(double v) {
//		System.out.println(v);
//		debug=true;

		int [] res = new int[maxDepth];
		Arrays.fill(indices,-1);
		Arrays.fill(res,-1);

		//first check direct match aainst unimod
		for(int i = 0; i < mods.size(); i++)
			if(Math.abs(v-mod_diffs.get(i)) < mod_tol) {
				for(int j = 0; j < allowed_list.size(); j++)
					if(allowed_list.get(j) == i) {
						res[0] = i;
						return res;
					}
				allowed_list.add(i);
				res[0] = i;
				return res;
			}
		for(int i = 0; i < vModNames.size(); i++){
			double vPriv = v - vModMasses.get(i); //privileged mass shift
			for(int j = 0; j < mods.size(); j++)
				if(Math.abs(vPriv-mod_diffs.get(j)) < mod_tol) {
					for(int k = 0; k < allowed_list.size(); k++)
						if(allowed_list.get(k) == j) {
							res[0] = i;
							res[1] = j;
							return res;
						}
					allowed_list.add(j);
					res[0] = i;
					res[1] = j;
					return res;
				}
		}
		if(debug) 
			for(int i = 0; i < allowed_list.size(); i++) {
				System.out.printf("%s - %.4f\n",mods.get(allowed_list.get(i)),mod_diffs.get(allowed_list.get(i)));
			}
		for(int i = 1; i <= maxDepth; i++) {
			cmaxDepth = i;
			found = false;
			dfs(0,v);
			if(found) {
				for(int j = 0; j < i; j++)
					res[j] = indices[j]; 
				break;
			}
		}
		if(!found) {
			mods.add(String.format("Unannotated mass-shift %.4f",v));
			mod_diffs.add(v);
			allowed_list.add(mods.size()-1);
			res[0] = mods.size()-1;
		}
		debug = false;
		return res;
	}
	
	public void addMod(String tag, double v) {
		for(int i = 0; i < mod_diffs.size(); i++) {
			if(Math.abs(v-mod_diffs.get(i)) < modEqual_tol) {
				mods.set(i, mods.get(i)+ "/"+tag);
				return;
			}
		}
		mods.add(tag);
		mod_diffs.add(v);
	}
	
	public void init(String varMods, String modSourcePath) throws Exception {
		mods = new ArrayList<String>(); //all modification names
		mod_diffs = new ArrayList<Double>(); //all modification mass shifts
		allowed_list = new ArrayList<Integer>(); //mods that can be added at any given step
		indices = new int[128];
		//userMods = new String("");
		vModNames = new ArrayList<String>();
		vModMasses = new ArrayList<Double>();

		userMods = varMods;

		if(!varMods.equals("None:0") && !varMods.equals("")) {
			List<String> tMods = Arrays.asList(varMods.split(","));
			for(int i = 0; i < tMods.size(); i++){
				List<String> md = Arrays.asList(tMods.get(i).trim().split(":"));
				String m = md.get(0).trim();
				double dMass = Double.parseDouble(md.get(1));
				vModNames.add(m);
				vModMasses.add(dMass);
			}
		}
		for (int i = 0; i < vModNames.size(); i++) {
			addMod(vModNames.get(i), vModMasses.get(i));
		}

		BufferedReader in;
		if (modSourcePath.equals("") || modSourcePath.toLowerCase().trim().equals("unimod")) {
			modSource = "unimod_20210623.txt";
			in = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(modSource)));
		} else if (modSourcePath.toLowerCase().trim().equals("common")) {
			modSource = "common_mods_20200813.txt";
			in = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(modSource)));
		} else if (modSourcePath.toLowerCase().trim().equals("glyco")) {
			modSource = "glyco_mods_20210127.txt";
			in = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(modSource)));
		} else {
			modSource = modSourcePath.trim();
			in = new BufferedReader(new FileReader(modSource));
		}

		String cline;
		//add isotopic peaks to modification list
		addMod("Isotopic peak error", -1*C13delta);
		addMod("First isotopic peak",C13delta);
		addMod("Second isotopic peak",2*C13delta);
		addMod("Third isotopic peak",3*C13delta);
		//addMod("Oxidation or Hydroxylation",15.994915);
		//make isotopic peaks statically accessible
		//add user-defined mods to modification list
        for (int i = 0; i < mods.size(); i++) {
			allowed_list.add(i);
		}

		while((cline = in.readLine())!= null) {
			String [] sp = cline.split("\\t | %");
			addMod(sp[0].trim(),Double.parseDouble(sp[1].trim()));
		}

		in.close();

		for(int i = 0; i < AAMasses.monoisotopic_masses.length; i++) {
			if(AAMasses.monoisotopic_masses[i] > 0) {
				addMod("Addition of " + (char)('A'+i),(double)AAMasses.monoisotopic_masses[i]);
				addMod("Deletion of " + (char)('A'+i),-1.0*(double)AAMasses.monoisotopic_masses[i]);
			}
		}
	}

	public void loadAnnotatedFile(File fin, double precursorTol, int precursorUnits) throws IOException {
		/* load annotated file */
		BufferedReader in = new BufferedReader(new FileReader(fin));
		this.headers = in.readLine().split("\t", -1);
		ArrayList<String> inLines = new ArrayList<>();
		String cline;
		while((cline = in.readLine())!= null)
			inLines.add(cline);

		/* populate class parameters */
		this.peaks = new double[3][inLines.size()];
		this.modMappings = new String[this.maxDepth][inLines.size()];

		for(int i = 0; i < inLines.size(); i++) {
			String [] sp = inLines.get(i).split("\t", -1);
			for (int j = 0; j < 3; j++)
				this.peaks[j][i] = Double.parseDouble(sp[j]);
			for (int j = 0; j < this.maxDepth; j++)
				this.modMappings[j][i] = sp[getColumn("mapped_mass_" + (j+1))];
		}

		/* set up indexer for fast access */
		this.fastLocator = new FastLocator(this.peaks, precursorTol, precursorUnits);

		in.close();
	}

	private int getColumn (String head) {
		for (int i = 0; i < headers.length; i++)
			if (headers[i].equals(head))
				return i;
		return -1;
	}

	public String[] getDeltaMassMappings(ArrayList<Float> massdiffs, ArrayList<Float> precursors, double precursorTol, int precursorUnits) {
		String[] massDiffAnnotations = new String[massdiffs.size()];

		/* make temp mass diff instance for sorting */
		ArrayList<MassDiffAnnotation> mdasUnimod = new ArrayList<>();
		for(int i = 0; i < this.mod_diffs.size(); i++)
			mdasUnimod.add(new MassDiffAnnotation(this.mods.get(i), this.mod_diffs.get(i)));
		Collections.sort(mdasUnimod);

		/* skip these entries */
		int zeroBin = this.fastLocator.getIndex(0.0000);
		/* get mod annotations for each delta mass */
		for (int i = 0; i < massdiffs.size(); i++) {
			if (this.fastLocator.getIndex(massdiffs.get(i)) == zeroBin) {
				massDiffAnnotations[i] = "";
				continue;
			}
			/* hold indices of possible PTMS and Unimod mods */
			ArrayList<Integer> cPtmsIs = new ArrayList<>();
			ArrayList<Integer> cUnimodIs = new ArrayList<>();
			/* set up search parameters */
			double tol = 0;
			if (precursorUnits == 0)
				tol = precursorTol;
			else if (precursorUnits == 1)
				tol = precursors.get(i) * precursorTol / 1000000;
			double lowerBound = massdiffs.get(i) - tol;
			double upperBound = massdiffs.get(i) + tol;
			/* zoom through Unimod annotations */
			for (int j = 0; j < mdasUnimod.size(); j++) {
				double cdiff = mdasUnimod.get(j).diff;
				if (cdiff > upperBound)
					break;
				else if (cdiff < lowerBound)
					continue;
				else {
					cUnimodIs.add(j);
				}
			}
			/* zoom through PTMS annotations */
			for (int j = 0; j < this.peaks[0].length; j++) {
				double cdiff = this.peaks[0][j];
				if (cdiff <= upperBound && cdiff >= lowerBound)
					cPtmsIs.add(j);
			}
			/* see if there are equivalent mods in PTMS and Unimod annotations and format output */
			DeltaMassAnnotation dma = new DeltaMassAnnotation();
			for (int j = 0; j < cUnimodIs.size(); j++) {
				dma.UnimodMasses.add(mdasUnimod.get(cUnimodIs.get(j)).diff);
				dma.UnimodAnnos.add("Mod1: " + mdasUnimod.get(cUnimodIs.get(j)).annotation);
			}
			for (int j = 0; j < cPtmsIs.size(); j++) {
				dma.PtmsMasses.add(this.peaks[0][cPtmsIs.get(j)]);
				String[] transposedAnnos = new String[this.maxDepth];
				for (int k = 0; k < this.maxDepth; k++)
					transposedAnnos[k] = this.modMappings[k][cPtmsIs.get(j)];
				dma.PtmsAnnos.add(formatMultiplePtmsMods(transposedAnnos));
			}

			massDiffAnnotations[i] = dma.toString();
		}

		return massDiffAnnotations;
	}

	public String[] getPeakApexMappings() {
		String[] singleColumnAnnotations = new String[this.modMappings[0].length];

		String[][] tStringAnnos = new String[this.modMappings[0].length][this.modMappings.length];
		for (int i  = 0; i < tStringAnnos.length; i++) {
			for (int j = 0; j < tStringAnnos[i].length; j++)
				tStringAnnos[i][j] = this.modMappings[j][i];
		}

		for (int i  = 0; i < tStringAnnos.length; i++) {
			StringBuffer sb = new StringBuffer();
			for (int j = 0; j < tStringAnnos[i].length; j++) {
				if (!tStringAnnos[i][j].equals(""))
					sb.append(tStringAnnos[i][j] + " + ");
			}
			if (sb.length() > 0)
				sb.setLength(sb.length() - 3);
			singleColumnAnnotations[i] = sb.toString();
		}

		return singleColumnAnnotations;
	}

	class MassDiffAnnotation implements Comparable<MassDiffAnnotation> {
		String annotation;
		double diff;

		MassDiffAnnotation(String anno, double dm) {
			this.annotation = anno;
			this.diff = dm;
		}

		public int compareTo(MassDiffAnnotation o) {
			return Double.valueOf(this.diff).compareTo(o.diff);
		}
	}

	class DeltaMassAnnotation {
		public ArrayList<String> UnimodAnnos;
		public ArrayList<String> PtmsAnnos;
		public ArrayList<Double> UnimodMasses;
		public ArrayList<Double> PtmsMasses;

		DeltaMassAnnotation() {
			this.UnimodAnnos = new ArrayList<>();
			this.PtmsAnnos =  new ArrayList<>();
			this.UnimodMasses = new ArrayList<>();
			this.PtmsMasses = new ArrayList<>();
		}

		@Override
		public String toString() {
			/* remove unannotated mass shfits from Unimod lists that probably shouldn't be there in the first place */
			for (int i = UnimodAnnos.size()-1; i >= 0; i--) {
				if (UnimodAnnos.get(i).contains("Unannotated")) {
					UnimodAnnos.remove(i);
					UnimodMasses.remove(i);
				}
			}

			StringBuffer sb = new StringBuffer();
			for (int i = 0; i < PtmsAnnos.size(); i++) {
				int umIndx = -1;
				for (int j = 0; j < UnimodAnnos.size(); j++) {
					if (PtmsAnnos.get(i).equals(UnimodAnnos.get(j))) {
						umIndx = j;
						break;
					}
				}
				if (umIndx > -1) {
					if (sb.length() > 0)
						sb.append(String.format("; %s (PeakApex: %.04f, Theoretical: %.04f)",
								PtmsAnnos.get(i), PtmsMasses.get(i), UnimodMasses.get(umIndx)));
					else
						sb.append(String.format("%s (PeakApex: %.04f, Theoretical: %.04f)",
								PtmsAnnos.get(i), PtmsMasses.get(i), UnimodMasses.get(umIndx)));
					UnimodAnnos.remove(umIndx);
					UnimodMasses.remove(umIndx);
				} else {
					if (sb.length() > 0)
						sb.append(String.format("; %s (PeakApex: %.04f)",
								PtmsAnnos.get(i), PtmsMasses.get(i)));
					else
						sb.append(String.format("%s (PeakApex: %.04f)",
								PtmsAnnos.get(i), PtmsMasses.get(i)));
				}
			}
			for (int i  = 0; i < UnimodAnnos.size(); i++) {
				if (sb.length() > 0)
					sb.append(String.format("; %s (Theoretical: %.04f)",
							UnimodAnnos.get(i), UnimodMasses.get(i)));
				else
					sb.append(String.format("%s (Theoretical: %.04f)",
							UnimodAnnos.get(i), UnimodMasses.get(i)));
			}
			return sb.toString();
		}
	}

	private String formatMultiplePtmsMods (String[] massAnnos) {
		StringBuffer cMods = new StringBuffer();
		for (int i = 0; i < massAnnos.length; i++) {
			if (massAnnos[i].equals(""))
				continue;
			if (cMods.length() > 0)
				cMods.append(String.format(", Mod2: %s", massAnnos[i]));
			else
				cMods.append(String.format("Mod1: %s", massAnnos[i]));
		}
		return cMods.toString();
	}
	
}
