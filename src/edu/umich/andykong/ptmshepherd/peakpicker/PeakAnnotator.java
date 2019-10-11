package edu.umich.andykong.ptmshepherd.peakpicker;

import java.io.*;
import java.util.*;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.AAMasses;

public class PeakAnnotator {
	
	ArrayList<String> mods;
	ArrayList<Double> mod_diffs;
	ArrayList<Integer> allowed_list;
	String userMods = "";
	ArrayList<String> vModNames;
	ArrayList<Double> vModMasses;

	static final double C13delta = 1.003355;
	static final double modEqual_tol = 0.001;
	//static final double mod_tol = 0.002;
	static final double mod_tol = Double.parseDouble(PTMShepherd.getParam(("precursor_tol")));
	
	int [] indices;
	static final int maxDepth = 3;
	static int cmaxDepth;
	static boolean found, debug;
	
	public void annotateTSV(File inTSV, File outTSV) throws Exception {
        //ArrayList<String> vModNames = new ArrayList<>();
        //ArrayList<Double> vModMasses = new ArrayList<>();

		//if(mods == null)
		//	init(varMods);
		
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
			out.print("\tPotential Modification "+i);
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
			int [] cres = checkMod(masses[i]);
			for(int j = 0; j < cres.length; j++) {
				if(cres[j] != -1) {
					res[i][j] = mods.get(cres[j]);
				} else
					res[i][j] = "";
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
	
	public int [] checkMod(double v) {
//		System.out.println(v);
//		debug=true;
		int [] res = new int[maxDepth];
		Arrays.fill(indices,-1);
		Arrays.fill(res,-1);
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
	
	public void init(String varMods) throws Exception {
		mods = new ArrayList<String>(); //all modification names
		mod_diffs = new ArrayList<Double>(); //all modification mass shifts
		allowed_list = new ArrayList<Integer>(); //mods that can be added at any given step
		indices = new int[128];
		userMods = new String("");
		vModNames = new ArrayList<String>();
		vModMasses = new ArrayList<Double>();

		userMods = varMods;

		if(!varMods.equals("")) {
			List<String> tMods = Arrays.asList(varMods.split(","));
			for(int i = 0; i < tMods.size(); i++){
				List<String> md = Arrays.asList(tMods.get(i).split(":"));
				String m = md.get(0);
				double dMass = Double.parseDouble(md.get(1));
				vModNames.add(m);
				vModMasses.add(dMass);
			}
		}

		BufferedReader in = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("unimod_20191002.txt")));
		String cline;
		//add isotopic peaks to modification list
		addMod("First isotopic peak",C13delta);
		addMod("Second isotopic peak",2*C13delta);
		addMod("Third isotopic peak",3*C13delta);
		//make isotopic peaks statically accessible
		for(int i = 0; i < 3; i++){
			allowed_list.add(i);
		}
		//add user-defined mods to modification list
        for (int i = 0; i < vModNames.size(); i++) {
            addMod(vModNames.get(i), vModMasses.get(i));
            allowed_list.add(i);
        };

		while((cline = in.readLine())!= null) {
			String [] sp = cline.split("\\t");
			addMod(sp[0],Double.parseDouble(sp[1]));
		}
		in.close();

		for(int i = 0; i < AAMasses.monoisotopic_masses.length; i++) {
			if(AAMasses.monoisotopic_masses[i] > 0) {
				addMod("Addition of " + (char)('A'+i),(double)AAMasses.monoisotopic_masses[i]);
				addMod("Deletion of " + (char)('A'+i),-1.0*(double)AAMasses.monoisotopic_masses[i]);
			}
		}
	}
	
}
