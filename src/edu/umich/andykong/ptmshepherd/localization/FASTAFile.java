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

package edu.umich.andykong.ptmshepherd.localization;

import java.io.*;
import java.util.*;

public class FASTAFile {
	
	public ArrayList<Protein> prots;
	public double [] norms;
	public int [] counts;
	
	public void countAA() {
		norms = new double[26];
		counts = new int[26];
		
		for(int i = 0; i < prots.size(); i++) {
			
		}
		
		long sum = 0;
		for(int i = 0; i < counts.length; i++)
			sum += counts[i];
		
	}
		
    public static String reverse(String s) {
        StringBuffer rs = new StringBuffer();
        for(int i = 0; i < s.length(); i++)
                rs.append(s.charAt(s.length()-i-1));
        return rs.toString();
    }
    
    public FASTAFile() {
    	prots = new ArrayList<Protein>();
    }
    
    public void placeFirst(String prefix) {
    	ArrayList<Protein> selected = new ArrayList<Protein>();
    	ArrayList<Protein> notselected = new ArrayList<Protein>();
    	for(int i = 0; i < prots.size(); i++) {
    		if(prots.get(i).seqname.startsWith(prefix))
    			selected.add(prots.get(i));
    		else
    			notselected.add(prots.get(i));
    	}
    	prots = selected;
    	prots.addAll(notselected);
    }
	
	public FASTAFile(String fn) throws Exception {
		StringBuffer cseq = new StringBuffer();
		String cname = "", cline;
		BufferedReader in = new BufferedReader(new FileReader(fn));
		
		prots = new ArrayList<Protein>();
		
		while((cline = in.readLine())!=null) {
			cline = cline.trim();
			if(cline.length() == 0)
				continue;
			if(cline.startsWith(">")) {
				if(cseq.length() > 0) {
					Protein np = new Protein();
					np.seqname = cname;
					np.seq = cseq.toString();
					prots.add(np);
				}
				cname = cline.substring(1);
				cseq = new StringBuffer();
			} else
				cseq.append(cline.trim());
		}
		if(cseq.length() > 0) {
			Protein np = new Protein();
			np.seqname = cname;
			np.seq = cseq.toString();
			prots.add(np);
		}
		in.close();

	}	
	
	public void writeOut(String f) throws Exception {
		PrintWriter out = new PrintWriter(f);
		for(int i = 0; i < prots.size(); i++) {
			Protein cp = prots.get(i);
			out.println(">"+cp.seqname);
			for(int j = 0; j < cp.seq.length(); j += 60) {
				out.println(cp.seq.substring(j,Math.min(j+60,cp.seq.length())));
			}
		}
		out.close();
	}
	
	public void validate() {
		String aas = "RHKDESTNQCGPAILMFWYV*";
		for(int i = 0; i < prots.size(); i++) {
			Protein cp = prots.get(i);
			boolean valid = true;
			for(int j = 0; j < cp.seq.length(); j++)
				if(aas.indexOf(cp.seq.charAt(j)) < 0)
					valid = false;
			if(!valid)
				System.out.println("Invalid protein: " + cp.seqname + "\n" + cp.seq);
		}
	}
	
}

class Protein {
	String seqname;
	String seq;
}
