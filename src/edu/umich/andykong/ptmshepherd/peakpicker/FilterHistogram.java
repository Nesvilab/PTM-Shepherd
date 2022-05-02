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

public class FilterHistogram {

	public static void filter(String input, String output, double lo, double hi) throws Exception {
		BufferedReader in = new BufferedReader(new FileReader(input));
		PrintWriter out = new PrintWriter(new FileWriter(output));
		String cline;
		out.println(in.readLine());
		while((cline = in.readLine()) != null) {
			String [] sp = cline.split("\t");
			double v = Double.parseDouble(sp[0]);
			if(v >= lo && v <= hi)
				out.println(cline);
		}
		in.close();
		out.close();
	}
	
	public static void main(String [] args) throws Exception {
		//filter("E:\\OpenPipeline\\combined.tsv","E:\\OpenPipeline\\deamidation.tsv",0.97,1.01);
		//filter("E:\\OpenPipeline\\combined.tsv","E:\\OpenPipeline\\formylation.tsv",27.98,28.02);
		filter("E:\\OpenPipeline\\combined.tsv","E:\\OpenPipeline\\phosphorylation.tsv",79.955,79.975);
	}
	
}
