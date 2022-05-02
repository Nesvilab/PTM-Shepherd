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

package edu.umich.andykong.ptmshepherd.utils;

import java.io.*;
import java.util.*;
import java.nio.file.Files;

public class PSMLocator {

	public static void main(String [] args) throws Exception {
		File [] ls = new File(".").listFiles();
		int PNNL = 1, JHU = 1;
		for(int i = 0; i < ls.length; i++) {
			if(!ls[i].getName().startsWith("TCGA"))
				continue;
			File psmTSV = ls[i].toPath().resolve("psm.tsv").toFile();
			System.out.print(ls[i].getName());
			if(!psmTSV.exists())
				System.out.println(" FAIL");
			else {
				System.out.println(" OK");
				File target = new File(ls[i].getName()+".tsv");
				if(!target.exists())
					Files.copy(psmTSV.toPath(), target.toPath());
				if(ls[i].getName().indexOf("JHU") >= 0) {
					System.err.printf("dataset = JHU_Set%d %s /scratch/datasets/CPTAC2_Ovarian\n", JHU++, target.getAbsolutePath());
				} else {
					System.err.printf("dataset = PNNL_Set%d %s /scratch/datasets/CPTAC2_Ovarian\n", PNNL++, target.getAbsolutePath());
				}
			}
		}
	}
	
}
