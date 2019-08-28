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
