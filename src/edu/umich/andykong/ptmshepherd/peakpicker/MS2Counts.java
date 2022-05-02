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
import java.util.*;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.MSFMGFFile;
import edu.umich.andykong.ptmshepherd.core.MZBINFile;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
import umich.ms.fileio.filetypes.mzxml.MZXMLFile;
import umich.ms.fileio.filetypes.thermo.ThermoRawFile;

public class MS2Counts {

	public static int countMS2Scans(File f, int threads) throws Exception {
		int count = 0;
		String ext = f.getName().substring(f.getName().lastIndexOf(".")+1).toLowerCase();
		if(ext.length() == 0)
			return 0;
		try {
			if (ext.equals("mzxml")) {
				MZXMLFile source = new MZXMLFile(f.getAbsolutePath());
				source.setNumThreadsForParsing(threads);
				source.setExcludeEmptyScans(false);
				ScanCollectionDefault scans = new ScanCollectionDefault();
				scans.setDataSource(source);
				scans.loadData(LCMSDataSubset.STRUCTURE_ONLY);
				TreeMap<Integer, IScan> num2scan = scans.getMapNum2scan();
				Set<Map.Entry<Integer, IScan>> scanEntries = num2scan.entrySet();
				for (Map.Entry<Integer, IScan> scanEntry : scanEntries) {
					IScan scan = scanEntry.getValue();
					if (scan.getMsLevel() == 2)
						count++;
				}
				scans.reset();
			} else if (ext.equals("mzml")) {
				MZMLFile source = new MZMLFile(f.getAbsolutePath());
				source.setNumThreadsForParsing(threads);
				source.setExcludeEmptyScans(false);
				ScanCollectionDefault scans = new ScanCollectionDefault();
				scans.setDataSource(source);
				scans.loadData(LCMSDataSubset.STRUCTURE_ONLY);
				TreeMap<Integer, IScan> num2scan = scans.getMapNum2scan();
				Set<Map.Entry<Integer, IScan>> scanEntries = num2scan.entrySet();
				for (Map.Entry<Integer, IScan> scanEntry : scanEntries) {
					IScan scan = scanEntry.getValue();
					if (scan.getMsLevel() == 2)
						count++;
				}
				scans.reset();
			} else if (ext.equals("raw")) {
				ThermoRawFile source = new ThermoRawFile(f.getAbsolutePath());
				source.setNumThreadsForParsing(threads);
				source.setExcludeEmptyScans(false);
				ScanCollectionDefault scans = new ScanCollectionDefault();
				scans.setDataSource(source);
				scans.loadData(LCMSDataSubset.MS2_WITH_SPECTRA);
				TreeMap<Integer, IScan> num2scan = scans.getMapNum2scan();
				Set<Map.Entry<Integer, IScan>> scanEntries = num2scan.entrySet();
				for (Map.Entry<Integer, IScan> scanEntry : scanEntries) {
					IScan scan = scanEntry.getValue();
					if (scan.getMsLevel() == 2)
						count++;
				}
				scans.reset();
				source.close();
			} else if (ext.equals("d")) {
				String tempMzbFp = f.getAbsolutePath().replaceFirst(".d$", ".mzBIN");
				File tempMzbF = new File(tempMzbFp);
				if (tempMzbF.exists()) {
					//String f = tempMzbFp;
					MZBINFile source = new MZBINFile(PTMShepherd.executorService, Integer.parseInt(PTMShepherd.getParam("threads")), tempMzbFp, true);
				} else {
					System.out.println("Cannot read .d files without associated .mzBIN");
					System.exit(1);
				}
			} else if (ext.equals("mzBIN")) {
				try {
					MZBINFile source = new MZBINFile(PTMShepherd.executorService, Integer.parseInt(PTMShepherd.getParam("threads")), f.getAbsolutePath(), true);
					for (Spectrum spec : source.specs) {
						if (spec.msLevel == 2) {
							count++;
						}
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
			} else if (ext.equals("mgf")) {
				try {
					MSFMGFFile source = new MSFMGFFile(PTMShepherd.executorService, Integer.parseInt(PTMShepherd.getParam("threads")), f.getAbsolutePath(), false);
					count = source.maxScan;
				} catch (Exception e) {
					e.printStackTrace();
				}
			} else {
				System.err.println("Unrecognized extension: " + ext);
			}
		} catch (Exception e) {
			System.err.println("Error: " + e.toString());
		}
		return count;
	}
	
	public static void main(String [] args) throws Exception {
//		System.out.println(countMS2Scans(new File("E:\\q01507.mzXML")));
	}
	
}
