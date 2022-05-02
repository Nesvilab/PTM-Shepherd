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

public class PrefixMover {

	static TreeMap<String, ArrayList<File>> prefixMap;
	
	public static void main(String[] args) throws Exception {
		prefixMap = new TreeMap<>();
		File baseDir = new File(args[0]);
		File [] ls = baseDir.listFiles();
		for(int i = 0; i < ls.length; i++) {
			if(ls[i].getName().endsWith(args[1])) {
				String fbase = ls[i].getName();
				fbase = fbase.substring(0,fbase.indexOf("."));
				String prefix = fbase.substring(0, fbase.lastIndexOf("_"));
				if(!prefixMap.containsKey(prefix))
					prefixMap.put(prefix, new ArrayList<>());
				prefixMap.get(prefix).add(ls[i]);
			}
		}
		for(String s : prefixMap.keySet()) {
			ArrayList<File> files = prefixMap.get(s);
			File target = baseDir.toPath().resolve(s).toFile();
			target.mkdirs();
			for(File f : files) {
				f.renameTo(target.toPath().resolve(f.getName()).toFile());
			}
			System.out.printf("%s -> %d files\n", s, files.size());
		}
	}
}
