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
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.util.*;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.w3c.dom.*;

import javax.xml.parsers.*;;

public class UnimodParser {

	public static void loadUniMod(String path) throws Exception {
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder(); 
		Document doc = db.parse(new File(path));
		NodeList nl = doc.getChildNodes();
		Node unimod = null;
		Node modifications = null;
		for(int i = 0; i < nl.getLength(); i++) {
			Node n = nl.item(i);
			if(n.getNodeName().equals("unimod"))
				unimod = n;
		}
		nl = unimod.getChildNodes();
		for(int i = 0; i < nl.getLength(); i++) {
			Node n = nl.item(i);
			if(n.getNodeName().equals("modifications"))
				modifications = n;
		}
		
		NodeList nds = modifications.getChildNodes();
		for(int i = 0; i < nds.getLength(); i++) {
			Node n = nds.item(i);
			if(!n.getNodeName().equals("modifications_row"))
				continue;
			Element e = (Element)n;
			System.out.println(e.getAttribute("full_name") +"\t"+ e.getAttribute("mono_mass"));
		}
	}
	
	public static void main(String [] args) throws Exception {
		loadUniMod("d:\\unimod_tables.xml");
	}
	
}
