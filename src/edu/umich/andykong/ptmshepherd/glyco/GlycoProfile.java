package edu.umich.andykong.ptmshepherd.glyco;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.specsimilarity.SimRTRecord;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;

public class GlycoProfile {
    public GlycoRecord [] records;
    public FastLocator locate;
    double [] masses;
    double peakTol;
    int precursorMassUnits;
    String[] capYStrs, diagIonStrs, remFragStrs;


    public GlycoProfile(double[] massRanges, int precMassUnits, double peakTol) {
        masses = Arrays.copyOf(massRanges, massRanges.length);
        this.peakTol = peakTol;
        locate = new FastLocator(masses, peakTol);
        records = new GlycoRecord[masses.length];
        for(int i = 0; i < masses.length; i++)
            records[i] = new GlycoRecord(masses[i], i);
        //get cap y, ox ions
        //String[] capYStrs, diagIonStrs, remFragStrs;
        //cap y ions
        if (PTMShepherd.getParam("cap_y_ions").length() > 0)
            capYStrs = PTMShepherd.getParam("cap_y_ions").split(",");
        else
            capYStrs = new String[0];
        //oxonium ions
        if (PTMShepherd.getParam("diag_ions").length() > 0)
            diagIonStrs = PTMShepherd.getParam("diag_ions").split(",");
        else
            diagIonStrs = new String[0];
        //remainder masses
        if (PTMShepherd.getParam("remainder_masses").length() > 0)
            remFragStrs = PTMShepherd.getParam("remainder_masses").split(",");
        else
            remFragStrs = new String[0];
    }

    public void writeProfile(String path) throws Exception {
        PrintWriter out = new PrintWriter(new FileWriter(path));
        //format header
        StringBuffer sb = new StringBuffer();
        sb.append("Peak\tMatched PSMs");
        for (int i = 0; i < capYStrs.length; i++)
            sb.append("\tY " + capYStrs[i] + " (PSMs)");
        for (int i = 0; i < diagIonStrs.length; i++)
            sb.append("\tox " + diagIonStrs[i] + " (PSMs)");
        for (int i = 0; i < remFragStrs.length; i++)
            sb.append("\tremainder "+remFragStrs[i] + " (PSMs)");
        for (int i = 0; i < capYStrs.length; i++)
            sb.append("\tY " + capYStrs[i] + " (% PSMs)");
        for (int i = 0; i < diagIonStrs.length; i++)
            sb.append("\tox " + diagIonStrs[i] + " (% PSMs)");
        for (int i = 0; i < remFragStrs.length; i++)
            sb.append("\tremainder "+remFragStrs[i] + " (% PSMs)");
        out.println(sb.toString());
        for(int i = 0; i < records.length; i++) {
            out.println(records[i].toString());
        }
        out.close();
    }

}
