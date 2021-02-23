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
    double [][] peaks;
    double peakTol;
    int precursorUnits;
    String[] capYStrs, diagIonStrs, remFragStrs;


    public GlycoProfile(double[][] peakVals, int precursorUnits, double peakTol) {
        this.masses = new double[peakVals[0].length];
        for (int i = 0; i < peakVals[0].length; i++) {
            this.masses[i] = peakVals[0][i];
            //System.out.println(this.masses[i]);
        }
        peaks = peakVals;
        this.peakTol = peakTol;
        this.precursorUnits = precursorUnits;
        locate = new FastLocator(peaks, peakTol, precursorUnits);
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
        sb.append("Peak\tPSMs");
        for (int i = 0; i < capYStrs.length; i++)
            sb.append("\tY_" + capYStrs[i] + "_(PSMs)");
        for (int i = 0; i < diagIonStrs.length; i++)
            sb.append("\tox_" + diagIonStrs[i] + "_(PSMs)");
        for (int i = 0; i < remFragStrs.length; i++)
            sb.append("\tremainder_"+remFragStrs[i] + "_(PSMs)");
        for (int i = 0; i < capYStrs.length; i++)
            sb.append("\tY_" + capYStrs[i] + "_(percent_of_PSMs)");
        for (int i = 0; i < diagIonStrs.length; i++)
            sb.append("\tox_" + diagIonStrs[i] + "_(percent_of_PSMs)");
        for (int i = 0; i < remFragStrs.length; i++)
            sb.append("\tremainder_"+remFragStrs[i] + "_(percent_of_PSMs)");
        out.println(sb.toString());
        for(int i = 0; i < records.length; i++) {
            out.println(records[i].toString());
        }
        out.close();
    }

}
