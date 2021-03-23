package edu.umich.andykong.ptmshepherd.glyco;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.specsimilarity.SimRTRecord;
import edu.umich.andykong.ptmshepherd.specsimilarity.Variance;

public class GlycoRecord {
    double mass;
    int originalOrder;
    int count;
    String [] capYStrs;
    String [] diagIonStrs;
    String [] remFragStrs;
    int [] capYCounts;
    int [] diagIonCounts;
    int [] remFragCounts;

    public GlycoRecord(double mass, int originalOrder) {
        this.mass = mass;
        this.originalOrder = originalOrder;
        count = 0;

        //get cap y, ox ions
        //cap y ions
        if (PTMShepherd.getParam("cap_y_ions").length() > 0)
            capYStrs = PTMShepherd.getParam("cap_y_ions").split(",| |/");
        else
            capYStrs = new String[0];
        this.capYCounts = new int[capYStrs.length];
        //oxonium ions
        if (PTMShepherd.getParam("diag_ions").length() > 0)
            diagIonStrs = PTMShepherd.getParam("diag_ions").split(",| |/");
        else
            diagIonStrs = new String[0];
        this.diagIonCounts = new int[diagIonStrs.length];
        //remainder masses
        if (PTMShepherd.getParam("remainder_masses").length() > 0)
            remFragStrs = PTMShepherd.getParam("remainder_masses").split(",| |/");
        else
            remFragStrs = new String[0];
        this.remFragCounts = new int[remFragStrs.length];
    }

    public void updateWithLine(String [] sp) {
        count++;
        double cInt;
        //get instances of Y ion identified
        int startCol = 8;   // 8 columns always padding left side
        int endCol = startCol + capYCounts.length;
        for(int i = startCol; i < endCol; i++) {
            cInt = Double.parseDouble(sp[i]);
            if (cInt > 0.0) {
                //System.out.println("*");
                //System.out.println(this.capYCounts[i-startCol]);
                this.capYCounts[i-startCol]++;
                //System.out.println(this.capYCounts[i-startCol]);
            }
        }
        //get instances of diagnostic ion identified
        startCol = endCol;
        endCol = startCol + diagIonCounts.length;
        //System.out.printf("%d\t%d\t11111\n",startCol,endCol);
        for(int i = startCol; i < endCol; i++) {
            cInt = Double.parseDouble(sp[i]);
            if (cInt > 0.0) {
                this.diagIonCounts[i-startCol]++;
            }
        }
        //get instances of remainder frag identified
        startCol = endCol;
        endCol = sp.length;
        for(int i = startCol; i < endCol; i++){
            cInt = Double.parseDouble(sp[i]); //is actually delta score, not intensity
            if (cInt > 0.0) {
                int t = i;
                t = t - startCol;
                t = t - (i - startCol)/2;
                this.remFragCounts[t]++;
                //this.remFragCounts[i-(startCol + (i - startCol) / 2)]++;
            }
            i++;
        }
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        //General stats
        sb.append(String.format("%.4f\t%d", mass, count));
        //append y ion stats
        for(int i = 0; i < capYCounts.length; i++)
            sb.append(String.format("\t%d", capYCounts[i]));
        //append diagnostic ion stats
        for(int i = 0; i < diagIonCounts.length; i++)
            sb.append(String.format("\t%d", diagIonCounts[i]));
        //append remainder ion stats
        for(int i = 0; i < remFragCounts.length; i++)
            sb.append(String.format("\t%d", remFragCounts[i]));
        //append y ion stats
        for(int i = 0; i < capYCounts.length; i++)
            sb.append(String.format("\t%.1f", (float)capYCounts[i] * 100.0 / (float)count));
        //append diagnostic ion stats
        for(int i = 0; i < diagIonCounts.length; i++)
            sb.append(String.format("\t%.1f", (float)diagIonCounts[i] * 100.0 / (float)count));
        //append remainder ion stats
        for(int i = 0; i < remFragCounts.length; i++)
            sb.append(String.format("\t%.1f", (float) remFragCounts[i] * 100.0 / (float) count));

        return sb.toString();
    }

}
