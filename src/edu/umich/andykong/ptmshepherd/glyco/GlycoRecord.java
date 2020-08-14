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
            capYStrs = PTMShepherd.getParam("cap_y_ions").split(",");
        else
            capYStrs = new String[0];
        capYCounts = new int[capYStrs.length];
        //oxonium ions
        if (PTMShepherd.getParam("diag_ions").length() > 0)
            diagIonStrs = PTMShepherd.getParam("diag_ions").split(",");
        else
            diagIonStrs = new String[0];
        diagIonCounts = new int[diagIonStrs.length];
        //remainder masses
        if (PTMShepherd.getParam("remainder_masses").length() > 0)
            remFragStrs = PTMShepherd.getParam("remainder_masses").split(",");
        else
            remFragStrs = new String[0];
        remFragCounts = new int[remFragStrs.length];
    }

    public void updateWithLine(String [] sp) {
        count++;
    }

}
