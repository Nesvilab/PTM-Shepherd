package edu.umich.andykong.ptmshepherd.glyco;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.FastLocator;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.specsimilarity.SimRTRecord;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;

public class GlycoProfile {
    public GlycoRecord [] records;
    public FastLocator locate;
    double [] masses;
    double peakTol;


    public GlycoProfile(double[] massRanges, double peakTol) {
        masses = Arrays.copyOf(massRanges, massRanges.length);
        this.peakTol = peakTol;
        locate = new FastLocator(masses, peakTol);
        records = new GlycoRecord[masses.length];

        for(int i = 0; i < masses.length; i++)
            records[i] = new GlycoRecord(masses[i], i);
    }


}
