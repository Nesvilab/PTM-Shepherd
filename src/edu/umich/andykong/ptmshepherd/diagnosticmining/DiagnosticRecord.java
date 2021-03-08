package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.core.Spectrum;

import java.util.HashMap;

/* This class holds the transformed peaklists of a single PSM
*  Each instance of this will be entered into a binary diagnostic output file */
public class DiagnosticRecord implements Comparable<DiagnosticRecord>  {
    public int scanNum;
    private float[][] immoniumPeaks;
    private float[][] capYPeaks;
    private HashMap<Character, float[][]> squigglePeaks;
    private String ionTypes;
    private String pepSeq;
    private String modifications;
    private float dmass;

    public int compareTo(DiagnosticRecord dr) {
        return -1 * Integer.compare(this.scanNum, dr.scanNum);
    };

    public DiagnosticRecord(Spectrum spec, String ionTypes, String pepSeq, String modifications, float dmass) {
        this.scanNum = spec.scanNum;
        this.ionTypes = ionTypes;
        this.pepSeq = pepSeq;
        this.modifications = modifications;
        this.dmass = dmass;
    }

    public void setImmoniumPeaks(float[][] immoniumPeaks) {
        this.immoniumPeaks = immoniumPeaks;
    }

    public void setCapYPeaks(float[][] capYPeaks) {
        this.capYPeaks = capYPeaks;
    }

    public void setSquigglePeaks(HashMap<Character, float[][]> squigglePeaks) {
        this.squigglePeaks = squigglePeaks;
    }
}
