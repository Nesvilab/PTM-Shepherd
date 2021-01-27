package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.core.Spectrum;

import java.util.ArrayList;

public class BinDiagMetric {
    double peakApex;
    double leftBound;
    double rightBound;
    public ArrayList<Spectrum> spectra;
    public ArrayList<String> specNames;
    public DiagnosticHisto immoniumIons;
    public DiagnosticHisto tildeIons; //split into b and y?
    public DiagnosticHisto capYIons;

    public BinDiagMetric(double[] peakBounds) {
        this.peakApex = peakBounds[0];
        this.leftBound = peakBounds[1];
        this.rightBound = peakBounds[2];
        spectra = new ArrayList<>();
    }

    public void addPSMs(String specName) {
        specNames.add(specName);
    }
}