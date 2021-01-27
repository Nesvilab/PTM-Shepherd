package edu.umich.andykong.ptmshepherd.diagnosticmining;

public class DiagnosticHisto {
    double [] bins;
    double min;
    double max;
    double buffersize = 5; //pads endings to mirror internal bins
    double binsPerDa;

    public DiagnosticHisto(double min, double max, float binWidth) {
        this.min = min;
        this.max = max;
        this.binsPerDa = (int) (1.0 / binWidth);
        System.out.println(this.binsPerDa);
        this.bins = new double [(int)(binsPerDa * (max - min + buffersize))]; //TODO
    }

    private int locateBin(double val) {
        int bini = -1;
        if (val < max && val > min)
            bini = (int)((val - min + buffersize) * binsPerDa);
        return bini;
    }
}
