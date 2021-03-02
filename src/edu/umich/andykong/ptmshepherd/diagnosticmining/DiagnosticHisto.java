package edu.umich.andykong.ptmshepherd.diagnosticmining;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class DiagnosticHisto {
    double [] bins;
    double min;
    double max;
    double buffersize = 1.0; //pads endings to mirror internal bins
    double binsPerDa;

    public DiagnosticHisto (double mn, double mx, double binWidth) {
        this.min = mn - buffersize;
        this.max = mx + buffersize;
        this.binsPerDa = (int) (1.0 / binWidth);
        this.bins = new double [(int)(this.binsPerDa * (max - min))]; //TODO test
    }

    private int locateBin(double val) {
        int bini = -1;
        bini = (int) ((val - this.min) * this.binsPerDa);
        return bini;
    }

    public void placeIon(double mz, double intensity) {
        this.bins[locateBin(mz)] += intensity;
    }

    public double binToMass(int bini) {
        return (((double) bini / (double) this.binsPerDa) + this.min);
    }

    public int findMax() {
        double max = -1;
        int maxi = -1;
        for (int i = 0; i < this.bins.length; i++) {
            if (this.bins[i] > max) {
                max = this.bins[i];
                maxi = i;
            }
        }
        return maxi;
    }

    public void printHisto(String fname) throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(new File(fname)));
        out.printf("mass\theight\n");
        for (int i = 0; i < this.bins.length; i++) {
            StringBuffer sb = new StringBuffer();
            for (int j = 0; j < 10000; j++) {
                if (i + j < this.bins.length) {
                    sb.append(String.format("%.04f\t%.04f\n", binToMass(i), this.bins[i]));
                    i++;
                } else {
                    continue;
                }
            }
            out.printf(sb.toString());
        }
        out.close();
    }
}
