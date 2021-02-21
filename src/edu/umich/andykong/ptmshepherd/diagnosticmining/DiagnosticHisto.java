package edu.umich.andykong.ptmshepherd.diagnosticmining;

public class DiagnosticHisto {
    double [] bins;
    double min;
    double max;
    double buffersize = 1.0; //pads endings to mirror internal bins
    double binsPerDa;

    public DiagnosticHisto (double mn, double mx, double binWidth) {
        this.min = mn - buffersize;
        this.max = mx + buffersize;
        //System.out.println(this.max);
        //System.out.println(this.min);
        this.binsPerDa = (int) (1.0 / binWidth);
        //System.out.println(1.0 / binWidth);
        //System.out.println((int) (1.0 / binWidth));
        //System.out.println((int)(this.binsPerDa * (this.max - this.min)));
        this.bins = new double [(int)(this.binsPerDa * (max - min))]; //TODO test
    }

    private int locateBin(double val) {
        int bini = -1;
        bini = (int) ((val - this.min) * this.binsPerDa);
        return bini;
    }

    public void placeIon(double mz, double intensity) {
        try {
            this.bins[locateBin(mz)] += intensity;
        } catch (Exception e) { //todo purge this
            System.out.printf("%.4f\t%.4f\t%d***\n", this.min, this.max, this.bins.length);
            System.out.printf("%.4f\t%d\t%.4f\n", mz, locateBin(mz), intensity);
        }
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
}
