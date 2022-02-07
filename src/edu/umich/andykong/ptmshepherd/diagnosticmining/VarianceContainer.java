package edu.umich.andykong.ptmshepherd.diagnosticmining;

// Follows Welford's online algorithm with M_2 modification for running value for numerical stability
// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
class VarianceContainer {
    private double sum;
    private int n;
    private double M2;
    private double mean;
    private double var;
    private double sampleVar;
    private boolean isValid = true;

    VarianceContainer() {
        this.sum = 0;
        this.M2 = 0;
        this.n = 0;
        this.mean = 0;
    }

    public void addVal(float val) {
        this.n++;
        double delta = val - this.mean;
        this.mean += delta / (double) this.n;
        double delta2 = val - this.mean;
        this.M2 += delta * delta2;
    }

    public void finalize() {
        if (this.n < 2)
            this.isValid = false;
        else {
            this.var = this.M2 / (double) this.n;
            this.sampleVar = this.M2 / (double) (this.n - 1);
        }
    }

    public double getSampleVar() { return this.sampleVar; }
    public double getPopulationVar() { return this.var; }
    public double getMean() { return this.mean; }
}