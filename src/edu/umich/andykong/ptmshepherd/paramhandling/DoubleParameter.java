package edu.umich.andykong.ptmshepherd.paramhandling;


public class DoubleParameter implements Parameter<Double> {
    private String key;
    private double value;
    private double min;
    private double max;

    public DoubleParameter(String key, double min, double max) {
        this.key = key;
        this.min = min;
        this.max = max;
    }

    @Override
    public Double getValue() {
        return this.value;
    }

    @Override
    public void setValue(Double value) {
        if (!isValid(value)) {
            throw new IllegalArgumentException("Value out of bounds");
        }
        this.value = value;
    }

    @Override
    public boolean isValid(Double value) {
        return ((value >= this.min) && (value <= this.max));
    }
}
