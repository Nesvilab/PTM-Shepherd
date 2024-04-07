package edu.umich.andykong.ptmshepherd.paramhandling;


public class IntegerParameter implements Parameter<Integer> {
    private String key;
    private int level;
    private int order;
    private int value;
    private int min;
    private int max;
    private int defaultValue;
    private String description;

    public IntegerParameter(String key, int level, int order, int min, int max, int value, String description) {
        this.key = key;
        this.level = level;
        this.order = order;
        this.min = min;
        this.max = max;
        this.value = value;
        this.defaultValue = value;
        this.description = description;
    }

    @Override
    public Integer getValue() {
        return value;
    }

    @Override
    public void setValue(Integer value) throws IllegalArgumentException {
        if (!isValid(value)) {
            throw new IllegalArgumentException(String.format(this.key + " received an invalid argument."));
        }
        this.value = value;
    }

    @Override
    public boolean isValid(Integer value) {
        return ((value >= min) && (value <= max));
    }


}

