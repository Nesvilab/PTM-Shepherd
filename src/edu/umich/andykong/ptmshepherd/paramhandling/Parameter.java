package edu.umich.andykong.ptmshepherd.paramhandling;

public interface Parameter<T> {
    T getValue();
    void setValue(T value) throws IllegalArgumentException;
    boolean isValid(T value);

}