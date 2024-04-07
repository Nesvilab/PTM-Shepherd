package edu.umich.andykong.ptmshepherd.paramhandling;

import java.util.HashMap;
import java.util.Map;

public class ParameterGroup {
    private String name;
    private Map<String, Parameter<?>> parameters;


    public ParameterGroup(String name) {
        this.name = name;
        this.parameters = new HashMap<>();
    }

    public void addParam(String key, Parameter<?> parameter) {
        parameters.put(key, parameter);
    }


    public Parameter<?> getParam(String key) {
        return parameters.get(key);
    }


    public <T> void setParamValue(String key, T value) {
        Parameter<T> parameter = (Parameter<T>) parameters.get(key);
        if (parameter == null)
            throw new IllegalArgumentException("Parameter not found: " + key);
        parameter.setValue(value);
    }



    public void printParameters() {
        System.out.println("Group: " + name);
        parameters.forEach((key, parameter) -> System.out.println(key + ": " + parameter.getValue()));
    }
}