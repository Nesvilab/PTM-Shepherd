package edu.umich.andykong.ptmshepherd.paramhandling;

import java.util.HashMap;
import java.util.Map;

public class ParamHandlerSingleton {
    private static ParamHandlerSingleton instance;
    private Map<String, ParameterGroup> groups;
    //todo make a map here that links to parameter keys to group parameters


    public static synchronized ParamHandlerSingleton getInstance() {
        if (instance == null) {
            instance = new ParamHandlerSingleton();
        }
        return instance;
    }


    public void addParamGroup(String key) {
        try {
            if (groups.containsKey(key)) {
                throw new IllegalArgumentException("Parameter group already exists: " + key);
            }
            groups.put(key, new ParameterGroup(key));
        } catch (IllegalArgumentException e) {
            return;
        }
    }


    @SuppressWarnings("unchecked")
    public ParameterGroup getParamGroup(String key) {
        ParameterGroup parameterGroup = groups.get(key);
        if (parameterGroup == null) {
            throw new IllegalArgumentException("Parameter group not found: " + key);
        }
        return parameterGroup;
    }

    public <T> void addParam(String group, String key, Parameter<T> param) {
        ParameterGroup parameterGroup = groups.get(group);
        parameterGroup.setParamValue(key, param);
    }

    public <T> void setParamValue(String group, String key, T value) {
        ParameterGroup parameterGroup = groups.get(group);
        parameterGroup.setParamValue(key, value);
    }


}