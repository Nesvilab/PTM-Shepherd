package edu.umich.andykong.ptmshepherd.core;

import java.util.HashMap;
import java.util.Map;

public class AAMutationRules {
    public static Map<Character, Character> mutateFrom = new HashMap<>();

    // These rules are taken from DIA-NN
    static {
        mutateFrom.put('G', 'L');
        mutateFrom.put('A', 'L');
        mutateFrom.put('V', 'L');
        mutateFrom.put('L', 'V');
        mutateFrom.put('I', 'V');
        mutateFrom.put('F', 'L');
        mutateFrom.put('M', 'L');
        mutateFrom.put('P', 'L');
        mutateFrom.put('W', 'L');
        mutateFrom.put('S', 'T');
        mutateFrom.put('C', 'S');
        mutateFrom.put('T', 'S');
        mutateFrom.put('Y', 'S');
        mutateFrom.put('H', 'S');
        mutateFrom.put('K', 'L');
        mutateFrom.put('R', 'L');
        mutateFrom.put('Q', 'N');
        mutateFrom.put('E', 'D');
        mutateFrom.put('N', 'Q');
        mutateFrom.put('D', 'E');
    }
}
