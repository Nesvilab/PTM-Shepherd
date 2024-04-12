package edu.umich.andykong.ptmshepherd.utils;

import java.util.ArrayList;
import java.util.List;

public class StringParsingUtils {
    public static String subString(String tar, String left, String right) {
        return tar.substring(tar.indexOf(left)+1, tar.indexOf(right));
    }

    public static List<String> splitStringByTab(String input) {
        List<String> result = new ArrayList<>();
        int start = 0;
        int end = input.indexOf('\t');
        while (end != -1) {
            result.add(input.substring(start, end));
            start = end + 1;
            end = input.indexOf('\t', start);
        }
        result.add(input.substring(start)); // Add the last segment
        return result;
    }

}
