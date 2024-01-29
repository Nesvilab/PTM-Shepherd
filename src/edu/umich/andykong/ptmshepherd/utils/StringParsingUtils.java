package edu.umich.andykong.ptmshepherd.utils;

public class StringParsingUtils {
    public static String subString(String tar, String left, String right) {
        return tar.substring(tar.indexOf(left)+1, tar.indexOf(right));
    }
}
