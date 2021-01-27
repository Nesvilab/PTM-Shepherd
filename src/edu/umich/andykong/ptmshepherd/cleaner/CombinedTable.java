package edu.umich.andykong.ptmshepherd.cleaner;

import com.google.common.base.Charsets;
import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CombinedTable {

    public static void writeCombinedTable(String dataset) throws IOException {
        List<String> subTabs = new ArrayList<>();
        List<List<String>> records = new ArrayList<>();
        subTabs.add("peaksummary.annotated.tsv");
        subTabs.add(dataset + ".simrtprofile.txt");
        subTabs.add(dataset + ".locprofile.txt");


        Set<Integer> observedLineCounts = new HashSet<>();
        int maxLineLen = 0;

        for (int subTabIdx = 0; subTabIdx < subTabs.size(); subTabIdx++) {
            String subTab = subTabs.get(subTabIdx);
            Path p = Paths.get(subTab);
            System.out.println("Reading file: " + p.toString());

            // read the whole file
            List<String> lines = readLines(p);

            // check for consistency, number of lines must be the same in all files
            observedLineCounts.add(lines.size());
            if (observedLineCounts.size() > 1) {
                throw new IllegalStateException(String.format("The number of lines [%d] in file "
                    + "[%s] does not match what we've seen in other files [%d]",
                    lines.size(), subTab, observedLineCounts.iterator().next()));
            }

            // update
            if  (subTabIdx == 0) {
                maxLineLen = Arrays.asList(lines.get(0).split("\t")).size();
                for (int lineIdx = 0; lineIdx < lines.size(); lineIdx++) {
                    String line = lines.get(lineIdx);
                    if (line == null) {
                        continue;
                    }
                    List<String> split = new ArrayList<String>(Arrays.asList(line.split("\t", -1)));
                    //if (split.size() > maxLineLen)
                    //    maxLineLen = split.size();
                    List<String> newLine = split.subList(0, 4); //peak values
                    newLine.add(split.get(split.size()-3)); //total PSMs
                    newLine.add(split.get(split.size()-4)); //total % in unmodified
                    newLine.add(split.get(split.size()-2)); //mod 1
                    newLine.add(split.get(split.size()-1)); //mod 2
                    records.add(newLine);
                    //records.add(newLine); //add first 4 fields to new file
                    //List<String> newLine2 = split.subList(maxLineLen-2, split.size()); //get mod annotationsm $

                    //int padding = 2 - newLine2.size();
                    //List<String> crec = records.get(lineIdx);
                    //crec.addAll(newLine2);
                    //for(int i = 0; i < padding; i++){
                    //    crec.add("");
                    //}
                    //crec.add(split.get(maxLineLen-1));
                    //records.set(lineIdx, crec);
                }
            } else if (subTabIdx == 2) {
                for (int lineIdx = 0; lineIdx < lines.size(); lineIdx++) {
                    String line = lines.get(lineIdx);
                    if (line == null) {
                        continue;
                    }
                    List<String> split = new ArrayList<String>(Arrays.asList(line.split("\t")));
                    List<String> newLine = split.subList(1, 7);
                    newLine.remove(1);
                    List<String> crec = records.get(lineIdx);
                    crec.addAll(newLine);
                    records.set(lineIdx, crec);
                }
            }
            else {
                for (int lineIdx = 0; lineIdx < lines.size(); lineIdx++) {
                    String line = lines.get(lineIdx);
                    if (line == null) {
                        continue;
                    }
                    List<String> split = new ArrayList<String>(Arrays.asList(line.split("\t")));
                    List<String> newLine = split.subList(2, split.size());
                    List<String> crec = records.get(lineIdx);
                    crec.addAll(newLine);
                    records.set(lineIdx, crec);
                }
            }

//
//            else {
//                for (int lineIdx = 0; lineIdx < lines.size(); lineIdx++) {
//                    String line = lines.get(lineIdx);
//                    if (line == null) {
//                        continue;
//                    }
//                    List<String> split = new ArrayList<String>(Arrays.asList(line.split("\t")));
//                    List<String> newLine = split.subList(1, 7);
//                    List<String> crec = records.get(lineIdx);
//                    crec.addAll(newLine);
//                    records.set(lineIdx, crec);
//
//                    String line = lines.get(lineIdx);
//                    if (line == null) {
//                        continue;
//                    }
//                    List<String> split = new ArrayList<String>(Arrays.asList(line.split("\t")));
//                    List<String> newLine = split.subList(2, split.size());
//                    String lineCut = String.join("\t", newLine);
//                    records.set(lineIdx, String.format("%s%s%s", records.get(lineIdx), "\t", lineCut));
//                }
//            }
        }


        String subTab = subTabs.get(0);
        Path p = Paths.get(subTab);
        List<String> lines = readLines(p);
        for (int lineIdx = 0; lineIdx < lines.size(); lineIdx++) {
            String line = lines.get(lineIdx);
            if (line == null) {
                continue;
            }
            List<String> split = new ArrayList<String>(Arrays.asList(line.split("\t")));
            List<String> newLine = split.subList(4, maxLineLen-4);
            List<String> crec = records.get(lineIdx);
            crec.addAll(newLine);
            records.set(lineIdx, crec);
        }

        String fnOut = dataset + ".profile.tsv";
        Path pOut = Paths.get(fnOut);
        System.out.println("Writing to file: " + pOut.toString());
        try (BufferedWriter bw = Files.newBufferedWriter(pOut)) {
            for (List record: records) {
                String outLine = String.join("\t", record);
                bw.write(outLine + System.lineSeparator());
            }
        }
        //records.forEach(record ->{System.out.println(record);});
    }

    private static List<String> readLines(Path path) throws IOException {
        List<String> lines = new ArrayList<>();
        try (InputStream is = Files.newInputStream(path)) {
            BufferedReader br = new BufferedReader(new InputStreamReader(is, "UTF-8"));
            String line;
            while ((line = br.readLine()) != null) {
                lines.add(line);
            }
        }
        return lines;
    }
}
