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
        List<String> records = new ArrayList<>();
        subTabs.add("peaksummary.annotated.tsv");
        subTabs.add(dataset + ".simrtprofile.txt");
        subTabs.add(dataset + ".locprofile.txt");


        Set<Integer> observedLineCounts = new HashSet<>();

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
            if (subTabIdx == 0) {
                records.addAll(lines);
            } else if (subTabIdx == 2) {
                for (int lineIdx = 0; lineIdx < lines.size(); lineIdx++) {
                    String line = lines.get(lineIdx);
                    if (line == null) {
                        continue;
                    }
                    String[] split = line.split("\t");
                    String lineCut = String.join("\t", Arrays.asList(split).subList(1, 7));
                    records.set(lineIdx, String.format("%s%s%s", records.get(lineIdx), "\t", lineCut));
                }
            }
            else {
                for (int lineIdx = 0; lineIdx < lines.size(); lineIdx++) {
                    String line = lines.get(lineIdx);
                    if (line == null) {
                        continue;
                    }
                    String[] split = line.split("\t");
                    String lineCut = String.join("\t", Arrays.asList(split).subList(1, split.length));
                    records.set(lineIdx, String.format("%s%s%s", records.get(lineIdx), "\t", lineCut));
                }
            }
        }


        String fnOut = dataset + ".profile.tsv";
        Path pOut = Paths.get(fnOut);
        System.out.println("Writing to file: " + pOut.toString());
        try (BufferedWriter bw = Files.newBufferedWriter(pOut)) {
            for (String record: records) {
                bw.write(record + System.lineSeparator());
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
