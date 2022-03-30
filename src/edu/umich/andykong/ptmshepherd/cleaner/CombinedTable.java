package edu.umich.andykong.ptmshepherd.cleaner;

import com.google.common.base.Charsets;
import java.io.*;
import java.lang.reflect.Array;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CombinedTable {
    public ArrayList<StringBuffer> data;
    public ArrayList<String> headers;
    public String fname;
    public String dataset;

    public CombinedTable(String ds) {
        this.dataset = ds;
        this.data = new ArrayList<>();
        this.headers = new ArrayList<>();
        this.fname = PTMShepherd.normFName(dataset + ".profile.tsv");
    }

    public void writeCombinedTable(int useIntensities, boolean calcIntensities) throws IOException {
        /* Process preaksummary.annotated.tsv file */

        /* Get headers that we're adding later */
        BufferedReader in = new BufferedReader(new FileReader(new File(
                PTMShepherd.normFName("peaksummary.annotated.tsv"))));
        String[] curHeaders = in.readLine().split("\t", -1);
        ArrayList<String> experiments = new ArrayList<>();
        for (int i = 0; i < curHeaders.length; i++) {
            if (curHeaders[i].endsWith("_percent_PSMs"))
                experiments.add(curHeaders[i].substring(0, curHeaders[i].indexOf("_percent_PSMs")));
        }

        /* Get cols to add to beginning from peaksummary.annotated.tsv */
        String[] colsToAdd = new String[]{"peak_apex", "peak_lower", "peak_upper", "PSMs", "peak_signal",
                    "percent_also_in_unmodified", "mapped_mass_1", "mapped_mass_2"};


        ArrayList<String> colsToAddLater = new ArrayList<>();

        /* Add headers */
        for (String col : colsToAdd)
            this.headers.add(col);

        /* Assure proper order for elements */
        for (String exp : experiments)
            colsToAddLater.add(exp + "_PSMs");
        for (String exp : experiments)
            colsToAddLater.add(exp + "_percent_PSMs");
        for (String exp : experiments)
            colsToAddLater.add(exp + "_peptides");
        for (String exp : experiments)
            colsToAddLater.add(exp + "_percent_also_in_unmodified");
        if (useIntensities == 1) {
            for (String exp : experiments)
                colsToAddLater.add(exp + "_intensity");
        }

        /* Read remaining lines and append to data */
        String cline;
        ArrayList<StringBuffer> linesToAddLater = new ArrayList<>();
        while ((cline = in.readLine()) != null) {
            String[] sp = cline.split("\t", -1);
            StringBuffer sb = new StringBuffer();
            StringBuffer sb2 = new StringBuffer();
            for (String col : colsToAdd)
                sb.append(sp[getColumn(col, curHeaders)] + "\t");
            for (String col : colsToAddLater)
                sb2.append(sp[getColumn(col, curHeaders)] + "\t");
            this.data.add(sb);
            linesToAddLater.add(sb2);
        }

        in.close();

        /* Process *.simrtprofile.txt */

        if (calcIntensities) {
            colsToAdd = new String[]{"similarity", "rt_shift", "int_log2fc"};
        } else {
            colsToAdd = new String[]{"similarity", "rt_shift"};
        }

        in = new BufferedReader(new FileReader(new File(
                PTMShepherd.normFName(dataset + ".simrtprofile.txt"))));

        /* Add headers */
        for (String col : colsToAdd)
            this.headers.add(col);
        /* Read remaining lines and append to data */
        curHeaders = in.readLine().split("\t", -1);
        int lineIndx = 0;
        while ((cline = in.readLine()) != null) {
            String[] sp = cline.split("\t", -1);
            StringBuffer sb = new StringBuffer();
            for (String col : colsToAdd)
                sb.append(sp[getColumn(col, curHeaders)] + "\t");
            this.data.get(lineIndx).append(sb);
            lineIndx++;
        }

        in.close();

        /* Process *.locprofile.txt */

        /* Get cols to add to beginning from *.locprofile */
        colsToAdd = new String[]{"localized_PSMs", "n-term_localization_rate",
                "AA1", "AA1_enrichment_score", "AA1_psm_count",
                "AA2", "AA2_enrichment_score", "AA2_psm_count",
                "AA3", "AA3_enrichment_score", "AA3_psm_count"};

        in = new BufferedReader(new FileReader(new File(
                PTMShepherd.normFName(dataset + ".locprofile.txt"))));

        /* Add headers */
        for (String col : colsToAdd)
            this.headers.add(col);

        /* Read remaining lines and append to data */
        curHeaders = in.readLine().split("\t", -1);
        lineIndx = 0;
        while ((cline = in.readLine()) != null) {
            String[] sp = cline.split("\t", -1);
            StringBuffer sb = new StringBuffer();
            for (String col : colsToAdd)
                sb.append(sp[getColumn(col, curHeaders)] + "\t");
            this.data.get(lineIndx).append(sb);
            lineIndx++;
        }

        in.close();

        /* Process remaining lines from peaksummary.annotated.tsv */

        /* Add remaining headers */
        for (String col : colsToAddLater)
            this.headers.add(col);

        /* Add remaining lines */
        for (int i = 0; i < this.data.size(); i++)
            this.data.get(i).append(linesToAddLater.get(i));

        /* Write output file */
        PrintWriter out = new PrintWriter(new FileWriter(this.fname));

        /* Print headers, remove trailing tab */
        StringBuffer headBuff = new StringBuffer();
        for (String head : this.headers)
            headBuff.append(head + "\t");
        out.println(headBuff.toString().substring(0, headBuff.length() - 1));

        /* Print all new lines */
        for (StringBuffer newLine : this.data)
            out.println(newLine.toString());//.substring(0, newLine.length()-1));

        out.close();
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

    public int getColumn(String head, String[] headers) {
        for(int i = 0; i < headers.length; i++)
            if(headers[i].equals(head))
                return i;
        return -1;
    }

}
