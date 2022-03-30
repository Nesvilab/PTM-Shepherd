package edu.umich.andykong.ptmshepherd.cleaner;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class CombinedExperimentsSummary {
    public ArrayList<String> headers;
    public ArrayList<StringBuffer> data;
    public File fname;
    public HashMap <String, ArrayList<String>> expToData;

    public CombinedExperimentsSummary(String fn) throws Exception {
        this(new File(fn));
    }

    public CombinedExperimentsSummary(File fn) throws Exception {
        this.fname = fn;
        this.headers = new ArrayList<>();
        this.data = new ArrayList<>();
    }

    /* Read global profile table to initialize data */
    public void initializeExperimentSummary(String fn, int useIntensity) throws IOException {
        /* Universal columns to be added */
        ArrayList<String> colsToAdd = new ArrayList<>(Arrays.asList("peak_apex", "peak_lower", "peak_upper",
                "PSMs", "percent_also_in_unmodified",
                "mapped_mass_1", "mapped_mass_2",
                "localized_PSMs", "n-term_localization_rate",
                "AA1", "AA1_enrichment_score", "AA1_psm_count",
                "AA2", "AA2_enrichment_score", "AA2_psm_count",
                "AA3", "AA3_enrichment_score", "AA3_psm_count",
                "peak_signal"));

        BufferedReader in = new BufferedReader(new FileReader(new File(fn)));

        /* Add exp-level values to file file */
        String[] curHeaders  = in.readLine().split("\t", -1);
        for (int i = 0; i < curHeaders.length; i++) {
            if ((curHeaders[i].contains("_PSMs") && !curHeaders[i].contains("_percent_PSMs"))||
                    curHeaders[i].contains("_percent_PSMs") ||
                    curHeaders[i].contains("_peptides") ||
                    curHeaders[i].contains("_percent_also_in_unmodified") ||
                    (curHeaders[i].contains("_intensity") && useIntensity == 1));
                colsToAdd.add(curHeaders[i]);
        }

        /* Add all new headers to new file */
        for (String col : colsToAdd)
            this.headers.add(col);

        /* Read remaining lines and add new lines to data */
        String cline;
        while((cline = in.readLine()) != null) {
            String[] sp = cline.split("\t", -1);
            StringBuffer sb = new StringBuffer();
            for (String col : colsToAdd)
                sb.append(sp[getColumn(col, curHeaders)] + "\t");
            this.data.add(sb);
        }

        in.close();
    }

    /* Add single experiment summary to the combined file */
    public void addExperimentSummary(String fn, String exp) throws IOException {
        /* These are the columns that will be added */
        String[] colsToAdd = new String[]{"localized_PSMs", "n-term_localization_rate",
                "AA1", "AA1_enrichment_score", "AA1_psm_count",
                "AA2", "AA2_enrichment_score", "AA2_psm_count",
                "AA3", "AA3_enrichment_score", "AA3_psm_count",
                "peak_signal"};

        BufferedReader in = new BufferedReader(new FileReader(new File(fn)));

        /* Add new headers to new file and tag them with experiment name */
        for (String col : colsToAdd)
            this.headers.add(exp + "_" + col);

        /* Read remaining lines and append to data */
        String[] curHeaders  = in.readLine().split("\t", -1);
        String cline;
        int lineIndx = 0;
        while((cline = in.readLine()) != null) {
            String[] sp = cline.split("\t", -1);
            StringBuffer sb = new StringBuffer();
            for (String col : colsToAdd)
                sb.append(sp[getColumn(col, curHeaders)] + "\t");
            this.data.get(lineIndx).append(sb);
            lineIndx++;
        }

        in.close();
    }

    /* Read localization profile table to initialize data */
    public void addLocalizationSummary(String fn, String ds) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(new File(fn)));

        /* Add exp-level values to file file */
        ArrayList<String> colsToAdd = new ArrayList<>();
        String[] curHeaders  = in.readLine().split("\t", -1);
        for (int i = 0; i < curHeaders.length; i++) {
            if (ds.equals("combined")) {
                if (curHeaders[i].endsWith("_enrichment"))
                    colsToAdd.add(curHeaders[i]);
            } else {
                if (curHeaders[i].contains("_enrichment"))
                    colsToAdd.add(curHeaders[i]);
            }
        }

        /* Add all new headers to new file */
        for (String col : colsToAdd)
            this.headers.add(ds+"_"+col);

        /* Read remaining lines and add new lines to data */
        String cline;
        int lineIndx = 0;
        while((cline = in.readLine()) != null) {
            String[] sp = cline.split("\t", -1);
            StringBuffer sb = new StringBuffer();
            for (String col : colsToAdd)
                sb.append(sp[getColumn(col, curHeaders)] + "\t");
            this.data.get(lineIndx).append(sb);
            lineIndx++;
        }

        in.close();
    }

    /* Read simrt profile table to initialize data */
    public void addSimilarityRTSummary(String fn, String ds, boolean calcIntensity) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(new File(fn)));

        String[] colsToAdd;
        if (calcIntensity) {
            colsToAdd = new String[]{"similarity", "similarity_(variance)", "rt_shift",
                    "rt_shift_(variance)", "int_log2fc", "int_log2fc_(variance)"};
        } else {
            colsToAdd = new String[]{"similarity", "similarity_(variance)", "rt_shift",
                    "rt_shift_(variance)"};
        }

        /* Add all new headers to new file */
        String[] curHeaders  = in.readLine().split("\t", -1);
        for (String col : colsToAdd) {
            if (ds.equals("combined"))
                this.headers.add(col);
            else
                this.headers.add(ds + "_" + col);
        }

        /* Read remaining lines and add new lines to data */
        String cline;
        int lineIndx = 0;
        while((cline = in.readLine()) != null) {
            String[] sp = cline.split("\t", -1);
            StringBuffer sb = new StringBuffer();
            for (String col : colsToAdd)
                sb.append(sp[getColumn(col, curHeaders)] + "\t");
            this.data.get(lineIndx).append(sb);
            lineIndx++;
        }

        in.close();
    }

    /* Print the file */
    public void printFile() throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(this.fname));

        /* Print headers, remove trailing tab */
        StringBuffer headBuff = new StringBuffer();
        for (String head : this.headers)
            headBuff.append(head + "\t");
        out.println(headBuff.toString().substring(0, headBuff.length()-1));

        /* Print all new lines, remove trailing tab */
        for (StringBuffer newLine : this.data)
            out.println(newLine.toString().substring(0, newLine.length()-1));

        out.close();

    }

    public int getColumn(String head, String[] headers) {
        for(int i = 0; i < headers.length; i++)
            if(headers[i].equals(head))
                return i;
        return -1;
    }
}
