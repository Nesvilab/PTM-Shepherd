package edu.umich.andykong.ptmshepherd.peakpicker;

import edu.umich.andykong.ptmshepherd.peakpicker.PeakAnnotator;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class ModSummary {
    HashMap<String, HashMap<String, Double>> psmCountsNorm;
    HashMap<String, HashMap<String, Double>> psmCounts;
    Set<String> datasets;
    Set<String> mods;
    LinkedHashMap<String, Double> modsMap;



    public ModSummary(File inPeaks, Set<String> dsets) throws Exception{
        datasets = dsets;

        ArrayList<String[]> inFile = new ArrayList<>();
        String cline;
        BufferedReader in = new BufferedReader(new FileReader(inPeaks));
        while((cline = in.readLine())!= null) {
            inFile.add(cline.split("\t"));
        }
        in.close();

        // get headers
        String [] headers = inFile.get(0);
        HashMap<String, Integer> headis = new HashMap<>();
        for(int i = 0; i < headers.length; i++) {
            headis.put(headers[i], i);
        }

        // get modifications and masses
        int modi = 0;
        mods = new HashSet<>();
        modsMap = new LinkedHashMap<>();

        //get modification names and map apexes
        for(int i = 1; i <= 2; i++){
            modi = headis.get("Potential Modification "+i);
            for(int j = 1; j < inFile.size(); j++) {
                //System.out.println(inFile.get(j));
                if(inFile.get(j).length > modi) {
                    String mod = inFile.get(j)[modi];
                    mods.add(mod);
                    for(int k = 0; k < PeakAnnotator.mods.size(); k++){
                        if(PeakAnnotator.mods.get(k).equals(mod)){
                            modsMap.put(mod, PeakAnnotator.mod_diffs.get(k));
                        }
                    }
                }
            }
        }
        modsMap.put("None", 0.0);

        // initialize empty hashmaps
        int mod1i = headis.get("Potential Modification 1");
        psmCounts = new HashMap<String, HashMap<String, Double>>();
        psmCountsNorm = new HashMap<String, HashMap<String, Double>>();
        for(String ds : datasets){
            psmCounts.put(ds, new HashMap<>());
            psmCountsNorm.put(ds, new HashMap<>());
            for (String mod : mods) {
                psmCounts.get(ds).put(mod, 0.0);
                psmCountsNorm.get(ds).put(mod, 0.0);
            }
            for(int j = 1; j < inFile.size(); j++){
                String line [] = inFile.get(j);
                if(line.length == mod1i){
                    psmCounts.get(ds).put("None", Double.parseDouble(line[headis.get(ds + " (PSMs)")]));
                    psmCountsNorm.get(ds).put("None", Double.parseDouble(line[headis.get(ds + " (PSMs/million)")]));
                }
            }
        }
        //populate hashmap
        for(int i = 1; i < inFile.size(); i++){
            for(int j = modi; j > modi-2; j--){
                String [] line = inFile.get(i);
                if(line.length > j){
                    for(String ds : datasets) {
                        int col = headis.get(ds + " (PSMs)");
                        double prevCount = psmCounts.get(ds).get(line[j]);
                        psmCounts.get(ds).put(line[j], Double.parseDouble(line[col]) + prevCount);
                        col = headis.get(ds + " (PSMs/million)");
                        prevCount = psmCountsNorm.get(ds).get(line[j]);
                        psmCountsNorm.get(ds).put(line[j], Double.parseDouble(line[col]) + prevCount);
                    }
                }
            }
        }
    }

    public void toFile(File modOut) throws Exception{
        //sort modifications by psm count for table
        LinkedHashMap<String, Integer> modsCount = new LinkedHashMap<>();
        mods.add("None"); //allow mods to be searched for
        for (String mod : mods) { //sum PSMs across DSs
            int nPsms = 0;
            for (String ds : datasets){
                nPsms += psmCounts.get(ds).get(mod);
            }
            modsCount.put(mod, nPsms);
        }
        List<Map.Entry<String, Integer>> entries = //sorts new LinkedHashMap based on summed values
                new ArrayList<Map.Entry<String, Integer>>(modsCount.entrySet());
        Collections.sort(entries, new Comparator<Map.Entry<String, Integer>>() {
            public int compare(Map.Entry<String, Integer> a, Map.Entry<String, Integer> b){
                return -a.getValue().compareTo(b.getValue());
            }
        });
        Map<String, Integer> sortedMap = new LinkedHashMap<String, Integer>();
        for (Map.Entry<String, Integer> entry : entries) {
            sortedMap.put(entry.getKey(), entry.getValue()); //psms sorted by sum
        }
        // write to file
        PrintWriter out = new PrintWriter(new FileWriter(modOut));
        out.print("Modification\tTheoretical Mass Shift");
        for(String ds : datasets){
            out.printf("\t%s (PSMs)\t%s (PSMs/million)", ds, ds);
        }
        out.println();
        // output sorted by sum of spectral counts
        for(String mod : sortedMap.keySet()){
            out.printf("%s\t%s", mod, modsMap.get(mod));
            for(String ds : datasets){
                out.printf("\t%.00f\t%.04f",psmCounts.get(ds).get(mod),psmCountsNorm.get(ds).get(mod));
            }
            out.println();
        }
        out.flush();
        out.close();
    }


}
