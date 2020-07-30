package edu.umich.andykong.ptmshepherd.glyco;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;

public class GlycoAnalysis {
    String dsName;
    File glycoFile;
    MXMLReader mr;
    HashMap<String, MXMLReader> multiMr;
    double ppmTol, condRatio, peakTol;
    int condPeaks;
    int specCol, pepCol, modpepCol, chargecol, deltaCol, rtCol, intCol, pmassCol;

    public GlycoAnalysis(String dsName) {
        this.dsName = dsName;
        glycoFile = new File(dsName+".rawglyco");
    }

    public void glycoPSMs(PSMFile pf, HashMap<String, File> mzMappings) throws Exception {
        //open up output file
        HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
        PrintWriter out = new PrintWriter(new FileWriter(glycoFile, true));
        //write header
        out.printf("%s\t%s\t%s\t%s\t%s\n", "Spectrum", "Peptide", "Mod_Peptide", "Mass_Shift", "Pep_Mass");
        //get necessary col indices
        specCol = pf.getColumn("Spectrum");
        pepCol = pf.getColumn("Peptide");
        modpepCol = pf.getColumn("Modified Peptide");
        //chargeCol = pf.getColumn("Charge");
        deltaCol = pf.dMassCol;
        rtCol = pf.getColumn("Retention");
        intCol = pf.getColumn("Intensity");

        //map PSMs to file
        for (int i = 0; i < pf.data.size(); i++) {
            String[] sp = pf.data.get(i).split("\t");
            String bn = sp[specCol].substring(0, sp[specCol].indexOf(".")); //fraction
            if (!mappings.containsKey(bn))
                mappings.put(bn, new ArrayList<>());
        }
        for (int i = 0; i < pf.data.size(); i++) {
            for (String fraction : mappings.keySet())
                mappings.get(fraction).add(i);
        }

        for (String cf : mappings.keySet()) { //for file in relevant spectral files
            long t1 = System.currentTimeMillis();
            mr = new MXMLReader(mzMappings.get(cf), Integer.parseInt(PTMShepherd.getParam("threads")));
            mr.readFully();
            ArrayList<Integer> clines = mappings.get(cf); //lines corr to curr spec file
            for (int i = 0; i < clines.size(); i++) {//for relevant line in curr spec file
                try {
                    //out.println(processLine());
                } catch (Exception e) {
                    e.printStackTrace();
                    System.out.println("Error in: " + pf.data.get(clines.get(i)));
                }
            }
        }
    }

    public void processLine(String line) {
        StringBuffer sb = new StringBuffer();
        String [] sp = line.split("\\t");
        String seq = sp[pepCol];
        float dmass = Float.parseFloat(sp[deltaCol]);

        String specName = sp[specCol];
        //String [] smods = sp[modCol].split(",");

        //sb.append(String.format("%s\t%s\t%s\t%.4f", specName,seq,sp[modCol],dmass));

        //StringBuffer sb = new StringBuffer();
        String [] s = line.split("\\t");
        //Spectrum spec = mr.getSpectrum(reNormName(specName))
        //findIonMasses();
        //localizeRemainderFragments();
    }

    public void findIonMasses(float pepMass) {

    }

    public String reNormName(String s) {
        String[] sp = s.split("\\.");
        int sn = Integer.parseInt(sp[1]);
        //with charge state
        //return String.format("%s.%d.%d.%s",sp[0],sn,sn,sp[3]);
        //without charge state
        return String.format("%s.%d.%d", sp[0], sn, sn);
    }

    public boolean isComplete() throws Exception {
        if(glycoFile.exists()) {
            RandomAccessFile raf = new RandomAccessFile(glycoFile, "r");
            raf.seek(Math.max(0, glycoFile.length() - 20));
            String cline;
            while((cline = raf.readLine())!=null)
                if(cline.equals("COMPLETE")) {
                    raf.close();
                    return true;
                }
            raf.close();
            glycoFile.delete();
        }
        return false;
    }

    public void complete() throws Exception {
        PrintWriter out = new PrintWriter(new FileWriter(glycoFile,true));
        out.println("COMPLETE");
        out.close();
    }



}
