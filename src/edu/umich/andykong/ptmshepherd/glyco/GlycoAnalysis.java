package edu.umich.andykong.ptmshepherd.glyco;

import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.util.HashMap;

public class GlycoAnalysis {
    String dsName;
    File glycoFile;
    MXMLReader mr;
    HashMap<String, MXMLReader> multiMr;
    double ppmTol, condRatio, peakTol;
    int condPeaks;
    int specCol, pepCol, modpepCol, chargecol, deltaCol, rtCol, intCol;

    public GlycoAnalysis(String dsName) {
        this.dsName = dsName;
        glycoFile = new File(dsName+".rawglyco");
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
