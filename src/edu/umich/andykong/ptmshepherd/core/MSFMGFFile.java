package edu.umich.andykong.ptmshepherd.core;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;
import java.util.concurrent.*;

public class MSFMGFFile {

    private static final int MSFMGFVersion = 20200824;

    public final File f;
    public final String runName;
    public ArrayList<Spectrum> specs;
    public ArrayList<Integer> scanNums;
    public int maxScan;

    private int runStartTime = 0;
    private String mzMLHeader = "";
    private String rawSHA1 = "";
    private String activationMethods = "";
    private String batmassVersion = "";
    private String timsLibVersion = "";

    //Two constructors, the first
    public MSFMGFFile(ExecutorService executorService, int nThread, String filePath, boolean loadScans) throws Exception {
        this(executorService, nThread, new File(filePath), loadScans);
    }

    public MSFMGFFile(ExecutorService executorService, int nThread, File f, boolean loadScans) throws Exception {
        this.f = f;
        this.runName = getBasename(f.getName());
        this.specs = new ArrayList<>();
        this.scanNums = new ArrayList<>();
        readMSFMGF(executorService, nThread, loadScans);
    }

    public static String getBasenameStat(String f) {
        String baseName;
        if (f.contains("_calibrated"))
            baseName = f.substring(0, f.indexOf("_calibrated"));
        else if (f.contains("_uncalibrated"))
            baseName = f.substring(0, f.indexOf("_uncalibrated"));
        else
            baseName = f.substring(0, f.lastIndexOf("."));
        return baseName;
    }

    //Parses MSFragger appended substring and file type out of name
    public String getBasename(String f) {
        String baseName;
        if (f.contains("_calibrated"))
            baseName = f.substring(0, f.indexOf("_calibrated"));
        else if (f.contains("_uncalibrated"))
            baseName = f.substring(0, f.indexOf("_uncalibrated"));
        else
            baseName = f.substring(0, f.lastIndexOf("."));
        return baseName;
    }

    private void readMSFMGF(ExecutorService executorService, int nThread, boolean loadScans) throws Exception {
        loadStructure(); //should load start indices here too //todo
        if (loadScans)
            loadMSFMGFScans(executorService, nThread);
        //todo
    }

    private void loadStructure() throws Exception {
        this.scanNums = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader(this.f));
        String cline;
        boolean readFlag = false;
        int bScan = 0;
        int eScan = 0;
        while ((cline = br.readLine()) != null) {
            cline = cline.trim();
            if (cline.equals("BEGIN IONS")) {
                bScan++;
                readFlag = true;
            } else if (cline.equals("END IONS")) {
                eScan++;
            }
            if (readFlag) {
                if (cline.substring(0, 4).equals("TITLE")) {
                    scanNums.add(getScanNum(cline));
                    readFlag = false;
                }
            }
        }
        this.maxScan = Math.max(bScan, eScan);
    }

    private int getScanNum(String line) {
        String scanStr = line.trim().split("=")[1].trim();
        int scanNum = Integer.parseInt(scanStr.split(".")[1]);
        return scanNum;
    }

    private void loadMSFMGFScans(ExecutorService executorService, int nThread) throws Exception {
        this.scanNums = new ArrayList<>();
        this.specs = new ArrayList<>();
        long t1 = System.currentTimeMillis();

        BufferedReader br = new BufferedReader(new FileReader(this.f));

        String cline;
        ArrayList<String> specBlock = new ArrayList<>();

        //long factor = this.maxScan / nThread;
        //int remainder = (int) (this.maxScan - (factor * nThread));
        final int BLOCKSIZE = 500; //number of scans to be parsed per thread (to cut down on thread creation overhead)
        int nBlocks = this.maxScan / (BLOCKSIZE); //number of jobs submitted to queue
        if (this.maxScan % BLOCKSIZE != 0) //if there are missing scans, add one more block
            nBlocks++;

        ArrayList<String> fileLines = new ArrayList<>();
        int[] scanStarts = new int[this.maxScan];

        int scanC = 0;
        int lineC = 0;
        while ((cline = br.readLine()) != null) {
            fileLines.add(cline);
            if (cline.equals("BEGIN IONS"))
                scanStarts[scanC] = lineC; //0 indexed
            else if (cline.equals("END IONS"))
                scanC++;
            lineC++;
        }

        System.out.println("Time: " + Long.toString(System.currentTimeMillis() - t1));

        ArrayList<Future> futureList = new ArrayList<>(nBlocks);

        for (int i = 0; i < nBlocks; i++) {
            int istart = i * BLOCKSIZE; //start scan 0
            int iend = (i + 1) * BLOCKSIZE; //end scan 500
            int ilinestart = scanStarts[istart]; //start at line corresponding to line 0 (0)
            int ilineend;
            if (iend > scanC) {
                iend = scanC;
                ilineend = fileLines.size() - 1; //last line
            } else
                ilineend = scanStarts[iend + 1] - 1; //end at line corresponding to line before scan 501

            futureList.add(executorService.submit(() -> processMSFMGFSpectraBlock(fileLines, ilinestart, ilineend)));
        }
//
//                ArrayList<String>
//                try {
//                    ArrayList<String> tSpecBlock = Arrays.copyOf(specBlock);
//                    ArrayList<Spectrum> specs = loadMSFMGFSpectra(specBlock, BLOCKSIZE);
//                    addSpecs(specs);
//                } catch (Exception e) {
//                    e.printStackTrace();
//                    System.exit(1);
//                }
        //    }));
        //}

        //System.out.println("Time: " + Long.toString(t1 - System.currentTimeMillis()));
        //System.exit(1);


        //loop through file with single thread to minimize hopping around file
//        while ((cline = br.readLine()) != null) {
//            specBlock.add(cline);
//            if (cline.equals("END IONS")) {
//                nScans++;
//                if (nScans >= BLOCKSIZE) {
//                    System.out.println(specBlock); //exists here
//                    futureList.add(executorService.submit(() -> {
//                        try {
//                            ArrayList<String> tSpecBlock = Arrays.copyOf(specBlock);
//                            ArrayList<Spectrum> specs = loadMSFMGFSpectra(specBlock, BLOCKSIZE);
//                            addSpecs(specs);
//                        } catch (Exception e) {
//                            e.printStackTrace();
//                            System.exit(1);
//                        }
//                    }));
//                    //specBlock.clear();
//                    nScans = 0;
//                }
//            }
//        }
//        //submit remainder to executor... I should probably make a separate function for this... //todo
//        futureList.add(executorService.submit(() -> {
//            try {
//                ArrayList<Spectrum> specs = loadMSFMGFSpectra(specBlock, BLOCKSIZE);
//                addSpecs(specs);
//            } catch (Exception e) {
//                e.printStackTrace();
//                System.exit(1);
//            }
//        }));
        System.out.println("Time: " + Long.toString(System.currentTimeMillis() - t1));
        for (Future future : futureList) { //checks to make sure all threads are done
            future.get();
        }
        System.out.println("Time: " + Long.toString(System.currentTimeMillis() - t1));
        //todo
    }

    //Processes a given block of the mgf file and adds it to the class variables
    //Note: the lines var is a shared resource, but I don't *think* sharing it between threads will slow it down
    //Note: will throw an non-fatal error if a scan is bad, but won't add it to the array so it will be hidden
    // from downstream processes. Can be made fatal under the catch block, but if that is an issue then it should be
    // made fatal earlier in the parsing.
    //TODO
    private void processMSFMGFSpectraBlock(ArrayList<String> lines, int istart, int iend) {
        ArrayList<Spectrum> spectra = new ArrayList<>();
        ArrayList<String> cspecLines = new ArrayList<>();
        String line;

        for (int i = istart; i < iend; i++) {
            line = lines.get(i);
            cspecLines.add(line);
            if (line.equals("END IONS")) {
                try {
                    spectra.add(loadMSFMGFSpectrum(cspecLines));
                } catch (Exception e) {
                    System.out.printf("Malformed scan: %s\n", cspecLines.get(1));
                    e.printStackTrace();
                    //System.exit(1);
                }
                cspecLines.clear();
            }
        }
        addSpecs(spectra);
    }

    //loads a BLOCKSIZE number of spectra to the class specs and scanNums vars
    //Note: doing this is a block of BLOCKSIZE reduces thread creation overhead and mitigates
    // threads competing for locks on the shared class variables
    private synchronized void addSpecs(ArrayList<Spectrum> specs) {
        this.specs.addAll(specs);
        for (Spectrum spec : specs)
            this.scanNums.add(spec.scanNum);
    }

    private int getScanNumFromTitle(String title) {
        int scanNum;
        String[] bits = title.trim().split("\\.");
        if (bits.length == 4)
            scanNum = Integer.parseInt(bits[2]);
        else
            scanNum = Integer.parseInt(bits[1]);
        return scanNum;
    }

    private String getScanNameFromTitle(int scanNum) {
        String scanname = String.format("%s.%d.%d", this.runName, scanNum, scanNum);
        return scanname;
    }

    //loads a single MGF given a series of lines from MGF file
    //Note: ignores any BEGIN and END IONS lines, so it doesn't handle exceptions where two spectra
    // are stuck together because one of these lines has been deleted from the mgf //todo
    private Spectrum loadMSFMGFSpectrum(ArrayList<String> block) throws Exception{
        String scanName = "";
        int scanNum = -1;
        int charge = 0;
        int msLevel = 2;
        double precursorMZ = -1;
        double retTime = -1;
        ArrayList<Float> peakmz = new ArrayList<>();
        ArrayList<Float> peakint = new ArrayList<>();
        String[] cinfo;
        String att, val;

        //this can be made more flexible using a lookup array with n-keywords for generic MGFs
        //where n == count(intersection of all lookup terms)
        try {
            for (String line : block) {
                if (line.equals("BEGIN IONS") || line.equals("END IONS"))
                        continue;
                if (line.contains("=")) {
                    cinfo = line.trim().split("=");
                    att = cinfo[0];
                    val = cinfo[1];
                    if (att.equals("TITLE")) {
                        scanNum = getScanNumFromTitle(val);
                        scanName = getScanNameFromTitle(scanNum);
                        //charge = Integer.parseInt(Character.toString((val.split("\\.")[val.length()-1]).charAt(0)));
                    } else if (att.equals("CHARGE")) {
                        charge = Integer.parseInt(val.substring(0,val.indexOf("+")));
                    } else if (att.equals("RTINSECONDS")) {
                        retTime = Double.parseDouble(val);
                    } else if (att.equals("PEPMASS")) {
                        precursorMZ = Double.parseDouble(val.split(" ")[0]);
                    }
                } else {
                    cinfo = line.trim().split(" ");
                    peakmz.add(Float.parseFloat(cinfo[0]));
                    peakint.add(Float.parseFloat(cinfo[0]));
                }
            }
        } catch (Exception e) {
            System.out.printf("Malformed scan: %s\n", block.get(1));
            e.printStackTrace();
            System.exit(1);
        }

        float[] peakmzArr = new float[peakmz.size()];
        for(int i = 0; i < peakmz.size(); i++)
            peakmzArr[i] = peakmz.get(i);
        float[] peakintArr = new float[peakint.size()];
        for(int i = 0; i < peakmz.size(); i++)
            peakintArr[i] = peakint.get(i);

        return new Spectrum(scanName, scanNum, charge, msLevel, precursorMZ, retTime, peakmzArr, peakintArr);
    }
}
