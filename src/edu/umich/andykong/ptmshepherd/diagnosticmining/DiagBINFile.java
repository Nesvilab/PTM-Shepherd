package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.core.Spectrum;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.zip.CRC32;

public class DiagBINFile {

    private static final int DiagBINFileVersion = 20210303;

    public final File f;
    public final String runName;
    public final ArrayList<DiagnosticRecord> spectra;

    private int indexWidth;
    private int mzScale;
    private byte[] diagBinIndex;
    private int nIonTypes;
    private ArrayList<Character> ionTypes;

    public DiagBINFile(ExecutorService executorService, int nThread, String filePath) throws Exception {
        this(executorService, nThread, new File(filePath));
    }

    public DiagBINFile(ExecutorService executorService, int nThread, File f) throws Exception {
        this.f = f;
        this.runName = getBasename(f.getName());
        this.spectra = new ArrayList<>();
        readDiagBin(executorService, nThread);
    }

    public DiagBINFile(ArrayList<DiagnosticRecord> diagnosticRecords, String filePath, String ionTypes) throws Exception {
        this(diagnosticRecords, new File(filePath), ionTypes);
    }

    public DiagBINFile(ArrayList<DiagnosticRecord> diagnosticRecords, File f, String ionTypes) throws Exception {
        this.f = f;
        this.runName = getBasename(f.getName());
        this.ionTypes = new ArrayList<>();
        for (int i = 0; i < ionTypes.length(); i++)
            this.ionTypes.add(ionTypes.toLowerCase().charAt(i));
        this.spectra = new ArrayList<>(diagnosticRecords);
        /* Ensure proper order of diagnostic record and ion types */
        Collections.sort(this.spectra);
        Collections.sort(this.ionTypes);
    }

    public void writeDiagBinFile() throws Exception {
        /*
         * Binary values are stored in network order (default for Java)
         * All values are signed
         * DiagBIN file contains 5 major sections
         * 1. DiagBIN header (20 bytes)
         * 2. spectra index  ((maxScanNumber + 1) * INDEXWidth) bytes)
         * 3. spectra peak data
         * 4. DiagBIN textual metadata
         * 5. DiagBIN footer (12 bytes)
         */

        /* Get max scan number */
        int maxScan = 0;
        for(DiagnosticRecord dr : this.spectra) {
            if (maxScan < dr.scanNum)
                maxScan = dr.scanNum;
        }

        /* Format DiagBIN header */
        this.indexWidth = 80;
        this.mzScale = 200000;
        /* The DiagBIN header consists of 28 bytes */
        byte[] diagBinHead = new byte[52];
        /* Bytes 0-4   : 'DGBN' - to indicate that this is a DGBN file
        *       4-8   : 0 - 4 byte integer indicating the current DGBN version number
        *       8-12 : 4 byte integer containing the scaling factor (MZscale) to perform integer <-> doubles for M/Z values (currently 200000)
        *       12-16 : 4 byte integer indicating the width (in bytes) of the spectra index (currently 80)
        *       16-20 : 4 byte integer containing the maximum scan ID as numbered by vendor library
        *       20-52 : 4 byte integers acting as bools to store ion types present in file
        */

        ByteBuffer.wrap(diagBinHead).asIntBuffer().put(0, 1297760847); //todo change to DGBN, just added 1 to MZBN
        ByteBuffer.wrap(diagBinHead).asIntBuffer().put(1, this.DiagBINFileVersion);
        ByteBuffer.wrap(diagBinHead).asIntBuffer().put(2, this.mzScale);
        ByteBuffer.wrap(diagBinHead).asIntBuffer().put(3, this.indexWidth);
        ByteBuffer.wrap(diagBinHead).asIntBuffer().put(4, maxScan);

        /* Convert iontypes to be printed to file */
        boolean[] boolsIonTypes = ionTypesToBoolArray(this.ionTypes);
        for (int i = 0; i < boolsIonTypes.length; i++)
            ByteBuffer.wrap(diagBinHead).asIntBuffer().put(i+5, boolsIonTypes[i]?1:0);

        // The scaled int data type is a m/z value stored in the mzBIN format as an integer
        // conversion between a floating point m/z value to an integer is defined as
        // <scaled int m/z> = round(float m/z / MZscale)
        // conversion in the unpacking process is simply a multiplication by MZscale

        /* Write DiagBIN header */
        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f), 1 << 24));
        CRC32 crc = new CRC32();
        dos.write(diagBinHead);
        crc.update(diagBinHead);
        //dos.write(this.diagBinIndex);
        //crc.update(this.diagBinIndex);
        dos.close();
    }

    public void readDiagBin(ExecutorService executorService, int nThread) throws Exception {
        readDiagBinStructure();
        loadDiagBinScans(executorService, nThread);
    }

    private void loadDiagBinScans(ExecutorService executorService, int nThread) throws Exception {
        ArrayList<Integer> scanNums = new ArrayList<>();
        int maxScans = this.diagBinIndex.length / this.indexWidth - 1;
        for (int i = 0; i <= maxScans; i++) {
            if (ByteBuffer.wrap(this.diagBinIndex).asIntBuffer().get((this.indexWidth / 4) * i + 2) > 0)
                scanNums.add(i);
        }
        loadDiagBinScans(executorService, nThread, scanNums);
    }

    private void loadDiagBinScans(ExecutorService executorService, int nThread, ArrayList<Integer> scanNums) throws Exception {
        FileChannel fc = FileChannel.open(f.toPath(), StandardOpenOption.READ);
        int factor = Math.max(4, 64 / nThread);
        List<Future> futureList = new ArrayList<>(factor * nThread);
        for (int i = 0; i < factor * nThread; i++) {
            int start = (int) (scanNums.size() * (long) i) / (factor * nThread);
            int end = (int) (scanNums.size() * (long) (i + 1)) / (factor * nThread);
            futureList.add(executorService.submit(() -> {
                ArrayList<DiagnosticRecord> spectra = new ArrayList<>();
                try {
                    for (int j = start; j < end; j++) {
                        //spectra.add(loadDiagBinSpectrum(fc, scanNums.get(j))); todo
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
                //addScans(scans); todo
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
        fc.close();
        Collections.sort(this.spectra);
    }

    private void loadDiagBinSpectrum(FileChannel fc, int scanNum) throws Exception {
        int currPos = (this.indexWidth / 4) * scanNum;

        //int nPeakArrays =  ByteBuffer.wrap(this.diagBinIndex).asIntBuffer(); todo
    }

    private void readDiagBinStructure() throws IOException {
        RandomAccessFile raf = new RandomAccessFile(f, "r");

        byte[] MZDGhead = new byte[52];
        raf.read(MZDGhead);

        /* Collect metadata from head of file */
        int diagBinId = ByteBuffer.wrap(MZDGhead).asIntBuffer().get(0);
        int diagBinVersion = ByteBuffer.wrap(MZDGhead).asIntBuffer().get(1);
        this.mzScale = ByteBuffer.wrap(MZDGhead).asIntBuffer().get(2);
        this.indexWidth = ByteBuffer.wrap(MZDGhead).asIntBuffer().get(3);
        int maxScan = ByteBuffer.wrap(MZDGhead).asIntBuffer().get(4);

        boolean[] boolIonTypes = new boolean[8];
        for (int i = 0; i < boolIonTypes.length; i++) {
            if (ByteBuffer.wrap(MZDGhead).asIntBuffer().get(i + 5) == 1)
                boolIonTypes[i] = true;
            else
                boolIonTypes[i] = false;
        }
        setIonTypes(boolIonTypes);

        /* Reads in the spectrum index for fast spec access */
        this.diagBinIndex = new byte[(maxScan + 1) * this.indexWidth];
        raf.read(this.diagBinIndex);

        raf.close();
    }

    private String getBasename(String f) {
        String baseName = f.substring(0,f.lastIndexOf("."));
        return baseName;
    }

    private static String getBasenameStat(String f) {
        String baseName = f.substring(0,f.lastIndexOf("."));
        return baseName;
    }

    /* Converts all 8 possible ion types into a 4 byte integer */
    private boolean[] ionTypesToBoolArray(ArrayList<Character> ionTypes) {
        boolean[] boolIonTypes = new boolean[8];

        /* Fill boolarray with trues for iontypes included in file */
        boolIonTypes[0] = true; /* Always true for oxonium ions */
        boolIonTypes[1] = true; /* Always true for capY ions */
        for (Character c : ionTypes) {
            if (c - 'a' <  13)
                boolIonTypes[Character.toLowerCase(c) - 'a' + 2] = true;
            else
                boolIonTypes[Character.toLowerCase(c) - 'x' + 5] = true;
        }

        return boolIonTypes;

        /* Convert boolean array to int
        int n = 0;
        for (boolean b : boolIonTypes)
            n = (n << 1) | (b ? 1 : 0);

        return n;
        */

    }

    /* Converts int to bool array containing iontypes */ //todo
    private void setIonTypes(boolean[] boolIonTypes) {
        this.nIonTypes = 2;
        this.ionTypes = new ArrayList<>();
        for (int i = 2; i < 5; i++) {
            if (boolIonTypes[i])
                this.ionTypes.add((char) ('a' + (i - 2)));
        }
        for (int i = 5; i < 8; i++) {
            if (boolIonTypes[i])
                this.ionTypes.add((char) ('x' + (i - 5)));
        }
    }
}
