package edu.umich.andykong.ptmshepherd.diagnosticmining;

import edu.umich.andykong.ptmshepherd.core.Spectrum;
import io.grpc.netty.shaded.io.netty.buffer.ByteBuf;

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

    public DiagBINFile(ExecutorService executorService, int nThread, String filePath, boolean loadScans) throws Exception {
        this(executorService, nThread, new File(filePath), loadScans);
    }

    public DiagBINFile(ExecutorService executorService, int nThread, File f, boolean loadScans) throws Exception {
        this.f = f;
        this.runName = getBasename(f.getName());
        this.spectra = new ArrayList<>();
        readDiagBin(executorService, nThread, loadScans);
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
         * 1. DiagBIN header (52 bytes)
         * 2. spectra index  ((maxScanNumber + 1) * indexWidth) bytes)
         * 3. spectra peak data
         * 4. DiagBIN footer (12 bytes)
         */

        /* Get max scan number */
        int maxScan = 0;
        for(DiagnosticRecord dr : this.spectra) {
            if (maxScan < dr.scanNum)
                maxScan = dr.scanNum;
        }

        this.indexWidth = 16;
        this.mzScale = 200000; // Unused, but will be useful for higher res data in near future

        /* The DiagBIN header consists of 52 bytes */
        byte[] diagBinHead = new byte[52];
        /* Bytes 0-4    :   'DGBN' - to indicate that this is a DGBN file
        *       4-8     :   0 - 4 byte integer indicating the current DGBN version number
        *       8-12    :   4 byte integer containing the scaling factor (MZscale) to perform integer <-> doubles for M/Z values (currently 200000)
        *       12-16   :   4 byte integer indicating the width (in bytes) of the spectra index (currently 80)
        *       16-20   :   4 byte integer containing the maximum scan ID as numbered by vendor library
        *       20-52   :   4 byte integers acting as bools to store ion types present in file
        */

        ByteBuffer.wrap(diagBinHead).asIntBuffer().put(0, 1145520718);
        ByteBuffer.wrap(diagBinHead).asIntBuffer().put(1, this.DiagBINFileVersion);
        ByteBuffer.wrap(diagBinHead).asIntBuffer().put(2, this.mzScale);
        ByteBuffer.wrap(diagBinHead).asIntBuffer().put(3, this.indexWidth);
        ByteBuffer.wrap(diagBinHead).asIntBuffer().put(4, maxScan);

        /* Convert iontypes to be printed to file */
        boolean[] boolsIonTypes = ionTypesToBoolArray(this.ionTypes);
        for (int i = 0; i < boolsIonTypes.length; i++)
            ByteBuffer.wrap(diagBinHead).asIntBuffer().put(i+5, boolsIonTypes[i]?1:0);

        /* Calculate DiagBIN index, index points to beginning of each spec's info */
        /* Bytes 0-8    :   8 byte long indicating where in the file the spectrum's info starts
         *       8-12   :   4 byte float indicating delta mass
         *       12-16  :   4 byte int indicating line width
         */
        this.diagBinIndex = new byte[(maxScan + 1) * this.indexWidth];

        long offset = diagBinHead.length + this.diagBinIndex.length;
        for (DiagnosticRecord dr : this.spectra) {
            /* Write metadata to index */
            int currIndexPos = (this.indexWidth / 4) * (dr.scanNum);
            ByteBuffer.wrap(this.diagBinIndex).asLongBuffer().put(currIndexPos / 2, offset);
            ByteBuffer.wrap(this.diagBinIndex).asFloatBuffer().put(currIndexPos + 2, dr.dmass);
            ByteBuffer.wrap(this.diagBinIndex).asIntBuffer().put(currIndexPos + 3, calculateDiagnosticRecordLength(dr));
            /* Add length of vals to offset to calc next spec's starting pos */
            int length = calculateDiagnosticRecordLength(dr);
            offset += length;
        }

        /* Write DiagBIN header and index*/
        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f), 1 << 24));
        CRC32 crc = new CRC32();
        dos.write(diagBinHead);
        crc.update(diagBinHead);
        dos.write(this.diagBinIndex);
        crc.update(this.diagBinIndex);

        /* Write diagnosticRecords information */
        for (DiagnosticRecord dr : this.spectra) {
            int length = calculateDiagnosticRecordLength(dr);
            byte[] cdata = new byte[length];
            int cpos = 0; // 32 bit / 4 byte base

            /* Write pep seq info */
            ByteBuffer.wrap(cdata).asIntBuffer().put(cpos, dr.pepSeq.length());
            cpos++;
            for (int i = 0; i < dr.pepSeq.length(); i++)
                ByteBuffer.wrap(cdata).asCharBuffer().put(cpos * 2 + i, dr.pepSeq.charAt(i));

            if (dr.pepSeq.length() % 2 != 0)
                cpos += (int) (dr.pepSeq.length() / 2) + 1;
            else
                cpos += dr.pepSeq.length() / 2;

            /* Write charge info */
            ByteBuffer.wrap(cdata).asIntBuffer().put(cpos, dr.charge);
            cpos++;

            /* Write modification info */
            ByteBuffer.wrap(cdata).asIntBuffer().put(cpos, dr.modifications.size());
            cpos++;
            for (Integer pos : dr.modifications.keySet()) {
                ByteBuffer.wrap(cdata).asIntBuffer().put(cpos, pos);
                cpos++;
                ByteBuffer.wrap(cdata).asFloatBuffer().put(cpos, dr.modifications.get(pos));
                cpos++;
            }

            /* Write immonium peaks */
            ByteBuffer.wrap(cdata).asIntBuffer().put(cpos, dr.immoniumPeaks.length);
            cpos++;
            for (int i = 0; i < dr.immoniumPeaks.length; i++) {
                ByteBuffer.wrap(cdata).asFloatBuffer().put(cpos + i * 2, dr.immoniumPeaks[i][0]);
                ByteBuffer.wrap(cdata).asFloatBuffer().put(cpos + i * 2 + 1, dr.immoniumPeaks[i][1]);
            }
            cpos += 2 * dr.immoniumPeaks.length;

            /* Write capY peaks */
            ByteBuffer.wrap(cdata).asIntBuffer().put(cpos, dr.capYPeaks.length);
            cpos++;
            for (int i = 0; i < dr.capYPeaks.length; i++) {
                ByteBuffer.wrap(cdata).asFloatBuffer().put(cpos + i * 2, dr.capYPeaks[i][0]);
                ByteBuffer.wrap(cdata).asFloatBuffer().put(cpos + i * 2 + 1, dr.capYPeaks[i][1]);
            }
            cpos += 2 * dr.capYPeaks.length;

            /* Write squiggle peaks */
            for (Character it : this.ionTypes) {
                ByteBuffer.wrap(cdata).asIntBuffer().put(cpos, dr.squigglePeaks.get(it).length);
                cpos++;
                for (int i = 0; i < dr.squigglePeaks.get(it).length; i++) {
                    ByteBuffer.wrap(cdata).asFloatBuffer().put(cpos + i * 2, dr.squigglePeaks.get(it)[i][0]);
                    ByteBuffer.wrap(cdata).asFloatBuffer().put(cpos + i * 2 + 1, dr.squigglePeaks.get(it)[i][1]);
                }
                cpos += 2 * dr.squigglePeaks.get(it).length;
            }
            /* Write spec to file */
            dos.write(cdata);
            crc.update(cdata);
        }

        /* The DiagBIN footer is composed of 12 bytes.  The first 8 bytes is a long value containing the CRC32
         * hash of all data in the DiagBIN file up to this point.  The last 4 bytes in a DiagBIN file is the byte
         * string 'NBGD', indicating the end of a DiagBin file. //todo
         */

        byte[] MZBINTail = new byte[12];
        ByteBuffer.wrap(MZBINTail).asLongBuffer().put(0, crc.getValue());
        ByteBuffer.wrap(MZBINTail).asIntBuffer().put(2, 1312966468);
        dos.write(MZBINTail);

        dos.close();
    }

    public void readDiagBin(ExecutorService executorService, int nThread, boolean loadScans) throws Exception {
        readDiagBinStructure();
        if (loadScans)
            loadDiagBinSpectra(executorService, nThread);
    }

    private void loadDiagBinSpectra(ExecutorService executorService, int nThread) throws Exception {
        ArrayList<Integer> scanNums = new ArrayList<>();
        int maxScans = this.diagBinIndex.length / this.indexWidth - 1;
        for (int i = 0; i <= maxScans; i++) {
            if (ByteBuffer.wrap(this.diagBinIndex).asIntBuffer().get((this.indexWidth / 4) * i + 3) > 0)
                scanNums.add(i);
        }
        loadDiagBinSpectra(executorService, nThread, scanNums);
    }

    private void loadDiagBinSpectra(ExecutorService executorService, int nThread, ArrayList<Integer> scanNums) throws Exception {
        FileChannel fc = FileChannel.open(f.toPath(), StandardOpenOption.READ);
        int factor = Math.max(4, 64 / nThread);
        List<Future> futureList = new ArrayList<>(factor * nThread);
        for (int i = 0; i < factor * nThread; i++) {
            int start = (int) (scanNums.size() * (long) i) / (factor * nThread);
            int end = (int) (scanNums.size() * (long) (i + 1)) / (factor * nThread);
            futureList.add(executorService.submit(() -> {
                ArrayList<DiagnosticRecord> scans = new ArrayList<>();
                try {
                    for (int j = start; j < end; j++) {
                        scans.add(loadDiagBinSpectrum(fc, scanNums.get(j)));
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
                addSpectra(scans);
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
        fc.close();
        Collections.sort(this.spectra);
        //System.out.println(this.spectra.size()); todo why am I missing one spec?
    }

    private DiagnosticRecord loadDiagBinSpectrum(FileChannel fc, int scanNum) throws Exception {
        /* Parse header data and find scan */
        int currIndexPos = (this.indexWidth / 4) * scanNum;
        long offset = ByteBuffer.wrap(this.diagBinIndex).asLongBuffer().get(currIndexPos / 2);
        float dmass = ByteBuffer.wrap(this.diagBinIndex).asFloatBuffer().get(currIndexPos + 2);
        int specLen = ByteBuffer.wrap(this.diagBinIndex).asIntBuffer().get(currIndexPos + 3);

        /* Get byte index and parse */
        byte[] specInfo = new byte[specLen];
        fc.read(ByteBuffer.wrap(specInfo), offset);

        /* Get pepSeq */
        ByteBuffer sibb = ByteBuffer.wrap(specInfo);

        int nChars = sibb.getInt();

        char[] pepSeqChars = new char[nChars];
        for (int i = 0; i < nChars; i++)
            pepSeqChars[i] = sibb.getChar();
        if (nChars % 2 != 0) // Remove padding on seq string
            sibb.getChar();
        String pepSeq = String.valueOf(pepSeqChars);

        /* Get charge */
        int charge = sibb.getInt();

        /* Get modification info */
        int nMods = sibb.getInt();
        float[] mods = new float[nChars + 1]; // 1 based index, 0 = n-term
        for (int i = 0; i < nMods; i++) {
            int modPos = sibb.getInt();
            float modMass = sibb.getFloat();
            mods[modPos] = modMass;
        }

        /* Get immonium peaks */
        int nImm = sibb.getInt();
        float[][] immoniumPeaks = new float[nImm][2];
        for (int i = 0; i < nImm; i++) {
            immoniumPeaks[i][0] = sibb.getFloat();
            immoniumPeaks[i][1] = sibb.getFloat();
        }

        /* Get capY peaks */
        int nCapY = sibb.getInt();
        float[][] capYPeaks = new float[nCapY][2];
        for (int i = 0; i < nCapY; i++) {
            capYPeaks[i][0] = sibb.getFloat();
            capYPeaks[i][1] = sibb.getFloat();
        }

        /* Get squiggle peaks */
        HashMap<Character, float[][]> squigglePeakSets = new HashMap<>();
        for (int i = 0; i < this.ionTypes.size(); i++) {
            int nSquiggleIons = sibb.getInt();
            float[][] squigglePeaks = new float[nSquiggleIons][2];
            for (int j = 0; j < nSquiggleIons; j++) {
                squigglePeaks[j][0] = sibb.getFloat();
                squigglePeaks[j][1] = sibb.getFloat();
            }
            squigglePeakSets.put(this.ionTypes.get(i), squigglePeaks);
        }

        return new DiagnosticRecord(scanNum, this.ionTypes, pepSeq, mods, dmass, charge, immoniumPeaks, capYPeaks, squigglePeakSets);
    }

    private synchronized void addSpectra(ArrayList<DiagnosticRecord> diagnosticRecords) {
        this.spectra.addAll(diagnosticRecords);
    }

    private void readDiagBinStructure() throws IOException {
        RandomAccessFile raf = new RandomAccessFile(f, "r");

        byte[] MZDGhead = new byte[52];
        raf.read(MZDGhead);

        /* Collect metadata from head of file */
        int diagBinId = ByteBuffer.wrap(MZDGhead).asIntBuffer().get(0); //todo add check
        int diagBinVersion = ByteBuffer.wrap(MZDGhead).asIntBuffer().get(1); //todo add check
        this.mzScale = ByteBuffer.wrap(MZDGhead).asIntBuffer().get(2);
        this.indexWidth = ByteBuffer.wrap(MZDGhead).asIntBuffer().get(3);
        int maxScan = ByteBuffer.wrap(MZDGhead).asIntBuffer().get(4);

        boolean[] boolIonTypes = new boolean[8]; //todo add check
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

    private int calculateDiagnosticRecordLength(DiagnosticRecord dr) {
        /* Calculates length of data in file */
        int length = 0;
        length += 4 + dr.pepSeq.length() * 2; // 4 Bs to store int nChars + 2 Bs for each char
        if (length % 4 != 0) // Pad pepSeq to make it 32-bit divisible
            length += 2;
        length += 4; // 4 Bs to store int charge
        length += 4 + dr.modifications.size() * 8; // 4 Bs to store int nMods then (4Bs int pos + 4 Bs float dmass)
        length += 4 + dr.immoniumPeaks.length * 8; // 4 Bs to store int nImmonPs then (4 Bs float int + 4Bs mzdub mz)
        length += 4 + dr.capYPeaks.length * 8; // 4 Bs to store int nCapyPs then (4 Bs float int + 4Bs mzdub mz)
        for (Character ionType : this.ionTypes)
            length += 4 + dr.squigglePeaks.get(ionType).length * 8; // 4 Bs to store int nCapyPs then (4 Bs float int + 4Bs mzdub mz)

        return length;
    }

    public LinkedHashMap<Integer, Float> getDmasses() {
        LinkedHashMap<Integer, Float> scanToDmass = new LinkedHashMap<>();
        int maxScans = this.diagBinIndex.length / this.indexWidth - 1;
        for (int i = 0; i <= maxScans; i++) {
            if (ByteBuffer.wrap(this.diagBinIndex).asIntBuffer().get((this.indexWidth / 4) * i + 3) > 0) {
                float dmass = ByteBuffer.wrap(this.diagBinIndex).asFloatBuffer().get((this.indexWidth / 4) * i + 2);
                scanToDmass.put(i, dmass);
            }
        }
        return scanToDmass;
    }

    //todo if this is slow, only load relevant scans
    public DiagnosticRecord getScan(int scanNum) {
        for (DiagnosticRecord dr : this.spectra) {
            if (dr.scanNum == scanNum) {
                return dr;
            }
        }
        return null;
    }
}
