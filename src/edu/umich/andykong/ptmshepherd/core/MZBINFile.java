package edu.umich.andykong.ptmshepherd.core;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.concurrent.*;
import java.util.zip.CRC32;

public class MZBINFile {

    private static final int MZBINVersion = 20190716;

    public final File f;
    public final String runName;
    public final ArrayList<Spectrum> specs;

    private String instrumentModel = "";
    private String instrumentSerial = "";
    private int runStartTime = 0;
    private String mzMLHeader = "";
    private String rawSHA1 = "";
    private String activationMethods = "";
    private String batmassVersion = "";
    private String timsLibVersion = "";

    private int INDEXWidth;
    private int MZScale;
    private byte[] MZBINIndex;


    public String getBasename(String f) {
        String baseName = f.substring(0,f.lastIndexOf("."));
        return baseName;
    }

    public static String getBasenameStat(String f) {
        String baseName = f.substring(0,f.lastIndexOf("."));
        return baseName;
    }

    public MZBINFile(ExecutorService executorService, int nThread, String filePath, boolean loadScans) throws Exception {
        this(executorService, nThread, new File(filePath), loadScans);
    }

    public MZBINFile(ExecutorService executorService, int nThread, File f, boolean loadScans) throws Exception {
        this.f = f;
        runName = getBasename(f.getName());
        specs = new ArrayList<>();
        readMZBIN(executorService, nThread, loadScans);
    }

    public MZBINFile(String filePath, Collection<Spectrum> specs, String batmassVersion, String timsLibVersion) {
        f = new File(filePath);
        runName = getBasename(f.getName());
        this.batmassVersion = batmassVersion;
        this.timsLibVersion = timsLibVersion;
        this.specs = new ArrayList<>(specs);
        Collections.sort(this.specs);
    }

    public MZBINFile(File filePath, Collection<Spectrum> specs, String batmassVersion, String timsLibVersion) {
        f = filePath;
        runName = getBasename(f.getName());
        this.batmassVersion = batmassVersion;
        this.timsLibVersion = timsLibVersion;
        this.specs = new ArrayList<>(specs);
        Collections.sort(this.specs);
    }

    public static Spectrum loadSpectrum(String path, int scanNum) throws Exception {
        String runName = getBasenameStat(path);

        RandomAccessFile raf = new RandomAccessFile(path, "r");

        byte[] MZBNhead = new byte[28];
        raf.read(MZBNhead);

        int MZBNid = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(0);
        int MZBNversion = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(1);
        long headerPos = ByteBuffer.wrap(MZBNhead).asLongBuffer().get(1);
        int maxScan = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(3);
        int MZSCALE = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(4);
        int INDEXwidth = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(5);
        maxScan = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(6);

        if (scanNum > maxScan)
            return null;

        byte[] scanHeader = new byte[INDEXwidth];
        raf.seek(28 + INDEXwidth * scanNum);
        raf.read(scanHeader);

        long offset = ByteBuffer.wrap(scanHeader).asLongBuffer().get(0);
        int nfrags = ByteBuffer.wrap(scanHeader).asIntBuffer().get(2);

        if (nfrags == 0)
            return null;

        float retentionTime = ByteBuffer.wrap(scanHeader).asFloatBuffer().get(3);
        float injectionTime = ByteBuffer.wrap(scanHeader).asFloatBuffer().get(4);
        byte msLevel = (byte) ByteBuffer.wrap(scanHeader).asIntBuffer().get(5);
        byte charge = (byte) ByteBuffer.wrap(scanHeader).asIntBuffer().get(6);
        double precursorTarget = __mzIntToDouble(ByteBuffer.wrap(scanHeader).asIntBuffer().get(7), MZSCALE);
        float precursorIntensity = ByteBuffer.wrap(scanHeader).asFloatBuffer().get(8);
        double precursorLow = __mzIntToDouble(ByteBuffer.wrap(scanHeader).asIntBuffer().get(9), MZSCALE);
        double precursorHigh = __mzIntToDouble(ByteBuffer.wrap(scanHeader).asIntBuffer().get(10), MZSCALE);
        int parentScan = ByteBuffer.wrap(scanHeader).asIntBuffer().get(11);
        double MZmin = __mzIntToDouble(ByteBuffer.wrap(scanHeader).asIntBuffer().get(12), MZSCALE);
        double MZmax = __mzIntToDouble(ByteBuffer.wrap(scanHeader).asIntBuffer().get(13), MZSCALE);
        double basePeakMZ = __mzIntToDouble(ByteBuffer.wrap(scanHeader).asIntBuffer().get(14), MZSCALE);
        float basePeakInt = ByteBuffer.wrap(scanHeader).asFloatBuffer().get(15);
        float TIC = ByteBuffer.wrap(scanHeader).asFloatBuffer().get(16);
        double selectedMass = __mzIntToDouble(ByteBuffer.wrap(scanHeader).asIntBuffer().get(17), MZSCALE);
        int activationMethod = ByteBuffer.wrap(scanHeader).asIntBuffer().get(18);
        float precursorIM = ByteBuffer.wrap(scanHeader).asFloatBuffer().get(19);

        byte[] rowData = new byte[8 * nfrags];
        raf.seek(offset + 8);
        raf.read(rowData);
        float[] peakMZ = new float[nfrags];
        float[] peakInt = new float[nfrags];
        for (int j = 0; j < nfrags; j++) {
            peakMZ[j] = (float) __mzIntToDouble(ByteBuffer.wrap(rowData).asIntBuffer().get(2 * j), MZSCALE);
            peakInt[j] = ByteBuffer.wrap(rowData).asFloatBuffer().get(2 * j + 1);
        }
        raf.close();

        return new Spectrum(runName + "." + scanNum + "." + scanNum, scanNum, charge, msLevel, selectedMass, retentionTime, peakMZ, peakInt);
    }

    public void writeMZBIN() throws Exception {

        // Binary values are stored in network order (default for Java)
        // All values are signed

        // mzBIN file contains 5 major sections
        // 1. mzBIN header (28 bytes)
        // 2. spectra index  ((maxScanNumber + 1) * INDEXWidth) bytes)
        // 3. spectra peak data
        // 4. mzBIN textual metadata
        // 5. mzBIN footer (12 bytes)

        // Compute maximum scan number (as numbered by the vendor library)

        int maxScan = 0;
        for (Spectrum spectrum : specs) {
            if (maxScan < spectrum.scanNum)
                maxScan = spectrum.scanNum;
        }

        // The mzBIN header consists of 28 bytes
        //
        // Bytes 0-4   : 'MZBN' - to indicate that this is a MZBN file
        //       4-8   : 0 - 4 byte integer indicating the current MZBN version number
        //       8-16  : 8 byte long indicating file offset of the text metadata store
        //       16-20 : 4 byte integer containing the scaling factor (MZscale) to perform integer <-> doubles for M/Z values (currently 200000)
        //       20-24 : 4 byte integer indicating the width (in bytes) of the spectra index (currently 80)
        //       24-28 : 4 byte integer containing the maximum scan ID as numbered by vendor library

        // The scaled int data type is a m/z value stored in the mzBIN format as an integer
        // conversion between a floating point m/z value to an integer is defined as
        // <scaled int m/z> = round(float m/z / MZscale)
        // conversion in the unpacking process is simply a multiplication by MZscale

        // MZscale is chosen so that the largest m/z value is less than 2^31-1
        // it is currently defaulted to 200,000

        INDEXWidth = 80;
        MZScale = 200000;

        byte[] MZBNHead = new byte[28];

        // spectra index contains (maxScan+1) rows of width INDEXWidth bytes starting at byte 29
        // the first 80 bytes are currently specified

        // Bytes  0-8   : peak data offset (long)
        //        8-12  : number of peaks (int)
        //       12-16  : retention time in minutes (float)
        //       16-20  : injection time in milliseconds (float)
        //       20-24  : MS level (int)
        //       24-28  : charge (int)
        //       28-32  : precursor target m/z (scaled int)
        //       32-36  : precursor intensity (float)
        //       36-40  : precursor window lower bound (scaled int)
        //       40-44  : precursor window upper bound (scaled int)
        //       44-48  : parent scan number (int)
        //       48-52  : lowest detected M/Z (scaled int)
        //       52-56  : highest detected M/Z (scaled int)
        //       56-60  : base peak M/Z (scaled int)
        //       60-64  : base peak intensity (float)
        //       64-68  : total ion current (float)
        //       68-72  : selected ion m/z (scaled int)
        //       72-76  : activation method (int)
        //       76-80  : precursor ion mobility (float)

        MZBINIndex = new byte[(maxScan + 1) * INDEXWidth];
        long offset = MZBNHead.length + MZBINIndex.length;
        for (Spectrum cs : specs) {
            int cp = (INDEXWidth / 4) * (cs.scanNum);
            ByteBuffer.wrap(MZBINIndex).asLongBuffer().put(cp / 2, offset);
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 2, cs.peakInt.length);
            ByteBuffer.wrap(MZBINIndex).asFloatBuffer().put(cp + 3, (float) cs.rt);
            // ByteBuffer.wrap(MZBINIndex).asFloatBuffer().put(cp+4, (float)cs.injectionTime);
            ByteBuffer.wrap(MZBINIndex).asFloatBuffer().put(cp + 4, 0f);
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 5, cs.msLevel);
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 6, cs.charge);
            // ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp+7, __mzDoubleToInt(cs.precursorTarget));
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 7, __mzDoubleToInt(0));
            //ByteBuffer.wrap(MZBINIndex).asFloatBuffer().put(cp + 8, cs.precursorIntensity);
            // ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp+9, __mzDoubleToInt(cs.precursorLow));
            // ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp+10, __mzDoubleToInt(cs.precursorHigh));
            // ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp+11, cs.parentScan);
            // ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp+12,__mzDoubleToInt(cs.MZmin));
            // ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp+13,__mzDoubleToInt(cs.MZmax));
            // ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp+14,__mzDoubleToInt(cs.basePeakMZ));
            // ByteBuffer.wrap(MZBINIndex).asFloatBuffer().put(cp+15, cs.basePeakInt);
            // ByteBuffer.wrap(MZBINIndex).asFloatBuffer().put(cp+16, cs.TIC);
            // ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp+17, __mzDoubleToInt(cs.selectedMass));
            // ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp+18, cs.activationMethod);
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 9, 0);
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 10, 0);
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 11, 0);
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 12, 0);
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 13, 0);
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 14, 0);
            ByteBuffer.wrap(MZBINIndex).asFloatBuffer().put(cp + 15, 0);
            ByteBuffer.wrap(MZBINIndex).asFloatBuffer().put(cp + 16, 0);
            //ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 17, __mzDoubleToInt(cs.precursorMz));
            ByteBuffer.wrap(MZBINIndex).asIntBuffer().put(cp + 18, 0);
            //ByteBuffer.wrap(MZBINIndex).asFloatBuffer().put(cp + 19, cs.precursorIM);
            offset += 8 * (cs.peakInt.length + 1);
        }

        ByteBuffer.wrap(MZBNHead).asIntBuffer().put(0, 1297760846);
        ByteBuffer.wrap(MZBNHead).asIntBuffer().put(1, MZBINVersion);
        ByteBuffer.wrap(MZBNHead).asLongBuffer().put(1, offset);
        ByteBuffer.wrap(MZBNHead).asIntBuffer().put(4, MZScale);
        ByteBuffer.wrap(MZBNHead).asIntBuffer().put(5, INDEXWidth);
        ByteBuffer.wrap(MZBNHead).asIntBuffer().put(6, maxScan);


        // start writing stuff
        DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f), 1 << 24));
        CRC32 crc = new CRC32();
        dos.write(MZBNHead);
        crc.update(MZBNHead);
        dos.write(MZBINIndex);
        crc.update(MZBINIndex);

        // the spectra peak lists are stored as records of 8 + 8*(number of peaks) bytes
        // Bytes 0-4 contain the vendor scan number and bytes 4-8 contain the number of peaks
        // This is followed by 8 byte peaks in the form of [ m/z (scaled int), intensity (float)]

        for (Spectrum cs : specs) {
            byte[] cdata = new byte[8 * (cs.peakInt.length + 1)];
            ByteBuffer.wrap(cdata).asIntBuffer().put(0, cs.scanNum);
            ByteBuffer.wrap(cdata).asIntBuffer().put(1, cs.peakInt.length);
            for (int j = 0; j < cs.peakInt.length; j++) {
                ByteBuffer.wrap(cdata).asIntBuffer().put(2 * j + 2, __mzDoubleToInt(cs.peakMZ[j]));
                ByteBuffer.wrap(cdata).asFloatBuffer().put(2 * j + 3, cs.peakInt[j]);
            }
            dos.write(cdata);
            crc.update(cdata);
        }

        // The textual metadata section is a dictionary for storing strings and other metadata
        // regarding the entire run.  It consists of 5 defined key-values but can be used to store
        // additional user/program annotations/metadata.

        // rawSHA1 (required) - SHA1 hash of RAW file

        // activationMethods (recommended) - activation methods
        // instModel (recommended) - instrument model
        // instSerial (recommended) - instrument serial number
        // instSerial (recommended) - instrument serial number
        // runStart (recommended) - run start time in seconds from UNIX epoch

        // mzMLheader (optional) - the contents of the mzML file from the start to the end of the <run ... > open tag

        // the textual metadata starts with a 4 byte integer representing the length of the metadata data section
        // the metadata data section is composed of key-value records represented as follows
        // for a key-value pair where the key is N bytes and value is M bytes in the UTF-8 character set
        //		N - 4 byte integer
        //		key - N bytes representing the key
        //		M - 4 byte integer
        //		value - M bytes representing the value
        //
        // see the method __packHeader for implementation details


        byte[] MZBNHeader = __packBINHeader();
        byte[] MZBNHeaderLength = new byte[4];
        ByteBuffer.wrap(MZBNHeaderLength).asIntBuffer().put(MZBNHeader.length);

        dos.write(MZBNHeaderLength);
        crc.update(MZBNHeaderLength);
        dos.write(MZBNHeader);
        crc.update(MZBNHeader);

        // The mzBIN footer is composed of 12 bytes.  The first 8 bytes is a long value containing the CRC32
        // hash of all data in the mzBIN file up to this point.  The last 4 bytes in a MZBN file is the byte
        // string 'NBZM', indicating the end of a MZBN file.

        byte[] MZBINTail = new byte[12];
        ByteBuffer.wrap(MZBINTail).asLongBuffer().put(0, crc.getValue());
        ByteBuffer.wrap(MZBINTail).asIntBuffer().put(2, 1312971341);
        dos.write(MZBINTail);
        dos.close();
    }

    public String getBatmassVersion() {
        return batmassVersion;
    }

    public String getTimslibVersion() {
        return timsLibVersion;
    }

    private static double __mzIntToDouble(int mz, int scale) {
        return ((double) mz) / scale;
    }

//	public void readMZBIN() throws Exception {
//		DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(f), 1 << 24));
//		CRC32 crc = new CRC32();
//
//		byte [] MZBINHead = new byte[28];
//		dis.read(MZBINHead);
//		crc.update(MZBINHead);
//
//		int MZBNid = ByteBuffer.wrap(MZBINHead).asIntBuffer().get(0);
//		int MZBNversion = ByteBuffer.wrap(MZBINHead).asIntBuffer().get(1);
//		long headerPos = ByteBuffer.wrap(MZBINHead).asLongBuffer().get(1);
//		int maxScan = ByteBuffer.wrap(MZBINHead).asIntBuffer().get(3);
//		MZScale = ByteBuffer.wrap(MZBINHead).asIntBuffer().get(4);
//		INDEXWidth = ByteBuffer.wrap(MZBINHead).asIntBuffer().get(5);
//		maxScan = ByteBuffer.wrap(MZBINHead).asIntBuffer().get(6);
//
//		if(MZBNid != 1297760846) {
//			System.err.println("MZBN header not found!");
//			System.exit(1);
//		}
//		if(MZBNversion != MZBINVersion) {
//			System.err.printf("Incompatible MZBN version %d (library is %d)!",MZBNversion, MZBINVersion);
//			System.exit(1);
//		}
//
//		MZBINIndex = new byte[(maxScan+1)*(INDEXWidth)];
//		dis.read(MZBINIndex);
//		crc.update(MZBINIndex);
//
//		specs = new ArrayList<>();
//
//		for(int i = 0; i <= maxScan; i++) {
//			int cp = (INDEXWidth /4)*i;
//			int nfrags = ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+2);
//			if(nfrags == 0)
//				continue;
//
//			float retentionTime = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp+3);
//			float injectionTime = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp+4);
//			byte msLevel = (byte) ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+5);
//			byte charge = (byte) ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+6);
//			double precursorTarget = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+7));
//			float precursorIntensity = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp+8);
//			double precursorLow = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+9));
//			double precursorHigh = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+10));
//			int parentScan = ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+11);
//			double MZmin = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+12));
//			double MZmax = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+13));
//			double basePeakMZ = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+14));
//			float basePeakInt = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp+15);
//			float TIC = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp+16);
//			double selectedMass = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+17));
//			int activationMethod = ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp+18);
//			specs.add(new Spectrum(nfrags, runName + "." + i + "." + i + "." + charge, i, charge, msLevel, retentionTime, selectedMass, precursorIntensity));
//		}
//
//		int offset = MZBINHead.length + MZBINIndex.length;
//		byte [] rowHeader = new byte[8];
//		for(int i = 0; i < specs.size(); i++) {
//			dis.read(rowHeader);
//			crc.update(rowHeader);
//			int scanNum = ByteBuffer.wrap(rowHeader).asIntBuffer().get(0);
//			int nfrags = ByteBuffer.wrap(rowHeader).asIntBuffer().get(1);
//
//			Spectrum cs = specs.get(i);
//			// verify offset
//			if(ByteBuffer.wrap(MZBINIndex).asLongBuffer().get((INDEXWidth /8)*scanNum) != offset) {
//				System.err.printf("Index corrupted! Scan %d at %d instead of %d in index.\n",scanNum,offset,
//						ByteBuffer.wrap(MZBINIndex).asIntBuffer().get((INDEXWidth /8)*scanNum));
//				System.exit(1);
//			}
//			// verify scan number and peaks
//			if(scanNum != cs.scanNum) {
//				System.err.printf("Data corruption! Scan %d: scanNum %d != scanNum %d.\n",i+1,scanNum,cs.scanNum);
//				System.exit(1);
//			}
//			if(nfrags != cs.peakInt.length) {
//				System.err.printf("Data corruption! Scan %d: fragments %d != fragments %d.\n",i+1,nfrags,cs.peakInt.length);
//				System.exit(1);
//			}
//
//			byte [] rowData = new byte[8*nfrags];
//			dis.read(rowData);
//			crc.update(rowData);
//
//			for(int j = 0; j < nfrags; j++) {
//				cs.peakMZ[j] = (float) __mzIntToDouble(ByteBuffer.wrap(rowData).asIntBuffer().get(2*j));
//				cs.peakInt[j] = ByteBuffer.wrap(rowData).asFloatBuffer().get(2*j+1);
//			}
//
//			offset += 8*(cs.peakInt.length+1);
//		}
//
//		// verify MZBN header offset
//		if(headerPos != offset) {
//			System.err.printf("Data corruption! MZBN metadata found at %d instead of %d\n",offset, headerPos);
//			System.exit(1);
//		}
//
//		byte [] MZBNheader_len = new byte[4];
//		dis.read(MZBNheader_len);
//		crc.update(MZBNheader_len);
//
//		byte [] MZBNheader = new byte[ByteBuffer.wrap(MZBNheader_len).asIntBuffer().get(0)];
//		dis.read(MZBNheader);
//		crc.update(MZBNheader);
//		__unpackBINHeader(MZBNheader);
//
//		byte [] MZBNtail = new byte[12];
//		dis.read(MZBNtail);
//
//		// verify CRC
//		long expectedCRC = ByteBuffer.wrap(MZBNtail).asLongBuffer().get(0);
//		if(crc.getValue() != expectedCRC) {
//			System.err.printf("Data corruption! MZBN CRC mismatch %d instead of expected %d\n",crc.getValue(),expectedCRC);
//			System.exit(1);
//		}
//		if(ByteBuffer.wrap(MZBNtail).asIntBuffer().get(2) != 1312971341) {
//			System.err.println("NBZM tail not found!");
//			System.exit(1);
//		}
//		dis.close();
//	}

    private static byte[] __packHeader(HashMap<String, String> header) throws Exception {
        ArrayList<String> strings = new ArrayList<>();
        for (String key : header.keySet()) {
            strings.add(key);
            strings.add(header.get(key));
        }
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        for (String s : strings) {
            ByteBuffer len = ByteBuffer.allocate(4);
            byte[] data = s.getBytes(StandardCharsets.UTF_8);
            len.asIntBuffer().put(data.length);
            bos.write(len.array());
            bos.write(data);
        }
        return bos.toByteArray();
    }

    private HashMap<String, String> __unpackHeader(byte[] header) throws Exception {
        ArrayList<String> strings = new ArrayList<>();
        HashMap<String, String> headers = new HashMap<>();
        ByteArrayInputStream bis = new ByteArrayInputStream(header);
        while (bis.available() > 0) {
            ByteBuffer len = ByteBuffer.allocate(4);
            bis.read(len.array());
            byte[] data = new byte[len.asIntBuffer().get(0)];
            bis.read(data);
            strings.add(new String(data, StandardCharsets.UTF_8));
        }
        for (int i = 0; i < strings.size(); i += 2) {
            headers.put(strings.get(i), strings.get(i + 1));
        }
        return headers;
    }

    private String getString(HashMap<String, String> header, String key) {
        if (header.containsKey(key))
            return header.get(key);
        return "";
    }

    private int getInt(HashMap<String, String> header, String key) {
        if (header.containsKey(key))
            return Integer.parseInt(header.get(key));
        return 0;
    }

    private void readMZBIN(ExecutorService executorService, int nThread, boolean loadScans) throws Exception {
        readMZBINStructure();
        if (loadScans) {
            loadMZBINScans(executorService, nThread);
        }
    }

    private synchronized void addScans(ArrayList<Spectrum> ispecs) {
        specs.addAll(ispecs);
    }

    private void loadMZBINScans(ExecutorService executorService, int nThread, ArrayList<Integer> scanNums) throws Exception {
        FileChannel fc = FileChannel.open(f.toPath(), StandardOpenOption.READ);
        int factor = Math.max(4, 64 / nThread);
        List<Future> futureList = new ArrayList<>(factor * nThread);
        for (int i = 0; i < factor * nThread; i++) {
            int start = (int) (scanNums.size() * (long) i) / (factor * nThread);
            int end = (int) (scanNums.size() * (long) (i + 1)) / (factor * nThread);
            futureList.add(executorService.submit(() -> {
                ArrayList<Spectrum> scans = new ArrayList<>();
                try {
                    for (int j = start; j < end; j++) {
                        scans.add(loadMZBINSpectrum(fc, scanNums.get(j)));
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
                addScans(scans);
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
        fc.close();
        Collections.sort(specs);
    }

    private void loadMZBINScans(ExecutorService executorService, int nThread) throws Exception {
        ArrayList<Integer> scanNums = new ArrayList<>();
        int maxScans = MZBINIndex.length / INDEXWidth - 1;
        for (int i = 0; i <= maxScans; i++) {
            if (ByteBuffer.wrap(MZBINIndex).asIntBuffer().get((INDEXWidth / 4) * i + 2) > 0)
                scanNums.add(i);
        }
        loadMZBINScans(executorService, nThread, scanNums);
    }

    private void readMZBINStructure() throws Exception {
        RandomAccessFile raf = new RandomAccessFile(f, "r");

        byte[] MZBNhead = new byte[28];
        raf.read(MZBNhead);

        int MZBNid = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(0);
        int MZBNversion = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(1);
        long headerPos = ByteBuffer.wrap(MZBNhead).asLongBuffer().get(1);
        int maxScan = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(3);
        MZScale = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(4);
        INDEXWidth = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(5);
        maxScan = ByteBuffer.wrap(MZBNhead).asIntBuffer().get(6);

        if (MZBNid != 1297760846) {
            System.err.println("MZBN header not found!");
            System.exit(1);
        }
        if (MZBNversion != MZBINVersion) {
            System.err.printf("Incompatible MZBN version %d (library is %d)!", MZBNversion, MZBINVersion);
            System.exit(1);
        }

        MZBINIndex = new byte[(maxScan + 1) * INDEXWidth];
        raf.read(MZBINIndex);

        raf.seek(headerPos);
        byte[] MZBNheader_len = new byte[4];
        raf.read(MZBNheader_len);

        byte[] MZBNheader = new byte[ByteBuffer.wrap(MZBNheader_len).asIntBuffer().get(0)];
        raf.read(MZBNheader);
        __unpackBINHeader(MZBNheader);

        raf.close();
    }

    private Spectrum loadMZBINSpectrum(FileChannel fc, int scanNum) throws Exception {
        int cp = (INDEXWidth / 4) * scanNum;
        int nfrags = ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 2);

        long offset = ByteBuffer.wrap(MZBINIndex).asLongBuffer().get(cp / 2);
        float retentionTime = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp + 3);
        float injectionTime = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp + 4);
        byte msLevel = (byte) ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 5);
        byte charge = (byte) ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 6);
        double precursorTarget = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 7));
        float precursorIntensity = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp + 8);
        double precursorLow = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 9));
        double precursorHigh = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 10));
        int parentScan = ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 11);
        double MZmin = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 12));
        double MZmax = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 13));
        double basePeakMZ = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 14));
        float basePeakInt = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp + 15);
        float TIC = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp + 16);
        double selectedMass = __mzIntToDouble(ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 17));
        int activationMethod = ByteBuffer.wrap(MZBINIndex).asIntBuffer().get(cp + 18);
        float precursorIM = ByteBuffer.wrap(MZBINIndex).asFloatBuffer().get(cp + 19);

        byte[] rowData = new byte[8 * nfrags];
        fc.read(ByteBuffer.wrap(rowData), offset + 8);

        float[] peakMZ = new float[nfrags];
        float[] peakInt = new float[nfrags];
        for (int j = 0; j < nfrags; j++) {
            peakMZ[j] = (float) __mzIntToDouble(ByteBuffer.wrap(rowData).asIntBuffer().get(2 * j), MZScale);
            peakInt[j] = ByteBuffer.wrap(rowData).asFloatBuffer().get(2 * j + 1);
        }

        return new Spectrum(runName + "." + scanNum + "." + scanNum, scanNum, charge, msLevel, selectedMass, retentionTime, peakMZ, peakInt);
    }

    private int __mzDoubleToInt(double mz) {
        if (Math.round(mz * MZScale) > Integer.MAX_VALUE) {
            System.err.printf("The mz %.4f exceeds allowed maximum value using currently scaling factor of %d. Check your data please.\n", mz, MZScale);
            System.exit(1);
        }
        return (int) (Math.round(mz * MZScale));
    }

    private double __mzIntToDouble(int mz) {
        return ((double) mz) / MZScale;
    }

    private byte[] __packBINHeader() throws Exception {
        HashMap<String, String> header = new HashMap<>();
        header.put("instModel", instrumentModel);
        header.put("instSerial", instrumentSerial);
        header.put("runStart", "" + runStartTime);
        header.put("rawSHA1", rawSHA1);
        header.put("mzMLHeader", mzMLHeader);
        header.put("activationMethods", activationMethods);
        header.put("batmassVersion", batmassVersion);
        header.put("timsLibVersion", timsLibVersion);
        return __packHeader(header);
    }

    private void __unpackBINHeader(byte[] headers) throws Exception {
        HashMap<String, String> header = __unpackHeader(headers);
        instrumentModel = getString(header, "instModel");
        instrumentSerial = getString(header, "instSerial");
        rawSHA1 = getString(header, "rawSHA1");
        mzMLHeader = getString(header, "mzMLHeader");
        runStartTime = getInt(header, "runStart");
        activationMethods = getString(header, "activationMethods");
        batmassVersion = getString(header, "batmassVersion");
        timsLibVersion = getString(header, "timsLibVersion");
    }
}
