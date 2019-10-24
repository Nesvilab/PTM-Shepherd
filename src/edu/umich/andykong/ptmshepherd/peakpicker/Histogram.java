package edu.umich.andykong.ptmshepherd.peakpicker;

import java.io.*;
import java.nio.*;
import java.util.*;

import org.apache.commons.math3.distribution.NormalDistribution;

public class Histogram {

	static double [] gweights;
	static double gfrac = 0.95;
	static int BUFSZ = 65536*16;
	
	public int start, end, expSize, binDivs;
	public double [] histo;
	
	public TreeMap<String,double []> merged;
	
	public Histogram(ArrayList<Float> vals, int expSize, int binDivs, int smoothBins) {
		double min = 1e100;
		double max = -1e100;
		
		this.expSize = expSize;
		this.binDivs = binDivs;
		
		for(int i = 0; i < vals.size(); i++) {
			if(vals.get(i) > max)
				max = vals.get(i);
			if(vals.get(i) < min)
				min = vals.get(i);
		}
		
		start = (int)(min-5);
		end = (int)(max + 5);
		
		histo = new double[(end-start)*binDivs];
		calcWeights(smoothBins);
		for(int i = 0; i < vals.size(); i++) {
			int cb = (int)(binDivs*(vals.get(i) + Math.random()/1000000 - 0.0000005 - start + 1.0 / binDivs));
			for(int j = cb - smoothBins/2; j <= (cb + smoothBins/2); j++) {
				histo[j] += gweights[j - (cb - smoothBins/2)];
			}
		}
	}
	
	public Histogram(int start, int end, int binDivs) {
		this.start = start;
		this.end = end;
		this.binDivs = binDivs;
		histo = new double[(end-start)*binDivs];
	}
	
	public Histogram() {
	}
	
	public void mergeHistogram(Histogram h, String saveName) {
		int offset = (h.start - this.start)*this.binDivs;
		double [] nv = new double[(end-start)*binDivs];
		double factor = h.expSize / 1000000.0;
		for(int i = 0; i < h.histo.length; i++)
			nv[offset+i] += (h.histo[i] / factor);
		if(saveName != null) {
			if(merged == null)
				merged = new TreeMap<>();
			merged.put(saveName, nv);
		}
		for(int i = 0; i < nv.length; i++)
			histo[i] += nv[i];
	}
	
	public double [] getOffsets() {
		double [] offsets = new double[(end-start)*binDivs];
		for(int i = 0; i < offsets.length; i++)
			offsets[i] = start + i*(1.0/binDivs);
		return offsets;
	}
	
	public void writeCombinedTSV(File f) throws Exception {
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(f),1<<24));
		out.print("BinCenter");
		String [] keys = new String[merged.size()];
		int cnt = 0;
		for(String cds : merged.keySet()) {
			keys[cnt++] = cds;
			out.print("\t"+cds);
		}
		out.println();
		int len = merged.get(keys[0]).length;
		for(int i = 0; i < len; i++) {
			out.printf("%.5f", start + i*(1.0/binDivs));
			for(int j = 0; j < keys.length; j++)
				out.printf("\t%.8f",merged.get(keys[j])[i]);
			out.println();
		}
		out.close();
	}
	
	public void writeHistogram(File f) throws Exception {
		DataOutputStream dos = new DataOutputStream(new FileOutputStream(f));
		dos.writeInt(start);
		dos.writeInt(end);
		dos.writeInt(binDivs);
		dos.writeInt(expSize);
		writeDouble(dos,histo);
		dos.close();
	}
	
	public static Histogram readHistogram(File f) throws Exception {
		Histogram h = new Histogram();
		DataInputStream dis = new DataInputStream(new FileInputStream(f));
		h.start = dis.readInt();
		h.end = dis.readInt();
		h.binDivs = dis.readInt();
		h.expSize = dis.readInt();
		h.histo = new double[(h.end-h.start)*h.binDivs];
		readDouble(dis, h.histo);		
		dis.close();
		return h;
	}
	
	public static Histogram readHistogramHeader(File f) throws Exception {
		Histogram h = new Histogram();
		DataInputStream dis = new DataInputStream(new FileInputStream(f));
		h.start = dis.readInt();
		h.end = dis.readInt();
		h.binDivs = dis.readInt();
		h.expSize = dis.readInt();
		dis.close();
		return h;
	}
	
    public static void readDouble(DataInputStream dis, double [] vals) throws Exception {
        byte [] cbuf = new byte[BUFSZ];
        int cpos = 0;
        while(cpos < vals.length) {
                int nread = dis.read(cbuf,0,8*Math.min(BUFSZ/8,(vals.length - cpos)));
                ByteBuffer bb = ByteBuffer.wrap(cbuf);
                DoubleBuffer fb = bb.asDoubleBuffer();
                fb.get(vals, cpos, nread / 8);
                cpos += (nread / 8);
        }
    }
	
	public static void writeDouble(DataOutputStream dos, double [] vals) throws Exception {
        ByteBuffer cbuf = ByteBuffer.allocate(BUFSZ);
        int cpos = 0;
        while(cpos < vals.length) {
                int nEle = Math.min(BUFSZ/8, vals.length - cpos);
                cbuf.asDoubleBuffer().put(vals,cpos,nEle);
                dos.write(cbuf.array(),0,8*nEle);
                cpos += nEle;
        }
	}

	public static void calcWeights(int Nbins) {
		NormalDistribution nd = new NormalDistribution();
		double a = (1-gfrac)/2;
		double lower = nd.inverseCumulativeProbability(a);
		double upper = nd.inverseCumulativeProbability(1-a);
		double width = (upper - lower) / Nbins;
		gweights = new double[Nbins];
		for(int i = 0; i < Nbins; i++) {
			double lo = nd.cumulativeProbability(lower+i*width);
			double hi = nd.cumulativeProbability(lower +(i+1)*width);
			gweights[i] = hi-lo;
			//System.out.println(gweights[i]);
		}
	}
}
