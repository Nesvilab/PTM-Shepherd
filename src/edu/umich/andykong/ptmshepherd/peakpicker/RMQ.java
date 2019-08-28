package edu.umich.andykong.ptmshepherd.peakpicker;

import java.util.Arrays;


public class RMQ {

	double [] vals;
	int BLOCKSZ;
	RMQLinearithmic lmql;
	
	public RMQ(double [] v) {
		vals = v;
		BLOCKSZ = (int)((Math.log(v.length)/Math.log(2))/4)+1;
		double [] blockMins = new double[vals.length / BLOCKSZ + 1];
		Arrays.fill(blockMins, 1e100);
		for(int i = 0; i < v.length; i++) 
			blockMins[i / BLOCKSZ] = Math.min(blockMins[i / BLOCKSZ], v[i]);
		lmql = new RMQLinearithmic(blockMins);
	}
	
	public double query(int l, int r) {
		double res = 1e100;
		int a = l / BLOCKSZ + 1;
		int b = r / BLOCKSZ - 1;
		for(int i = l; i <= r && i < a*BLOCKSZ; i++)
			res = Math.min(res,vals[i]);
		for(int i = Math.max(l,(b+1)*BLOCKSZ); i <= r; i++)
			res = Math.min(res,vals[i]);
		if(b >= a)
			res = Math.min(res,lmql.query(a, b));
		return res;
	}
	
}
