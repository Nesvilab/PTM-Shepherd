package edu.umich.andykong.ptmshepherd.core;

public class AAMasses {
	
	public static final float protMass = (float)1.00727647;
	public static final float monoisotopic_nterm_mass = (float) 1.0078250321;
	//public static final double monoisotopic_nterm_masses = []{,1.0078250321,}
	//public static final float average_nterm_mass = (float) 1.0078250321;
	public static final float monoisotopic_cterm_mass = (float) 17.0027396542;
	//public static final float average_cterm_mass = (float) 17.0027396542;
	//public static final double monisotopic_cterm_mass = []{18.03386,17.0027396542,
	public static final float averagineIsotopeMass = 1.00235f;

	public static double getMonoisotopicNeutralMass(String pep) {
		double cmass = AAMasses.monoisotopic_nterm_mass + AAMasses.monoisotopic_cterm_mass;
		for(int i = 0; i < pep.length(); i++)
			cmass += AAMasses.monoisotopic_masses[pep.charAt(i)-'A'];
		return cmass;
	}
	
	public static final float monoisotopic_masses[] = {
	  (float) 71.03711, // A
	  (float) 0, // B
	  (float) 103.00919, // C
	  (float) 115.02694, // D
	  (float) 129.04259, // E
	  (float) 147.06841, // F
	  (float) 57.02146, // G
	  (float) 137.05891, // H
	  (float) 113.08406, // I
	  (float) 0, // J
	  (float) 128.09496, // K
	  (float) 113.08406, // L
	  (float) 131.04049, // M
	  (float) 114.04293, // N
	  (float) 114.07931, // O
	  (float) 97.05276, // P
	  (float) 128.05858, // Q
	  (float) 156.10111, // R
	  (float) 87.03203, // S
	  (float) 101.04768, // T
	  (float) 0, // U
	  (float) 99.06841, // V
	  (float) 186.07931, // W
	  (float) 0, // X
	  (float) 163.06333, //  Y
	  (float) 0, // Z
	  (float) 10000.00000 //[
	};

	public static final float ionTypeShifts[] = {
		(float) -27.99591562, //a
		(float) 0.0, //b
		(float) 17.0265491, //c
		(float) 43.98982924, //x
		(float) AAMasses.monoisotopic_cterm_mass + AAMasses.monoisotopic_nterm_mass, //y
		(float) 1.991840615 //z
	};
}
