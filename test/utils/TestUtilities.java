package utils;

import edu.umich.andykong.ptmshepherd.PSMFile;


public class TestUtilities {

    public static int countGlycoPSMs(PSMFile psmFile, double glycoFDR) {
        int totalPSMs = 0;
        int glycoPSMs = 0;
        int glycoPSMsBelowThreshold = 0;

        for (int i = 0; i < psmFile.psms.size(); i++) {
            totalPSMs++;
            if (psmFile.psms.get(i).glycanAssignmentResult != null && psmFile.psms.get(i).glycanAssignmentResult.foundGlycan) {
                glycoPSMs++;
                if (psmFile.psms.get(i).glycanAssignmentResult.glycanQval <= glycoFDR) {
                    glycoPSMsBelowThreshold++;
                }
            }
        }

        System.out.printf("%s: Total PSMs: %d, glycoPSMs: %d, Glyco PSMs below FDR threshold (%.2f): %d%n", psmFile.fname, totalPSMs, glycoPSMs, glycoFDR, glycoPSMsBelowThreshold);
        return glycoPSMsBelowThreshold;
    }

}
