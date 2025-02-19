package psm;

import edu.umich.andykong.ptmshepherd.PSM;
import edu.umich.andykong.ptmshepherd.PSMFile;
import org.junit.jupiter.api.Test;

import java.io.File;

public class PSMFileTest {

    private static final double tol = 0.0001;

    @Test
    public void parsePSMFileTest() {
        File testFile = new File("test-resources/test_psms.tsv");
        PSMFile psmFile = null;
        try {
            psmFile = new PSMFile(testFile, 0);
        } catch (Exception e) {
            e.printStackTrace();
        }

        // are headers found?
        assert psmFile.specCol != -1;
        assert psmFile.peptideCol != -1;
        assert psmFile.dMassCol != -1;
        assert psmFile.assignedModCol != -1;

        PSM psm1 = psmFile.psms.get(0);
        assert psm1.getScanNum() == 382;
        assert psm1.getDMass() + -0.0045 < tol;
        assert psm1.getAssignedMods().isEmpty();

        PSM psm4 = psmFile.psms.get(3);
        assert psm4.getScanNum() == 2903;
        assert psm4.getDMass() - 0.0048 - 2512.8455 < tol;
        assert psm4.getAssignedMods().isEmpty();
    }

    @Test
    public void parsePSMFileRemoveDeltaTest() {
        File testFile = new File("test-resources/test_psms_remove-delta.tsv");
        PSMFile psmFile = null;
        try {
            psmFile = new PSMFile(testFile, 1);
        } catch (Exception e) {
            e.printStackTrace();
        }

        assert psmFile != null;
        PSM psm1 = psmFile.psms.get(0);
        assert psm1.getScanNum() == 382;
        assert psm1.getDMass() + -0.0045 < tol;
        assert psm1.getAssignedMods().isEmpty();

        PSM psm4 = psmFile.psms.get(3);
        assert psm4.getScanNum() == 2903;
        assert psm4.getDMass() - 0.0048 - 2512.8455 < tol;
        assert psm4.getOriginalDeltaMass() - 0.0048 < tol;
        assert psm4.getAssignedMods().isEmpty();
        assert psm4.getOriginalAssignedMods().get(3) - 2512.8455 < tol;
    }

    @Test
    public void parsePSMFileKeepDeltaTest() {
        File testFile = new File("test-resources/test_psms_keep-delta.tsv");
        PSMFile psmFile = null;
        try {
            psmFile = new PSMFile(testFile, 2);
        } catch (Exception e) {
            e.printStackTrace();
        }

        assert psmFile != null;
        PSM psm1 = psmFile.psms.get(0);
        assert psm1.getScanNum() == 382;
        assert psm1.getDMass() + -0.0045 < tol;
        assert psm1.getAssignedMods().isEmpty();

        PSM psm4 = psmFile.psms.get(3);
        assert psm4.getScanNum() == 2903;
        assert psm4.getDMass() - 0.0048 - 2512.8455 < tol;
        assert psm4.getOriginalDeltaMass() - 0.0048 - 2512.8455 < tol;
        assert psm4.getAssignedMods().isEmpty();
        assert psm4.getOriginalAssignedMods().get(3) - 2512.8455 < tol;
    }

    @Test
    public void updateDeltaMassTest() throws Exception {
        File testFile = new File("test-resources/test_psms.tsv");
        PSMFile psmFile = new PSMFile(testFile, 0);

        PSM psm = psmFile.psms.get(3);
        double prevCalcMass = Double.parseDouble(psm.spLine.get(psmFile.peptideCalcMassCol));
        psm.updateDeltaMass(2512, 0, 0, psmFile.peptideCalcMassCol, psmFile.calcMZcol, psmFile.dMassCol, psmFile.assignedModCol);
        assert psm.getDMass() - 0.8503 < tol;
        assert Double.parseDouble(psm.spLine.get(psmFile.peptideCalcMassCol)) - (prevCalcMass + 2512) < tol;

    }
}
