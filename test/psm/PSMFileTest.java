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
        float prevCalcMass = psm.getCalcPepMass();
        float prevDMass = psm.getDMass();

        psm.updateDeltaMass(2512.8455f, 0, 0, psmFile.peptideCalcMassCol, psmFile.calcMZcol, psmFile.dMassCol, psmFile.assignedModCol);
        assert psm.getDMass() - 0.8503 < tol;
        assert psm.getDMass() - prevDMass < tol;
        assert psm.getCalcPepMass() - (prevCalcMass + 2512.8455f) < tol;
        assert Float.parseFloat(psm.spLine.get(psmFile.peptideCalcMassCol)) - (prevCalcMass + 2512.8455) < tol;

        psm.updateDeltaMass(2512.8455f, 1, 2512.8455f, psmFile.peptideCalcMassCol, psmFile.calcMZcol, psmFile.dMassCol, psmFile.assignedModCol);
        assert psm.getDMass() - 0.8503 < tol;
        assert psm.getCalcPepMass() - (prevCalcMass + 2512.8455f) < tol;
    }

    @Test
    public void editModifiedPeptideTest0() throws Exception {
        File testFile = new File("test-resources/test_psms.tsv");
        PSMFile psmFile = new PSMFile(testFile, 0);

        PSM psm = psmFile.psms.get(3);
        assert psm.getOriginalModifiedPeptide().equals("NFNDSSTK");
        assert psm.getModifiedPeptide().equals("NFNDSSTK");
        assert psm.spLine.get(psmFile.modPeptideCol).isEmpty();

        psm.editModifiedPeptide(3, 2512.8455, 0, psmFile.modPeptideCol);
        assert psm.getModifiedPeptide().equals("NFN[2627]DSSTK");
        assert psm.spLine.get(psmFile.modPeptideCol).equals("NFN[2627]DSSTK");

        // test re-run of the same file with previous mod
        psm.editModifiedPeptide(3, 2512.8455, 1, psmFile.modPeptideCol);
        assert psm.getModifiedPeptide().equals("NFN[2627]DSSTK");
        assert psm.spLine.get(psmFile.modPeptideCol).equals("NFN[2627]DSSTK");

        // test re-run of the same file with a new mod
        psm.editModifiedPeptide(3, 2000, 1, psmFile.modPeptideCol);
        assert psm.getModifiedPeptide().equals("NFN[2115]DSSTK");
        assert psm.spLine.get(psmFile.modPeptideCol).equals("NFN[2115]DSSTK");
    }

    @Test
    public void editModifiedPeptideTest1() throws Exception {
        File testFile = new File("test-resources/test_psms_remove-delta.tsv");
        PSMFile psmFile = new PSMFile(testFile, 1);

        PSM psm = psmFile.psms.get(3);
        assert psm.getOriginalModifiedPeptide().equals("NFN[2627]DSSTK");
        assert psm.getModifiedPeptide().equals("NFNDSSTK");
        assert psm.spLine.get(psmFile.modPeptideCol).equals("NFN[2627]DSSTK");

        psm.editModifiedPeptide(3, 2512.8455, 1, psmFile.modPeptideCol);
        assert psm.getModifiedPeptide().equals("NFN[2627]DSSTK");
        assert psm.spLine.get(psmFile.modPeptideCol).equals("NFN[2627]DSSTK");
    }

    @Test
    public void editModifiedPeptideTest2() throws Exception {
        File testFile = new File("test-resources/test_psms_keep-delta.tsv");
        PSMFile psmFile = new PSMFile(testFile, 2);

        PSM psm = psmFile.psms.get(3);
        assert psm.getOriginalModifiedPeptide().equals("NFN[2627]DSSTK");
        assert psm.getModifiedPeptide().equals("NFNDSSTK");
        assert psm.spLine.get(psmFile.modPeptideCol).equals("NFN[2627]DSSTK");

        psm.editModifiedPeptide(3, 2512.8455, 1, psmFile.modPeptideCol);
        assert psm.getModifiedPeptide().equals("NFN[2627]DSSTK");
        assert psm.spLine.get(psmFile.modPeptideCol).equals("NFN[2627]DSSTK");
    }
}
