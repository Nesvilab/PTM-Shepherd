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
        PSMFile psmFile = new PSMFile(testFile, 0);

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
        PSMFile psmFile = new PSMFile(testFile, 1);

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
        PSMFile psmFile = new PSMFile(testFile, 2);

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
    public void updateDeltaMassTest() {
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
    // test editing modified peptide for mass-diff-to-varmod = 0
    public void editModifiedPeptideTest0() {
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
        assert psm.getModifiedPeptide().equals("NFN[2114]DSSTK");
        assert psm.spLine.get(psmFile.modPeptideCol).equals("NFN[2114]DSSTK");

        // test another PSM with a different mod location and subsequent AA
        PSM psm2 = psmFile.psms.get(4);
        assert psm2.getOriginalModifiedPeptide().equals("SNATKPQCPK");
        assert psm2.getModifiedPeptide().equals("SNATKPQCPK");

        psm2.editModifiedPeptide(2, 1864.6342, 0, psmFile.modPeptideCol);
        assert psm2.getModifiedPeptide().equals("SN[1979]ATKPQCPK");
        assert psm2.spLine.get(psmFile.modPeptideCol).equals("SN[1979]ATKPQCPK");

        // test editing a PSM with a glycan modification at the end of the peptide
        PSM psm3 = psmFile.psms.get(7);
        assert psm3.getOriginalModifiedPeptide().equals("KETLHKQYHLVKSHTN");
        assert psm3.getModifiedPeptide().equals("KETLHKQYHLVKSHTN");
        assert psm3.spLine.get(psmFile.modPeptideCol).isEmpty();

        psm3.editModifiedPeptide(16, 203.0794, 0, psmFile.modPeptideCol);
        assert psm3.getModifiedPeptide().equals("KETLHKQYHLVKSHTN[317]");
        assert psm3.spLine.get(psmFile.modPeptideCol).equals("KETLHKQYHLVKSHTN[317]");

        // test adding a glycan to a PSM with a variable modification
        PSM psm4 = psmFile.psms.get(8);
        assert psm4.getOriginalModifiedPeptide().equals("HKDDCERM[147]NITVKN");
        assert psm4.getModifiedPeptide().equals("HKDDCERM[147]NITVKN");
        assert psm4.spLine.get(psmFile.modPeptideCol).equals("HKDDCERM[147]NITVKN");

        psm4.editModifiedPeptide(9, 1038.3751, 0, psmFile.modPeptideCol);
        assert psm4.getModifiedPeptide().equals("HKDDCERM[147]N[1152]ITVKN");
        assert psm4.getOriginalModifiedPeptide().equals("HKDDCERM[147]NITVKN");
        assert psm4.spLine.get(psmFile.modPeptideCol).equals("HKDDCERM[147]N[1152]ITVKN");
    }

    @Test
    // test editing modified peptide for mass-diff-to-varmod = 1
    public void editModifiedPeptideTest1() {
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
    // test editing modified peptide for mass-diff-to-varmod = 2
    public void editModifiedPeptideTest2() {
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
