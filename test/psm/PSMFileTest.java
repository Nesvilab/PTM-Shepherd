package psm;

import edu.umich.andykong.ptmshepherd.PSM;
import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.glyco.GlycoParams;
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
        assert Math.abs(psm1.getDMass() + 0.0045) < tol;
        assert psm1.getAssignedMods().isEmpty();

        PSM psm4 = psmFile.psms.get(3);
        assert psm4.getScanNum() == 2903;
        assert Math.abs(psm4.getDMass() - 0.0048 - 2512.8455) < tol;
        assert Math.abs(psm4.getOriginalDeltaMass() - 0.0048) < tol;
        assert psm4.getAssignedMods().isEmpty();
        assert psm4.getOriginalAssignedMods().get(0).mass - 2512.8455 < tol;

        // N-term modification from varmod
        PSM psmNvarmod = psmFile.psms.get(6);
        assert psmNvarmod.getScanNum() == 7005;
        assert Math.abs(psmNvarmod.getDMass() + 0.0007) < tol;
        assert Math.abs(psmNvarmod.getAssignedMods().get(0).mass - 42.0106) < tol;
        assert Math.abs(psmNvarmod.getOriginalAssignedMods().get(0).mass - 42.0106) < tol;

        // N-term modification from offset
        PSM psmNoffset = psmFile.psms.get(7);
        assert psmNoffset.getScanNum() == 7516;
        assert Math.abs(psmNoffset.getDMass() - 0.0007 - 27.9949) < tol;
        assert Math.abs(psmNoffset.getOriginalDeltaMass() - 0.0007) < tol;
        assert psmNoffset.getAssignedMods().isEmpty();
        assert Math.abs(psmNoffset.getOriginalAssignedMods().get(0).mass - 27.9949) < tol;

        // C-term modification from offset
        PSM psmCoffset = psmFile.psms.get(8);
        assert psmCoffset.getScanNum() == 118089;
        assert Math.abs(psmCoffset.getDMass() - 0.0059 + 0.98401) < tol;
        assert Math.abs(psmCoffset.getOriginalDeltaMass() - 0.0059) < tol;
        assert psmCoffset.getAssignedMods().isEmpty();
        assert Math.abs(psmCoffset.getOriginalAssignedMods().get(0).mass + 0.98401) < tol;

        PSM psmFixedAndOffsetCys = psmFile.psms.get(9);
        assert psmFixedAndOffsetCys.getScanNum() == 38977;
        assert Math.abs(psmFixedAndOffsetCys.getDMass() - 57.02146 + 9.0367 + 0.0052) < tol;
        assert Math.abs(psmFixedAndOffsetCys.getOriginalDeltaMass() + 0.0052) < tol;
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
        assert Math.abs(psm4.getDMass() - 0.0048 - 2512.8455) < tol;
        assert Math.abs(psm4.getOriginalDeltaMass() - 0.0048 - 2512.8455) < tol;
        assert psm4.getAssignedMods().isEmpty();
        assert Math.abs(psm4.getOriginalAssignedMods().get(0).mass - 2512.8455) < tol;
    }

    @Test
    public void updateDeltaMassTest() {
        File testFile = new File("test-resources/test_psms.tsv");
        PSMFile psmFile = new PSMFile(testFile, 0);

        PSM psm = psmFile.psms.get(3);
        float prevCalcMass = psm.getCalcPepMass();
        float prevDMass = psm.getDMass();
        GlycoParams params = new GlycoParams("", "", "");
        params.removeGlycanDeltaMass = true;
        params.writeGlycansToAssignedMods = true;

        psmFile.writeGlycanToAssignedMod(psm, "HexNAc(2)Hex(13)", false, params);
        assert Math.abs(psm.getDMass() - 0.0049) < tol;
        assert Math.abs(psm.getOriginalDeltaMass() - prevDMass) < tol;
        assert Math.abs(psm.getCalcPepMass() - (prevCalcMass + 2512.8455f)) < tol;
        assert Math.abs(Float.parseFloat(psm.spLine.get(psmFile.peptideCalcMassCol)) - (prevCalcMass + 2512.8455)) < tol;

        // simulate re-run of the same file with previous mod
        psm.spLine.set(psmFile.peptideCalcMassCol, String.format("%.4f", psm.getCalcPepMass()));    // reset calc mass column as if reading from new psm file
        psm.initializeMods(1, psmFile.dMassCol, psmFile.assignedModCol, psmFile.modPeptideCol, psmFile.msfraggerLocalizationCol, psmFile.peptideCalcMassCol);
        psmFile.massdiffToVarmod = 1;
        psmFile.writeGlycanToAssignedMod(psm, "HexNAc(2)Hex(13)", false, params);
        assert Math.abs(psm.getDMass() - 0.0049) < tol;
        assert Math.abs(psm.getCalcPepMass() - (prevCalcMass + 2512.8455f)) < tol;
        assert psm.getAssignedMods().size() == 1;
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

        PSM nVarmodPsm = psmFile.psms.get(6);
        nVarmodPsm.editModifiedPeptide(0, 42.0106, 1, psmFile.modPeptideCol);
        assert nVarmodPsm.getOriginalModifiedPeptide().equals("n[43]AKPAQGAK");
        assert nVarmodPsm.getModifiedPeptide().equals("n[43]AKPAQGAK");
        assert nVarmodPsm.spLine.get(psmFile.modPeptideCol).equals("n[43]AKPAQGAK");

        PSM nOffsetPsm = psmFile.psms.get(7);
        nOffsetPsm.editModifiedPeptide(0, 27.9949, 1, psmFile.modPeptideCol);
        assert nOffsetPsm.getOriginalModifiedPeptide().equals("n[29]AKHHAISAK");
        assert nOffsetPsm.getModifiedPeptide().equals("n[29]AKHHAISAK");
        assert nOffsetPsm.spLine.get(psmFile.modPeptideCol).equals("n[29]AKHHAISAK");

        PSM cOffsetPsm = psmFile.psms.get(8);
        cOffsetPsm.editModifiedPeptide(19, -0.98401, 1, psmFile.modPeptideCol);
        assert cOffsetPsm.getOriginalModifiedPeptide().equals("ILTEAEIDAHLVALAERDc[17]");
        assert cOffsetPsm.getModifiedPeptide().equals("ILTEAEIDAHLVALAERDc[17]");
        assert cOffsetPsm.spLine.get(psmFile.modPeptideCol).equals("ILTEAEIDAHLVALAERDc[17]");
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
