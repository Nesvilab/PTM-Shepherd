package psm;

import edu.umich.andykong.ptmshepherd.Mod;
import edu.umich.andykong.ptmshepherd.PSM;
import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.glyco.GlycoParams;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

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

        // Cys disulfide special handling
        PSM psmCysDisulfide = psmFile.psms.get(10);
        assert psmCysDisulfide.getScanNum() == 8213;
        assert Math.abs(psmCysDisulfide.getDMass() + 2.01565 - psmCysDisulfide.getOriginalDeltaMass()) < tol;
        for (int i=0; i < psmCysDisulfide.getOriginalAssignedMods().size(); i++) {
            assert psmCysDisulfide.getAssignedMods().get(i).equals(psmCysDisulfide.getOriginalAssignedMods().get(i));
        }
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

        // Cys disulfide special handling
        PSM psmCysDisulfide = psmFile.psms.get(6);
        assert psmCysDisulfide.getScanNum() == 8213;
        assert Math.abs(psmCysDisulfide.getDMass() - psmCysDisulfide.getOriginalDeltaMass() - (2 * 57.02146)) < tol;
        for (int i=0; i < psmCysDisulfide.getOriginalAssignedMods().size(); i++) {
            assert psmCysDisulfide.getAssignedMods().get(i).equals(psmCysDisulfide.getOriginalAssignedMods().get(i));
        }
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
    // test that mergeGlycoTable adds exactly 3 columns to every PSM, including PSMs sharing a scan number
    public void mergeGlycoTableMultiRankTest() throws IOException {
        // Copy to a temp file so mergeGlycoTable's save() doesn't overwrite the test resource
        Path tempFile = Files.createTempFile("test_psms_multirank", ".tsv");
        try {
            Files.copy(Paths.get("test-resources/test_psms.tsv"), tempFile, StandardCopyOption.REPLACE_EXISTING);
            PSMFile psmFile = new PSMFile(tempFile.toFile(), 0);

            int initialHeaderCount = psmFile.headers.length;
            // PSMs at indices 9, 10, 11 all share the same spectrum (scan 7005) — the multi-rank case
            assert psmFile.psms.get(9).getScanNum() == 7005;
            assert psmFile.psms.get(10).getScanNum() == 7005;
            assert psmFile.psms.get(11).getScanNum() == 7005;

            // Record each PSM's column count before the call (TSV rows may have trailing fields stripped)
            int[] initialSpLineSizes = new int[psmFile.psms.size()];
            for (int i = 0; i < psmFile.psms.size(); i++) {
                initialSpLineSizes[i] = psmFile.psms.get(i).spLine.size();
            }

            // Leave all psm.glycanAssignmentResult null (no glycan found) — simplest valid input
            GlycoParams glycoParams = new GlycoParams("", "", "");
            psmFile.mergeGlycoTable("test_dataset", glycoParams);

            // Headers should gain exactly 3 glyco columns
            assert psmFile.headers.length == initialHeaderCount + 3 :
                    String.format("Expected %d headers, got %d", initialHeaderCount + 3, psmFile.headers.length);

            // Every PSM must have gained exactly 3 columns — including the multi-rank scan 7005 PSMs at
            // indices 9, 10, 11 which share the same spectrum value and previously triggered the bug where
            // rank-1 got 0 extra columns and rank-2 got 6 extra columns via the scanToLineMap overwrite.
            for (int i = 0; i < psmFile.psms.size(); i++) {
                PSM psm = psmFile.psms.get(i);
                int gained = psm.spLine.size() - initialSpLineSizes[i];
                assert gained == 3 :
                        String.format("PSM lineNum=%d gained %d columns, expected 3", psm.lineNum, gained);
            }
        } finally {
            Files.deleteIfExists(tempFile);
        }
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
