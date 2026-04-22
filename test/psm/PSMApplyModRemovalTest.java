package psm;

import edu.umich.andykong.ptmshepherd.Mod;
import edu.umich.andykong.ptmshepherd.PSM;
import edu.umich.andykong.ptmshepherd.PSMFile;
import org.junit.jupiter.api.Test;

import java.io.File;

import static org.junit.jupiter.api.Assertions.*;

public class PSMApplyModRemovalTest {

    private static final double tol = 0.0001;

    /**
     * Test removing a single assigned mod (oxidation on Met) from a PSM with one mod.
     * Uses PSM index 2: peptide GVSHHGHGGMSGSHR with 10M(15.9949) and modified peptide GVSHHGHGGM[147]SGSHR
     */
    @Test
    public void testRemoveSingleMod() {
        File testFile = new File("test-resources/test_psms.tsv");
        PSMFile psmFile = new PSMFile(testFile, 0);

        PSM psm = psmFile.psms.get(2);
        assertEquals("GVSHHGHGGMSGSHR", psm.getPeptide());
        assertEquals(1, psm.getAssignedMods().size());
        assertEquals(10, psm.getAssignedMods().get(0).position);
        assertEquals(15.9949, psm.getAssignedMods().get(0).mass, tol);
        assertEquals("GVSHHGHGGM[147]SGSHR", psm.getModifiedPeptide());

        double prevDMass = psm.getDMass();
        double prevCalcPepMass = psm.getCalcPepMass();
        Mod modToRemove = psm.getAssignedMods().get(0);

        psm.applyModRemoval(modToRemove, psmFile.dMassCol, psmFile.assignedModCol,
                psmFile.modPeptideCol, psmFile.peptideCalcMassCol, psmFile.calcMZcol);

        // assignedMods should be empty after removal
        assertTrue(psm.getAssignedMods().isEmpty());

        // modified peptide should have the [147] bracket removed
        assertEquals("GVSHHGHGGMSGSHR", psm.getModifiedPeptide());

        // delta mass should increase by mod mass
        assertEquals(prevDMass + 15.9949, psm.getDMass(), tol);

        // calc pep mass should decrease by mod mass
        assertEquals(prevCalcPepMass - 15.9949, psm.getCalcPepMass(), tol);

        // spLine columns should be updated
        assertEquals(String.format("%.4f", psm.getDMass()), psm.spLine.get(psmFile.dMassCol));
        assertEquals(String.format("%.4f", psm.getCalcPepMass()), psm.spLine.get(psmFile.peptideCalcMassCol));
        assertEquals("", psm.spLine.get(psmFile.assignedModCol));  // no remaining mods
        assertEquals("GVSHHGHGGMSGSHR", psm.spLine.get(psmFile.modPeptideCol));
    }

    /**
     * Test removing a mod from a PSM that also has a variable mod shown in the modified peptide.
     * Uses PSM index 8: peptide HKDDCERMNITVKN with 8M(15.9949) and modified peptide HKDDCERM[147]NITVKN.
     * This is a realistic glyco scenario: the oxidized Met mod is incorrectly assigned and should be removed.
     */
    @Test
    public void testRemoveModWithBracketInModifiedPeptide() {
        File testFile = new File("test-resources/test_psms.tsv");
        PSMFile psmFile = new PSMFile(testFile, 0);

        PSM psm = psmFile.psms.get(8);
        assertEquals("HKDDCERMNITVKN", psm.getPeptide());
        assertEquals(1, psm.getAssignedMods().size());
        assertEquals(8, psm.getAssignedMods().get(0).position);
        assertEquals("HKDDCERM[147]NITVKN", psm.getModifiedPeptide());

        double prevDMass = psm.getDMass();
        double prevCalcPepMass = psm.getCalcPepMass();
        Mod modToRemove = psm.getAssignedMods().get(0);

        psm.applyModRemoval(modToRemove, psmFile.dMassCol, psmFile.assignedModCol,
                psmFile.modPeptideCol, psmFile.peptideCalcMassCol, psmFile.calcMZcol);

        // assignedMods should be empty
        assertTrue(psm.getAssignedMods().isEmpty());

        // modified peptide should no longer have [147]
        assertEquals("HKDDCERMNITVKN", psm.getModifiedPeptide());

        // mass adjustments
        assertEquals(prevDMass + 15.9949, psm.getDMass(), tol);
        assertEquals(prevCalcPepMass - 15.9949, psm.getCalcPepMass(), tol);
    }

    /**
     * Test that spLine columns (delta mass, calc pep mass, calc M/Z, assigned mods, modified peptide)
     * are all properly updated after mod removal.
     */
    @Test
    public void testSpLineColumnsUpdated() {
        File testFile = new File("test-resources/test_psms.tsv");
        PSMFile psmFile = new PSMFile(testFile, 0);

        PSM psm = psmFile.psms.get(2);  // charge 2, has 10M(15.9949)
        assertEquals(2, psm.getCharge());

        Mod modToRemove = psm.getAssignedMods().get(0);
        psm.applyModRemoval(modToRemove, psmFile.dMassCol, psmFile.assignedModCol,
                psmFile.modPeptideCol, psmFile.peptideCalcMassCol, psmFile.calcMZcol);

        // verify all spLine columns match the PSM state
        assertEquals(String.format("%.4f", psm.getDMass()), psm.spLine.get(psmFile.dMassCol));
        assertEquals(String.format("%.4f", psm.getCalcPepMass()), psm.spLine.get(psmFile.peptideCalcMassCol));
        assertEquals("", psm.spLine.get(psmFile.assignedModCol));
        assertEquals("GVSHHGHGGMSGSHR", psm.spLine.get(psmFile.modPeptideCol));

        // verify calc M/Z is consistent with calc pep mass and charge
        double expectedMZ = (psm.getCalcPepMass() + psm.getCharge() * 1.00728) / psm.getCharge();
        double actualMZ = Double.parseDouble(psm.spLine.get(psmFile.calcMZcol));
        assertEquals(expectedMZ, actualMZ, 0.01);  // wider tolerance for float rounding in neutralMassToMZ
    }

    /**
     * Test that removing a mod with no bracket in the modified peptide leaves modified peptide unchanged.
     * This happens for fixed mods (e.g., Cys carbamidomethylation) which are in Assigned Modifications
     * but not shown as brackets in Modified Peptide.
     */
    @Test
    public void testRemoveModWithoutBracketLeavesModifiedPeptideUnchanged() {
        File testFile = new File("test-resources/test_psms.tsv");
        PSMFile psmFile = new PSMFile(testFile, 0);

        // PSM index 5: KSCCTSEKPSCCSNGK with fixed Cys mods (no brackets in modified peptide)
        PSM psm = psmFile.psms.get(5);
        assertEquals("KSCCTSEKPSCCSNGK", psm.getPeptide());
        assertEquals("KSCCTSEKPSCCSNGK", psm.getModifiedPeptide());  // no brackets for fixed mods
        assertEquals(4, psm.getAssignedMods().size());

        double prevDMass = psm.getDMass();
        double prevCalcPepMass = psm.getCalcPepMass();

        // Remove one fixed Cys mod (position 3, mass 57.0214)
        Mod fixedMod = null;
        for (Mod mod : psm.getAssignedMods()) {
            if (mod.position == 3) {
                fixedMod = mod;
                break;
            }
        }
        assertNotNull(fixedMod);

        psm.applyModRemoval(fixedMod, psmFile.dMassCol, psmFile.assignedModCol,
                psmFile.modPeptideCol, psmFile.peptideCalcMassCol, psmFile.calcMZcol);

        // mod removed from assignedMods list
        assertEquals(3, psm.getAssignedMods().size());
        for (Mod mod : psm.getAssignedMods()) {
            assertFalse(mod.position == 3 && Math.abs(mod.mass - 57.0214) < tol,
                    "Removed mod should not be in assignedMods");
        }

        // modified peptide should be unchanged (no bracket was present to remove)
        assertEquals("KSCCTSEKPSCCSNGK", psm.getModifiedPeptide());

        // masses still adjusted
        assertEquals(prevDMass + 57.0214, psm.getDMass(), 0.001);
        assertEquals(prevCalcPepMass - 57.0214, psm.getCalcPepMass(), 0.001);
    }

    /**
     * Test removing a mod using the remove-delta test file (massdiffToVarmod=1).
     * PSM index 3 has the delta mass placed as a variable mod at N3: NFN[2627]DSSTK.
     * After massdiffToVarmod processing, the mod was removed from assignedMods and modified peptide.
     * We test removing the original assigned mod from the re-parsed state.
     */
    @Test
    public void testRemoveModAfterMassdiffToVarmod() {
        File testFile = new File("test-resources/test_psms_remove-delta.tsv");
        PSMFile psmFile = new PSMFile(testFile, 1);

        // PSM index 6: N-term acetylation. peptide AKPAQGAK with N-term(42.0106), modified peptide n[43]AKPAQGAK
        PSM psm = psmFile.psms.get(6);
        assertEquals("AKPAQGAK", psm.getPeptide());
        assertEquals(1, psm.getAssignedMods().size());
        assertEquals(0, psm.getAssignedMods().get(0).position);  // N-term
        assertEquals(42.0106, psm.getAssignedMods().get(0).mass, tol);

        double prevDMass = psm.getDMass();
        double prevCalcPepMass = psm.getCalcPepMass();
        Mod modToRemove = psm.getAssignedMods().get(0);

        // N-term removal is a no-op for modified peptide (returns unchanged)
        psm.applyModRemoval(modToRemove, psmFile.dMassCol, psmFile.assignedModCol,
                psmFile.modPeptideCol, psmFile.peptideCalcMassCol, psmFile.calcMZcol);

        // mod removed from assignedMods
        assertTrue(psm.getAssignedMods().isEmpty());

        // masses adjusted
        assertEquals(prevDMass + 42.0106, psm.getDMass(), tol);
        assertEquals(prevCalcPepMass - 42.0106, psm.getCalcPepMass(), tol);
    }
}
