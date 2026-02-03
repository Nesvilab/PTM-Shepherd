package utils;

import edu.umich.andykong.ptmshepherd.utils.Peptide;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

public class PeptideTest {

    @Test
    void generateDecoy() {
        String seq = "PEPTIDE";
        double[] mods = new double[]{1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0};
        Random rng = new Random(1);
        String expectedSeq = "PTIPDEE";
        double[] expectedMods = new double[]{1.0, 2.0, 0.0, 0.0, 0.0 ,0.0, 3.0};

        Peptide decoy =  Peptide.generateDecoy(seq, mods, rng, "shuffled");
        assertEquals(expectedSeq, decoy.pepSeq);
        assertArrayEquals(expectedMods, decoy.mods);

        decoy = Peptide.generateDecoy(seq, mods, rng, "mutated");
        expectedSeq = "PDPTIEE";
        assertEquals(expectedSeq, decoy.pepSeq);

        decoy = Peptide.generateDecoy(seq, mods, rng, "mono-mutated");
        expectedSeq = "PEPTIDD";
        assertEquals(expectedSeq, decoy.pepSeq);
        assertEquals(6, decoy.mutatedResidue);
    }

    @Test
    void calculatePeptideFragments() {
        String seq = "PEPT";
        double[] mods = new double[]{0.0, 0.0, 0.0, 0.0};
        Peptide pep = new Peptide(seq, mods);

        ArrayList<Float> sitePepFrags = pep.calculatePeptideFragments("b", 1);
        ArrayList<Float> expectedFrags = new ArrayList<>(Arrays.asList(227.1027f, 324.1554f));
        assertEquals(expectedFrags.size(), sitePepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++)
            assertEquals(expectedFrags.get(i), sitePepFrags.get(i), 0.0001, "Mismatch at index " + i);

        sitePepFrags = pep.calculatePeptideFragments("y", 1);
        expectedFrags = new ArrayList<>(Arrays.asList(120.0656f, 217.1183f, 346.1609f));
        assertEquals(expectedFrags.size(), sitePepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++) {
            assertEquals(expectedFrags.get(i), sitePepFrags.get(i), 0.0001, "Mismatch at index " + i);
        }

        seq = "GDRGEIGPPGPR";
        mods = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0004, 0.0, 0.0, 0.0, 0.0};
        pep = new Peptide(seq, mods);

        sitePepFrags = pep.calculatePeptideFragments("b", 1);

        System.out.println(sitePepFrags.toString());

    }

    @Test
    void calculatePeptideFragmentsBetween() {
        String seq = "PEPTIDE";
        double[] mods = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        Peptide pep = new Peptide(seq, mods);

        ArrayList<Float> expectedFrags = new ArrayList<>(Arrays.asList(324.1554f, 425.20306f, 376.17145f, 477.21912f));
        ArrayList<Float> sitePepFrags = pep.calculatePeptideFragmentsBetween("by", 2, 4, 1);

        assertEquals(expectedFrags.size(), sitePepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++) {
            assertEquals(expectedFrags.get(i), sitePepFrags.get(i), 0.0001, "Mismatch at index " + i);
        }

        expectedFrags = new ArrayList<>(Arrays.asList(227.102633f, 324.155397f));
        sitePepFrags = pep.calculatePeptideFragmentsBetween("b", 1, 3, 1);

        assertEquals(expectedFrags.size(), sitePepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++) {
            assertEquals(expectedFrags.get(i), sitePepFrags.get(i), 0.0001, "Mismatch at index " + i);
        }

        expectedFrags = new ArrayList<>(Arrays.asList(477.219120f, 574.271884f));
        sitePepFrags = pep.calculatePeptideFragmentsBetween("y", 1, 3, 1);

        assertEquals(expectedFrags.size(), sitePepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++) {
            assertEquals(expectedFrags.get(i), sitePepFrags.get(i), 0.0001, "Mismatch at index " + i);
        }
    }

    @Test
    void calculateComplementaryPeptideFragments() {
        String seq = "PEPT";
        double[] mods = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        Peptide pep = new Peptide(seq, mods);

        // Simulate 1.0f dmass on position P3
        ArrayList<Float> expectedFrags = new ArrayList<>(Arrays.asList(228.102633f, 324.155397f));
        ArrayList<Float> sitePepFrags = pep.calculateComplementaryFragments("b", 1.0f, 2, 1);

        assertEquals(expectedFrags.size(), sitePepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++) {
            assertEquals(expectedFrags.get(i), sitePepFrags.get(i), 0.0001, "Mismatch at index " + i);
        }

        // Simulate 1.0f dmass on position P3
        expectedFrags = new ArrayList<>(Arrays.asList(121.065520f, 217.118284f, 346.160877f));
        sitePepFrags = pep.calculateComplementaryFragments("y", 1.0f, 2, 1);

        assertEquals(expectedFrags.size(), sitePepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++) {
            assertEquals(expectedFrags.get(i), sitePepFrags.get(i), 0.0001, "Mismatch at index " + i);
        }

    }

    @Test
    void AddModTest() {
        String seq = "PEP";
        double[] mods = new double[]{1.0, 0.0, 0.0};
        double[] expectedMods = new double[]{0.0, 0.0, 0.0};
        Peptide pep =  new Peptide(seq, mods);

        pep.addMod(-1*1.0, 0);
        assertArrayEquals(pep.mods, expectedMods);

        pep.addMod(1.0, 0);
        expectedMods = new double[]{1.0, 0.0, 0.0};
        assertArrayEquals(pep.mods, expectedMods);
    }
}
