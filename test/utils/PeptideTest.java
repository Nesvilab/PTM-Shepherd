package utils;

import edu.umich.andykong.ptmshepherd.utils.Peptide;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

public class PeptideTest {

    @Test
    void generateDecoy() {
        String seq = "PEPTIDE";
        float[] mods = new float[]{1.0f, 0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 3.0f};
        Random rng = new Random(1);
        String expectedSeq = "PTIPDEE";
        float[] expectedMods = new float[]{1.0f, 2.0f, 0.0f, 0.0f, 0.0f ,0.0f, 3.0f};

        Peptide decoy =  Peptide.generateDecoy(seq, mods, rng, "shuffled");
        assertEquals(expectedSeq, decoy.pepSeq);
        assertArrayEquals(expectedMods, decoy.mods);

        decoy =  Peptide.generateDecoy(seq, mods, rng, "mutated");
        expectedSeq = "PDPTIEE";
        assertEquals(expectedSeq, decoy.pepSeq);
    }
}
