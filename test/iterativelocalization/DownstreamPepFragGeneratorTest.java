package iterativelocalization;

import edu.umich.andykong.ptmshepherd.iterativelocalization.DownstreamPepFragGenerator;
import edu.umich.andykong.ptmshepherd.utils.Peptide;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class DownstreamPepFragGeneratorTest {

    @Test
    void calculatePeptideFragments() {
        String seq = "PEPT";
        float[] mods = new float[]{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        Peptide pep = new Peptide(seq, mods);

        // Only b2, b3
        ArrayList<Float> pepFrags = DownstreamPepFragGenerator.calculatePeptideFragments(pep, "b", 0);
        ArrayList<Float> expectedFrags = new ArrayList<>(Arrays.asList(227.1027f, 324.1554f));
        assertEquals(expectedFrags.size(), pepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++)
            assertEquals(expectedFrags.get(i), pepFrags.get(i), 0.0001, "Mismatch at index " + i);

        // Only b2, b3, so empty
        pepFrags = DownstreamPepFragGenerator.calculatePeptideFragments(pep, "y", 0);
        expectedFrags = new ArrayList<>(Arrays.asList());
        assertEquals(expectedFrags.size(), pepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++)
            assertEquals(expectedFrags.get(i), pepFrags.get(i), 0.0001, "Mismatch at index " + i);


        // y1, y2, y3, so empty
        pepFrags = DownstreamPepFragGenerator.calculatePeptideFragments(pep, "b", 3);
        expectedFrags = new ArrayList<>(Arrays.asList());
        assertEquals(expectedFrags.size(), pepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++)
            assertEquals(expectedFrags.get(i), pepFrags.get(i), 0.0001, "Mismatch at index " + i);

        // y1, y2, y3
        pepFrags = DownstreamPepFragGenerator.calculatePeptideFragments(pep, "y", 3);
        expectedFrags = new ArrayList<>(Arrays.asList(120.06556f, 217.11833f, 346.16092f));
        assertEquals(expectedFrags.size(), pepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++)
            assertEquals(expectedFrags.get(i), pepFrags.get(i), 0.0001, "Mismatch at index " + i);


        // b3, y2, y3
        pepFrags = DownstreamPepFragGenerator.calculatePeptideFragments(pep, "by", 2);
        expectedFrags = new ArrayList<>(Arrays.asList(324.15544f, 217.11833f, 346.16092f));
        assertEquals(expectedFrags.size(), pepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++)
            assertEquals(expectedFrags.get(i), pepFrags.get(i), 0.0001, "Mismatch at index " + i);

        // b2, b3, y3
        pepFrags = DownstreamPepFragGenerator.calculatePeptideFragments(pep, "by", 1);
        expectedFrags = new ArrayList<>(Arrays.asList(227.1027f, 324.1554f, 346.16092f));
        assertEquals(expectedFrags.size(), pepFrags.size());
        for (int i = 0; i < expectedFrags.size(); i++)
            assertEquals(expectedFrags.get(i), pepFrags.get(i), 0.0001, "Mismatch at index " + i);


    }
}
