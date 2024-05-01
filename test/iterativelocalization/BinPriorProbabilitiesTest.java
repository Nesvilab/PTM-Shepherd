package iterativelocalization;

import edu.umich.andykong.ptmshepherd.iterativelocalization.BinPriorProbabilities;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class BinPriorProbabilitiesTest {
    BinPriorProbabilities bpp;

    @BeforeEach
    public void setUp() {
        // Initialize test data
        this.bpp = new BinPriorProbabilities();
    }

    @Test
    public void computePriorProbs() {
        String seq = "PEPTI";
        boolean[] allowedPoses = new boolean[]{true, true, true, true, true};

        double[] expectedPriorProbs = new double[]{0.2, 0.2, 0.2, 0.2, 0.2};
        assertArrayEquals(bpp.computePriorProbs(seq, allowedPoses), expectedPriorProbs, 0.0001);
    }

    @Test
    public void computeUniformPriorProbs() {
        String seq = "PEPTI";
        boolean[] allowedPoses = new boolean[]{true, true, true, true, true};

        double[] expectedPriorProbs = new double[]{0.2, 0.2, 0.2, 0.2, 0.2};
        assertArrayEquals(bpp.computeUniformPriorProbs(seq, allowedPoses), expectedPriorProbs, 0.0001);
    }
}
