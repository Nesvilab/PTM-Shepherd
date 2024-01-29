package iterativelocalization;

import edu.umich.andykong.ptmshepherd.iterativelocalization.PoissonBinomialLikelihood;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class PoissonBinomialLikelihoodTest {

    @Test
    void calculateProbXGreaterThanY() {
        // Standard example
        double[] xIonProbs = {-0.1, 0.5, 0.5, -0.9, -0.7};
        double[] yIonProbs = {0.1, -0.5, 0.5, -0.9, -0.7};
        double expected = 0.7;
        assertEquals(expected, PoissonBinomialLikelihood.calculateProbXGreaterThanY(xIonProbs, yIonProbs));

        // Only matched ions for X
        xIonProbs = new double[]{0.5, 0.5};
        yIonProbs = new double[]{-0.5, -0.5};
        expected = 0.875;
        assertEquals(expected, PoissonBinomialLikelihood.calculateProbXGreaterThanY(xIonProbs, yIonProbs));

        // Only matched ions for Y
        xIonProbs = new double[]{-0.5, -0.5};
        yIonProbs = new double[]{0.5, 0.5};
        expected = 0.125;
        assertEquals(expected, PoissonBinomialLikelihood.calculateProbXGreaterThanY(xIonProbs, yIonProbs));

        // No matched ions for X or Y
        xIonProbs = new double[]{-0.5, -0.5};
        yIonProbs = new double[]{-0.5, -0.5};
        expected = 0.5;
        assertEquals(expected, PoissonBinomialLikelihood.calculateProbXGreaterThanY(xIonProbs, yIonProbs));

        // Future-proofing for making unmatched ions 0 instead of negative
        xIonProbs = new double[]{0.5, 0.5};
        yIonProbs = new double[]{0.0, 0.0};
        expected = 0.875;
        assertEquals(expected, PoissonBinomialLikelihood.calculateProbXGreaterThanY(xIonProbs, yIonProbs));
    }

    @Test
    void calculateProbXMax() {
        // Standard example
        // Two sites, one matches two ions with probability 0.5 and 0.5 of being true, the second matches two ions
        // with probability 0.5 and 0.1 of being true. 0.5 ions are the same ion, so they cancel, leaving only a
        // comparison between the 0.5 and the 0.1 ions.
        double[][] ionProbs = new double[][]{
                {-0.1, 0.5, 0.5, -0.9, -0.7},
                {0.1, -0.5, 0.5, -0.9, -0.7},
        };
        double[] expected = new double[]{0.7, 0.3};
        assertArrayEquals(expected, PoissonBinomialLikelihood.calculateProbXMax(ionProbs), 0.001);


        // Standard example but with no unique ions in row 1
        ionProbs = new double[][]{
                {-0.1, 0.5, 0.5, -0.9, -0.7, -0.5},
                {0.1, -0.5, 0.5, -0.9, -0.7, -0.5},
                {-0.1, 0.5, 0.5, -0.9, -0.7, 0.5},
        };
        expected = new double[]{0.2244, 0.0981, 0.6775};
        assertArrayEquals(expected, PoissonBinomialLikelihood.calculateProbXMax(ionProbs), 0.001);

        // Standard example but with no unique ions in row 3
        ionProbs = new double[][]{
                {0.1, -0.5, 0.5, -0.9, -0.7, -0.5},
                {-0.1, 0.5, 0.5, -0.9, -0.7, 0.5},
                {-0.1, 0.5, 0.5, -0.9, -0.7, -0.5}
        };
        expected = new double[]{0.0981, 0.6775, 0.2244};
        assertArrayEquals(expected, PoissonBinomialLikelihood.calculateProbXMax(ionProbs), 0.001);

        // Literally only wrote this for the code coverage
        PoissonBinomialLikelihood.printMatrix(ionProbs);
    }
}
