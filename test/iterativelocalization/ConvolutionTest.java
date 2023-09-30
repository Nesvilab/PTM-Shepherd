package iterativelocalization;

import edu.umich.andykong.ptmshepherd.iterativelocalization.Convolution;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class ConvolutionTest {

    @Test
    void convolve() {
        // Single convolution method
        double[] rv1 = {0.5, 0.5};
        double[] rv2 = {0.5, 0.5};
        double[] expected = {0.25, 0.50, 0.25};
        assertArrayEquals(expected, Convolution.convolve(rv1, rv2), 0.001);

        // Empty X
        rv1 = new double[]{};
        rv2 = new double[]{0.5, 0.5};
        expected = new double[]{0.5, 0.5};
        assertArrayEquals(expected, Convolution.convolve(rv1, rv2));

        // Empty Y
        rv1 = new double[]{0.1, 0.9};
        rv2 = new double[]{};
        expected = new double[]{0.1, 0.9};
        assertArrayEquals(expected, Convolution.convolve(rv1, rv2));

        // Empty X and Empty Y
        rv1 = new double[]{};
        rv2 = new double[]{};
        expected = new double[]{1.0};
        assertArrayEquals(expected, Convolution.convolve(rv1, rv2));

        // Multiple convolution method
        double[][] randomVariables = {
                {0.5, 0.5},
                {0.5, 0.5},
                {0.5, 0.5}
        };
        expected = new double[]{0.125, 0.375, 0.375, 0.125};
        assertArrayEquals(expected, Convolution.convolve(randomVariables), 0.001);

        // Multiple convolution method with empty variable
        randomVariables = new double[][]{
                {},
                {0.5, 0.5},
                {0.5, 0.5}
        };
        expected = new double[]{0.25, 0.5, 0.25};
        assertArrayEquals(expected, Convolution.convolve(randomVariables), 0.001);

        // Multiple convolution method with empty variable
        randomVariables = new double[][]{
                {0.5, 0.5},
                {0.5, 0.5},
                {}
        };
        expected = new double[]{0.25, 0.5, 0.25};
        assertArrayEquals(expected, Convolution.convolve(randomVariables), 0.001);
    }

    @Test
    void calculatePXGreaterThanY() {
        double[] pmfX = {0.5, 0.5};
        double[] pmfY = {0.5, 0.5};
        double expected = 0.25;
        assertEquals(expected, Convolution.calculatePXGreaterThanY(pmfX, pmfY), 0.001);

        pmfX = new double[]{};
        pmfY = new double[]{0.5, 0.5};
        expected = 0.5;
        assertEquals(expected, Convolution.calculatePXGreaterThanY(pmfX, pmfY), 0.001);

        pmfX = new double[]{0.5, 0.5};
        pmfY = new double[]{};
        expected = 0.5;
        assertEquals(expected, Convolution.calculatePXGreaterThanY(pmfX, pmfY), 0.001);

        pmfX = new double[]{};
        pmfY = new double[]{};
        expected = 0.0;
        assertEquals(expected, Convolution.calculatePXGreaterThanY(pmfX, pmfY), 0.001);
    }

    @Test
    void calculateXYWinProbability() {
        double[] pmfX = {0.5, 0.5};
        double[] pmfY = {0.5, 0.5};
        double[] expected = new double[]{0.25, 0.5, 0.25};
        assertArrayEquals(expected, Convolution.calculateXYWinProbability(pmfX, pmfY), 0.001);

        pmfX = new double[]{0.5, 0.5};
        pmfY = new double[]{};
        expected = new double[]{0.5, 0.5, 0.0};
        assertArrayEquals(expected, Convolution.calculateXYWinProbability(pmfX, pmfY), 0.001);

        pmfX = new double[]{};
        pmfY = new double[]{0.5, 0.5};
        expected = new double[]{0, 0.5, 0.5};
        assertArrayEquals(expected, Convolution.calculateXYWinProbability(pmfX, pmfY), 0.001);

        pmfX = new double[]{};
        pmfY = new double[]{};
        expected = new double[]{0, 1.0, 0.0};
        assertArrayEquals(expected, Convolution.calculateXYWinProbability(pmfX, pmfY), 0.001);

        pmfX = new double[]{0.0, 0.1666, 0.1666, 0.1666, 0.1666, 0.1666, 0.1666};
        pmfY = new double[]{0.0, 0.1666, 0.1666, 0.1666, 0.1666, 0.1666, 0.1666};
        expected = new double[]{0.4166, 0.1666, 0.4166};
        assertArrayEquals(expected, Convolution.calculateXYWinProbability(pmfX, pmfY), 0.001);
    }
}