package edu.umich.andykong.ptmshepherd.iterativelocalization;

/**
 * This class contains static functions for computing convolutions of two random categorical variables.
 */
public class Convolution {
     /**
     * This function can compute the convolution of >=2 random variables by sequentially convoluting them. This is
     * possible because convolutions are commutative. The input is an array of size i x j, with i variables having
     * a probability p that outcome j is the result.  The result is the probability mass function (PMF)
     * of the sum of outcomes.
     *
     * Example inputs would be:
     * coin : {
     *      {0.5, 0.5}
     *      }
     * die : {
     *      {0.0, 0.166, 0.166, 0.166, 0.166, 0.166, 0.166} First index is chance of rolling a 0
     *      }
     * die + coin : {
     *      {0.0, 0.166, 0.166, 0.166, 0.166, 0.166, 0.166},
     *      {0.5, 0.5}
     *      }
     *
     * @param randomVariables   matrix of random variables
     * @return PMF of sum of random variable outcomes
     */
    public static double[] convolve(double[][] randomVariables) {
        double[] result = randomVariables[0];  // Start with the first random variable

        for (int i = 1; i < randomVariables.length; i++) {
            result = convolve(result, randomVariables[i]);  // Convolve with each subsequent random variable
        }

        return result;
    }

    /**
     * Effective inner function for convolve(double[][]), performing one of n-1 convolutions. Can be called
     * separately as well.
     * @param rv1   PMF for one random variable in format {0.5, 0.5}
     * @param rv2   PMF for another random variable in format {0.5, 0.5}
     * @return PMF of sum of rv1 and rv2
     */
    public static double[] convolve(double[] rv1, double[] rv2) {
        // Handle cases where there are missing values, make them 100% chance of equal to 0
        if (rv1.length == 0)
            rv1 = new double[]{1.0};
        if (rv2.length == 0)
            rv2 = new double[]{1.0};

        int n = rv1.length + rv2.length - 1;
        double[] result = new double[n];

        for (int i = 0; i < rv1.length; i++) {
            for (int j = 0; j < rv2.length; j++) {
                result[i + j] += rv1[i] * rv2[j];
            }
        }

        return result;
    }

    public static double calculatePXGreaterThanY(double[] pmfX, double[] pmfY) {
        int lenX = pmfX.length;
        int lenY = pmfY.length;

        // Deal with cases where one is empty
        if (pmfX.length == 0 && pmfY.length == 0)
            return 0.0;
        if (pmfX.length == 0)
            return 1.0 - pmfY[0];
        if (pmfY.length == 0)
            return 1.0 - pmfX[0];

        // Reverse and pad pmfY
        double[] pmfNegY = new double[lenY * 2 - 1];
        for (int i = 0; i < lenY; i++) {
            pmfNegY[i] = pmfY[lenY - i - 1];
        }

        // Compute the convolution of pmfX and pmfNegY
        double[] pmfZ = new double[lenX + pmfNegY.length - 1];
        for (int i = 0; i < pmfZ.length; i++) {
            for (int j = Math.max(0, i + 1 - pmfNegY.length); j < Math.min(lenX, i + 1); j++) {
                pmfZ[i] += pmfX[j] * pmfNegY[i - j];
            }
        }

        // Index corresponding to z = 0
        int zeroIndex = lenY - 1;

        // Sum the probabilities for z > 0
        double pXGreaterThanY = 0;
        for (int i = zeroIndex + 1; i < pmfZ.length; i++) {
            pXGreaterThanY += pmfZ[i];
        }

        return pXGreaterThanY;
    }

    /**
     * Computed the probably mass function (PMF) of Z, where Z = X - Y. This is the PMF of X > Y.
     * @param pmfX  PMF of X
     * @param pmfY  PMF of Y
     * @return length 3 array where [0] is X win, [1] is equal, and [2] is Y win
     */
    public static double[] calculateXYWinProbability(double[] pmfX, double[] pmfY) {
        double [] pXYWin = new double[3];

        // Deal with cases where there is an empty array
        if (pmfX.length == 0 && pmfY.length == 0) {
            pXYWin[0] = 0.0; // No chance of X win
            pXYWin[1] = 1.0; // They can only tie at 0 with probability 1
            pXYWin[2] = 0.0; // No chance of Y win
            return pXYWin;
        } else if (pmfX.length == 0) {
            pXYWin[0] = 0.0; // No chance of X win
            double pYWin = 0.0;
            for (int i = 1; i < pmfY.length; i++)
                pYWin += pmfY[i];
            pXYWin[1] = 1.0 - pYWin; // Chance of tie
            pXYWin[2] = pYWin; // Chcance of Y win
            return pXYWin;
        } else if (pmfY.length == 0) {
            pXYWin[2] = 0.0; // No chance of Y win
            double pXWin = 0.0;
            for (int i = 1; i < pmfX.length; i++)
                pXWin += pmfX[i];
            pXYWin[1] = 1.0 - pXWin; // Chance of tie
            pXYWin[0] = pXWin; // Chcance of X win
            return pXYWin;
        }

        // Reverse and pad pmfY
        double[] pmfNegY = new double[pmfY.length * 2 - 1];
        for (int i = 0; i < pmfY.length; i++)
            pmfNegY[i] = pmfY[pmfY.length - i - 1];

        // Compute the convolution of pmfX and pmfNegY
        double[] pmfZ = new double[pmfX.length + pmfNegY.length - 1];
        for (int i = 0; i < pmfZ.length; i++) {
            for (int j = Math.max(0, i + 1 - pmfNegY.length); j < Math.min(pmfX.length, i + 1); j++) {
                pmfZ[i] += pmfX[j] * pmfNegY[i - j];
            }
        }

        // Index corresponding to z = 0
        int zeroIndex = pmfY.length - 1;
        // Get probability that z = 0
        double pXEqualsY = pmfZ[zeroIndex];

        // Sum the probabilities for z > 0
        double pXGreaterThanY = 0;
        for (int i = zeroIndex + 1; i < pmfZ.length; i++)
            pXGreaterThanY += pmfZ[i];

        // Sum probabilities for z < 0
        double pYGreaterThanX = 0;
        for (int i = zeroIndex - 1; i >= 0; i--)
            pYGreaterThanX += pmfZ[i];

        pXYWin = new double[]{pXGreaterThanY, pXEqualsY, pYGreaterThanX};

        return pXYWin;
    }
}