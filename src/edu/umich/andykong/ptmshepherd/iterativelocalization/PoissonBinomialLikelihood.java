package edu.umich.andykong.ptmshepherd.iterativelocalization;

import java.util.ArrayList;
import java.util.Arrays;

import static edu.umich.andykong.ptmshepherd.iterativelocalization.Convolution.calculateXYWinProbability;
import static edu.umich.andykong.ptmshepherd.iterativelocalization.Convolution.convolve;

public class PoissonBinomialLikelihood {
    /**
     * Calculate the probability that site X has the maximum evidence among all sites on the peptide.
     * The steps are as follows:
     * 1) Compute P(X > Y) for every site on the peptide. Ties here are handled by splitting the P(X=Y) between X and Y.
     * The assumption is that even if one site cannot be determined to be the winner, the modification can only be on
     * one site, so there is 50-50 chance that it is on either one. This assumption also simplifies the rest of the
     * computation. Rather than doing n^2 - n convolutions (minus n because there are no X~X comparisons) for every
     * P(X>Y) and P(Y>X), all P(Y>X) are the complement 1-P(X>Y), so can be computed directly.
     * 2) Calculate the probability that each site is the maximum by multiply the probabilities that it beats all the
     * others.
     * 3) Step 2 will not sum to 1 because there are missing probabilities from ties where the most likely site is
     * ambiguous. Make the probability distribution sum to 1 by distributing the remaining weight across the sites.
     *
     * There are two places where we make assumptions about ties. In step 1, ties are distributed equally. This is the
     * most conservative approach, and it cuts the number of convolutions we need to do in half by making the bottom
     * half of the matrix the complement of the top. Two alternative are to model the ties directly in another matrix,
     * but for large peptides with lots of sites the dependencies between them might become cumbersome in the final
     * modeling step. One alternative is to normalize by the sum of P(max). This is much more aggressive, and the best
     * method will need to be tested. In step 3, ties are also distributed equally. NOTE: I think that the only
     * probability that is being missed here (if we distribute ties in step 1) is the probability of MULTIPLE ties.
     * There are two alternatives here:
     * as above, normalize by sum of P(max). Alternatively, compute the PMF of each site being the max, and distribute
     * the probability that it is tied with n other sites by weighing it by p/n. This will sum to more than 1 at the
     * peptide level because there are redundant probabilities calculated among ties, and I'm not sure how to resolve
     * them. We already compute the pairwise ties in the first step, so it is possible to store them in a separate
     * matrix, but it will double the number of convolutions in the first step. Completely ignoring ties at both
     * levels will lead to 0s for sites where there are no matched ions, so even very strong priors will be zeroed out
     * which is undesirable behavior.
     *
     * The current implementation is probably the most conservative and the fastest. If the model can be relaxed, we can
     * either store the ties and distribute their weights among the sites (slower) or normalize by sum of P(max)
     * (faster).
     *
     *
     * Or we could skip this whole thing and just sample it.
     *
     * @param ionProbs
     * @return
     */
    public static double[] calculateProbXMax(double[][] ionProbs) {
        // Initialize array to hold individual X > Y comparisons for all sites
        double[][] probXGreaterThanY = new double[ionProbs.length][ionProbs.length];

        // Fill upper half of array excluding diagonal with prob X > Y
        for (int i = 0; i < ionProbs.length; i++) {
            for (int j = i + 1; j < ionProbs.length; j++) {
                probXGreaterThanY[i][j] = calculateProbXGreaterThanY(ionProbs[i], ionProbs[j]);
            }
        }

        // Fill diagonal
        for (int i = 0; i < probXGreaterThanY.length; i++)
            probXGreaterThanY[i][i] = 1.0;
        // Bottom half is the complement of the top half
        for (int i = 0; i < probXGreaterThanY.length; i++) {
            for(int j = i + 1; j < ionProbs.length; j++) {
                probXGreaterThanY[j][i] = 1.0 - probXGreaterThanY[i][j];
            }
        }

        double[] siteProbs = new double[ionProbs.length];

        // Find probability that each one is the maximum
        double sumProbs = 0.0;
        for (int i = 0; i < probXGreaterThanY.length; i++) {
            double pMax = 1.0;
            for (int j = 0; j < probXGreaterThanY.length; j++)
                pMax *= probXGreaterThanY[i][j];
            siteProbs[i] = pMax;
            sumProbs += pMax;
        }

        // These probabilities will not add to 1 because we do not enumerate all possible ways there can be ties
        // If no winner can be found, balance the remaining probability between sites
        double remainingProbPerSite = (1.0 - sumProbs) / ionProbs.length;
        for (int i = 0; i < siteProbs.length; i++) {
            siteProbs[i] += remainingProbPerSite;
        }

        return siteProbs;

    }

    /**
     * Currently unused. Helper function to strip out self comparison in calculateProbMax matrix and format
     * convolution input.
     * @param probMatrixRow row of the matrix from calculateProbMax
     * @param index site of modification to be analyzed, also row of matrix
     * @return Poisson Binomial convolution input with stripped out self comparison
     */
    private static double[][] transformEPmfParams(double[] probMatrixRow, int index) {
        double[][] ePMFParams = new double[probMatrixRow.length-1][2];
        int cComparison = 0;
        for (int i = 0; i < probMatrixRow.length; i++) {
            if (i == index) // Skip diagonal
                continue;
            ePMFParams[cComparison][0] = 1.0 - probMatrixRow[i];
            ePMFParams[cComparison][1] = probMatrixRow[i];
            cComparison++;
        }

        return ePMFParams;
    }


    /**
     * Calculates the probability that there is more evidence for site X than site Y. X is distributed
     * Poisson Binomial (ionProbsX) and Y is distribtued Poisson Binomial (ionProbsY).
     *
     * Determining these distributions can be intractable for large datasets. Currently, the PMFs are determined by
     * taking the convolutions of the ions probabilities of X and Y. For small repeated calculations, this is likely to
     * be the second fasted algorithm for computation according to Biscarri et al Computational Statistics & Data
     * Analysis (2018). The refined normal approximation will be much faster, but can have large errors with small n's
     * and mean probability close to 0 or 1. It may be possible to implement that for a subset of well-behaved site
     * comparisons if this is too slow. This algorithm is also sped up by ignoring shared ions because their
     * contribution to the PMF of Z is degenerate. The PMF of the distribution Z = X - Y is also distributed
     * Poisson Binomial and can be solved by convolving PMF(X) and -PMF(Y), and where Z is positive is X > Y. Ties
     * (Z = 0) are handled by splitting the mass between X wins and Y wins since only one site can have the PTM.
     *
     * Currently, input is positive is matched and negative if unmatched.
     * Arrays will always be of the same length.
     * @param xIonProbs array of n matched ions' intensities -> probability mappings for X
     * @param yIonProbs array of n matched ions' intensities -> probability mappings for Y
     */
    public static double calculateProbXGreaterThanY(double[] xIonProbs, double[] yIonProbs) {
        // Iterate through arrays once to determine size of output arrays
        int nMatchedIonsX = 0;
        int nMatchedIonsY = 0;
        int nSharedIons = 0;
        for (int i = 0; i < xIonProbs.length; i++) {
            if (xIonProbs[i] > 0.0)
                nMatchedIonsX++;
            if (yIonProbs[i] > 0.0)
                nMatchedIonsY++;
            if (xIonProbs[i] > 0.0 && yIonProbs[i] > 0.0)
                nSharedIons++;
        }

        // Set up arrays to be convolved
        double[][] xBinPoiParams = new double[nMatchedIonsX - nSharedIons][2];
        double[][] yBinPoiParams = new double[nMatchedIonsY - nSharedIons][2];

        // Construct Binomial Poisson parameters for X and Y
        // Dimensions are [nIons][2], [i][0] value holds chance of 0 (fail); [i][1] holds chance of 1 (success)
        // I am very proud of this logic... it does max 2 checks for redundancy and a match to either X or Y
        int xIonPoint = 0;
        int yIonPoint = 0;
        for (int i = 0; i < xIonProbs.length; i++) {
            if (xIonProbs[i] == yIonProbs[i]) // if both matched or both unmatched
                continue;
            if (xIonProbs[i] > 0.0) {
                xBinPoiParams[xIonPoint][0] = 1.0 - xIonProbs[i]; // chance of fail
                xBinPoiParams[xIonPoint][1] = xIonProbs[i]; // chance of success
                xIonPoint++;
            } else {
                yBinPoiParams[yIonPoint][0] = 1.0 - yIonProbs[i]; // chance of fail
                yBinPoiParams[yIonPoint][1] = yIonProbs[i]; // chance of success
                yIonPoint++;
            }
        }

        // Handle cases where there are no ions matched by inserting dummy variable that there is 1 probability of 0
        if (xBinPoiParams.length == 0) {
            xBinPoiParams = new double[][]{
                    {1.0, 0.0}
            };
        }
        if (yBinPoiParams.length == 0) {
            yBinPoiParams = new double[][]{
                    {1.0, 0.0}
            };
        }

        // Compute PMF(X) and PMF(Y)
        double[] xPMF = convolve(xBinPoiParams);
        double[] yPMF = convolve(yBinPoiParams);

        // Compute probability of X > Y, X == Y, and X < Y
        double[] result = calculateXYWinProbability(xPMF, yPMF);

        // Compute probability of X > Y
        double pXGreaterThanY = result[0] + result[1]/2.0; // wins + half of ties.

        return pXGreaterThanY;
    }

    public static void printMatrix(double[][] mat) {
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat.length; j++) {
                System.out.printf("%.2f ", mat[i][j]);
            }
            System.out.println();
        }
    }
}
