package edu.umich.andykong.ptmshepherd.glyco;

import java.util.*;

import edu.umich.andykong.ptmshepherd.PTMShepherd;
import org.apache.commons.math3.linear.*;

/**
 * The Scoring class implements Linear Discriminant Analysis (LDA) to score and
 * filter inferred glycopeptides based on multiple features.
 *
 * This class calculates discriminant functions to separate target PSMs from decoy PSMs,
 * and applies an FDR-based threshold to filter glycopeptides.
 */
public class ScoreLDA {
    /** Collection of feature vectors for decoy PSMs */
    public List<double[]> decoyData;

    /** Collection of feature vectors for target PSMs */
    public List<double[]> targetData;

    /** Score threshold for accepting glycopeptides */
    public double thresholdScore;

    /** Coefficients of the LDA discriminant function */
    private double[] coefficients;

    /** Intercept term of the LDA discriminant function */
    private double intercept;

    /**
     * Initializes the scoring system by creating new data collections.
     * Must be called before any scoring operations are performed.
     */
    public ScoreLDA() {
        decoyData = new ArrayList<>();
        targetData = new ArrayList<>();
    }

    /**
     * Calculates the LDA model and determines the score threshold based on the
     * specified false discovery rate (FDR).
     *
     * @param fdrCutoff The maximum acceptable FDR
     */
    public void runLDA(ArrayList<GlycanAssignmentResult> results, double fdrCutoff) {
        // Calculate mean vectors for target and decoy datasets
        double[] decoyMean = calculateMean(decoyData);
        double[] targetMean = calculateMean(targetData);

        // Calculate covariance matrices
        RealMatrix decoyCov = calculateCovariance(decoyData, decoyMean);
        RealMatrix targetCov = calculateCovariance(targetData, targetMean);

        // Calculate pooled covariance matrix
        RealMatrix pooledCov = calculatePooledCovariance(targetCov, decoyCov,
                targetData.size(), decoyData.size());

        // Check for singularity and regularize if necessary
        double det = new LUDecomposition(pooledCov).getDeterminant();
        if (Double.isNaN(det) || det == 0) {
            PTMShepherd.print("Warning: Pooled covariance matrix is singular or invalid. Regularizing...");
            pooledCov = regularizeCovarianceMatrix(pooledCov, 0.001);
        }

        // Invert the covariance matrix
        LUDecomposition lu = new LUDecomposition(pooledCov);
        RealMatrix invPooledCov = lu.getSolver().getInverse();

        // Calculate LDA model parameters
        coefficients = calculateCoefficients(targetMean, decoyMean, invPooledCov);
        intercept = calculateIntercept(targetMean, decoyMean, coefficients);

        // Output model coefficients
        PTMShepherd.print("Coefficients: " + java.util.Arrays.toString(coefficients));

        // Calculate scores and determine threshold
        for (GlycanAssignmentResult result: results) {
            result.glycanScore = calculateScore(result.featureVec);
        }
        thresholdScore = getCutOff(results, fdrCutoff);
    }


    //--------------------------------
    // SCORING CALCULATION METHODS
    //--------------------------------

    /**
     * Calculates the LDA score for a data point.
     *
     * @param dataPoint Feature vector to score
     * @return The calculated score
     */
    public double calculateScore(double[] dataPoint) {
        return dotProduct(dataPoint, coefficients) + intercept;
    }

    /**
     * Determines the score threshold based on the specified FDR.
     *
     * @param fdrCutOff Maximum acceptable FDR
     * @param results List of GlycanAssignmentResults containing target and decoy scores
     * @return The score threshold
     */
    public static double getCutOff(List<GlycanAssignmentResult> results, double fdrCutOff) {
        // Sort target scores in descending order
        results.sort(Comparator.comparingDouble((GlycanAssignmentResult result) -> result.glycanScore).reversed());

        // Get decoy indices in the combined sorted list
        List<Integer> decoyIndexes = getDecoyIndexes(results);

        // Calculate FDR at each point
        int decoyCount = decoyIndexes.size();
        int targetCount = results.size() - decoyCount;
        double currentMinQ = Double.MAX_VALUE;
        double scoreThreshold = Double.NaN;
        boolean foundThreshold = false;
        for (int i = results.size() - 1; i >= 0; i--) {
            if (results.get(i).isDecoyGlycan) {
                decoyCount--;
            } else {
                targetCount--;
            }
            double fdr = Math.min((double) (decoyCount + 1) / targetCount, currentMinQ);        // q = (d+1)/t recommended per 10.1021/acs.jproteome.6b00144
            if (fdr < currentMinQ) {
                currentMinQ = fdr;
            }
            if (!foundThreshold) {
                if (fdr < fdrCutOff) {
                    scoreThreshold = results.get(i).glycanScore;
                    PTMShepherd.print(String.format("Found LDA score threshold: %.2f with %d decoys, %d targets for %.4f estimated FDR (%d total inputs)",
                            results.get(i).glycanScore, decoyCount, targetCount, fdr, results.size()));
                    foundThreshold = true;
                }
            }
            results.get(i).glycanQval = fdr;
        }
        return scoreThreshold;
    }

    //--------------------------------
    // HELPER METHODS
    //--------------------------------

    /**
     * Gets the indices of decoy scores in the combined sorted score list.
     */
    private static List<Integer> getDecoyIndexes(List<GlycanAssignmentResult> results) {
        List<Integer> decoyIndexes = new ArrayList<>();
        for (int i = 0; i < results.size(); i++) {
            if (results.get(i).isDecoyGlycan) {
                decoyIndexes.add(i);
            }
        }
        return decoyIndexes;
    }

    /**
     * Finds the index of the last value in a list that is below a threshold.
     */
    public static Integer getLastValueBelowThreshold(List<Double> list, double threshold) {
        for (int i = list.size() - 1; i >= 0; i--) {
            if (list.get(i) <= threshold) {
                return i;
            }
        }
        return 0; // Default if no value found
    }

    /**
     * Calculates the mean vector for a dataset.
     */
    private static double[] calculateMean(List<double[]> data) {
        int featureCount = data.get(0).length;
        double[] mean = new double[featureCount];

        // Sum values for each feature
        for (double[] dataPoint : data) {
            for (int i = 0; i < featureCount; i++) {
                mean[i] += dataPoint[i];
            }
        }

        // Divide by count to get mean
        for (int i = 0; i < featureCount; i++) {
            mean[i] /= data.size();
        }

        return mean;
    }

    /**
     * Calculates the covariance matrix for a dataset.
     */
    private static RealMatrix calculateCovariance(List<double[]> data, double[] mean) {
        int featureCount = data.get(0).length;
        RealMatrix covariance = new Array2DRowRealMatrix(featureCount, featureCount);

        // Calculate covariance sums
        for (double[] dataPoint : data) {
            for (int i = 0; i < featureCount; i++) {
                for (int j = 0; j < featureCount; j++) {
                    double currentValue = covariance.getEntry(i, j);
                    double update = (dataPoint[i] - mean[i]) * (dataPoint[j] - mean[j]);
                    covariance.setEntry(i, j, currentValue + update);
                }
            }
        }

        // Normalize by degrees of freedom
        for (int i = 0; i < featureCount; i++) {
            for (int j = 0; j < featureCount; j++) {
                double normalizedValue = covariance.getEntry(i, j) / (data.size() - 1);
                covariance.setEntry(i, j, normalizedValue);
            }
        }

        return covariance;
    }

    /**
     * Calculates the pooled covariance matrix from two class-specific matrices.
     */
    private static RealMatrix calculatePooledCovariance(RealMatrix cov1, RealMatrix cov2, int n1, int n2) {
        double weight1 = (double) (n1 - 1) / (n1 + n2 - 2);
        double weight2 = (double) (n2 - 1) / (n1 + n2 - 2);
        return cov1.scalarMultiply(weight1).add(cov2.scalarMultiply(weight2));
    }

    /**
     * Regularizes a covariance matrix to ensure it's invertible.
     */
    private static RealMatrix regularizeCovarianceMatrix(RealMatrix covarianceMatrix, double alpha) {
        int size = covarianceMatrix.getRowDimension();
        RealMatrix identity = MatrixUtils.createRealIdentityMatrix(size);
        return covarianceMatrix.scalarMultiply(1 - alpha).add(identity.scalarMultiply(alpha));
    }

    /**
     * Calculates the coefficients for the LDA discriminant function.
     */
    private static double[] calculateCoefficients(double[] mean1, double[] mean2, RealMatrix invPooledCov) {
        double[] coefficients = new double[mean1.length];
        double[] difference = new double[mean1.length];

        // Calculate mean difference vector
        for (int i = 0; i < mean1.length; i++) {
            difference[i] = mean1[i] - mean2[i];
        }

        // Calculate coefficients using matrix operations
        ArrayRealVector diffRealVector = new ArrayRealVector(difference);
        for (int i = 0; i < mean1.length; i++) {
            coefficients[i] = invPooledCov.getRowVector(i).dotProduct(diffRealVector);
        }

        return coefficients;
    }

    /**
     * Calculates the intercept term for the LDA discriminant function.
     */
    private static double calculateIntercept(double[] mean1, double[] mean2, double[] coefficients) {
        return -0.5 * (dotProduct(mean1, coefficients) + dotProduct(mean2, coefficients));
    }

    /**
     * Calculates the dot product of two vectors.
     */
    private static double dotProduct(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("Arrays must have the same length");
        }

        double result = 0;
        for (int i = 0; i < a.length; i++) {
            result += a[i] * b[i];
        }

        return result;
    }
}



