package edu.umich.andykong.ptmshepherd.glyco;

import java.util.*;

import edu.umich.andykong.ptmshepherd.PTMShepherd;

/**
 * Neural network scorer for glycan composition assignment. Trains a feedforward MLP
 * to separate target (correct) from decoy (incorrect) glycan assignments using the
 * same feature vectors as the LDA scorer, but captures nonlinear patterns in the data.
 *
 * Architecture: Input -> Hidden1 (ReLU, dropout) -> Hidden2 (ReLU, dropout) -> Sigmoid output
 * Training: Mini-batch Adam optimizer with L2 weight decay, label smoothing, and early stopping.
 * Features are z-score standardized before training and inference.
 */
public class ScoreNN {

    /** Feature vectors from target glycan assignments (training positive class) */
    public List<double[]> targetData;

    /** Feature vectors from decoy glycan assignments (training negative class) */
    public List<double[]> decoyData;

    // Network dimensions
    private int inputSize;
    private int hidden1Size;
    private int hidden2Size;

    // Weights and biases
    private double[][] w1;   // [inputSize][hidden1Size]
    private double[] b1;     // [hidden1Size]
    private double[][] w2;   // [hidden1Size][hidden2Size]
    private double[] b2;     // [hidden2Size]
    private double[] w3;     // [hidden2Size] -> single output
    private double b3;

    // Feature standardization parameters (computed from training data)
    private double[] featureMean;
    private double[] featureStd;

    // Hyperparameters
    private static final int DEFAULT_HIDDEN1 = 32;
    private static final int DEFAULT_HIDDEN2 = 16;
    private static final int SMALL_HIDDEN1 = 16;
    private static final int SMALL_HIDDEN2 = 8;
    private static final int SMALL_INPUT_THRESHOLD = 5;
    private static final double LEARNING_RATE = 0.001;
    private static final double WEIGHT_DECAY = 1e-4;
    private static final double LABEL_SMOOTHING = 0.1;
    private static final double DROPOUT_RATE = 0.3;
    private static final int BATCH_SIZE = 64;
    private static final int MAX_EPOCHS = 200;
    private static final int PATIENCE = 20;
    private static final double GRAD_CLIP_NORM = 5.0;
    private static final int MIN_TRAINING_SAMPLES = 10;

    // Adam optimizer constants
    private static final double BETA1 = 0.9;
    private static final double BETA2 = 0.999;
    private static final double ADAM_EPS = 1e-8;

    // Adam first/second moment estimates
    private double[][] mW1, vW1, mW2, vW2;
    private double[] mB1, vB1, mB2, vB2, mW3, vW3;
    private double mB3, vB3;
    private int adamT;

    private final Random randomGenerator;

    public ScoreNN(Random randomGenerator) {
        targetData = new ArrayList<>();
        decoyData = new ArrayList<>();
        this.randomGenerator = randomGenerator;
    }

    /**
     * Train the neural network on collected target/decoy data and score all results.
     * Same interface as ScoreLDA.runLDA().
     *
     * @param results   all glycan assignment results to score
     * @param ldaHeader feature name header string for logging
     * @param targetProp proportion of top-scoring targets to use for training
     */
    public void runNN(ArrayList<GlycanAssignmentResult> results, String ldaHeader, double targetProp) {
        filterTargetData(targetProp);

        if (targetData.size() < MIN_TRAINING_SAMPLES || decoyData.size() < MIN_TRAINING_SAMPLES) {
            PTMShepherd.print("\tWarning: insufficient data for NN scoring (targets=" + targetData.size() +
                    ", decoys=" + decoyData.size() + "). Using summed scores instead.");
            return;
        }

        inputSize = targetData.get(0).length;
        hidden1Size = inputSize < SMALL_INPUT_THRESHOLD ? SMALL_HIDDEN1 : DEFAULT_HIDDEN1;
        hidden2Size = inputSize < SMALL_INPUT_THRESHOLD ? SMALL_HIDDEN2 : DEFAULT_HIDDEN2;

        computeStandardization();

        double[][] X = buildFeatureMatrix();
        double[] y = buildLabels();

        initWeights();
        initAdam();
        int epochsTrained = train(X, y);

        PTMShepherd.print(String.format("\tNN trained: %d features, [%d, %d] hidden, %d epochs, %d targets, %d decoys",
                inputSize, hidden1Size, hidden2Size, epochsTrained, targetData.size(), decoyData.size()));
        PTMShepherd.print("\tNN feature names: " + ldaHeader);

        // Apply scores to all candidates (same pattern as ScoreLDA.runLDA)
        for (GlycanAssignmentResult result : results) {
            if (result.foundGlycan) {
                for (GlycanCandidateResult candidate : result.allCandidates) {
                    candidate.nnScore = score(candidate.featureVec);
                    candidate.glycanScore = candidate.nnScore;
                }
                result.glycanScore = result.bestCandidate.nnScore;
            }
        }
    }

    /**
     * Score a single feature vector using the trained network (no dropout at inference).
     *
     * @param features raw (unstandardized) feature vector
     * @return sigmoid output in [0,1]; higher = more target-like
     */
    public double score(double[] features) {
        double[] x = standardize(features);
        double[] h1 = forwardDense(x, w1, b1, true);
        double[] h2 = forwardDense(h1, w2, b2, true);
        return sigmoid(dot(h2, w3) + b3);
    }

    // ================================================================
    // Training loop
    // ================================================================

    /**
     * Train the network using mini-batch Adam with early stopping.
     *
     * @return number of epochs trained
     */
    private int train(double[][] X, double[] y) {
        int n = X.length;
        int[] indices = new int[n];
        for (int i = 0; i < n; i++) indices[i] = i;

        double bestLoss = Double.MAX_VALUE;
        int patienceCounter = 0;
        int epoch;

        // Best weights (for early stopping restore)
        double[][] bestW1 = null, bestW2 = null;
        double[] bestB1 = null, bestB2 = null, bestW3out = null;
        double bestB3 = 0;

        for (epoch = 0; epoch < MAX_EPOCHS; epoch++) {
            shuffle(indices);
            double epochLoss = 0;

            for (int batchStart = 0; batchStart < n; batchStart += BATCH_SIZE) {
                int batchEnd = Math.min(batchStart + BATCH_SIZE, n);
                int bs = batchEnd - batchStart;

                // Gradient accumulators
                double[][] gW1 = new double[inputSize][hidden1Size];
                double[] gB1 = new double[hidden1Size];
                double[][] gW2 = new double[hidden1Size][hidden2Size];
                double[] gB2 = new double[hidden2Size];
                double[] gW3 = new double[hidden2Size];
                double gB3local = 0;
                double batchLoss = 0;

                for (int bi = batchStart; bi < batchEnd; bi++) {
                    int idx = indices[bi];
                    double[] x = X[idx];
                    double label = y[idx];

                    // ---- Forward pass with dropout ----

                    // Hidden layer 1
                    double[] z1 = new double[hidden1Size];
                    double[] a1 = new double[hidden1Size];
                    boolean[] drop1 = new boolean[hidden1Size];
                    for (int j = 0; j < hidden1Size; j++) {
                        double sum = b1[j];
                        for (int k = 0; k < inputSize; k++) sum += x[k] * w1[k][j];
                        z1[j] = sum;
                        drop1[j] = randomGenerator.nextDouble() < DROPOUT_RATE;
                        a1[j] = drop1[j] ? 0.0 : Math.max(0, z1[j]) / (1.0 - DROPOUT_RATE);
                    }

                    // Hidden layer 2
                    double[] z2 = new double[hidden2Size];
                    double[] a2 = new double[hidden2Size];
                    boolean[] drop2 = new boolean[hidden2Size];
                    for (int j = 0; j < hidden2Size; j++) {
                        double sum = b2[j];
                        for (int k = 0; k < hidden1Size; k++) sum += a1[k] * w2[k][j];
                        z2[j] = sum;
                        drop2[j] = randomGenerator.nextDouble() < DROPOUT_RATE;
                        a2[j] = drop2[j] ? 0.0 : Math.max(0, z2[j]) / (1.0 - DROPOUT_RATE);
                    }

                    // Output
                    double z3 = b3;
                    for (int k = 0; k < hidden2Size; k++) z3 += a2[k] * w3[k];
                    double pred = sigmoid(z3);

                    // Binary cross-entropy loss
                    double clamped = Math.max(1e-7, Math.min(1.0 - 1e-7, pred));
                    batchLoss += -(label * Math.log(clamped) + (1.0 - label) * Math.log(1.0 - clamped));

                    // ---- Backward pass ----

                    // Output gradient (BCE + sigmoid combined)
                    double dz3 = pred - label;

                    // w3, b3 gradients
                    for (int k = 0; k < hidden2Size; k++) gW3[k] += dz3 * a2[k];
                    gB3local += dz3;

                    // Backprop through hidden layer 2
                    double[] dz2 = new double[hidden2Size];
                    for (int j = 0; j < hidden2Size; j++) {
                        if (drop2[j]) {
                            dz2[j] = 0;
                        } else {
                            double da2j = dz3 * w3[j];
                            dz2[j] = z2[j] > 0 ? da2j / (1.0 - DROPOUT_RATE) : 0;
                        }
                    }
                    for (int j = 0; j < hidden2Size; j++) {
                        for (int k = 0; k < hidden1Size; k++) gW2[k][j] += dz2[j] * a1[k];
                        gB2[j] += dz2[j];
                    }

                    // Backprop through hidden layer 1
                    double[] da1 = new double[hidden1Size];
                    for (int j = 0; j < hidden1Size; j++) {
                        for (int k = 0; k < hidden2Size; k++) da1[j] += dz2[k] * w2[j][k];
                    }
                    double[] dz1 = new double[hidden1Size];
                    for (int j = 0; j < hidden1Size; j++) {
                        if (drop1[j]) {
                            dz1[j] = 0;
                        } else {
                            dz1[j] = z1[j] > 0 ? da1[j] / (1.0 - DROPOUT_RATE) : 0;
                        }
                    }
                    for (int j = 0; j < hidden1Size; j++) {
                        for (int k = 0; k < inputSize; k++) gW1[k][j] += dz1[j] * x[k];
                        gB1[j] += dz1[j];
                    }
                }

                // Average gradients over batch and add L2 weight decay
                double invBs = 1.0 / bs;
                for (int i = 0; i < inputSize; i++)
                    for (int j = 0; j < hidden1Size; j++)
                        gW1[i][j] = gW1[i][j] * invBs + WEIGHT_DECAY * w1[i][j];
                for (int j = 0; j < hidden1Size; j++) gB1[j] *= invBs;

                for (int i = 0; i < hidden1Size; i++)
                    for (int j = 0; j < hidden2Size; j++)
                        gW2[i][j] = gW2[i][j] * invBs + WEIGHT_DECAY * w2[i][j];
                for (int j = 0; j < hidden2Size; j++) gB2[j] *= invBs;

                for (int j = 0; j < hidden2Size; j++)
                    gW3[j] = gW3[j] * invBs + WEIGHT_DECAY * w3[j];
                gB3local *= invBs;

                // Gradient clipping (global norm)
                double clipScale = computeClipScale(gW1, gB1, gW2, gB2, gW3, gB3local);
                if (clipScale < 1.0) {
                    scale2D(gW1, clipScale);
                    scale1D(gB1, clipScale);
                    scale2D(gW2, clipScale);
                    scale1D(gB2, clipScale);
                    scale1D(gW3, clipScale);
                    gB3local *= clipScale;
                }

                // Adam parameter update
                adamUpdate(gW1, gB1, gW2, gB2, gW3, gB3local);

                epochLoss += batchLoss;
            }

            epochLoss /= n;

            // Early stopping with best-weight checkpointing
            if (epochLoss < bestLoss - 1e-6) {
                bestLoss = epochLoss;
                patienceCounter = 0;
                bestW1 = deepCopy2D(w1);
                bestB1 = b1.clone();
                bestW2 = deepCopy2D(w2);
                bestB2 = b2.clone();
                bestW3out = w3.clone();
                bestB3 = b3;
            } else {
                patienceCounter++;
                if (patienceCounter >= PATIENCE) {
                    epoch++;
                    break;
                }
            }
        }

        // Restore best weights
        if (bestW1 != null) {
            w1 = bestW1;
            b1 = bestB1;
            w2 = bestW2;
            b2 = bestB2;
            w3 = bestW3out;
            b3 = bestB3;
        }

        return epoch;
    }

    // ================================================================
    // Network layer operations
    // ================================================================

    /**
     * Forward pass through a dense layer: output = relu(input * weights + bias) or linear if !relu.
     */
    private static double[] forwardDense(double[] input, double[][] weights, double[] bias, boolean relu) {
        int outSize = bias.length;
        double[] output = new double[outSize];
        for (int j = 0; j < outSize; j++) {
            double sum = bias[j];
            for (int k = 0; k < input.length; k++) sum += input[k] * weights[k][j];
            output[j] = relu ? Math.max(0, sum) : sum;
        }
        return output;
    }

    private static double sigmoid(double x) {
        if (x >= 0) {
            return 1.0 / (1.0 + Math.exp(-x));
        } else {
            double ex = Math.exp(x);
            return ex / (1.0 + ex);
        }
    }

    private static double dot(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) sum += a[i] * b[i];
        return sum;
    }

    // ================================================================
    // Weight initialization (Xavier/Glorot)
    // ================================================================

    private void initWeights() {
        w1 = xavierInit(inputSize, hidden1Size);
        b1 = new double[hidden1Size];
        w2 = xavierInit(hidden1Size, hidden2Size);
        b2 = new double[hidden2Size];
        w3 = xavierInit1D(hidden2Size);
        b3 = 0;
    }

    private double[][] xavierInit(int fanIn, int fanOut) {
        double scale = Math.sqrt(2.0 / (fanIn + fanOut));
        double[][] w = new double[fanIn][fanOut];
        for (int i = 0; i < fanIn; i++)
            for (int j = 0; j < fanOut; j++)
                w[i][j] = randomGenerator.nextGaussian() * scale;
        return w;
    }

    private double[] xavierInit1D(int fanIn) {
        double scale = Math.sqrt(2.0 / (fanIn + 1));
        double[] w = new double[fanIn];
        for (int i = 0; i < fanIn; i++) w[i] = randomGenerator.nextGaussian() * scale;
        return w;
    }

    // ================================================================
    // Adam optimizer
    // ================================================================

    private void initAdam() {
        mW1 = new double[inputSize][hidden1Size];
        vW1 = new double[inputSize][hidden1Size];
        mB1 = new double[hidden1Size];
        vB1 = new double[hidden1Size];
        mW2 = new double[hidden1Size][hidden2Size];
        vW2 = new double[hidden1Size][hidden2Size];
        mB2 = new double[hidden2Size];
        vB2 = new double[hidden2Size];
        mW3 = new double[hidden2Size];
        vW3 = new double[hidden2Size];
        mB3 = 0;
        vB3 = 0;
        adamT = 0;
    }

    private void adamUpdate(double[][] gW1, double[] gB1, double[][] gW2, double[] gB2, double[] gW3, double gB3val) {
        adamT++;
        double bc1 = 1.0 - Math.pow(BETA1, adamT);
        double bc2 = 1.0 - Math.pow(BETA2, adamT);

        // w1
        for (int i = 0; i < inputSize; i++) {
            for (int j = 0; j < hidden1Size; j++) {
                mW1[i][j] = BETA1 * mW1[i][j] + (1 - BETA1) * gW1[i][j];
                vW1[i][j] = BETA2 * vW1[i][j] + (1 - BETA2) * gW1[i][j] * gW1[i][j];
                w1[i][j] -= LEARNING_RATE * (mW1[i][j] / bc1) / (Math.sqrt(vW1[i][j] / bc2) + ADAM_EPS);
            }
        }
        // b1
        for (int j = 0; j < hidden1Size; j++) {
            mB1[j] = BETA1 * mB1[j] + (1 - BETA1) * gB1[j];
            vB1[j] = BETA2 * vB1[j] + (1 - BETA2) * gB1[j] * gB1[j];
            b1[j] -= LEARNING_RATE * (mB1[j] / bc1) / (Math.sqrt(vB1[j] / bc2) + ADAM_EPS);
        }
        // w2
        for (int i = 0; i < hidden1Size; i++) {
            for (int j = 0; j < hidden2Size; j++) {
                mW2[i][j] = BETA1 * mW2[i][j] + (1 - BETA1) * gW2[i][j];
                vW2[i][j] = BETA2 * vW2[i][j] + (1 - BETA2) * gW2[i][j] * gW2[i][j];
                w2[i][j] -= LEARNING_RATE * (mW2[i][j] / bc1) / (Math.sqrt(vW2[i][j] / bc2) + ADAM_EPS);
            }
        }
        // b2
        for (int j = 0; j < hidden2Size; j++) {
            mB2[j] = BETA1 * mB2[j] + (1 - BETA1) * gB2[j];
            vB2[j] = BETA2 * vB2[j] + (1 - BETA2) * gB2[j] * gB2[j];
            b2[j] -= LEARNING_RATE * (mB2[j] / bc1) / (Math.sqrt(vB2[j] / bc2) + ADAM_EPS);
        }
        // w3
        for (int j = 0; j < hidden2Size; j++) {
            this.mW3[j] = BETA1 * this.mW3[j] + (1 - BETA1) * gW3[j];
            this.vW3[j] = BETA2 * this.vW3[j] + (1 - BETA2) * gW3[j] * gW3[j];
            w3[j] -= LEARNING_RATE * (this.mW3[j] / bc1) / (Math.sqrt(this.vW3[j] / bc2) + ADAM_EPS);
        }
        // b3
        this.mB3 = BETA1 * this.mB3 + (1 - BETA1) * gB3val;
        this.vB3 = BETA2 * this.vB3 + (1 - BETA2) * gB3val * gB3val;
        b3 -= LEARNING_RATE * (this.mB3 / bc1) / (Math.sqrt(this.vB3 / bc2) + ADAM_EPS);
    }

    // ================================================================
    // Data preparation
    // ================================================================

    /**
     * Filter target data to keep only top-scoring fraction (same approach as ScoreLDA).
     * Uses sum of feature values as a proxy quality metric since no model score exists yet.
     */
    private void filterTargetData(double targetProp) {
        Map<Integer, Double> scored = new HashMap<>();
        for (int i = 0; i < targetData.size(); i++) {
            double sum = 0;
            for (double v : targetData.get(i)) sum += v;
            scored.put(i, sum);
        }
        List<Map.Entry<Integer, Double>> sorted = new ArrayList<>(scored.entrySet());
        sorted.sort(Map.Entry.comparingByValue(Comparator.reverseOrder()));

        int cutoff = (int) Math.ceil(targetData.size() * targetProp);
        List<double[]> filtered = new ArrayList<>();
        for (int i = 0; i < cutoff && i < sorted.size(); i++) {
            filtered.add(targetData.get(sorted.get(i).getKey()));
        }
        targetData = filtered;
    }

    /**
     * Compute z-score standardization parameters from all training data.
     */
    private void computeStandardization() {
        int nFeatures = targetData.get(0).length;
        featureMean = new double[nFeatures];
        featureStd = new double[nFeatures];

        List<double[]> allData = new ArrayList<>(targetData);
        allData.addAll(decoyData);
        int n = allData.size();

        for (double[] x : allData) {
            for (int i = 0; i < nFeatures; i++) featureMean[i] += x[i];
        }
        for (int i = 0; i < nFeatures; i++) featureMean[i] /= n;

        for (double[] x : allData) {
            for (int i = 0; i < nFeatures; i++) {
                double d = x[i] - featureMean[i];
                featureStd[i] += d * d;
            }
        }
        for (int i = 0; i < nFeatures; i++) {
            featureStd[i] = Math.sqrt(featureStd[i] / n);
            if (featureStd[i] < 1e-10) featureStd[i] = 1.0; // constant feature: no scaling
        }
    }

    private double[] standardize(double[] x) {
        double[] result = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            result[i] = (x[i] - featureMean[i]) / featureStd[i];
        }
        return result;
    }

    /**
     * Build standardized feature matrix from target + decoy data (targets first).
     */
    private double[][] buildFeatureMatrix() {
        int n = targetData.size() + decoyData.size();
        double[][] X = new double[n][];
        int idx = 0;
        for (double[] x : targetData) X[idx++] = standardize(x);
        for (double[] x : decoyData) X[idx++] = standardize(x);
        return X;
    }

    /**
     * Build label array with label smoothing: targets get (1 - smoothing), decoys get smoothing.
     */
    private double[] buildLabels() {
        int n = targetData.size() + decoyData.size();
        double[] y = new double[n];
        double posLabel = 1.0 - LABEL_SMOOTHING;
        double negLabel = LABEL_SMOOTHING;
        for (int i = 0; i < targetData.size(); i++) y[i] = posLabel;
        for (int i = targetData.size(); i < n; i++) y[i] = negLabel;
        return y;
    }

    // ================================================================
    // Utilities
    // ================================================================

    private void shuffle(int[] arr) {
        for (int i = arr.length - 1; i > 0; i--) {
            int j = randomGenerator.nextInt(i + 1);
            int tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }

    /**
     * Compute global gradient norm and return scale factor for clipping.
     * Returns 1.0 if no clipping needed, < 1.0 if gradients should be scaled down.
     */
    private static double computeClipScale(double[][] gW1, double[] gB1, double[][] gW2, double[] gB2,
                                            double[] gW3, double gB3) {
        double norm = 0;
        for (double[] row : gW1) for (double v : row) norm += v * v;
        for (double v : gB1) norm += v * v;
        for (double[] row : gW2) for (double v : row) norm += v * v;
        for (double v : gB2) norm += v * v;
        for (double v : gW3) norm += v * v;
        norm += gB3 * gB3;
        norm = Math.sqrt(norm);
        return norm > GRAD_CLIP_NORM ? GRAD_CLIP_NORM / norm : 1.0;
    }

    private static void scale2D(double[][] arr, double s) {
        for (double[] row : arr)
            for (int j = 0; j < row.length; j++)
                row[j] *= s;
    }

    private static void scale1D(double[] arr, double s) {
        for (int i = 0; i < arr.length; i++) arr[i] *= s;
    }

    private static double[][] deepCopy2D(double[][] src) {
        double[][] copy = new double[src.length][];
        for (int i = 0; i < src.length; i++) copy[i] = src[i].clone();
        return copy;
    }
}
