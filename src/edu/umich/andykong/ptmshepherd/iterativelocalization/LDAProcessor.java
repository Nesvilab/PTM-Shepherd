package edu.umich.andykong.ptmshepherd.iterativelocalization;

import org.apache.commons.math3.linear.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class LDAProcessor {
    private MatchedIonTable data;
    private Map<Boolean, RealVector> meanVectors;
    public RealMatrix projectedData;
    public List<RealVector> eigenVectors;
    public RealMatrix eigenMatrix;

    public LDAProcessor(MatchedIonTable data) {
        this.data = data;
    }

    public void solveLDA(ExecutorService executorService) throws Exception {
        this.meanVectors = computeMeanVectors();
        RealMatrix sW = computeWithinClassScatterMatrix(executorService);
        RealMatrix sB = computeBetweenClassScatterMatrix(executorService);
        List<RealVector> eigenVectors = solveEigenProblem(sW, sB, 1);
        this.projectedData = projectTrainData(eigenVectors);
    }


    private Map<Boolean, RealVector> computeMeanVectors() {
        Map<Boolean, ArrayList<MatchedIonTable.MatchedIonRow>> classifiedData = classifyData();
        Map<Boolean, RealVector> meanVectors = new HashMap<>();

        for (Map.Entry<Boolean, ArrayList<MatchedIonTable.MatchedIonRow>> entry : classifiedData.entrySet()) {
            Boolean classLabel = entry.getKey();
            ArrayList<MatchedIonTable.MatchedIonRow> rows = entry.getValue();

            // Initialize sum vectors
            RealVector sumVector = new ArrayRealVector(2); // Assuming 2 features for simplicity
            for (MatchedIonTable.MatchedIonRow row : rows) {
                RealVector rowVector = new ArrayRealVector(new double[]{row.massError, row.intensity});
                sumVector = sumVector.add(rowVector);
            }

            // Calculate mean by dividing sumVector by the number of rows for the class
            RealVector meanVector = sumVector.mapDivide(rows.size());
            meanVectors.put(classLabel, meanVector);
        }

        return meanVectors;
    }


    private RealVector computeOverallMean() {
        if (this.data.rows.isEmpty()) {
            return new ArrayRealVector(new double[]{0.0, 0.0});
        }

        double sumMassError = 0.0;
        double sumIntensity = 0.0;

        for (MatchedIonTable.MatchedIonRow row : this.data.rows) {
            sumMassError += row.massError;
            sumIntensity += row.intensity;
        }

        double meanMassError = sumMassError / this.data.rows.size();
        double meanIntensity = sumIntensity / this.data.rows.size();

        return new ArrayRealVector(new double[]{meanMassError, meanIntensity});
    }


    private Map<Boolean, ArrayList<MatchedIonTable.MatchedIonRow>> classifyData() {
        Map<Boolean, ArrayList<MatchedIonTable.MatchedIonRow>> classifiedData = new HashMap<>();
        classifiedData.put(true, new ArrayList<>());
        classifiedData.put(false, new ArrayList<>());

        for (MatchedIonTable.MatchedIonRow row : this.data.rows) {
            classifiedData.get(row.isDecoy).add(row);
        }

        return classifiedData;
    }


    private RealMatrix computeWithinClassScatterMatrix(ExecutorService executorService) throws ExecutionException, InterruptedException {
        // Use Apache Commons Math RealMatrix for matrix representation
        RealMatrix sW = MatrixUtils.createRealMatrix(2, 2); // Assuming 2 features

        List<Future<RealMatrix>> futures = new ArrayList<>();

        for (Boolean classLabel : meanVectors.keySet()) {
            Future<RealMatrix> future = executorService.submit(() -> {
                RealMatrix scatterMatrix = MatrixUtils.createRealMatrix(2, 2);
                RealVector meanVector = meanVectors.get(classLabel);
                for (MatchedIonTable.MatchedIonRow row : this.data.rows) {
                    if (row.isDecoy == classLabel) {
                        RealVector x = new ArrayRealVector(new double[]{row.massError, row.intensity});
                        RealVector xMinusMean = x.subtract(meanVector);
                        RealMatrix outerProduct = xMinusMean.outerProduct(xMinusMean);
                        scatterMatrix = scatterMatrix.add(outerProduct);
                    }
                }
                return scatterMatrix;
            });
            futures.add(future);
        }

        for (Future<RealMatrix> future : futures) {
            sW = sW.add(future.get());
        }

        return sW;
    }


    private RealMatrix computeBetweenClassScatterMatrix(ExecutorService executorService) throws Exception {
        // Step 1: Compute the overall mean of all samples
        RealVector overallMean = computeOverallMean();

        // Step 2: Prepare tasks for parallel execution
        List<Callable<RealMatrix>> tasks = new ArrayList<>();
        for (Map.Entry<Boolean, RealVector> entry : meanVectors.entrySet()) {
            Boolean classLabel = entry.getKey();
            RealVector classMean = entry.getValue();

            tasks.add(() -> {
                // Count the number of samples in this class
                long count = data.stream().filter(row -> row.isDecoy == classLabel).count();

                // Compute the matrix (m_i - m) * (m_i - m)^T
                RealVector meanDiff = classMean.subtract(overallMean);
                RealMatrix sbi = meanDiff.outerProduct(meanDiff).scalarMultiply(count);

                return sbi;
            });
        }

        // Step 3: Execute tasks and aggregate results
        RealMatrix sB = new Array2DRowRealMatrix(meanVectors.get(true).getDimension(), meanVectors.get(true).getDimension());
        for (Future<RealMatrix> future : executorService.invokeAll(tasks)) {
            RealMatrix sbi = future.get(); // This call blocks until the computation is complete
            sB = sB.add(sbi);
        }

        return sB;
    }


    /**
     * Solves the eigenvalue problem for the LDA transformation matrix.
     *
     * @param sW The within-class scatter matrix.
     * @param sB The between-class scatter matrix.
     * @return The matrix of eigenvectors that form the transformation matrix for LDA.
     */
    private List<RealVector> solveEigenProblem(RealMatrix sW, RealMatrix sB, int nEigenVectors) throws IllegalArgumentException {
        // Calculate the inverse (or pseudo-inverse) of sW
        RealMatrix sWInverse = MatrixUtils.inverse(sW);

        // Compute the matrix for which we need to find eigenvalues and eigenvectors: sW^-1 * sB
        RealMatrix targetMatrix = sWInverse.multiply(sB);

        // Perform eigen decomposition on the target matrix
        EigenDecomposition eigenDecomposition = new EigenDecomposition(targetMatrix);

        // Sort eigenvectors by the magnitude of their corresponding eigenvalues in descending order
        List<Integer> sortedIndices = IntStream.range(0, eigenDecomposition.getRealEigenvalues().length)
                .boxed()
                .sorted((i, j) -> Double.compare(eigenDecomposition.getRealEigenvalue(j), eigenDecomposition.getRealEigenvalue(i)))
                .collect(Collectors.toList());

        // Retrieve and store the top nEigenVectors eigenvectors based on sorted eigenvalues
        List<RealVector> eigvenVectors = sortedIndices.stream()
                .limit(nEigenVectors)
                .map(i -> new ArrayRealVector(eigenDecomposition.getEigenvector(i).toArray()))
                .collect(Collectors.toList());


        return eigvenVectors;
    }


    /**
     * Projects the data onto the new axes defined by the discriminants.
     *
     * @param eigenVectors The list of eigenvectors that form the transformation matrix for LDA.
     * @return A RealMatrix containing the projected data.
     */
    private RealMatrix projectTrainData(List<RealVector> eigenVectors) {
        // Convert discriminants to a matrix where each column is a discriminant
        RealMatrix eigenMatrix = new Array2DRowRealMatrix(eigenVectors.get(0).getDimension(), eigenVectors.size());
        for (int i = 0; i < eigenVectors.size(); i++) {
            eigenMatrix.setColumnVector(i, eigenVectors.get(i));
        }

        this.eigenMatrix = eigenMatrix; // TODO break out eigen matrix into new function

        // Convert data to a matrix where each row is a data point
        List<RealVector> dataPoints = data.stream()
                .map(row -> new ArrayRealVector(new double[]{row.massError, row.intensity}))
                .collect(Collectors.toList());

        RealMatrix dataMatrix = new Array2DRowRealMatrix(dataPoints.size(), dataPoints.get(0).getDimension());
        for (int i = 0; i < dataPoints.size(); i++) {
            dataMatrix.setRowVector(i, dataPoints.get(i));
        }

        // Project the data onto the discriminants
        RealMatrix projectedData = dataMatrix.multiply(eigenMatrix);

        return projectedData;
    }

    /**
     * Projects a new data point onto the subspace defined by the LDA eigenvectors.
     *
     * @param massError The mass error feature of the new data point.
     * @param intensity The intensity feature of the new data point.
     * @return The projected data point in the subspace defined by the eigenvectors.
     */
    public RealMatrix projectData(double massError, double intensity) {
        // Create a row matrix from the new data point
        double[][] dataPointArray = {{massError, intensity}};
        RealMatrix dataPointMatrix = new Array2DRowRealMatrix(dataPointArray);

        // Project the data point onto the eigenvector directions
        RealMatrix projectedData = dataPointMatrix.multiply(this.eigenMatrix);

        return projectedData;
    }
}
