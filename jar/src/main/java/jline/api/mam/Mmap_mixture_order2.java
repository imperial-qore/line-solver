package jline.api.mam;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.Pair;

public final class Mmap_mixture_order2 {
    private Mmap_mixture_order2() {}

    /**
     * Creates a second-order MMAP mixture from a collection of MMAPs.
     * This function reduces higher-order MMAPs to order 2 and creates a mixture.
     *
     * @param mmaps   List of MMAPs to mix
     * @param weights Mixing weights (must sum to 1)
     * @return Second-order MMAP mixture
     */
    public static MatrixCell mmap_mixture_order2(List<MatrixCell> mmaps, double[] weights) {
        if (mmaps.isEmpty()) {
            throw new IllegalArgumentException("Must provide at least one MMAP");
        }
        if (weights.length != mmaps.size()) {
            throw new IllegalArgumentException("Number of weights must match number of MMAPs");
        }
        double sum = 0.0;
        for (double w : weights) sum += w;
        if (Math.abs(sum - 1.0) >= 1e-10) {
            throw new IllegalArgumentException("Weights must sum to 1");
        }

        // First, reduce each MMAP to order 2 if necessary
        List<MatrixCell> order2Mmaps = new ArrayList<MatrixCell>(mmaps.size());
        for (MatrixCell mmap : mmaps) {
            if (mmap.get(0).getNumRows() <= 2) {
                order2Mmaps.add(mmap);
            } else {
                order2Mmaps.add(reduceToOrder2(mmap));
            }
        }

        // Create mixture using the regular mixture function
        Map<Integer, MatrixCell> mapsMap = new HashMap<Integer, MatrixCell>();
        for (int index = 0; index < order2Mmaps.size(); index++) {
            mapsMap.put(index, order2Mmaps.get(index));
        }
        Matrix alpha = new Matrix(1, weights.length);
        for (int index = 0; index < weights.length; index++) {
            alpha.set(0, index, weights[index]);
        }
        return Mmap_mixture.mmap_mixture(alpha, mapsMap);
    }

    /**
     * Reduces an MMAP to order 2 using moment matching.
     */
    private static MatrixCell reduceToOrder2(MatrixCell mmap) {
        int originalOrder = mmap.get(0).getNumRows();
        int numClasses = mmap.size() - 2;

        if (originalOrder <= 2) {
            return mmap;
        }

        // Extract key characteristics from original MMAP
        double M1 = Map_moment.map_moment(mmap, 1);
        double M2 = Map_moment.map_moment(mmap, 2);
        double M3 = Map_moment.map_moment(mmap, 3);
        Matrix classProbs = Mmap_pc.mmap_pc(mmap);

        // Create order-2 approximation using AMAP2 fitting
        Pair<MatrixCell, List<MatrixCell>> baseAmap2 = Amap2_fit_gamma.amap2_fit_gamma(M1, M2, M3, 0.0);
        MatrixCell baseOrder2 = baseAmap2.getFirst();
        if (baseOrder2 == null) {
            throw new RuntimeException("Failed to create order-2 approximation");
        }

        // Convert to multi-class MMAP by distributing classes
        MatrixCell order2Mmap = new MatrixCell(numClasses + 2);
        order2Mmap.set(0, baseOrder2.get(0)); // D0
        order2Mmap.set(1, new Matrix(2, 2)); // Total arrivals matrix

        // Distribute class arrivals
        for (int c = 0; c < numClasses; c++) {
            Matrix classMatrix = new Matrix(2, 2);
            double classWeight = (c < classProbs.length()) ? classProbs.get(0, c) : 1.0 / numClasses;

            // Scale base arrival matrix by class probability
            Matrix D1 = baseOrder2.get(1);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    double value = D1.get(i, j) * classWeight;
                    classMatrix.set(i, j, value);
                    order2Mmap.get(1).set(i, j, order2Mmap.get(1).get(i, j) + value); // Add to total arrivals
                }
            }

            order2Mmap.set(c + 2, classMatrix);
        }

        return order2Mmap;
    }

    /**
     * Functional callback used internally for weight-combination enumeration.
     */
    public interface WeightsCallback {
        void apply(double[] weights);
    }

    /**
     * Creates a second-order MMAP mixture with automatic weight selection.
     * Weights are chosen to minimize the approximation error.
     *
     * @param mmaps                 List of MMAPs to mix
     * @param targetCharacteristics Target moments to match (M1, M2, M3)
     * @return Pair of (optimized mixture, optimal weights)
     */
    public static Pair<MatrixCell, double[]> mmap_mixture_order2_optimal(
            List<MatrixCell> mmaps,
            double[] targetCharacteristics) {

        if (mmaps.isEmpty()) {
            throw new IllegalArgumentException("Must provide at least one MMAP");
        }
        if (targetCharacteristics.length < 3) {
            throw new IllegalArgumentException("Must provide at least 3 target characteristics (M1, M2, M3)");
        }

        final int n = mmaps.size();

        // Compute characteristics of each component MMAP
        final double[][] componentCharacteristics = new double[n][3];
        for (int idx = 0; idx < n; idx++) {
            MatrixCell mmap = mmaps.get(idx);
            componentCharacteristics[idx][0] = Map_moment.map_moment(mmap, 1);
            componentCharacteristics[idx][1] = Map_moment.map_moment(mmap, 2);
            componentCharacteristics[idx][2] = Map_moment.map_moment(mmap, 3);
        }

        // Find optimal weights using simple optimization
        final double[][] bestWeightsHolder = new double[1][];
        bestWeightsHolder[0] = new double[n];
        for (int i = 0; i < n; i++) bestWeightsHolder[0][i] = 1.0 / n; // uniform start
        final double[] bestErrorHolder = new double[]{Double.MAX_VALUE};
        final double[] target = targetCharacteristics;

        // Grid search for optimal weights (simplified approach)
        int gridSize = 20;
        generateWeightCombinations(n, gridSize, new WeightsCallback() {
            @Override
            public void apply(double[] weights) {
                // Compute mixture characteristics
                double[] mixtureChar = new double[3];
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < 3; j++) {
                        mixtureChar[j] += weights[i] * componentCharacteristics[i][j];
                    }
                }

                // Compute error
                double error = 0.0;
                for (int j = 0; j < 3; j++) {
                    double relativeError = Math.abs(mixtureChar[j] - target[j]) /
                            Math.max(target[j], 1e-6);
                    error += relativeError * relativeError;
                }

                if (error < bestErrorHolder[0]) {
                    bestErrorHolder[0] = error;
                    bestWeightsHolder[0] = weights.clone();
                }
            }
        });

        double[] bestWeights = bestWeightsHolder[0];

        // Create mixture with optimal weights
        MatrixCell mixture = mmap_mixture_order2(mmaps, bestWeights);

        return new Pair<MatrixCell, double[]>(mixture, bestWeights);
    }

    /**
     * Generate weight combinations for optimization.
     */
    private static void generateWeightCombinations(int n, int gridSize, WeightsCallback callback) {
        double[] weights = new double[n];
        generateRecursive(0, 1.0, n, gridSize, weights, callback);
    }

    private static void generateRecursive(int index, double remainingWeight, int n, int gridSize,
                                          double[] weights, WeightsCallback callback) {
        if (index == n - 1) {
            weights[index] = remainingWeight;
            callback.apply(weights);
            return;
        }

        double step = remainingWeight / (gridSize - index);
        for (int i = 0; i <= gridSize - index; i++) {
            weights[index] = i * step;
            generateRecursive(index + 1, remainingWeight - weights[index], n, gridSize, weights, callback);
        }
    }
}
