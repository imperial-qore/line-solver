package jline.lib.perm;

import jline.util.matrix.Matrix;

/**
 * Heuristic approximation to the permanent of a positive matrix.
 *
 * This implementation uses Sinkhorn scaling to make the matrix approximately
 * doubly stochastic, then applies a mean-field approximation with van der Waerden
 * bounds combined with a Gurvits-like capacity bound.
 *
 * The algorithm:
 * 1. Uses Sinkhorn scaling to make A approximately doubly stochastic
 * 2. Applies a mean-field approximation with van der Waerden bounds
 * 3. Combines with a Gurvits-like capacity bound
 * 4. Scales the result back to the original matrix scale
 *
 * This is a heuristic approximation suitable for large matrices where exact
 * computation is computationally prohibitive. The approximation quality depends
 * on the structure of the input matrix.
 */
public class HeuristicPermanent extends PermSolver {

    private final double tolerance;
    private final int maxIterations;

    /**
     * Constructor with all parameters.
     *
     * @param matrix         The matrix for which to compute the permanent (must be strictly positive)
     * @param tolerance      Convergence threshold for Sinkhorn scaling
     * @param maxIterations  Maximum number of Sinkhorn iterations
     * @param solve          Whether to automatically run solve() after construction
     * @throws IllegalArgumentException if matrix contains non-positive elements
     */
    public HeuristicPermanent(Matrix matrix, double tolerance, int maxIterations, boolean solve) {
        super(matrix);
        this.tolerance = tolerance;
        this.maxIterations = maxIterations;

        // Validate that all matrix elements are non-negative
        // We'll add a small epsilon to zeros during computation
        for (int i = 0; i < matrix.getNumRows(); i++) {
            for (int j = 0; j < matrix.getNumCols(); j++) {
                if (matrix.get(i, j) < 0.0) {
                    throw new IllegalArgumentException("Matrix must be non-negative. Found negative element at ("
                            + i + ", " + j + "): " + matrix.get(i, j));
                }
            }
        }

        if (solve) {
            solve();
        }
    }

    /**
     * Constructor with just matrix and solve parameters.
     * Uses default tolerance (1e-10) and maxIterations (1000).
     */
    public HeuristicPermanent(Matrix matrix, boolean solve) {
        this(matrix, 1e-10, 1000, solve);
    }

    /**
     * Constructor with just matrix parameter.
     * Uses default tolerance (1e-10), maxIterations (1000), and solve (false).
     */
    public HeuristicPermanent(Matrix matrix) {
        this(matrix, 1e-10, 1000, false);
    }

    @Override
    public void compute() {
        value = computeHeuristicPermanent();
    }

    /**
     * Computes the heuristic permanent approximation using Sinkhorn scaling
     * and mean-field approximation.
     *
     * @return The approximate permanent value
     */
    private double computeHeuristicPermanent() {
        // A zero used to be replaced by 1e-15 here. That is not invertible: it
        // changes the permanent by n!*eps, which is O(1) by n=18. Refuse instead.
        if (n > 0) {
            PermSupport.requireFullSupport(matrix, "heur");
        }
        Matrix workingMatrix = matrix.copy();

        // Sinkhorn scaling to make matrix approximately doubly stochastic
        SinkhornResult sinkhorn = sinkhornScaling(workingMatrix);
        Matrix B = sinkhorn.B;
        double[] r = sinkhorn.r;
        double[] c = sinkhorn.c;

        // Compute approximate permanent of doubly stochastic matrix B
        // Mean-field approximation: product of row sums / n^n * n!
        double[] rowSums = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += B.get(i, j);
            }
            rowSums[i] = sum;
        }

        double rowProd = 1.0;
        for (int i = 0; i < n; i++) {
            rowProd *= rowSums[i];
        }

        double pMeanfield = factorial(n) * (rowProd / Math.pow(n, n));

        // Gurvits-like capacity bound (optional refinement)
        double logSumRowSums = 0.0;
        for (int i = 0; i < n; i++) {
            logSumRowSums += Math.log(rowSums[i]);
        }
        double cap = Math.exp(logSumRowSums / n);
        double pGurvits = factorial(n) * Math.pow(cap / n, n);

        // Combine (simple average)
        double pEst = 0.5 * (pMeanfield + pGurvits);

        // Undo scaling
        double scaleFactor = 1.0;
        for (int i = 0; i < n; i++) {
            scaleFactor *= (1.0 / r[i]);
        }
        for (int j = 0; j < n; j++) {
            scaleFactor *= (1.0 / c[j]);
        }

        pEst *= scaleFactor;

        return pEst;
    }

    /**
     * Helper container for sinkhornScaling result (a 3-tuple).
     */
    private static final class SinkhornResult {
        final Matrix B;
        final double[] r;
        final double[] c;

        SinkhornResult(Matrix B, double[] r, double[] c) {
            this.B = B;
            this.r = r;
            this.c = c;
        }
    }

    /**
     * Performs Sinkhorn scaling to make the matrix approximately doubly stochastic.
     */
    private SinkhornResult sinkhornScaling(Matrix inputMatrix) {
        Matrix B = inputMatrix.copy();
        double[] r = new double[n];
        double[] c = new double[n];
        for (int i = 0; i < n; i++) { r[i] = 1.0; c[i] = 1.0; }

        boolean converged = false;
        double lastError = Double.POSITIVE_INFINITY;
        for (int iter = 0; iter < maxIterations; iter++) {
            // Update row scaling: r = 1 / (B * c)
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    sum += B.get(i, j) * c[j];
                }
                r[i] = 1.0 / sum;
            }

            // Update column scaling: c = 1 / (B' * r)
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int i = 0; i < n; i++) {
                    sum += B.get(i, j) * r[i];
                }
                c[j] = 1.0 / sum;
            }

            // Check convergence: max(abs(r .* (B * c) - 1)) < tol
            double maxDiff = 0.0;
            for (int i = 0; i < n; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < n; j++) {
                    rowSum += B.get(i, j) * c[j];
                }
                double diff = Math.abs(r[i] * rowSum - 1.0);
                if (diff > maxDiff) {
                    maxDiff = diff;
                }
            }

            if (maxDiff < tolerance) {
                converged = true;
                lastError = maxDiff;
                break;
            }
            lastError = maxDiff;
        }
        if (!converged) {
            throw new IllegalArgumentException(
                    "The Sinkhorn scaling did not converge to a doubly stochastic matrix in "
                    + maxIterations + " sweeps (margin error " + lastError
                    + " against a tolerance of " + tolerance + "). The estimate below assumes"
                    + " convergence, so no value is returned. The usual cause is a matrix"
                    + " without total support.");
        }

        // Apply scaling to matrix: B = diag(r) * B * diag(c)
        Matrix scaledB = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                scaledB.set(i, j, r[i] * B.get(i, j) * c[j]);
            }
        }

        return new SinkhornResult(scaledB, r, c);
    }

    /**
     * Computes n! (factorial of n).
     * Exact up to n = 170, the largest factorial representable as a double, which
     * keeps the estimate identical to the MATLAB twin perm_heur.m; beyond that
     * Stirling's approximation avoids the overflow.
     */
    private double factorial(int n) {
        if (n <= 1) return 1.0;
        if (n <= 170) {
            // For representable n, compute exactly
            double result = 1.0;
            for (int i = 2; i <= n; i++) {
                result *= i;
            }
            return result;
        } else {
            // For large n, use Stirling's approximation: n! ~ sqrt(2*pi*n) * (n/e)^n
            return Math.sqrt(2.0 * Math.PI * n) * Math.pow(n / Math.E, n);
        }
    }
}
