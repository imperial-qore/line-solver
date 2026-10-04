package jline.lib.perm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import jline.util.matrix.Matrix;

/**
 * Sampling method to approximate the permanent using the Huber-Law bound.
 *
 * The Huber-Law algorithm uses acceptance-rejection sampling with a carefully constructed
 * bound to estimate the permanent. The method rescales the matrix to be doubly stochastic
 * and then uses the Huber-Law bound for efficient sampling.
 */
public class HuberLawSampler extends PermSolver {
    private final double delta;
    private final double alpha2;
    private final double epsilon;
    private final String mode;
    private final int numberOfSamples;
    private final long maximumTime; // milliseconds

    /** Draw budget of the classic mode; exceeding it throws rather than spins. */
    private static final int MAX_DRAWS = 1000000;
    /** Sweep budget of the doubly stochastic rescaling; exceeding it throws. */
    private static final int MAX_SINKHORN = 10000;

    private final Random random = new Random();
    private Matrix C; // Rescaled doubly stochastic matrix
    private double rescalingConstant = 1.0;

    private List<Integer> sampleAccepted = Collections.emptyList();
    private List<Long> sampleTime = Collections.emptyList();
    private List<Double> permStep = Collections.emptyList();

    public HuberLawSampler(Matrix matrix) {
        this(matrix, 0.1, 0.000001, 0.1, "classic", 1000, 30000L, false);
    }

    public HuberLawSampler(Matrix matrix, double delta, double alpha2, double epsilon, String mode,
                           int numberOfSamples, long maximumTime, boolean solve) {
        super(matrix);
        if (matrix.getNumRows() > 0) {
            PermSupport.requireFullSupport(matrix, "huberlaw");
        }
        this.delta = delta;
        this.alpha2 = alpha2;
        this.epsilon = epsilon;
        this.mode = mode;
        this.numberOfSamples = numberOfSamples;
        this.maximumTime = maximumTime;
        if (solve) solve();
    }

    public List<Integer> getSampleAccepted() { return sampleAccepted; }
    public List<Long> getSampleTime() { return sampleTime; }
    public List<Double> getPermStep() { return permStep; }

    @Override
    public void compute() {
        if ("classic".equals(mode)) value = samplingClassic();
        else if ("time".equals(mode)) value = samplingTime();
        else if ("sample".equals(mode)) value = samplingSample();
        else value = samplingClassic();
    }

    private double samplingClassic() {
        rescale();
        int K = (int) (14.0 * Math.pow(delta, -2) * Math.log(2.0 / epsilon));
        long startTime = System.currentTimeMillis();
        List<Integer> acceptedList = new ArrayList<Integer>();
        List<Long> timeList = new ArrayList<Long>();
        int acceptedCount = 0;
        // Bounded independently of the scaling: with perm(A) = 0 the acceptance
        // probability is 0 and this loop would never terminate. A cap that
        // RETURNS a number would be a workaround, so it throws.
        while (acceptedCount < K) {
            if (acceptedList.size() >= MAX_DRAWS) {
                throw new IllegalArgumentException("Only " + acceptedCount + " of the " + K
                        + " required acceptances were obtained in " + acceptedList.size()
                        + " draws. The acceptance probability is too low for this budget;"
                        + " relax delta or use the exact engine.");
            }
            int[] sigma = sample();
            int isAccepted = anyEquals(sigma, n) ? 0 : 1;
            acceptedList.add(isAccepted);
            timeList.add(System.currentTimeMillis() - startTime);
            if (isAccepted == 1) acceptedCount++;
        }
        sampleAccepted = acceptedList;
        sampleTime = timeList;
        computePermStep();
        return sumInts(sampleAccepted) / (double) sampleAccepted.size() * rescalingConstant;
    }

    private double samplingTime() {
        rescale();
        long startTime = System.currentTimeMillis();
        List<Integer> acceptedList = new ArrayList<Integer>();
        List<Long> timeList = new ArrayList<Long>();
        while (System.currentTimeMillis() - startTime < maximumTime) {
            int[] sigma = sample();
            int isAccepted = anyEquals(sigma, n) ? 0 : 1;
            acceptedList.add(isAccepted);
            timeList.add(System.currentTimeMillis() - startTime);
        }
        sampleAccepted = acceptedList;
        sampleTime = timeList;
        computePermStep();
        return sampleAccepted.isEmpty() ? 0.0
                : sumInts(sampleAccepted) / (double) sampleAccepted.size() * rescalingConstant;
    }

    private double samplingSample() {
        rescale();
        long startTime = System.currentTimeMillis();
        List<Integer> acceptedList = new ArrayList<Integer>();
        List<Long> timeList = new ArrayList<Long>();
        while (acceptedList.size() < numberOfSamples) {
            int[] sigma = sample();
            int isAccepted = anyEquals(sigma, n) ? 0 : 1;
            acceptedList.add(isAccepted);
            timeList.add(System.currentTimeMillis() - startTime);
        }
        sampleAccepted = acceptedList;
        sampleTime = timeList;
        computePermStep();
        return sumInts(sampleAccepted) / (double) sampleAccepted.size() * rescalingConstant;
    }

    private int[] sample() {
        Matrix M = C.copy();
        int[] sigma = new int[n];
        for (int j = 0; j < n; j++) {
            double[] rowSums = new double[n];
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int k = 0; k < n; k++) s += M.get(i, k);
                rowSums[i] = s;
            }
            double ub = 1.0;
            for (int i = 0; i < n; i++) ub *= h(rowSums[i]);
            ub /= Math.exp(n);
            double[] p = precomputing(M, j);
            double[] normalizedP = new double[n];
            double probSum = 0.0;
            for (int i = 0; i < n; i++) {
                normalizedP[i] = p[i] / ub;
                probSum += normalizedP[i];
            }
            double[] prob = new double[n + 1];
            System.arraycopy(normalizedP, 0, prob, 0, n);
            prob[n] = 1.0 - probSum;
            if (prob[n] < 0) {
                double posSum = 0.0;
                for (int i = 0; i < n; i++) posSum += prob[i];
                for (int i = 0; i < n; i++) prob[i] /= posSum;
                prob[n] = 0.0;
            }
            double randVal = random.nextDouble();
            double cumSum = 0.0;
            int selectedI = n;
            for (int i = 0; i < prob.length; i++) {
                cumSum += prob[i];
                if (randVal <= cumSum) { selectedI = i; break; }
            }
            if (selectedI == n) {
                int[] rejected = new int[n];
                for (int i = 0; i < n; i++) rejected[i] = n;
                return rejected;
            }
            sigma[j] = selectedI;
            Matrix newMatrix = Matrix.zeros(n, n);
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    if ((row == selectedI && col == j) || (row != selectedI && col != j)) {
                        newMatrix.set(row, col, M.get(row, col));
                    }
                }
            }
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) M.set(row, col, newMatrix.get(row, col));
            }
        }
        return sigma;
    }

    private double h(double r) {
        if (r >= 1.0) {
            return r + 0.5 * Math.log(Math.max(r, 1.0)) + Math.E - 1;
        } else {
            return 1 + (Math.E - 1) * r;
        }
    }

    private double[] precomputing(Matrix M, int j) {
        double[] c = new double[n];
        double[] r = new double[n];
        for (int i = 0; i < n; i++) {
            c[i] = M.get(i, j);
            double sum = 0.0;
            for (int k = 0; k < n; k++) sum += M.get(i, k);
            r[i] = sum - c[i];
        }
        double[] hr = new double[n];
        double hrProduct = 1.0;
        for (int i = 0; i < n; i++) { hr[i] = h(r[i]); hrProduct *= hr[i]; }
        double[] result = new double[n];
        double expFactor = Math.exp(n - 1.0);
        for (int i = 0; i < n; i++) {
            result[i] = hrProduct / hr[i] * c[i] / expFactor;
        }
        return result;
    }

    private void rescale() {
        double[][] logMatrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            // The matrix is strictly positive here (requireFullSupport in the
            // constructor), so no log floor is needed.
            for (int j = 0; j < n; j++) logMatrix[i][j] = Math.log(matrix.get(i, j));
        }
        Matrix logMat = new Matrix(logMatrix);
        // A real maximum-weight assignment, not the row-by-row greedy that used
        // to stand in for it: the weight is alpha3, which sets the flooring
        // level alpha1 of the Huber-Law bound.
        int[] assignment = PermSupport.maxWeightAssignment(logMat);

        double maxElement = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (matrix.get(i, j) > maxElement) maxElement = matrix.get(i, j);
            }
        }
        double[][] MScaledArr = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) MScaledArr[i][j] = matrix.get(i, j) / maxElement;
        }
        Matrix MScaled = new Matrix(MScaledArr);

        // alpha3 is a permanent lower bound of the scaled matrix, the one floored below
        double alpha3 = 1.0;
        for (int i = 0; i < assignment.length; i++) alpha3 *= MScaled.get(i, assignment[i]);

        double alpha1 = alpha3 * delta / 3 / factorial(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (MScaled.get(i, j) < alpha1) MScaled.set(i, j, alpha1);
            }
        }

        Matrix[] dsResult = makeDoublyStochastic(MScaled);
        Matrix doublyStochastic = dsResult[0];
        Matrix X = dsResult[1];
        Matrix Y = dsResult[2];

        Matrix Z = Matrix.eye(n);
        for (int i = 0; i < n; i++) {
            double maxVal = 0.0;
            for (int j = 0; j < n; j++) if (doublyStochastic.get(i, j) > maxVal) maxVal = doublyStochastic.get(i, j);
            Z.set(i, i, 1.0 / maxVal);
        }
        C = Z.mult(doublyStochastic);

        double[] rowSums = new double[n];
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) s += C.get(i, j);
            rowSums[i] = s;
        }
        double hProduct = 1.0;
        for (int i = 0; i < n; i++) hProduct *= h(rowSums[i]) / Math.E;

        double diagonalProduct = 1.0;
        for (int i = 0; i < n; i++) diagonalProduct *= X.get(i, i) * Y.get(i, i) * Z.get(i, i);
        rescalingConstant = hProduct / diagonalProduct * Math.pow(maxElement, n);
    }

    private double factorial(int n) {
        double result = 1.0;
        for (int i = 1; i <= n; i++) result *= i;
        return result;
    }

    private Matrix[] makeDoublyStochastic(Matrix M) {
        Matrix result = M.copy();
        Matrix X = Matrix.eye(n);
        Matrix Y = Matrix.eye(n);
        double maxRowError = Double.POSITIVE_INFINITY;
        double maxColError = Double.POSITIVE_INFINITY;
        // Capped: a row that sums to zero leaves maxRowError at 1 forever and
        // the guarded normalization below skips it, so this loop used to spin
        // without terminating. A cap that RETURNS is a workaround; this throws.
        int sweeps = 0;
        while (maxRowError > alpha2 || maxColError > alpha2) {
            if (++sweeps > MAX_SINKHORN) {
                throw new IllegalArgumentException(
                        "The doubly stochastic rescaling did not converge in " + MAX_SINKHORN
                        + " sweeps (row error " + maxRowError + ", column error " + maxColError
                        + " against a tolerance of " + alpha2 + "). The usual cause is a matrix"
                        + " without total support.");
            }
            double[] colSums = new double[n];
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += result.get(i, j);
                colSums[j] = s;
            }
            for (int j = 0; j < n; j++) {
                if (colSums[j] > 0) {
                    for (int i = 0; i < n; i++) {
                        result.set(i, j, result.get(i, j) / colSums[j]);
                        Y.set(i, j, Y.get(i, j) / colSums[j]);
                    }
                }
            }
            double[] rowSums = new double[n];
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int j = 0; j < n; j++) s += result.get(i, j);
                rowSums[i] = s;
            }
            for (int i = 0; i < n; i++) {
                if (rowSums[i] > 0) {
                    for (int j = 0; j < n; j++) {
                        result.set(i, j, result.get(i, j) / rowSums[i]);
                        X.set(i, j, X.get(i, j) / rowSums[i]);
                    }
                }
            }
            double[] newColSums = new double[n];
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += result.get(i, j);
                newColSums[j] = s;
            }
            double[] newRowSums = new double[n];
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int j = 0; j < n; j++) s += result.get(i, j);
                newRowSums[i] = s;
            }
            maxColError = 0.0;
            for (double v : newColSums) maxColError = Math.max(maxColError, Math.abs(v - 1.0));
            maxRowError = 0.0;
            for (double v : newRowSums) maxRowError = Math.max(maxRowError, Math.abs(v - 1.0));
        }
        return new Matrix[]{result, X, Y};
    }

    private void computePermStep() {
        List<Double> steps = new ArrayList<Double>();
        int cumAccepted = 0;
        for (int i = 0; i < sampleAccepted.size(); i++) {
            cumAccepted += sampleAccepted.get(i);
            steps.add(rescalingConstant * cumAccepted / (i + 1));
        }
        permStep = steps;
    }

    private static int sumInts(List<Integer> list) {
        int s = 0;
        for (int v : list) s += v;
        return s;
    }

    private static boolean anyEquals(int[] arr, int v) {
        for (int x : arr) if (x == v) return true;
        return false;
    }
}
