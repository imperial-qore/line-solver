/**
 * @file BMAP/PH/N/N bufferless retrial queue analysis
 *
 * Implements analysis of BMAP/PH/N/N bufferless retrial queueing systems
 * with admission control. Uses the QBD-based algorithm from:
 * Dudin et al., "Analysis of BMAP/PH/N-Type Queueing System with Flexible
 * Retrials Admission Control", Mathematics 2025, 13(9), 1434.
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.ops.DConvertMatrixStruct;

import jline.GlobalConstants;
import jline.io.InputOutput;
import jline.lang.constant.RetrialPolicy;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixEntry;

public final class Qsys_bmapphnn_retrial {

    /** Default relative orbit-truncation error target. */
    public static final double DEFAULT_TAIL_TOLERANCE = 1e-6;

    /** Default cap on the total generator dimension explored by the adaptive refinement. */
    public static final double DEFAULT_MAX_DIM = 2e5;

    /** Default cap on the per-level block size V*d. */
    public static final int DEFAULT_MAX_BLOCK_SIZE = 5000;

    private Qsys_bmapphnn_retrial() {}

    public static QsysRetrialResult qsys_bmapphnn_retrial(
            Matrix[] D, Matrix beta, Matrix S, int N,
            double alpha, double gamma, double p, int R) {
        return qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, R, -1, 1e-10, false);
    }

    public static QsysRetrialResult qsys_bmapphnn_retrial(
            Matrix[] D, Matrix beta, Matrix S, int N,
            double alpha, double gamma, double p, int R, int maxLevel) {
        return qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, R, maxLevel, 1e-10, false);
    }

    public static QsysRetrialResult qsys_bmapphnn_retrial(
            Matrix[] D, Matrix beta, Matrix S, int N,
            double alpha, double gamma, double p, int R, int maxLevel, double tolerance, boolean verbose) {
        int[] Rarr = new int[D[0].getNumRows()];
        for (int i = 0; i < Rarr.length; i++) Rarr[i] = R;
        return qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, Rarr, maxLevel, tolerance, verbose);
    }

    public static QsysRetrialResult qsys_bmapphnn_retrial(
            Matrix[] D, Matrix beta, Matrix S, int N,
            double alpha, double gamma, double p, int R, int maxLevel, double tolerance, boolean verbose,
            double tailTolerance, double maxDim, int maxBlockSize) {
        int[] Rarr = new int[D[0].getNumRows()];
        for (int i = 0; i < Rarr.length; i++) Rarr[i] = R;
        return qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, Rarr, maxLevel, tolerance, verbose,
                tailTolerance, maxDim, maxBlockSize);
    }

    /**
     * As above with a scalar admission threshold and an explicit retrial policy.
     *
     * @param retrialPolicy {@link RetrialPolicy#LINEAR} or {@link RetrialPolicy#CONSTANT}
     * @return the retrial-queue result
     */
    public static QsysRetrialResult qsys_bmapphnn_retrial(
            Matrix[] D, Matrix beta, Matrix S, int N,
            double alpha, double gamma, double p, int R, int maxLevel, double tolerance, boolean verbose,
            double tailTolerance, double maxDim, int maxBlockSize, int retrialPolicy) {
        int[] Rarr = new int[D[0].getNumRows()];
        for (int i = 0; i < Rarr.length; i++) Rarr[i] = R;
        return qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, Rarr, maxLevel, tolerance, verbose,
                tailTolerance, maxDim, maxBlockSize, retrialPolicy);
    }

    public static QsysRetrialResult qsys_bmapphnn_retrial(
            Matrix[] D, Matrix beta, Matrix S, int N,
            double alpha, double gamma, double p, int[] R) {
        return qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, R, -1, 1e-10, false);
    }

    public static QsysRetrialResult qsys_bmapphnn_retrial(
            Matrix[] D, Matrix beta, Matrix S, int N,
            double alpha, double gamma, double p, int[] R, int maxLevel) {
        return qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, R, maxLevel, 1e-10, false);
    }

    public static QsysRetrialResult qsys_bmapphnn_retrial(
            Matrix[] D, Matrix beta, Matrix S, int N,
            double alpha, double gamma, double p, int[] R,
            int maxLevel, double tolerance, boolean verbose) {
        return qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, R, maxLevel, tolerance, verbose,
                DEFAULT_TAIL_TOLERANCE, DEFAULT_MAX_DIM, DEFAULT_MAX_BLOCK_SIZE);
    }

    /**
     * Analyzes a BMAP/PH/N/N bufferless retrial queue with admission control.
     *
     * @param maxLevel       fixed orbit truncation level; when non-positive the level is chosen
     *                       adaptively from the residual tail mass
     * @param tolerance      convergence tolerance of the underlying linear algebra
     * @param verbose        print progress messages
     * @param tailTolerance  relative orbit-truncation error target for the adaptive refinement
     * @param maxDim         cap on the total generator dimension explored by the refinement
     * @param maxBlockSize   cap on the per-level block size V*d
     */
    public static QsysRetrialResult qsys_bmapphnn_retrial(
            Matrix[] D, Matrix beta, Matrix S, int N,
            double alpha, double gamma, double p, int[] R,
            int maxLevel, double tolerance, boolean verbose,
            double tailTolerance, double maxDim, int maxBlockSize) {
        return qsys_bmapphnn_retrial(D, beta, S, N, alpha, gamma, p, R, maxLevel, tolerance, verbose,
                tailTolerance, maxDim, maxBlockSize, RetrialPolicy.LINEAR);
    }

    /**
     * As above, with an explicit retrial policy.
     *
     * @param retrialPolicy {@link RetrialPolicy#LINEAR} (each orbiting customer retries at
     *                      its own rate, aggregate rate n*nu) or {@link RetrialPolicy#CONSTANT}
     *                      (the orbit retries as a whole at rate nu when non-empty)
     * @return the retrial-queue result
     */
    public static QsysRetrialResult qsys_bmapphnn_retrial(
            Matrix[] D, Matrix beta, Matrix S, int N,
            double alpha, double gamma, double p, int[] R,
            int maxLevel, double tolerance, boolean verbose,
            double tailTolerance, double maxDim, int maxBlockSize, int retrialPolicy) {

        // see _kb/03-api-layer.md for rationale
        validateRetrialInputs(D, beta, S, N, alpha, gamma, p, R);

        // BMAP parameters
        int K = D.length - 1;
        int V = D[0].getNumRows();
        int M = S.getNumRows();

        // Compute generator of fundamental process: D^(1) = sum(D_k)
        Matrix D1_gen = new Matrix(V, V);
        for (int k = 0; k < D.length; k++) {
            for (int i = 0; i < V; i++) {
                for (int j = 0; j < V; j++) {
                    D1_gen.set(i, j, D1_gen.get(i, j) + D[k].get(i, j));
                }
            }
        }

        // Stationary distribution of fundamental process
        double[] theta = computeStationaryVector(D1_gen);

        // Mean arrival rate: lambda = theta * sum(k * D_k) * e
        Matrix sumKDk = new Matrix(V, V);
        for (int k = 1; k < D.length; k++) {
            for (int i = 0; i < V; i++) {
                for (int j = 0; j < V; j++) {
                    sumKDk.set(i, j, sumKDk.get(i, j) + ((double) k) * D[k].get(i, j));
                }
            }
        }
        double lambda = 0.0;
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                lambda += theta[j] * sumKDk.get(j, i);
            }
        }

        // PH service parameters
        Matrix S0 = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < M; j++) {
                rowSum += S.get(i, j);
            }
            S0.set(i, 0, -rowSum);
        }

        // Mean service time: b1 = beta * (-S)^{-1} * ones
        Matrix negSinv = S.scale(-1.0).inv();
        Matrix onesM = Matrix.ones(M, 1);
        double b1 = beta.mult(negSinv).mult(onesM).get(0, 0);

        // Compute T_n values
        int[] T = new int[N + 1];
        for (int n = 0; n <= N; n++) {
            T[n] = (int) CombinatoricsUtils.binomialCoefficient(n + M - 1, M - 1);
        }
        int d = 0;
        for (int n = 0; n <= N; n++) {
            d += T[n];
        }

        // Build state mapping
        int[][][] stateMap = buildStateMap(N, M, T);

        // Determine truncation level
        double rho = lambda * b1 / N;

        if (verbose) {
            System.out.println("Solving BMAP/PH/N/N retrial queue...");
            System.out.println("  V=" + V + ", M=" + M + ", N=" + N + ", K=" + K);
            System.out.println("  d=" + d + ", block size Vd=" + (V * d));
            System.out.println("  lambda=" + String.format("%.4f", lambda) + ", mu=" + String.format("%.4f", 1.0 / b1));
            System.out.println("  Offered load rho=" + String.format("%.4f", rho));
        }

        // Context for helper functions
        RetrialContext ctx = new RetrialContext(D, beta, S, S0, M, N, V, K, d, T, R, alpha, gamma, p, stateMap, retrialPolicy);

        int Vd = V * d;

        // see _kb/03-api-layer.md for rationale
        if (Vd > maxBlockSize) {
            InputOutput.line_error("qsys_bmapphnn_retrial", String.format(
                    "Per-level block size V*d = %d exceeds MaxBlockSize = %d. "
                    + "The service distribution has %d phases and the station has %d servers, which yields "
                    + "%d service configurations. Reduce the phase-type order (e.g. fit the service "
                    + "distribution with fewer phases), reduce the number of servers, or raise "
                    + "'MaxBlockSize' if the memory cost is acceptable.", Vd, maxBlockSize, M, N, d));
        }

        // see _kb/03-api-layer.md for rationale
        int truncLevel;
        double truncError;
        Matrix piMatrix;
        if (maxLevel > 0) {
            truncLevel = maxLevel;
            piMatrix = solveAtLevel(ctx, truncLevel, Vd, verbose);
            truncError = orbitTruncationError(piMatrix, truncLevel);
        } else {
            truncLevel = Math.max(100, (int) Math.ceil(50.0 / (1.0 - Math.min(rho, 0.99))));
            truncError = Double.POSITIVE_INFINITY;
            boolean converged = false;
            piMatrix = null;
            while (true) {
                piMatrix = solveAtLevel(ctx, truncLevel, Vd, verbose);
                truncError = orbitTruncationError(piMatrix, truncLevel);
                if (truncError <= tailTolerance) {
                    converged = true;
                    break;
                }
                int nextLevel = 2 * truncLevel;
                if (((double) nextLevel + 1.0) * Vd > maxDim) {
                    break;
                }
                if (verbose) {
                    System.out.println("  Truncation error " + String.format("%.3e", truncError)
                            + " > " + String.format("%.3e", tailTolerance) + ", refining to level " + nextLevel);
                }
                truncLevel = nextLevel;
            }
            if (!converged) {
                InputOutput.line_warning("qsys_bmapphnn_retrial", String.format(
                        "Orbit truncation did not reach the requested accuracy: "
                        + "residual %.3e > TailTolerance %.3e at level %d (dimension cap MaxDim = %d). "
                        + "Orbit measures are underestimated; raise 'MaxDim' or set 'MaxLevel' explicitly.",
                        truncError, tailTolerance, truncLevel, (long) Math.round(maxDim)));
            }
        }

        if (verbose) {
            System.out.println("  Truncation level: " + truncLevel
                    + " (residual " + String.format("%.3e", truncError) + ")");
        }

        // Compute performance measures

        double L_orbit = 0.0;
        for (int i = 1; i <= truncLevel; i++) {
            double levelProb = 0.0;
            for (int j = 0; j < Vd; j++) {
                levelProb += piMatrix.get(i, j);
            }
            L_orbit += ((double) i) * levelProb;
        }

        double N_server = 0.0;
        for (int i = 0; i <= truncLevel; i++) {
            for (int nu = 0; nu < V; nu++) {
                for (int n = 0; n <= N; n++) {
                    int offset = nu * d + getBlockOffset(T, n);
                    for (int t = 0; t < T[n]; t++) {
                        int idx = offset + t;
                        if (idx < Vd) {
                            N_server += ((double) n) * piMatrix.get(i, idx);
                        }
                    }
                }
            }
        }

        double P_idle = 0.0;
        for (int i = 0; i <= truncLevel; i++) {
            for (int nu = 0; nu < V; nu++) {
                int offset = nu * d;
                P_idle += piMatrix.get(i, offset);
            }
        }

        double P_empty_orbit = 0.0;
        for (int j = 0; j < Vd; j++) {
            P_empty_orbit += piMatrix.get(0, j);
        }

        double P_empty = 0.0;
        for (int nu = 0; nu < V; nu++) {
            int offset = nu * d;
            P_empty += piMatrix.get(0, offset);
        }

        if (verbose) {
            System.out.println("Solution complete.");
        }

        return new QsysRetrialResult(
                L_orbit,
                N_server,
                L_orbit + N_server,
                N_server / N,
                N_server / b1,
                P_idle,
                P_empty_orbit,
                P_empty,
                piMatrix,
                truncLevel,
                truncError,
                "LINE:qsys_bmapphnn_retrial");
    }

    // ========== Helper classes and functions ==========

    private static class RetrialContext {
        final Matrix[] D;
        final Matrix beta;
        final Matrix S;
        final Matrix S0;
        final int M;
        final int N;
        final int V;
        final int K;
        final int d;
        final int[] T;
        final int[] R;
        final double alpha;
        final double gamma;
        final double p;
        final int[][][] stateMap;
        final int retrialPolicy;

        RetrialContext(Matrix[] D, Matrix beta, Matrix S, Matrix S0, int M, int N, int V, int K,
                       int d, int[] T, int[] R, double alpha, double gamma, double p, int[][][] stateMap) {
            this(D, beta, S, S0, M, N, V, K, d, T, R, alpha, gamma, p, stateMap, RetrialPolicy.LINEAR);
        }

        RetrialContext(Matrix[] D, Matrix beta, Matrix S, Matrix S0, int M, int N, int V, int K,
                       int d, int[] T, int[] R, double alpha, double gamma, double p, int[][][] stateMap,
                       int retrialPolicy) {
            this.retrialPolicy = retrialPolicy;
            this.D = D;
            this.beta = beta;
            this.S = S;
            this.S0 = S0;
            this.M = M;
            this.N = N;
            this.V = V;
            this.K = K;
            this.d = d;
            this.T = T;
            this.R = R;
            this.alpha = alpha;
            this.gamma = gamma;
            this.p = p;
            this.stateMap = stateMap;
        }

        /**
         * Same context with the retrial rate and the orbit impatience rate replaced.
         * Used to isolate the two orbit terms, which scale differently with the orbit
         * size: impatience is always proportional to it, the retrial rate is weighted
         * by the policy.
         */
        RetrialContext withRates(double newAlpha, double newGamma) {
            return new RetrialContext(this.D, this.beta, this.S, this.S0, this.M, this.N, this.V,
                    this.K, this.d, this.T, this.R, newAlpha, newGamma, this.p, this.stateMap,
                    this.retrialPolicy);
        }
    }

    /**
     * Rejects inputs that are not a well-formed BMAP/PH pair. The generator build
     * is driven entirely by these matrices, so a NaN, an Inf, or a non-square
     * block would otherwise be written into the generator and only surface as a
     * meaningless stationary vector or as an out-of-memory failure.
     */
    private static void validateRetrialInputs(Matrix[] D, Matrix beta, Matrix S, int N,
                                              double alpha, double gamma, double p, int[] R) {
        String mfilename = "qsys_bmapphnn_retrial";

        if (D == null || D.length == 0) {
            InputOutput.line_error(mfilename,
                    "BMAP arrival representation D must be a non-empty cell array {D0,D1,...}.");
        }
        if (D.length < 2) {
            InputOutput.line_error(mfilename,
                    "BMAP arrival representation D must contain at least {D0,D1}. "
                    + "A single-element representation is a non-Markovian distribution (e.g. Det or a trace) "
                    + "and is not admissible in the matrix-analytic retrial engine.");
        }
        int V = D[0].getNumRows();
        for (int k = 0; k < D.length; k++) {
            Matrix Dk = D[k];
            if (Dk == null || Dk.getNumRows() != Dk.getNumCols() || Dk.getNumRows() != V) {
                InputOutput.line_error(mfilename, String.format(
                        "BMAP matrix D{%d} must be a %dx%d numeric matrix.", k + 1, V, V));
            }
            if (Dk.hasNaN() || Dk.hasInfinite()) {
                InputOutput.line_error(mfilename, String.format(
                        "BMAP matrix D{%d} contains NaN or Inf entries. The arrival "
                        + "process is disabled or not phase-type representable.", k + 1));
            }
            if (k >= 1 && minEntry(Dk) < -GlobalConstants.FineTol) {
                InputOutput.line_error(mfilename, String.format(
                        "BMAP arrival matrix D{%d} must be non-negative.", k + 1));
            }
        }
        for (int i = 0; i < V; i++) {
            if (D[0].get(i, i) > GlobalConstants.FineTol) {
                InputOutput.line_error(mfilename, "BMAP matrix D0 must have non-positive diagonal entries.");
            }
        }
        for (int i = 0; i < V; i++) {
            double rowSum = 0.0;
            for (int k = 0; k < D.length; k++) {
                for (int j = 0; j < V; j++) {
                    rowSum += D[k].get(i, j);
                }
            }
            if (Math.abs(rowSum) > Math.sqrt(GlobalConstants.FineTol)) {
                InputOutput.line_error(mfilename,
                        "BMAP matrices are inconsistent: sum_k D_k must have zero row sums.");
            }
        }

        int nbeta = (beta == null) ? 0 : beta.getNumRows() * beta.getNumCols();
        if (beta == null || nbeta == 0 || beta.hasNaN() || beta.hasInfinite()) {
            InputOutput.line_error(mfilename,
                    "Phase-type service vector beta is empty or contains NaN/Inf. The service "
                    + "distribution is disabled or not phase-type representable.");
        }
        int M = (S == null) ? 0 : S.getNumRows();
        if (S == null || S.getNumCols() != M || nbeta != M) {
            InputOutput.line_error(mfilename, String.format(
                    "Phase-type service subgenerator S must be square and conformant with beta (%d phases).", nbeta));
        }
        if (S.hasNaN() || S.hasInfinite()) {
            InputOutput.line_error(mfilename,
                    "Phase-type service subgenerator S contains NaN or Inf entries. The service "
                    + "distribution is disabled or not phase-type representable.");
        }
        for (int i = 0; i < M; i++) {
            if (S.get(i, i) >= 0) {
                InputOutput.line_error(mfilename,
                        "Phase-type service subgenerator S must have strictly negative diagonal entries.");
            }
        }
        double betaSum = 0.0;
        for (int i = 0; i < nbeta; i++) {
            double bi = beta.get(i);
            betaSum += bi;
            if (bi < -GlobalConstants.FineTol) {
                InputOutput.line_error(mfilename,
                        "Phase-type service vector beta must be non-negative and sum to one.");
            }
        }
        if (Math.abs(betaSum - 1.0) > Math.sqrt(GlobalConstants.FineTol)) {
            InputOutput.line_error(mfilename,
                    "Phase-type service vector beta must be non-negative and sum to one.");
        }
        for (int i = 0; i < M; i++) {
            double exitRate = 0.0;
            for (int j = 0; j < M; j++) {
                exitRate -= S.get(i, j);
            }
            if (exitRate < -GlobalConstants.FineTol) {
                InputOutput.line_error(mfilename,
                        "Phase-type service subgenerator S must have non-negative exit rates.");
            }
        }

        if (N < 1) {
            InputOutput.line_error(mfilename, "Number of servers N must be a positive integer.");
        }
        if (!isFiniteValue(alpha) || alpha < 0) {
            InputOutput.line_error(mfilename, "Retrial rate alpha must be a finite non-negative scalar.");
        }
        if (!isFiniteValue(gamma) || gamma < 0) {
            InputOutput.line_error(mfilename, "Orbit impatience rate gamma must be a finite non-negative scalar.");
        }
        if (!isFiniteValue(p) || p < 0 || p > 1) {
            InputOutput.line_error(mfilename, "Batch rejection probability p must lie in [0,1].");
        }
        if (R == null) {
            InputOutput.line_error(mfilename, String.format("Admission threshold R must lie in [0,%d].", N));
        }
        for (int i = 0; i < R.length; i++) {
            if (R[i] < 0 || R[i] > N) {
                InputOutput.line_error(mfilename, String.format("Admission threshold R must lie in [0,%d].", N));
            }
        }
    }

    private static boolean isFiniteValue(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }

    private static double minEntry(Matrix A) {
        double min = 0.0;
        Iterator<MatrixEntry> it = A.nonZeroIterator();
        while (it.hasNext()) {
            MatrixEntry e = it.next();
            if (e.value < min) {
                min = e.value;
            }
        }
        return min;
    }

    /**
     * Relative contribution that the truncated tail would add to the mean orbit
     * length. Truncation reflects the probability flow that would leave the top
     * level back into it, so the mass sitting at the top level bounds the error.
     */
    private static double orbitTruncationError(Matrix pi, int truncLevel) {
        int nCols = pi.getNumCols();
        double lOrbit = 0.0;
        double lastMass = 0.0;
        for (int i = 0; i <= truncLevel; i++) {
            double mass = 0.0;
            for (int j = 0; j < nCols; j++) {
                mass += pi.get(i, j);
            }
            lOrbit += ((double) i) * mass;
            if (i == truncLevel) {
                lastMass = mass;
            }
        }
        return truncLevel * lastMass / Math.max(lOrbit, Double.MIN_NORMAL);
    }

    /**
     * Builds the level-truncated generator and solves pi*Q = 0, pi*e = 1.
     *
     * <p>The level blocks are level-homogeneous apart from the orbit terms, which
     * are linear in the level index: the diagonal block is Qdiag0 + i*Qdiag1
     * (retrial and impatience departures from an orbit of size i), the
     * subdiagonal block is i*Qsub1 (one of the i orbiting customers succeeds or
     * abandons) and the k-th superdiagonal block is level-independent. Building
     * those four shapes once and replicating them keeps the assembly linear in
     * the truncation level, which the adaptive refinement relies on.</p>
     */
    private static Matrix solveAtLevel(RetrialContext ctx, int truncLevel, int Vd, boolean verbose) {
        long totalDimLong = ((long) truncLevel + 1L) * (long) Vd;
        if (totalDimLong > Integer.MAX_VALUE) {
            InputOutput.line_error("qsys_bmapphnn_retrial", String.format(
                    "Generator dimension %d exceeds the addressable matrix size; lower the truncation level.",
                    totalDimLong));
        }
        int totalDim = (int) totalDimLong;

        if (verbose) {
            System.out.println("Total matrix dimension: " + totalDim + " x " + totalDim);
        }

        // see _kb/03-api-layer.md for rationale
        RetrialContext ctxGamma = ctx.withRates(0.0, ctx.gamma);   // impatience only
        RetrialContext ctxAlpha = ctx.withRates(ctx.alpha, 0.0);   // retrials only

        Matrix Qdiag0 = buildGeneratorLevel(ctx, 0, 0);              // diagonal block, empty orbit
        Matrix Qdiag1G = buildGeneratorLevel(ctxGamma, 1, 1).add(-1.0, Qdiag0); // per-customer impatience
        Matrix Qdiag1A = buildGeneratorLevel(ctxAlpha, 1, 1).add(-1.0, Qdiag0); // one-unit retrial
        Matrix QsubG = buildGeneratorLevel(ctxGamma, 1, 0);          // subdiagonal, impatience part
        Matrix QsubA = buildGeneratorLevel(ctxAlpha, 1, 0);          // subdiagonal, retrial part
        Matrix[] Qsup = new Matrix[ctx.K + 1];
        for (int k = 1; k <= ctx.K; k++) {
            Qsup[k] = buildGeneratorLevel(ctx, 0, k);
        }

        // see _kb/03-api-layer.md for rationale
        BlockPattern diag = mergePattern(new Matrix[] { Qdiag0, Qdiag1G, Qdiag1A }, true, Vd);
        BlockPattern sub = mergePattern(new Matrix[] { null, QsubG, QsubA }, false, Vd);
        BlockPattern[] sup = new BlockPattern[ctx.K + 1];
        for (int k = 1; k <= ctx.K; k++) {
            sup[k] = mergePattern(new Matrix[] { Qsup[k] }, false, Vd);
        }

        long nnzEstimate = (long) (truncLevel + 1) * diag.size()
                + (long) truncLevel * sub.size() + totalDim;
        for (int k = 1; k <= ctx.K; k++) {
            nnzEstimate += (long) Math.max(0, truncLevel - k + 1) * sup[k].size();
        }
        if (nnzEstimate > Integer.MAX_VALUE) {
            InputOutput.line_error("qsys_bmapphnn_retrial", String.format(
                    "Generator has %d non-zeros, exceeding the addressable matrix size; lower the truncation level.",
                    nnzEstimate));
        }

        DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(totalDim, totalDim, (int) nnzEstimate);
        int lastCol = totalDim - 1;

        for (int i = 0; i <= truncLevel; i++) {
            int rowBase = i * Vd;

            // Retrial weight of this level: the orbit size under LINEAR, one whenever
            // the orbit is non-empty under CONSTANT.
            double retrialWeight;
            if (ctx.retrialPolicy == RetrialPolicy.CONSTANT) {
                retrialWeight = (i >= 1) ? 1.0 : 0.0;
            } else {
                retrialWeight = i;
            }

            // Row sums of the full level-i block row, used to force zero row sums.
            double[] rowSum = new double[Vd];
            for (int r = 0; r < Vd; r++) {
                rowSum[r] = diag.rowSum0[r] + ((double) i) * diag.rowSum1[r]
                        + retrialWeight * diag.rowSum2[r];
                if (i >= 1) {
                    rowSum[r] += ((double) i) * sub.rowSum1[r] + retrialWeight * sub.rowSum2[r];
                }
            }
            for (int k = 1; k <= ctx.K; k++) {
                if (i + k <= truncLevel) {
                    for (int r = 0; r < Vd; r++) {
                        rowSum[r] += sup[k].rowSum0[r];
                    }
                }
            }

            // Diagonal block
            for (int e = 0; e < diag.size(); e++) {
                int r = diag.rows[e];
                int c = diag.cols[e];
                double v = diag.vals0[e] + ((double) i) * diag.vals1[e]
                        + retrialWeight * diag.vals2[e];
                if (r == c) {
                    v -= rowSum[r];
                }
                int gc = rowBase + c;
                if (v != 0.0 && gc != lastCol) {
                    triplet.addItem(rowBase + r, gc, v);
                }
            }

            // Subdiagonal block, rate proportional to the orbit size
            if (i >= 1) {
                int colBase = (i - 1) * Vd;
                for (int e = 0; e < sub.size(); e++) {
                    double v = ((double) i) * sub.vals1[e] + retrialWeight * sub.vals2[e];
                    int gc = colBase + sub.cols[e];
                    if (v != 0.0 && gc != lastCol) {
                        triplet.addItem(rowBase + sub.rows[e], gc, v);
                    }
                }
            }

            // Superdiagonal blocks, level independent
            for (int k = 1; k <= ctx.K; k++) {
                if (i + k > truncLevel) {
                    continue;
                }
                int colBase = (i + k) * Vd;
                for (int e = 0; e < sup[k].size(); e++) {
                    double v = sup[k].vals0[e];
                    int gc = colBase + sup[k].cols[e];
                    if (v != 0.0 && gc != lastCol) {
                        triplet.addItem(rowBase + sup[k].rows[e], gc, v);
                    }
                }
            }
        }

        // Replace the last column with ones for normalization
        for (int i = 0; i < totalDim; i++) {
            triplet.addItem(i, lastCol, 1.0);
        }

        DMatrixSparseCSC Qcsc = DConvertMatrixStruct.convert(triplet, (DMatrixSparseCSC) null);
        Matrix Q = new Matrix((org.ejml.data.DMatrix) Qcsc);

        // Solve pi * Q = 0, pi * e = 1
        if (verbose) {
            System.out.println("Solving linear system...");
        }

        Matrix Qt = Q.transpose();
        Matrix b = new Matrix(totalDim, 1, 1);
        b.set(totalDim - 1, 0, 1.0);
        Matrix piCol = new Matrix(totalDim, 1, totalDim);

        boolean ok = Matrix.solveDirect(Qt, b, piCol);
        if (!ok || piCol.hasNaN() || piCol.hasInfinite()) {
            // The sparse factorization failed; fall back to the checked dense solve,
            // which is affordable only for a moderate truncation level.
            if (totalDim > DEFAULT_MAX_BLOCK_SIZE) {
                InputOutput.line_error("qsys_bmapphnn_retrial", String.format(
                        "The generator of dimension %d is singular to working precision at truncation level %d.",
                        totalDim, truncLevel));
            }
            piCol = new Matrix(totalDim, 1, totalDim);
            Matrix.solveSafe(Qt.copy().toDense(), b.copy().toDense(), piCol);
        }

        // Read the solution once: element-wise access on a sparse column vector is
        // a per-call column scan, which dominates the solve at high truncation levels.
        double[] piVec = new double[totalDim];
        if (piCol.isSparse()) {
            Iterator<MatrixEntry> it = piCol.nonZeroIterator();
            while (it.hasNext()) {
                MatrixEntry e = it.next();
                piVec[e.row] = e.value;
            }
        } else {
            for (int i = 0; i < totalDim; i++) {
                piVec[i] = piCol.get(i, 0);
            }
        }

        // Reshape to level structure
        DMatrixRMaj piDense = new DMatrixRMaj(truncLevel + 1, Vd);
        boolean hasNeg = false;
        double piSum = 0.0;
        for (int i = 0; i <= truncLevel; i++) {
            for (int j = 0; j < Vd; j++) {
                double v = piVec[i * Vd + j];
                if (v < -1e-8) {
                    hasNeg = true;
                }
                if (v < 0) {
                    v = 0.0;
                }
                piDense.set(i, j, v);
                piSum += v;
            }
        }
        if (hasNeg) {
            InputOutput.line_warning("qsys_bmapphnn_retrial", "Negative probabilities detected, clipping to zero");
        }
        if (piSum > 0) {
            for (int i = 0; i <= truncLevel; i++) {
                for (int j = 0; j < Vd; j++) {
                    piDense.set(i, j, piDense.get(i, j) / piSum);
                }
            }
        }
        return new Matrix(piDense);
    }

    /**
     * Sparsity pattern of a level block, split into the three ways an entry can
     * depend on the orbit size i: a level-independent part (vals0), a part
     * proportional to i (vals1, orbit impatience, since each orbiting customer
     * abandons on its own clock), and the retrial part (vals2), whose weight is
     * i under RetrialPolicy.LINEAR and 1 when the orbit is non-empty under
     * RetrialPolicy.CONSTANT.
     */
    private static class BlockPattern {
        final int[] rows;
        final int[] cols;
        final double[] vals0;
        final double[] vals1;
        final double[] vals2;
        final double[] rowSum0;
        final double[] rowSum1;
        final double[] rowSum2;

        BlockPattern(int[] rows, int[] cols, double[] vals0, double[] vals1, double[] vals2, int n) {
            this.rows = rows;
            this.cols = cols;
            this.vals0 = vals0;
            this.vals1 = vals1;
            this.vals2 = vals2;
            this.rowSum0 = new double[n];
            this.rowSum1 = new double[n];
            this.rowSum2 = new double[n];
            for (int e = 0; e < rows.length; e++) {
                this.rowSum0[rows[e]] += vals0[e];
                this.rowSum1[rows[e]] += vals1[e];
                this.rowSum2[rows[e]] += vals2[e];
            }
        }

        int size() {
            return rows.length;
        }
    }

    /**
     * Merges up to three level-block shapes into a single sparsity pattern, so that
     * the block at orbit level i is vals0 + w1*vals1 + w2*vals2 without duplicated
     * triplet entries.
     *
     * @param mats        the shapes, in slot order (constant, impatience, retrial);
     *                    a null slot contributes nothing
     * @param reserveDiag reserve every diagonal position, needed for the block that
     *                    receives the row-sum correction
     * @param n           the block size
     * @return the merged pattern
     */
    private static BlockPattern mergePattern(Matrix[] mats, boolean reserveDiag, int n) {
        Map<Integer, Integer> keyToIdx = new HashMap<Integer, Integer>();
        ArrayList<int[]> idx = new ArrayList<int[]>();
        ArrayList<double[]> vals = new ArrayList<double[]>();

        if (reserveDiag) {
            for (int r = 0; r < n; r++) {
                keyToIdx.put(Integer.valueOf(r * n + r), Integer.valueOf(idx.size()));
                idx.add(new int[] { r, r });
                vals.add(new double[3]);
            }
        }

        for (int which = 0; which < mats.length; which++) {
            Matrix A = mats[which];
            if (A == null) {
                continue;
            }
            Iterator<MatrixEntry> it = A.nonZeroIterator();
            while (it.hasNext()) {
                MatrixEntry e = it.next();
                if (e.value == 0.0) {
                    continue;
                }
                Integer key = Integer.valueOf(e.row * n + e.col);
                Integer at = keyToIdx.get(key);
                if (at == null) {
                    keyToIdx.put(key, Integer.valueOf(idx.size()));
                    double[] triple = new double[3];
                    triple[which] = e.value;
                    idx.add(new int[] { e.row, e.col });
                    vals.add(triple);
                } else {
                    vals.get(at.intValue())[which] += e.value;
                }
            }
        }

        int m = idx.size();
        int[] rows = new int[m];
        int[] cols = new int[m];
        double[] v0 = new double[m];
        double[] v1 = new double[m];
        double[] v2 = new double[m];
        for (int e = 0; e < m; e++) {
            rows[e] = idx.get(e)[0];
            cols[e] = idx.get(e)[1];
            v0[e] = vals.get(e)[0];
            v1[e] = vals.get(e)[1];
            v2[e] = vals.get(e)[2];
        }
        return new BlockPattern(rows, cols, v0, v1, v2, n);
    }

    private static double[] computeStationaryVector(Matrix Q) {
        int n = Q.getNumRows();
        Matrix A = Q.transpose();
        for (int j = 0; j < n; j++) {
            A.set(n - 1, j, 1.0);
        }
        Matrix b = new Matrix(n, 1);
        b.set(n - 1, 0, 1.0);
        Matrix x = new Matrix(n, 1);
        Matrix.solveSafe(A, b, x);
        double[] theta = new double[n];
        for (int i = 0; i < n; i++) {
            theta[i] = x.get(i, 0);
        }
        return theta;
    }

    private static int[][][] buildStateMap(int N, int M, int[] T) {
        int[][][] result = new int[N + 1][][];
        for (int n = 0; n <= N; n++) {
            result[n] = generateCompositions(n, M);
        }
        return result;
    }

    private static int[][] generateCompositions(int n, int M) {
        if (M == 1) {
            return new int[][] { new int[] { n } };
        }
        ArrayList<int[]> result = new ArrayList<int[]>();
        for (int m1 = n; m1 >= 0; m1--) {
            int[][] subComps = generateCompositions(n - m1, M - 1);
            for (int[] sub : subComps) {
                int[] comp = new int[M];
                comp[0] = m1;
                System.arraycopy(sub, 0, comp, 1, M - 1);
                result.add(comp);
            }
        }
        return result.toArray(new int[result.size()][]);
    }

    private static int getBlockOffset(int[] T, int n) {
        int offset = 0;
        for (int i = 0; i < n; i++) {
            offset += T[i];
        }
        return offset;
    }

    private static boolean intArrayEquals(int[] a, int[] b) {
        if (a.length != b.length) return false;
        for (int i = 0; i < a.length; i++) if (a[i] != b[i]) return false;
        return true;
    }

    private static Matrix computeL(RetrialContext ctx, int n) {
        if (n == 0) return null;
        int rows = ctx.T[n];
        int cols = ctx.T[n - 1];
        Matrix L = new Matrix(rows, cols);
        int[][] compsN = ctx.stateMap[n];
        int[][] compsNm1 = ctx.stateMap[n - 1];

        for (int i = 0; i < compsN.length; i++) {
            int[] m = compsN[i];
            for (int l = 0; l < ctx.M; l++) {
                if (m[l] > 0) {
                    int[] mPrime = m.clone();
                    mPrime[l] = mPrime[l] - 1;
                    for (int j = 0; j < compsNm1.length; j++) {
                        if (intArrayEquals(compsNm1[j], mPrime)) {
                            L.set(i, j, L.get(i, j) + ((double) m[l]) * ctx.S0.get(l, 0));
                            break;
                        }
                    }
                }
            }
        }
        return L;
    }

    private static Matrix computeA(RetrialContext ctx, int n) {
        if (n == 0) return null;
        int sz = ctx.T[n];
        Matrix A = new Matrix(sz, sz);
        int[][] comps = ctx.stateMap[n];

        for (int i = 0; i < comps.length; i++) {
            int[] m = comps[i];
            for (int l = 0; l < ctx.M; l++) {
                if (m[l] > 0) {
                    for (int lPrime = 0; lPrime < ctx.M; lPrime++) {
                        if (lPrime != l && ctx.S.get(l, lPrime) > 0) {
                            int[] mPrime = m.clone();
                            mPrime[l] = mPrime[l] - 1;
                            mPrime[lPrime] = mPrime[lPrime] + 1;
                            for (int j = 0; j < comps.length; j++) {
                                if (intArrayEquals(comps[j], mPrime)) {
                                    A.set(i, j, A.get(i, j) + ((double) m[l]) * ctx.S.get(l, lPrime));
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        return A;
    }

    private static Matrix computeP(RetrialContext ctx, int n) {
        if (n >= ctx.N) return null;
        int rows = ctx.T[n];
        int cols = ctx.T[n + 1];
        Matrix P = new Matrix(rows, cols);
        int[][] compsN = ctx.stateMap[n];
        int[][] compsNp1 = ctx.stateMap[n + 1];

        for (int i = 0; i < compsN.length; i++) {
            int[] m = compsN[i];
            for (int l = 0; l < ctx.M; l++) {
                double betaL = (ctx.beta.getNumCols() > 1) ? ctx.beta.get(0, l) : ctx.beta.get(l, 0);
                if (betaL > 0) {
                    int[] mPrime = m.clone();
                    mPrime[l] = mPrime[l] + 1;
                    for (int j = 0; j < compsNp1.length; j++) {
                        if (intArrayEquals(compsNp1[j], mPrime)) {
                            P.set(i, j, P.get(i, j) + betaL);
                            break;
                        }
                    }
                }
            }
        }
        return P;
    }

    private static Matrix computeDelta(RetrialContext ctx, int n) {
        if (n == 0) return null;
        int sz = ctx.T[n];
        Matrix Delta = new Matrix(sz, sz);
        int[][] comps = ctx.stateMap[n];

        for (int i = 0; i < comps.length; i++) {
            int[] m = comps[i];
            double total = 0.0;
            for (int l = 0; l < ctx.M; l++) {
                total += ((double) m[l]) * (-ctx.S.get(l, l));
            }
            Delta.set(i, i, total);
        }
        return Delta;
    }

    private static Matrix computeGamma(RetrialContext ctx, int nu) {
        Matrix Gam = new Matrix(ctx.d, ctx.d);
        int offset = 0;
        for (int n = 0; n <= ctx.N; n++) {
            if (n > ctx.R[nu]) {
                for (int t = 0; t < ctx.T[n]; t++) {
                    Gam.set(offset + t, offset + t, 1.0);
                }
            }
            offset += ctx.T[n];
        }
        return Gam;
    }

    private static Matrix computeG_nn(RetrialContext ctx, int n, int nu, int nuPrime) {
        int sz = ctx.T[n];
        Matrix G = new Matrix(sz, sz);
        if (n <= ctx.N - ctx.K) {
            return G;
        }
        double total = 0.0;
        for (int k = (ctx.N - n + 1); k <= ctx.K; k++) {
            if (k >= 1 && k <= ctx.K) {
                total += ctx.D[k].get(nu, nuPrime);
            }
        }
        for (int i = 0; i < sz; i++) {
            G.set(i, i, ctx.p * total);
        }
        return G;
    }

    private static Matrix computeB(RetrialContext ctx, int nu) {
        Matrix B = new Matrix(ctx.d, ctx.d);

        Matrix[] L = new Matrix[ctx.N + 1];
        Matrix[] A = new Matrix[ctx.N + 1];
        Matrix[] Pmat = new Matrix[ctx.N + 1];
        Matrix[] Delta = new Matrix[ctx.N + 1];

        for (int n = 0; n <= ctx.N; n++) {
            L[n] = computeL(ctx, n);
            A[n] = computeA(ctx, n);
            Pmat[n] = computeP(ctx, n);
            Delta[n] = computeDelta(ctx, n);
        }

        for (int n = 0; n <= ctx.N; n++) {
            int rowStart = getBlockOffset(ctx.T, n);

            Matrix G_nn = computeG_nn(ctx, n, nu, nu);
            if (n == 0) {
                B.set(rowStart, rowStart, G_nn.get(0, 0));
            } else {
                int sz = ctx.T[n];
                for (int i = 0; i < sz; i++) {
                    for (int j = 0; j < sz; j++) {
                        B.set(rowStart + i, rowStart + j, A[n].get(i, j) + Delta[n].get(i, j) + G_nn.get(i, j));
                    }
                }
            }

            if (n >= 1) {
                int colStart = getBlockOffset(ctx.T, n - 1);
                Matrix Ln = L[n];
                for (int i = 0; i < Ln.getNumRows(); i++) {
                    for (int j = 0; j < Ln.getNumCols(); j++) {
                        B.set(rowStart + i, colStart + j, Ln.get(i, j));
                    }
                }
            }

            for (int k = 1; k <= ctx.K; k++) {
                if (n + k <= ctx.N) {
                    int colStart = getBlockOffset(ctx.T, n + k);
                    double D_k_nu_nu = ctx.D[k].get(nu, nu);

                    Matrix Pprod = Matrix.eye(ctx.T[n]);
                    for (int jj = n; jj < n + k; jj++) {
                        if (jj < ctx.N) {
                            Pprod = Pprod.mult(Pmat[jj]);
                        }
                    }
                    Matrix block = Pprod.scale(D_k_nu_nu);
                    for (int i = 0; i < block.getNumRows(); i++) {
                        for (int j = 0; j < block.getNumCols(); j++) {
                            B.set(rowStart + i, colStart + j, block.get(i, j));
                        }
                    }
                }
            }
        }
        return B;
    }

    private static Matrix computeBbar(RetrialContext ctx, int nu) {
        Matrix Bbar = new Matrix(ctx.d, ctx.d);
        int upper = Math.min(ctx.R[nu], ctx.N - 1);
        for (int n = 0; n <= upper; n++) {
            int rowStart = getBlockOffset(ctx.T, n);
            int colStart = getBlockOffset(ctx.T, n + 1);
            Matrix P_n = computeP(ctx, n);
            if (P_n != null) {
                for (int i = 0; i < P_n.getNumRows(); i++) {
                    for (int j = 0; j < P_n.getNumCols(); j++) {
                        Bbar.set(rowStart + i, colStart + j, P_n.get(i, j));
                    }
                }
            }
        }
        return Bbar;
    }

    private static Matrix computeBtilde(RetrialContext ctx, int nu, int nuPrime) {
        Matrix Btilde = new Matrix(ctx.d, ctx.d);

        Matrix[] Pmat = new Matrix[ctx.N + 1];
        for (int n = 0; n <= ctx.N; n++) {
            Pmat[n] = computeP(ctx, n);
        }

        for (int n = 0; n <= ctx.N; n++) {
            int rowStart = getBlockOffset(ctx.T, n);

            Matrix G_nn = computeG_nn(ctx, n, nu, nuPrime);
            for (int i = 0; i < G_nn.getNumRows(); i++) {
                for (int j = 0; j < G_nn.getNumCols(); j++) {
                    Btilde.set(rowStart + i, rowStart + j, G_nn.get(i, j));
                }
            }

            for (int k = 1; k <= ctx.K; k++) {
                if (n + k <= ctx.N) {
                    int colStart = getBlockOffset(ctx.T, n + k);
                    double D_k_nu_nuPrime = ctx.D[k].get(nu, nuPrime);

                    Matrix Pprod = Matrix.eye(ctx.T[n]);
                    for (int jj = n; jj < n + k; jj++) {
                        if (jj < ctx.N) {
                            Pprod = Pprod.mult(Pmat[jj]);
                        }
                    }
                    Matrix block = Pprod.scale(D_k_nu_nuPrime);
                    for (int i = 0; i < block.getNumRows(); i++) {
                        for (int j = 0; j < block.getNumCols(); j++) {
                            Btilde.set(rowStart + i, colStart + j, block.get(i, j));
                        }
                    }
                }
            }
        }
        return Btilde;
    }

    private static Matrix computeC(RetrialContext ctx, int n, int k, int nu, int nuPrime) {
        if (n < ctx.N - ctx.K + k) {
            return new Matrix(ctx.T[n], ctx.T[ctx.N]);
        } else if (n < ctx.N) {
            int batchSize = ctx.N - n + k;
            if (batchSize < 1 || batchSize > ctx.K) {
                return new Matrix(ctx.T[n], ctx.T[ctx.N]);
            }
            double D_batch = ctx.D[batchSize].get(nu, nuPrime);

            Matrix[] Pmat = new Matrix[ctx.N + 1];
            for (int nn = 0; nn <= ctx.N; nn++) {
                Pmat[nn] = computeP(ctx, nn);
            }

            Matrix Pprod = Matrix.eye(ctx.T[n]);
            for (int jj = n; jj < ctx.N; jj++) {
                if (Pmat[jj] != null) {
                    Pprod = Pprod.mult(Pmat[jj]);
                }
            }
            return Pprod.scale((1.0 - ctx.p) * D_batch);
        } else {
            // n == N
            if (k < 1 || k > ctx.K) {
                return new Matrix(ctx.T[ctx.N], ctx.T[ctx.N]);
            }
            double D_k = ctx.D[k].get(nu, nuPrime);
            Matrix C = Matrix.eye(ctx.T[ctx.N]);
            return C.scale((1.0 - ctx.p) * D_k);
        }
    }

    private static Matrix buildGeneratorLevel(RetrialContext ctx, int i, int j) {
        int Vd = ctx.V * ctx.d;
        Matrix Q = new Matrix(Vd, Vd);

        if (j < Math.max(0, i - 1) || j > i + ctx.K) {
            return Q;
        }

        Matrix[] B = new Matrix[ctx.V];
        Matrix[] Bbar = new Matrix[ctx.V];
        Matrix[] Gam = new Matrix[ctx.V];
        for (int nu = 0; nu < ctx.V; nu++) {
            B[nu] = computeB(ctx, nu);
            Bbar[nu] = computeBbar(ctx, nu);
            Gam[nu] = computeGamma(ctx, nu);
        }

        if (i == j) {
            for (int nu = 0; nu < ctx.V; nu++) {
                int rowStart = nu * ctx.d;
                for (int nuPrime = 0; nuPrime < ctx.V; nuPrime++) {
                    int colStart = nuPrime * ctx.d;
                    if (nu == nuPrime) {
                        double D0_nu_nu = ctx.D[0].get(nu, nu);
                        for (int r = 0; r < ctx.d; r++) {
                            for (int c = 0; c < ctx.d; c++) {
                                double value = B[nu].get(r, c);
                                if (r == c) {
                                    value += D0_nu_nu - ((double) i) * (ctx.gamma + ctx.alpha) + ((double) i) * ctx.alpha * Gam[nu].get(r, c);
                                } else {
                                    value += ((double) i) * ctx.alpha * Gam[nu].get(r, c);
                                }
                                Q.set(rowStart + r, colStart + c, value);
                            }
                        }
                    } else {
                        Matrix Btilde = computeBtilde(ctx, nu, nuPrime);
                        double D0_nu_nuPrime = ctx.D[0].get(nu, nuPrime);
                        for (int r = 0; r < ctx.d; r++) {
                            for (int c = 0; c < ctx.d; c++) {
                                double value = Btilde.get(r, c);
                                if (r == c) {
                                    value += D0_nu_nuPrime;
                                }
                                Q.set(rowStart + r, colStart + c, value);
                            }
                        }
                    }
                }
            }
        } else if (j == i - 1 && i >= 1) {
            for (int nu = 0; nu < ctx.V; nu++) {
                int rowStart = nu * ctx.d;
                int colStart = rowStart;
                for (int r = 0; r < ctx.d; r++) {
                    for (int c = 0; c < ctx.d; c++) {
                        double value = ((double) i) * ctx.alpha * Bbar[nu].get(r, c);
                        if (r == c) {
                            value += ((double) i) * ctx.gamma;
                        }
                        Q.set(rowStart + r, colStart + c, value);
                    }
                }
            }
        } else if (j > i && j <= i + ctx.K) {
            int kk = j - i;
            for (int nu = 0; nu < ctx.V; nu++) {
                int rowStart = nu * ctx.d;
                for (int nuPrime = 0; nuPrime < ctx.V; nuPrime++) {
                    int colStart = nuPrime * ctx.d;
                    Matrix block = new Matrix(ctx.d, ctx.d);
                    for (int n = 0; n <= ctx.N; n++) {
                        Matrix C_nk = computeC(ctx, n, kk, nu, nuPrime);
                        if (C_nk != null && C_nk.getNumCols() == ctx.T[ctx.N]) {
                            int nRowStart = getBlockOffset(ctx.T, n);
                            int NColStart = getBlockOffset(ctx.T, ctx.N);
                            for (int r = 0; r < ctx.T[n]; r++) {
                                for (int c = 0; c < ctx.T[ctx.N]; c++) {
                                    double v = C_nk.get(r, c);
                                    if (v != 0.0) {
                                        block.set(nRowStart + r, NColStart + c, block.get(nRowStart + r, NColStart + c) + v);
                                    }
                                }
                            }
                        }
                    }
                    for (int r = 0; r < ctx.d; r++) {
                        for (int c = 0; c < ctx.d; c++) {
                            Q.set(rowStart + r, colStart + c, block.get(r, c));
                        }
                    }
                }
            }
        }

        return Q;
    }
}
