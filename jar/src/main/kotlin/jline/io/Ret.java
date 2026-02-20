/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.lang.Network;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

import java.util.*;
import java.util.function.Function;


/**
 * Container class for return types used throughout the LINE queueing network solver library.
 * 
 * <p>This class provides a centralized collection of data structures used as return types
 * from various queueing network analysis algorithms. Each inner class represents a specific
 * return type tailored to the output requirements of different solvers and methods.</p>
 * 
 * <p>The return types are organized into several categories:
 * <ul>
 *   <li><b>MVA (Mean Value Analysis) variants:</b> pfqnMVA, pfqnMVALD, pfqnMVALDMX, etc.</li>
 *   <li><b>Normalizing constant methods:</b> pfqnNc, pfqnNcXQ, pfqnNcComplex, etc.</li>
 *   <li><b>Approximation methods:</b> pfqnAMVA, LinearizerResult, etc.</li>
 *   <li><b>Cache analysis:</b> cacheMVA, cacheXiFp, cacheGamma, etc.</li>
 *   <li><b>Distribution fitting:</b> mamAPH2Fit, mamMMAPMixtureFit, etc.</li>
 *   <li><b>Sampling and simulation:</b> ctmcSimulation, SampleResult, etc.</li>
 *   <li><b>Solver results:</b> DistributionResult, ProbabilityResult, etc.</li>
 * </ul>
 * </p>
 * 
 * <p>These return types enable type-safe communication between different components
 * of the solver library and provide clear interfaces for algorithm outputs.</p>
 * 
 * @since 1.0
 */
public class Ret {

    /**
     * Result type for event-based state space exploration functions.
     * 
     * <p>This class encapsulates the output of afterEvent* functions which compute
     * the possible states reachable after an event occurs in the queueing network,
     * along with the rates and probabilities of transitions.</p>
     */
    public static class EventResult {

        public final Matrix outspace;
        public final Matrix outrate;
        public final Matrix outprob;


        public EventResult(Matrix outspace, Matrix outrate, Matrix outprob) {
            this.outspace = outspace;
            this.outrate = outrate;
            this.outprob = outprob;
        }

    }

    /**
     * Index key for caching intermediate results in the pfqn_gld (Generalized Local Balance) algorithm.
     * 
     * <p>This class provides a three-dimensional index (a, b, c) used as a key
     * for storing and retrieving intermediate values computed during the
     * generalized local balance algorithm for product-form queueing networks.</p>
     */
    public static class pfqnGldIndex {
        int a;
        int b;
        int c;

        public pfqnGldIndex(int a, int b, int c) {
            this.a = a;
            this.b = b;
            this.c = c;
        }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof pfqnGldIndex)) {
                return false;
            }
            pfqnGldIndex other = (pfqnGldIndex) obj;
            return (this.a == other.a) && (this.b == other.b)
                    && (this.c == other.c);
        }

        @Override
        public int hashCode() {
            return (a + b + c) * a * b * c;
        }
    }

    /**
     * Result type for the MVA (Mean Value Analysis) algorithm.
     * 
     * <p>This class encapsulates the key performance metrics computed by the MVA algorithm
     * for product-form queueing networks. MVA is an iterative algorithm that computes
     * steady-state performance measures without explicitly solving the underlying Markov chain.</p>
     * 
 * @see jline.api.PFQN#pfqn_mva
     */
    public static class pfqnMVA {
        public Matrix X;
        public Matrix Q;
        public Matrix U;
        public Matrix R;
        public double lGN;

        public pfqnMVA(Matrix XN, Matrix QN, Matrix UN, Matrix RN, double lGN) {
            this.X = XN;
            this.Q = QN;
            this.U = UN;
            this.R = RN;
            this.lGN = lGN;
        }
    }

    /**
     * Extended result type for MVA with load-dependent stations and additional metrics.
     * 
     * <p>This class extends the basic MVA results with additional information needed
     * for load-dependent queueing stations, including state-dependent service rates
     * and probability distributions over the number of customers at each station.</p>
     * 
     * @see jline.api.PFQN#pfqn_mvaldmx
     */
    public static class pfqnMVALDMX {
        public Matrix X;
        public Matrix Q;
        public Matrix U;
        public Matrix R;
        public double lG;
        public Matrix Pc;

        public pfqnMVALDMX(Matrix XN, Matrix QN, Matrix UN, Matrix RN, double lG, Matrix Pc) {
            this.X = XN;
            this.Q = QN;
            this.U = UN;
            this.R = RN;
            this.lG = lG;
            this.Pc = Pc;
        }
    }

    /**
     * Data structure to hold extended results from the MVA computation, particularly
     * focusing on error corrections. It includes the matrices EC, E, Eprime, and Lo.
     */
    public static class pfqnMVALDMXEC {
        public Matrix EC;
        public Matrix E;
        public Matrix Eprime;
        public Matrix Lo;

        public pfqnMVALDMXEC(Matrix EC, Matrix E, Matrix Eprime, Matrix Lo) {
            this.EC = EC;
            this.E = E;
            this.Eprime = Eprime;
            this.Lo = Lo;
        }
    }

    /**
     * Data structure for storing results from an MVA computation using the LD method.
     * Contains performance metrics such as throughput (X), queue length (Q), utilization (U),
     * response time (R), a list of logarithmic normalization constants (lG), a numerical
     * stability flag (isNumStable), and a probability matrix (pi).
     */
    public static class pfqnMVALD {
        public Matrix X;
        public Matrix Q;
        public Matrix U;
        public Matrix R;
        public List<Double> lG;
        public boolean isNumStable;
        public Matrix pi;

        public pfqnMVALD(Matrix XN, Matrix QN, Matrix UN, Matrix RN, List<Double> lG, boolean isNumStable, Matrix pi) {
            this.X = XN;
            this.Q = QN;
            this.U = UN;
            this.R = RN;
            this.lG = lG;
            this.isNumStable = isNumStable;
            this.pi = pi;
        }
    }

    /**
     * Data structure for storing results from the AMVA MS (Approximate Mean Value Analysis
     * Multiservice) method. Includes queue length (Q), utilization (U), response time (R),
     * a matrix of class-think times (C), throughput (X), and the total number of iterations (totiter).
     */
    public static class pfqnAMVAMS {
        public Matrix Q;
        public Matrix U;
        public Matrix R;
        public Matrix C;
        public Matrix X;
        public int totiter;

        public pfqnAMVAMS(Matrix Q, Matrix U, Matrix R, Matrix C, Matrix X, int totiter) {
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.C = C;
            this.X = X;
            this.totiter = totiter;
        }
    }

    /**
     * Result type for the Schmidt AMVA algorithm.
     *
     * <p>This class encapsulates the output of pfqn_schmidt and pfqn_schmidt_ext functions
     * which implement the Schmidt population recursion algorithm for analyzing closed
     * queueing networks with multi-class workloads and class-dependent FCFS scheduling.</p>
     */
    public static class pfqnAMVASchmidt {
        /** System throughput vector (1 x R) */
        public Matrix X;
        /** Queue length matrix (M x R) */
        public Matrix Q;
        /** Utilization matrix (M x R) */
        public Matrix U;
        /** Waiting time matrix (M x R) */
        public Matrix C;
        /** State probabilities per station */
        public List<Map<Object, Double>> P;

        public pfqnAMVASchmidt(Matrix XN, Matrix QN, Matrix UN, Matrix CN, List<Map<Object, Double>> PN) {
            this.X = XN;
            this.Q = QN;
            this.U = UN;
            this.C = CN;
            this.P = PN;
        }
    }

    /**
     * Unified result type for linearizer approximation methods.
     * 
     * <p>The linearizer is an iterative approximation technique for solving queueing
     * networks with non-product-form features. This class provides a flexible structure
     * that accommodates results from different linearizer variants including:
     * <ul>
     *   <li>Basic linearizer with queue lengths and waiting times</li>
     *   <li>Multi-server linearizer with blocking probabilities</li>
     *   <li>Forward MVA linearizer without iteration counts</li>
     * </ul>
     * </p>
     * 
     * <p>Fields may be null depending on the specific linearizer variant used.</p>
     */
    public static class LinearizerResult {
        public Matrix Q;   // Queue length
        public Matrix W;   // Wait time
        public Matrix T;   // Throughput (T in some algorithms)
        public Matrix X;   // Throughput (X in other algorithms) 
        public Matrix P;   // Probability matrices (optional)
        public Matrix PB;  // Blocking probabilities (optional)
        public Integer iter; // Number of iterations (optional)

        // Full constructor for all fields
        public LinearizerResult(Matrix Q, Matrix W, Matrix T, Matrix X, Matrix P, Matrix PB, Integer iter) {
            this.Q = Q;
            this.W = W;
            this.T = T;
            this.X = X;
            this.P = P;
            this.PB = PB;
            this.iter = iter;
        }

        // Constructor for basic linearizer (Q, W, T, iter)
        public LinearizerResult(Matrix Q, Matrix W, Matrix T, int iter) {
            this(Q, W, T, null, null, null, iter);
        }

        // Constructor for forward MVA (Q, W, T) without iterations
        public LinearizerResult(Matrix Q, Matrix W, Matrix T) {
            this(Q, W, T, null, null, null, null);
        }

        // Constructor for multi-server with probabilities (Q, W, T, P, PB, iter)
        public LinearizerResult(Matrix Q, Matrix W, Matrix T, Matrix P, Matrix PB, int iter) {
            this(Q, W, T, null, P, PB, iter);
        }

        // Constructor for multi-server forward MVA (Q, W, T, P, PB) without iterations
        public LinearizerResult(Matrix Q, Matrix W, Matrix T, Matrix P, Matrix PB) {
            this(Q, W, T, null, P, PB, null);
        }

        // Constructor for core with X instead of T (Q, W, X, P, PB, iter)
        public static LinearizerResult withX(Matrix Q, Matrix W, Matrix X, Matrix P, Matrix PB, int iter) {
            return new LinearizerResult(Q, W, null, X, P, PB, iter);
        }
    }


    /**
     * Data structure for storing estimated intermediate results from the MS linearizer method.
     * Includes arrays of matrices for queue length estimates (Q_1), probabilities (P_1), and a
     * matrix of blocking probabilities (PB_1).
     */
    public static class pfqnLinearizerMSEstimate {
        public Matrix[] Q_1;
        public Matrix[] P_1;
        public Matrix PB_1;

        public pfqnLinearizerMSEstimate(Matrix[] Q_1, Matrix[] P_1, Matrix PB_1) {
            this.Q_1 = Q_1;
            this.P_1 = P_1;
            this.PB_1 = PB_1;
        }
    }


    /**
     * Data structure for storing results from the AMVA (Approximate Mean Value Analysis) method.
     * Contains queue length (Q), utilization (U), response time (R), throughput (T),
     * a matrix for class-think times (C), throughput (X), and the total number of iterations (totiter).
     */
    public static class pfqnAMVA {
        public Matrix Q;
        public Matrix U;
        public Matrix R;
        public Matrix T;
        public Matrix C;
        public Matrix X;
        public int totiter;

        public pfqnAMVA(Matrix Q, Matrix U, Matrix R, Matrix T, Matrix C, Matrix X, int totiter) {
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.T = T;
            this.C = C;
            this.X = X;
            this.totiter = totiter;
        }
    }


    /**
     * Data structure for storing intermediate estimates from the linearizer method.
     * Contains arrays of matrices for queue length estimates (Q_1) and a matrix for
     * the throughput estimates (T_1).
     */
    public static class pfqnLinearizerEstimate {
        public Matrix[] Q_1;
        public Matrix T_1;

        public pfqnLinearizerEstimate(Matrix[] Q_1, Matrix T_1) {
            this.Q_1 = Q_1;
            this.T_1 = T_1;
        }
    }


    /**
     * Data structure for storing results from a fixed-point iteration method.
     * Contains a matrix of utilization (u) and a matrix of differences (d).
     */
    static public class pfqnLeFpi {
        public final Matrix u;
        public final Matrix d;

        public pfqnLeFpi(Matrix u, Matrix d) {
            this.u = u;
            this.d = d;
        }
    }

    /**
     * Data structure for storing results from a fixed-point iteration method with normalization.
     * Contains a matrix of utilization (u), a normalization constant (v), and a matrix of differences (d).
     */
    static public class pfqnLeFpiZ {
        public final Matrix u;
        public final double v;
        public final Matrix d;

        public pfqnLeFpiZ(Matrix u, double v, Matrix d) {
            this.u = u;
            this.v = v;
            this.d = d;
        }
    }

    /**
     * Data structure for storing results the COMOM method.
     * Contains a logarithmic normalization constant (lG) and a basis matrix (lGbasis).
     */
    public static class pfqnComomrm {
        public double lG;
        public Matrix lGbasis;

        public pfqnComomrm(double lG, Matrix lGbasis) {
            this.lG = lG;
            this.lGbasis = lGbasis;
        }
    }

    /**
     * Data structure for storing sanitized input parameters for a normalizing constant calculation.
     * Includes arrival rates (lambda), service demand matrix (L), population vector (N),
     * think times (Z), and a logarithmic remainder (lGremaind).
     */
    public static class pfqnNcSanitize {
        public Matrix lambda;
        public Matrix L;
        public Matrix N;
        public Matrix Z;
        public double lGremaind;

        public pfqnNcSanitize(Matrix lambda, Matrix L, Matrix N, Matrix Z, double lGremaind) {
            this.lambda = lambda;
            this.L = L;
            this.N = N;
            this.Z = Z;
            this.lGremaind = lGremaind;
        }
    }

    /**
     * Data structure for storing results from a normalizing constant calculation involving throughputs and
     * queue lengths.
     * Contains normalization constants (G, lG), throughput (X), queue length (Q), and method description.
     */
    public static class pfqnNcXQ {
        public Double G;
        public Double lG;
        public Matrix X;
        public Matrix Q;
        public String method;

        public pfqnNcXQ(Double G, Double lG, Matrix X, Matrix Q, String method) {
            this.G = G;
            this.lG = lG;
            this.X = new Matrix(X);
            this.Q = new Matrix(Q);
            this.method = method;
        }

        public pfqnNcXQ(Double lG, Matrix X, Matrix Q, String method) {
            this.G = FastMath.exp(lG);
            this.lG = lG;
            this.X = new Matrix(X);
            this.Q = new Matrix(Q);
            this.method = method;
        }
    }

    /**
     * Result type for normalizing constant computations in product-form queueing networks.
     * 
     * <p>The normalizing constant G(N) is a fundamental quantity in product-form
     * queueing network analysis. It appears in the denominator of the steady-state
     * probability distribution and is used to compute performance measures.</p>
     * 
     * <p>This class stores both the normalizing constant G and its natural logarithm lG
     * to handle numerical issues with very large or small values.</p>
     */
    public static class pfqnNc {
        public Double G;
        public Double lG;
        public String method;

        public pfqnNc(Double G, Double lG) {
            this.G = G;
            this.lG = lG;
            this.method = null;
        }

        public pfqnNc(Double G, Double lG, String method) {
            this.lG = lG;
            this.G = G;
            this.method = method;
        }
    }

    /**
     * Data structure for storing results from the RD method.
     * Contains a logarithmic normalization constant (lG) and an optional coefficient (Cgamma).
     */
    public static class pfqnRd {
        public Double lG;
        public Double Cgamma;

        public pfqnRd(Double lG, Double Cgamma) {
            this.lG = lG;
            this.Cgamma = Cgamma;
        }

        public pfqnRd(Double lG) {
            this.lG = lG;
            this.Cgamma = null;
        }

    }

    /**
     * Data structure for storing complex results from a normalizing constant calculation.
     * Contains complex normalization constants (G, lG) and an optional method description.
     */
    public static class pfqnNcComplex {
        public Complex G;
        public Complex lG;
        public String method;

        public pfqnNcComplex(Complex G, Complex lG) {
            this.G = G;
            this.lG = lG;
            this.method = null;
        }
    }

    /**
     * Data structure for storing results from a FNC (Fitting Normalizing Constants) calculation.
     * Contains matrices of mean service rates (mu) and coefficient of variation (c).
     */
    public static class pfqnFnc {
        public Matrix mu;
        public Matrix c;

        public pfqnFnc(Matrix mu, Matrix c) {
            this.mu = mu;
            this.c = c;
        }
    }

    /**
     * Function class for the integrand used in cubature calculations.
     * This function calculates the integrand value given a vector of inputs,
     * based on the service demand matrix (L), total population (Nt), and a
     * vector of coefficients (beta).
     */
    public static class pfqnCUB implements Function<double[], Double> {
        double[][] L;
        double Nt;
        double[] beta;

        public pfqnCUB(Matrix L, double Nt, Matrix beta) {
            this.L = L.toArray2D();
            this.Nt = Nt;
            this.beta = beta.toArray1D();
        }

        @Override
        public Double apply(double[] x) {
            double sumX = Arrays.stream(x).sum();
            double oneMinusSumX = 1.0 - sumX;

            double[] extendedX = Arrays.copyOf(x, x.length + 1);
            extendedX[x.length] = oneMinusSumX;

            // Transpose extendedX and multiply with Lv
            double[] result = new double[L[0].length];
            Arrays.fill(result, 0.0);
            for (int j = 0; j < L[0].length; j++) {
                for (int i = 0; i < extendedX.length; i++) {
                    result[j] += extendedX[i] * L[i][j];
                }
            }
            for (int i = 0; i < result.length; i++) {
                result[i] = beta[i] * FastMath.log(result[i]);
            }
            double sumhx = Arrays.stream(result).sum();

            return FastMath.exp(Nt * sumhx);
        }
    }

    /**
     * Data structure for storing results from the load-dependent COMOM method.
     * Contains a normalization constant (GN), a logarithmic normalization constant (lG),
     * and a matrix of probabilities (prob).
     */
    public static class pfqnComomrmLd {
        public double GN;
        public double lG;
        public Matrix prob;

        public pfqnComomrmLd(double GN, double lG, Matrix prob) {
            this.GN = GN;
            this.lG = lG;
            this.prob = prob;
        }
    }

    /**
     * Data structure for storing results from the procomom2 method.
     * Contains marginal state probabilities (pk), logarithmic normalization constant (lG),
     * normalization constant (G), transition matrices (T), and auxiliary matrices (F, B).
     */
    public static class pfqnProcomom2 {
        public Matrix pk;
        public double lG;
        public double G;
        public List<Matrix> T;
        public Matrix F;
        public Matrix B;

        public pfqnProcomom2(Matrix pk, double lG, double G, List<Matrix> T, Matrix F, Matrix B) {
            this.pk = pk;
            this.lG = lG;
            this.G = G;
            this.T = T;
            this.F = F;
            this.B = B;
        }
    }


    /**
     * Data structure for storing linearizer estimtate results from a queueing network analysis.
     * Contains queue length (Q), wait time (W), throughput (X), probability matrices (P, PB),
     * and the number of iterations performed (iter).
     */
    public static class pfqnEstimate {
        public final MatrixCell Q_1;
        public final MatrixCell P_1;
        public final Matrix PB_1;
        public final Matrix T_1;

        public pfqnEstimate(MatrixCell Q_1, MatrixCell P_1, Matrix PB_1, Matrix T_1) {
            this.Q_1 = Q_1;
            this.P_1 = P_1;
            this.PB_1 = PB_1;
            this.T_1 = T_1;
        }
    }


    /**
     * Represents the return type for cache miss rate computations with the SPM method.
     * This class encapsulates various metrics related to the cache miss rates and probabilities.
     */
    public static class cacheMissSpm {
        double M;
        double[] MU;
        double[] MI;
        double[] pi0;
        double lE;

        public cacheMissSpm(double M, double[] MU, double[] MI, double[] pi0, double lE) {
            this.M = M;
            this.MU = MU;
            this.MI = MI;
            this.pi0 = pi0;
            this.lE = lE;
        }

        public double getM() {
            return M;
        }

        public double[] getMI() {
            return MI;
        }

        public double[] getMU() {
            return MU;
        }

        public double[] getPi0() {
            return pi0;
        }

        public double getlE() {
            return lE;
        }
    }
    
    /**
     * @deprecated Use cacheMissSpm instead
     */
    @Deprecated
    public static class cacheMissRayInt extends cacheMissSpm {
        public cacheMissRayInt(double M, double[] MU, double[] MI, double[] pi0, double lE) {
            super(M, MU, MI, pi0, lE);
        }
    }

    /**
     * Result type for the cache characteristic time (xi) fixed-point algorithm.
     * 
     * <p>This algorithm computes the characteristic times (xi) for cache replacement
     * policies using a fixed-point iteration. The characteristic time represents the
     * average time between consecutive misses for each item in the cache.</p>
     * 
     * <p>The results include:
     * <ul>
     *   <li>xi: Characteristic times for each item</li>
     *   <li>pi0: Steady-state probability that each item is not in cache</li>
     *   <li>pij: Joint probabilities for cache states</li>
     *   <li>it: Number of iterations until convergence</li>
     * </ul>
     * </p>
     */
    public static class cacheXiFp {
        public Matrix xi;    // Matrix representing the xi terms
        public Matrix pi0;   // Matrix representing the initial state probabilities
        public Matrix pij;   // Matrix representing the cache state probabilities
        public int it;       // Integer representing the number of iterations

        /**
         * Creates a new cacheXiFp result object.
         *
         * @param xi  Matrix of characteristic times (one per item)
         * @param pi0 Matrix of miss probabilities (probability each item is not in cache)
         * @param pij Matrix of joint cache state probabilities
         * @param it  Number of fixed-point iterations performed
         */
        public cacheXiFp(Matrix xi, Matrix pi0, Matrix pij, int it) {
            this.xi = xi;
            this.pi0 = pi0;
            this.pij = pij;
            this.it = it;
        }
    }

    /**
     * Represents the return type for the cache gamma linear program computations.
     * This class encapsulates the computed gamma matrix along with the number of users, items, and levels.
     */
    public static class cacheGamma {
        public Matrix gamma; // Matrix representing the gamma values
        public int u;        // Integer representing the number of users
        public int n;        // Integer representing the number of items
        public int h;        // Integer representing the number of levels

        /**
         * Constructor for initializing the cacheGammaLpReturn object.
         *
         * @param gamma - Matrix representing the gamma values.
         * @param u     - Integer representing the number of users.
         * @param n     - Integer representing the number of items.
         * @param h     - Integer representing the number of levels.
         */
        public cacheGamma(Matrix gamma, int u, int n, int h) {
            this.gamma = gamma;
            this.u = u;
            this.n = n;
            this.h = h;
        }
    }

    /**
     * Represents the return type for the cache MVA (Mean Value Analysis) computations.
     * This class encapsulates the results including various probability matrices and other related metrics.
     */
    public static class cacheMVA {
        public Matrix pi;    // Matrix representing the cache state probabilities
        public Matrix pi0;   // Matrix representing the initial state probabilities
        public Matrix pij;   // Matrix representing the cache state probabilities for specific items and levels
        public Matrix x;     // Matrix representing additional metric data
        public Matrix u;     // Matrix representing another computed metric
        public int E;        // Integer representing an additional metric or counter

        /**
         * Constructor for initializing the cacheMVAReturn object.
         *
         * @param pi  - Matrix representing the cache state probabilities.
         * @param pi0 - Matrix representing the initial state probabilities.
         * @param pij - Matrix representing the cache state probabilities for specific items and levels.
         * @param x   - Matrix representing additional metric data.
         * @param u   - Matrix representing another computed metric.
         * @param E   - Integer representing an additional metric or counter.
         */
        public cacheMVA(Matrix pi, Matrix pi0, Matrix pij, Matrix x, Matrix u, int E) {
            this.pi = pi;
            this.pi0 = pi0;
            this.pij = pij;
            this.x = x;
            this.u = u;
            this.E = E;
        }
    }

    /**
     * Represents the return type for cache ray method.
     * This class encapsulates the results including the partition function (Z),
     * its logarithm (lZ), and the xi terms.
     */
    public static class cacheSpm {
        public double Z;   // Double representing the partition function value
        public double lZ;  // Double representing the logarithm of the partition function
        public Matrix xi;  // Matrix representing the xi terms

        /**
         * Constructor for initializing the cacheSpmReturn object.
         *
         * @param Z  - Double representing the partition function value.
         * @param lZ - Double representing the logarithm of the partition function.
         * @param xi - Matrix representing the xi terms.
         */
        public cacheSpm(double Z, double lZ, Matrix xi) {
            this.Z = Z;
            this.lZ = lZ;
            this.xi = xi;
        }

    }
    
    /**
     * @deprecated Use cacheSpm instead
     */
    @Deprecated
    public static class cacheRayInt extends cacheSpm {
        public cacheRayInt(double Z, double lZ, Matrix xi) {
            super(Z, lZ, xi);
        }
    }

    /**
     * Represents the return type for cache importance sampling method.
     * This class encapsulates the estimated normalizing constant E and its logarithm lE.
     */
    public static class cacheIs {
        public double E;   // Estimated normalizing constant
        public double lE;  // Logarithm of the normalizing constant

        /**
         * Constructor for initializing the cacheIs object.
         *
         * @param E  - Estimated normalizing constant.
         * @param lE - Logarithm of the normalizing constant.
         */
        public cacheIs(double E, double lE) {
            this.E = E;
            this.lE = lE;
        }
    }

    public static class FJAuxClassKey {

        private final int x;
        private final int y;

        public FJAuxClassKey(int x, int y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof FJAuxClassKey)) return false;
            FJAuxClassKey key = (FJAuxClassKey) o;
            return x == key.x && y == key.y;
        }

        @Override
        public int hashCode() {
            return Objects.hash(x, y);
        }
    }

    public static class FJsortForks {
        public Matrix outerForks;
        public Matrix parentForks;

        public FJsortForks(Matrix outerForks, Matrix parentForks) {
            this.outerForks = outerForks;
            this.parentForks = parentForks;
        }
    }

    public static class FJApprox {
        public Network nonfjmodel;
        public Matrix fjclassmap;
        public Matrix fjforkmap;
        public Map<Integer, Integer> fj_auxiliary_delays;
        public Map<Integer, Integer> fanout;

        public FJApprox(Network nonfjmodel, Matrix fjclassmap, Matrix fjforkmap, Map<Integer, Integer> fj_auxiliary_delays, Map<Integer, Integer> fanout) {
            this.nonfjmodel = nonfjmodel;
            this.fjclassmap = fjclassmap;
            this.fjforkmap = fjforkmap;
            this.fj_auxiliary_delays = fj_auxiliary_delays;
            this.fanout = fanout;
        }
    }

    /**
     * Constructs a lossnErlangFPReturn object with the specified queue-length, loss probability,
     * blocking probability, and iteration count.
     */
    public static class lossnErlangFP {
        public Matrix qLen;
        public Matrix lossProb;
        public Matrix blockProb;
        public int niter;

        /*
         * @param q The mean queue-length for each route.
         * @param l The loss probability for each route.
         * @param b The blocking probability for each link.
         * @param n The number of iterations taken to converge.
         */

        public lossnErlangFP(Matrix q, Matrix l, Matrix b, int n) {
            qLen = q;
            lossProb = l;
            blockProb = b;
            niter = n;
        }
    }

    /**
     * Class representing the return type for the fitting of a 2-phase APH (Acyclic Phase-Type) distribution.
     */
    public static class mamAPH2Fit {
        public MatrixCell APH;
        public Map<Integer, MatrixCell> APHS;

        public mamAPH2Fit() {
            APH = new MatrixCell();
            APHS = new HashMap<>();
        }
    }

    /**
     * Constructor initializing the sample data, number of types, and type indices.
     */
    public static class mamMMAPSample {
        double[] samples; // inter-arrival times
        int ntypes; // number of event types
        int[] types; // event type for each sample
        int[] states; // state sequence (state at the time of each event)

        /*
         * @param s  the array of sampled inter-arrival times
         * @param nt the number of different types of arrivals
         * @param t  the array indicating the type of each sampled arrival
         */
        public mamMMAPSample(double[] s, int nt, int[] t) {
            samples = s;
            types = t;
            ntypes = nt;
            states = null;
        }

        /*
         * @param s  the array of sampled inter-arrival times
         * @param nt the number of different types of arrivals
         * @param t  the array indicating the type of each sampled arrival
         * @param sts the array of state indices at each event
         */
        public mamMMAPSample(double[] s, int nt, int[] t, int[] sts) {
            samples = s;
            types = t;
            ntypes = nt;
            states = sts;
        }

        public double[] getSamples() {
            return samples;
        }

        public int[] getTypes() {
            return types;
        }

        public int getNtypes() {
            return ntypes;
        }

        public int[] getStates() {
            return states;
        }
    }

    /**
     * Class representing the return type for fitting a mixture model to a MMAP.
     */
    public static class mamMMAPMixtureFit {
        public MatrixCell MMAP;
        public Map<Integer[], MatrixCell> PHs;

        public mamMMAPMixtureFit() {
            MMAP = new MatrixCell();
            PHs = new HashMap<>();
        }
    }

    public static class ctmcSimulation {
        public int[] states;
        public double[] sojournTimes;
    }

    /**
     * Data structure to return the results of non-product-form queueing network approximation.
     */
    public static class npfqnNonexpApprox {
        public Matrix ST;
        public Matrix gamma;
        public Matrix nservers;
        public Matrix rho;
        public Matrix scva;
        public Matrix scvs;
        public Matrix eta;

        public npfqnNonexpApprox(Matrix ST, Matrix gamma, Matrix nservers, Matrix rho, Matrix scva, Matrix scvs, Matrix eta) {
            this.ST = ST;
            this.gamma = gamma;
            this.nservers = nservers;
            this.rho = rho;
            this.scva = scva;
            this.scvs = scvs;
            this.eta = eta;
        }
    }

    /**
     * A class to store the results of queueing system analysis.
     */
    public static class qsys {
        public static double W;
        public static double rho;

        /**
         * Constructs a qsysReturn object.
         *
         * @param W   Average waiting time.
         * @param rho Utilization.
         */
        public qsys(double W, double rho) {
            qsys.W = W;
            qsys.rho = rho;
        }
    }

    /**
     * Return type for multi-class queueing system analysis with priorities.
     * Contains per-class waiting times and overall system utilization.
     */
    public static class qsys_prio {
        public static Matrix W;   // Vector of waiting times per class
        public static double rho; // Overall utilization

        /**
         * Constructs a qsys_prio return object.
         *
         * @param W   Vector of average waiting times per priority class.
         * @param rho Overall system utilization.
         */
        public qsys_prio(Matrix W, double rho) {
            qsys_prio.W = W;
            qsys_prio.rho = rho;
        }
    }

    /**
     * A unified return type for demand-related methods, supporting both simple demands (D, Z)
     * and comprehensive chain demands with optional chain-specific parameters.
     */
    public static class snGetDemands {
        // Basic demand parameters (always present)
        public Matrix D;
        public Matrix Z;
        
        // Chain-specific parameters (optional, only present for chain operations)
        public Matrix Dchain;
        public Matrix STchain;
        public Matrix Vchain;
        public Matrix alpha;
        public Matrix Nchain;
        public Matrix SCVchain;
        public Matrix refstatchain;
        
        // Constructor for simple demands (D, Z only)
        public snGetDemands(Matrix D, Matrix Z) {
            this.D = D;
            this.Z = Z;
        }
        
        // Constructor for chain demands (all parameters)
        public snGetDemands(Matrix D, Matrix Z, Matrix Dchain, Matrix STchain, Matrix Vchain, 
                           Matrix alpha, Matrix Nchain, Matrix SCVchain, Matrix refstatchain) {
            this.D = D;
            this.Z = Z;
            this.Dchain = Dchain;
            this.STchain = STchain;
            this.Vchain = Vchain;
            this.alpha = alpha;
            this.Nchain = Nchain;
            this.SCVchain = SCVchain;
            this.refstatchain = refstatchain;
        }
        
        // Constructor for chain demands only (legacy snGetDemandsChain compatibility)
        public snGetDemands(Matrix Dchain, Matrix STchain, Matrix Vchain, Matrix alpha,
                           Matrix Nchain, Matrix SCVchain, Matrix refstatchain) {
            this.Dchain = Dchain;
            this.STchain = STchain;
            this.Vchain = Vchain;
            this.alpha = alpha;
            this.Nchain = Nchain;
            this.SCVchain = SCVchain;
            this.refstatchain = refstatchain;
        }
    }

    /**
     * A unified return type for methods returning product form parameters.
     * Supports both class-level and chain-level product form parameters.
     */
    public static class snGetProductFormParams {
        public Matrix lambda;
        public Matrix D;
        public Matrix N;
        public Matrix Z;
        public Matrix mu;
        public Matrix S;
        public Matrix V;
    }


    /**
     * A return type for the snDeaggregateChainResults method, encapsulating multiple chain-related matrix results.
     */
    public static class snDeaggregateChainResults {
        public Matrix Q;
        public Matrix U;
        public Matrix R;
        public Matrix T;
        public Matrix C;
        public Matrix X;

        public snDeaggregateChainResults(Matrix q, Matrix u, Matrix r, Matrix t,
                                         Matrix c, Matrix x) {
            Q = q;
            U = u;
            R = r;
            T = t;
            C = c;
            X = x;
        }
    }


    /**
     * A class to represent the return type of the map2_fit function, holding the transition matrices and possibly other fitting results.
     */
    public static class mamMAPFitReturn {
        public MatrixCell MAP;
        public double error;

        public mamMAPFitReturn() {
            MAP = new MatrixCell();
            error = 0;
        }
    }

    public static class Eigs {
        public Matrix values; // vector of eigenvalues
        public Matrix vectors; // eigenvectors

        public Eigs(Matrix s, Matrix V) {
            this.values = s;
            this.vectors = V;
        }
    }

    /**
     * Result class for SVD (Singular Value Decomposition)
     */
    public static class SVD {
        public Matrix u;  // left singular vectors (m x m)
        public Matrix s;  // singular values as column vector (min(m,n) x 1)
        public Matrix v;  // right singular vectors (n x n)

        public SVD(Matrix u, Matrix s, Matrix v) {
            this.u = u;
            this.s = s;
            this.v = v;
        }
    }

    public static class SpectralDecomposition {
        public Matrix spectrum;
        public MatrixCell projectors;

        public SpectralDecomposition(Matrix s, MatrixCell p) {
            this.spectrum = s;
            this.projectors = p;
        }
    }

    /**
     * Result class for getHashOrAdd method
     */
    public static class getHashOrAddResult {
        public Matrix hashid;
        public jline.lang.NetworkStruct sn;

        public getHashOrAddResult(Matrix hashid, jline.lang.NetworkStruct sn) {
            this.hashid = hashid;
            this.sn = sn;
        }
    }

    /**
     * Result class for afterEventHashedOrAdd method
     */
    public static class afterEventHashedOrAddResult {
        public Matrix outhash;
        public Matrix outrate;
        public Matrix outprob;
        public jline.lang.NetworkStruct sn;

        public afterEventHashedOrAddResult(Matrix outhash, Matrix outrate, Matrix outprob, jline.lang.NetworkStruct sn) {
            this.outhash = outhash;
            this.outrate = outrate;
            this.outprob = outprob;
            this.sn = sn;
        }
    }

    /**
     * Result class for reachableSpaceGenerator method
     */
    public static class reachableSpaceGeneratorResult {
        public Matrix SSq;
        public Matrix SSh;
        public jline.lang.NetworkStruct sn;

        public reachableSpaceGeneratorResult(Matrix SSq, Matrix SSh, jline.lang.NetworkStruct sn) {
            this.SSq = SSq;
            this.SSh = SSh;
            this.sn = sn;
        }
    }

    /**
     * Unified result type for distribution computations in queueing network solvers.
     * 
     * <p>This class provides a standardized interface for returning cumulative distribution
     * functions (CDFs) from various distribution analyses including:
     * <ul>
     *   <li>Response time distributions</li>
     *   <li>Passage time distributions</li>
     *   <li>Queue length distributions</li>
     *   <li>Transient distributions at specific time points</li>
     * </ul>
     * </p>
     * 
     * <p>The CDF data is organized as a 2D structure indexed by [station][class],
     * where each element contains a matrix of [quantile, value] pairs representing
     * the empirical CDF.</p>
     */
    public static class DistributionResult {
        
        /**
         * The CDF data organized as a cell array [stations x classes].
         * Each cell contains a matrix with [quantile, value] pairs.
         */
        public List<List<Matrix>> cdfData;
        
        /**
         * Number of stations in the model.
         */
        public int numStations;
        
        /**
         * Number of job classes in the model.
         */
        public int numClasses;
        
        /**
         * Type of distribution (e.g., "response_time", "passage_time").
         */
        public String distributionType;
        
        /**
         * Indicates whether this is a transient distribution.
         */
        public boolean isTransient;
        
        /**
         * Time points for transient distributions.
         */
        public Matrix timePoints;
        
        /**
         * Runtime for computing the distribution.
         */
        public double runtime;
        
        /**
         * Default constructor.
         */
        public DistributionResult() {
            this.cdfData = new ArrayList<>();
            this.numStations = 0;
            this.numClasses = 0;
            this.distributionType = "";
            this.isTransient = false;
            this.timePoints = new Matrix(0, 0);
            this.runtime = 0.0;
        }
        
        /**
         * Constructor for steady-state distributions.
         * 
         * @param numStations number of stations
         * @param numClasses number of classes
         * @param distributionType type of distribution
         */
        public DistributionResult(int numStations, int numClasses, String distributionType) {
            this();
            this.numStations = numStations;
            this.numClasses = numClasses;
            this.distributionType = distributionType;
            this.isTransient = false;
            
            // Initialize CDF data structure
            initializeCdfData();
        }
        
        /**
         * Constructor for transient distributions.
         * 
         * @param numStations number of stations
         * @param numClasses number of classes
         * @param distributionType type of distribution
         * @param timePoints time points for transient analysis
         */
        public DistributionResult(int numStations, int numClasses, String distributionType, Matrix timePoints) {
            this(numStations, numClasses, distributionType);
            this.isTransient = true;
            this.timePoints = timePoints.copy();
        }
        
        /**
         * Initializes the CDF data structure.
         */
        private void initializeCdfData() {
            cdfData = new ArrayList<>(numStations);
            for (int i = 0; i < numStations; i++) {
                List<Matrix> stationCdfs = new ArrayList<>(numClasses);
                for (int j = 0; j < numClasses; j++) {
                    stationCdfs.add(new Matrix(0, 0)); // Empty matrix initially
                }
                cdfData.add(stationCdfs);
            }
        }
        
        /**
         * Sets the CDF for a specific station and class.
         * 
         * @param station station index
         * @param jobClass class index
         * @param cdf CDF matrix with [quantile, value] pairs
         */
        public void setCdf(int station, int jobClass, Matrix cdf) {
            if (station < 0 || station >= numStations) {
                throw new IndexOutOfBoundsException("Station index out of bounds: " + station);
            }
            if (jobClass < 0 || jobClass >= numClasses) {
                throw new IndexOutOfBoundsException("Class index out of bounds: " + jobClass);
            }
            
            cdfData.get(station).set(jobClass, cdf.copy());
        }
        
        /**
         * Gets the CDF for a specific station and class.
         * 
         * @param station station index
         * @param jobClass class index
         * @return CDF matrix with [quantile, value] pairs
         */
        public Matrix getCdf(int station, int jobClass) {
            if (station < 0 || station >= numStations) {
                throw new IndexOutOfBoundsException("Station index out of bounds: " + station);
            }
            if (jobClass < 0 || jobClass >= numClasses) {
                throw new IndexOutOfBoundsException("Class index out of bounds: " + jobClass);
            }
            
            return cdfData.get(station).get(jobClass).copy();
        }
        
        /**
         * Sets an exponential approximation CDF based on mean value.
         * This is used as a default approximation when exact CDFs are not available.
         * 
         * @param station station index
         * @param jobClass class index
         * @param meanValue mean value for the exponential distribution
         */
        public void setExponentialApproximation(int station, int jobClass, double meanValue) {
            if (meanValue <= 0) {
                // Set degenerate distribution at 0
                Matrix cdf = new Matrix(1, 2);
                cdf.set(0, 0, 1.0); // quantile
                cdf.set(0, 1, 0.0); // value
                setCdf(station, jobClass, cdf);
            } else {
                // Create exponential CDF approximation
                int numPoints = 100;
                Matrix cdf = new Matrix(numPoints, 2);
                double lambda = 1.0 / meanValue;
                
                for (int i = 0; i < numPoints; i++) {
                    double quantile = 0.001 + (0.999 - 0.001) * i / (numPoints - 1);
                    double value = -Math.log(1 - quantile) / lambda;
                    cdf.set(i, 0, quantile);
                    cdf.set(i, 1, value);
                }
                
                setCdf(station, jobClass, cdf);
            }
        }
        
        /**
         * Checks if the CDF data is available for a specific station and class.
         * 
         * @param station station index
         * @param jobClass class index
         * @return true if CDF data is available
         */
        public boolean hasCdf(int station, int jobClass) {
            if (station < 0 || station >= numStations || jobClass < 0 || jobClass >= numClasses) {
                return false;
            }
            
            Matrix cdf = cdfData.get(station).get(jobClass);
            return cdf != null && !cdf.isEmpty();
        }
        
        /**
         * Gets the complete CDF data structure.
         * 
         * @return the CDF data as a list of lists of matrices
         */
        public List<List<Matrix>> getAllCdfData() {
            List<List<Matrix>> result = new ArrayList<>(numStations);
            for (int i = 0; i < numStations; i++) {
                List<Matrix> stationCdfs = new ArrayList<>(numClasses);
                for (int j = 0; j < numClasses; j++) {
                    stationCdfs.add(cdfData.get(i).get(j).copy());
                }
                result.add(stationCdfs);
            }
            return result;
        }
    }

    /**
     * Unified result type for probability computations in queueing network solvers.
     * 
     * <p>This class provides a standardized interface for returning probability values
     * from various probability analyses including:
     * <ul>
     *   <li>State probabilities (marginal or joint)</li>
     *   <li>Normalizing constants</li>
     *   <li>Aggregated state probabilities</li>
     *   <li>Node-specific marginal probabilities</li>
     * </ul>
     * </p>
     * 
     * <p>The result can represent scalar probabilities, probability vectors,
     * or probability matrices depending on the specific analysis performed.</p>
     */
    public static class ProbabilityResult {
        
        /**
         * The probability value(s). Can be a scalar, vector, or matrix depending on the analysis type.
         */
        public Matrix probability;
        
        /**
         * The logarithm of the normalizing constant (for methods that compute it).
         */
        public double logNormalizingConstant;
        
        /**
         * Indicates whether the result is for an aggregated state space.
         */
        public boolean isAggregated;
        
        /**
         * The node/station index for marginal probabilities (null for system-wide results).
         */
        public Integer nodeIndex;
        
        /**
         * The state specification for which the probability was computed.
         */
        public Matrix state;
        
        /**
         * Default constructor.
         */
        public ProbabilityResult() {
            this.probability = new Matrix(0, 0);
            this.logNormalizingConstant = Double.NaN;
            this.isAggregated = false;
            this.nodeIndex = null;
            this.state = null;
        }
        
        /**
         * Constructor for scalar probability results.
         * 
         * @param probability the probability value
         */
        public ProbabilityResult(double probability) {
            this();
            this.probability = new Matrix(1, 1);
            this.probability.set(0, 0, probability);
        }
        
        /**
         * Constructor for matrix probability results.
         * 
         * @param probability the probability matrix
         */
        public ProbabilityResult(Matrix probability) {
            this();
            this.probability = probability.copy();
        }
        
        /**
         * Constructor for normalizing constant results.
         * 
         * @param logNormalizingConstant the log of the normalizing constant
         */
        public ProbabilityResult(double logNormalizingConstant, boolean isLogNormConst) {
            this();
            this.logNormalizingConstant = logNormalizingConstant;
        }
        
        /**
         * Gets the probability as a scalar value.
         * 
         * @return the scalar probability value
         * @throws IllegalStateException if the result is not a scalar
         */
        public double getScalarProbability() {
            if (probability.getNumRows() != 1 || probability.getNumCols() != 1) {
                throw new IllegalStateException("Probability result is not a scalar");
            }
            return probability.get(0, 0);
        }
        
        /**
         * Gets the probability matrix.
         * 
         * @return the probability matrix
         */
        public Matrix getProbabilityMatrix() {
            return probability.copy();
        }
        
        /**
         * Checks if this result contains a valid probability value.
         * 
         * @return true if the probability is valid
         */
        public boolean hasValidProbability() {
            return probability != null && !probability.isEmpty() && 
                   !Double.isNaN(probability.get(0, 0));
        }
        
        /**
         * Checks if this result contains a valid normalizing constant.
         * 
         * @return true if the normalizing constant is valid
         */
        public boolean hasValidNormalizingConstant() {
            return !Double.isNaN(logNormalizingConstant);
        }
    }

    /**
     * Unified result type for sampling and simulation in queueing network solvers.
     * 
     * <p>This class provides a standardized interface for returning sampled trajectories
     * from various sampling methods including:
     * <ul>
     *   <li>Discrete-event simulation</li>
     *   <li>Perfect sampling</li>
     *   <li>Time-averaged sampling</li>
     *   <li>Single-node or system-wide sampling</li>
     * </ul>
     * </p>
     * 
     * <p>The state trajectory can be either a single matrix (for node-specific sampling)
     * or a list of matrices (for system-wide sampling), allowing flexible representation
     * of different sampling scenarios.</p>
     */
    public static class SampleResult {
        
        /**
         * Handle or identifier for the sampling method used.
         */
        public String handle;
        
        /**
         * Time points for the sampled trajectory.
         */
        public Matrix t;
        
        /**
         * State trajectory. For single node sampling, this is a matrix.
         * For system sampling, this is a list of matrices (one per node).
         */
        public Object state; // Can be Matrix or List<Matrix>
        
        /**
         * Event sequence corresponding to the state trajectory.
         */
        public Matrix event;
        
        /**
         * Indicates whether the result uses aggregated states.
         */
        public boolean isAggregate;
        
        /**
         * The node index for single-node sampling (null for system-wide sampling).
         */
        public Integer nodeIndex;
        
        /**
         * Number of samples generated.
         */
        public int numSamples;
        
        /**
         * Default constructor.
         */
        public SampleResult() {
            this.handle = "";
            this.t = new Matrix(0, 0);
            this.state = new Matrix(0, 0);
            this.event = new Matrix(0, 0);
            this.isAggregate = false;
            this.nodeIndex = null;
            this.numSamples = 0;
        }
        
        /**
         * Constructor for single-node sampling results.
         * 
         * @param handle the sampling method handle
         * @param t time points
         * @param state state trajectory matrix
         * @param event event sequence
         * @param isAggregate whether states are aggregated
         * @param nodeIndex the node index
         * @param numSamples number of samples
         */
        public SampleResult(String handle, Matrix t, Matrix state, Matrix event, 
                           boolean isAggregate, Integer nodeIndex, int numSamples) {
            this.handle = handle;
            this.t = t.copy();
            this.state = state.copy();
            this.event = event.copy();
            this.isAggregate = isAggregate;
            this.nodeIndex = nodeIndex;
            this.numSamples = numSamples;
        }
        
        /**
         * Constructor for system-wide sampling results.
         * 
         * @param handle the sampling method handle
         * @param t time points
         * @param systemState list of state trajectories (one per node)
         * @param event event sequence
         * @param isAggregate whether states are aggregated
         * @param numSamples number of samples
         */
        public SampleResult(String handle, Matrix t, List<Matrix> systemState, Matrix event, 
                           boolean isAggregate, int numSamples) {
            this.handle = handle;
            this.t = t.copy();
            this.state = new ArrayList<>(systemState);
            this.event = event.copy();
            this.isAggregate = isAggregate;
            this.nodeIndex = null; // System-wide sampling
            this.numSamples = numSamples;
        }
        
        /**
         * Gets the state trajectory as a matrix (for single-node sampling).
         * 
         * @return the state trajectory matrix
         * @throws IllegalStateException if this is a system-wide sampling result
         */
        public Matrix getStateMatrix() {
            if (!(state instanceof Matrix)) {
                throw new IllegalStateException("This is a system-wide sampling result, use getSystemStateList()");
            }
            return ((Matrix) state).copy();
        }
        
        /**
         * Gets the system state trajectories as a list of matrices (for system-wide sampling).
         * 
         * @return the list of state trajectory matrices
         * @throws IllegalStateException if this is a single-node sampling result
         */
        @SuppressWarnings("unchecked")
        public List<Matrix> getSystemStateList() {
            if (!(state instanceof List)) {
                throw new IllegalStateException("This is a single-node sampling result, use getStateMatrix()");
            }
            return new ArrayList<>((List<Matrix>) state);
        }
        
        /**
         * Checks if this is a system-wide sampling result.
         * 
         * @return true if this is a system-wide sampling result
         */
        public boolean isSystemWide() {
            return nodeIndex == null;
        }
        
        /**
         * Checks if this result contains valid sample data.
         * 
         * @return true if the result contains valid samples
         */
        public boolean hasValidSamples() {
            return t != null && !t.isEmpty() && state != null && numSamples > 0;
        }
    }

    /**
     * Result type for the Method of Moments (MoM) exact algorithm.
     * 
     * <p>The Method of Moments is an exact algorithm for computing the normalizing
     * constant and performance measures in product-form queueing networks. Unlike
     * floating-point algorithms, MoM uses exact rational arithmetic (BigFraction)
     * to avoid numerical errors, making it suitable for networks with extreme
     * parameter values or when exact results are required.</p>
     * 
     * <p>The algorithm computes normalizing constants recursively and derives
     * performance measures (throughputs and queue lengths) from these constants.</p>
     */
    public static class pfqnMom {
        public Matrix X;   // Throughputs
        public Matrix Q;   // Queue lengths
        public org.apache.commons.math3.fraction.BigFraction G;  // Exact normalizing constant
        public double lG;  // Natural logarithm of G
        public org.apache.commons.math3.fraction.BigFraction[] g;   // Final normalizing constants
        public org.apache.commons.math3.fraction.BigFraction[] g_1; // Pre-final normalizing constants

        public pfqnMom(Matrix X, Matrix Q, org.apache.commons.math3.fraction.BigFraction G, double lG,
                      org.apache.commons.math3.fraction.BigFraction[] g,
                      org.apache.commons.math3.fraction.BigFraction[] g_1) {
            this.X = X;
            this.Q = Q;
            this.G = G;
            this.lG = lG;
            this.g = g;
            this.g_1 = g_1;
        }
    }

    /**
     * Result type for the Akyildiz-Bolch (A/B) linearizer method for load-dependent multi-server BCMP networks.
     * 
     * <p>This class encapsulates the output of the pfqn_ab function which implements the 
     * Akyildiz-Bolch linearizer approach for analyzing multi-server queueing networks with 
     * processor sharing (PS) scheduling. The method uses weight functions and marginal 
     * probability estimation for iterative convergence.</p>
     */
    public static class pfqnAB {
        /** Queue length matrix (M x R) */
        public Matrix QN;
        /** Utilization matrix (M x R) */
        public Matrix UN;
        /** Response time matrix (M x R) */
        public Matrix RN;
        /** Throughput matrix (M x R) */
        public Matrix TN;
        /** Waiting time matrix (M x R) */
        public Matrix CN;
        /** System throughput vector (1 x R) */
        public Matrix XN;
        /** Method identifier */
        public String method;
        /** Number of iterations */
        public int iter;
        /** Runtime in seconds */
        public double runtime;

        public pfqnAB(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix CN, Matrix XN, 
                      String method, int iter, double runtime) {
            this.QN = QN;
            this.UN = UN;
            this.RN = RN;
            this.TN = TN;
            this.CN = CN;
            this.XN = XN;
            this.method = method;
            this.iter = iter;
            this.runtime = runtime;
        }
    }

    /**
     * Result type for the Schmidt method for load-dependent MVA with multi-server stations.
     * 
     * <p>This class encapsulates the output of the pfqn_schmidt function which implements 
     * the Schmidt population recursion algorithm for analyzing closed queueing networks 
     * with load-dependent service and multi-server stations. Supports INF, PS, and FCFS 
     * scheduling strategies with both single and multi-server configurations.</p>
     */
    public static class pfqnSchmidt {
        /** Queue length matrix (M x R) */
        public Matrix QN;
        /** Utilization matrix (M x R) */
        public Matrix UN;
        /** Response time matrix (M x R) */
        public Matrix RN;
        /** Throughput matrix (M x R) */
        public Matrix TN;
        /** Waiting time matrix (M x R) */
        public Matrix CN;
        /** System throughput vector (1 x R) */
        public Matrix XN;
        /** State probabilities per station */
        public List<Matrix> PN;
        /** Method identifier */
        public String method;
        /** Number of iterations (states enumerated) */
        public int iter;
        /** Runtime in seconds */
        public double runtime;

        public pfqnSchmidt(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix CN, Matrix XN,
                          List<Matrix> PN, String method, int iter, double runtime) {
            this.QN = QN;
            this.UN = UN;
            this.RN = RN;
            this.TN = TN;
            this.CN = CN;
            this.XN = XN;
            this.PN = PN;
            this.method = method;
            this.iter = iter;
            this.runtime = runtime;
        }
    }

    /**
     * Result type for Harel et al. throughput bounds for closed queueing networks.
     *
     * <p>This class encapsulates the output of the pfqn_harel_bounds function which
     * implements the throughput bounds from Harel, Namn, and Sturm (1999)
     * "Simple bounds for closed queueing networks", Queueing Systems 31:125-135.</p>
     *
     * <p>The bounds provide:
     * <ul>
     *   <li>LB: Lower bound on throughput (tighter than Zahorjan's balanced job bound)</li>
     *   <li>UB(n): Upper bounds on throughput for various n values</li>
     * </ul>
     * </p>
     */
    public static class pfqnHarelBounds {
        /** Lower bound on throughput */
        public double LB;
        /** Upper bounds UB(n) for n=2,3,...,maxN */
        public double[] UB;
        /** Exact throughput values TH(n) for n=1,2,...,maxN (used to compute UB) */
        public double[] TH;
        /** Population size N */
        public int N;
        /** Number of queues k */
        public int k;
        /** Maximum n for which UB(n) was computed */
        public int maxN;

        public pfqnHarelBounds(double LB, double[] UB, double[] TH, int N, int k, int maxN) {
            this.LB = LB;
            this.UB = UB;
            this.TH = TH;
            this.N = N;
            this.k = k;
            this.maxN = maxN;
        }
    }

    /**
     * Result type for the snToAG conversion.
     */
    public static class snToAG {
        public Matrix[][] R;
        public Matrix AP;
        public Matrix processMap;
        public List<ActionMapEntry> actionMap;
        public int[] N;

        public snToAG(Matrix[][] R, Matrix AP, Matrix processMap, List<ActionMapEntry> actionMap, int[] N) {
            this.R = R;
            this.AP = AP;
            this.processMap = processMap;
            this.actionMap = actionMap;
            this.N = N;
        }
    }

    /**
     * Entry for the action map in snToAG.
     */
    public static class ActionMapEntry {
        public int fromStation;
        public int fromClass;
        public int toStation;
        public int toClass;
        public double prob;
        public boolean isNegative;

        public ActionMapEntry(int fromStation, int fromClass, int toStation, int toClass, double prob, boolean isNegative) {
            this.fromStation = fromStation;
            this.fromClass = fromClass;
            this.toStation = toStation;
            this.toClass = toClass;
            this.prob = prob;
            this.isNegative = isNegative;
        }
    }

}
