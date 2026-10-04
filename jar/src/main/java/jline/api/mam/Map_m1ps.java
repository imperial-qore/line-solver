/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mam;

import org.apache.commons.math3.special.Gamma;

import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;

/**
 * Sojourn time distribution of the MAP/M/1 processor-sharing queue.
 *
 * <p>Port of {@code matlab/src/api/mam/map_compute_R.m},
 * {@code map_m1ps_h_recursive.m} and {@code map_m1ps_sojourn.m}, twin of the
 * native Python {@code api.mam.mapm1ps}. The processor-sharing discipline
 * shares the server equally, so with n customers present each is served at rate
 * mu/n; the sojourn law is therefore not the FCFS waiting law and is not
 * available in closed form.
 *
 * <p>Implements Theorem 1 of Masuyama, H. and Takine, T., "Sojourn time
 * distribution in a MAP/M/1 processor-sharing queue", Operations Research
 * Letters 31(6), 2003, 406-412:
 * <pre>
 *   Wbar(x) = (1/lambda) sum_n pi_0 R^n D sum_k p_k(x) h_{n,k}
 * </pre>
 * with p_k(x) the Poisson(( theta + mu) x) mass, R the minimal nonnegative
 * solution of D + R(C - mu I) + mu R^2 = 0, and h_{n,k} the uniformized
 * coefficients of {@link #map_m1ps_h_recursive}.
 */
public final class Map_m1ps {

    private Map_m1ps() {
    }

    /** Result of {@link #map_m1ps_sojourn}. */
    public static final class SojournResult {
        /** 1 x numel(x), Pr[W &gt; x]. */
        public Matrix Wbar;
        /**
         * (N+1) x numel(x), the CCDF conditional on the arrival finding n
         * customers already in the system, n = 0..N. The phase law at such an
         * arrival epoch is pi_0 R^n D normalized, so
         * sum_n P(N=n) WbarN(n,.) = Wbar with P(N=n) = pi_0 R^n D e / lambda.
         */
        public Matrix WbarN;
        /** Queue-length truncation order actually used. */
        public int nEpsilon;
    }

    /**
     * Minimal nonnegative solution R of D + R(C - mu I) + mu R^2 = 0, by the
     * fixed point R &lt;- -D (C - mu I + mu R)^-1 started at R = 0.
     */
    public static Matrix map_compute_R(Matrix C, Matrix D, double mu) {
        int M = C.getNumRows();
        Matrix R = new Matrix(M, M);
        int maxIter = 1000;
        double tol = 1e-10;
        for (int iter = 0; iter < maxIter; iter++) {
            Matrix Rold = R.copy();
            Matrix A = new Matrix(M, M);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    A.set(i, j, C.get(i, j) - (i == j ? mu : 0.0) + mu * R.get(i, j));
                }
            }
            Matrix Rn = D.mult(A.inv());
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    Rn.set(i, j, -Rn.get(i, j));
                }
            }
            R = Rn;
            double dev = 0.0;
            for (int i = 0; i < M; i++) {
                double rowdev = 0.0;
                for (int j = 0; j < M; j++) {
                    rowdev += Math.abs(R.get(i, j) - Rold.get(i, j));
                }
                dev = Math.max(dev, rowdev);
            }
            if (dev < tol) {
                break;
            }
        }
        // numerical error can push an entry a hair below zero; R is nonnegative
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                if (R.get(i, j) < 0) {
                    R.set(i, j, 0.0);
                }
            }
        }
        return R;
    }

    /**
     * The h_{n,k} vectors of Theorem 1, returned as h[n][k], each M x 1:
     * <pre>
     *   h_{n,0}   = e
     *   h_{n,k+1} = [ n mu/(n+1) h_{n-1,k} + (theta I + C) h_{n,k} + D h_{n+1,k} ]
     *               / (theta + mu)
     * </pre>
     * with h_{-1,k} = 0 and the level dimension truncated at N, i.e.
     * h_{N+1,k} = 0.
     */
    public static double[][][] map_m1ps_h_recursive(Matrix C, Matrix D, double mu, int N, int K) {
        int M = C.getNumRows();
        double theta = 0.0;
        for (int i = 0; i < M; i++) {
            theta = Math.max(theta, Math.abs(C.get(i, i)));
        }
        return hRecursive(C, D, mu, N, K, theta);
    }

    private static double[][][] hRecursive(Matrix C, Matrix D, double mu, int N, int K,
                                           double theta) {
        int M = C.getNumRows();
        double thetaPlusMu = theta + mu;
        double[][] tIC = new double[M][M];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                tIC[i][j] = (i == j ? theta : 0.0) + C.get(i, j);
            }
        }
        double[][][] h = new double[N + 1][K + 1][M];
        for (int n = 0; n <= N; n++) {
            for (int i = 0; i < M; i++) {
                h[n][0][i] = 1.0;
            }
        }
        for (int k = 0; k < K; k++) {
            for (int n = 0; n <= N; n++) {
                for (int i = 0; i < M; i++) {
                    double v = 0.0;
                    for (int j = 0; j < M; j++) {
                        v += tIC[i][j] * h[n][k][j];
                    }
                    if (n > 0) {
                        v += (n * mu / (n + 1.0)) * h[n - 1][k][i];
                    }
                    if (n < N) {
                        for (int j = 0; j < M; j++) {
                            v += D.get(i, j) * h[n + 1][k][j];
                        }
                    }
                    h[n][k + 1][i] = v / thetaPlusMu;
                }
            }
        }
        return h;
    }

    /** map_m1ps_sojourn at the reference defaults epsilon = 1e-11, epsilon' = 1e-10. */
    public static SojournResult map_m1ps_sojourn(Matrix C, Matrix D, double mu, double[] x) {
        return map_m1ps_sojourn(C, D, mu, x, 1e-11, 1e-10);
    }

    /**
     * @param C            M x M MAP transitions without an arrival
     * @param D            M x M MAP transitions with an arrival
     * @param mu           service rate, mu &gt; 0
     * @param x            time points at which Pr[W &gt; x] is evaluated
     * @param epsilon      queue-length truncation tolerance
     * @param epsilonPrime uniformization truncation tolerance
     */
    public static SojournResult map_m1ps_sojourn(Matrix C, Matrix D, double mu, double[] x,
                                                 double epsilon, double epsilonPrime) {
        int M = C.getNumRows();
        if (C.getNumCols() != M || D.getNumRows() != M || D.getNumCols() != M) {
            line_error("map_m1ps_sojourn", "C and D must be square matrices of equal size.");
        }
        if (!(mu > 0)) {
            line_error("map_m1ps_sojourn", "the service rate mu must be positive.");
        }
        for (int i = 0; i < x.length; i++) {
            if (x[i] < 0) {
                line_error("map_m1ps_sojourn", "the time points must be nonnegative.");
            }
        }

        // Stationary phase law of the MAP: the LEFT null vector of C+D.
        Matrix pi = Map_prob.map_prob(C, D);
        double lambda = 0.0;
        for (int i = 0; i < M; i++) {
            double rowsum = 0.0;
            for (int j = 0; j < M; j++) {
                rowsum += D.get(i, j);
            }
            lambda += pi.get(i) * rowsum;
        }
        double rho = lambda / mu;
        if (rho >= 1.0) {
            line_error("map_m1ps_sojourn", "the system is unstable (rho = " + rho + " >= 1).");
        }

        Matrix R = map_compute_R(C, D, mu);
        // pi_0 = pi (I - R)
        double[] pi0 = new double[M];
        for (int j = 0; j < M; j++) {
            double s = pi.get(j);
            for (int i = 0; i < M; i++) {
                s -= pi.get(i) * R.get(i, j);
            }
            pi0[j] = s;
        }

        // N(epsilon): the smallest N with (1/lambda) sum_{n<=N} pi_0 R^n D e > 1 - epsilon
        int nEpsilon = 100;
        double cum = 0.0;
        Matrix Rpow = Matrix.eye(M);
        for (int n = 0; n <= 1000; n++) {
            double[] w = rowTimes(pi0, Rpow);
            double[] wd = rowTimes(w, D);
            double s = 0.0;
            for (int j = 0; j < M; j++) {
                s += wd[j];
            }
            cum += s / lambda;
            if (cum > 1 - epsilon) {
                nEpsilon = n;
                break;
            }
            Rpow = Rpow.mult(R);
        }

        double theta = 0.0;
        for (int i = 0; i < M; i++) {
            theta = Math.max(theta, Math.abs(C.get(i, i)));
        }
        double thetaPlusMu = theta + mu;

        Matrix Wbar = new Matrix(1, x.length);
        Matrix WbarN = new Matrix(nEpsilon + 1, x.length);

        for (int idx = 0; idx < x.length; idx++) {
            double xv = x[idx];
            double meanVal = thetaPlusMu * xv;
            int kLo;
            int kMax;
            if (meanVal > 0) {
                kLo = (int) Math.max(0.0, Math.floor(meanVal - 10.0 * Math.sqrt(meanVal)));
                kMax = (int) Math.ceil(meanVal + 10.0 * Math.sqrt(meanVal));
                while (poissonMass(kLo, kMax, meanVal) < 1 - epsilonPrime && kMax < 10000) {
                    kMax += 10;
                }
            } else {
                kLo = 0;
                kMax = 0;
            }

            double[][][] h = hRecursive(C, D, mu, nEpsilon, kMax, theta);

            Matrix Rp = Matrix.eye(M);
            double acc = 0.0;
            for (int n = 0; n <= nEpsilon; n++) {
                double[] weight = rowTimes(rowTimes(pi0, Rp), D);
                double[] sumK = new double[M];
                for (int k = kLo; k <= kMax; k++) {
                    double p = poissonPmf(k, meanVal);
                    if (p == 0.0) {
                        continue;
                    }
                    for (int i = 0; i < M; i++) {
                        sumK[i] += p * h[n][k][i];
                    }
                }
                double dot = 0.0;
                double wsum = 0.0;
                for (int i = 0; i < M; i++) {
                    dot += weight[i] * sumK[i];
                    wsum += weight[i];
                }
                acc += dot / lambda;
                WbarN.set(n, idx, (wsum > 0) ? dot / wsum : 0.0);
                Rp = Rp.mult(R);
            }
            Wbar.set(0, idx, acc);
        }

        SojournResult out = new SojournResult();
        out.Wbar = Wbar;
        out.WbarN = WbarN;
        out.nEpsilon = nEpsilon;
        return out;
    }

    /** Row vector times matrix. */
    private static double[] rowTimes(double[] v, Matrix A) {
        int m = A.getNumRows();
        int n = A.getNumCols();
        double[] out = new double[n];
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int i = 0; i < m; i++) {
                s += v[i] * A.get(i, j);
            }
            out[j] = s;
        }
        return out;
    }

    /** Poisson pmf via log-gamma, so a large k does not overflow the factorial. */
    private static double poissonPmf(int k, double mean) {
        if (mean <= 0.0) {
            return (k == 0) ? 1.0 : 0.0;
        }
        return Math.exp(k * Math.log(mean) - Gamma.logGamma(k + 1.0) - mean);
    }

    private static double poissonMass(int lo, int hi, double mean) {
        double s = 0.0;
        for (int k = lo; k <= hi; k++) {
            s += poissonPmf(k, mean);
        }
        return s;
    }
}
