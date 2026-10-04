/**
 * @file Loss Network Analysis via Monte Carlo Importance-Sampling Summation
 *
 * @since LINE 3.0
 */
package jline.api.lossn;

import java.util.Random;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

/**
 * Monte Carlo importance-sampling summation for loss networks, after
 * Ross and Wang, "Monte Carlo Summation Applied to Product-Form Loss
 * Networks", Probability in the Engineering and Informational Sciences,
 * 6 (1992), 323-348.
 *
 * Estimates the product-form normalization constant g(C) and per-class
 * blocking probabilities with confidence intervals. Unlike the Erlang
 * fixed point, no link-independence assumption is made and the estimator
 * is consistent.
 */
public final class Lossn_mci {
    private Lossn_mci() {}

    /**
     * Monte Carlo importance-sampling summation for loss networks.
     *
     * @param nuVec   Offered load per class (1xR).
     * @param Amat    Circuit requirement of link j for class r (JxR).
     * @param cVec    Link capacity (Jx1).
     * @param samples Number of Monte Carlo samples.
     * @param gammaVec Importance-sampling parameters (1xR), or null for the
     *                 Section 3.4 heuristic.
     * @param seed    RNG seed; use a negative value for a nondeterministic seed.
     * @param alpha   Confidence-interval significance level (e.g. 0.05).
     */
    public static Ret.lossnMCI lossn_mci(Matrix nuVec, Matrix Amat, Matrix cVec,
                                         int samples, Matrix gammaVec, long seed, double alpha) {
        double[] nu = nuVec.toArray1D();
        double[][] A = Amat.toArray2D();
        double[] C = cVec.toArray1D();
        int R = nu.length;
        int J = C.length;
        if (A.length != J || A[0].length != R) {
            throw new RuntimeException("A must be J x R (" + J + " x " + R + ")");
        }
        int S = samples;
        Random rng = (seed < 0) ? new Random() : new Random(seed);

        // Per-class maximum feasible occupancy N_r = min_j floor(C_j/A_jr)
        int[] N = new int[R];
        for (int k = 0; k < R; k++) {
            double best = Double.POSITIVE_INFINITY;
            boolean any = false;
            for (int j = 0; j < J; j++) {
                if (A[j][k] > 0) {
                    any = true;
                    best = FastMath.min(best, FastMath.floor(C[j] / A[j][k]));
                }
            }
            N[k] = any ? (int) best : 0;
        }

        // Importance-sampling parameters gamma (Section 3.4 heuristic)
        double[] gamma = new double[R];
        if (gammaVec != null) {
            gamma = gammaVec.toArray1D();
        } else {
            double delta = 0.0;
            for (int j = 0; j < J; j++) {
                double load = 0.0;
                for (int k = 0; k < R; k++) {
                    load += A[j][k] * nu[k];
                }
                delta = FastMath.max(delta, load / C[j]);
            }
            double base = FastMath.max(1.0 - 0.15 * (1.0 - delta), 1e-6);
            for (int k = 0; k < R; k++) {
                double bk = 0.0;
                for (int j = 0; j < J; j++) {
                    bk = FastMath.max(bk, A[j][k]);
                }
                gamma[k] = nu[k] * FastMath.pow(base, bk);
            }
        }
        for (int k = 0; k < R; k++) {
            gamma[k] = FastMath.max(gamma[k], 1e-300);
        }

        // Normalization constant c of the importance distribution (log space)
        // and per-class sampling CDFs.
        double log_c = 0.0;
        double[][] cdf = new double[R][];
        for (int k = 0; k < R; k++) {
            double[] logpmf = new double[N[k] + 1];
            for (int l = 0; l <= N[k]; l++) {
                logpmf[l] = l * FastMath.log(gamma[k]) - Maths.factln(l);
            }
            double lse = logsumexp(logpmf);
            log_c += lse;
            double[] ck = new double[N[k] + 1];
            double acc = 0.0;
            for (int l = 0; l <= N[k]; l++) {
                acc += FastMath.exp(logpmf[l] - lse);
                ck[l] = acc;
            }
            ck[N[k]] = 1.0;
            cdf[k] = ck;
        }

        double[] logratio = new double[R];
        for (int k = 0; k < R; k++) {
            logratio[k] = FastMath.log(nu[k]) - FastMath.log(gamma[k]);
        }

        // see _kb/03-api-layer.md for rationale (lossn/ section)
        int[][] V = new int[S][R];
        double[] logAlpha = new double[S];
        boolean[] inOmega = new boolean[S];
        boolean[][] inOmegaK = new boolean[S][R];
        double maxLaO = Double.NEGATIVE_INFINITY;
        boolean anyOmega = false;

        for (int s = 0; s < S; s++) {
            double la = 0.0;
            for (int k = 0; k < R; k++) {
                double u = rng.nextDouble();
                double[] ck = cdf[k];
                int v = 0;
                while (v < ck.length && u > ck[v]) {
                    v++;
                }
                if (v >= ck.length) {
                    v = ck.length - 1;
                }
                V[s][k] = v;
                la += v * logratio[k];
            }
            logAlpha[s] = la;
            // Feasibility: A*n <= C and A*n <= C - A_.k
            boolean feas = true;
            double[] AV = new double[J];
            for (int j = 0; j < J; j++) {
                double val = 0.0;
                for (int k = 0; k < R; k++) {
                    val += A[j][k] * V[s][k];
                }
                AV[j] = val;
                if (val > C[j]) {
                    feas = false;
                }
            }
            inOmega[s] = feas;
            for (int k = 0; k < R; k++) {
                boolean feasK = true;
                for (int j = 0; j < J; j++) {
                    if (AV[j] > C[j] - A[j][k]) {
                        feasK = false;
                        break;
                    }
                }
                inOmegaK[s][k] = feasK;
            }
            if (feas) {
                anyOmega = true;
                if (la > maxLaO) {
                    maxLaO = la;
                }
            }
        }

        // Normalization constant estimate g(C) = (c/S) sum alpha 1(Omega)
        double lG;
        if (!anyOmega) {
            lG = Double.NEGATIVE_INFINITY;
        } else {
            double acc = 0.0;
            for (int s = 0; s < S; s++) {
                if (inOmega[s]) {
                    acc += FastMath.exp(logAlpha[s] - maxLaO);
                }
            }
            lG = log_c + maxLaO + FastMath.log(acc) - FastMath.log(S);
        }

        // Ratio estimators with shifted weights for stability
        double M = anyOmega ? maxLaO : 0.0;
        double[] w = new double[S];
        double[] Z = new double[S];
        double sumZ = 0.0;
        for (int s = 0; s < S; s++) {
            w[s] = FastMath.exp(logAlpha[s] - M);
            Z[s] = inOmega[s] ? w[s] : 0.0;
            sumZ += Z[s];
        }
        double meanZ = sumZ / S;
        double varZ = sampleVar(Z, meanZ, S);

        double crit = FastMath.sqrt(2.0) * Erf.erfInv(1.0 - alpha);
        double[] accept = new double[R];
        double[] loss = new double[R];
        double[] qLen = new double[R];
        double[][] accCI = new double[R][2];
        double[][] lossCI = new double[R][2];

        for (int k = 0; k < R; k++) {
            double sumY = 0.0;
            double[] Y = new double[S];
            for (int s = 0; s < S; s++) {
                Y[s] = inOmegaK[s][k] ? w[s] : 0.0;
                sumY += Y[s];
            }
            double meanY = sumY / S;
            double phi;
            double half;
            if (meanZ <= 0) {
                phi = Double.NaN;
                half = Double.NaN;
            } else {
                phi = meanY / meanZ;
                double varY = sampleVar(Y, meanY, S);
                double cov = 0.0;
                for (int s = 0; s < S; s++) {
                    cov += (Y[s] - meanY) * (Z[s] - meanZ);
                }
                cov /= (S - 1);
                double sig2 = (varY - 2 * phi * cov + phi * phi * varZ) / (S * meanZ * meanZ);
                sig2 = FastMath.max(sig2, 0.0);
                half = crit * FastMath.sqrt(sig2);
            }
            accept[k] = phi;
            loss[k] = 1.0 - phi;
            qLen[k] = nu[k] * phi;
            accCI[k][0] = phi - half;
            accCI[k][1] = phi + half;
            lossCI[k][0] = 1.0 - (phi + half);
            lossCI[k][1] = 1.0 - (phi - half);
        }

        return new Ret.lossnMCI(new Matrix(qLen), new Matrix(loss), lG,
                new Matrix(accCI), new Matrix(lossCI), S);
    }

    private static double logsumexp(double[] x) {
        double m = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < x.length; i++) {
            if (x[i] > m) {
                m = x[i];
            }
        }
        if (Double.isInfinite(m)) {
            return m;
        }
        double s = 0.0;
        for (int i = 0; i < x.length; i++) {
            s += FastMath.exp(x[i] - m);
        }
        return m + FastMath.log(s);
    }

    private static double sampleVar(double[] x, double mean, int S) {
        double acc = 0.0;
        for (int s = 0; s < S; s++) {
            double d = x[s] - mean;
            acc += d * d;
        }
        return acc / (S - 1);
    }
}
