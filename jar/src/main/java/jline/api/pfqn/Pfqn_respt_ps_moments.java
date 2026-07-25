/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.special.Gamma;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.ops.DConvertMatrixStruct;

/**
 * Sojourn-time moments at the processor-sharing station of the closed
 * terminal-driven system of Mitra and Morrison (1983): a bank of terminals in
 * series with a single processor-sharing CPU, with class-dependent exponential
 * think times (mean {@code Z[r]}) and class-dependent exponential service times
 * (mean {@code S[r]}), and {@code N[r]} jobs of class r cycling between the two.
 *
 * <p>Two routes to the moments are implemented, both from that paper.</p>
 *
 * <p>{@code exact} solves the linear system {@code c'[A - q_J I] = -pi'B} of
 * Proposition 3 on the state space {@code {n : 0 <= n <= K}}, K being the
 * population vector with the tagged class decremented by one. The moments are
 * then {@code E[W_J] = sum_n c(n)} and
 * {@code (q_J/2) E[W_J^2] = sum_n (n'1+1) c(n)}. Exact to solver precision, at
 * the cost of a linear solve of dimension {@code prod_r (K[r]+1)}.</p>
 *
 * <p>{@code asymptotic} evaluates the two leading terms of the asymptotic
 * expansion in inverse powers of the large parameter
 * {@code Nexp = max_r Z[r]/S[r]}, {@code E[W_J^2] ~ c0 + c1/Nexp}, of
 * Proposition 6. The cost is a linear system of dimension R, the number of
 * classes, and is therefore independent of the populations. Note that the
 * expansion parameter is the think-to-service ratio and NOT the population, so a
 * model with short think times is expanded in a small parameter no matter how
 * many jobs it holds.</p>
 *
 * <p>{@code auto} takes the exact route when the state space has at most
 * {@link #AUTO_MAX} states and the asymptotic route otherwise.</p>
 *
 * <p>The asymptotic route requires the normal-usage condition alpha &gt; 0.
 * Where it fails and the exact route is not affordable, the entry of W and W2 is
 * NaN and the result records "unavailable"; asking for "asymptotic" explicitly
 * in that regime throws rather than returning a blank.</p>
 *
 * <p>Reference: D. Mitra, J. A. Morrison, "Asymptotic Expansions of Moments of
 * the Waiting Time in Closed and Open Processor-Sharing Systems with Multiple
 * Job Classes", Adv. Appl. Prob. 15(4):813-839, 1983, Propositions 3 and 6.</p>
 *
 * <p>Port of MATLAB pfqn_respt_ps_moments.m.</p>
 */
public class Pfqn_respt_ps_moments {

    /** State-space size below which the auto route goes exact. */
    public static final int AUTO_MAX = 4096;
    /** Hard bound on an explicitly requested exact solve. */
    public static final int EXACT_MAX = 65536;

    private Pfqn_respt_ps_moments() {
    }

    /**
     * Sojourn-time moments at the PS station, choosing the route automatically.
     *
     * @param S per-class mean service times at the PS station, positive
     * @param N per-class populations, non-negative integers
     * @param Z per-class mean think times, positive where N &gt; 0
     * @return per-class moments and the route taken
     */
    public static PfqnResptPsResult pfqn_respt_ps_moments(double[] S, double[] N, double[] Z) {
        return pfqn_respt_ps_moments(S, N, Z, "auto");
    }

    /**
     * Sojourn-time moments at the PS station.
     *
     * @param S      per-class mean service times at the PS station, positive
     * @param N      per-class populations, non-negative integers
     * @param Z      per-class mean think times, positive where N &gt; 0
     * @param method "auto", "exact" or "asymptotic"
     * @return per-class moments and the route taken
     */
    public static PfqnResptPsResult pfqn_respt_ps_moments(double[] S, double[] N, double[] Z,
                                                          String method) {
        String mth = method == null ? "auto" : method.trim().toLowerCase();
        if (!mth.equals("auto") && !mth.equals("exact") && !mth.equals("asymptotic")) {
            throw new RuntimeException("pfqn_respt_ps_moments: method must be one of auto, exact, asymptotic");
        }
        int R = S.length;
        if (N.length != R || Z.length != R) {
            throw new RuntimeException("pfqn_respt_ps_moments: S, N and Z must have the same number of classes");
        }
        for (int r = 0; r < R; r++) {
            if (!isFinite(S[r]) || S[r] <= 0) {
                throw new RuntimeException("pfqn_respt_ps_moments: S must be finite and positive");
            }
            if (!isFinite(N[r]) || N[r] < 0 || N[r] != Math.rint(N[r])) {
                throw new RuntimeException("pfqn_respt_ps_moments: N must contain non-negative integers");
            }
            if (N[r] > 0 && (!isFinite(Z[r]) || Z[r] <= 0)) {
                throw new RuntimeException("pfqn_respt_ps_moments: Z must be finite and positive for every populated class");
            }
        }

        double[] W = new double[R];
        double[] W2 = new double[R];
        double[] c0 = new double[R];
        double[] c1 = new double[R];
        double[] alphaOut = new double[R];
        double[] nstates = new double[R];
        String[] mout = new String[R];
        for (int r = 0; r < R; r++) {
            W[r] = Double.NaN;
            W2[r] = Double.NaN;
            c0[r] = Double.NaN;
            c1[r] = Double.NaN;
            alphaOut[r] = Double.NaN;
            nstates[r] = Double.NaN;
            mout[r] = "none";
        }

        int na = 0;
        for (int r = 0; r < R; r++) {
            if (N[r] > 0) {
                na++;
            }
        }
        if (na == 0) {
            return new PfqnResptPsResult(W, W2, mout, c0, c1, alphaOut, nstates, Double.NaN);
        }
        int[] act = new int[na];
        int k = 0;
        for (int r = 0; r < R; r++) {
            if (N[r] > 0) {
                act[k++] = r;
            }
        }
        double[] qa = new double[na];
        double[] pa = new double[na];
        double expansionParam = 0.0;
        for (int j = 0; j < na; j++) {
            qa[j] = 1.0 / S[act[j]];
            pa[j] = 1.0 / Z[act[j]];
            expansionParam = Math.max(expansionParam, qa[j] / pa[j]);
        }

        for (int jj = 0; jj < na; jj++) {
            int J = act[jj];
            int[] K = new int[na];
            long ns = 1;
            for (int j = 0; j < na; j++) {
                K[j] = (int) Math.rint(N[act[j]]);
                if (j == jj) {
                    K[j] -= 1;
                }
                ns *= (K[j] + 1);
            }
            double lamSum = 0.0;
            for (int j = 0; j < na; j++) {
                lamSum += pa[j] * K[j] / qa[j];
            }
            double alpha = 1.0 - lamSum;
            alphaOut[J] = alpha;
            nstates[J] = ns;
            boolean useExact = mth.equals("exact") || (mth.equals("auto") && ns <= AUTO_MAX);
            if (useExact) {
                if (ns > EXACT_MAX) {
                    throw new RuntimeException(String.format(
                            "pfqn_respt_ps_moments: the exact route needs a linear solve of dimension %d, "
                                    + "above the bound of %d; use method = 'asymptotic'", ns, EXACT_MAX));
                }
                double[] wm = exactMoments(pa, qa, K, jj, (int) ns);
                W[J] = wm[0];
                W2[J] = wm[1];
                mout[J] = "exact";
                continue;
            }
            if (alpha <= 0) {
                if (mth.equals("asymptotic")) {
                    throw new RuntimeException(String.format(
                            "pfqn_respt_ps_moments: the asymptotic expansion needs normal usage alpha > 0, "
                                    + "but class %d gives alpha = %.6f", J, alpha));
                }
                mout[J] = "unavailable";
                continue;
            }
            double[] wm = asymptoticMoments(pa, qa, K, jj);
            W[J] = wm[0];
            W2[J] = wm[1];
            c0[J] = wm[2];
            c1[J] = wm[3];
            mout[J] = "asymptotic";
        }
        return new PfqnResptPsResult(W, W2, mout, c0, c1, alphaOut, nstates, expansionParam);
    }

    /**
     * Proposition 3: the moments follow from c, the solution of
     * c'[A - q_J I] = -pi'B, with A the generator-like operator of equation (26)
     * and B the diagonal operator B(n,n) = n'1+1.
     */
    private static double[] exactMoments(double[] p, double[] q, int[] K, int J, int ns) {
        int R = K.length;
        int[] dims = new int[R];
        int[] stride = new int[R];
        int acc = 1;
        for (int j = 0; j < R; j++) {
            dims[j] = K[j] + 1;
            stride[j] = acc;
            acc *= dims[j];
        }
        int[][] states = new int[ns][R];
        int[] tot = new int[ns];
        for (int i = 0; i < ns; i++) {
            int res = i;
            int s = 0;
            for (int j = 0; j < R; j++) {
                states[i][j] = res % dims[j];
                res /= dims[j];
                s += states[i][j];
            }
            tot[i] = s;
        }

        // stationary law (15), in logs so that large populations do not overflow
        double[] logpi = new double[ns];
        double maxlog = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < ns; i++) {
            double v = Gamma.logGamma(tot[i] + 1.0);
            for (int j = 0; j < R; j++) {
                int nj = states[i][j];
                v += Gamma.logGamma(K[j] + 1.0) - Gamma.logGamma(nj + 1.0)
                        - Gamma.logGamma(K[j] - nj + 1.0);
                double rj = p[j] / q[j];
                if (rj > 0) {
                    v += nj * Math.log(rj);
                } else if (nj > 0) {
                    v = Double.NEGATIVE_INFINITY;
                }
            }
            logpi[i] = v;
            if (v > maxlog) {
                maxlog = v;
            }
        }
        double sumpi = 0.0;
        double[] pin = new double[ns];
        for (int i = 0; i < ns; i++) {
            pin[i] = Math.exp(logpi[i] - maxlog);
            sumpi += pin[i];
        }
        for (int i = 0; i < ns; i++) {
            pin[i] /= sumpi;
        }

        DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(ns, ns, ns * (2 * R + 1));
        for (int i = 0; i < ns; i++) {
            double diag = -q[J];
            for (int j = 0; j < R; j++) {
                int nj = states[i][j];
                if (nj >= 1) {
                    triplet.addItem(i - stride[j], i, p[j] * (K[j] - nj + 1) * tot[i]);
                }
                if (nj <= K[j] - 1) {
                    triplet.addItem(i + stride[j], i, (nj + 1) * q[j]);
                }
                diag -= p[j] * (K[j] - nj) * (tot[i] + 1) + nj * q[j];
            }
            triplet.addItem(i, i, diag);
        }
        DMatrixSparseCSC csc = DConvertMatrixStruct.convert(triplet, (DMatrixSparseCSC) null);
        Matrix A = new Matrix((org.ejml.data.DMatrix) csc);

        Matrix rhs = new Matrix(new DMatrixRMaj(ns, 1));
        for (int i = 0; i < ns; i++) {
            rhs.set(i, 0, -(tot[i] + 1.0) * pin[i]);
        }
        Matrix c = new Matrix(new DMatrixRMaj(ns, 1));
        Matrix.solve(A.transpose(), rhs, c);

        double W = 0.0;
        double sw2 = 0.0;
        for (int i = 0; i < ns; i++) {
            W += c.get(i, 0);
            sw2 += (tot[i] + 1.0) * c.get(i, 0);
        }
        return new double[]{W, 2.0 / q[J] * sw2};
    }

    /**
     * Proposition 6: the two leading terms of the expansion in 1/Nexp. Equation
     * numbers below are those of Mitra and Morrison (1983).
     */
    private static double[] asymptoticMoments(double[] p, double[] q, int[] K, int J) {
        int R = K.length;
        double[] lam = new double[R];
        double lamSum = 0.0;
        double Nexp = 0.0;
        for (int j = 0; j < R; j++) {
            lam[j] = p[j] * K[j];
            lamSum += lam[j] / q[j];
            Nexp = Math.max(Nexp, q[j] / p[j]);                       // (50)
        }
        double alpha = 1.0 - lamSum;
        double[] Gam = new double[R];
        double[] beta = new double[R];
        double bg2 = 0.0;
        for (int j = 0; j < R; j++) {
            Gam[j] = Nexp * p[j] / q[j];                              // (51)
            beta[j] = K[j] / Nexp;                                    // (51)
            bg2 += beta[j] * Gam[j] * Gam[j];
        }
        double qJ = q[J];

        double den = 1.0;
        double numer = 1.0;
        for (int j = 0; j < R; j++) {
            den -= lam[j] / (q[j] + qJ);
            numer -= lam[j] * (q[j] - qJ) / (q[j] * (q[j] + qJ));
        }
        double F10 = (-1.0 / (alpha * alpha * qJ)) * numer / den;      // (110)
        double c0 = -2.0 / qJ * F10;

        double[] f1 = new double[R];
        double[] S2j = new double[R];
        for (int j = 0; j < R; j++) {
            f1[j] = lam[j] / (q[j] + qJ) * (F10 - 2.0 / (alpha * alpha * q[j]));      // (113iii)
            S2j[j] = 6.0 / Math.pow(alpha, 4) * (alpha * beta[j] * Gam[j] * Gam[j]
                    + 2.0 * bg2 * beta[j] * Gam[j]);                                  // (113i)
        }

        Matrix Amat = new Matrix(new DMatrixRMaj(R, R));
        Matrix rhs = new Matrix(new DMatrixRMaj(R, 1));
        for (int j = 0; j < R; j++) {
            double diag = 1.0;
            double b = -f1[j];
            for (int s = 0; s < R; s++) {
                double d = q[j] + q[s] + qJ;
                double s2js = 3.0 / Math.pow(alpha, 3) * (beta[j] * Gam[j]) * (beta[s] * Gam[s]); // (113ii)
                b += s2js / d;
                diag -= lam[s] / d;
                if (s != j) {
                    Amat.set(j, s, -lam[j] / d);
                }
            }
            diag -= lam[j] / (2 * q[j] + qJ);
            Amat.set(j, j, diag);
            rhs.set(j, 0, b);
        }
        Matrix F2 = new Matrix(new DMatrixRMaj(R, 1));
        Matrix.solve(Amat, rhs, F2);                                                  // (112)

        double f10 = -3.0 / (Math.pow(alpha, 3) * qJ) * bg2;                          // (98)
        double acc = 0.0;
        for (int j = 0; j < R; j++) {
            acc += (2.0 * Gam[j] * q[j] * F2.get(j, 0) + S2j[j]) / (q[j] + qJ);
        }
        double F20 = (acc - f10) / den;                                               // (111)
        double c1 = -2.0 / qJ * F20 + c0 / (alpha * alpha) * bg2;                      // (114ii)

        double W = 1.0 / (alpha * qJ) * (1.0 - 2.0 / Nexp * bg2 / (alpha * alpha));    // (68)
        double W2 = c0 + c1 / Nexp;                                                    // (114i)
        return new double[]{W, W2, c0, c1};
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }
}
