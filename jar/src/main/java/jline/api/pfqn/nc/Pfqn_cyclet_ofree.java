/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.List;
import java.util.function.UnaryOperator;

import org.apache.commons.math3.complex.Complex;

import jline.api.lti.Laplace_invert;

/**
 * Exact passage-time density, distribution and moments along an OVERTAKE-FREE
 * PATH of a closed single-chain tree-like product-form network.
 *
 * <p>Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time
 * Distributions in Large Markov Chains", 2002, Sec. 7.1, Theorems 1 and 2, after
 * P. G. Harrison, J. Appl. Prob. 27, 1990 and H. Duduna, Adv. Appl. Prob. 14,
 * 1982. The underlying sojourn-time result for overtake-free paths is F. Kelly
 * and P. Pollett, Adv. Appl. Prob. 15, 1983.
 *
 * <p>THE ONE FACT THAT MAKES ALL THREE ROUTES WORK. Conditional on the path,
 *
 * <pre>
 *     T | z  =  sum_{j in z} Erlang(u_{z_j} + 1, mu_{z_j})
 * </pre>
 *
 * with u distributed as the network's equilibrium population vector AT N-1 (the
 * arrival theorem). Hence the transform of Theorem 1 collapses to
 *
 * <pre>
 *     L(s|z) = prod_{j in z} mu_j/(s+mu_j) * G(y(s), N-1) / G(x, N-1)
 * </pre>
 *
 * where x_i = v_i/mu_i and y_i(s) = x_i mu_i/(s+mu_i) on the path, x_i off it:
 * one Buzen convolution per value of s.
 *
 * <p>MOMENTS ARE NEVER TAKEN FROM THE DENSITY. They come from running the same
 * Buzen convolution in the ring of truncated power series in s, so they are
 * exact to machine precision, are unaffected by the time grid, and stay valid
 * when the rates coincide and Theorem 2 does not apply.
 *
 * <p>NOTE ON THE PAPER. The inner sum of Theorem 2 reads (v_j t)^(c-i)/(c-i)!
 * and that is CORRECT as printed, however odd the visit ratio looks against a
 * time: substituting the service rate instead returns negative densities.
 * Verified against a direct mixture-of-Erlangs oracle to 1e-15, and at the
 * paper's own N = 18 example against the transform route to 1e-11.
 *
 * <p>Node indices are 0-based here and 1-based in the MATLAB reference.
 */
public final class Pfqn_cyclet_ofree {

    private Pfqn_cyclet_ofree() {
    }

    /** Default separation below which two path rates count as coincident. */
    public static final double DEFAULT_TOL = 1e-8;

    public static PfqnCycletResult pfqn_cyclet_ofree(double[] v, double[] mu, int N,
                                                     List<int[]> paths, double[] tset) {
        return pfqn_cyclet_ofree(v, mu, N, paths, tset, "auto", 3, null, "euler", DEFAULT_TOL);
    }

    /**
     * @param method "auto" (default) uses "exact" when the path rates are
     *               separated and "lt" otherwise; "exact" is Theorem 2 in closed
     *               form and REQUIRES DISTINCT RATES on the path, since its
     *               partial fractions divide by prod_{i!=j}(mu_i - mu_j)
     */
    public static PfqnCycletResult pfqn_cyclet_ofree(double[] v, double[] mu, int N,
                                                     List<int[]> paths, double[] tset,
                                                     String method, int nmom, double[] pathprob,
                                                     String ltiMethod, double tol) {
        int M = v.length;
        if (mu.length != M) {
            throw new RuntimeException(
                    "pfqn_cyclet_ofree: v and mu must name the same number of nodes");
        }
        for (int i = 0; i < M; i++) {
            if (!(mu[i] > 0.0)) {
                throw new RuntimeException(
                        "pfqn_cyclet_ofree: every service rate must be positive");
            }
        }
        if (N < 1) {
            throw new RuntimeException("pfqn_cyclet_ofree: the population N must be positive");
        }
        if (paths == null || paths.isEmpty()) {
            throw new RuntimeException("pfqn_cyclet_ofree: no path was given");
        }
        double[] pp = pathprob;
        if (pp == null) {
            pp = new double[paths.size()];
            for (int i = 0; i < pp.length; i++) {
                pp[i] = 1.0 / paths.size();
            }
        }
        if (pp.length != paths.size()) {
            throw new RuntimeException(
                    "pfqn_cyclet_ofree: pathprob must carry one probability per path");
        }

        double[] x = new double[M];
        for (int i = 0; i < M; i++) {
            x[i] = v[i] / mu[i];
        }
        double[] gfull = buzen(x, N - 1);
        double Gn1 = gfull[N - 1];
        if (!(Gn1 > 0.0)) {
            throw new RuntimeException("pfqn_cyclet_ofree: the network normalizing constant at "
                    + "population N-1 vanished; check v and mu");
        }

        double[] f = new double[tset.length];
        double[] F = new double[tset.length];
        double[] mom = new double[nmom];
        List<String> methods = new ArrayList<String>();

        for (int ip = 0; ip < paths.size(); ip++) {
            int[] z = paths.get(ip);
            if (z.length == 0) {
                throw new RuntimeException("pfqn_cyclet_ofree: an overtake-free path must contain "
                        + "at least the root node");
            }
            for (int a = 0; a < z.length; a++) {
                if (z[a] < 0 || z[a] >= M) {
                    throw new RuntimeException(
                            "pfqn_cyclet_ofree: a path node is outside the network");
                }
                for (int b = a + 1; b < z.length; b++) {
                    if (z[a] == z[b]) {
                        throw new RuntimeException(
                                "pfqn_cyclet_ofree: a path must have distinct nodes");
                    }
                }
            }
            int m = z.length;
            boolean[] onpath = new boolean[M];
            for (int j : z) {
                onpath[j] = true;
            }

            String mth = method;
            if ("auto".equalsIgnoreCase(mth)) {
                if (m == 1) {
                    mth = "exact";
                } else {
                    double sep = Double.POSITIVE_INFINITY;
                    double mx = 0.0;
                    for (int a = 0; a < m; a++) {
                        mx = Math.max(mx, mu[z[a]]);
                        for (int b = a + 1; b < m; b++) {
                            sep = Math.min(sep, Math.abs(mu[z[a]] - mu[z[b]]));
                        }
                    }
                    mth = (sep > tol * mx) ? "exact" : "lt";
                }
            }

            double[] fi = new double[tset.length];
            double[] Fi = new double[tset.length];
            if ("exact".equalsIgnoreCase(mth)) {
                thm2(v, mu, N, z, onpath, x, Gn1, tset, fi, Fi);
            } else if ("lt".equalsIgnoreCase(mth)) {
                final int[] zf = z;
                final double[] xf = x;
                final double gf = Gn1;
                final double[] muf = mu;
                final int Mf = M;
                final int Nf = N;
                UnaryOperator<Complex> L = new UnaryOperator<Complex>() {
                    @Override
                    public Complex apply(Complex s) {
                        Complex[] y = new Complex[Mf];
                        for (int i = 0; i < Mf; i++) {
                            y[i] = new Complex(xf[i], 0.0);
                        }
                        for (int j : zf) {
                            Complex mj = new Complex(muf[j], 0.0);
                            y[j] = y[j].multiply(mj.divide(s.add(mj)));
                        }
                        Complex acc = buzenComplex(y, Nf - 1)[Nf - 1].divide(gf);
                        for (int j : zf) {
                            Complex mj = new Complex(muf[j], 0.0);
                            acc = acc.multiply(mj.divide(s.add(mj)));
                        }
                        return acc;
                    }
                };
                fi = Laplace_invert.laplace_invert_pdf(L, tset, ltiMethod, 0);
                Fi = Laplace_invert.laplace_invert_cdf(L, tset, ltiMethod, 0);
            } else {
                throw new RuntimeException("pfqn_cyclet_ofree: unknown method '" + method
                        + "', expected auto, exact or lt");
            }

            double[] momi = seriesMoments(v, mu, N, z, onpath, x, Gn1, nmom);
            for (int i = 0; i < tset.length; i++) {
                f[i] += pp[ip] * fi[i];
                F[i] += pp[ip] * Fi[i];
            }
            for (int q = 0; q < nmom; q++) {
                mom[q] += pp[ip] * momi[q];
            }
            methods.add(mth);
        }

        for (int i = 0; i < tset.length; i++) {
            f[i] = Math.max(0.0, f[i]);
            F[i] = Math.min(1.0, Math.max(0.0, F[i]));
        }
        return new PfqnCycletResult(f, F, mom, methods, Math.log(Gn1));
    }

    /**
     * Buzen's convolution: g[k] = G at population k for the node set y.
     *
     * <p>This is the k(y,a,b) recursion of Sec. 7.1 with the node index rolled
     * up: k(y,a,b) = k(y,a-1,b) + y_a k(y,a,b-1), k(y,a,0) = 1, k(y,0,b>0) = 0.
     */
    private static double[] buzen(double[] y, int n) {
        double[] g = new double[n + 1];
        g[0] = 1.0;
        for (int i = 0; i < y.length; i++) {
            for (int k = 1; k <= n; k++) {
                g[k] = g[k] + y[i] * g[k - 1];
            }
        }
        return g;
    }

    private static Complex[] buzenComplex(Complex[] y, int n) {
        Complex[] g = new Complex[n + 1];
        for (int k = 0; k <= n; k++) {
            g[k] = new Complex(0.0, 0.0);
        }
        g[0] = new Complex(1.0, 0.0);
        for (int i = 0; i < y.length; i++) {
            for (int k = 1; k <= n; k++) {
                g[k] = g[k].add(y[i].multiply(g[k - 1]));
            }
        }
        return g;
    }

    /**
     * Theorem 2 in closed form. The density is a finite sum of terms
     * t^k exp(-mu_j t), so its integral is an incomplete gamma and the CDF comes
     * out in closed form too rather than by quadrature.
     */
    private static void thm2(double[] v, double[] mu, int N, int[] z, boolean[] onpath, double[] x,
                             double Gn1, double[] tset, double[] f, double[] F) {
        int M = v.length;
        int m = z.length;
        int nOff = 0;
        for (int i = 0; i < M; i++) {
            if (!onpath[i]) {
                nOff++;
            }
        }
        double[] xoff = new double[nOff];
        int p = 0;
        for (int i = 0; i < M; i++) {
            if (!onpath[i]) {
                xoff[p++] = x[i];
            }
        }
        double[] Gm = buzen(xoff, N - 1);

        double[][] coef = new double[m][N];
        for (int j = 0; j < m; j++) {
            double den = 1.0;
            for (int i = 0; i < m; i++) {
                if (i != j) {
                    den *= (mu[z[i]] - mu[z[j]]);
                }
            }
            if (den == 0.0) {
                throw new RuntimeException("pfqn_cyclet_ofree: Theorem 2 needs distinct service "
                        + "rates on the path; two coincide. Use method 'lt'");
            }
            double[] w = new double[m - 1];
            int q = 0;
            for (int i = 0; i < m; i++) {
                if (i != j) {
                    w[q++] = (v[z[i]] - v[z[j]]) / (mu[z[i]] - mu[z[j]]);
                }
            }
            double[] K = buzen(w, N - 1);
            for (int c = 0; c < N; c++) {
                double Gmc = Gm[N - 1 - c];
                if (Gmc == 0.0) {
                    continue;
                }
                for (int i = 0; i <= c; i++) {
                    coef[j][c - i] += Gmc * K[i] / den;
                }
            }
        }

        double pref = 1.0 / Gn1;
        for (int j = 0; j < m; j++) {
            pref *= mu[z[j]];
        }
        for (int j = 0; j < m; j++) {
            double mj = mu[z[j]];
            double vj = v[z[j]];
            for (int k = 0; k < N; k++) {
                double c = coef[j][k];
                if (c == 0.0) {
                    continue;
                }
                double vk = Math.pow(vj, k);
                double kf = factorial(k);
                for (int it = 0; it < tset.length; it++) {
                    double t = tset[it];
                    if (t < 0.0) {
                        continue;
                    }
                    f[it] += pref * c * vk * Math.pow(t, k) / kf * Math.exp(-mj * t);
                    // int_0^t s^k exp(-mu s) ds = k!/mu^(k+1) P(k+1, mu t)
                    F[it] += pref * c * vk / Math.pow(mj, k + 1) * gammaincInt(k, mj * t);
                }
            }
        }
    }

    /**
     * The regularized lower incomplete gamma P(k+1, x) for INTEGER shape, which
     * is all this class needs: the finite sum is exact for the integer shapes
     * the Erlang terms produce, so no series or continued fraction is required.
     */
    private static double gammaincInt(int k, double x) {
        if (!(x > 0.0)) {
            return 0.0;
        }
        double term = Math.exp(-x);
        double acc = term;
        for (int j = 1; j <= k; j++) {
            term *= x / j;
            acc += term;
        }
        return Math.min(1.0, Math.max(0.0, 1.0 - acc));
    }

    private static double factorial(int k) {
        double r = 1.0;
        for (int i = 2; i <= k; i++) {
            r *= i;
        }
        return r;
    }

    /**
     * The same Buzen convolution run in the ring of truncated power series in s.
     * Every operation in the recursion is an addition or a multiplication, so
     * the series ring carries it unchanged, and E[T^q] = (-1)^q q! [s^q] L(s).
     */
    private static double[] seriesMoments(double[] v, double[] mu, int N, int[] z,
                                          boolean[] onpath, double[] x, double Gn1, int nmom) {
        int M = v.length;
        int K = nmom;
        double[][] Y = new double[M][K + 1];
        for (int i = 0; i < M; i++) {
            if (onpath[i]) {
                double p = 1.0;
                for (int k = 0; k <= K; k++) {
                    Y[i][k] = x[i] * p;
                    p *= (-1.0 / mu[i]);
                }
            } else {
                Y[i][0] = x[i];
            }
        }
        double[][] G = new double[N][K + 1];
        G[0][0] = 1.0;
        for (int i = 0; i < M; i++) {
            for (int nn = 1; nn < N; nn++) {
                double[] add = seriesMul(Y[i], G[nn - 1], K);
                for (int k = 0; k <= K; k++) {
                    G[nn][k] += add[k];
                }
            }
        }
        double[] L = new double[K + 1];
        for (int k = 0; k <= K; k++) {
            L[k] = G[N - 1][k] / Gn1;
        }
        for (int j : z) {
            double[] e = new double[K + 1];
            double p = 1.0;
            for (int k = 0; k <= K; k++) {
                e[k] = p;
                p *= (-1.0 / mu[j]);
            }
            L = seriesMul(L, e, K);
        }
        double[] mom = new double[nmom];
        for (int q = 1; q <= nmom; q++) {
            mom[q - 1] = ((q % 2 == 1) ? -1.0 : 1.0) * factorial(q) * L[q];
        }
        return mom;
    }

    private static double[] seriesMul(double[] a, double[] b, int K) {
        double[] c = new double[K + 1];
        for (int i = 0; i <= K; i++) {
            if (a[i] == 0.0) {
                continue;
            }
            for (int j = 0; i + j <= K; j++) {
                c[i + j] += a[i] * b[j];
            }
        }
        return c;
    }
}
