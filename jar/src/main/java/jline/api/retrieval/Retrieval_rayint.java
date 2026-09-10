/**
 * @file Retrieval_rayint.java
 * @brief Ray (WKB) asymptotic expansion of the list-based cache normalizing constant.
 *
 * Port of matlab/src/api/retrieval/retrieval_rayint.m. Approximates the constant
 * that Cache_erec computes exactly, in the SAME normalization, so the two are
 * interchangeable:
 *
 *   E(m,n) = E(m,n-1) + sum_j gamma_{n,j} m_j E(m-1_j,n-1),   E(0,0)=1
 *
 * Writing E = prod_j m_j! * Et, the relaxation Et ~ H exp(phi/eps) with n = y/eps
 * and m_j = x_j/eps gives the eikonal e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_j},
 * whose rays carry the constants xi_j = e^{-phi_j}. With S(v) = 1 + sum_j
 * gamma_j(v) xi_j,
 *
 *   x_j    = int_0^y gamma_j(v) xi_j / S(v) dv          (the saddle conditions)
 *   phi    = int_0^y log S(v) dv - sum_j x_j log xi_j
 *   H      = (2 pi)^{-h/2} sqrt(S(y)/S(0)) / sqrt(prod_j xi_j * det A)
 *   A_{ik} = d x_i / d xi_k
 *
 * and E ~ prod_j m_j! * eps^{h/2} H exp(phi/eps).
 *
 * DISCRETE (gamma as an n x h matrix). The ray integrals are the sums they
 * discretize and the expansion collapses to the Laplace form
 *
 *   Et ~ (2 pi)^{-h/2} exp(sum_k log D_k - sum_j m_j log xi_j) / sqrt(det Sigma)
 *
 * with D_k = 1 + sum_j gamma_{k,j} xi_j, sum_k gamma_{k,j} xi_j / D_k = m_j and
 * Sigma = A * diag(xi) the Hessian in log xi. This is the more accurate of the
 * two forms; the sqrt(S(y)/S(0)) factor is exactly the Euler-Maclaurin term
 * relating sum_k to int dv and is already accounted for.
 *
 * CONTINUUM (gamma as a Profile on v in [0,1]). Composite Simpson quadrature on
 * the profile itself, the form written in the note. Costs roughly a factor two
 * in accuracy but does not need the n rows.
 *
 * ACCURACY. The relative error is O(1/n) at fixed occupancy but is governed by
 * the smallest occupancy rather than by n, tracking
 * 0.14 * (1/min_j m_j + 1/(n - sum_j m_j)), so a per cent needs every m_j and
 * n - sum_j m_j above about 15 and a part in a thousand needs them above about
 * 150. Returned in Result.relerrEst. Lists with m_j = 0 contribute nothing and
 * are dropped before the saddle is solved.
 *
 * This is the no-fetch (q=0) case, i.e. the same quantity as Cache_erec. The
 * delayed-hit extension carrying the fetch coordinates is NOT implemented: its
 * eikonal is known but its amplitude has not been derived.
 *
 * @since LINE 3.0
 */
package jline.api.retrieval;

public final class Retrieval_rayint {
    private Retrieval_rayint() {}

    /** Access-factor profile gamma_j(v) for v in [0,1]; returns a v.length x h matrix. */
    public interface Profile {
        double[][] at(double[] v);
    }

    /** Expansion result. */
    public static final class Result {
        /** Normalizing constant, same normalization as Cache_erec (may overflow; use logE). */
        public double E;
        /** Natural logarithm of E, safe for large n. */
        public double logE;
        /** Saddle point xi_j, one entry per list (0 for a list of zero capacity). */
        public double[] xi;
        /** log det of the Hessian in log xi. */
        public double logdetSigma;
        /** Estimated relative error, 0.14*(1/min_j m_j + 1/(n - sum_j m_j)). */
        public double relerrEst;
        /** "saddle", "rayint" or "boundary". */
        public String method;
        /** Newton iterations used. */
        public int iter;
    }

    /**
     * Discrete (saddle) form.
     *
     * @param gamma access factors gamma[k][j], n x h
     * @param m     cache list capacities, length h
     * @return the expansion result
     */
    public static Result retrieval_rayint(double[][] gamma, double[] m) {
        if (gamma == null || gamma.length == 0) {
            throw new RuntimeException("retrieval_rayint: the access factors must be a non-empty n x h matrix.");
        }
        return run(gamma, null, gamma.length, m, 0);
    }

    /**
     * Continuum (ray-integral) form.
     *
     * @param gfun  access-factor profile on v in [0,1]
     * @param m     cache list capacities, length h
     * @param n     number of items
     * @return the expansion result
     */
    public static Result retrieval_rayint(Profile gfun, double[] m, int n) {
        return retrieval_rayint(gfun, m, n, 4097);
    }

    /**
     * Continuum (ray-integral) form with an explicit quadrature resolution.
     *
     * @param gfun  access-factor profile on v in [0,1]
     * @param m     cache list capacities, length h
     * @param n     number of items
     * @param nquad composite Simpson nodes (forced odd, at least 5)
     * @return the expansion result
     */
    public static Result retrieval_rayint(Profile gfun, double[] m, int n, int nquad) {
        if (gfun == null) {
            throw new RuntimeException("retrieval_rayint: the access-factor profile must not be null.");
        }
        int nq = Math.max(5, 2 * (nquad / 2) + 1);
        return run(null, gfun, n, m, nq);
    }

    // ---------------------------------------------------------------------

    private static Result run(double[][] gammaMat, Profile gfun, int n, double[] m, int nquad) {
        if (m == null || m.length == 0) {
            throw new RuntimeException("retrieval_rayint: the capacity vector must not be empty.");
        }
        int h = m.length;
        double msum = 0;
        for (int j = 0; j < h; j++) {
            if (m[j] < 0 || Math.abs(m[j] - Math.rint(m[j])) > 0) {
                throw new RuntimeException("retrieval_rayint: list capacities must be non-negative integers.");
            }
            msum += m[j];
        }
        if (gammaMat != null && gammaMat[0].length != h) {
            throw new RuntimeException("retrieval_rayint: the capacity vector must have one entry per cache list ("
                    + h + " given, " + gammaMat[0].length + " expected).");
        }

        Result out = new Result();
        out.xi = new double[h];
        out.logdetSigma = Double.NaN;
        out.relerrEst = Double.NaN;

        if (msum > n) {
            out.E = 0.0; out.logE = Double.NEGATIVE_INFINITY; out.method = "boundary"; return out;
        }
        if (msum == 0) {
            out.E = 1.0; out.logE = 0.0; out.method = "boundary"; return out;
        }
        if (msum == n) {
            throw new RuntimeException("retrieval_rayint: the expansion requires sum(m) < n; at sum(m) = n the "
                    + "saddle point escapes to infinity. Use Cache_erec for a full cache.");
        }

        // lists of zero capacity contribute nothing and would make the saddle singular
        int hk = 0;
        for (int j = 0; j < h; j++) { if (m[j] > 0) hk++; }
        int[] keep = new int[hk];
        double[] mk = new double[hk];
        int t = 0;
        for (int j = 0; j < h; j++) { if (m[j] > 0) { keep[t] = j; mk[t] = m[j]; t++; } }

        double[][] G;
        double[] w;
        double[] tgt;
        if (gammaMat != null) {
            G = new double[n][hk];
            for (int k = 0; k < n; k++) {
                for (int a = 0; a < hk; a++) { G[k][a] = gammaMat[k][keep[a]]; }
            }
            w = new double[n];
            for (int k = 0; k < n; k++) { w[k] = 1.0; }
            tgt = mk;
        } else {
            double[] v = new double[nquad];
            for (int k = 0; k < nquad; k++) { v[k] = ((double) k) / (nquad - 1); }
            double[][] Gfull = gfun.at(v);
            if (Gfull.length != nquad) {
                throw new RuntimeException("retrieval_rayint: the access-factor profile must return one row per "
                        + "evaluation point.");
            }
            if (Gfull[0].length != h) {
                throw new RuntimeException("retrieval_rayint: the profile must return one column per cache list ("
                        + h + " expected, " + Gfull[0].length + " given).");
            }
            G = new double[nquad][hk];
            for (int k = 0; k < nquad; k++) {
                for (int a = 0; a < hk; a++) { G[k][a] = Gfull[k][keep[a]]; }
            }
            w = new double[nquad];
            double step = 1.0 / (3.0 * (nquad - 1));
            for (int k = 0; k < nquad; k++) {
                double c;
                if (k == 0 || k == nquad - 1) { c = 1.0; }
                else if (k % 2 == 1) { c = 4.0; }
                else { c = 2.0; }
                w[k] = c * step;
            }
            tgt = new double[hk];
            for (int a = 0; a < hk; a++) { tgt[a] = mk[a] / n; }
        }

        int[] iters = new int[1];
        double[] xi = saddle(G, w, tgt, iters);

        int K = G.length;
        double[] S = new double[K];
        for (int k = 0; k < K; k++) {
            double s = 1.0;
            for (int a = 0; a < hk; a++) { s += G[k][a] * xi[a]; }
            S[k] = s;
        }
        double[][] Sig = hessian(G, w, xi);
        double logdet = logdetChol(Sig);

        double logEt;
        if (gammaMat != null) {
            double phi = 0;
            for (int k = 0; k < K; k++) { phi += Math.log(S[k]); }
            for (int a = 0; a < hk; a++) { phi -= mk[a] * Math.log(xi[a]); }
            logEt = -0.5 * hk * Math.log(2 * Math.PI) + phi - 0.5 * logdet;
            out.method = "saddle";
        } else {
            double phi = 0;
            for (int k = 0; k < K; k++) { phi += w[k] * Math.log(S[k]); }
            for (int a = 0; a < hk; a++) { phi -= (mk[a] / n) * Math.log(xi[a]); }
            double s0 = S[0];
            double sy = S[K - 1];
            logEt = -0.5 * hk * Math.log(n) - 0.5 * hk * Math.log(2 * Math.PI) + n * phi
                    - 0.5 * logdet + 0.5 * Math.log(sy / s0);
            out.method = "rayint";
        }

        double logfact = 0;
        for (int j = 0; j < h; j++) { logfact += logGamma(m[j] + 1.0); }
        out.logE = logEt + logfact;
        out.E = Math.exp(out.logE);
        for (int a = 0; a < hk; a++) { out.xi[keep[a]] = xi[a]; }
        out.logdetSigma = logdet;
        out.iter = iters[0];
        double mmin = Double.POSITIVE_INFINITY;
        for (int a = 0; a < hk; a++) { if (mk[a] < mmin) mmin = mk[a]; }
        out.relerrEst = 0.14 * (1.0 / mmin + 1.0 / (n - msum));
        return out;
    }

    /**
     * Newton on theta = log xi for sum_k w_k gamma_{k,j} xi_j / (1 + sum_l gamma_{k,l} xi_l) = tgt_j.
     * The objective sum_k w_k log(1 + sum_l gamma_{k,l} e^{theta_l}) - tgt'*theta is strictly
     * convex, so the root is unique and damped Newton converges globally.
     */
    private static double[] saddle(double[][] G, double[] w, double[] tgt, int[] iters) {
        int K = G.length;
        int h = tgt.length;
        double wsum = 0;
        for (int k = 0; k < K; k++) { wsum += w[k]; }
        double tsum = 0;
        for (int a = 0; a < h; a++) { tsum += tgt[a]; }
        double slack = Math.max(1.0 - tsum / wsum, 1e-9);
        double[] th = new double[h];
        for (int a = 0; a < h; a++) {
            double gb = 0;
            for (int k = 0; k < K; k++) { gb += w[k] * G[k][a]; }
            th[a] = Math.log(Math.max(tgt[a], 1e-12) / Math.max(gb * slack, 1e-12));
        }
        double[] xi = new double[h];
        double tmax = 0;
        for (int a = 0; a < h; a++) { tmax = Math.max(tmax, Math.abs(tgt[a])); }
        int it = 0;
        for (it = 1; it <= 200; it++) {
            for (int a = 0; a < h; a++) { xi[a] = Math.exp(th[a]); }
            double[] g = new double[h];
            for (int k = 0; k < K; k++) {
                double s = 1.0;
                for (int a = 0; a < h; a++) { s += G[k][a] * xi[a]; }
                for (int a = 0; a < h; a++) { g[a] += w[k] * G[k][a] * xi[a] / s; }
            }
            double gmax = 0;
            for (int a = 0; a < h; a++) { g[a] -= tgt[a]; gmax = Math.max(gmax, Math.abs(g[a])); }
            if (gmax <= 1e-12 * Math.max(1.0, tmax)) { break; }
            double[][] H = hessian(G, w, xi);
            double[] d = solve(H, g);
            double dmax = 0;
            for (int a = 0; a < h; a++) { d[a] = -d[a]; dmax = Math.max(dmax, Math.abs(d[a])); }
            double step = 1.0;
            while (step * dmax > 2.0) { step /= 2.0; }
            for (int a = 0; a < h; a++) { th[a] += step * d[a]; }
        }
        iters[0] = it;
        for (int a = 0; a < h; a++) { xi[a] = Math.exp(th[a]); }
        return xi;
    }

    /** Hessian in theta = log xi; equals A*diag(xi) with A_{ik} = d x_i / d xi_k. */
    private static double[][] hessian(double[][] G, double[] w, double[] xi) {
        int K = G.length;
        int h = xi.length;
        double[][] H = new double[h][h];
        for (int k = 0; k < K; k++) {
            double s = 1.0;
            for (int a = 0; a < h; a++) { s += G[k][a] * xi[a]; }
            for (int a = 0; a < h; a++) {
                double aa = G[k][a] * xi[a] / s;
                H[a][a] += w[k] * aa;
                for (int b = 0; b < h; b++) {
                    double ab = G[k][b] * xi[b] / s;
                    H[a][b] -= w[k] * aa * ab;
                }
            }
        }
        return H;
    }

    /** Gaussian elimination with partial pivoting; h is the number of cache lists, so tiny. */
    private static double[] solve(double[][] A, double[] b) {
        int h = b.length;
        double[][] M = new double[h][h + 1];
        for (int i = 0; i < h; i++) {
            System.arraycopy(A[i], 0, M[i], 0, h);
            M[i][h] = b[i];
        }
        for (int c = 0; c < h; c++) {
            int piv = c;
            for (int i = c + 1; i < h; i++) { if (Math.abs(M[i][c]) > Math.abs(M[piv][c])) { piv = i; } }
            if (Math.abs(M[piv][c]) < 1e-300) {
                throw new RuntimeException("retrieval_rayint: the saddle-point Hessian is singular; the ray map "
                        + "is degenerate here.");
            }
            double[] tmp = M[c]; M[c] = M[piv]; M[piv] = tmp;
            for (int i = c + 1; i < h; i++) {
                double f = M[i][c] / M[c][c];
                for (int j = c; j <= h; j++) { M[i][j] -= f * M[c][j]; }
            }
        }
        double[] x = new double[h];
        for (int i = h - 1; i >= 0; i--) {
            double s = M[i][h];
            for (int j = i + 1; j < h; j++) { s -= M[i][j] * x[j]; }
            x[i] = s / M[i][i];
        }
        return x;
    }

    /** log det via Cholesky; the Hessian of a strictly convex objective is positive definite. */
    private static double logdetChol(double[][] A) {
        int h = A.length;
        double[][] L = new double[h][h];
        double ld = 0;
        for (int i = 0; i < h; i++) {
            for (int j = 0; j <= i; j++) {
                double s = 0.5 * (A[i][j] + A[j][i]);
                for (int k = 0; k < j; k++) { s -= L[i][k] * L[j][k]; }
                if (i == j) {
                    if (s <= 0) {
                        throw new RuntimeException("retrieval_rayint: the saddle-point Hessian is not positive "
                                + "definite; the ray map is singular here.");
                    }
                    L[i][j] = Math.sqrt(s);
                    ld += 2.0 * Math.log(L[i][j]);
                } else {
                    L[i][j] = s / L[j][j];
                }
            }
        }
        return ld;
    }

    /** Lanczos log-Gamma, so that the m_j! normalization does not overflow. */
    private static double logGamma(double x) {
        double[] c = {676.5203681218851, -1259.1392167224028, 771.32342877765313,
                -176.61502916214059, 12.507343278686905, -0.13857109526572012,
                9.9843695780195716e-6, 1.5056327351493116e-7};
        if (x < 0.5) {
            return Math.log(Math.PI / Math.sin(Math.PI * x)) - logGamma(1 - x);
        }
        double z = x - 1;
        double a = 0.99999999999980993;
        double tt = z + 7.5;
        for (int i = 0; i < c.length; i++) { a += c[i] / (z + i + 1); }
        return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(tt) - tt + Math.log(a);
    }
}
