/**
 * @file Cache_spm_size.java
 * @brief Ray (WKB) asymptotic expansion of the cost-capped cache normalizing constant.
 *
 * Port of matlab/src/api/cache/cache_spm_size.m. Approximates what
 * Cache_erec.cache_erec(gamma, m, sigma, k) computes exactly, in the SAME
 * normalization, so the two are interchangeable. This is the item-size extension
 * of Retrieval_rayint, which carries the size-free expansion; call that one when
 * there are no storage costs.
 *
 * Writing E = prod_j m_j! * H, the size-free recursion
 *
 *   E(m,n) = E(m,n-1) + sum_j gamma_{n,j} m_j E(m-1_j,n-1),  E(0,0)=1
 *
 * relaxes to H ~ exp(phi/eps) with n = y/eps, m_j = x_j/eps, whose eikonal
 * e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_j} carries the ray constants
 * xi_j = e^{-phi_j}. With per-item storage costs sigma_i and per-list cost caps
 * k_j the recursion gains the cost coordinate,
 *
 *   E(m,k) = E_i(m,k) + sum_j m_j gamma_ij E_i(m-1_j, k-sigma_i 1_j),
 *
 * so the shift 1_j becomes e_j(y) = (1_j, s(y) 1_j) in the enlarged space
 * X = (x,kappa) and the eikonal picks up the size tilt
 *
 *   e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_{x_j} - s(y) phi_{kappa_j}},
 *
 * with the second family of ray constants zeta_j = e^{-phi_{kappa_j}}. The rays
 * integrate to the discrete saddle point of the product generating function
 *
 *   sum_{m,k} H(m,k) prod_j z_j^{m_j} w_j^{k_j}
 *       = prod_i ( 1 + sum_j gamma_ij z_j w_j^{sigma_i} ),
 *
 * namely, with D_i = 1 + sum_j gamma_ij xi_j zeta_j^{sigma_i} and
 * Psi = sum_i log D_i,
 *
 *   m_j = sum_i gamma_ij xi_j zeta_j^{sigma_i} / D_i,
 *   k_j = sum_i sigma_i gamma_ij xi_j zeta_j^{sigma_i} / D_i,
 *   log H(m,k) ~ Psi - sum_j m_j log xi_j - sum_j k_j log zeta_j
 *                - (d/2) log(2 pi) - (1/2) log det grad^2 Psi,
 *
 * where d is the number of saddle coordinates and, with
 * pi_ij = gamma_ij xi_j zeta_j^{sigma_i} / D_i and
 * Q^i_{jl} = delta_{jl} pi_ij - pi_ij pi_il,
 *
 *   grad^2 Psi = sum_i [1; sigma_i] [1; sigma_i]' (x) Q^i .
 *
 * Setting zeta_j = 1 recovers the size-free expansion exactly.
 *
 * CAPS ARE CUMULATIVE. Cache_erec sums over the states of cost AT MOST k_j, so
 * this class does the same by default (ATMOST). The shadow price
 * eta_j = log zeta_j <= 0 obeys complementary slackness: a list whose
 * unconstrained mean cost already meets its cap is SLACK, keeps zeta_j = 1, and
 * drops out of the saddle, which then degenerates continuously to the size-free
 * expansion; a list whose cap BINDS sits at eta_j < 0, and the states below the
 * boundary decay geometrically with ratio zeta_j, contributing the amplitude
 * factor 1/(1-zeta_j). EXACT gives instead the constant resolving the cost
 * exactly at k_j, which is the raw Laplace formula above with no such factor.
 *
 * SIZE DIVERSITY IS REQUIRED. The Hessian integrand
 * [1;sigma_i][1;sigma_i]' (x) Q^i has rank h, not 2h, so grad^2 Psi is
 * nonsingular only if the sizes actually vary. This is not an artefact: with a
 * single item size the cost of list j is sigma*m_j identically and the cap
 * carries no information. That case is detected and answered exactly rather than
 * passed to a singular saddle. If the sizes share a common divisor the cost
 * lives on a sublattice; the sizes and caps are divided through by their gcd,
 * which is an exact reduction and removes the corresponding lattice factor.
 *
 * OCCUPANCY. Result.pij is the saddle occupancy
 * pi_il = gamma_il xi_l zeta_l^{sigma_i} / D_i and Result.K its per-list cost.
 * These are EXACT-COST quantities: the saddle conditions are sum_i pi_ij = m_j
 * and sum_i sigma_i pi_ij = k_j, so Result.K equals the cap exactly on every
 * binding list. Under cumulative caps the true mean cost is strictly below the
 * cap; use Cache_cost and Cache_prob_erec for that.
 *
 * ACCURACY. The expansion is O(1/n) at fixed occupancy. With a well-separated
 * cap the observed error in log E is around 1e-2 at n = 200 and halves at each
 * doubling of n. It degrades as a binding zeta_j approaches 1, i.e. in the
 * transition between the binding and slack regimes, where the geometric
 * resummation 1/(1-zeta_j) is no longer sharp; Result.zeta and Result.binding
 * report where the saddle sits.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

public final class Cache_spm_size {
    private Cache_spm_size() {}

    /** Sum over states of cost AT MOST k, matching Cache_erec. */
    public static final String ATMOST = "atmost";
    /** Resolve the cost exactly at k. */
    public static final String EXACT = "exact";

    /** Expansion result. */
    public static final class Result {
        /** Normalizing constant, same normalization as Cache_erec (may overflow; use logE). */
        public double E;
        /** Natural logarithm of E, safe for large n. */
        public double logE;
        /** Saddle point xi_j, one entry per list (0 for a list of zero capacity). */
        public double[] xi;
        /** Cost tilt zeta_j on the original size lattice (1 for a slack or absent list). */
        public double[] zeta;
        /** Whether each list's cost cap binds. */
        public boolean[] binding;
        /** Occupancy pi, n x (h+1), column 0 the miss probability. */
        public double[][] pij;
        /** Mean storage cost held by each list. */
        public double[] K;
        /** The exponent Psi - m.log xi - k.log zeta. */
        public double phi;
        /** log det of the Hessian in (log xi, log zeta), restricted to the free coordinates. */
        public double logdetSigma;
        /** gcd of the item sizes, divided out as an exact lattice reduction. */
        public int span;
        /** "spm-size", "spm", "uniform-size", "lattice" or "boundary". */
        public String method;
        /** Size-free error baseline 0.14*(1/min_j m_j + 1/(n - sum_j m_j)); see ACCURACY. */
        public double relerrEst;
        /** Newton iterations used. */
        public int iter;
    }

    /**
     * Cost-capped ray expansion with cumulative caps (matches Cache_erec).
     *
     * @param gamma access factors gamma[i][j], n x h
     * @param m     cache list capacities, length h
     * @param sigma item storage costs, length n, positive integers
     * @param k     per-list storage cost caps, length h
     * @return the expansion result
     */
    public static Result cache_spm_size(double[][] gamma, double[] m, double[] sigma, double[] k) {
        return cache_spm_size(gamma, m, sigma, k, ATMOST);
    }

    /**
     * Cost-capped ray expansion.
     *
     * @param gamma    access factors gamma[i][j], n x h
     * @param m        cache list capacities, length h
     * @param sigma    item storage costs, length n, positive integers
     * @param k        per-list storage cost caps, length h
     * @param costmode ATMOST (matches Cache_erec) or EXACT
     * @return the expansion result
     */
    public static Result cache_spm_size(double[][] gamma, double[] m, double[] sigma,
                                           double[] k, String costmode) {
        if (gamma == null || gamma.length == 0) {
            throw new RuntimeException("cache_spm_size: the access factors must be a non-empty "
                    + "n x h matrix.");
        }
        int n0 = gamma.length;
        int h0 = gamma[0].length;
        if (m == null || m.length != h0) {
            throw new RuntimeException("cache_spm_size: the capacity vector must have one entry "
                    + "per cache list (" + (m == null ? 0 : m.length) + " given, " + h0 + " expected).");
        }
        for (int j = 0; j < h0; j++) {
            if (m[j] < 0 || Math.abs(m[j] - Math.rint(m[j])) > 0) {
                throw new RuntimeException("cache_spm_size: list capacities must be non-negative "
                        + "integers.");
            }
        }
        if (sigma == null || k == null || sigma.length == 0 || k.length == 0) {
            throw new RuntimeException("cache_spm_size: the item sizes and the cost caps are both "
                    + "required. Use Retrieval_rayint for the size-free expansion.");
        }
        if (sigma.length != n0) {
            throw new RuntimeException("cache_spm_size: the item size vector must have one entry "
                    + "per item.");
        }
        if (k.length != h0) {
            throw new RuntimeException("cache_spm_size: the cost cap vector must have one entry "
                    + "per cache list.");
        }
        for (int i = 0; i < n0; i++) {
            if (sigma[i] <= 0 || Math.abs(sigma[i] - Math.rint(sigma[i])) > 0) {
                throw new RuntimeException("cache_spm_size: item sizes must be positive integers.");
            }
        }
        for (int j = 0; j < h0; j++) {
            if (Math.abs(k[j] - Math.rint(k[j])) > 0) {
                throw new RuntimeException("cache_spm_size: storage cost caps must be integers.");
            }
        }
        String mode = (costmode == null) ? ATMOST : costmode.toLowerCase();
        if (!ATMOST.equals(mode) && !EXACT.equals(mode)) {
            throw new RuntimeException("cache_spm_size: the cost mode must be '" + ATMOST + "' or '"
                    + EXACT + "' ('" + costmode + "' given).");
        }
        boolean exact = EXACT.equals(mode);
        boolean capped = true;   // cleared below when a single item size makes the cap uninformative

        Result out = new Result();
        out.xi = new double[h0];
        out.zeta = new double[h0];
        out.binding = new boolean[h0];
        out.K = new double[h0];
        out.span = 1;
        out.phi = Double.NaN;
        out.logdetSigma = Double.NaN;
        out.relerrEst = Double.NaN;
        out.method = "";
        for (int j = 0; j < h0; j++) { out.zeta[j] = 1.0; }
        out.pij = new double[n0][h0 + 1];
        for (int i = 0; i < n0; i++) { out.pij[i][0] = 1.0; }

        // --- boundaries, matching Cache_erec ---
        double msum = 0;
        for (int j = 0; j < h0; j++) { msum += m[j]; }
        boolean negcap = false;
        for (int j = 0; j < h0; j++) { if (k[j] < 0) { negcap = true; } }
        if (msum > n0 || negcap) {
            out.E = 0.0; out.logE = Double.NEGATIVE_INFINITY; out.method = "boundary"; return out;
        }
        if (msum == 0) {
            boolean poscap = false;
            for (int j = 0; j < h0; j++) { if (k[j] > 0) { poscap = true; } }
            out.method = "boundary";
            if (exact && poscap) {
                out.E = 0.0; out.logE = Double.NEGATIVE_INFINITY;
            } else {
                out.E = 1.0; out.logE = 0.0;
            }
            return out;
        }

        // --- items that can never be cached and lists of zero capacity drop out ---
        int n = 0;
        int[] alive = new int[n0];
        for (int i = 0; i < n0; i++) {
            double s = 0;
            for (int j = 0; j < h0; j++) { s += gamma[i][j]; }
            if (s > 0) { alive[n] = i; n++; }
        }
        int hk = 0;
        for (int j = 0; j < h0; j++) { if (m[j] > 0) { hk++; } }
        int[] keep = new int[hk];
        double[] mk = new double[hk];
        double[] kk = new double[hk];
        int t = 0;
        for (int j = 0; j < h0; j++) { if (m[j] > 0) { keep[t] = j; mk[t] = m[j]; kk[t] = k[j]; t++; } }
        double[][] G = new double[n][hk];
        double[] sg = new double[n];
        for (int a = 0; a < n; a++) {
            sg[a] = sigma[alive[a]];
            for (int b = 0; b < hk; b++) { G[a][b] = gamma[alive[a]][keep[b]]; }
        }
        double mksum = 0;
        for (int b = 0; b < hk; b++) { mksum += mk[b]; }
        if (mksum > n) {
            out.E = 0.0; out.logE = Double.NEGATIVE_INFINITY; out.method = "boundary"; return out;
        }
        if (mksum == n) {
            throw new RuntimeException("cache_spm_size: the expansion requires sum(m) < n; at "
                    + "sum(m) = n the saddle point escapes to infinity. Use Cache_erec for a full cache.");
        }

        // --- exact reductions on the cost lattice ---
        long span = 0;
        for (int a = 0; a < n; a++) { span = gcd(span, (long) Math.rint(sg[a])); }
        out.span = (int) span;
        if (exact) {
            for (int b = 0; b < hk; b++) {
                if (Math.rint(kk[b]) % span != 0) {   // unreachable off the sublattice
                    out.E = 0.0; out.logE = Double.NEGATIVE_INFINITY; out.method = "lattice"; return out;
                }
            }
        }
        for (int a = 0; a < n; a++) { sg[a] = sg[a] / span; }
        for (int b = 0; b < hk; b++) { kk[b] = Math.floor(kk[b] / span); }
        // per-list feasibility: the m_j cheapest (dearest) reachable items bound the cost
        for (int b = 0; b < hk; b++) {
            int cnt = 0;
            for (int a = 0; a < n; a++) { if (G[a][b] > 0) { cnt++; } }
            int mj = (int) Math.rint(mk[b]);
            if (cnt < mj) {
                out.E = 0.0; out.logE = Double.NEGATIVE_INFINITY; out.method = "boundary"; return out;
            }
            double[] srt = new double[cnt];
            int c = 0;
            for (int a = 0; a < n; a++) { if (G[a][b] > 0) { srt[c] = sg[a]; c++; } }
            java.util.Arrays.sort(srt);
            double lo = 0;
            for (int a = 0; a < mj; a++) { lo += srt[a]; }
            if (lo > kk[b]) {
                out.E = 0.0; out.logE = Double.NEGATIVE_INFINITY; out.method = "boundary"; return out;
            }
            if (exact) {
                double hi = 0;
                for (int a = cnt - mj; a < cnt; a++) { hi += srt[a]; }
                if (hi < kk[b]) {
                    out.E = 0.0; out.logE = Double.NEGATIVE_INFINITY; out.method = "boundary"; return out;
                }
            }
        }
        // a single item size makes the cost of list j equal to sigma*m_j identically,
        // so the cap carries no information and the 2h saddle is singular (rank h)
        boolean uniform = true;
        for (int a = 1; a < n; a++) { if (sg[a] != sg[0]) { uniform = false; } }
        if (uniform) {
            boolean feasible = true;
            for (int b = 0; b < hk; b++) {
                double cost = sg[0] * mk[b];
                if (exact ? (cost != kk[b]) : (cost > kk[b])) { feasible = false; }
            }
            if (!feasible) {
                out.E = 0.0; out.logE = Double.NEGATIVE_INFINITY; out.method = "uniform-size"; return out;
            }
            capped = false;                     // fall through to the size-free expansion
            out.method = "uniform-size";
        }

        // --- saddle point ---
        double[] th;
        double[] et = new double[hk];
        boolean[] bind = new boolean[hk];
        int[] iters = new int[1];
        double[][] P;
        double[] D;
        double phi;
        double logdet;
        double logH;
        if (capped) {
            double[][] state = new double[2][];
            saddleCost(G, mk, sg, kk, !exact, state, bind, iters);
            th = state[0];
            et = state[1];
            D = new double[n];
            P = occupancy(G, sg, th, et, D);
            int nb = 0;
            for (int b = 0; b < hk; b++) { if (bind[b]) { nb++; } }
            int[] ix = new int[nb];
            int c = 0;
            for (int b = 0; b < hk; b++) { if (bind[b]) { ix[c] = b; c++; } }
            phi = 0;
            for (int a = 0; a < n; a++) { phi += Math.log(D[a]); }
            for (int b = 0; b < hk; b++) { phi -= mk[b] * th[b]; }
            for (int b = 0; b < nb; b++) { phi -= kk[ix[b]] * et[ix[b]]; }
            logdet = logdetChol(hessian(P, sg, ix));
            logH = phi - 0.5 * (hk + nb) * Math.log(2 * Math.PI) - 0.5 * logdet;
            if (!exact) {
                for (int b = 0; b < nb; b++) {   // geometric resummation below the cap
                    logH -= Math.log1p(-Math.exp(et[ix[b]]));
                }
            }
            if (out.method.isEmpty()) { out.method = "spm-size"; }
        } else {
            double[] zeros = new double[n];
            th = saddleFree(G, mk, iters);
            D = new double[n];
            P = occupancy(G, zeros, th, new double[hk], D);
            phi = 0;
            for (int a = 0; a < n; a++) { phi += Math.log(D[a]); }
            for (int b = 0; b < hk; b++) { phi -= mk[b] * th[b]; }
            logdet = logdetChol(hessian(P, zeros, new int[0]));
            logH = phi - 0.5 * hk * Math.log(2 * Math.PI) - 0.5 * logdet;
            if (out.method.isEmpty()) { out.method = "spm"; }
        }

        double logfact = 0;
        for (int j = 0; j < h0; j++) { logfact += logGamma(m[j] + 1.0); }
        out.logE = logH + logfact;              // back to the Cache_erec normalization
        out.E = Math.exp(out.logE);

        // --- ray quantities, reported on the original item and list indexing ---
        for (int b = 0; b < hk; b++) {
            out.xi[keep[b]] = Math.exp(th[b]);
            out.zeta[keep[b]] = Math.exp(et[b] / span);
            out.binding[keep[b]] = bind[b];
        }
        for (int a = 0; a < n; a++) {
            double hit = 0;
            for (int b = 0; b < hk; b++) { out.pij[alive[a]][1 + keep[b]] = P[a][b]; hit += P[a][b]; }
            out.pij[alive[a]][0] = 1.0 - hit;
        }
        for (int j = 0; j < h0; j++) {
            double s = 0;
            for (int i = 0; i < n0; i++) { s += sigma[i] * out.pij[i][1 + j]; }
            out.K[j] = s;
        }
        out.phi = phi;
        out.logdetSigma = logdet;
        out.iter = iters[0];
        double mmin = Double.POSITIVE_INFINITY;
        for (int b = 0; b < hk; b++) { if (mk[b] < mmin) { mmin = mk[b]; } }
        out.relerrEst = 0.14 * (1.0 / mmin + 1.0 / (n - mksum));
        return out;
    }

    // ---------------------------------------------------------------------
    // The saddle machinery. h is the list count, so every linear algebra step
    // below is on a matrix of order at most 2h and a dense solve is right.

    private static long gcd(long a, long b) {
        long x = Math.abs(a);
        long y = Math.abs(b);
        while (y != 0) { long r = x % y; x = y; y = r; }
        return x;
    }

    /**
     * pi_ij = gamma_ij xi_j zeta_j^{sigma_i} / D_i, D_i = 1 + sum_j (that numerator).
     *
     * @param dOut filled with the D_i when non-null, so the caller need not recompute them
     */
    private static double[][] occupancy(double[][] G, double[] sg, double[] th, double[] et,
                                        double[] dOut) {
        int n = G.length;
        int h = th.length;
        double[][] P = new double[n][h];
        for (int a = 0; a < n; a++) {
            double s = 1.0;
            for (int b = 0; b < h; b++) {
                P[a][b] = G[a][b] * Math.exp(th[b] + sg[a] * et[b]);
                s += P[a][b];
            }
            if (dOut != null) { dOut[a] = s; }
            for (int b = 0; b < h; b++) { P[a][b] /= s; }
        }
        return P;
    }

    private static double[] theta0(double[][] G, double[] tgt) {
        int n = G.length;
        int h = tgt.length;
        double tsum = 0;
        for (int b = 0; b < h; b++) { tsum += tgt[b]; }
        double slack = Math.max(1.0 - tsum / n, 1e-9);
        double[] th = new double[h];
        for (int b = 0; b < h; b++) {
            double gb = 0;
            for (int a = 0; a < n; a++) { gb += G[a][b]; }
            th[b] = Math.log(Math.max(tgt[b], 1e-12) / Math.max(gb * slack, 1e-12));
        }
        return th;
    }

    /** The convex dual f = sum_i log D_i - m.theta - k.eta, over the binding eta only. */
    private static double obj(double[][] G, double[] sg, double[] th, double[] et,
                              double[] tgtm, double[] tgtk, boolean[] bind) {
        double[] D = new double[G.length];
        occupancy(G, sg, th, et, D);
        double f = 0;
        for (int a = 0; a < G.length; a++) { f += Math.log(D[a]); }
        for (int b = 0; b < th.length; b++) {
            f -= tgtm[b] * th[b];
            if (bind[b]) { f -= tgtk[b] * et[b]; }
        }
        return f;
    }

    /**
     * grad^2 Psi in (theta, eta), restricted to the free eta coordinates ix. Coordinate
     * u &lt; h is theta_u with weight 1; coordinate h+b is eta_{ix[b]} with weight sigma_i.
     * The block is then sum_i w_u(i) w_v(i) Q^i_{ju,jv} with
     * Q^i_{jl} = delta_{jl} pi_ij - pi_ij pi_il.
     */
    private static double[][] hessian(double[][] P, double[] sg, int[] ix) {
        int n = P.length;
        int h = (n == 0) ? 0 : P[0].length;
        int nb = ix.length;
        int p = h + nb;
        int[] coord = new int[p];
        for (int u = 0; u < h; u++) { coord[u] = u; }
        for (int b = 0; b < nb; b++) { coord[h + b] = ix[b]; }
        double[][] H = new double[p][p];
        double[] w = new double[p];
        for (int a = 0; a < n; a++) {
            for (int u = 0; u < h; u++) { w[u] = 1.0; }
            for (int b = 0; b < nb; b++) { w[h + b] = sg[a]; }
            for (int u = 0; u < p; u++) {
                int ju = coord[u];
                for (int v = 0; v < p; v++) {
                    int jv = coord[v];
                    double q = (ju == jv ? P[a][ju] : 0.0) - P[a][ju] * P[a][jv];
                    H[u][v] += w[u] * w[v] * q;
                }
            }
        }
        return H;
    }

    /** Gaussian elimination with partial pivoting; the order is at most 2h, so tiny. */
    private static double[] solve(double[][] A, double[] b) {
        int p = b.length;
        double[][] M = new double[p][p + 1];
        for (int i = 0; i < p; i++) {
            System.arraycopy(A[i], 0, M[i], 0, p);
            M[i][p] = b[i];
        }
        for (int c = 0; c < p; c++) {
            int piv = c;
            for (int i = c + 1; i < p; i++) { if (Math.abs(M[i][c]) > Math.abs(M[piv][c])) { piv = i; } }
            if (Math.abs(M[piv][c]) < 1e-300) {
                throw new RuntimeException("cache_spm_size: the saddle-point Newton step is not "
                        + "finite. With item sizes this is the rank-h degeneracy of the size-tilted "
                        + "Hessian: the sizes must genuinely vary for the cost coordinate to carry "
                        + "information.");
            }
            double[] tmp = M[c]; M[c] = M[piv]; M[piv] = tmp;
            for (int i = c + 1; i < p; i++) {
                double f = M[i][c] / M[c][c];
                for (int j = c; j <= p; j++) { M[i][j] -= f * M[c][j]; }
            }
        }
        double[] x = new double[p];
        for (int i = p - 1; i >= 0; i--) {
            double s = M[i][p];
            for (int j = i + 1; j < p; j++) { s -= M[i][j] * x[j]; }
            x[i] = s / M[i][i];
        }
        return x;
    }

    /** log det via Cholesky; the Hessian of a strictly convex objective is positive definite. */
    private static double logdetChol(double[][] A) {
        int p = A.length;
        double[][] L = new double[p][p];
        double ld = 0;
        for (int i = 0; i < p; i++) {
            for (int j = 0; j <= i; j++) {
                double s = 0.5 * (A[i][j] + A[j][i]);
                for (int c = 0; c < j; c++) { s -= L[i][c] * L[j][c]; }
                if (i == j) {
                    if (s <= 0) {
                        throw new RuntimeException("cache_spm_size: the saddle-point Hessian is "
                                + "not positive definite; the ray map is singular here. With item "
                                + "sizes this happens when the sizes do not vary over the items the "
                                + "cache can hold, in which case the cost cap carries no information.");
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

    /** Keeps the tilts within a factor e^2 per iteration. */
    private static double damp(double[] d) {
        double dmax = 0;
        for (int i = 0; i < d.length; i++) { dmax = Math.max(dmax, Math.abs(d[i])); }
        double step = 1.0;
        while (step * dmax > 2.0) { step /= 2.0; }
        return step;
    }

    /** Size-free saddle: Newton on theta = log xi for sum_i gamma_ij xi_j / D_i = m_j. */
    private static double[] saddleFree(double[][] G, double[] tgt, int[] iters) {
        int n = G.length;
        int h = tgt.length;
        double[] zeros = new double[n];
        double[] zh = new double[h];
        double[] th = theta0(G, tgt);
        double tmax = 1.0;
        for (int b = 0; b < h; b++) { tmax = Math.max(tmax, Math.abs(tgt[b])); }
        int it = 0;
        for (it = 1; it <= 200; it++) {
            double[][] P = occupancy(G, zeros, th, zh, null);
            double[] g = new double[h];
            double gmax = 0;
            for (int b = 0; b < h; b++) {
                double s = 0;
                for (int a = 0; a < n; a++) { s += P[a][b]; }
                g[b] = s - tgt[b];
                gmax = Math.max(gmax, Math.abs(g[b]));
            }
            if (gmax <= 1e-12 * tmax) { break; }
            double[] d = solve(hessian(P, zeros, new int[0]), g);
            for (int b = 0; b < h; b++) { d[b] = -d[b]; }
            double step = damp(d);
            for (int b = 0; b < h; b++) { th[b] += step * d[b]; }
        }
        iters[0] = it;
        return th;
    }

    /**
     * Cost-constrained saddle. Minimises the convex dual
     * f(theta,eta) = sum_i log D_i - m.theta - k.eta over eta &lt;= 0 when the caps
     * are cumulative, so that complementary slackness selects the binding lists;
     * over all of R^{2h} when the cost is resolved exactly.
     */
    private static void saddleCost(double[][] G, double[] tgtm, double[] sg, double[] tgtk,
                                   boolean cumulative, double[][] state, boolean[] bind, int[] iters) {
        int n = G.length;
        int h = tgtm.length;
        double[] th = theta0(G, tgtm);
        double[] et = new double[h];
        for (int b = 0; b < h; b++) { bind[b] = true; }
        double tol = 1.0;
        for (int b = 0; b < h; b++) {
            tol = Math.max(tol, Math.abs(tgtm[b]));
            tol = Math.max(tol, Math.abs(tgtk[b]));
        }
        tol *= 1e-12;
        int it = 0;
        for (it = 1; it <= 200; it++) {
            double[][] P = occupancy(G, sg, th, et, null);
            double[] gth = new double[h];
            double[] get = new double[h];
            for (int b = 0; b < h; b++) {
                double s1 = 0;
                double s2 = 0;
                for (int a = 0; a < n; a++) { s1 += P[a][b]; s2 += sg[a] * P[a][b]; }
                gth[b] = s1 - tgtm[b];
                get[b] = s2 - tgtk[b];
            }
            if (cumulative) {
                // at eta_j = 0 the cap binds when the mean cost exceeds it
                for (int b = 0; b < h; b++) { bind[b] = (et[b] < 0) || (get[b] > 0); }
            }
            int nb = 0;
            for (int b = 0; b < h; b++) { if (bind[b]) { nb++; } }
            int[] ix = new int[nb];
            int c = 0;
            for (int b = 0; b < h; b++) { if (bind[b]) { ix[c] = b; c++; } }
            double[] g = new double[h + nb];
            double gmax = 0;
            for (int b = 0; b < h; b++) { g[b] = gth[b]; gmax = Math.max(gmax, Math.abs(g[b])); }
            for (int b = 0; b < nb; b++) {
                g[h + b] = get[ix[b]];
                gmax = Math.max(gmax, Math.abs(g[h + b]));
            }
            if (gmax <= tol) { break; }
            double[] d = solve(hessian(P, sg, ix), g);
            for (int b = 0; b < d.length; b++) { d[b] = -d[b]; }
            double step = damp(d);
            double fcur = obj(G, sg, th, et, tgtm, tgtk, bind);
            // Backtrack until the dual decreases. The slack is essential, not cosmetic:
            // Newton reaches the floating-point floor of f in a handful of steps, and a
            // strict test then rejects every step and halves to zero without converging.
            double ftol = 1e-12 * (1.0 + Math.abs(fcur));
            double[] thn = new double[h];
            double[] etn = new double[h];
            for (int ls = 0; ls < 40; ls++) {
                for (int b = 0; b < h; b++) { thn[b] = th[b] + step * d[b]; etn[b] = et[b]; }
                for (int b = 0; b < nb; b++) { etn[ix[b]] = et[ix[b]] + step * d[h + b]; }
                if (cumulative) {
                    for (int b = 0; b < h; b++) { if (etn[b] > 0) { etn[b] = 0.0; } }
                }
                if (obj(G, sg, thn, etn, tgtm, tgtk, bind) <= fcur + ftol) { break; }
                step /= 2.0;
            }
            double moved = 0;
            for (int b = 0; b < h; b++) {
                moved = Math.max(moved, Math.abs(thn[b] - th[b]));
                moved = Math.max(moved, Math.abs(etn[b] - et[b]));
            }
            System.arraycopy(thn, 0, th, 0, h);
            System.arraycopy(etn, 0, et, 0, h);
            if (moved <= 1e-13) { break; }   // the iterate can no longer move: at the floor
        }
        iters[0] = it;
        if (cumulative) {
            for (int b = 0; b < h; b++) { bind[b] = et[b] < 0; }
        }
        state[0] = th;
        state[1] = et;
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
