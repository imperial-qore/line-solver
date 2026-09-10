/**
 * Mean Busy Period of Order n for a Subnetwork of a Multichain Network
 *
 * @since LINE 3.0
 */
package jline.api.pfqn;

import jline.util.matrix.Matrix;

/**
 * Multichain generalization of {@link Pfqn_busyp}.
 *
 * H. Daduna, "Busy Periods for Subnetworks in Stochastic Networks: Mean Value
 * Analysis", J. ACM 35(3), 1988, states Theorems 1 and 3 for a single chain and
 * notes in Section 5 that they carry over to the whole product-form class. The
 * proof uses only that the stationary law is product form and that the busy
 * period is Keilson's mean ergodic sojourn time on a level set, neither of which
 * is single-chain, so replacing the scalar population by a per-chain vector m
 * gives, for a closed network,
 *
 * <pre>
 *             sum_{m : |m| &gt;= n}   G_I(m) H(N-m)
 *   b(n,I) = --------------------------------------------------
 *             sum_{m : |m| = n-1}  G_I(m) sum_r A_r(I) H(N-m-e_r)
 * </pre>
 *
 * with G_I and H the normalizing constants of the subnetwork and of its
 * complement at a population VECTOR and A_r(I) the chain-r arrival flow into I.
 * The denominator is the exact chain-r flow across the cut: a chain-r departure
 * from the complement at population k occurs at rate alpha_ir H(k-e_r)/H(k),
 * and the H(k) cancels the state weight. At R=1 the inner sum holds the single
 * term m=n-1 and H(N-m-e_1)=H(N-n), so it collapses to Theorem 1 exactly.
 *
 * THE OPEN CASE NEEDS NO LATTICE. In an open product-form network the stations
 * are independent and the total occupancy of a node depends on the AGGREGATE
 * load sum_r alpha_ir/mu_ir alone, since summing the station function over the
 * compositions of t collapses the multinomial to (sum_r rho_ir)^t. It is
 * therefore reduced here to the single-chain routine on aggregated demands.
 *
 * PER CLASS: with jobclass = r the level set becomes {m_r >= n}, the jobs of
 * chain r alone. Only chain-r arrivals move that level, so the flow sum loses its
 * sum over r and the same two lattices serve every class.
 *
 * A MIXED MODEL keeps the closed lattice with its OPEN dimensions TRUNCATED. The
 * closed chains are conserved between the subnetwork and its complement, the open
 * ones are not: the complement's open count is free, so its open dimensions are
 * summed out and no e_r shift applies to an open chain, removing one job from an
 * unbounded dimension leaving the same sum. The truncation grows until the answer
 * stops moving, and is the only approximation in that branch. For a per-class
 * query on a CLOSED chain there is an exact shortcut: marginalizing the open
 * chains leaves a closed network with the demands deflated by 1/(1-rho_i^open),
 * where the RATES and not the visits must be deflated, since A_r is built from
 * the visit ratios.
 */
public final class Pfqn_busyp_multiclass {
    private Pfqn_busyp_multiclass() {}

    /**
     * Mean busy period of order n for the subnetwork, multichain.
     *
     * @param alpha  relative arrival rates (JxR), one column per chain
     * @param mu     service rates (JxR), the chain-r rate at node j
     * @param P      routing matrices, one per chain (length 1 = shared by all)
     * @param N      population per chain, infinite entries for an open chain
     * @param subnet zero-based node indexes forming the subnetwork
     * @param n      busy period orders, counting the jobs of every chain
     * @param gamma  external arrival rates (JxR), null for a closed network
     * @param phi    dimensionless load-dependent scaling (JxK), null = single server
     * @param tol    relative tolerance of the open-network tail truncation
     * @param jobclass zero-based chain whose own jobs are counted, -1 for every chain
     * @return mean busy period duration for each requested order
     */
    public static double[] pfqn_busyp_multiclass(Matrix alpha, Matrix mu, Matrix[] P,
                                                 double[] N, int[] subnet, int[] n,
                                                 Matrix gamma, Matrix phi, double tol,
                                                 int jobclass) {
        final int J = alpha.getNumRows();
        final int R = alpha.getNumCols();
        if (N.length != R)
            throw new IllegalArgumentException(
                    "pfqn_busyp_multiclass: the population vector must have one entry per chain");
        boolean isClosed = true, isOpen = true;
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(N[r])) {
                isClosed = false;
            } else {
                isOpen = false;
            }
        }
        final boolean isMixed = !isClosed && !isOpen;

        int[] target = Pfqn_busyp_multiclass.sortedUnique(subnet);
        if (target.length == 0)
            throw new IllegalArgumentException("pfqn_busyp_multiclass: the subnetwork must be non-empty");
        if (isClosed && target.length >= J)
            throw new IllegalArgumentException(
                    "pfqn_busyp_multiclass: in a closed network the subnetwork must be a proper subset");
        if (target[target.length - 1] >= J)
            throw new IllegalArgumentException("pfqn_busyp_multiclass: the subnetwork indexes are out of range");
        int[] compl = complement(J, target);

        int totalN = 0;
        for (int r = 0; r < R; r++)
            if (!Double.isInfinite(N[r])) totalN += (int) Math.round(N[r]);
        Matrix scaling = phi;
        if (scaling == null) {
            scaling = new Matrix(J, Math.max(1, totalN));
            for (int i = 0; i < J; i++)
                for (int k = 0; k < scaling.getNumCols(); k++) scaling.set(i, k, 1.0);
        }

        // demands L(i,r) = alpha(i,r)/mu(i,r), zero where chain r does not visit node i
        double[][] L = new double[J][R];
        for (int i = 0; i < J; i++)
            for (int r = 0; r < R; r++)
                L[i][r] = (alpha.get(i, r) > 0 && mu.get(i, r) > 0)
                        ? alpha.get(i, r) / mu.get(i, r) : 0.0;

        // A_r(I): the chain-r rate at which jobs enter the subnetwork from outside it
        double[] A = new double[R];
        double inflow = 0.0;
        for (int r = 0; r < R; r++) {
            Matrix Pr = (P.length == 1) ? P[0] : P[r];
            for (int i = 0; i < compl.length; i++)
                for (int j = 0; j < target.length; j++)
                    A[r] += alpha.get(compl[i], r) * Pr.get(compl[i], target[j]);
            if (gamma != null)
                for (int j = 0; j < target.length; j++) A[r] += gamma.get(target[j], r);
            inflow += A[r];
        }
        if (inflow <= 0)
            throw new IllegalArgumentException(
                    "pfqn_busyp_multiclass: no job ever enters the subnetwork, its busy period is undefined");
        if (jobclass >= R)
            throw new IllegalArgumentException("pfqn_busyp_multiclass: the job class index is out of range");
        if (jobclass >= 0 && A[jobclass] <= 0)
            throw new IllegalArgumentException(
                    "pfqn_busyp_multiclass: no job of that class ever enters the subnetwork");

        if (isOpen && !isMixed) {
            // exact reduction to the single-chain routine on a per-station scalar; the
            // synthetic problem carries no routing, the whole inflow riding on gamma
            Matrix scalar = new Matrix(1, J);
            Matrix gsyn = new Matrix(1, J);
            double[] rho = new double[J];
            for (int i = 0; i < J; i++)
                for (int r = 0; r < R; r++) rho[i] += L[i][r];
            if (jobclass < 0) {
                // the total occupancy depends on the AGGREGATE load alone
                for (int i = 0; i < J; i++) scalar.set(0, i, rho[i]);
                gsyn.set(0, target[0], inflow);
            } else {
                // the class-r marginal is geometric in rho_ir/(1-rho_i+rho_ir), NOT in
                // rho_ir: the other classes inflate the queue the class-r jobs sit in.
                // That collapse assumes a load-INDEPENDENT station.
                for (int t = 0; t < target.length; t++)
                    for (int k = 0; k < scaling.getNumCols(); k++)
                        if (scaling.get(target[t], k) != 1.0)
                            throw new IllegalArgumentException(
                                    "pfqn_busyp_multiclass: a per-class busy period of an open "
                                    + "subnetwork requires load-independent stations");
                for (int i = 0; i < J; i++) {
                    final double den = 1.0 - rho[i] + L[i][jobclass];
                    scalar.set(0, i, den > 0 ? L[i][jobclass] / den : 0.0);
                }
                gsyn.set(0, target[0], A[jobclass]);
            }
            return Pfqn_busyp.pfqn_busyp(scalar, scaling, new Matrix(J, J),
                                         Double.POSITIVE_INFINITY, target, n, gsyn, tol);
        }

        if (!isMixed) {
            int bound = (jobclass < 0) ? totalN : (int) Math.round(N[jobclass]);
            for (int t = 0; t < n.length; t++)
                if (n[t] < 1 || n[t] > bound)
                    throw new IllegalArgumentException(
                            "pfqn_busyp_multiclass: the busy period order must be an integer in "
                            + "1..sum(N), or in 1..N(r) for the busy period of class r alone");
            int[] pop = new int[R];
            boolean[] openChain = new boolean[R];
            for (int r = 0; r < R; r++) pop[r] = (int) Math.round(N[r]);
            return latticeBusyp(L, scaling, target, compl, pop, n, A, jobclass, openChain);
        }

        // A mixed model grows the truncation of its open dimensions until the answer
        // stops moving; that truncation is the only approximation in this branch.
        int nmax = 1;
        for (int t = 0; t < n.length; t++) nmax = Math.max(nmax, n[t]);
        int trunc = 8 + 2 * nmax;
        double[] prev = null;
        boolean[] openChain = new boolean[R];
        for (int r = 0; r < R; r++) openChain[r] = Double.isInfinite(N[r]);
        while (true) {
            int[] pop = new int[R];
            for (int r = 0; r < R; r++)
                pop[r] = openChain[r] ? trunc : (int) Math.round(N[r]);
            double[] b = latticeBusyp(L, scaling, target, compl, pop, n, A, jobclass, openChain);
            if (prev != null) {
                boolean settled = true;
                for (int t = 0; t < b.length; t++)
                    if (Math.abs(b[t] - prev[t]) > 1e-10 * Math.abs(b[t])) settled = false;
                if (settled) return b;
            }
            prev = b;
            trunc *= 2;
            if (trunc > 4096)
                throw new IllegalArgumentException(
                        "pfqn_busyp_multiclass: the mixed busy period did not converge, so some "
                        + "station of the subnetwork is nearly saturated");
        }
    }

    /**
     * The lattice evaluation shared by the closed and the mixed branch. `pop` bounds
     * every chain: the population of a closed one, the truncation of an open one;
     * `openChain` names the dimensions that are NOT conserved, whose complement
     * counts are summed out rather than read at N-m.
     */
    private static double[] latticeBusyp(double[][] L, Matrix scaling, int[] target,
                                         int[] compl, int[] pop, int[] n, double[] A,
                                         int jobclass, boolean[] openChain) {
        final int R = pop.length;
        final int[] stride = new int[R];
        stride[0] = 1;
        for (int r = 1; r < R; r++) stride[r] = stride[r - 1] * (pop[r - 1] + 1);
        int size = 1;
        for (int r = 0; r < R; r++) size *= pop[r] + 1;
        int[][] mvec = new int[size][R];
        for (int idx = 0; idx < size; idx++)
            for (int r = 0; r < R; r++) mvec[idx][r] = (idx / stride[r]) % (pop[r] + 1);

        double[] lG = lgvec(select(L, target), select(scaling, target), mvec, stride, pop);
        double[] lH = lgvec(select(L, compl), select(scaling, compl), mvec, stride, pop);

        // Hbar sums the complement over its unconserved dimensions, so it is indexed
        // by the CLOSED components alone; with no open chain it is lH itself.
        boolean anyOpen = false;
        for (int r = 0; r < R; r++) anyOpen |= openChain[r];
        double[] lHbar = lH;
        if (anyOpen) {
            lHbar = new double[size];
            java.util.Arrays.fill(lHbar, Double.NEGATIVE_INFINITY);
            for (int idx = 0; idx < size; idx++) {
                int j = 0;
                for (int r = 0; r < R; r++)
                    if (!openChain[r]) j += mvec[idx][r] * stride[r];
                lHbar[j] = lse(new double[]{lHbar[j], lH[idx]});
            }
        }

        double[] b = new double[n.length];
        for (int t = 0; t < n.length; t++) {
            final int nt = n[t];
            java.util.List<Double> num = new java.util.ArrayList<Double>();
            java.util.List<Double> den = new java.util.ArrayList<Double>();
            for (int idx = 0; idx < size; idx++) {
                // the level set is |m| for the aggregate busy period and m_r for the
                // class-r one; only chain-r arrivals move m_r
                int level = 0;
                if (jobclass < 0)
                    for (int r = 0; r < R; r++) level += mvec[idx][r];
                else
                    level = mvec[idx][jobclass];
                if (level >= nt) {
                    int j = 0;
                    for (int r = 0; r < R; r++)
                        if (!openChain[r]) j += (pop[r] - mvec[idx][r]) * stride[r];
                    num.add(lG[idx] + lHbar[j]);
                }
                if (level + 1 != nt) continue;
                java.util.List<Double> terms = new java.util.ArrayList<Double>();
                for (int r = 0; r < R; r++) {
                    if (jobclass >= 0 && r != jobclass) continue;
                    if (A[r] <= 0) continue;
                    if (!openChain[r] && pop[r] == mvec[idx][r]) continue;
                    int j = 0;
                    for (int s = 0; s < R; s++) {
                        if (openChain[s]) continue;
                        // a closed chain conserves jobs, so the departing one is removed
                        j += (pop[s] - mvec[idx][s] - (s == r ? 1 : 0)) * stride[s];
                    }
                    terms.add(Math.log(A[r]) + lHbar[j]);
                }
                if (!terms.isEmpty()) den.add(lG[idx] + lse(toArray(terms)));
            }
            b[t] = Math.exp(lse(toArray(num)) - lse(toArray(den)));
        }
        return b;
    }

    /** Closed-network form with the default tolerance and a shared routing matrix. */
    public static double[] pfqn_busyp_multiclass(Matrix alpha, Matrix mu, Matrix P,
                                                 double[] N, int[] subnet, int[] n) {
        return pfqn_busyp_multiclass(alpha, mu, new Matrix[]{P}, N, subnet, n, null, null,
                                     Pfqn_busyp.DEFAULT_TOL, -1);
    }

    /**
     * Log normalizing constants over the whole lattice of a set of nodes.
     *
     * The station function is X_i(m) = multinomial(|m|; m) prod_r L(i,r)^m_r /
     * prod_{k=1}^{|m|} phi_i(k), which at R=1 is the prod_k alpha_i/mu_i(k) of the
     * single-chain routine and at phi(k)=k the infinite-server form. A node whose
     * scaling row is all ones takes the Buzen recursion, O(R) per lattice point;
     * any other node needs the full sub-lattice convolution.
     */
    private static double[] lgvec(double[][] L, double[][] phi, int[][] mvec, int[] stride,
                                  int[] pop) {
        final int nodes = L.length;
        final int R = pop.length;
        final int size = mvec.length;
        double[] lg = new double[size];
        java.util.Arrays.fill(lg, Double.NEGATIVE_INFINITY);
        lg[0] = 0.0;
        int totalN = 0;
        for (int r = 0; r < R; r++) totalN += pop[r];
        for (int i = 0; i < nodes; i++) {
            boolean isLI = true;
            for (int k = 0; k < Math.min(phi[i].length, Math.max(1, totalN)); k++)
                if (phi[i][k] != 1.0) isLI = false;
            if (isLI) {
                double[] lgnew = lg.clone();
                for (int idx = 0; idx < size; idx++) {
                    double acc = lgnew[idx];
                    for (int r = 0; r < R; r++)
                        if (mvec[idx][r] > 0 && L[i][r] > 0)
                            acc = lse(new double[]{acc, Math.log(L[i][r]) + lgnew[idx - stride[r]]});
                    lgnew[idx] = acc;
                }
                lg = lgnew;
            } else {
                double[] lX = station(L[i], phi[i], mvec);
                double[] lgnew = new double[size];
                java.util.Arrays.fill(lgnew, Double.NEGATIVE_INFINITY);
                for (int a = 0; a < size; a++) {
                    if (Double.isInfinite(lg[a]) && lg[a] < 0) continue;
                    for (int c = 0; c < size; c++) {
                        if (Double.isInfinite(lX[c]) && lX[c] < 0) continue;
                        boolean fits = true;
                        int j = 0;
                        for (int r = 0; r < R && fits; r++) {
                            final int s = mvec[a][r] + mvec[c][r];
                            if (s > pop[r]) fits = false;
                            else j += s * stride[r];
                        }
                        if (fits) lgnew[j] = lse(new double[]{lgnew[j], lg[a] + lX[c]});
                    }
                }
                lg = lgnew;
            }
        }
        return lg;
    }

    /** log X_i(m) over the lattice for one node. */
    private static double[] station(double[] Li, double[] phii, int[][] mvec) {
        final int size = mvec.length, R = mvec[0].length;
        double[] out = new double[size];
        for (int idx = 0; idx < size; idx++) {
            int tot = 0;
            for (int r = 0; r < R; r++) tot += mvec[idx][r];
            double v = factln(tot);
            boolean ok = true;
            for (int r = 0; r < R && ok; r++) {
                if (mvec[idx][r] == 0) continue;
                if (Li[r] <= 0) ok = false;
                else v += -factln(mvec[idx][r]) + mvec[idx][r] * Math.log(Li[r]);
            }
            if (!ok) {
                out[idx] = Double.NEGATIVE_INFINITY;
                continue;
            }
            for (int k = 1; k <= tot; k++) v -= Math.log(phii[Math.min(k, phii.length) - 1]);
            out[idx] = v;
        }
        return out;
    }

    /** log(k!) through the log-gamma function. */
    private static double factln(int k) {
        if (k < 2) return 0.0;
        return org.apache.commons.math3.special.Gamma.logGamma(k + 1.0);
    }

    /** log-sum-exp, stable when every entry is -infinity. */
    private static double lse(double[] v) {
        double m = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < v.length; i++)
            if (v[i] > m) m = v[i];
        if (Double.isInfinite(m) || Double.isNaN(m)) return m;
        double s = 0.0;
        for (int i = 0; i < v.length; i++) s += Math.exp(v[i] - m);
        return m + Math.log(s);
    }

    private static double[] toArray(java.util.List<Double> v) {
        double[] out = new double[v.size()];
        for (int i = 0; i < out.length; i++) out[i] = v.get(i);
        return out;
    }

    private static double[][] select(double[][] v, int[] idx) {
        double[][] out = new double[idx.length][];
        for (int i = 0; i < idx.length; i++) out[i] = v[idx[i]];
        return out;
    }

    private static double[][] select(Matrix m, int[] idx) {
        double[][] out = new double[idx.length][m.getNumCols()];
        for (int i = 0; i < idx.length; i++)
            for (int k = 0; k < m.getNumCols(); k++) out[i][k] = m.get(idx[i], k);
        return out;
    }

    private static int[] sortedUnique(int[] v) {
        int[] copy = v.clone();
        java.util.Arrays.sort(copy);
        int m = 0;
        for (int i = 0; i < copy.length; i++)
            if (i == 0 || copy[i] != copy[i - 1]) copy[m++] = copy[i];
        return java.util.Arrays.copyOf(copy, m);
    }

    private static int[] complement(int J, int[] subnet) {
        boolean[] in = new boolean[J];
        for (int i = 0; i < subnet.length; i++) in[subnet[i]] = true;
        int count = 0;
        for (int j = 0; j < J; j++)
            if (!in[j]) count++;
        int[] out = new int[count];
        int t = 0;
        for (int j = 0; j < J; j++)
            if (!in[j]) out[t++] = j;
        return out;
    }
}
