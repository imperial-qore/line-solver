/**
 * Busy Period of a Subnetwork through Normalizing-Constant Point Evaluations
 *
 * @since LINE 3.0
 */
package jline.api.pfqn;

import jline.api.pfqn.nc.Pfqn_clw;
import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Busy period of a subnetwork evaluated from point values of the normalizing
 * constant rather than from the whole population ladder.
 *
 * WHAT THIS BUYS OVER {@link Pfqn_busyp} / {@link Pfqn_busyp_multiclass}. Those
 * walk the whole ladder (the whole lattice, multichain) because the numerator
 * sums over {|m| &gt;= n}. The complement of that set is the SHELLS |m| &lt;= n-1,
 * and summing the product form over the WHOLE lattice is the full network's own
 * normalizing constant, G_I and H convolving to it:
 *
 * <pre>
 *   sum_{m : |m| &gt;= n} G_I(m) H(N-m) = G(N) - sum_{m : |m| &lt;= n-1} G_I(m) H(N-m)
 * </pre>
 *
 * so order n needs only the n lowest shells plus ONE evaluation of G(N). The
 * ordinary busy period n=1 collapses to three constants,
 *
 * <pre>
 *   b(1,I) = [G(N) - H(N)] / sum_r A_r(I) H(N-e_r)
 * </pre>
 *
 * all point evaluations at or near the full population, which is what the
 * normalizing-constant methods are built for. This routine calls CLW
 * (Choudhury-Leung-Whitt, J. ACM 42, 1995, numerical inversion of the generating
 * function); any method returning lG(N) can take its place. The cost stops
 * depending on N: O(shells up to n-1) plus O(n*R) constant evaluations, against
 * O(lattice) for the ladder routines.
 *
 * THE OPEN CASE NEEDS NO INVERSION. The subnetwork's constant sequence has
 * generating function g(z) = prod_{i in I} f_i(z) and the tail is g(1) minus a
 * partial sum, with f_i(1) = 1/(1-rho_i) at a single server and exp(rho_i) at an
 * infinite one. That removes the tail TRUNCATION of the ladder routine, not just
 * its cost: the tail is exact.
 *
 * ACCURACY: the numerator is a difference of two nearly equal quantities when the
 * level set is unlikely, so the relative error grows with n, measured 6e-12 at
 * n=1 against 2.7e-08 at n=N on a three-station closed model at N=20. Cost grows
 * with n too, so the routine is most accurate where it is fastest.
 *
 * SCOPE: CLW's generating function covers single-server and infinite-server
 * stations, so a general load-dependent scaling belongs to
 * {@link Pfqn_busyp_multiclass}. The identity is for the AGGREGATE level set: a
 * per-class one has complement {m_r &lt;= n-1}, the whole lattice in the other
 * chains, which buys nothing.
 */
public final class Pfqn_busyp_clw {
    private Pfqn_busyp_clw() {}

    /**
     * Mean busy period of order n for the subnetwork, via NC point evaluations.
     *
     * @param alpha   relative arrival rates (JxR), one column per chain
     * @param mu      service rates (JxR), the chain-r rate at node j
     * @param P       routing matrices, one per chain (length 1 = shared by all)
     * @param N       population per chain, infinite entries for an open chain
     * @param subnet  zero-based node indexes forming the subnetwork
     * @param n       busy period orders, counting the jobs of every chain
     * @param gamma   external arrival rates (JxR), null for a closed network
     * @param isdelay infinite-server nodes, null meaning all single servers
     * @param method  method name of the normalizing-constant method ("clw")
     * @return mean busy period duration for each requested order
     */
    public static double[] pfqn_busyp_clw(Matrix alpha, Matrix mu, Matrix[] P, double[] N,
                                          int[] subnet, int[] n, Matrix gamma,
                                          boolean[] isdelay, String method) {
        final int J = alpha.getNumRows();
        final int R = alpha.getNumCols();
        boolean isClosed = true, isOpen = true;
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(N[r])) isClosed = false;
            else isOpen = false;
        }
        if (!isClosed && !isOpen)
            throw new IllegalArgumentException(
                    "pfqn_busyp_clw: a mixed model needs the lattice routine pfqn_busyp_multiclass");
        final boolean[] delay = (isdelay == null) ? new boolean[J] : isdelay;
        final String tok = (method == null) ? "clw" : method;

        int[] target = sortedUnique(subnet);
        if (target.length == 0)
            throw new IllegalArgumentException("pfqn_busyp_clw: the subnetwork must be non-empty");
        if (isClosed && target.length >= J)
            throw new IllegalArgumentException(
                    "pfqn_busyp_clw: in a closed network the subnetwork must be a proper subset");
        int[] compl = complement(J, target);

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
                    "pfqn_busyp_clw: no job ever enters the subnetwork, its busy period is undefined");

        int nmax = 1;
        for (int t = 0; t < n.length; t++) nmax = Math.max(nmax, n[t]);

        if (isOpen) {
            // g_I(1) in closed form, so the tail is exact rather than truncated
            double[] rho = new double[target.length];
            double lg1 = 0.0;
            for (int t = 0; t < target.length; t++) {
                for (int r = 0; r < R; r++) rho[t] += L[target[t]][r];
                if (delay[target[t]]) {
                    lg1 += rho[t];
                } else {
                    if (rho[t] >= 1)
                        throw new IllegalArgumentException(
                                "pfqn_busyp_clw: the subnetwork is not stable, its busy period is infinite");
                    lg1 += -Math.log1p(-rho[t]);
                }
            }
            boolean[] delayI = new boolean[target.length];
            for (int t = 0; t < target.length; t++) delayI[t] = delay[target[t]];
            double[] lseq = openCoefficients(rho, delayI, nmax);
            double[] b = new double[n.length];
            for (int t = 0; t < n.length; t++) {
                double[] head = new double[n[t]];
                System.arraycopy(lseq, 0, head, 0, n[t]);
                final double tail = lg1 + Math.log1p(-Math.exp(lse(head) - lg1));
                b[t] = Math.exp(tail - lseq[n[t] - 1] - Math.log(inflow));
            }
            return b;
        }

        int[] pop = new int[R];
        for (int r = 0; r < R; r++) pop[r] = (int) Math.round(N[r]);
        int[] bound = new int[R];
        for (int r = 0; r < R; r++) bound[r] = Math.min(pop[r], Math.max(nmax - 1, 0));
        final int[] stride = new int[R];
        stride[0] = 1;
        for (int r = 1; r < R; r++) stride[r] = stride[r - 1] * (bound[r - 1] + 1);
        int size = 1;
        for (int r = 0; r < R; r++) size *= bound[r] + 1;
        int[][] mvec = new int[size][R];
        for (int idx = 0; idx < size; idx++)
            for (int r = 0; r < R; r++) mvec[idx][r] = (idx / stride[r]) % (bound[r] + 1);

        // the low shells of the subnetwork, the only lattice this routine walks
        double[] lGlow = new double[size];
        for (int idx = 0; idx < size; idx++)
            lGlow[idx] = stationSum(L, delay, target, 0, mvec[idx]);

        final double lGfull = logNc(L, delay, allNodes(J), pop, R, tok);

        double[] b = new double[n.length];
        for (int t = 0; t < n.length; t++) {
            java.util.List<Double> corr = new java.util.ArrayList<Double>();
            java.util.List<Double> den = new java.util.ArrayList<Double>();
            for (int idx = 0; idx < size; idx++) {
                int level = 0;
                for (int r = 0; r < R; r++) level += mvec[idx][r];
                if (level <= n[t] - 1) {
                    // numerator: the full constant minus the shells the level set excludes
                    int[] left = new int[R];
                    for (int r = 0; r < R; r++) left[r] = pop[r] - mvec[idx][r];
                    corr.add(lGlow[idx] + logNc(L, delay, compl, left, R, tok));
                }
                if (level != n[t] - 1) continue;
                // denominator: the flow out of the shell |m| = n-1
                java.util.List<Double> terms = new java.util.ArrayList<Double>();
                for (int r = 0; r < R; r++) {
                    if (A[r] <= 0) continue;
                    int[] left = new int[R];
                    boolean ok = true;
                    for (int s = 0; s < R; s++) {
                        left[s] = pop[s] - mvec[idx][s] - (s == r ? 1 : 0);
                        if (left[s] < 0) ok = false;
                    }
                    if (!ok) continue;
                    terms.add(Math.log(A[r]) + logNc(L, delay, compl, left, R, tok));
                }
                if (!terms.isEmpty()) den.add(lGlow[idx] + lse(toArray(terms)));
            }
            final double num = lGfull + Math.log1p(-Math.exp(lse(toArray(corr)) - lGfull));
            b[t] = Math.exp(num - lse(toArray(den)));
        }
        return b;
    }

    /** Closed-network form with a shared routing matrix and the CLW default. */
    public static double[] pfqn_busyp_clw(Matrix alpha, Matrix mu, Matrix P, double[] N,
                                          int[] subnet, int[] n) {
        return pfqn_busyp_clw(alpha, mu, new Matrix[]{P}, N, subnet, n, null, null, "clw");
    }

    /**
     * log G(k) of a set of nodes: the single servers go to the normalizing-constant
     * method and the infinite servers into its aggregate think time.
     */
    private static double logNc(double[][] L, boolean[] delay, int[] nodes, int[] k, int R,
                                String method) {
        boolean allZero = true;
        for (int r = 0; r < R; r++) {
            if (k[r] < 0) return Double.NEGATIVE_INFINITY;
            if (k[r] != 0) allZero = false;
        }
        if (allZero) return 0.0;
        java.util.List<Integer> queues = new java.util.ArrayList<Integer>();
        double[] Z = new double[R];
        for (int t = 0; t < nodes.length; t++) {
            if (delay[nodes[t]]) {
                for (int r = 0; r < R; r++) Z[r] += L[nodes[t]][r];
            } else {
                queues.add(nodes[t]);
            }
        }
        if (queues.isEmpty()) {
            // only infinite servers left: G(k) = prod_r Z_r^k_r / k_r!
            double out = 0.0;
            for (int r = 0; r < R; r++) {
                if (k[r] == 0) continue;
                if (Z[r] <= 0) return Double.NEGATIVE_INFINITY;
                out += k[r] * Math.log(Z[r]) - factln(k[r]);
            }
            return out;
        }
        if (!"clw".equalsIgnoreCase(method))
            throw new IllegalArgumentException(
                    "pfqn_busyp_clw: only the clw method is wired here; the point evaluation is a "
                    + "plug-in, so another token needs its own call rather than a silent substitution");
        Matrix Lq = new Matrix(queues.size(), R);
        for (int i = 0; i < queues.size(); i++)
            for (int r = 0; r < R; r++) Lq.set(i, r, L[queues.get(i)][r]);
        Matrix Nk = new Matrix(1, R);
        Matrix Zm = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            Nk.set(0, r, k[r]);
            Zm.set(0, r, Z[r]);
        }
        Ret.pfqnNc res = Pfqn_clw.pfqn_clw(Lq, Nk, Zm);
        return res.lG;
    }

    /**
     * log G_I(m) by direct enumeration of the splits of m across the subnetwork
     * nodes, cheap because m is bounded by n-1 and n is small wherever this wins.
     */
    private static double stationSum(double[][] L, boolean[] delay, int[] nodes, int from,
                                     int[] m) {
        final int R = m.length;
        if (from >= nodes.length) {
            for (int r = 0; r < R; r++)
                if (m[r] > 0) return Double.NEGATIVE_INFINITY;
            return 0.0;
        }
        int total = 1;
        for (int r = 0; r < R; r++) total *= m[r] + 1;
        java.util.List<Double> acc = new java.util.ArrayList<Double>();
        int[] head = new int[R];
        for (int idx = 0; idx < total; idx++) {
            int t = idx;
            for (int r = 0; r < R; r++) {
                head[r] = t % (m[r] + 1);
                t /= m[r] + 1;
            }
            final double lterm = nodeTerm(L[nodes[from]], delay[nodes[from]], head);
            if (Double.isInfinite(lterm) && lterm < 0) continue;
            int[] rest = new int[R];
            for (int r = 0; r < R; r++) rest[r] = m[r] - head[r];
            final double lrest = stationSum(L, delay, nodes, from + 1, rest);
            if (Double.isInfinite(lrest) && lrest < 0) continue;
            acc.add(lterm + lrest);
        }
        return lse(toArray(acc));
    }

    /** log X_i(k): multinomial at a single server, 1/prod k_r! at an infinite one. */
    private static double nodeTerm(double[] Li, boolean isdelay, int[] k) {
        int tot = 0;
        for (int r = 0; r < k.length; r++) tot += k[r];
        if (tot == 0) return 0.0;
        double v = isdelay ? 0.0 : factln(tot);
        for (int r = 0; r < k.length; r++) {
            if (k[r] == 0) continue;
            if (Li[r] <= 0) return Double.NEGATIVE_INFINITY;
            v += -factln(k[r]) + k[r] * Math.log(Li[r]);
        }
        return v;
    }

    /** log G_I(0..kmax) of an OPEN subnetwork, convolving the per-node series. */
    private static double[] openCoefficients(double[] rho, boolean[] delay, int kmax) {
        double[] lg = new double[kmax + 1];
        java.util.Arrays.fill(lg, Double.NEGATIVE_INFINITY);
        lg[0] = 0.0;
        for (int i = 0; i < rho.length; i++) {
            double[] li = new double[kmax + 1];
            double acc = 0.0;
            for (int k = 1; k <= kmax; k++) {
                // 1/(1-rho z) has coefficients rho^k; exp(rho z) has rho^k/k!
                acc += Math.log(rho[i]) - (delay[i] ? Math.log(k) : 0.0);
                li[k] = acc;
            }
            double[] lgnew = new double[kmax + 1];
            java.util.Arrays.fill(lgnew, Double.NEGATIVE_INFINITY);
            for (int m = 0; m <= kmax; m++) {
                double[] terms = new double[m + 1];
                for (int k = 0; k <= m; k++) terms[k] = lg[m - k] + li[k];
                lgnew[m] = lse(terms);
            }
            lg = lgnew;
        }
        return lg;
    }

    private static double factln(int k) {
        if (k < 2) return 0.0;
        return org.apache.commons.math3.special.Gamma.logGamma(k + 1.0);
    }

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

    private static int[] allNodes(int J) {
        int[] out = new int[J];
        for (int i = 0; i < J; i++) out[i] = i;
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
