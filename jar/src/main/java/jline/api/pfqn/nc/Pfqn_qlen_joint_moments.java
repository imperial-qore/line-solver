package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;

import jline.api.moment.Moment_joint_binomial_from_tail;
import jline.api.moment.Moment_joint_central_from_raw;
import jline.api.moment.Moment_joint_cumulant_from_raw;
import jline.api.moment.Moment_joint_factorial_from_binomial;
import jline.api.moment.Moment_joint_raw_from_factorial;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCOptions;
import jline.util.matrix.Matrix;

/**
 * Joint moments of the queue-length vector of a closed product-form network,
 * obtained from normalizing constants.
 *
 * <p>The coordinates are (station,class) pairs. Two pairs sharing a class give
 * the cross-station covariance of that class; two pairs sharing a station give
 * the cross-class covariance at that station, which is what a class-oriented
 * method of moments ({@link Pfqn_comomrm} and its relatives) is positioned to
 * deliver. Two exact routes reach the joint survival array, and both end in the
 * same conversion, the tail edge of the house of moments
 * ({@link jline.api.moment}) followed by the joint central-moment and cumulant
 * conversions:
 *
 * <ul>
 * <li>SINGLE CLASS (R = 1), route "tail". The survival probabilities are ratios
 * of normalizing constants of the network itself,
 * P(n_i &gt;= k_i for all i) = (prod_i L_i^k_i) G(N - sum_i k_i) / G(N), which
 * holds because a load-independent single-class station has the geometric
 * occupancy L_i^n. Only N+1 constants of the ORIGINAL model are needed, which
 * is why any normalizing-constant algorithm serves it.</li>
 * <li>MULTICLASS, route "pmf". The geometric factorization fails, since a
 * multiclass load-independent station carries the multinomial occupancy
 * f_i(n_i) = (|n_i|)! prod_r L_(i,r)^n_(i,r) / n_(i,r)!. What holds instead is
 * the joint law of the selected stations in terms of the COMPLEMENTARY network,
 * the model with those stations deleted and the think times kept,
 * P(n_i = m_i, i in S) = prod_i f_i(m_i) G_(S^c)(N - sum_i m_i) / G(N). The
 * survival array is the reverse cumulative sum of that array, exactly, since
 * the box covers the support.</li>
 * </ul>
 *
 * <p>Neither the factorial nor the raw moments have a one-constant closed form;
 * the survival array is the queue-length functional that does. The
 * normalizing-constant algorithm is INJECTED rather than called at a fixed
 * site: the whole set of populations is known before any evaluation, so it is
 * emitted in one batch and an algorithm that produces several constants in one
 * pass serves it without recomputation. A precomputed table, which the MATLAB
 * and native-Python versions accept directly, is passed here by wrapping it in
 * the same functional interface.
 *
 * <p>Reference:
 * M. Reiser and S. S. Lavenberg. Mean-value analysis of closed multichain
 * queuing networks. Journal of the ACM, 27(2):313-322, 1980.
 *
 * @since LINE 3.0
 */
public final class Pfqn_qlen_joint_moments {
    private Pfqn_qlen_joint_moments() {}

    /**
     * Joint queue-length moments with the default options and no injected
     * source of normalizing constants.
     *
     * @param L service demand matrix of the queueing stations (M x R)
     * @param N population vector (1 x R)
     * @param Z think time vector (1 x R), may be null
     * @param pairs P x 2 matrix of 0-based (station,class) pairs
     * @return the moment arrays over the selected coordinates
     */
    public static Ret.pfqnQlenMoments pfqn_qlen_joint_moments(Matrix L, Matrix N, Matrix Z,
                                                              int[][] pairs) {
        return pfqn_qlen_joint_moments(L, N, Z, pairs, "auto", null, null);
    }

    /**
     * Joint queue-length moments of a closed product-form network.
     *
     * @param L service demand matrix of the QUEUEING stations (M x R). Delay
     *          stations belong in Z, their marginals following a different law.
     *          Load-dependent and multiserver stations are out of scope for
     *          both routes
     * @param N population vector (1 x R)
     * @param Z think time vector (1 x R), may be null for a closed model with
     *          no delay
     * @param pairs P x 2 matrix of 0-based (station,class) pairs, one per
     *              dimension of the returned arrays; null selects every class
     *              of every station
     * @param route "auto" (tail when R = 1, pmf otherwise), "tail" or "pmf"
     * @param lgSource where log G comes from, or null to always call
     *                 {@link Pfqn_nc}. It is invoked ONCE per network with the
     *                 demand matrix and a P x R matrix of populations, and must
     *                 return a P x 1 matrix of log G with NaN where it cannot
     *                 serve; those rows are filled in by Pfqn_nc
     * @param options options passed to Pfqn_nc; null uses the exact method,
     *                since an approximate normalizing constant would silently
     *                make the moments approximate
     * @return the moment arrays over the selected coordinates
     */
    public static Ret.pfqnQlenMoments pfqn_qlen_joint_moments(Matrix L, Matrix N, Matrix Z,
                                                              int[][] pairs, String route,
                                                              BiFunction<Matrix, Matrix, Matrix> lgSource,
                                                              SolverOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (N.length() != R) {
            throw new IllegalArgumentException("pfqn_qlen_joint_moments: N must have one entry "
                    + "per class.");
        }
        Matrix Zv = new Matrix(1, R);
        if (Z != null) {
            for (int r = 0; r < R; r++) {
                Zv.set(0, r, Z.get(r));
            }
        }
        if (pairs == null) {
            pairs = new int[M * R][2];
            int t = 0;
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    pairs[t][0] = i;
                    pairs[t][1] = r;
                    t++;
                }
            }
        }
        int d = pairs.length;
        if (d == 0) {
            throw new IllegalArgumentException("pfqn_qlen_joint_moments: at least one "
                    + "(station,class) pair is required.");
        }
        for (int j = 0; j < d; j++) {
            if (pairs[j].length != 2 || pairs[j][0] < 0 || pairs[j][0] >= M
                    || pairs[j][1] < 0 || pairs[j][1] >= R) {
                throw new IllegalArgumentException("pfqn_qlen_joint_moments: a (station,class) "
                        + "pair is out of range.");
            }
            for (int l = 0; l < j; l++) {
                if (pairs[l][0] == pairs[j][0] && pairs[l][1] == pairs[j][1]) {
                    throw new IllegalArgumentException("pfqn_qlen_joint_moments: the "
                            + "(station,class) pairs must be distinct.");
                }
            }
        }
        if (route == null || "auto".equals(route)) {
            route = (R == 1) ? "tail" : "pmf";
        }
        if (!"tail".equals(route) && !"pmf".equals(route)) {
            throw new IllegalArgumentException("pfqn_qlen_joint_moments: the route must be "
                    + "auto, tail or pmf.");
        }
        if ("tail".equals(route) && R > 1) {
            throw new IllegalArgumentException("pfqn_qlen_joint_moments: the tail route needs "
                    + "the geometric occupancy of a single-class load-independent station; with "
                    + "several classes the multinomial factor breaks the survival identity, so "
                    + "use the pmf route.");
        }
        if (options == null) {
            options = new NCOptions();
            options.method = "exact";
        }

        int[] Nv = new int[R];
        for (int r = 0; r < R; r++) {
            Nv[r] = (int) Math.round(N.get(r));
        }
        int[] dims = new int[d];
        for (int j = 0; j < d; j++) {
            dims[j] = Nv[pairs[j][1]] + 1;
        }

        double[] tail;
        int[] counters = new int[3];   // points, served, evals
        if ("tail".equals(route)) {
            tail = tailRoute(L, Nv, Zv, pairs, dims, lgSource, options, counters);
        } else {
            int[][] outDims = new int[1][];
            tail = pmfRoute(L, Nv, Zv, pairs, lgSource, options, counters, outDims);
            dims = outDims[0];
        }

        double[] b = Moment_joint_binomial_from_tail.moment_joint_binomial_from_tail(tail, dims);
        double[] f = Moment_joint_factorial_from_binomial
                .moment_joint_factorial_from_binomial(b, dims);
        double[] m = Moment_joint_raw_from_factorial.moment_joint_raw_from_factorial(f, dims);
        double[] mc = Moment_joint_central_from_raw.moment_joint_central_from_raw(m, dims);
        double[] kap = Moment_joint_cumulant_from_raw.moment_joint_cumulant_from_raw(m, dims);

        int[] stride = strides(dims);
        Matrix mean = new Matrix(1, d);
        Matrix cov = new Matrix(d, d);
        for (int j = 0; j < d; j++) {
            mean.set(0, j, m[stride[j]]);
            for (int l = 0; l < d; l++) {
                cov.set(j, l, kap[stride[j] + stride[l]]);
            }
        }
        return new Ret.pfqnQlenMoments(tail, b, f, m, mc, kap, dims, mean, cov, route,
                counters[0], counters[1], counters[2]);
    }

    /**
     * Row-major strides of a multi-dimensional array.
     *
     * @param dims extents of the array
     * @return the stride of each dimension
     */
    private static int[] strides(int[] dims) {
        int d = dims.length;
        int[] stride = new int[d];
        stride[d - 1] = 1;
        for (int i = d - 2; i >= 0; i--) {
            stride[i] = stride[i + 1] * dims[i + 1];
        }
        return stride;
    }

    /**
     * Single-class route: the survival array is a ratio of normalizing
     * constants of the original network.
     *
     * @param L demand matrix
     * @param Nv population vector
     * @param Z think time vector
     * @param pairs selected coordinates
     * @param dims extents of the survival array
     * @param lgSource injected source of log G, possibly null
     * @param options options for Pfqn_nc
     * @param counters accumulator for points, served and evals
     * @return the flattened survival array
     */
    private static double[] tailRoute(Matrix L, int[] Nv, Matrix Z, int[][] pairs, int[] dims,
                                      BiFunction<Matrix, Matrix, Matrix> lgSource,
                                      SolverOptions options, int[] counters) {
        int d = dims.length;
        int nel = 1;
        for (int j = 0; j < d; j++) {
            nel *= dims[j];
        }
        List<int[]> need = new ArrayList<int[]>();
        int[] a = new int[d];
        for (int ia = 0; ia < nel; ia++) {
            int s = 0;
            for (int j = 0; j < d; j++) {
                s += a[j];
            }
            if (Nv[0] - s >= 0) {
                need.add(new int[]{Nv[0] - s});
            }
            odometer(a, dims);
        }
        need.add(new int[]{Nv[0]});
        Map<String, Double> index = batchLg(L, dedup(need), Z, lgSource, options, counters);
        double lgN = index.get(key(new int[]{Nv[0]}));

        double[] tail = new double[nel];
        a = new int[d];
        for (int ia = 0; ia < nel; ia++) {
            int s = 0;
            for (int j = 0; j < d; j++) {
                s += a[j];
            }
            if (Nv[0] - s >= 0) {
                double acc = 0.0;
                boolean ok = true;
                for (int j = 0; j < d; j++) {
                    if (a[j] > 0) {
                        double lij = L.get(pairs[j][0], pairs[j][1]);
                        if (lij <= 0) {
                            ok = false;
                            break;
                        }
                        acc += a[j] * Math.log(lij);
                    }
                }
                if (ok) {
                    tail[ia] = Math.exp(acc + index.get(key(new int[]{Nv[0] - s})) - lgN);
                }
            }
            odometer(a, dims);
        }
        return tail;
    }

    /**
     * Multiclass route: the joint law of the selected stations follows from the
     * complementary network, and the survival array is its reverse cumulative
     * sum.
     *
     * @param L demand matrix
     * @param Nv population vector
     * @param Z think time vector
     * @param pairs selected coordinates
     * @param lgSource injected source of log G, possibly null
     * @param options options for Pfqn_nc
     * @param counters accumulator for points, served and evals
     * @param outDims receives the extents of the survival array
     * @return the flattened survival array
     */
    private static double[] pmfRoute(Matrix L, int[] Nv, Matrix Z, int[][] pairs,
                                     BiFunction<Matrix, Matrix, Matrix> lgSource,
                                     SolverOptions options, int[] counters, int[][] outDims) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        int d = pairs.length;
        boolean[] isSel = new boolean[M];
        for (int j = 0; j < d; j++) {
            isSel[pairs[j][0]] = true;
        }
        List<Integer> stations = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (isSel[i]) {
                stations.add(i);
            }
        }
        int ns = stations.size();
        int dc = ns * R;
        int[][] coords = new int[dc][2];
        int t = 0;
        for (int si = 0; si < ns; si++) {
            for (int r = 0; r < R; r++) {
                coords[t][0] = stations.get(si);
                coords[t][1] = r;
                t++;
            }
        }
        int[] cdims = new int[dc];
        int ncel = 1;
        for (int j = 0; j < dc; j++) {
            cdims[j] = Nv[coords[j][1]] + 1;
            ncel *= cdims[j];
        }
        Matrix Lsub = new Matrix(M - ns, R);
        int rowc = 0;
        for (int i = 0; i < M; i++) {
            if (!isSel[i]) {
                for (int r = 0; r < R; r++) {
                    Lsub.set(rowc, r, L.get(i, r));
                }
                rowc++;
            }
        }

        List<int[]> need = new ArrayList<int[]>();
        int[] a = new int[dc];
        for (int ia = 0; ia < ncel; ia++) {
            int[] n = shift(Nv, coords, a, R);
            if (n != null) {
                need.add(n);
            }
            odometer(a, cdims);
        }
        int[][] pops = dedup(need);
        Map<String, Double> indexc;
        if (M - ns == 0) {
            indexc = new LinkedHashMap<String, Double>();
            for (int p = 0; p < pops.length; p++) {
                indexc.put(key(pops[p]), delayLg(Z, pops[p]));
            }
            counters[0] += pops.length;
            counters[1] += pops.length;
        } else {
            indexc = batchLg(Lsub, pops, Z, lgSource, options, counters);
        }
        List<int[]> full = new ArrayList<int[]>();
        full.add(Nv);
        Map<String, Double> indexN = batchLg(L, dedup(full), Z, null, options, counters);
        double lgN = indexN.get(key(Nv));

        int[] dims = new int[d];
        int nel = 1;
        for (int j = 0; j < d; j++) {
            dims[j] = Nv[pairs[j][1]] + 1;
            nel *= dims[j];
        }
        int[] stride = strides(dims);
        double[] marg = new double[nel];
        a = new int[dc];
        for (int ia = 0; ia < ncel; ia++) {
            int[] n = shift(Nv, coords, a, R);
            if (n != null) {
                Double gc = indexc.get(key(n));
                if (gc != null && !Double.isInfinite(gc) && !Double.isNaN(gc)) {
                    double acc = gc - lgN;
                    boolean ok = true;
                    for (int si = 0; si < ns && ok; si++) {
                        int station = stations.get(si);
                        int tot = 0;
                        for (int j = 0; j < dc; j++) {
                            if (coords[j][0] == station) {
                                tot += a[j];
                            }
                        }
                        acc += logFactorial(tot);
                        for (int j = 0; j < dc; j++) {
                            if (coords[j][0] != station || a[j] == 0) {
                                continue;
                            }
                            double lij = L.get(station, coords[j][1]);
                            if (lij <= 0) {
                                ok = false;
                                break;
                            }
                            acc += a[j] * Math.log(lij) - logFactorial(a[j]);
                        }
                    }
                    if (ok) {
                        int pos = 0;
                        for (int j = 0; j < d; j++) {
                            for (int jc = 0; jc < dc; jc++) {
                                if (coords[jc][0] == pairs[j][0] && coords[jc][1] == pairs[j][1]) {
                                    pos += a[jc] * stride[j];
                                }
                            }
                        }
                        marg[pos] += Math.exp(acc);
                    }
                }
            }
            odometer(a, cdims);
        }

        // reverse cumulative sum along every dimension turns the law into the
        // survival array, exactly, because the box covers the support
        double[] tail = marg.clone();
        for (int mode = 0; mode < d; mode++) {
            int n = dims[mode];
            int st = stride[mode];
            int outer = nel / (n * st);
            for (int o = 0; o < outer; o++) {
                for (int s = 0; s < st; s++) {
                    int base = o * n * st + s;
                    for (int i = n - 2; i >= 0; i--) {
                        tail[base + i * st] += tail[base + (i + 1) * st];
                    }
                }
            }
        }
        outDims[0] = dims;
        return tail;
    }

    /**
     * Population left after removing the orders of a multi-index.
     *
     * @param Nv population vector
     * @param coords the (station,class) coordinate of every dimension
     * @param a the multi-index
     * @param R number of classes
     * @return the shifted population, or null if it is infeasible
     */
    private static int[] shift(int[] Nv, int[][] coords, int[] a, int R) {
        int[] n = new int[R];
        System.arraycopy(Nv, 0, n, 0, R);
        for (int j = 0; j < a.length; j++) {
            n[coords[j][1]] -= a[j];
        }
        for (int r = 0; r < R; r++) {
            if (n[r] < 0) {
                return null;
            }
        }
        return n;
    }

    /**
     * Advance a multi-index, last dimension fastest (row-major order).
     *
     * @param a the multi-index, updated in place
     * @param dims extents
     */
    private static void odometer(int[] a, int[] dims) {
        for (int l = a.length - 1; l >= 0; l--) {
            a[l]++;
            if (a[l] < dims[l]) {
                return;
            }
            a[l] = 0;
        }
    }

    /**
     * Deduplicate a list of population vectors, preserving order.
     *
     * @param rows the populations
     * @return the distinct populations
     */
    private static int[][] dedup(List<int[]> rows) {
        Map<String, int[]> seen = new LinkedHashMap<String, int[]>();
        for (int p = 0; p < rows.size(); p++) {
            seen.put(key(rows.get(p)), rows.get(p));
        }
        int[][] out = new int[seen.size()][];
        int t = 0;
        for (Map.Entry<String, int[]> e : seen.entrySet()) {
            out[t++] = e.getValue();
        }
        return out;
    }

    /**
     * Canonical key of a population vector.
     *
     * @param n the population
     * @return the key
     */
    private static String key(int[] n) {
        StringBuilder sb = new StringBuilder();
        for (int r = 0; r < n.length; r++) {
            sb.append(n[r]).append(',');
        }
        return sb.toString();
    }

    /**
     * Log normalizing constant of a pure-delay network,
     * prod_r Z_r^n_r / n_r!.
     *
     * @param Z think time vector
     * @param n population vector
     * @return the logarithm, negative infinity when a class with a positive
     *         population has no think time
     */
    private static double delayLg(Matrix Z, int[] n) {
        double acc = 0.0;
        for (int r = 0; r < n.length; r++) {
            if (n[r] == 0) {
                continue;
            }
            double z = Z.get(r);
            if (z <= 0) {
                return Double.NEGATIVE_INFINITY;
            }
            acc += n[r] * Math.log(z) - logFactorial(n[r]);
        }
        return acc;
    }

    /**
     * Logarithm of a factorial, by summation, exact for the small orders a
     * population box reaches.
     *
     * @param n the argument
     * @return log(n!)
     */
    private static double logFactorial(int n) {
        double acc = 0.0;
        for (int t = 2; t <= n; t++) {
            acc += Math.log(t);
        }
        return acc;
    }

    /**
     * Evaluate log G at a batch of populations, honouring the injected source.
     *
     * @param Lsub demand matrix of the network whose constants are wanted
     * @param pops distinct populations
     * @param Z think time vector
     * @param lgSource injected source, possibly null
     * @param options options for Pfqn_nc
     * @param counters accumulator for points, served and evals
     * @return map from population key to log G
     */
    private static Map<String, Double> batchLg(Matrix Lsub, int[][] pops, Matrix Z,
                                               BiFunction<Matrix, Matrix, Matrix> lgSource,
                                               SolverOptions options, int[] counters) {
        int P = pops.length;
        int R = Z.getNumCols();
        double[] lg = new double[P];
        for (int p = 0; p < P; p++) {
            lg[p] = Double.NaN;
        }
        counters[0] += P;
        if (lgSource != null) {
            Matrix req = new Matrix(P, R);
            for (int p = 0; p < P; p++) {
                for (int r = 0; r < R; r++) {
                    req.set(p, r, pops[p][r]);
                }
            }
            Matrix got = lgSource.apply(Lsub, req);
            if (got == null || got.length() != P) {
                throw new IllegalArgumentException("pfqn_qlen_joint_moments: the lgSource must "
                        + "return one value per requested population.");
            }
            for (int p = 0; p < P; p++) {
                lg[p] = got.get(p);
                if (!Double.isNaN(lg[p])) {
                    counters[1]++;
                }
            }
        }
        Matrix lambda = new Matrix(1, R);
        Map<String, Double> index = new LinkedHashMap<String, Double>();
        for (int p = 0; p < P; p++) {
            if (Double.isNaN(lg[p])) {
                Matrix n = new Matrix(1, R);
                for (int r = 0; r < R; r++) {
                    n.set(0, r, pops[p][r]);
                }
                lg[p] = Pfqn_nc.pfqn_nc(lambda, Lsub, n, Z, options).lG;
                counters[2]++;
            }
            index.put(key(pops[p]), lg[p]);
        }
        return index;
    }
}
