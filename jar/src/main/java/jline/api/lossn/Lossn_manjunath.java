/**
 * @file Exact analysis of loss networks by the Manjunath-Sikdar transform
 *
 * @since LINE 3.0
 */
package jline.api.lossn;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Exact normalization constant, carried load and blocking of a loss network by
 * the transform technique of Manjunath and Sikdar.
 *
 * Calls on route r arrive Poisson at rate nu_r with unit mean holding time, so
 * nu_r is the offered load, and a call is admitted only while every constraint
 * holds, sum_r A(j,r) n_r &lt;= C(j). The admissible set is coordinate convex, so
 * Kelly's truncation theorem gives the truncated product form
 * p(n) = nu^n / n! / g(C) and every metric is a ratio of normalization
 * constants:
 *
 * <pre>
 *     g(C)   = sum_{A n &lt;= C} prod_r nu_r^{n_r} / n_r!
 *     E[n_r] = nu_r g(C - A e_r) / g(C)
 *     Loss_r = 1 - g(C - A e_r) / g(C)
 * </pre>
 *
 * because a class r call is blocked exactly when the state cannot absorb one
 * more unit of its own requirement vector.
 *
 * WHY IT IS A COEFFICIENT COMPUTATION AND NOT A QUADRATURE. Writing each
 * indicator as a contour integral turns g(C) into a J-fold integral over the
 * unit circle whose integrand factorizes into the per-route z-transforms. Inside
 * the circle the only pole in z_j sits at the origin with order C_j+1, so each
 * integration is a residue, i.e. a Taylor coefficient. The routine therefore
 * never evaluates an integral: it builds the generating function as a
 * multivariate power series truncated at degree C_j in z_j, one
 * shift-and-accumulate convolution per route, and discharges each '&lt;='
 * constraint by summing the coefficients of degrees 0..C_j along that dimension.
 * Truncation is exact because A is nonnegative, so a monomial above degree C_j
 * can never contribute to an extracted coefficient.
 *
 * THE ELIMINATION ORDER IS THE MEMORY BOUND. Contour integrations are
 * interleaved with the product rather than deferred: variable z_j is created
 * when the first route with A(j,r) != 0 is multiplied in and integrated out
 * immediately after the last one. Peak memory is therefore the product of
 * (C_j+1) over the SIMULTANEOUSLY LIVE links, an induced width of the route-link
 * incidence, not over all J links. That product is bounded by
 * DEFAULT_MAX_LIVE_STATES and a region above it is refused by name rather than
 * allowed to exhaust the heap: the algorithm is exact but not unconditionally
 * cheap, and Lossn_mci answers the same question at any size.
 *
 * Unlike Lossn_erlangfp this is exact rather than a reduced-load approximation,
 * and unlike Lossn_mci it carries no sampling error, which is what matters for
 * rare blocking: a loss probability of 1e-4 recovered from a sampled throughput
 * is dominated by the estimator variance.
 *
 * Reference: D. Manjunath and B. Sikdar, Integral Expressions for the Numerical
 * Evaluation of Product Form Expressions Over Irregular Multidimensional Integer
 * Spaces.
 */
public final class Lossn_manjunath {
    private Lossn_manjunath() {}

    /**
     * Cap on the product of (C_j+1) over the simultaneously live links, i.e. on
     * the number of series coefficients held at once. 2^26 coefficients is half a
     * gigabyte of double, which no region a FiniteCapacityRegion can express
     * reaches by accident. Raise it deliberately.
     */
    public static final long DEFAULT_MAX_LIVE_STATES = 1L << 26;

    /**
     * The reduced admission rule: rows that constrain nothing dropped, every
     * remaining row divided by the gcd of its entries and its right-hand side.
     */
    private static final class Rule {
        long[][] A;
        long[] C;
    }

    /** gcd on nonnegative long, with gcd(0, a) = a as MATLAB's gcd has it. */
    private static long gcd(long a, long b) {
        while (b != 0) {
            long t = a % b;
            a = b;
            b = t;
        }
        return a < 0 ? -a : a;
    }

    /**
     * Integer view of (A, C), refusing a fractional or negative entry by name.
     *
     * Both the dropping and the gcd reduction are exact. The second is what makes
     * the cost tractable when the class sizes share a factor: a row (4, 8) &lt;= 20
     * becomes (1, 2) &lt;= 5 and its dimension shrinks from 21 coefficients to 6.
     */
    private static Rule integralize(double[][] A, double[] C) {
        int J = C.length;
        int R = (J > 0) ? A[0].length : 0;
        long[][] rows = new long[J][];
        long[] rhs = new long[J];
        int kept = 0;
        for (int j = 0; j < J; j++) {
            double c = C[j];
            if (c < 0.0 || Math.abs(c - Math.rint(c)) > 1e-9) {
                throw new RuntimeException("lossn_manjunath: C must contain nonnegative integers -- the "
                        + "residue argument counts whole units of capacity. Use lossn_mci, which "
                        + "compares in real arithmetic, or lossn_erlangfp");
            }
            long[] row = new long[R];
            boolean any = false;
            for (int r = 0; r < R; r++) {
                double a = A[j][r];
                if (a < 0.0 || Math.abs(a - Math.rint(a)) > 1e-9) {
                    throw new RuntimeException("lossn_manjunath: A must contain nonnegative integers -- "
                            + "the residue argument counts whole units of capacity. Use lossn_mci, "
                            + "which compares in real arithmetic, or lossn_erlangfp");
                }
                row[r] = Math.round(a);
                if (row[r] != 0) {
                    any = true;
                }
            }
            // A row of zeros bounds nothing and is dropped, not carried as a
            // one-coefficient dimension: keeping it would leave first/last
            // undefined for that link.
            if (!any) {
                continue;
            }
            long g = Math.round(c);
            for (int r = 0; r < R; r++) {
                if (row[r] != 0) {
                    g = gcd(g, row[r]);
                }
            }
            if (g > 1) {
                for (int r = 0; r < R; r++) {
                    row[r] /= g;
                }
                rhs[kept] = Math.round(c) / g;   // floor, both nonnegative
            } else {
                rhs[kept] = Math.round(c);
            }
            rows[kept] = row;
            kept++;
        }
        Rule rule = new Rule();
        rule.A = new long[kept][];
        rule.C = new long[kept];
        for (int j = 0; j < kept; j++) {
            rule.A[j] = rows[j];
            rule.C[j] = rhs[j];
        }
        return rule;
    }

    /** Peak live coefficient count, threaded out of the series by reference. */
    private static final class Peak {
        long value;
    }

    /**
     * Coefficient-domain evaluation of the J-fold contour integral at the
     * right-hand side C, which is the full rule for g(C) and the rule shifted by
     * one class requirement for g(C - A e_r).
     *
     * The series lives in a flat array indexed column-major over the live links,
     * stride[k] = prod_{i&lt;k} curdim[i], with curdim[j] == 1 while link j is not
     * live and C[j]+1 while it is. f[r] holds the scaled terms nu_r^n / n! for
     * n = 0..nmaxFull[r].
     */
    private static double series(double[][] f, long[][] A, long[] C, long[] nmaxFull,
                                 long maxLiveStates, Peak peak) {
        int R = f.length;
        int J = C.length;

        // The elimination order: link j is created at its first route and summed
        // out after its last, so only an induced width of links is ever live.
        int[] first = new int[J];
        int[] last = new int[J];
        for (int j = 0; j < J; j++) {
            boolean seen = false;
            for (int r = 0; r < R; r++) {
                if (A[j][r] == 0) {
                    continue;
                }
                if (!seen) {
                    first[j] = r;
                    seen = true;
                }
                last[j] = r;
            }
            if (!seen) {
                throw new RuntimeException("lossn_manjunath: a constraint row with no nonzero entry "
                        + "reached the series; the rule was not reduced");
            }
        }

        int[] curdim = new int[J];
        for (int j = 0; j < J; j++) {
            curdim[j] = 1;
        }
        double[] ser = new double[1];
        ser[0] = 1.0;

        for (int r = 0; r < R; r++) {
            // 1. Create the links whose first route is this one, keeping the
            // existing content at degree zero in the new variable.
            for (int j = 0; j < J; j++) {
                if (first[j] != r) {
                    continue;
                }
                int newdim = (int) C[j] + 1;
                long pre = 1;
                long post = 1;
                for (int k = 0; k < j; k++) {
                    pre *= curdim[k];
                }
                for (int k = j + 1; k < J; k++) {
                    post *= curdim[k];
                }
                if (pre * post > maxLiveStates / newdim) {
                    throw new RuntimeException("lossn_manjunath: the exact transform would hold more than "
                            + maxLiveStates + " series coefficients at once. Peak memory is the "
                            + "product of (C_j+1) over the links live at the same time, so a wide "
                            + "constraint row with a large capacity is what costs; raise "
                            + "DEFAULT_MAX_LIVE_STATES deliberately, or use lossn_mci, which is "
                            + "unbiased at any size");
                }
                double[] grown = new double[(int) (pre * newdim * post)];
                for (int q = 0; q < post; q++) {
                    for (int p = 0; p < pre; p++) {
                        grown[(int) (p + q * pre * newdim)] = ser[(int) (p + q * pre)];
                    }
                }
                ser = grown;
                curdim[j] = newdim;
                if (ser.length > peak.value) {
                    peak.value = ser.length;
                }
            }

            // 2. Multiply in route r.
            boolean constrained = false;
            for (int j = 0; j < J; j++) {
                if (A[j][r] != 0) {
                    constrained = true;
                }
            }

            if (!constrained) {
                // Bounded by no link, so its z-transform is a constant: the whole
                // truncated sequence sums into the series. For a route absent from
                // every row this is the factor exp(nu_r) the caller folded into
                // lGfree, hence a multiply by 1.
                double s = 0.0;
                for (long n = 0; n <= nmaxFull[r] && n < f[r].length; n++) {
                    s += f[r][(int) n];
                }
                for (int i = 0; i < ser.length; i++) {
                    ser[i] *= s;
                }
                continue;
            }

            // The degree of route r is capped by every row it appears in,
            // evaluated at THIS right-hand side: the shifted series for
            // g(C - A e_r) admits strictly fewer calls than g(C).
            long nmax = Math.min(nmaxFull[r], f[r].length - 1);
            for (int j = 0; j < J; j++) {
                if (A[j][r] > 0) {
                    nmax = Math.min(nmax, C[j] / A[j][r]);
                }
            }

            long[] stride = new long[J];
            stride[0] = 1;
            for (int k = 1; k < J; k++) {
                stride[k] = stride[k - 1] * curdim[k - 1];
            }
            int P = ser.length;
            double[] next = new double[P];
            int[] sub = new int[J];
            for (long n = 0; n <= nmax; n++) {
                double c = f[r][(int) n];
                if (c == 0.0) {
                    continue;
                }
                if (n == 0) {
                    for (int i = 0; i < P; i++) {
                        next[i] += c * ser[i];
                    }
                    continue;
                }
                // Shift by n requirement vectors, dropping the coefficients the
                // shift would push past the capacity. Those monomials can never
                // contribute to an extracted coefficient, which is exactly why the
                // truncation is exact and not an approximation.
                for (int j = 0; j < J; j++) {
                    sub[j] = 0;
                }
                boolean anyok = false;
                for (int i = 0; i < P; i++) {
                    boolean ok = true;
                    long tgt = 0;
                    for (int j = 0; j < J && ok; j++) {
                        long d = sub[j] + A[j][r] * n;
                        if (d > C[j]) {
                            ok = false;
                        } else {
                            tgt += d * stride[j];
                        }
                    }
                    if (ok) {
                        next[(int) tgt] += c * ser[i];
                        anyok = true;
                    }
                    // Odometer over the live grid, in the same column-major order
                    // the strides encode.
                    for (int j = 0; j < J; j++) {
                        if (++sub[j] < curdim[j]) {
                            break;
                        }
                        sub[j] = 0;
                    }
                }
                // The shift only grows with n, so once nothing fits nothing will.
                if (!anyok) {
                    break;
                }
            }
            ser = next;

            // 3. Integrate out the links whose last route was this one. The
            // multiplier (z^{C+1}-1)/(z-1) of a '<=' constraint turns the residue
            // into the partial sum of the coefficients of degrees 0..C_j, which is
            // the sum along that dimension.
            for (int j = 0; j < J; j++) {
                if (last[j] != r) {
                    continue;
                }
                long pre = 1;
                long post = 1;
                for (int k = 0; k < j; k++) {
                    pre *= curdim[k];
                }
                for (int k = j + 1; k < J; k++) {
                    post *= curdim[k];
                }
                int dj = curdim[j];
                double[] summed = new double[(int) (pre * post)];
                for (int q = 0; q < post; q++) {
                    for (int d = 0; d < dj; d++) {
                        for (int p = 0; p < pre; p++) {
                            summed[(int) (p + q * pre)] += ser[(int) (p + d * pre + q * pre * dj)];
                        }
                    }
                }
                ser = summed;
                curdim[j] = 1;
            }
        }

        if (ser.length != 1) {
            throw new RuntimeException("lossn_manjunath: a link was never integrated out; the elimination "
                    + "order is inconsistent with the constraint rows");
        }
        return ser[0];
    }

    /**
     * Exact normalization constant, carried load and blocking of a loss network.
     *
     * @param nuVec Offered load of route r, nonnegative (1xR).
     * @param Amat  Capacity requirement of link j for route r (JxR nonnegative integers).
     * @param cVec  Available capacity of link j (Jx1 nonnegative integers).
     */
    public static Ret.lossnManjunath lossn_manjunath(Matrix nuVec, Matrix Amat, Matrix cVec) {
        return lossn_manjunath(nuVec, Amat, cVec, DEFAULT_MAX_LIVE_STATES);
    }

    /**
     * Exact normalization constant, carried load and blocking of a loss network.
     *
     * @param nuVec          Offered load of route r, nonnegative (1xR).
     * @param Amat           Capacity requirement of link j for route r (JxR).
     * @param cVec           Available capacity of link j (Jx1).
     * @param maxLiveStates  Cap on the simultaneously live series coefficients.
     */
    public static Ret.lossnManjunath lossn_manjunath(Matrix nuVec, Matrix Amat, Matrix cVec,
                                       long maxLiveStates) {
        double[] nu = nuVec.toArray1D();
        double[] C = cVec.toArray1D();
        int R = nu.length;
        int Jin = C.length;
        double[][] A = (Jin > 0) ? Amat.toArray2D() : new double[0][0];
        if (Jin > 0 && (A.length != Jin || A[0].length != R)) {
            throw new RuntimeException("lossn_manjunath: A must be " + Jin + "x" + R + " (J x R)");
        }
        for (int r = 0; r < R; r++) {
            if (nu[r] < 0.0) {
                throw new RuntimeException("lossn_manjunath: nu must be nonnegative");
            }
        }

        double[] qLen = new double[R];
        double[] loss = new double[R];

        Rule rule = integralize(A, C);
        int J = rule.C.length;

        // A route absent from every remaining row never blocks: its marginal is an
        // untruncated Poisson, so it carries its full offered load and factors
        // exp(nu_r) out of g(C).
        boolean[] free = new boolean[R];
        for (int r = 0; r < R; r++) {
            free[r] = true;
        }
        for (int j = 0; j < J; j++) {
            for (int r = 0; r < R; r++) {
                if (rule.A[j][r] != 0) {
                    free[r] = false;
                }
            }
        }
        double lGfree = 0.0;
        boolean allFree = true;
        for (int r = 0; r < R; r++) {
            if (free[r]) {
                qLen[r] = nu[r];
                lGfree += nu[r];
            } else {
                allFree = false;
            }
        }
        if (J == 0 || allFree) {
            return new Ret.lossnManjunath(new Matrix(qLen), new Matrix(loss), lGfree, 1, 0L);
        }

        // Per-route truncation, and the terms f_r(n) = nu_r^n / n!. Each sequence
        // is built in the LOG domain and shifted by its own maximum before
        // exponentiating, so its largest entry is exactly 1; the scale cancels in
        // every ratio below and is added back in log g. Shifting before the
        // exponential rather than after is what keeps a heavy route in range:
        // nu^n/n! peaks near n = nu at roughly exp(nu)/sqrt(2 pi nu), so forming
        // the terms first and dividing by their maximum afterwards overflows to
        // infinity for a load above about 700, and infinity is not finite, so a
        // rescaling step guarded on finiteness would then decline to run and the
        // infinity would reach g, returning NaN for every metric. Measured at
        // nu = C = 900.
        long[] nmax = new long[R];
        double[][] f = new double[R][];
        double logscale = 0.0;
        for (int r = 0; r < R; r++) {
            if (free[r]) {
                f[r] = new double[1];
                f[r][0] = 1.0;
                continue;
            }
            long v = Long.MAX_VALUE;
            for (int j = 0; j < J; j++) {
                if (rule.A[j][r] > 0) {
                    v = Math.min(v, rule.C[j] / rule.A[j][r]);
                }
            }
            nmax[r] = v;
            int len = (int) v + 1;
            if (nu[r] == 0.0) {
                // No offered load: only the empty term survives, and log(0) has no
                // shift to take.
                f[r] = new double[len];
                f[r][0] = 1.0;
                continue;
            }
            double lnu = FastMath.log(nu[r]);
            double[] lf = new double[len];
            double m = Double.NEGATIVE_INFINITY;
            for (int n = 0; n < len; n++) {
                lf[n] = n * lnu - Gamma.logGamma(n + 1.0);
                if (lf[n] > m) {
                    m = lf[n];
                }
            }
            f[r] = new double[len];
            for (int n = 0; n < len; n++) {
                f[r][n] = FastMath.exp(lf[n] - m);
            }
            logscale += m;
        }

        Peak peak = new Peak();
        double G = series(f, rule.A, rule.C, nmax, maxLiveStates, peak);
        if (G <= 0.0) {
            throw new RuntimeException("lossn_manjunath: the admissible set is empty -- no state satisfies "
                    + "A n <= C, so the loss network has no stationary distribution");
        }
        double lG = FastMath.log(G) + logscale + lGfree;

        for (int r = 0; r < R; r++) {
            if (free[r]) {
                continue;
            }
            long[] Cr = new long[J];
            boolean overflows = false;
            for (int j = 0; j < J; j++) {
                Cr[j] = rule.C[j] - rule.A[j][r];
                if (Cr[j] < 0) {
                    overflows = true;
                }
            }
            if (overflows) {
                // A single class r call already exceeds a capacity, so the route is
                // blocked in every state, including the empty one.
                loss[r] = 1.0;
                qLen[r] = 0.0;
                continue;
            }
            double Gr = series(f, rule.A, Cr, nmax, maxLiveStates, peak);
            // The scale cancels here, which is what lets the terms be normalized.
            double ratio = Gr / G;
            qLen[r] = nu[r] * ratio;
            loss[r] = 1.0 - ratio;
        }

        return new Ret.lossnManjunath(new Matrix(qLen), new Matrix(loss), lG, 1, peak.value);
    }
}
