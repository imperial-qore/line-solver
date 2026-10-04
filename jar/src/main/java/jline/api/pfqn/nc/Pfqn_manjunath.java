/**
 * @file Exact normalizing constant of a product-form network over a state space
 *       cut by linear integer constraints (Manjunath-Sikdar)
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Exact normalizing constant of a closed multiclass product-form network whose
 * state space carries arbitrary linear integer constraints.
 *
 * This is the queueing-network half of the transform technique of Manjunath and
 * Sikdar, of which {@link jline.api.lossn.Lossn_manjunath} is the loss-network
 * half. The two solve the same problem -- sum a product form over an irregular
 * integer state space -- from opposite ends of the paper: Lossn_manjunath
 * implements Section 2.2, a set of '&lt;=' rows over the Poisson terms nu^n/n!,
 * while this routine implements Section 3 together with Section 5.3, a MIXED set
 * of '=', '&lt;=' and '&gt;' rows over the BCMP terms, where the population
 * constraint of a closed network is itself one of the equalities.
 *
 * THE MODEL. M queueing stations (FCFS, PS or LCFS; rows of L) and Mz delay
 * stations (rows of Z) serve R closed classes with populations N. Writing n_ir
 * for the class r jobs at station i and n_i = sum_r n_ir, the BCMP product form
 * of Baskett-Chandy-Muntz-Palacios is
 *
 * <pre>
 *   p(n) = (1/G) prod_{i queueing} n_i! prod_r L_ir^{n_ir}/n_ir!
 *                prod_{i delay}         prod_r Z_ir^{n_ir}/n_ir!
 * </pre>
 *
 * Every state obeys the R population equalities sum_i n_ir = N_r; on top of
 * those the caller may impose any number of further rows
 *
 * <pre>
 *   sum_{i,r} A(j, i + S*r) n_ir  {=, &lt;=, &gt;}  b(j),   S = M + Mz,
 * </pre>
 *
 * i.e. A acts on the (M+Mz)-by-R occupancy read column by column with the
 * queueing stations first. With no extra rows the routine returns exactly the
 * normalizing constant of {@link Pfqn_ca}, which is the parity oracle used by
 * the tests; with extra rows it answers a question no other routine in the pfqn
 * family can, the convolution and MVA recursions having nowhere to carry a
 * second constraint.
 *
 * WHY THE GENERATING FUNCTION IS A PRODUCT, AND WHERE THE n_i! GOES. Marking
 * class r by z_r and row j by y_j, and abbreviating the monomial one class r job
 * at station i contributes as u_ir = z_r prod_j y_j^{A(j, i + S*r)}, the sum over
 * the occupancies of a single QUEUEING station is, by the multinomial theorem,
 *
 * <pre>
 *   sum_{n_i.} n_i! prod_r (L_ir u_ir)^{n_ir}/n_ir!
 *     = sum_k (sum_r L_ir u_ir)^k = 1 / (1 - sum_r L_ir u_ir),
 * </pre>
 *
 * so the n_i! that couples the classes at a queueing station is exactly what
 * turns its factor from an exponential into a geometric one. The paper reaches
 * the same place through the Euler integral n! = int_0^inf e^-t t^n dt (Eqns
 * 16-18), which is that geometric series evaluated; the closed form is used here
 * because there is then no quadrature to discretize. A DELAY station has no n_i!
 * and keeps its exponential. Hence
 *
 * <pre>
 *   F(z,y) = prod_i 1/(1 - sum_r L_ir u_ir) prod_k prod_r exp(Z_kr u_kr)
 * </pre>
 *
 * and G is read off F as a coefficient: degree exactly N_r in z_r, and for row j
 * the degree its sense dictates -- exactly b_j for '=', the sum of degrees
 * 0..b_j for '&lt;=' (the multiplier (y^{b+1}-1)/(y-1) of Eqn 5, whose residue is
 * that partial sum), and the complement of the latter for '&gt;' (Eqn 6).
 *
 * WHY IT IS A COEFFICIENT COMPUTATION AND NOT A QUADRATURE. The contour
 * integrals of Eqn 9 all have their only pole at the origin, of order one more
 * than the right-hand side, so each is a residue and hence a Taylor coefficient.
 * The routine therefore never integrates: it carries F as a multivariate power
 * series truncated at degree N_r in z_r and b_j in y_j. Truncation is exact
 * because A is nonnegative -- no monomial above a cut can be brought back down
 * by a later factor.
 *
 * Each queueing station is applied by SOLVING (1 - sum_r L_ir u_ir) x = ser
 * rather than by expanding the geometric series, which keeps the cost at one
 * pass. Every monomial of the operator carries z_r to a strictly higher power, so
 * sweeping the lattice in increasing flat index lets each coefficient read only
 * coefficients already final: a Gauss-Seidel sweep whose result is the exact
 * solve, not an iterate. A delay station has no such recurrence and is convolved
 * with exp term by term, which is where its extra factor of the population in
 * the cost comes from.
 *
 * THE ELIMINATION ORDER IS THE MEMORY BOUND. Variable y_j is created when the
 * first station its row touches is multiplied in and discharged immediately
 * after the last, so peak memory is prod_r (N_r+1) times the product of (b_j+1)
 * over the SIMULTANEOUSLY LIVE rows, not over all rows. A row constraining one
 * station therefore costs essentially nothing. The class axes are live
 * throughout, so prod_r (N_r+1) is a floor -- the same lattice Pfqn_ca walks.
 *
 * SCOPE. Load-dependent and multiserver stations are NOT covered: their
 * per-station term is not geometric, and while the paper admits an arbitrary
 * f_i(n_i) in the single-class case (Section 2), the multiclass n_i! coupling
 * used above then breaks. Use Pfqn_gld or Pfqn_conwayms for those. A and b must
 * be integer valued and A nonnegative, since the residue argument counts whole
 * units; a fractional entry is refused rather than rounded.
 *
 * Reference: D. Manjunath and B. Sikdar, Integral Expressions for the Numerical
 * Evaluation of Product Form Expressions Over Irregular Multidimensional Integer
 * Spaces. Sections 3 and 5.3.
 */
public final class Pfqn_manjunath {
    private Pfqn_manjunath() {}

    /** Exact normalizing constant of an unconstrained closed network. */
    public static Ret.pfqnManjunath pfqn_manjunath(Matrix L, Matrix N) {
        return pfqn_manjunath(L, N, null, null, null, null);
    }

    /** Exact normalizing constant of an unconstrained closed network with delay. */
    public static Ret.pfqnManjunath pfqn_manjunath(Matrix L, Matrix N, Matrix Z) {
        return pfqn_manjunath(L, N, Z, null, null, null);
    }

    /**
     * Exact normalizing constant of a constrained closed product-form network.
     *
     * @param L     service demand of class r at queueing station i (MxR)
     * @param N     population of class r (1xR nonnegative integers)
     * @param Z     think time of class r at delay station k (MzxR), may be null
     * @param A     extra constraint coefficients on the occupancy read column by
     *              column (Jx((M+Mz)*R) nonnegative integers), may be null
     * @param b     extra constraint right-hand sides (Jx1 integers), may be null
     * @param sense one character per row, 'E' (=), 'L' (&lt;=) or 'G' (&gt;);
     *              null means all 'L'
     */
    public static Ret.pfqnManjunath pfqn_manjunath(Matrix L, Matrix N, Matrix Z,
                                                   Matrix A, Matrix b, String sense) {
        return pfqn_manjunath(L, N, Z, A, b, sense, false);
    }

    /**
     * As above, additionally returning the per-class decomposition when
     * {@code stats} is set. That requires the ONE configuration in which the
     * truncated product form is the EXACT stationary law: a single queueing
     * station inside the region and a SINGLE DELAY STATION OUTSIDE IT. Anything
     * else is refused by name rather than answered wrongly.
     *
     * <p>WHY THE CONFIGURATION IS NOT A CONVENIENCE. With one queueing station
     * the state is the queue occupancy alone (the delay holds the complement)
     * and every transition moves one job of one class by one unit, so the chain
     * is a multidimensional birth-death process. That process is reversible, and
     * Kelly's truncation theorem then applies verbatim: restricting it to the
     * coordinate-convex set A n &lt;= b and renormalizing gives exactly the
     * truncated product form. Add a second queueing station and the
     * delay -&gt; q1 -&gt; q2 -&gt; delay cycle destroys reversibility; truncation no
     * longer preserves the product form, measured at 131% relative error on the
     * stationary law of a 2-class, N = [2 2] instance. G and lG stay correct as
     * a sum over the admissible set in every configuration; only the metrics are
     * withheld.</p>
     *
     * <p>Everything follows from two ratios of normalizing constants, both taken
     * in the log domain so the internal rescaling cancels without being
     * reconstructed: X_r = G(N - e_r ; b - A(:,qcol_r)) / G(N ; b), and
     * P(n_qr = k) = G(N ; b plus the row n_qr = k) / G(N ; b). The first is the
     * loss network's g(C - A e_r) in another guise.</p>
     */
    public static Ret.pfqnManjunath pfqn_manjunath(Matrix L, Matrix N, Matrix Z,
                                                   Matrix A, Matrix b, String sense,
                                                   boolean stats) {
        double[] Nd = N.toArray1D();
        int R = Nd.length;

        double[][] Lm = toRows(L, R, "L");
        double[][] Zm = toRows(Z, R, "Z");
        int M = Lm.length;
        int Mz = Zm.length;
        int S = M + Mz;

        double[] bd = (b == null) ? new double[0] : b.toArray1D();
        int J = bd.length;
        double[][] Am;
        if (J == 0) {
            Am = new double[0][S * R];
        } else {
            Am = A.toArray2D();
            if (Am.length != J || Am[0].length != S * R) {
                throw new RuntimeException("pfqn_manjunath: A must be " + J + "x" + (S * R)
                        + " (J x (M+Mz)*R), acting on the occupancy read column by column");
            }
        }
        char[] sn = new char[J];
        if (sense == null) {
            for (int j = 0; j < J; j++) {
                sn[j] = 'L';
            }
        } else {
            String su = sense.toUpperCase();
            if (su.length() != J) {
                throw new RuntimeException("pfqn_manjunath: sense must have one character per "
                        + "row of b");
            }
            for (int j = 0; j < J; j++) {
                sn[j] = su.charAt(j);
                if (sn[j] != 'E' && sn[j] != 'L' && sn[j] != 'G') {
                    throw new RuntimeException("pfqn_manjunath: sense must contain only 'E' (=), "
                            + "'L' (<=) or 'G' (>)");
                }
            }
        }

        int[] Ni = new int[R];
        for (int r = 0; r < R; r++) {
            if (Nd[r] < 0.0) {
                return empty(0.0, Double.NEGATIVE_INFINITY, 0L, stats, R);
            }
            if (Math.abs(Nd[r] - Math.rint(Nd[r])) > 1e-9) {
                throw new RuntimeException("pfqn_manjunath: N must contain nonnegative integers");
            }
            Ni[r] = (int) Math.round(Nd[r]);
        }
        long[][] Ai = new long[J][];
        long[] bi = new long[J];
        for (int j = 0; j < J; j++) {
            Ai[j] = new long[S * R];
            for (int c = 0; c < S * R; c++) {
                double a = Am[j][c];
                if (a < 0.0 || Math.abs(a - Math.rint(a)) > 1e-9) {
                    throw new RuntimeException("pfqn_manjunath: A must contain nonnegative "
                            + "integers; the residue argument counts whole units");
                }
                Ai[j][c] = Math.round(a);
            }
            if (Math.abs(bd[j] - Math.rint(bd[j])) > 1e-9) {
                throw new RuntimeException("pfqn_manjunath: b must contain integers; the residue "
                        + "argument counts whole units");
            }
            bi[j] = Math.round(bd[j]);
        }

        // A row of zeros constrains nothing, so it is decided here rather than
        // carried as a one-coefficient dimension: 0 = b, 0 <= b and 0 > b are
        // each settled by the sign of b alone. Same for a negative right-hand
        // side, which no nonnegative combination can meet ('E', 'L') or can fail
        // to beat ('G').
        boolean[] keep = new boolean[J];
        int kept = 0;
        for (int j = 0; j < J; j++) {
            boolean trivial = true;
            for (int c = 0; c < S * R; c++) {
                if (Ai[j][c] != 0) {
                    trivial = false;
                    break;
                }
            }
            keep[j] = true;
            if (sn[j] == 'E') {
                if (bi[j] < 0 || (trivial && bi[j] != 0)) {
                    return empty(0.0, Double.NEGATIVE_INFINITY, 0L, stats, R);
                }
                if (trivial) {
                    keep[j] = false;
                }
            } else if (sn[j] == 'L') {
                if (bi[j] < 0) {
                    return empty(0.0, Double.NEGATIVE_INFINITY, 0L, stats, R);
                }
                if (trivial) {
                    keep[j] = false;
                }
            } else {
                if (bi[j] < 0) {
                    keep[j] = false;                 // 0 > negative always holds
                } else if (trivial) {
                    return empty(0.0, Double.NEGATIVE_INFINITY, 0L, stats, R);
                }
            }
            if (keep[j]) {
                kept++;
            }
        }
        long[][] Ak = new long[kept][];
        long[] bk = new long[kept];
        char[] sk = new char[kept];
        int w = 0;
        for (int j = 0; j < J; j++) {
            if (keep[j]) {
                Ak[w] = Ai[j];
                bk[w] = bi[j];
                sk[w] = sn[j];
                w++;
            }
        }
        J = kept;

        if (S == 0) {
            // No station: the only state is empty, admissible when every class is too.
            boolean isEmpty = true;
            for (int r = 0; r < R; r++) {
                if (Ni[r] != 0) {
                    isEmpty = false;
                }
            }
            return isEmpty ? empty(1.0, 0.0, 1L, stats, R)
                           : empty(0.0, Double.NEGATIVE_INFINITY, 0L, stats, R);
        }

        // Every monomial that survives the extraction has total degree sum(N) in
        // the demands, so a common rescaling of L and Z moves lG by a known
        // amount and nothing else. The exponent is chosen as in Pfqn_ca, from the
        // largest single term the network can produce, so the series is centred
        // near 1.
        int Nt = 0;
        for (int r = 0; r < R; r++) {
            Nt += Ni[r];
        }
        double cscale = 1.0;
        if (Nt > 0) {
            double lGest = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < M; i++) {
                double t = 0.0;
                boolean ok = true;
                for (int r = 0; r < R && ok; r++) {
                    if (Ni[r] > 0) {
                        if (Lm[i][r] > 0.0) {
                            t += Ni[r] * FastMath.log(Lm[i][r]);
                        } else {
                            ok = false;
                        }
                    }
                }
                if (ok && t > lGest) {
                    lGest = t;
                }
            }
            if (Mz > 0) {
                double t = 0.0;
                boolean ok = true;
                for (int r = 0; r < R && ok; r++) {
                    if (Ni[r] > 0) {
                        double zs = 0.0;
                        for (int k = 0; k < Mz; k++) {
                            zs += Zm[k][r];
                        }
                        if (zs > 0.0) {
                            t += Ni[r] * FastMath.log(zs) - Gamma.logGamma(Ni[r] + 1.0);
                        } else {
                            ok = false;
                        }
                    }
                }
                if (ok && t > lGest) {
                    lGest = t;
                }
            }
            if (!Double.isInfinite(lGest) && !Double.isNaN(lGest)) {
                cscale = Math.scalb(1.0, matlabRound(lGest / (Nt * FastMath.log(2.0))));
            }
        }
        double[][] Ls = new double[M][R];
        double[][] Zs = new double[Mz][R];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Ls[i][r] = Lm[i][r] / cscale;
            }
        }
        for (int k = 0; k < Mz; k++) {
            for (int r = 0; r < R; r++) {
                Zs[k][r] = Zm[k][r] / cscale;
            }
        }

        // A '>' row is the complement of a '<=' row at the same right-hand side,
        // which is how the paper discharges it (Eqn 6). With several such rows
        // the product of the complements expands by inclusion-exclusion, so the
        // series is evaluated once per subset of them, with the subset's rows
        // re-entered as '<=' and the rest dropped. Exact, and the only place the
        // cost is exponential -- in the number of '>' rows, normally zero.
        int nGt = 0;
        for (int j = 0; j < J; j++) {
            if (sk[j] == 'G') {
                nGt++;
            }
        }
        int[] gt = new int[nGt];
        int g = 0;
        for (int j = 0; j < J; j++) {
            if (sk[j] == 'G') {
                gt[g++] = j;
            }
        }
        double Gs = 0.0;
        long peak = 0L;
        for (int mask = 0; mask < (1 << nGt); mask++) {
            boolean[] on = new boolean[J];
            for (int j = 0; j < J; j++) {
                on[j] = (sk[j] != 'G');
            }
            for (int t = 0; t < nGt; t++) {
                if ((mask & (1 << t)) != 0) {
                    on[gt[t]] = true;
                }
            }
            int cnt = 0;
            for (int j = 0; j < J; j++) {
                if (on[j]) {
                    cnt++;
                }
            }
            long[][] As = new long[cnt][];
            long[] bs = new long[cnt];
            char[] ss = new char[cnt];
            int q = 0;
            for (int j = 0; j < J; j++) {
                if (on[j]) {
                    As[q] = Ak[j];
                    bs[q] = bk[j];
                    ss[q] = (sk[j] == 'G') ? 'L' : sk[j];
                    q++;
                }
            }
            Peak pk = new Peak();
            double gv = series(Ls, Zs, Ni, As, bs, ss, pk);
            int bits = Integer.bitCount(mask);
            Gs += ((bits % 2 == 0) ? 1.0 : -1.0) * gv;
            if (pk.value > peak) {
                peak = pk.value;
            }
        }

        if (Gs <= 0.0) {
            // Either the admissible set is empty or the '>' complements cancelled it.
            return empty(0.0, Double.NEGATIVE_INFINITY, peak, stats, R);
        }
        double lG = FastMath.log(Gs) + Nt * FastMath.log(cscale);
        Ret.pfqnManjunath out = new Ret.pfqnManjunath(FastMath.exp(lG), lG, peak);
        if (stats) {
            fillStats(out, L, N, Z, Ak, bk, sk, lG, M, Mz, S, R);
        }
        return out;
    }


    /** Early exit, carrying an all-zero decomposition when one was asked for. */
    private static Ret.pfqnManjunath empty(double G, double lG, long peak, boolean stats, int R) {
        Ret.pfqnManjunath out = new Ret.pfqnManjunath(G, lG, peak);
        if (stats) {
            out.Q = new Matrix(1, R);
            out.X = new Matrix(1, R);
            out.U = new Matrix(1, R);
            out.think = new Matrix(1, R);
            out.blocked = new Matrix(1, R);
            out.delay = new Matrix(1, R);
        }
        return out;
    }

    /**
     * Per-class decomposition, for the one configuration in which the truncated
     * product form is the exact stationary law: a single queueing station inside
     * the region and a single delay station outside it.
     *
     * <p>A refused admission is a DELETED transition, so a blocked job never
     * leaves the delay, and because the think time is exponential a held job is
     * indistinguishable from one still thinking. The delay population carries
     * both and Little's law separates them. This is NOT the WAITQ rule of
     * SolverSSA/SolverCTMC/JMT, which moves a refused job out of the delay into
     * a per-region FIFO counted at no station.</p>
     */
    private static void fillStats(Ret.pfqnManjunath out, Matrix L, Matrix N, Matrix Z,
                                  long[][] A, long[] b, char[] sense, double lG,
                                  int M, int Mz, int S, int R) {
        if (Mz != 1) {
            throw new RuntimeException("pfqn_manjunath: the per-class decomposition needs "
                    + "exactly one delay station, got " + Mz + ". Pass Z as a 1xR row of "
                    + "think times");
        }
        if (M != 1) {
            throw new RuntimeException("pfqn_manjunath: the per-class decomposition needs "
                    + "exactly one queueing station, got " + M + ". With two or more the "
                    + "delay->q1->q2->delay cycle makes the chain irreversible, Kelly "
                    + "truncation no longer holds, and the truncated product form is not the "
                    + "stationary law (measured at 131% error). G and lG are still returned "
                    + "and still correct as a sum over the admissible set");
        }
        int J = b.length;
        // The delay must sit OUTSIDE the region: its columns are S*r + (S-1).
        for (int r = 0; r < R; r++) {
            int dcol = S * r + (S - 1);
            for (int j = 0; j < J; j++) {
                if (A[j][dcol] != 0) {
                    throw new RuntimeException("pfqn_manjunath: constraint row(s) reference the "
                            + "delay station in class " + (r + 1) + " (column " + dcol + "). The "
                            + "delay must lie OUTSIDE the finite capacity region, because the "
                            + "decomposition charges every held job to it");
                }
            }
        }

        Matrix Am = new Matrix(J, S * R);
        for (int j = 0; j < J; j++) {
            for (int c = 0; c < S * R; c++) {
                Am.set(j, c, (double) A[j][c]);
            }
        }
        String sn = new String(sense);

        // Ratios of normalizing constants are taken in the LOG domain, so the
        // internal power-of-two rescaling cancels without being reconstructed.
        Matrix Q = new Matrix(1, R);
        Matrix X = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            int qcol = S * r;               // column of (queueing station, class r)

            // Throughput. One class r job removed from the queue leaves a state of
            // population N - e_r whose admission rule is shifted by that job's own
            // requirement column, exactly as the loss network's g(C - A e_r).
            if (N.get(0, r) >= 1.0) {
                Matrix Nr = N.copy();
                Nr.set(0, r, N.get(0, r) - 1.0);
                Matrix br = new Matrix(J, 1);
                for (int j = 0; j < J; j++) {
                    br.set(j, 0, (double) (b[j] - A[j][qcol]));
                }
                Ret.pfqnManjunath sub = pfqn_manjunath(L, Nr, Z, Am, br, sn, false);
                if (!Double.isInfinite(sub.lG)) {
                    X.set(0, r, FastMath.exp(sub.lG - lG));
                }
            }

            // Mean queue length from the marginal law. An '=' row is discharged
            // by picking a single coefficient, so each call returns the mass of
            // exactly that occupancy.
            int Nr = (int) Math.round(N.get(0, r));
            double q = 0.0;
            for (int k = 1; k <= Nr; k++) {
                Matrix Ak = new Matrix(J + 1, S * R);
                for (int j = 0; j < J; j++) {
                    for (int c = 0; c < S * R; c++) {
                        Ak.set(j, c, (double) A[j][c]);
                    }
                }
                Ak.set(J, qcol, 1.0);
                Matrix bk = new Matrix(J + 1, 1);
                for (int j = 0; j < J; j++) {
                    bk.set(j, 0, (double) b[j]);
                }
                bk.set(J, 0, (double) k);
                Ret.pfqnManjunath sub = pfqn_manjunath(L, N, Z, Ak, bk, sn + "E", false);
                if (!Double.isInfinite(sub.lG)) {
                    q += k * FastMath.exp(sub.lG - lG);
                }
            }
            Q.set(0, r, q);
        }

        Matrix U = new Matrix(1, R);
        Matrix think = new Matrix(1, R);
        Matrix delay = new Matrix(1, R);
        Matrix blocked = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            U.set(0, r, X.get(0, r) * L.get(0, r));       // one server, one visit
            think.set(0, r, X.get(0, r) * Z.get(0, r));   // Little's law at the delay
            delay.set(0, r, N.get(0, r) - Q.get(0, r));   // all that is not at the queue
            blocked.set(0, r, delay.get(0, r) - think.get(0, r));
        }
        out.Q = Q;
        out.X = X;
        out.U = U;
        out.think = think;
        out.blocked = blocked;
        out.delay = delay;
    }

    /** Peak live coefficient count, threaded out of the series by reference. */
    private static final class Peak {
        long value;
    }

    /**
     * MATLAB's round: half away from zero. Java's Math.round is half UP, so it
     * disagrees on negative halves, and the reference is MATLAB.
     */
    private static int matlabRound(double x) {
        return (int) (x >= 0.0 ? Math.floor(x + 0.5) : Math.ceil(x - 0.5));
    }

    private static double[][] toRows(Matrix X, int R, String name) {
        if (X == null || X.getNumRows() == 0 || X.getNumCols() == 0) {
            return new double[0][R];
        }
        if (X.getNumCols() != R) {
            throw new RuntimeException("pfqn_manjunath: " + name + " must have " + R
                    + " columns, one per class");
        }
        double[][] out = new double[X.getNumRows()][R];
        for (int i = 0; i < X.getNumRows(); i++) {
            for (int r = 0; r < R; r++) {
                out[i][r] = X.get(i, r);
            }
        }
        return out;
    }

    /**
     * Coefficient-domain evaluation of the multiple contour integral of Eqn 9 for
     * a set of '=' and '<=' rows. The series lives in a flat array indexed
     * column-major, its first R axes the class markers z_r of extent N_r+1
     * throughout and its remaining J axes the row markers y_j, of extent 1 while
     * row j is not live and b_j+1 while it is.
     */
    private static double series(double[][] L, double[][] Z, int[] N, long[][] A, long[] b,
                                 char[] sense, Peak peak) {
        int M = L.length;
        int Mz = Z.length;
        int S = M + Mz;
        int R = N.length;
        int J = b.length;
        int D = R + J;

        // Row j is created at the first station it touches and discharged after
        // the last, so only an induced width of rows is ever live. A row that
        // reached here touches at least one station, the trivial ones having been
        // decided already.
        int[] first = new int[J];
        int[] last = new int[J];
        for (int j = 0; j < J; j++) {
            boolean seen = false;
            for (int i = 0; i < S; i++) {
                boolean touch = false;
                for (int r = 0; r < R; r++) {
                    if (A[j][i + S * r] != 0) {
                        touch = true;
                    }
                }
                if (!touch) {
                    continue;
                }
                if (!seen) {
                    first[j] = i;
                    seen = true;
                }
                last[j] = i;
            }
            if (!seen) {
                throw new RuntimeException("pfqn_manjunath: a constraint row with no nonzero "
                        + "entry reached the series; the rule was not reduced");
            }
        }

        int[] dims = new int[D];
        for (int r = 0; r < R; r++) {
            dims[r] = N[r] + 1;
        }
        for (int j = 0; j < J; j++) {
            dims[R + j] = 1;
        }
        int P = 1;
        for (int d = 0; d < D; d++) {
            P *= dims[d];
        }
        double[] ser = new double[P];
        ser[0] = 1.0;
        if (P > peak.value) {
            peak.value = P;
        }

        for (int i = 0; i < S; i++) {
            for (int j = 0; j < J; j++) {
                if (first[j] == i) {
                    int newdim = (int) b[j] + 1;
                    ser = expand(ser, dims, R + j, newdim);
                    dims[R + j] = newdim;
                    if (ser.length > peak.value) {
                        peak.value = ser.length;
                    }
                }
            }

            int[] stride = strides(dims);
            // The monomial one class r job at station i contributes: z_r gains one
            // degree and y_j gains A(j, i + S*r). A row not live at this station
            // has a zero entry here by construction of first/last, so a dead axis
            // is never shifted.
            long[][] delta = new long[R][D];
            int[] off = new int[R];
            boolean[] fits = new boolean[R];
            for (int r = 0; r < R; r++) {
                delta[r][r] = 1;
                for (int j = 0; j < J; j++) {
                    delta[r][R + j] = A[j][i + S * r];
                }
                fits[r] = true;
                long o = 0;
                for (int d = 0; d < D; d++) {
                    if (delta[r][d] > dims[d] - 1) {
                        fits[r] = false;    // a single job already breaks the cut
                    }
                    o += delta[r][d] * stride[d];
                }
                off[r] = (int) o;
            }

            if (i < M) {
                // Queueing station: solve (1 - sum_r L_ir u_ir) x = ser in place.
                // Every monomial of the operator raises the total class degree by
                // one, so p - off[r] is always a strictly smaller flat index and a
                // sweep in increasing flat index reads only final coefficients:
                // the sweep IS the solve.
                double[] coef = L[i];
                int[] sub = new int[D];
                for (int p = 0; p < ser.length; p++) {
                    for (int r = 0; r < R; r++) {
                        if (coef[r] == 0.0 || !fits[r]) {
                            continue;
                        }
                        boolean ok = true;
                        for (int d = 0; d < D && ok; d++) {
                            if (sub[d] < delta[r][d]) {
                                ok = false;
                            }
                        }
                        if (ok) {
                            ser[p] += coef[r] * ser[p - off[r]];
                        }
                    }
                    for (int d = 0; d < D; d++) {
                        if (++sub[d] < dims[d]) {
                            break;
                        }
                        sub[d] = 0;
                    }
                }
            } else {
                // Delay station: no n_i! coupling, so the factor is a product of
                // exponentials, one per class, each convolved in term by term.
                // There is no first-order recurrence to exploit here, which is why
                // the delay costs a factor of the population that the queueing
                // station does not.
                double[] coef = Z[i - M];
                for (int r = 0; r < R; r++) {
                    if (coef[r] == 0.0 || !fits[r] || N[r] == 0) {
                        continue;
                    }
                    double[] nxt = ser.clone();
                    double[] term = ser;
                    for (int n = 1; n <= N[r]; n++) {
                        double[] shifted = new double[ser.length];
                        int[] sub = new int[D];
                        for (int p = 0; p < ser.length; p++) {
                            boolean ok = true;
                            for (int d = 0; d < D && ok; d++) {
                                if (sub[d] < delta[r][d]) {
                                    ok = false;
                                }
                            }
                            if (ok) {
                                shifted[p] = term[p - off[r]];
                            }
                            for (int d = 0; d < D; d++) {
                                if (++sub[d] < dims[d]) {
                                    break;
                                }
                                sub[d] = 0;
                            }
                        }
                        double c = coef[r] / n;
                        boolean any = false;
                        for (int p = 0; p < shifted.length; p++) {
                            shifted[p] = c * shifted[p];
                            if (shifted[p] != 0.0) {
                                any = true;
                            }
                        }
                        term = shifted;
                        if (!any) {
                            break;
                        }
                        for (int p = 0; p < nxt.length; p++) {
                            nxt[p] += term[p];
                        }
                    }
                    ser = nxt;
                }
            }

            for (int j = 0; j < J; j++) {
                if (last[j] == i) {
                    ser = reduce(ser, dims, R + j, sense[j], (int) b[j]);
                    dims[R + j] = 1;
                }
            }
        }

        int expect = 1;
        for (int r = 0; r < R; r++) {
            expect *= N[r] + 1;
        }
        if (ser.length != expect) {
            throw new RuntimeException("pfqn_manjunath: a constraint row was never discharged; "
                    + "the elimination order is inconsistent");
        }
        // The closed network's own equalities: degree exactly N_r in every class.
        int flat = 0;
        int st = 1;
        for (int r = 0; r < R; r++) {
            flat += N[r] * st;
            st *= N[r] + 1;
        }
        return ser[flat];
    }

    private static int[] strides(int[] dims) {
        int[] stride = new int[dims.length];
        stride[0] = 1;
        for (int d = 1; d < dims.length; d++) {
            stride[d] = stride[d - 1] * dims[d - 1];
        }
        return stride;
    }

    /**
     * Create marker k, keeping the existing content at degree zero: nothing
     * multiplied in so far carries any power of it.
     */
    private static double[] expand(double[] ser, int[] dims, int k, int newdim) {
        int pre = 1;
        int post = 1;
        for (int d = 0; d < k; d++) {
            pre *= dims[d];
        }
        for (int d = k + 1; d < dims.length; d++) {
            post *= dims[d];
        }
        double[] grown = new double[pre * newdim * post];
        for (int q = 0; q < post; q++) {
            for (int p = 0; p < pre; p++) {
                grown[p + q * pre * newdim] = ser[p + q * pre];
            }
        }
        return grown;
    }

    /**
     * Discharge marker k. The multiplier (y^{b+1}-1)/(y-1) of a '&lt;=' row turns
     * its residue into the partial sum of the coefficients of degrees 0..b, and
     * the multiplier 1/y^{b+1} of an '=' row picks the coefficient of degree b.
     */
    private static double[] reduce(double[] ser, int[] dims, int k, char sense, int rhs) {
        int pre = 1;
        int post = 1;
        for (int d = 0; d < k; d++) {
            pre *= dims[d];
        }
        for (int d = k + 1; d < dims.length; d++) {
            post *= dims[d];
        }
        int dk = dims[k];
        double[] out = new double[pre * post];
        if (sense == 'E') {
            for (int q = 0; q < post; q++) {
                for (int p = 0; p < pre; p++) {
                    out[p + q * pre] = ser[p + rhs * pre + q * pre * dk];
                }
            }
        } else {
            for (int q = 0; q < post; q++) {
                for (int d = 0; d < dk; d++) {
                    for (int p = 0; p < pre; p++) {
                        out[p + q * pre] += ser[p + d * pre + q * pre * dk];
                    }
                }
            }
        }
        return out;
    }
}
