package jline.api.spn;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NonNegativeConstraint;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;

import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.NodeType;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Linear-programming bounds on the mean marking and the throughputs of a
 * stochastic timed Petri net.
 *
 * <p>The stationary chain is relaxed to a MOMENT POLYTOPE: the uniformized
 * evolution equation is written for E[X_p], E[X_p^2] and E[X_p1 X_p2], which
 * gives linear equalities among the mean marking x, the enabling probabilities
 * q and the products y(p,t) = E[X_p e_t]; behavioural and probabilistic
 * inequalities are added on top; and every reported measure is then obtained by
 * minimising and maximising its linear form over that polytope. Any stationary
 * point of the true chain satisfies every row, so the two optima BRACKET the
 * exact value whatever the polytope leaves out.</p>
 *
 * <p>This is the Petri-net sibling of the QRF bounds in SolverBA: same
 * technique, a different index space, and a LINEAR objective, so there is no
 * stationary point to escape from and the answer is a property of the model
 * alone.</p>
 *
 * <p>VARIABLES, over place levels l = 0..L-1 and modes e = 0..E-1: x(l) the
 * mean tokens, q(e) the probability that mode e is enabled, th(e) its
 * throughput, u(e) the state-equation firing counts, and y(l,e) = E[X_l e_e]
 * on the Markovian side only. u is EXISTENTIAL and is not reported: E[X] is a
 * convex combination of reachable markings, each of which is m0 + C h for some
 * nonnegative integer h, so the mean satisfies m0 + C u for some nonnegative
 * real u.</p>
 *
 * <p>LEVELS ARE (place, class) PAIRS, PLACE-MAJOR, level pp*R + k, the same
 * coordinates {@link Spn_mdd}, {@link Spn_sinvariants} and {@code Spn_conv}
 * use. MODES are (transition, mode) pairs in node order.</p>
 *
 * <p>THE TOKEN COUNTS CREATED BY A FIRING ARE DETERMINISTIC IN LINE, which
 * removes a whole branch of the reference: it allows sigma_(t,p)(n) to be
 * random and splits the covariance family into an independent case (its eq. 7)
 * and a selective one (its eq. 8). setFiringOutcome takes an integer weight, so
 * E[sigma^2] = sigma^2 and E[sigma_p1 sigma_p2] = sigma_p1 sigma_p2 hold
 * exactly and eq. (7) is the correct form. Eq. (8) has no LINE model behind it
 * and is deliberately absent.</p>
 *
 * <p>LIVENESS IS OFF BY DEFAULT, AND THAT IS DELIBERATE. The reference's two
 * liveness rows (sum_t q_t &gt;= 1 and x_p &lt;= sum_t y_(p,t)) hold only on a
 * live net, and liveness is not something this class can cheaply certify -- an
 * inhibitor arc alone is enough to deadlock a net that looks well formed. A
 * bound that silently assumed it would be wrong rather than loose on exactly
 * the models where a bound is most wanted, so the rows are opt-in.</p>
 *
 * <p>WHAT THE ROWS ARE WORTH, MEASURED. They are the whole of the lower side.
 * On the reference's own Table 2 (its Fig. 2b production line, five rate
 * vectors) assumelive reproduces its published l.b. column to four decimals --
 * 1.1653 against 1.165, 1.8288 against 1.829, 1.5814 against 1.581, 1.3592
 * against 1.359, 1.3497 against 1.350 -- while without them the Markovian lower
 * bound collapses onto the OPERATIONAL one on four of the five. The upper side
 * needs neither row and matches the published u.b.2 either way.</p>
 *
 * <p>THE APACHE SIMPLEX SOLVER ASSUMES ONLY NON-NEGATIVITY, so every variable
 * upper bound (q &lt;= 1, x &lt;= B, y &lt;= B) is materialised as an explicit
 * LEQ row, the same workaround {@code Mapqn_bnd_lr_pf} carries. The constraint
 * set is built ONCE and reused across the one solve per reported cell; only the
 * objective vector changes.</p>
 *
 * <p>Reference: Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets
 * Using Linear Programming Approach", IEEE Trans. Software Engineering 24(11),
 * 1998, 1014-1030. The constraint families are its Table 1, p. 1022; the
 * bracket statement is its Theorem 3, p. 1021.</p>
 *
 * <p>MATLAB twin: {@code spn_lpbnd.m}. Python twin: {@code api/spn/lpbnd.py}.</p>
 */
public class Spn_lpbnd {

    private Spn_lpbnd() {}

    /** Options of the relaxation. */
    public static class SpnLpOptions {
        /**
         * true (default) uses the second-moment, covariance and Little's law
         * families, which need exponential firing times; false drops them and
         * the whole y block, leaving the operational bound, which needs only a
         * mean firing time and so admits any phase-type law.
         */
        public boolean markovian = true;
        /** true adds the two liveness rows, valid only on a live net. */
        public boolean assumelive = false;
        /** Initial tokens per place level, place-major; null reads the model. */
        public double[] init = null;
        /** Slack added to the inequality sides. */
        public double tol = 0.0;
        /** Print the polytope size. */
        public boolean verbose = false;
    }

    /** One (transition, mode) pair over place-major levels. */
    public static class SpnLpMode {
        public int trans;
        public int mode;
        public double[] enab;
        public double[] inhib;
        public double[] fire;
        public double rate;
    }

    /** The brackets, each a 2 x n array with row 0 the minimum. */
    public static class SpnLpBounds {
        public int[] places;
        public String[] levelname;
        public List<SpnLpMode> modes;
        public double[][] tokens;
        public double[][] placeTput;
        public double[][] modeTput;
        public double[][] modeUtil;
        public double[] bound;
        public int nplacelevels;
        public int nclasses;
        public boolean markovian;
        public int nvars;
        public int nrows;
    }

    /** Bracket the mean tokens and the throughputs, with the default options. */
    public static SpnLpBounds spn_lpbnd(NetworkStruct sn) {
        return spn_lpbnd(sn, new SpnLpOptions());
    }

    /**
     * Bracket the mean tokens and the throughputs of a stochastic Petri net.
     *
     * @param sn a NetworkStruct holding Places and Transitions
     * @param opt the relaxation options; null takes the defaults
     * @return the brackets, per place level and per mode
     */
    public static SpnLpBounds spn_lpbnd(NetworkStruct sn, SpnLpOptions opt) {
        if (opt == null) {
            opt = new SpnLpOptions();
        }
        List<Integer> places = new ArrayList<Integer>();
        boolean hasTransition = false;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Place) {
                places.add(Integer.valueOf(i));
            } else if (sn.nodetype.get(i) == NodeType.Transition) {
                hasTransition = true;
            }
        }
        if (places.isEmpty() || !hasTransition) {
            throw new RuntimeException("spn_lpbnd: the model holds no Place or no Transition node");
        }
        int R = sn.nclasses;
        int P = places.size();
        int L = P * R;

        List<SpnLpMode> md = modes(sn, places, R, opt.markovian);
        int E = md.size();

        // Arc tables as E x L: PI consumes, SG creates, ETA inhibits (infinite
        // where there is no inhibitor arc), NET = SG - PI is the incidence.
        double[][] pi = new double[E][];
        double[][] sg = new double[E][];
        double[][] eta = new double[E][];
        double[][] net = new double[E][L];
        double[] mu = new double[E];
        for (int e = 0; e < E; e++) {
            pi[e] = md.get(e).enab;
            sg[e] = md.get(e).fire;
            eta[e] = md.get(e).inhib;
            mu[e] = md.get(e).rate;
            for (int l = 0; l < L; l++) {
                net[e][l] = sg[e][l] - pi[e][l];
            }
        }

        // Per-level a priori bounds and the conserved sums, both off the same
        // minimal-support P-invariant basis. The reference writes its "cycle
        // population" family for UNWEIGHTED cycles; Spn_sinvariants returns the
        // weighted invariants S m = V, which are equally linear and strictly
        // tighter, so those are what is emitted.
        double[] init = opt.init != null ? opt.init : initMarking(sn, places, R);
        Spn_sinvariants.SpnInvariants inv = Spn_sinvariants.spn_sinvariants(sn, init);
        long[][] S = inv.S;
        long[] V = inv.V;
        long[] m0 = inv.m0;
        double[] B = levelBounds(S, V, L);

        // ---- variable layout
        int ix = 0;
        int iq = L;
        int ith = L + E;
        int iu = L + 2 * E;
        int iy = L + 3 * E;
        int nv = opt.markovian ? L + 3 * E + L * E : L + 3 * E;

        List<LinearConstraint> rows = new ArrayList<LinearConstraint>();
        Builder bd = new Builder(nv, rows, opt.tol);

        // Apache's SimplexSolver assumes only non-negativity, so every variable
        // upper bound is an explicit row. See the class note.
        for (int e = 0; e < E; e++) {
            bd.start().add(iq + e, 1.0).le(1.0);
        }
        for (int l = 0; l < L; l++) {
            if (!Double.isInfinite(B[l])) {
                bd.start().add(ix + l, 1.0).le(B[l]);
                if (opt.markovian) {
                    for (int e = 0; e < E; e++) {
                        bd.start().add(iy + l * E + e, 1.0).le(B[l]);
                    }
                }
            }
        }

        // ---- (1) throughput: th_e = mu_e q_e
        for (int e = 0; e < E; e++) {
            bd.start().add(ith + e, 1.0).add(iq + e, -mu[e]).eq(0.0);
        }

        // ---- (2) flow balance: tokens are created at a level at the rate they
        // are consumed there. Holds for any stable net, Markovian or not.
        for (int l = 0; l < L; l++) {
            Row r = bd.start();
            for (int e = 0; e < E; e++) {
                r.add(iq + e, mu[e] * net[e][l]);
            }
            r.eq(0.0);
        }

        // ---- (3)+(4) second moment and population covariance, from the
        // stationarity of E[X_l1 X_l2] under the uniformized chain. Table 1
        // writes the q side as four sums over set intersections; since the
        // memberships are exactly "sigma > 0" and "pi > 0", those collapse to
        //     -(sigma_1 - pi_1)(sigma_2 - pi_2) = -net_1 net_2
        // per mode, which also makes the l1 == l2 case reduce to the
        // second-moment family with no separate derivation.
        if (opt.markovian) {
            for (int l1 = 0; l1 < L; l1++) {
                for (int l2 = l1; l2 < L; l2++) {
                    Row r = bd.start();
                    for (int e = 0; e < E; e++) {
                        // at l1 == l2 the two y terms address the same column
                        // and the builder sums them, which is the factor of two
                        // the reference's (6) carries
                        r.add(iy + l1 * E + e, mu[e] * net[e][l2]);
                        r.add(iy + l2 * E + e, mu[e] * net[e][l1]);
                        r.add(iq + e, mu[e] * net[e][l1] * net[e][l2]);
                    }
                    r.eq(0.0);
                }
            }
        }

        // ---- (5) liveness, only when the caller vouches for it
        if (opt.assumelive) {
            Row r = bd.start();
            for (int e = 0; e < E; e++) {
                r.add(iq + e, 1.0);
            }
            r.ge(1.0);
            if (opt.markovian) {
                for (int l = 0; l < L; l++) {
                    Row rr = bd.start().add(ix + l, 1.0);
                    for (int e = 0; e < E; e++) {
                        rr.add(iy + l * E + e, -1.0);
                    }
                    rr.le(0.0);
                }
            }
        }

        // ---- (6) conflicting transitions: a mode that consumes no more and is
        // inhibited no sooner is enabled whenever the other is
        for (int e1 = 0; e1 < E; e1++) {
            for (int e2 = 0; e2 < E; e2++) {
                if (e1 == e2) {
                    continue;
                }
                boolean dominated = true;
                for (int l = 0; l < L; l++) {
                    if (!(pi[e1][l] <= pi[e2][l] && eta[e1][l] >= eta[e2][l])) {
                        dominated = false;
                        break;
                    }
                }
                if (dominated) {
                    bd.start().add(iq + e1, 1.0).add(iq + e2, -1.0).ge(0.0);
                }
            }
        }

        // ---- (7) boundedness, per level; and (8) cycle population, as the
        // weighted invariant equalities and their y companions
        if (opt.markovian) {
            for (int l = 0; l < L; l++) {
                if (Double.isInfinite(B[l])) {
                    continue;
                }
                for (int e = 0; e < E; e++) {
                    bd.start().add(iy + l * E + e, 1.0).add(iq + e, -B[l]).le(0.0);
                    bd.start().add(ix + l, 1.0).add(iy + l * E + e, -1.0)
                            .add(iq + e, B[l]).le(B[l]);
                    if (B[l] > 0) {
                        bd.start().add(ix + l, 1.0 - 1.0 / B[l]).add(iy + l * E + e, -1.0)
                                .add(iq + e, 1.0).ge(0.0);
                    }
                }
            }
        }
        for (int i = 0; i < S.length; i++) {
            Row r = bd.start();
            for (int l = 0; l < L; l++) {
                r.add(ix + l, S[i][l]);
            }
            r.eq(V[i]);
            if (opt.markovian) {
                for (int e = 0; e < E; e++) {
                    Row rr = bd.start();
                    for (int l = 0; l < L; l++) {
                        rr.add(iy + l * E + e, S[i][l]);
                    }
                    rr.add(iq + e, -(double) V[i]).eq(0.0);
                }
            }
        }

        // ---- (9) reachable marking: the mean lies in the state-equation cone
        for (int l = 0; l < L; l++) {
            Row r = bd.start().add(ix + l, 1.0);
            for (int e = 0; e < E; e++) {
                r.add(iu + e, -net[e][l]);
            }
            r.eq(m0[l]);
        }

        // ---- (10) sample-path comparisons
        if (opt.markovian) {
            double mutot = 0;
            for (int e = 0; e < E; e++) {
                mutot += mu[e];
            }
            for (int l = 0; l < L; l++) {
                for (int e = 0; e < E; e++) {
                    bd.start().add(iy + l * E + e, 1.0).add(ix + l, -1.0).le(0.0);
                    if (pi[e][l] > 0) {
                        bd.start().add(iy + l * E + e, 1.0).add(iq + e, -pi[e][l]).ge(0.0);
                    }
                    if (!Double.isInfinite(eta[e][l])) {
                        bd.start().add(iy + l * E + e, 1.0)
                                .add(iq + e, -(eta[e][l] - 1.0)).le(0.0);
                    }
                }
                Row r = bd.start().add(ix + l, mutot);
                for (int e = 0; e < E; e++) {
                    r.add(iy + l * E + e, -mu[e]);
                }
                r.ge(0.0);
            }
            for (int e = 0; e < E; e++) {
                int ent = -1;
                boolean single = true;
                boolean inhibited = false;
                for (int l = 0; l < L; l++) {
                    if (pi[e][l] > 0) {
                        if (ent >= 0) {
                            single = false;
                        }
                        ent = l;
                    }
                    if (!Double.isInfinite(eta[e][l])) {
                        inhibited = true;
                    }
                }
                if (single && ent >= 0 && !inhibited) {
                    bd.start().add(ix + ent, 1.0).add(iy + ent * E + e, -1.0)
                            .le(pi[e][ent] - 1.0);
                }
            }
        }

        // ---- (11) enabling bounds, from Chernoff's inequality on the marking
        for (int e = 0; e < E; e++) {
            List<Integer> ent = new ArrayList<Integer>();
            List<Integer> inh = new ArrayList<Integer>();
            for (int l = 0; l < L; l++) {
                if (pi[e][l] > 0) {
                    ent.add(Integer.valueOf(l));
                }
                if (!Double.isInfinite(eta[e][l])) {
                    inh.add(Integer.valueOf(l));
                }
            }
            int d = ent.size() + inh.size();
            if (d == 0) {
                continue;
            }
            boolean entBounded = !ent.isEmpty();
            for (int j = 0; j < ent.size(); j++) {
                if (Double.isInfinite(B[ent.get(j).intValue()])) {
                    entBounded = false;
                }
            }
            if (entBounded) {
                Row r = bd.start().add(iq + e, 1.0);
                double rhs = 1.0;
                boolean ok = true;
                for (int j = 0; j < ent.size() && ok; j++) {
                    int l = ent.get(j).intValue();
                    double den = B[l] - pi[e][l] + 1.0;
                    if (den <= 0) {
                        ok = false;
                        break;
                    }
                    r.add(ix + l, -1.0 / den);
                    rhs -= B[l] / den;
                }
                if (ok) {
                    for (int j = 0; j < inh.size(); j++) {
                        int l = inh.get(j).intValue();
                        r.add(ix + l, 1.0 / eta[e][l]);
                    }
                    r.ge(rhs);
                } else {
                    r.discard();
                }
            }
            // Upper side. Every term of the sum over input levels carries a
            // "min" operator, and Table 1's convention is that either operand
            // may be taken; each choice is a valid row and the whole set is the
            // tightest linear relaxation, so all of them are emitted while the
            // count stays small.
            boolean ok = true;
            for (int j = 0; j < inh.size(); j++) {
                int l = inh.get(j).intValue();
                if (Double.isInfinite(B[l]) || B[l] - eta[e][l] + 1.0 <= 0) {
                    ok = false;
                }
            }
            if (!ok) {
                continue;
            }
            int nc = ent.size();
            int[] combos = nc <= 4 ? new int[1 << nc] : new int[] {0, (1 << nc) - 1};
            if (nc <= 4) {
                for (int c = 0; c < combos.length; c++) {
                    combos[c] = c;
                }
            }
            for (int ci = 0; ci < combos.length; ci++) {
                int c = combos[ci];
                Row r = bd.start().add(iq + e, d);
                double rhs = 0.0;
                for (int j = 0; j < inh.size(); j++) {
                    int l = inh.get(j).intValue();
                    double den = B[l] - eta[e][l] + 1.0;
                    rhs += B[l] / den;
                    r.add(ix + l, 1.0 / den);
                }
                for (int j = 0; j < nc; j++) {
                    int l = ent.get(j).intValue();
                    if (((c >> j) & 1) == 0) {
                        r.add(ix + l, -1.0 / pi[e][l]);
                    } else {
                        rhs += 1.0;
                    }
                }
                r.le(rhs);
            }
        }

        // ---- (12) Little's law at each level: the mean sojourn time of a
        // token is at least the mean minimum firing time of the modes that can
        // remove it
        if (opt.markovian) {
            for (int l = 0; l < L; l++) {
                double out = 0;
                for (int e = 0; e < E; e++) {
                    if (pi[e][l] > 0) {
                        out += mu[e];
                    }
                }
                if (out <= 0) {
                    continue;
                }
                Row r = bd.start().add(ix + l, out);
                for (int e = 0; e < E; e++) {
                    r.add(iq + e, -mu[e] * sg[e][l]);
                }
                r.ge(0.0);
            }
        }

        if (opt.verbose) {
            System.out.printf("%nSPN -> LP: %d place levels, %d modes, %d variables, %d rows%n",
                    L, E, nv, rows.size());
        }

        LinearConstraintSet cs = new LinearConstraintSet(rows);
        SpnLpBounds out = new SpnLpBounds();
        out.tokens = new double[2][L];
        out.placeTput = new double[2][L];
        for (int l = 0; l < L; l++) {
            double[] c = new double[nv];
            c[ix + l] = 1.0;
            out.tokens[0][l] = solve(cs, c, nv, true);
            out.tokens[1][l] = solve(cs, c, nv, false);
            double[] ct = new double[nv];
            boolean any = false;
            for (int e = 0; e < E; e++) {
                if (pi[e][l] > 0) {
                    ct[iq + e] = mu[e] * pi[e][l];
                    any = true;
                }
            }
            out.placeTput[0][l] = any ? solve(cs, ct, nv, true) : 0.0;
            out.placeTput[1][l] = any ? solve(cs, ct, nv, false) : 0.0;
        }
        out.modeTput = new double[2][E];
        out.modeUtil = new double[2][E];
        for (int e = 0; e < E; e++) {
            double[] c = new double[nv];
            c[ith + e] = 1.0;
            out.modeTput[0][e] = solve(cs, c, nv, true);
            out.modeTput[1][e] = solve(cs, c, nv, false);
            double[] cq = new double[nv];
            cq[iq + e] = 1.0;
            out.modeUtil[0][e] = solve(cs, cq, nv, true);
            out.modeUtil[1][e] = solve(cs, cq, nv, false);
        }

        out.places = new int[P];
        out.levelname = new String[L];
        for (int pp = 0; pp < P; pp++) {
            out.places[pp] = places.get(pp).intValue();
            for (int k = 0; k < R; k++) {
                out.levelname[pp * R + k] = sn.nodenames.get(out.places[pp]) + "."
                        + sn.classnames.get(k);
            }
        }
        out.modes = md;
        out.bound = B;
        out.nplacelevels = L;
        out.nclasses = R;
        out.markovian = opt.markovian;
        out.nvars = nv;
        out.nrows = rows.size();
        return out;
    }

    // -----------------------------------------------------------------------
    /**
     * Row accumulator. Repeated indices are SUMMED rather than overwritten,
     * matching the MATLAB sparse() and Python coo-to-csr assemblies, which is
     * what supplies the factor of two the l1 == l2 covariance row carries.
     */
    private static final class Builder {
        private final int nv;
        private final List<LinearConstraint> rows;
        private final double tol;

        Builder(int nv, List<LinearConstraint> rows, double tol) {
            this.nv = nv;
            this.rows = rows;
            this.tol = tol;
        }

        Row start() {
            return new Row(this);
        }
    }

    /** One row under construction. */
    private static final class Row {
        private final Builder bd;
        private double[] coeff;

        Row(Builder bd) {
            this.bd = bd;
            this.coeff = new double[bd.nv];
        }

        Row add(int j, double v) {
            if (coeff != null && v != 0) {
                coeff[j] += v;
            }
            return this;
        }

        /** Abandon the row without emitting it. */
        void discard() {
            coeff = null;
        }

        void eq(double rhs) {
            emit(Relationship.EQ, rhs);
        }

        void le(double rhs) {
            emit(Relationship.LEQ, rhs + bd.tol);
        }

        void ge(double rhs) {
            emit(Relationship.GEQ, rhs - bd.tol);
        }

        private void emit(Relationship rel, double rhs) {
            if (coeff == null) {
                return;
            }
            boolean any = false;
            for (int j = 0; j < coeff.length && !any; j++) {
                any = coeff[j] != 0;
            }
            if (any) {
                bd.rows.add(new LinearConstraint(coeff, rel, rhs));
            }
            coeff = null;
        }
    }

    /**
     * One LP. An infeasible or unbounded program returns NaN rather than
     * throwing, matching the MATLAB reference, whose own note is that the exit
     * flag is the wrong success predicate and the finiteness of the answer is
     * the right one.
     */
    private static double solve(LinearConstraintSet cs, double[] c, int nv, boolean minimize) {
        try {
            SimplexSolver solver = new SimplexSolver();
            PointValuePair sol = solver.optimize(
                    new LinearObjectiveFunction(c, 0.0), cs,
                    new NonNegativeConstraint(true),
                    minimize ? GoalType.MINIMIZE : GoalType.MAXIMIZE);
            if (sol == null) {
                return Double.NaN;
            }
            double v = sol.getValue();
            return Double.isFinite(v) ? v : Double.NaN;
        } catch (RuntimeException ex) {
            return Double.NaN;
        }
    }

    /**
     * The (transition, mode) table over place-major levels. {@link Spn_mdd}
     * builds the same table, but only as a step of reachable-set construction,
     * which is the cost this bound exists to avoid.
     */
    private static List<SpnLpMode> modes(NetworkStruct sn, List<Integer> places, int R,
                                         boolean markovian) {
        List<SpnLpMode> md = new ArrayList<SpnLpMode>();
        for (int ind = 0; ind < sn.nodetype.size(); ind++) {
            if (sn.nodetype.get(ind) != NodeType.Transition) {
                continue;
            }
            Node node = sn.nodes.get(ind);
            NodeParam param = sn.nodeparam.get(node);
            if (!(param instanceof TransitionNodeParam)) {
                continue;
            }
            TransitionNodeParam tp = (TransitionNodeParam) param;
            for (int m = 0; m < tp.nmodes; m++) {
                if (tp.timing != null && tp.timing.size() > m
                        && tp.timing.get(m) == TimingStrategy.IMMEDIATE) {
                    throw new RuntimeException("spn_lpbnd: mode " + (m + 1) + " of node "
                            + (ind + 1) + " is IMMEDIATE; the moment relaxation is written for a "
                            + "net whose transitions all have finite rates, so vanishing states "
                            + "must be eliminated first");
                }
                if (tp.firingdep != null && tp.firingdep.size() > m
                        && tp.firingdep.get(m) != null) {
                    throw new RuntimeException("spn_lpbnd: mode " + (m + 1) + " of node "
                            + (ind + 1) + " has a marking-dependent firing rate; the "
                            + "uniformization step needs one rate per mode");
                }
                double srv = (tp.nmodeservers != null && tp.nmodeservers.getNumElements() > m)
                        ? tp.nmodeservers.get(m) : 1.0;
                if (srv != 1.0) {
                    throw new RuntimeException("spn_lpbnd: mode " + (m + 1) + " of node "
                            + (ind + 1) + " has " + srv + " servers; the relaxation is derived "
                            + "under single-server semantics, where the firing rate is mu*q. Its "
                            + "infinite-server form needs the K-fold transition expansion of the "
                            + "reference's Section 7, which is not implemented");
                }
                // A PHASE-TYPE FIRING LAW IS WHERE THE TWO VARIANTS PART. The
                // mean of a (D0,D1) pair is pie*(-D0)^-1*1 and is all the
                // operational bound needs; the Markovian one needs the marking
                // alone to be the state, which a multi-phase mode breaks. A law
                // that is not phase-type at all reaches sn as its own parameter
                // list, with no mean recoverable without a per-distribution
                // table, so BOTH variants refuse it.
                MatrixCell proc = firingProc(tp, m);
                if (proc == null || proc.size() < 2) {
                    throw new RuntimeException("spn_lpbnd: mode " + (m + 1) + " of node "
                            + (ind + 1) + " has a firing law that is not phase-type; sn carries "
                            + "its own parameters rather than a (D0,D1) pair, so neither the "
                            + "Markovian nor the operational bound can read a mean firing rate "
                            + "from it. Use a phase-type law, or an exact solver");
                }
                Matrix d0 = proc.get(0);
                Matrix d1 = proc.get(1);
                int nph = d0.getNumRows();
                if (markovian && nph > 1) {
                    throw new RuntimeException("spn_lpbnd: mode " + (m + 1) + " of node "
                            + (ind + 1) + " has a phase-type firing time; the relaxation is "
                            + "written over the marking alone, and a phase-type mode needs the "
                            + "state-machine expansion of the reference's Section 7, which is not "
                            + "implemented. Use the operational bound, which needs only the mean");
                }
                double rate;
                if (nph == 1) {
                    rate = d1.get(0, 0);
                } else {
                    rate = 1.0 / phMean(d0, firingPie(tp, m, nph));
                }
                if (!(rate > 0) || !Double.isFinite(rate)) {
                    throw new RuntimeException("spn_lpbnd: mode " + (m + 1) + " of node "
                            + (ind + 1) + " has mean firing rate " + rate + "; a bound needs a "
                            + "finite positive one");
                }
                SpnLpMode e = new SpnLpMode();
                e.trans = ind;
                e.mode = m;
                e.enab = arcvec(tp.enabling, m, places, R, 0.0);
                e.inhib = arcvec(tp.inhibiting, m, places, R, Double.POSITIVE_INFINITY);
                e.fire = arcvec(tp.firing, m, places, R, 0.0);
                e.rate = rate;
                md.add(e);
            }
        }
        if (md.isEmpty()) {
            throw new RuntimeException("spn_lpbnd: the net has no firing mode");
        }
        return md;
    }

    /** Mean time to absorption of a phase-type law, pie*(-D0)^-1*1. */
    private static double phMean(Matrix d0, double[] pie) {
        int n = d0.getNumRows();
        double[][] a = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a[i][j] = -d0.get(i, j);
            }
            a[i][n] = 1.0;
        }
        // Gaussian elimination with partial pivoting on (-D0) z = 1
        for (int c = 0; c < n; c++) {
            int piv = c;
            for (int r = c + 1; r < n; r++) {
                if (Math.abs(a[r][c]) > Math.abs(a[piv][c])) {
                    piv = r;
                }
            }
            double[] t = a[c];
            a[c] = a[piv];
            a[piv] = t;
            if (a[c][c] == 0) {
                return Double.NaN;
            }
            for (int r = 0; r < n; r++) {
                if (r == c) {
                    continue;
                }
                double f = a[r][c] / a[c][c];
                for (int j = c; j <= n; j++) {
                    a[r][j] -= f * a[c][j];
                }
            }
        }
        double mean = 0;
        for (int i = 0; i < n; i++) {
            mean += pie[i] * a[i][n] / a[i][i];
        }
        return mean;
    }

    private static MatrixCell firingProc(TransitionNodeParam tp, int m) {
        if (tp.firingproc == null) {
            return null;
        }
        int i = 0;
        for (Map.Entry<jline.lang.Mode, MatrixCell> en : tp.firingproc.entrySet()) {
            if (i == m) {
                return en.getValue();
            }
            i++;
        }
        return null;
    }

    private static double[] firingPie(TransitionNodeParam tp, int m, int nph) {
        double[] pie = new double[nph];
        Matrix pv = null;
        if (tp.firingpie != null) {
            int i = 0;
            for (Map.Entry<jline.lang.Mode, Matrix> en : tp.firingpie.entrySet()) {
                if (i == m) {
                    pv = en.getValue();
                    break;
                }
                i++;
            }
        }
        if (pv == null || pv.getNumElements() == 0) {
            pie[0] = 1.0;
            return pie;
        }
        double s = 0;
        for (int a = 0; a < nph && a < pv.getNumElements(); a++) {
            pie[a] = pv.get(a);
            s += pie[a];
        }
        if (s <= 0) {
            pie[0] = 1.0;
            return pie;
        }
        for (int a = 0; a < nph; a++) {
            pie[a] /= s;
        }
        return pie;
    }

    /**
     * (nnodes x nclasses) arc matrix -&gt; place-major length-(P*R) level
     * vector. The row index is a NODE index; fillval is 0 for enabling and
     * firing and infinite for inhibiting, where infinity means "no arc" and 0
     * would mean "inhibited always".
     */
    private static double[] arcvec(List<Matrix> mats, int m, List<Integer> places, int R,
                                   double fillval) {
        int P = places.size();
        double[] v = new double[P * R];
        for (int i = 0; i < v.length; i++) {
            v[i] = fillval;
        }
        if (mats == null || mats.size() <= m || mats.get(m) == null) {
            return v;
        }
        Matrix mat = mats.get(m);
        for (int pp = 0; pp < P; pp++) {
            for (int k = 0; k < R; k++) {
                double x = mat.get(places.get(pp).intValue(), k);
                if (!Double.isInfinite(fillval) && fillval == 0) {
                    x = Math.max(0.0, x);
                }
                v[pp * R + k] = x;
            }
        }
        return v;
    }

    /** Declared tokens per (place, class), or null when nothing is set. */
    private static double[] initMarking(NetworkStruct sn, List<Integer> places, int R) {
        int P = places.size();
        double[] init = new double[P * R];
        boolean anyset = false;
        for (int pp = 0; pp < P; pp++) {
            Node node = sn.nodes.get(places.get(pp).intValue());
            if (!(node instanceof StatefulNode)) {
                continue;
            }
            Matrix st = ((StatefulNode) node).getState();
            if (st == null || st.getNumElements() == 0) {
                continue;
            }
            for (int k = 0; k < R && k < st.getNumElements(); k++) {
                init[pp * R + k] = st.get(k);
                anyset = anyset || st.get(k) > 0;
            }
        }
        // null lets Spn_sinvariants take the reference-station default
        return anyset ? init : null;
    }

    /**
     * Tightest a priori bound per level: an invariant with weight w on level l
     * and conserved value V caps that level at floor(V/w).
     */
    private static double[] levelBounds(long[][] S, long[] V, int L) {
        double[] B = new double[L];
        for (int l = 0; l < L; l++) {
            B[l] = Double.POSITIVE_INFINITY;
        }
        for (int i = 0; i < S.length; i++) {
            for (int l = 0; l < L; l++) {
                if (S[i][l] > 0) {
                    B[l] = Math.min(B[l], Math.floor((double) V[i] / (double) S[i][l]));
                }
            }
        }
        return B;
    }
}
