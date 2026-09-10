package jline.api.spn;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.api.mc.Ctmc_solve;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

/**
 * Product form of a stochastic Petri net: decide whether one exists and derive
 * the per-level factors g_l that {@link jline.api.mdd.Mdd_rec} and
 * {@link Spn_metrics} take as input.
 *
 * <p>THIS IS THE PART THE MDD-REC PAPER DECLARES OUT OF SCOPE (FGCS Sec. 3.2).
 * Every other function in api/spn receives the g_l already formed; this one
 * derives them from the net, which is what lets a solver reach them.</p>
 *
 * <p><b>The theory, in one paragraph.</b> Write I(t), O(t) for the input and
 * output vectors of mode t and lambda_t for its rate constant. Henderson-Taylor
 * and Coleman-Henderson-Taylor show that a net whose firing rate has the form
 * r_t(m) = lambda_t psi(m - I(t)) / psi(m) for m &gt;= I(t) has invariant measure
 * pi(m) = psi(m) prod_l y_l^(m_l) whenever the positive vector y satisfies
 * COMPLEX BALANCE: reading the distinct vectors appearing as some I(t) or O(t)
 * as the COMPLEXES of the net, the flow into every complex equals the flow out
 * of it, sum over t with O(t)=v of lambda_t y^(I(t)) equals
 * (sum over t with I(t)=v of lambda_t) times y^v.</p>
 *
 * <p>Two choices of psi are realisable in LINE's own rate law, and they are the
 * two this class tests for. <b>psi = 1</b> gives r_t = lambda_t, the rate of a
 * SINGLE-SERVER mode, and pi(m) = prod_l y_l^(m_l), so g_l(k) = y_l^k.
 * <b>psi = prod_l 1/m_l!</b> gives r_t = lambda_t prod_l m_l!/(m_l-I_l)!, MASS
 * ACTION, reached through Transition.setFiringRateDependence or, for a mode
 * drawing one token from one place, by infinite-server semantics; then
 * pi(m) = prod_l y_l^(m_l)/m_l!, so g_l(k) = y_l^k/k!.</p>
 *
 * <p>Which one holds is not guessed from the model API: the effective rate LINE
 * would use, lambda_t min(enabling degree, servers) g(m), is EVALUATED at every
 * reachable marking and compared against both laws. A net that matches neither
 * under one common psi is refused by name, never approximated.</p>
 *
 * <p><b>Solving for y.</b> Complex balance reads A_lambda Psi(y) = 0 with
 * A_lambda the Laplacian of the weighted digraph on complexes and
 * Psi(y)_v = y^v. That Laplacian is the TRANSPOSED GENERATOR of a Markov chain
 * that hops from complex to complex at the rate of the mode joining them, so its
 * kernel on one linkage class is that chain's stationary distribution and
 * {@link Ctmc_solve} returns it -- strictly positive exactly when the class is
 * strongly connected, which is weak reversibility. With that positive vector
 * kappa in hand y follows from the LINEAR system in x = log y,
 * (v - v0) x = log kappa_v - log kappa_v0 for v, v0 in the same linkage class,
 * solved in minimum norm. Feinberg's Deficiency Zero Theorem says this system is
 * consistent for every choice of rate constants when the net is weakly
 * reversible and its deficiency c - l - s is zero, which is why those two
 * numbers are reported; but consistency is CHECKED rather than assumed, so a net
 * of positive deficiency whose particular rates still admit a complex-balanced
 * point is accepted on the evidence.</p>
 *
 * <p><b>The gauge, and why the minimum-norm solution is the canonical one.</b>
 * Complex balance fixes y only up to y -&gt; y .* exp(u) for any u orthogonal to
 * the stoichiometric subspace S. Such a shift multiplies pi(m) by exp(u'm),
 * which is CONSTANT on one compatibility class, so every reported measure is
 * invariant under it -- but the normalising constant G itself is not, it scales
 * by that constant. A gauge must therefore be FIXED, or the four codebases would
 * report four different G on the same net. The one fixed here is x in the row
 * space of the constraint matrix, i.e. the minimum-norm solution, reached in a
 * form that is unique whichever least-squares primitive a codebase carries:
 * solve (rows rows^T) w = rhs and set x = rows^T w. Any two solutions w of that
 * system give the SAME rows^T w, so the answer does not depend on how the
 * rank-deficient solve breaks its tie.</p>
 *
 * <p>References: J. L. Coleman, W. Henderson, P. G. Taylor, "Product form
 * equilibrium distributions and a convolution algorithm for stochastic Petri
 * nets", Performance Evaluation 26(3), 1996. M. Feinberg, "Complex balancing in
 * general kinetic systems", Arch. Rational Mech. Anal. 49, 1972.
 * D. F. Anderson, G. Craciun, T. G. Kurtz, "Product-form stationary
 * distributions for deficiency zero chemical reaction networks", Bull. Math.
 * Biol. 72, 2010.</p>
 *
 * @see jline.api.mdd.Mdd_rec
 * @see Spn_metrics
 * @see Spn_mdd
 */
public class Spn_pf {

    private Spn_pf() {
    }

    /** Options of the product-form derivation. */
    public static class SpnPfOptions {
        /** Per-place-level token bound, passed to spn_mdd; null infers it. */
        public double[] bound;
        /** Relative tolerance of the rate-law and complex-balance checks. */
        public double tol = 1e-9;
        public boolean verbose = false;
    }

    /** The product form, and the certificate that it is one. */
    public static class SpnPfResult {
        /** g[l][k] = g_l(k), ready for Mdd_rec. */
        public double[][] g;
        /** Positive vector solving complex balance. */
        public double[] y;
        /** "geometric" or "massaction", the psi that was found. */
        public String kind;
        /** The distinct complexes, one row each. */
        public double[][] complexes;
        public int deficiency;
        public int linkage;
        public int srank;
        public boolean weaklyReversible;
        /** Relative complex-balance residual at y. */
        public double residual;
        public Spn_mdd.SpnResult spn;
        /** Filled by the caller that runs Spn_metrics on this product form. */
        public Spn_metrics.SpnMetricsResult metrics;
    }

    /** Derive the product form with the default options. */
    public static SpnPfResult spn_pf(Network model) {
        return spn_pf(model, new SpnPfOptions());
    }

    /**
     * Derive the product form of a stochastic Petri net.
     *
     * @param model   a Network holding Places and Transitions
     * @param options tolerances and bounds; null takes the defaults
     * @return the per-level factors and the certificate
     */
    public static SpnPfResult spn_pf(Network model, SpnPfOptions options) {
        if (options == null) {
            options = new SpnPfOptions();
        }
        Spn_mdd.SpnOptions mddopt = new Spn_mdd.SpnOptions();
        mddopt.descriptor = false;
        mddopt.bound = options.bound;
        Spn_mdd.SpnResult spn = Spn_mdd.spn_mdd(model, mddopt);
        Spn_mdd.SpnInfo info = spn.info;

        int L = info.nplacelevels;
        List<Spn_mdd.SpnMode> md = info.modes;
        int E = md.size();
        if (E == 0) {
            throw new RuntimeException("spn_pf: the net has no timed mode");
        }

        // ---- a queueing place holds an embedded server, not a token container
        NetworkStruct sn = model.getStruct();
        for (int pp = 0; pp < info.places.length; pp++) {
            int ist = (int) sn.nodeToStation.get(info.places[pp]);   // linear: the row/column
                                                          // orientation is not fixed
            if (ist >= 0 && sn.sched.get(sn.stations.get(ist)) != SchedStrategy.INF) {
                throw new RuntimeException("spn_pf: place " + info.placenames[pp]
                        + " is a QUEUEING place (scheduling "
                        + sn.sched.get(sn.stations.get(ist))
                        + "): its embedded service is state that the marking does not carry, so "
                        + "the net is not the token-container Petri net this product form is "
                        + "written for");
            }
        }

        // ---- the rate constants and the structural vectors
        double[] lambda = new double[E];
        double[][] Iv = new double[E][L];
        double[][] Ov = new double[E][L];
        for (int e = 0; e < E; e++) {
            Spn_mdd.SpnMode mde = md.get(e);
            lambda[e] = mde.D1[0][0];
            for (int l = 0; l < L; l++) {
                Iv[e][l] = mde.enab[l];
                Ov[e][l] = mde.fire[l];
                if (!Double.isInfinite(mde.inhib[l])) {
                    throw new RuntimeException("spn_pf: mode " + (mde.mode + 1) + " of node "
                            + (mde.trans + 1) + " has an inhibitor arc. An inhibitor zeroes the "
                            + "firing rate on markings that still satisfy m >= I(t), so the rate "
                            + "is not lambda*psi(m-I)/psi(m) on any psi and the net has no "
                            + "product form of this kind");
                }
            }
            boolean consumes = false;
            for (int l = 0; l < L; l++) {
                if (mde.enab[l] > 0) {
                    consumes = true;
                    break;
                }
            }
            if (mde.srv != 1 && !consumes) {
                throw new RuntimeException("spn_pf: mode " + (mde.mode + 1) + " of node "
                        + (mde.trans + 1) + " has " + mde.srv + " servers but consumes from no "
                        + "place, so its enabling degree is unbounded and its firing rate "
                        + "undefined");
            }
            if (!(lambda[e] > 0)) {
                throw new RuntimeException("spn_pf: mode " + (mde.mode + 1) + " of node "
                        + (mde.trans + 1) + " has a non-positive firing rate");
            }
        }

        // ---- which psi does LINE's own rate law follow on this net?
        int[][] states = info.mdd.enumerate();
        String kind = rateLaw(states, md, lambda, info, options.tol);

        // ---- complexes and the weighted digraph on them
        List<double[]> clist = new ArrayList<double[]>();
        Map<String, Integer> seen = new HashMap<String, Integer>();
        int[] src = new int[E];
        int[] dst = new int[E];
        for (int e = 0; e < E; e++) {
            src[e] = complexIndex(Iv[e], clist, seen);
        }
        for (int e = 0; e < E; e++) {
            dst[e] = complexIndex(Ov[e], clist, seen);
        }
        int c = clist.size();
        double[][] C = new double[c][];
        for (int v = 0; v < c; v++) {
            C[v] = clist.get(v);
        }

        // ---- Laplacian of the complex graph: A[j][i] is the rate of the arc i -> j
        double[][] A = new double[c][c];
        for (int e = 0; e < E; e++) {
            if (src[e] == dst[e]) {
                continue;                                  // a mode that moves nothing
            }
            A[dst[e]][src[e]] += lambda[e];
            A[src[e]][src[e]] -= lambda[e];
        }

        // ---- linkage classes, weak reversibility, deficiency
        int[] lclass = new int[c];
        int nlink = linkage(c, src, dst, lclass);
        boolean wr = weaklyReversible(c, src, dst, lclass, nlink);
        Matrix netm = new Matrix(E, L);
        for (int e = 0; e < E; e++) {
            for (int l = 0; l < L; l++) {
                netm.set(e, l, Ov[e][l] - Iv[e][l]);
            }
        }
        int srank = netm.rank();
        int deficiency = c - nlink - srank;

        // ---- kappa: the positive balance flow on each linkage class
        double[] kappa = new double[c];
        for (int b = 0; b < nlink; b++) {
            List<Integer> idx = new ArrayList<Integer>();
            for (int v = 0; v < c; v++) {
                if (lclass[v] == b) {
                    idx.add(Integer.valueOf(v));
                }
            }
            if (idx.size() == 1) {
                kappa[idx.get(0).intValue()] = 1.0;
                continue;
            }
            // The block is the transposed generator of the complex-hopping chain,
            // so its kernel is that chain's stationary law.
            int nb = idx.size();
            Matrix Qb = new Matrix(nb, nb);
            for (int i = 0; i < nb; i++) {
                for (int j = 0; j < nb; j++) {
                    Qb.set(i, j, A[idx.get(j).intValue()][idx.get(i).intValue()]);
                }
            }
            Matrix pb = Ctmc_solve.ctmc_solve(Qb);
            double mx = 0;
            for (int i = 0; i < nb; i++) {
                mx = Math.max(mx, pb.get(i));
            }
            for (int i = 0; i < nb; i++) {
                double v = pb.get(i);
                if (!(v > 0)) {
                    throw new RuntimeException("spn_pf: linkage class " + (b + 1) + " of the "
                            + "complex graph carries no flow through complex "
                            + (idx.get(i).intValue() + 1) + ", so the net admits no positive "
                            + "complex-balanced point. A weakly reversible net has a strictly "
                            + "positive balance flow on every linkage class; this one is "
                            + wrText(wr));
                }
                kappa[idx.get(i).intValue()] = v / mx;
            }
        }

        // ---- x = log y from the linear system on each linkage class
        List<double[]> rows = new ArrayList<double[]>();
        List<Double> rhs = new ArrayList<Double>();
        for (int b = 0; b < nlink; b++) {
            int v0 = -1;
            for (int v = 0; v < c; v++) {
                if (lclass[v] == b) {
                    if (v0 < 0) {
                        v0 = v;
                    } else {
                        double[] row = new double[L];
                        for (int l = 0; l < L; l++) {
                            row[l] = C[v][l] - C[v0][l];
                        }
                        rows.add(row);
                        rhs.add(Double.valueOf(Math.log(kappa[v]) - Math.log(kappa[v0])));
                    }
                }
            }
        }
        double[] x = new double[L];
        if (!rows.isEmpty()) {
            Matrix Ar = new Matrix(rows.size(), L);
            Matrix br = new Matrix(rows.size(), 1);
            for (int i = 0; i < rows.size(); i++) {
                for (int l = 0; l < L; l++) {
                    Ar.set(i, l, rows.get(i)[l]);
                }
                br.set(i, 0, rhs.get(i).doubleValue());
            }
            // minimum norm through the row space; see the gauge note above
            Matrix AAt = Ar.mult(Ar.transpose());
            Matrix xm = Ar.transpose().mult(AAt.pinv().mult(br));
            for (int l = 0; l < L; l++) {
                x[l] = xm.get(l, 0);
            }
            double res = 0;
            double scale = 1;
            for (int i = 0; i < rows.size(); i++) {
                double s = 0;
                for (int l = 0; l < L; l++) {
                    s += rows.get(i)[l] * x[l];
                }
                res = Math.max(res, Math.abs(s - rhs.get(i).doubleValue()));
                scale = Math.max(scale, Math.abs(rhs.get(i).doubleValue()));
            }
            if (res > options.tol * scale) {
                throw new RuntimeException(String.format("spn_pf: the complex-balance equations "
                        + "are inconsistent (residual %.3e): this net has no product form of the "
                        + "tested kind at these rates. Its deficiency is %d and it is %s; the "
                        + "Deficiency Zero Theorem guarantees a solution only at deficiency 0 "
                        + "with weak reversibility", res, deficiency, wrText(wr)));
            }
        }
        double[] y = new double[L];
        for (int l = 0; l < L; l++) {
            y[l] = Math.exp(x[l]);
        }

        // ---- verify complex balance itself, which is what makes pi stationary
        double[] psi = new double[c];
        for (int v = 0; v < c; v++) {
            double p = 1;
            for (int l = 0; l < L; l++) {
                p *= Math.pow(y[l], C[v][l]);
            }
            psi[v] = p;
        }
        double resb = 0;
        double scaleb = Double.MIN_VALUE;
        for (int v = 0; v < c; v++) {
            double s = 0;
            double sa = 0;
            for (int u = 0; u < c; u++) {
                s += A[v][u] * psi[u];
                sa += Math.abs(A[v][u]) * psi[u];
            }
            resb = Math.max(resb, Math.abs(s));
            scaleb = Math.max(scaleb, sa);
        }
        if (resb > options.tol * scaleb) {
            throw new RuntimeException(String.format("spn_pf: complex balance fails at the "
                    + "computed point (relative residual %.3e), so the product form would not be "
                    + "stationary", resb / scaleb));
        }

        // ---- the per-level factors, tabulated over the reachable domain
        double[][] g = new double[L][];
        for (int l = 0; l < L; l++) {
            int d = spn.mdds.domain[l];
            g[l] = new double[d];
            double fact = 1;
            for (int k = 0; k < d; k++) {
                if (k > 0) {
                    fact *= k;
                }
                double v = Math.pow(y[l], k);
                g[l][k] = "massaction".equals(kind) ? v / fact : v;
            }
        }

        SpnPfResult out = new SpnPfResult();
        out.g = g;
        out.y = y;
        out.kind = kind;
        out.complexes = C;
        out.deficiency = deficiency;
        out.linkage = nlink;
        out.srank = srank;
        out.weaklyReversible = wr;
        out.residual = resb / scaleb;
        out.spn = spn;
        if (options.verbose) {
            System.out.format("%nSPN product form: %s, %d complexes, %d linkage classes, rank %d, "
                    + "deficiency %d, %s%n", kind, c, nlink, srank, deficiency, wrText(wr));
        }
        return out;
    }

    /**
     * Which psi reproduces the rate LINE would actually use, everywhere. Both
     * candidates are tried on every mode at every reachable marking; a net
     * matching neither, or matching different ones on different modes, has no
     * product form of this family and is refused by name.
     */
    private static String rateLaw(int[][] states, List<Spn_mdd.SpnMode> md, double[] lambda,
                                  Spn_mdd.SpnInfo info, double tol) {
        int E = md.size();
        int L = info.nplacelevels;
        boolean okgeo = true;
        boolean okma = true;
        double[] m = new double[L];
        for (int s = 0; s < states.length; s++) {
            for (int l = 0; l < L; l++) {
                m[l] = states[s][l];
            }
            Matrix marc = null;
            for (int e = 0; e < E; e++) {
                Spn_mdd.SpnMode mde = md.get(e);
                boolean enabled = true;
                for (int l = 0; l < L; l++) {
                    if (m[l] < mde.enab[l]) {
                        enabled = false;
                        break;
                    }
                }
                if (!enabled) {
                    continue;
                }
                double actual = lambda[e] * servers(m, mde, L);
                if (mde.dep != null) {
                    if (marc == null) {
                        marc = arcMatrix(m, info);
                    }
                    actual = actual * mde.dep.apply(marc).doubleValue();
                }
                double ma = lambda[e] * massAction(m, mde.enab, L);
                okgeo = okgeo && close(actual, lambda[e], tol);
                okma = okma && close(actual, ma, tol);
                if (!okgeo && !okma) {
                    throw new RuntimeException(String.format("spn_pf: mode %d of node %d fires at "
                            + "rate %g in a reachable marking, which is neither its rate constant "
                            + "(single-server, psi = 1) nor its mass-action rate %g "
                            + "(psi = prod 1/m!). LINE's rate law on this mode is "
                            + "lambda*min(enabling degree, servers)*g(m), and no psi puts that in "
                            + "the form lambda*psi(m-I)/psi(m)",
                            mde.mode + 1, mde.trans + 1, actual, ma));
                }
            }
        }
        return okgeo ? "geometric" : "massaction";
    }

    /** min(enabling degree, servers): the sets of tokens firing at once. */
    private static double servers(double[] m, Spn_mdd.SpnMode mde, int L) {
        double deg = Double.POSITIVE_INFINITY;
        for (int l = 0; l < L; l++) {
            if (mde.enab[l] > 0) {
                deg = Math.min(deg, Math.floor(m[l] / mde.enab[l]));
            }
        }
        if (Double.isInfinite(deg)) {
            deg = 1;                                       // consumes nothing: always one set
        }
        return Math.min(deg, mde.srv);
    }

    /** prod_l m_l!/(m_l - I_l)!, the ordered ways to pick the input tokens. */
    private static double massAction(double[] m, double[] enab, int L) {
        double r = 1;
        for (int l = 0; l < L; l++) {
            for (int j = 0; j < (int) enab[l]; j++) {
                r *= (m[l] - j);
            }
        }
        return r;
    }

    private static boolean close(double a, double b, double tol) {
        return Math.abs(a - b) <= tol * Math.max(1.0, Math.max(Math.abs(a), Math.abs(b)));
    }

    /**
     * Place-major level vector to the (nnodes x nclasses) marking matrix that a
     * Transition.setFiringRateDependence handle is written against.
     */
    private static Matrix arcMatrix(double[] m, Spn_mdd.SpnInfo info) {
        int R = info.nclasses;
        Matrix M = new Matrix(info.nnodes, R);
        for (int pp = 0; pp < info.places.length; pp++) {
            for (int k = 0; k < R; k++) {
                M.set(info.places[pp], k, m[pp * R + k]);
            }
        }
        return M;
    }

    /**
     * Index of one complex, appending it in first-seen order so that the complex
     * indices agree with the MATLAB, python and C++ twins.
     */
    private static int complexIndex(double[] v, List<double[]> clist, Map<String, Integer> seen) {
        StringBuilder sb = new StringBuilder();
        for (int l = 0; l < v.length; l++) {
            sb.append(v[l]).append(',');
        }
        String key = sb.toString();
        Integer at = seen.get(key);
        if (at != null) {
            return at.intValue();
        }
        double[] copy = new double[v.length];
        System.arraycopy(v, 0, copy, 0, v.length);
        clist.add(copy);
        seen.put(key, Integer.valueOf(clist.size() - 1));
        return clist.size() - 1;
    }

    /** Connected components of the UNDIRECTED complex graph. */
    private static int linkage(int c, int[] src, int[] dst, int[] lclass) {
        for (int v = 0; v < c; v++) {
            lclass[v] = -1;
        }
        int nlink = 0;
        for (int v = 0; v < c; v++) {
            if (lclass[v] >= 0) {
                continue;
            }
            List<Integer> stack = new ArrayList<Integer>();
            stack.add(Integer.valueOf(v));
            lclass[v] = nlink;
            while (!stack.isEmpty()) {
                int u = stack.remove(stack.size() - 1).intValue();
                for (int e = 0; e < src.length; e++) {
                    int w = -1;
                    if (src[e] == u) {
                        w = dst[e];
                    } else if (dst[e] == u) {
                        w = src[e];
                    }
                    if (w >= 0 && lclass[w] < 0) {
                        lclass[w] = nlink;
                        stack.add(Integer.valueOf(w));
                    }
                }
            }
            nlink++;
        }
        return nlink;
    }

    /** Every linkage class strongly connected in the DIRECTED complex graph. */
    private static boolean weaklyReversible(int c, int[] src, int[] dst, int[] lclass, int nlink) {
        for (int b = 0; b < nlink; b++) {
            int v0 = -1;
            for (int v = 0; v < c; v++) {
                if (lclass[v] == b) {
                    v0 = v;
                    break;
                }
            }
            boolean[] fwd = reach(v0, src, dst, c);
            boolean[] bwd = reach(v0, dst, src, c);
            for (int v = 0; v < c; v++) {
                if (lclass[v] == b && (!fwd[v] || !bwd[v])) {
                    return false;
                }
            }
        }
        return true;
    }

    private static boolean[] reach(int v0, int[] from, int[] to, int c) {
        boolean[] seen = new boolean[c];
        seen[v0] = true;
        List<Integer> stack = new ArrayList<Integer>();
        stack.add(Integer.valueOf(v0));
        while (!stack.isEmpty()) {
            int u = stack.remove(stack.size() - 1).intValue();
            for (int e = 0; e < from.length; e++) {
                if (from[e] == u && !seen[to[e]]) {
                    seen[to[e]] = true;
                    stack.add(Integer.valueOf(to[e]));
                }
            }
        }
        return seen;
    }

    private static String wrText(boolean wr) {
        return wr ? "weakly reversible" : "not weakly reversible";
    }
}
