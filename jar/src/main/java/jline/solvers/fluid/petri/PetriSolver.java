/**
 * @file Fluid analysis of a stochastic Petri net: one simultaneous algebraic solve per active set.
 *
 * Port of {@code matlab/src/solvers/FLD/solver_fluid_petri.m}.
 *
 * @since LINE 3.0
 */
package jline.solvers.fluid.petri;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.LSODAExt;
import jline.solvers.fluid.moments.FluidLyapunov;
import jline.util.matrix.Matrix;

/**
 * Fluid analysis of a stochastic Petri net.
 *
 * <p>A net's fluid limit is NOT the queueing drift with places in place of stations. Its immediate
 * transitions have no rate at all: their limit is a FLOW, an algebraic unknown pinned by the
 * constraint that the input place they bind holds no mass. The marking, the diffusion covariance,
 * those flows, the multi-server phase latches and the capacity gates therefore solve
 * SIMULTANEOUSLY, as one algebraic system per active set, rather than by integrating an ODE to its
 * fixed point.
 *
 * <p>THE ACTIVE SET IS WHAT ITERATES. A pass solves the system for a fixed choice of which
 * immediate modes fire, which coordinate each pins, and which capacities bind; the answer says
 * whether that choice was right (a negative flow, a negative marking, a violated or a released
 * capacity), and the next pass makes ONE move. The moves are ordered so that a failure of each
 * invalidates the next.
 *
 * <p>Reference: {@code matlab/src/solvers/FLD/solver_fluid_petri.m}, and see
 * {@code _kb/06-solver-catalog.md} for the closures and their measured accuracy.
 */
public final class PetriSolver {

    private static final double FINE_TOL = GlobalConstants.FineTol;
    private static final double COARSE_TOL = GlobalConstants.CoarseTol;

    private final NetworkStruct sn;
    private final SolverOptions options;

    public PetriSolver(NetworkStruct sn, SolverOptions options) {
        this.sn = sn;
        this.options = options;
    }

    /** What the Petri route computes that the station table has no column for. */
    public static final class PetriReport {
        public Matrix marking;
        public Matrix markingVar;
        public List<String> modeLabel = new ArrayList<String>();
        public double[] modeFlow;
        public double[] immediateFlow;
        public List<String> invariantLabel = new ArrayList<String>();
        public double[] invariantValue;
        public double[] invariantError;
        public List<String> capacityLabel = new ArrayList<String>();
        public int[] capacityActive;
        public double[] capacityFraction;
        public int[] pinned;
        public Matrix Sigma;
    }

    /** Everything {@link #solve} returns. */
    public static final class Result {
        public Matrix QN;
        public Matrix UN;
        public Matrix RN;
        public Matrix TN;
        /** The seed trajectory: (npoints) times and (npoints x nstate) states. */
        public double[] t;
        public double[][] xvecT;
        /** Per-station-class trajectories, [station][class][point]. */
        public double[][][] QNt;
        public double[][][] UNt;
        public double[][][] TNt;
        public double[] x;
        public int iters;
        public double resnorm;
        public boolean converged;
        public double runtime;
        public Matrix Sigma;
        public Matrix QVar;
        public Matrix QStd;
        public PetriReport petri;
        public List<String> warnings = new ArrayList<String>();
        public PetriTerms terms;
    }

    /**
     * Solve the net.
     *
     * @return the station table, the marking, the covariance and the report
     */
    public Result solve() {
        final long t0 = System.nanoTime();
        final int M = sn.nstations;
        final int K = sn.nclasses;

        PetriApplicable.Verdict v = PetriApplicable.applicable(sn, options);
        if (!v.ok) {
            throw new RuntimeException(
                    "The fluid Petri route cannot solve this model: " + v.reason + ".");
        }

        final PetriTerms terms = PetriTerms.build(sn, options);
        final int n = terms.nstate;
        final int npair = terms.npair;

        int maxstate = 100;
        if (options != null && options.config != null && options.config.dae_maxstate > 0) {
            maxstate = options.config.dae_maxstate;
        }
        if (n > maxstate) {
            throw new RuntimeException("The fluid Petri route solves a " + n
                    + "-unknown algebraic system with a finite-difference Jacobian, above the "
                    + "limit of " + maxstate + " set by options.config.dae_maxstate. Raise that "
                    + "limit, or use SolverSSA for a net of this size.");
        }

        final PetriConservation cons = PetriConservation.conservation(terms);
        if (cons.leak > 1e-7) {
            throw new RuntimeException("The conserved directions and the jump matrix disagree: the "
                    + "largest leak per unit rate is " + cons.leak + ", where it must be zero.");
        }
        final PetriConstraints con = PetriConstraints.constraints(sn, terms);
        final int ncon = con.b.length;

        // ---- seed ------------------------------------------------------------
        Seed sd = seed(terms);
        double[] x = sd.xt[sd.xt.length - 1].clone();

        PetriImmediate imm = PetriImmediate.immediate(terms, x, null);
        final int nimm = imm.n;

        // THE VARIANCE IS SEEDED POSITIVE: sigma2 = 0 is where min() has no
        // derivative, and a saturated net's first-order fixed point sits there.
        double[] s2 = new double[npair];
        boolean[] ondiag = new boolean[npair];
        for (int p = 0; p < npair; p++) {
            ondiag[p] = (terms.covPairs[p][0] == terms.covPairs[p][1]);
            if (ondiag[p]) {
                s2[p] = Math.max(FINE_TOL, x[terms.covPairs[p][0]]);
            }
        }

        // THE NEWTON TOLERANCE IS THE FINE ONE BY DEFAULT: the conservation rows
        // are LINEAR in the unknowns, so the residual norm IS the token-count
        // error, and stopping at 1e-4 would leave a closed net holding 4.0001
        // tokens.
        double tol = FINE_TOL;
        if (options != null && Double.isFinite(options.tol) && options.tol < tol) {
            tol = options.tol;
        }
        int newtonMax = 50;
        if (options != null && options.iter_max > 0) {
            newtonMax = Math.max(newtonMax, options.iter_max);
        }

        // ---- steady state ----------------------------------------------------
        List<Integer> active = new ArrayList<Integer>();
        double[] phi = new double[nimm];
        final int nlatch = terms.latchMode.size();
        double[] mu = new double[nlatch];
        double[] zeta = new double[0];
        int iters = 0;
        final int asetMax = Math.max(6, 2 * (ncon + nimm) + 2);
        boolean converged = false;
        double resnorm = Double.POSITIVE_INFINITY;
        Best best = null;

        for (int sweep = 0; sweep < asetMax; sweep++) {
            final Ctx ctx = context(terms, cons, con, imm, active);
            double[] u0 = concat(x, s2, phi, mu, ones(active.size()));
            double[] lb = new double[u0.length];
            Arrays.fill(lb, Double.NEGATIVE_INFINITY);
            for (int k = 0; k < nimm; k++) {
                lb[n + npair + k] = 0.0;
            }
            for (int k = 0; k < active.size(); k++) {
                lb[n + npair + nimm + nlatch + k] = 0.0;
            }
            for (int p = 0; p < npair; p++) {
                if (ondiag[p]) {
                    lb[n + p] = 0.0;
                }
            }

            NewtonResult nr = newton(ctx, u0, tol, newtonMax, lb);
            iters += nr.iterations;
            resnorm = nr.resnorm;
            converged = nr.converged;
            Unpacked up = unpack(nr.u, terms, imm, active.size());
            x = up.x;
            s2 = up.s2;
            phi = up.phi;
            mu = up.mu;
            zeta = up.zeta;

            if (best == null || resnorm < best.resnorm) {
                best = new Best(x.clone(), s2.clone(), phi.clone(), mu.clone(), zeta.clone(),
                        new ArrayList<Integer>(active), imm, resnorm, converged);
            }

            // The active-set moves, in the order that a failure of each
            // invalidates the next: a negative flow means the mode does not fire
            // at all, a negative marking means the wrong coordinate was pinned,
            // and only then is it worth asking which capacity rows bind.
            boolean moved = false;
            if (nimm > 0) {
                double scale = 1.0;
                for (double p : phi) {
                    scale = Math.max(scale, Math.abs(p));
                }
                double thr = -Math.max(tol, 1e-10) * scale;
                for (int k = 0; k < nimm && !moved; k++) {
                    if (phi[k] < thr) {
                        imm.active[k] = false;
                        imm = PetriImmediate.immediate(terms, x, imm);
                        moved = true;
                    }
                }
            }
            if (!moved && nimm > 0) {
                double thr = -Math.max(tol, 1e-10);
                for (int s = 0; s < terms.nm && !moved; s++) {
                    if (x[s] >= thr) {
                        continue;
                    }
                    for (int k = 0; k < nimm; k++) {
                        PetriMode md = terms.modes.get(terms.immIdx.get(k));
                        if (imm.active[k] && md.arcSlot.contains(Integer.valueOf(s))
                                && imm.bind[k] != s) {
                            imm.bind[k] = s;
                            imm = PetriImmediate.immediate(terms, x, imm);
                            moved = true;
                            break;
                        }
                    }
                }
            }
            if (!moved && ncon > 0) {
                List<Integer> over = new ArrayList<Integer>();
                for (int c = 0; c < ncon; c++) {
                    double val = 0.0;
                    for (int s = 0; s < n; s++) {
                        val += con.A.get(c, s) * x[s];
                    }
                    if (val > con.b[c] + Math.max(1e-9, tol) && !active.contains(c)) {
                        over.add(c);
                    }
                }
                if (!over.isEmpty()) {
                    active.addAll(over);
                    moved = true;
                } else {
                    List<Integer> keep = new ArrayList<Integer>();
                    for (int i = 0; i < active.size(); i++) {
                        if (!(i < zeta.length && zeta[i] > 1 + Math.max(1e-9, tol))) {
                            keep.add(active.get(i));
                        }
                    }
                    if (keep.size() != active.size()) {
                        active = keep;
                        moved = true;
                    }
                }
            }
            if (!moved) {
                break;
            }
        }

        Result out = new Result();
        if (best != null && !converged && best.converged) {
            x = best.x;
            s2 = best.s2;
            phi = best.phi;
            mu = best.mu;
            zeta = best.zeta;
            active = best.active;
            imm = best.imm;
            converged = true;
            resnorm = best.resnorm;
        }

        if (!converged) {
            out.warnings.add(String.format("The simultaneous closure solve stopped at residual "
                    + "%.3e after %d Newton steps without reaching %.3e. The reported point is "
                    + "the last iterate.", resnorm, iters, tol));
        }
        // A NEGATIVE MARKING IS NOT ROUNDING: the fixed point wanted mass a place
        // cannot supply and no immediate mode could be rebound to pin it, so the
        // answer is outside the model's own state space and is reported as such.
        for (int s = 0; s < terms.nm; s++) {
            if (x[s] < -Math.max(tol, 1e-10)) {
                out.warnings.add(String.format("The fixed point holds %g tokens at %s, which is "
                        + "negative: no immediate transition could be rebound to pin that place "
                        + "at zero. Use SolverCTMC, SolverSSA or SolverLDES for this net.",
                        x[s], terms.namesNode.get(terms.coordNode[s])));
                break;
            }
        }

        final Ctx ctx = context(terms, cons, con, imm, active);
        Residual fin = residual(concat(x, s2, phi, mu, zeta), ctx, false);

        // ---- the trajectory --------------------------------------------------
        // The immediate flow ALONG the reported path is not the steady-state
        // one: on the seed path it is the large-finite-rate approximation the
        // seed itself integrated, so the throughput table and the trajectory it
        // is read off describe the same run.
        out.t = sd.t;
        out.xvecT = sd.xt;
        double[][] phit = seedFlows(terms, sd.xt, sd.lam);

        out.QN = new Matrix(M, K);
        out.UN = new Matrix(M, K);
        out.RN = new Matrix(M, K);
        out.TN = new Matrix(M, K);
        metrics(terms, x, fin.r, M, K, out.QN, out.UN, out.RN, out.TN);
        metricsT(terms, out.xvecT, s2, phit, M, K, out);

        out.x = x;
        out.iters = iters;
        out.resnorm = resnorm;
        out.converged = converged;
        out.Sigma = fin.Sigma;
        out.QVar = qvar(terms, fin.Sigma, M, K);
        out.QStd = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                out.QStd.set(i, k, Math.sqrt(Math.max(0.0, out.QVar.get(i, k))));
            }
        }
        out.petri = report(terms, cons, con, imm, active, x, fin.r, phi, zeta, fin.Sigma);
        out.terms = terms;
        out.runtime = (System.nanoTime() - t0) / 1e9;
        return out;
    }

    // ===================== the algebraic system ==============================

    /** Everything the residual needs that does not change within one active set. */
    private static final class Ctx {
        PetriTerms terms;
        PetriImmediate imm;
        int[] active;
        PetriConstraints con;
        Matrix C;
        double[] N;
        Matrix Dp;
        Matrix Dn;
    }

    /**
     * THE CONSERVATION ROWS A BINDING CAP BREAKS ARE DROPPED: a capped place loses the tokens that
     * do not fit, so a conserved quantity supported on it is not conserved while the cap binds, and
     * keeping its row would state an equation the drift contradicts -- a singular Newton system
     * rather than an inaccuracy.
     */
    private static Ctx context(PetriTerms terms, PetriConservation cons, PetriConstraints con,
                               PetriImmediate imm, List<Integer> active) {
        Ctx ctx = new Ctx();
        ctx.terms = terms;
        ctx.imm = imm;
        ctx.active = new int[active.size()];
        for (int i = 0; i < active.size(); i++) {
            ctx.active[i] = active.get(i);
        }
        ctx.con = con;
        Matrix C = cons.C;
        double[] N = cons.N;
        if (ctx.active.length > 0 && C.getNumRows() > 0) {
            boolean[] hit = new boolean[terms.nstate];
            for (int c : ctx.active) {
                for (int s = 0; s < terms.nstate; s++) {
                    if (con.cover[c][s]) {
                        hit[s] = true;
                    }
                }
            }
            List<Integer> keep = new ArrayList<Integer>();
            for (int r = 0; r < C.getNumRows(); r++) {
                boolean drop = false;
                for (int s = 0; s < terms.nstate && !drop; s++) {
                    if (hit[s] && C.get(r, s) != 0.0) {
                        drop = true;
                    }
                }
                if (!drop) {
                    keep.add(r);
                }
            }
            Matrix Ck = new Matrix(keep.size(), terms.nstate);
            double[] Nk = new double[keep.size()];
            for (int i = 0; i < keep.size(); i++) {
                for (int s = 0; s < terms.nstate; s++) {
                    Ck.set(i, s, C.get(keep.get(i), s));
                }
                Nk[i] = N[keep.get(i)];
            }
            C = Ck;
            N = Nk;
        }
        ctx.C = C;
        ctx.N = N;
        ctx.Dp = new Matrix(terms.D.getNumRows(), terms.D.getNumCols());
        ctx.Dn = new Matrix(terms.D.getNumRows(), terms.D.getNumCols());
        for (int i = 0; i < terms.D.getNumRows(); i++) {
            for (int j = 0; j < terms.D.getNumCols(); j++) {
                double d = terms.D.get(i, j);
                ctx.Dp.set(i, j, Math.max(d, 0.0));
                ctx.Dn.set(i, j, Math.min(d, 0.0));
            }
        }
        return ctx;
    }

    /** The unknown vector, split into its blocks. */
    private static final class Unpacked {
        double[] x;
        double[] s2;
        double[] phi;
        double[] mu;
        double[] zeta;
    }

    private static Unpacked unpack(double[] u, PetriTerms terms, PetriImmediate imm, int na) {
        final int n = terms.nstate;
        final int npair = terms.npair;
        final int ni = imm.n;
        final int nl = terms.latchMode.size();
        Unpacked up = new Unpacked();
        up.x = Arrays.copyOfRange(u, 0, n);
        up.s2 = Arrays.copyOfRange(u, n, n + npair);
        up.phi = Arrays.copyOfRange(u, n + npair, n + npair + ni);
        for (int i = 0; i < up.phi.length; i++) {
            up.phi[i] = Math.max(0.0, up.phi[i]);
        }
        up.mu = Arrays.copyOfRange(u, n + npair + ni, n + npair + ni + nl);
        up.zeta = Arrays.copyOfRange(u, n + npair + ni + nl, n + npair + ni + nl + na);
        for (int i = 0; i < up.zeta.length; i++) {
            up.zeta[i] = Math.max(0.0, up.zeta[i]);
        }
        return up;
    }

    /**
     * The reduction of the fluctuation onto the manifold the fast and clamped directions leave
     * free.
     *
     * <p>AN IMMEDIATE PIN REDUCES OBLIQUELY, ALONG THE FAST REACTION ITSELF, and this is the one
     * place where the orthogonal projector the queueing twin uses is WRONG rather than merely
     * different: a slow event depositing into a pinned place is answered instantly by the immediate
     * transition, so its effective jump is its own plus the immediate flow it triggers -- the token
     * is forwarded, not lost. An orthogonal projection deletes the deposit and destroys mass in the
     * diffusion.
     *
     * <pre>    P = I - Cf * G * E_B,   G = G0 * (E_B Cf G0)^-1</pre>
     *
     * <p>A CAPACITY CAP REDUCES ORTHOGONALLY -- mass that does not fit is genuinely lost, so there
     * is nothing to forward it to. A SERVER LATCH REDUCES ORTHOGONALLY TOO, on the LINEARISED row
     * {@code [-dtheta_j/dm, 1 over the phase block]}, which is why this projector depends on the
     * iterate.
     *
     * @return the reduction, or null when it is the identity
     */
    private static Matrix clampTangent(PetriTerms terms, PetriImmediate imm, PetriConstraints con,
                                       int[] active, PetriSystem.Theta th) {
        final int[] idx = terms.covIdx;
        final int nc = idx.length;
        Matrix T = Matrix.eye(nc);

        List<Integer> actk = new ArrayList<Integer>();
        for (int k = 0; k < imm.n; k++) {
            if (imm.active[k] && imm.bind[k] >= 0) {
                actk.add(k);
            }
        }
        final int[] B = imm.pins;
        if (!actk.isEmpty() && B.length > 0) {
            Matrix Cf = new Matrix(nc, actk.size());
            for (int a = 0; a < actk.size(); a++) {
                double[] cv = terms.modes.get(terms.immIdx.get(actk.get(a))).cvec;
                for (int i = 0; i < nc; i++) {
                    Cf.set(i, a, cv[idx[i]]);
                }
            }
            Matrix G0 = new Matrix(actk.size(), B.length);
            for (int jb = 0; jb < B.length; jb++) {
                List<Integer> grp = new ArrayList<Integer>();
                double wsum = 0.0;
                for (int q = 0; q < actk.size(); q++) {
                    if (imm.bind[actk.get(q)] == B[jb]) {
                        grp.add(q);
                        wsum += terms.modes.get(terms.immIdx.get(actk.get(q))).weight;
                    }
                }
                if (grp.isEmpty()) {
                    continue;
                }
                boolean uniform = !(wsum > 0);
                for (int q : grp) {
                    double w = uniform ? 1.0
                            : terms.modes.get(terms.immIdx.get(actk.get(q))).weight;
                    G0.set(q, jb, w / (uniform ? grp.size() : wsum));
                }
            }
            Matrix EB = new Matrix(B.length, nc);
            for (int jb = 0; jb < B.length; jb++) {
                EB.set(jb, B[jb], 1.0);
            }
            Matrix Mb = EB.mult(Cf).mult(G0);
            Matrix corr = Cf.mult(G0.mult(pinv(Mb))).mult(EB);
            for (int i = 0; i < nc; i++) {
                for (int j = 0; j < nc; j++) {
                    T.set(i, j, T.get(i, j) - corr.get(i, j));
                }
            }
        }

        List<double[]> R = new ArrayList<double[]>();
        for (int c : active) {
            double[] row = new double[nc];
            for (int i = 0; i < nc; i++) {
                row[i] = con.A.get(c, idx[i]);
            }
            R.add(row);
        }
        if (th != null) {
            for (int j : terms.latchMode) {
                double[] row = new double[nc];
                PetriMode md = terms.modes.get(j);
                for (int z : md.zblk) {
                    int at = indexOf(idx, z);
                    if (at >= 0) {
                        row[at] = 1.0;
                    }
                }
                if (th.dslot[j] != null) {
                    for (int q = 0; q < th.dslot[j].length; q++) {
                        int at = indexOf(idx, th.dslot[j][q]);
                        if (at >= 0) {
                            row[at] -= th.dval[j][q];
                        }
                    }
                }
                R.add(row);
            }
        }
        if (!R.isEmpty()) {
            boolean any = false;
            for (double[] row : R) {
                for (double v : row) {
                    if (Math.abs(v) > 1e-14) {
                        any = true;
                        break;
                    }
                }
            }
            if (any) {
                Matrix Rm = new Matrix(R.size(), nc);
                for (int i = 0; i < R.size(); i++) {
                    for (int j = 0; j < nc; j++) {
                        Rm.set(i, j, R.get(i)[j]);
                    }
                }
                Matrix RRt = Rm.mult(Rm.transpose());
                Matrix P = Matrix.eye(nc);
                Matrix corr = Rm.transpose().mult(pinv(RRt)).mult(Rm);
                for (int i = 0; i < nc; i++) {
                    for (int j = 0; j < nc; j++) {
                        P.set(i, j, P.get(i, j) - corr.get(i, j));
                    }
                }
                T = P.mult(T);
            }
        }
        double worst = 0.0;
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                worst = Math.max(worst, Math.abs(T.get(i, j) - (i == j ? 1.0 : 0.0)));
            }
        }
        return (worst <= 1e-14) ? null : T;
    }

    /**
     * One Lyapunov solve, over the marking coordinates.
     *
     * <p>THE DIFFUSION COUNTS THE STOCHASTIC EVENTS ONLY: an immediate flow is not a Poisson stream
     * with an intensity but the limit of an infinitely fast one whose fluctuation is slaved, and
     * its pinned coordinate is projected out.
     */
    private static Matrix sigmaOf(PetriTerms terms, Matrix A, double[] r, Matrix clampT) {
        final int[] idx = terms.covIdx;
        final int nc = idx.length;
        final int ns = terms.stochCol.length;
        Matrix Dc = new Matrix(nc, ns);
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < ns; j++) {
                Dc.set(i, j, terms.D.get(idx[i], terms.stochCol[j]));
            }
        }
        Matrix Am = new Matrix(nc, nc);
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                Am.set(i, j, A.get(idx[i], idx[j]));
            }
        }
        if (clampT != null) {
            // BOTH the jump directions and the generator are reduced: reducing
            // Dc alone would fix the subspace but leave the generator's
            // orthogonal component on it, which is not the reduced dynamics when
            // the reduction is oblique.
            Dc = clampT.mult(Dc);
            Am = clampT.mult(Am);
        }
        Matrix Q = new Matrix(nc, nc);
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                double v = 0.0;
                for (int e = 0; e < ns; e++) {
                    v += Dc.get(i, e) * r[terms.stochCol[e]] * Dc.get(j, e);
                }
                Q.set(i, j, v);
            }
        }
        FluidLyapunov.Result lr = FluidLyapunov.solve(Am, Q, Dc, Math.sqrt(2.220446049250313e-16));
        Matrix Sigma = new Matrix(terms.nstate, terms.nstate);
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                Sigma.set(idx[i], idx[j], lr.sigma.get(i, j));
            }
        }
        return Sigma;
    }

    /** The stacked residual, the rates it was evaluated at, and the covariance. */
    private static final class Residual {
        double[] G;
        double[] r;
        PetriSystem.Theta th;
        Matrix Sigma;
    }

    /**
     * The coupled algebraic system, stacked.
     *
     * <p>Returns a null {@code G} when the closure cannot be evaluated at this iterate, so the line
     * search can back off; the first evaluation of a pass runs with {@code quiet=false}, where a
     * genuine failure surfaces.
     */
    private static Residual residual(double[] u, Ctx ctx, boolean quiet) {
        final PetriTerms terms = ctx.terms;
        final PetriImmediate imm = ctx.imm;
        final int n = terms.nstate;
        Unpacked up = unpack(u, terms, imm, ctx.active.length);

        Residual res = new Residual();
        res.th = PetriSystem.theta(terms, up.x, up.s2);
        res.r = PetriSystem.rates(terms, up.x, up.phi, up.mu, res.th);

        // the deposit gate of every binding capacity, as a product of fractions
        double[] gain = new double[n];
        Arrays.fill(gain, 1.0);
        for (int k = 0; k < ctx.active.length; k++) {
            for (int s = 0; s < n; s++) {
                if (ctx.con.cover[ctx.active[k]][s]) {
                    gain[s] *= up.zeta[k];
                }
            }
        }
        double[] drift = new double[n];
        for (int s = 0; s < n; s++) {
            double neg = 0.0;
            double pos = 0.0;
            for (int e = 0; e < terms.nev; e++) {
                neg += ctx.Dn.get(s, e) * res.r[e];
                pos += ctx.Dp.get(s, e) * res.r[e];
            }
            drift[s] = neg + gain[s] * pos;
        }

        try {
            Matrix A = PetriSystem.jacobian(terms, res.th);
            res.Sigma = sigmaOf(terms, A, res.r,
                    clampTangent(terms, imm, ctx.con, ctx.active, res.th));
        } catch (RuntimeException ex) {
            if (!quiet) {
                throw ex;
            }
            res.G = null;
            return res;
        }

        List<Double> G = new ArrayList<Double>();
        for (double d : drift) {
            G.add(d);
        }
        for (int r = 0; r < ctx.C.getNumRows(); r++) {
            double v = 0.0;
            for (int s = 0; s < n; s++) {
                v += ctx.C.get(r, s) * up.x[s];
            }
            G.add(v - ctx.N[r]);
        }
        for (int p = 0; p < terms.npair; p++) {
            G.add(up.s2[p] - res.Sigma.get(terms.covPairs[p][0], terms.covPairs[p][1]));
        }
        for (PetriImmediate.Row row : imm.rows) {
            if (row.kind == PetriImmediate.PIN) {
                G.add(up.x[row.a]);
            } else if (row.kind == PetriImmediate.RATIO) {
                G.add(up.phi[row.a] * row.wb - up.phi[row.b] * row.wa);
            } else {
                G.add(up.phi[row.a]);
            }
        }
        for (int j : terms.latchMode) {
            double s = 0.0;
            for (int z : terms.modes.get(j).zblk) {
                s += up.x[z];
            }
            G.add(s - res.th.theta[j]);
        }
        for (int k = 0; k < ctx.active.length; k++) {
            double v = 0.0;
            for (int s = 0; s < n; s++) {
                v += ctx.con.A.get(ctx.active[k], s) * up.x[s];
            }
            G.add(v - ctx.con.b[ctx.active[k]]);
        }
        res.G = new double[G.size()];
        for (int i = 0; i < G.size(); i++) {
            res.G[i] = G.get(i);
        }
        return res;
    }

    // ===================== the Newton solve ==================================

    private static final class NewtonResult {
        double[] u;
        int iterations;
        boolean converged;
        double resnorm;
    }

    /**
     * An iterate projected onto its feasible box, the lower bound only.
     *
     * <p>The bound is a VECTOR, one entry per unknown, never a count: a scalar "nfree" is
     * indistinguishable from a one-unknown bound vector, which is how the MATLAB twin crashed on
     * the simplest net in the tree.
     */
    private static double[] project(double[] u, double[] lb) {
        if (lb == null || lb.length == 0) {
            return u;
        }
        if (lb.length != u.length) {
            throw new RuntimeException("The bound vector has " + lb.length + " entries for "
                    + u.length + " unknowns.");
        }
        double[] w = u.clone();
        for (int i = 0; i < w.length; i++) {
            if (Double.isFinite(lb[i])) {
                w[i] = Math.max(lb[i], w[i]);
            }
        }
        return w;
    }

    /** Damped projected Newton with a finite-difference Jacobian. */
    private static NewtonResult newton(Ctx ctx, double[] u0, double tol, int maxit, double[] lb) {
        NewtonResult out = new NewtonResult();
        double[] u = project(u0.clone(), lb);
        Residual res = residual(u, ctx, false);
        if (res.G == null) {
            out.u = u;
            out.iterations = 0;
            out.converged = false;
            out.resnorm = Double.POSITIVE_INFINITY;
            return out;
        }
        double[] G = res.G;
        double resnorm = infNorm(G);
        int it = 0;
        for (it = 1; it <= maxit; it++) {
            if (resnorm <= tol) {
                out.u = u;
                out.iterations = it - 1;
                out.converged = true;
                out.resnorm = resnorm;
                return out;
            }
            Matrix J = fdjac(ctx, u, G);
            Matrix rhs = new Matrix(G.length, 1);
            for (int i = 0; i < G.length; i++) {
                rhs.set(i, 0, -G[i]);
            }
            Matrix duM;
            try {
                duM = lstsq(J, rhs);
            } catch (RuntimeException ex) {
                break;
            }
            double lam = 1.0;
            boolean improved = false;
            for (int b = 0; b < 30; b++) {
                double[] un = new double[u.length];
                for (int i = 0; i < u.length; i++) {
                    un[i] = u[i] + lam * duM.get(i, 0);
                }
                un = project(un, lb);
                Residual rn = residual(un, ctx, true);
                if (rn.G != null) {
                    double v = infNorm(rn.G);
                    if (v < resnorm) {
                        u = un;
                        G = rn.G;
                        resnorm = v;
                        improved = true;
                        break;
                    }
                }
                lam *= 0.5;
            }
            if (!improved) {
                break;
            }
        }
        out.u = u;
        out.iterations = it;
        out.converged = resnorm <= tol;
        out.resnorm = resnorm;
        return out;
    }

    /** Forward-difference Jacobian of the residual. */
    private static Matrix fdjac(Ctx ctx, double[] u, double[] G) {
        final int n = u.length;
        final int m = G.length;
        Matrix J = new Matrix(m, n);
        for (int k = 0; k < n; k++) {
            double h = 1e-7 * Math.max(1.0, Math.abs(u[k]));
            double[] up = u.clone();
            up[k] += h;
            Residual rp = residual(up, ctx, true);
            if (rp.G != null) {
                for (int i = 0; i < m; i++) {
                    J.set(i, k, (rp.G[i] - G[i]) / h);
                }
                continue;
            }
            up[k] = u[k] - h;
            Residual rm = residual(up, ctx, true);
            if (rm.G == null) {
                continue;
            }
            for (int i = 0; i < m; i++) {
                J.set(i, k, (G[i] - rm.G[i]) / h);
            }
        }
        return J;
    }

    // ===================== the seed ==========================================

    private static final class Seed {
        double[] t;
        double[][] xt;
        double lam;
    }

    /**
     * The first-order drift: the same rates at zero variance, with an immediate mode firing at LAM
     * times its enabling degree and its firing weight, and the server latch relaxed at LAM towards
     * the enabling degree instead of solved. Both are approximations of an algebraic constraint by
     * a fast reaction, and both are confined to the seed.
     */
    private static double[] seedDrift(PetriTerms terms, double[] xin, double lam) {
        double[] x = xin.clone();
        for (int i = 0; i < x.length; i++) {
            x[i] = Math.max(0.0, x[i]);
        }
        PetriSystem.Theta th = PetriSystem.theta(terms, x, new double[Math.max(terms.npair, 1)]);
        double[] phi = new double[terms.immIdx.size()];
        for (int k = 0; k < phi.length; k++) {
            int j = terms.immIdx.get(k);
            phi[k] = lam * terms.modes.get(j).weight * th.theta[j];
        }
        double[] mu = new double[terms.latchMode.size()];
        for (int q = 0; q < mu.length; q++) {
            int j = terms.latchMode.get(q);
            double s = 0.0;
            for (int z : terms.modes.get(j).zblk) {
                s += x[z];
            }
            mu[q] = lam * (th.theta[j] - s);
        }
        double[] r = PetriSystem.rates(terms, x, phi, mu, th);
        double[] d = new double[terms.nstate];
        for (int s = 0; s < terms.nstate; s++) {
            double v = 0.0;
            for (int e = 0; e < terms.nev; e++) {
                v += terms.D.get(s, e) * r[e];
            }
            d[s] = v;
        }
        return d;
    }

    /**
     * One first-order trajectory from the initial marking.
     *
     * <p>Newton needs a point in the basin, not an answer. The immediate modes get a large FINITE
     * rate here and only here, scaled to the model's own timescale rather than taken from
     * {@code GlobalConstants.Immediate}: 1e8 against a rate of order one is a stiffness the seed
     * does not need, and the answer does not depend on the seed's accuracy.
     */
    private Seed seed(final PetriTerms terms) {
        double rmax = 0.0;
        for (int j : terms.timedIdx) {
            double s = 0.0;
            for (double v : terms.modes.get(j).d1) {
                s += v;
            }
            rmax = Math.max(rmax, s);
        }
        for (int e = 0; e < terms.nev; e++) {
            if (terms.evKind[e] == 3) {
                rmax = Math.max(rmax, terms.rateBase[e]);
            }
        }
        if (!(rmax > 0)) {
            rmax = 1.0;
        }
        final double lam = Math.min(1e8, 1e4 * rmax);

        double mass = 0.0;
        for (int s = 0; s < terms.nm; s++) {
            mass += terms.x0[s];
        }
        double T = 50.0 * (mass + 1.0) / rmax;

        final int n = terms.nstate;
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return n;
            }

            @Override
            public void computeDerivatives(double tt, double[] y, double[] yDot) {
                double[] d = seedDrift(terms, y, lam);
                System.arraycopy(d, 0, yDot, 0, n);
            }
        };

        List<double[]> path = new ArrayList<double[]>();
        List<Double> times = new ArrayList<Double>();
        double[] xend = terms.x0.clone();
        for (int attempt = 0; attempt < 6; attempt++) {
            path.clear();
            times.clear();
            // The trajectory is sampled on a fixed grid rather than at the
            // integrator's own steps: the reported path only has to show the
            // approach, and a uniform grid keeps the transient table the same
            // shape whatever the stiffness made the stepper do.
            final int npoints = 200;
            double[] y = terms.x0.clone();
            path.add(y.clone());
            times.add(0.0);
            // The stepper's own tuning, as SolverFluid constructs it: Adams for a
            // non-stiff stretch, BDF for the stiff one the large finite immediate
            // rate creates, with the step count high enough that the seed is not
            // cut short by it.
            LSODAExt ode1 = new LSODAExt(0.0, 0.0, COARSE_TOL, FINE_TOL, 12, 5, 10000000);
            boolean failed = false;
            for (int p = 1; p < npoints; p++) {
                double ta = T * (p - 1) / (npoints - 1.0);
                double tb = T * p / (npoints - 1.0);
                double[] yOut = new double[n];
                try {
                    ode1.integrate(ode, ta, y, tb, yOut);
                } catch (RuntimeException ex) {
                    failed = true;
                    break;
                }
                y = yOut.clone();
                path.add(y.clone());
                times.add(tb);
            }
            xend = y;
            if (failed) {
                T *= 4.0;
                continue;
            }
            double[] d = seedDrift(terms, xend, lam);
            if (infNorm(d) <= COARSE_TOL * Math.max(1.0, rmax * (mass + 1.0))) {
                break;
            }
            T *= 4.0;
        }
        Seed sd = new Seed();
        sd.t = new double[times.size()];
        for (int i = 0; i < times.size(); i++) {
            sd.t[i] = times.get(i);
        }
        sd.xt = path.toArray(new double[0][]);
        sd.lam = lam;
        return sd;
    }

    /**
     * The immediate flows the SEED trajectory carried, point by point -- the approximation the seed
     * integrated, and therefore the one its throughput table has to be read at.
     */
    private static double[][] seedFlows(PetriTerms terms, double[][] xt, double lam) {
        final int nt = xt.length;
        final int ni = terms.immIdx.size();
        double[][] phit = new double[nt][ni];
        if (ni == 0) {
            return phit;
        }
        for (int a = 0; a < nt; a++) {
            PetriSystem.Theta th =
                    PetriSystem.theta(terms, xt[a], new double[Math.max(terms.npair, 1)]);
            for (int k = 0; k < ni; k++) {
                int j = terms.immIdx.get(k);
                phit[a][k] = lam * terms.modes.get(j).weight * th.theta[j];
            }
        }
        return phit;
    }

    // ===================== the metrics =======================================

    /**
     * The station table, in the conventions the exact SPN engines report.
     *
     * <p>A PLACE IS AN INF STATION: its queue length is its mean token count and its utilization is
     * the same number. Its throughput is the rate at which TOKENS leave it, so a consuming mode
     * contributes its firing rate times the arc multiplicity -- SolverCTMC's convention, and the
     * one Little's law needs. SolverSSA's NRM sums the UNWEIGHTED propensity, so the two disagree
     * wherever an input arc has multiplicity above one; the exact engine is the reference.
     */
    private static void metrics(PetriTerms terms, double[] x, double[] r, int M, int K,
                                Matrix QN, Matrix UN, Matrix RN, Matrix TN) {
        for (int s = 0; s < terms.nm; s++) {
            int ist = terms.coordStation[s];
            int k = terms.coordClass[s];
            if (ist < 0 || ist >= M) {
                continue;
            }
            QN.set(ist, k, QN.get(ist, k) + x[s]);
            UN.set(ist, k, QN.get(ist, k));
        }
        for (Map.Entry<Integer, List<Integer>> en : terms.consumers.entrySet()) {
            int key = en.getKey();
            int ist = key / K;
            int k = key % K;
            List<Double> w = terms.consumerW.get(key);
            double v = 0.0;
            for (int i = 0; i < en.getValue().size(); i++) {
                v += w.get(i) * r[en.getValue().get(i)];
            }
            TN.set(ist, k, TN.get(ist, k) + v);
        }
        for (Map.Entry<Integer, List<Integer>> en : terms.producers.entrySet()) {
            int key = en.getKey();
            int ist = key / K;
            int k = key % K;
            double v = 0.0;
            for (int e : en.getValue()) {
                v += r[e];
            }
            TN.set(ist, k, TN.get(ist, k) + v);
        }
        // TN is zero only to the integrator's accuracy.
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                if (TN.get(i, k) > GlobalConstants.Zero) {
                    RN.set(i, k, QN.get(i, k) / TN.get(i, k));
                }
            }
        }
    }

    /** The same reader along a trajectory. */
    private static void metricsT(PetriTerms terms, double[][] xt, double[] s2, double[][] phit,
                                 int M, int K, Result out) {
        final int nt = xt.length;
        out.QNt = new double[M][K][nt];
        out.UNt = new double[M][K][nt];
        out.TNt = new double[M][K][nt];
        double[][] Rt = new double[nt][];
        for (int a = 0; a < nt; a++) {
            // the latch flow moves no marking, so it cannot change a place throughput
            PetriSystem.Theta th = PetriSystem.theta(terms, xt[a], s2);
            Rt[a] = PetriSystem.rates(terms, xt[a], phit.length > 0 ? phit[a] : null, null, th);
        }
        for (int s = 0; s < terms.nm; s++) {
            int ist = terms.coordStation[s];
            int k = terms.coordClass[s];
            if (ist < 0 || ist >= M) {
                continue;
            }
            for (int a = 0; a < nt; a++) {
                out.QNt[ist][k][a] += xt[a][s];
                out.UNt[ist][k][a] = out.QNt[ist][k][a];
            }
        }
        for (Map.Entry<Integer, List<Integer>> en : terms.consumers.entrySet()) {
            int key = en.getKey();
            int ist = key / K;
            int k = key % K;
            List<Double> w = terms.consumerW.get(key);
            for (int a = 0; a < nt; a++) {
                for (int i = 0; i < en.getValue().size(); i++) {
                    out.TNt[ist][k][a] += w.get(i) * Rt[a][en.getValue().get(i)];
                }
            }
        }
        for (Map.Entry<Integer, List<Integer>> en : terms.producers.entrySet()) {
            int key = en.getKey();
            int ist = key / K;
            int k = key % K;
            for (int a = 0; a < nt; a++) {
                for (int e : en.getValue()) {
                    out.TNt[ist][k][a] += Rt[a][e];
                }
            }
        }
    }

    /** The variance of each station-class token count. */
    private static Matrix qvar(PetriTerms terms, Matrix Sigma, int M, int K) {
        Matrix QVar = new Matrix(M, K);
        for (int s = 0; s < terms.nm; s++) {
            int ist = terms.coordStation[s];
            if (ist < 0 || ist >= M) {
                continue;
            }
            QVar.set(ist, terms.coordClass[s], Math.max(0.0, Sigma.get(s, s)));
        }
        return QVar;
    }

    /** What the Petri route computes that the station table has no column for. */
    private static PetriReport report(PetriTerms terms, PetriConservation cons,
                                      PetriConstraints con, PetriImmediate imm,
                                      List<Integer> active, double[] x, double[] r, double[] phi,
                                      double[] zeta, Matrix Sigma) {
        PetriReport rep = new PetriReport();
        rep.marking = new Matrix(terms.I, terms.K);
        rep.markingVar = new Matrix(terms.I, terms.K);
        for (int s = 0; s < terms.nm; s++) {
            rep.marking.set(terms.coordNode[s], terms.coordClass[s], x[s]);
            rep.markingVar.set(terms.coordNode[s], terms.coordClass[s],
                    Math.max(0.0, Sigma.get(s, s)));
        }
        for (PetriMode md : terms.modes) {
            rep.modeLabel.add(md.label);
        }
        rep.modeFlow = new double[terms.modes.size()];
        for (int j = 0; j < terms.modes.size(); j++) {
            double v = 0.0;
            for (int e = 0; e < terms.nev; e++) {
                if (terms.evMode[e] == j && (terms.evKind[e] == 1 || terms.evKind[e] == 4)) {
                    v += r[e];
                }
            }
            rep.modeFlow[j] = v;
        }
        rep.immediateFlow = phi.clone();
        rep.invariantLabel = cons.label;
        rep.invariantValue = cons.N.clone();
        rep.invariantError = new double[cons.C.getNumRows()];
        for (int c = 0; c < cons.C.getNumRows(); c++) {
            double v = 0.0;
            for (int s = 0; s < terms.nstate; s++) {
                v += cons.C.get(c, s) * x[s];
            }
            rep.invariantError[c] = v - cons.N[c];
        }
        rep.capacityLabel = con.label;
        rep.capacityActive = new int[active.size()];
        for (int i = 0; i < active.size(); i++) {
            rep.capacityActive[i] = active.get(i);
        }
        rep.capacityFraction = zeta.clone();
        rep.pinned = imm.pins.clone();
        rep.Sigma = Sigma;
        return rep;
    }

    // ===================== small linear algebra ==============================

    /**
     * Moore-Penrose pseudo-inverse.
     *
     * {@link Matrix#pinv()} is used rather than a regularised normal-equations solve, and the
     * difference is not cosmetic: the oblique reduction of an immediate pin must ANNIHILATE the
     * pinned coordinate exactly, and a Tikhonov term of 1e-12 leaves that entry at 1e-12 instead
     * of zero. The reduced generator then keeps a marginal eigenvalue at -1e-13 and the Lyapunov
     * solve refuses the fixed point as non-hyperbolic, which is a failure to reduce reported as a
     * property of the model. Matrix.pinv uses MATLAB's own rank tolerance and gives the exact
     * inverse on the well-conditioned systems this raises.
     */
    private static Matrix pinv(Matrix A) {
        return A.pinv();
    }

    /** Least-squares solve of {@code J du = rhs}, minimum-norm where J is rank deficient. */
    private static Matrix lstsq(Matrix J, Matrix rhs) {
        return J.pinv().mult(rhs);
    }

    private static double infNorm(double[] v) {
        double m = 0.0;
        for (double d : v) {
            m = Math.max(m, Math.abs(d));
        }
        return m;
    }

    private static double[] ones(int n) {
        double[] v = new double[n];
        Arrays.fill(v, 1.0);
        return v;
    }

    private static int indexOf(int[] a, int v) {
        for (int i = 0; i < a.length; i++) {
            if (a[i] == v) {
                return i;
            }
        }
        return -1;
    }

    private static double[] concat(double[]... parts) {
        int n = 0;
        for (double[] p : parts) {
            n += p.length;
        }
        double[] out = new double[n];
        int at = 0;
        for (double[] p : parts) {
            System.arraycopy(p, 0, out, at, p.length);
            at += p.length;
        }
        return out;
    }

    /** The best iterate seen across the active-set sweeps. */
    private static final class Best {
        final double[] x;
        final double[] s2;
        final double[] phi;
        final double[] mu;
        final double[] zeta;
        final List<Integer> active;
        final PetriImmediate imm;
        final double resnorm;
        final boolean converged;

        Best(double[] x, double[] s2, double[] phi, double[] mu, double[] zeta,
             List<Integer> active, PetriImmediate imm, double resnorm, boolean converged) {
            this.x = x;
            this.s2 = s2;
            this.phi = phi;
            this.mu = mu;
            this.zeta = zeta;
            this.active = active;
            this.imm = imm;
            this.resnorm = resnorm;
            this.converged = converged;
        }
    }
}
