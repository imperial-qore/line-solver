/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.tr;

import java.util.ArrayList;
import java.util.List;

import jline.api.pfqn.nc.Pfqn_bk;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.io.Ret;
import jline.lang.JobClass;
import jline.lang.ModelAdapter;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.ServiceStation;
import jline.lang.nodes.Station;
import jline.lang.processes.Distribution;
import jline.lang.processes.DistributionScaling;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * LOAD CONCEALMENT as a model transformation, and the first ITERATED strategy
 * of {@link TransformSolve}.
 *
 * <p>Birman and Kogan (Stochastic Models 8(3):543-563, 1992), Algorithm 2. The
 * saddle point analysis of their Corollary 2 shows that chain {@code l} may be
 * solved on its own provided every station is slowed by the residual capacity
 * the other chains leave it, {@code A_i = 1 - sum_{k!=l} L(i,k) X_k}, so that
 * chain {@code l} sees the concealed demand {@code L(i,l)/A_i}.
 *
 * <p>WHAT THIS ADDS OVER THE KERNEL. {@code Pfqn_bk.pfqn_bklc} solves each
 * single-chain subproblem on a DEMAND VECTOR with the inner solve hard-wired to
 * MVA or the uniform expansion. Here the subproblem is a real single-class
 * Network, so the inner solve is the CALLER'S OWN solver, which is what makes
 * the concealment approximation measurable rather than merely asserted.
 *
 * <p>IT IS NOT A STRICTLY BETTER LC. The kernel sees only {@code L}; the chain
 * aggregation refits the chain service law to two moments. On a product-form
 * model they coincide; off it they are different approximations.
 *
 * <p>THE TOLERANCE IS FIXED AT 1e-10 and is deliberately not
 * {@code options.iter_tol}: a looser one stops the sweep at a different
 * iteration in each codebase.
 *
 * <p>Mirrors MATLAB {@code solver_tr_lc_analyzer.m} and python
 * {@code transform_driver._lc_strategy}.
 */
public class LcStrategy implements TransformSolve.Strategy {

    /** The fixed convergence tolerance; a parity requirement, not a knob. */
    public static final double TOL = 1e-10;

    private static final double FINE_TOL = 1e-12;

    /** Carries the concealment state across the sweep. */
    public static final class LcContext extends TransformSolve.Expanded {
        final NetworkStruct snOrig;
        final Matrix alpha;
        final ModelAdapter.DeaggInfo deagg;
        /** Queueing demands per chain, delay rows zeroed. */
        final Matrix L;
        final double[] Z;
        final double[] N;
        final boolean[] isDelay;
        /** The unconcealed service law of every station of every subproblem. */
        final Distribution[][] baseService;
        final int M;
        final int R;
        /** True when every chain holds one class, so the lift is a re-indexing. */
        final boolean identity;
        final int korig;
        double[] X;
        double[] Xold;

        LcContext(List<Network> submodels, NetworkStruct snOrig, Matrix alpha,
                  ModelAdapter.DeaggInfo deagg, Matrix L, double[] Z, double[] N, boolean[] isDelay,
                  Distribution[][] baseService, int M, int R, boolean identity, int korig,
                  double[] X) {
            super(submodels, true);
            this.snOrig = snOrig;
            this.alpha = alpha;
            this.deagg = deagg;
            this.L = L;
            this.Z = Z;
            this.N = N;
            this.isDelay = isDelay;
            this.baseService = baseService;
            this.M = M;
            this.R = R;
            this.identity = identity;
            this.korig = korig;
            this.X = X;
            this.Xold = new double[R];
            for (int r = 0; r < R; r++) {
                this.Xold[r] = -1.0;
            }
        }
    }

    @Override
    public TransformSolve.Expanded expand(Network model, NetworkStruct sn, SolverOptions options) {
        ModelAdapter.AggregateChainResult agg = ModelAdapter.aggregateChains(model, "");
        Network chainModel = agg.getChainModel();
        NetworkStruct snChain = chainModel.getStruct(true);
        int R = snChain.nclasses;
        int M = snChain.nstations;

        Ret.snGetDemands dem = SnGetDemandsChain.snGetDemandsChain(snChain);
        Matrix Lc = dem.Dchain != null ? dem.Dchain : dem.STchain;
        boolean[] isDelay = new boolean[M];
        List<Station> stations = chainModel.getStations();
        for (int i = 0; i < M; i++) {
            isDelay[i] = snChain.sched.get(stations.get(i)) == SchedStrategy.INF;
        }
        // The concealment slows QUEUEING stations only: a delay holds no queue,
        // so zeroing its rows makes A come out as exactly 1 there.
        Matrix L = new Matrix(M, R);
        double[] Z = new double[R];
        double[] N = new double[R];
        for (int r = 0; r < R; r++) {
            N[r] = dem.Nchain.get(r);
            for (int i = 0; i < M; i++) {
                if (isDelay[i]) {
                    Z[r] += Lc.get(i, r);
                } else {
                    L.set(i, r, Lc.get(i, r));
                }
            }
        }

        // One single-class model per chain, built ONCE and re-concealed in place.
        List<JobClass> origClasses = new ArrayList<JobClass>(chainModel.getClasses());
        List<Network> submodels = new ArrayList<Network>();
        Distribution[][] baseService = new Distribution[R][M];
        for (int l = 0; l < R; l++) {
            Network m = chainModel;
            for (int k = 0; k < R; k++) {
                if (k != l) {
                    m = ModelAdapter.removeClass(m, origClasses.get(k));
                }
            }
            submodels.add(m);
            JobClass only = m.getClasses().get(0);
            List<Station> st = m.getStations();
            for (int i = 0; i < M; i++) {
                baseService[l][i] = ((ServiceStation) st.get(i)).getServiceProcess(only);
            }
        }

        LcContext ctx = new LcContext(submodels, sn, agg.getAlpha(), agg.getDeaggInfo(),
                L, Z, N, isDelay, baseService, M, R,
                sn.nchains >= sn.nclasses, sn.nclasses, seed(L, N, Z));
        concealAll(ctx);
        return ctx;
    }

    /**
     * Step 1 of Algorithm 2, reproducing the kernel's own seed exactly: the
     * saddle point utilizations of Corollary 1, the {@code N/(Z+sum L)} fallback
     * and the {@code 1/max L} capacity clamp.
     */
    static double[] seed(Matrix L, double[] N, double[] Z) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double[] X = new double[R];
        try {
            Matrix Nm = new Matrix(1, R);
            Matrix Zm = new Matrix(1, R);
            for (int r = 0; r < R; r++) {
                Nm.set(0, r, N[r]);
                Zm.set(0, r, Z[r]);
            }
            Ret.pfqnNc bk = Pfqn_bk.pfqn_bk(L, Nm, Zm);
            if (bk != null && bk.X != null && bk.X.length() == R) {
                for (int r = 0; r < R; r++) {
                    X[r] = bk.X.get(r);
                }
            }
        } catch (Exception e) {
            // the seed is a starting point, not an answer: an unusable saddle
            // point falls through to the closed-form bound below
            X = new double[R];
        }
        for (int r = 0; r < R; r++) {
            if (!Double.isFinite(X[r]) || X[r] < 0) {
                X[r] = 0;
            }
            if (X[r] == 0 && N[r] > 0) {
                double sum = Z[r];
                for (int i = 0; i < M; i++) {
                    sum += L.get(i, r);
                }
                X[r] = sum > 0 ? N[r] / sum : 0;
            }
            double cap = 0;
            for (int i = 0; i < M; i++) {
                cap = Math.max(cap, L.get(i, r));
            }
            if (cap > 0) {
                X[r] = Math.min(X[r], 1.0 / cap);
            }
        }
        return X;
    }

    /** Re-conceal every subproblem against the current throughput vector. */
    private static void concealAll(LcContext ctx) {
        for (int l = 0; l < ctx.R; l++) {
            Network m = ctx.submodels.get(l);
            JobClass only = m.getClasses().get(0);
            List<Station> st = m.getStations();
            for (int i = 0; i < ctx.M; i++) {
                if (ctx.isDelay[i]) {
                    continue;
                }
                double busy = 0;
                for (int k = 0; k < ctx.R; k++) {
                    if (k != l) {
                        busy += ctx.L.get(i, k) * ctx.X[k];
                    }
                }
                double a = Math.max(1.0 - busy, FINE_TOL);
                // scaleRate multiplies the RATE, so a factor of A_i divides the
                // mean by A_i: exactly the concealed demand L(i,l)/A_i.
                ((ServiceStation) st.get(i)).setService(only,
                        DistributionScaling.scaleRate(ctx.baseService[l][i], a));
            }
            m.refreshStruct(true);
        }
    }

    private static double chainTput(TransformSolve.Inner r) {
        if (r == null || r.X == null || r.X.length() == 0) {
            return 0;
        }
        double x = r.X.get(0);
        return (Double.isFinite(x) && x >= 0) ? x : 0;
    }

    @Override
    public TransformSolve.Expanded couple(TransformSolve.Expanded ctx, List<TransformSolve.Inner> res, int e) {
        // GAUSS-SEIDEL: chain e's throughput is published the moment it is
        // known, so chain e+1 of this same sweep already sees it.
        LcContext c = (LcContext) ctx;
        c.X[e] = chainTput(res.get(e));
        concealAll(c);
        return c;
    }

    @Override
    public boolean converged(TransformSolve.Expanded ctx, List<TransformSolve.Inner> res, int it) {
        LcContext c = (LcContext) ctx;
        double diff = 0;
        double scale = 1;
        for (int r = 0; r < c.R; r++) {
            diff = Math.max(diff, Math.abs(c.X[r] - c.Xold[r]));
            scale = Math.max(scale, Math.abs(c.X[r]));
        }
        System.arraycopy(c.X, 0, c.Xold, 0, c.R);
        return diff <= TOL * scale;
    }

    @Override
    public TransformSolve.Lifted lift(TransformSolve.Expanded ctx, List<TransformSolve.Inner> res) {
        LcContext c = (LcContext) ctx;
        int M = c.M;
        int R = c.R;
        Matrix Q = new Matrix(M, R), U = new Matrix(M, R), Rr = new Matrix(M, R), T = new Matrix(M, R);
        Matrix X = new Matrix(1, R);
        for (int l = 0; l < R; l++) {
            TransformSolve.Inner r = res.get(l);
            for (int i = 0; i < M; i++) {
                Q.set(i, l, at(r.Q, i));
                U.set(i, l, at(r.U, i));
                Rr.set(i, l, at(r.R, i));
                T.set(i, l, at(r.T, i));
            }
            X.set(0, l, chainTput(r));
        }

        if (c.identity) {
            // ONE CLASS PER CHAIN: the chain answer already IS the class answer,
            // and the aggregation returned no deaggregation tables, so the lift
            // is a re-indexing of chains onto their classes.
            int K = c.korig;
            Matrix Qo = new Matrix(M, K), Uo = new Matrix(M, K);
            Matrix Ro = new Matrix(M, K), To = new Matrix(M, K), Xo = new Matrix(1, K);
            for (int l = 0; l < R; l++) {
                int k = (int) c.snOrig.inchain.get(l).get(0);
                for (int i = 0; i < M; i++) {
                    Qo.set(i, k, Q.get(i, l));
                    Uo.set(i, k, U.get(i, l));
                    Ro.set(i, k, Rr.get(i, l));
                    To.set(i, k, T.get(i, l));
                }
                Xo.set(0, k, X.get(0, l));
            }
            Matrix Co = new Matrix(1, K);
            for (int k = 0; k < K; k++) {
                double s = 0;
                for (int i = 0; i < M; i++) {
                    s += Ro.get(i, k);
                }
                Co.set(0, k, s);
            }
            return new TransformSolve.Lifted(Qo, Uo, Ro, To, Co, Xo);
        }

        Ret.snDeaggregateChainResults d = SnDeaggregateChainResults.snDeaggregateChainResults(
                c.snOrig, c.deagg.Lchain, new Matrix(0, 0), c.deagg.STchain, c.deagg.Vchain,
                c.alpha, Q, U, Rr, T, new Matrix(0, 0), X);
        return new TransformSolve.Lifted(d.Q, d.U, d.R, d.T, d.C, d.X);
    }

    private static double at(Matrix m, int i) {
        if (m == null || i >= m.getNumRows() || m.getNumCols() == 0) {
            return 0;
        }
        return m.get(i, 0);
    }
}
