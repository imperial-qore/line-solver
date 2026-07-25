/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.handlers;

import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.QueueNodeParam;
import jline.lang.nodes.Node;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import jline.api.pfqn.nc.Pfqn_ncoi;
import jline.api.pfqn.ld.Pfqn_oi_fnc;
import jline.api.pfqn.ld.Pfqn_oi_insvc;
import jline.io.Ret;

import java.util.ArrayList;
import java.util.List;
import java.util.function.ToDoubleFunction;

/**
 * Exact normalizing-constant analysis of a closed queueing network that mixes
 * order-independent (OI) stations with ordinary BCMP product-form stations.
 * Port of matlab/src/solvers/NC/solver_nc_oi_analyzer.m.
 *
 * Supported stations: OI (SchedStrategy.OI / PAS with empty swap graph),
 * analyzed by the balanced-fairness rank rate; and any BCMP product-form
 * station (infinite server / delay, PS, LCFS-PR, SIRO, class-independent-rate
 * FCFS, single- or multi-server), analyzed by the load-dependent BCMP weight
 * table W_i(n) = (sum n)!/prod(n_r!) * prod_r D_{i,r}^{n_r} / prod_k min(k,c).
 * The full-network normalizing-constant table G(P) is assembled by convolution
 * (OI + delay via {@link Pfqn_ncoi}, then each BCMP queue folded in); the exact
 * per-class mean queue length at any station follows from the OI functional-
 * server (FNC) identity {@link Pfqn_oi_fnc}: E[n_{i,r}] = G^{+}_{i,r}/G - 1.
 * General per-class visits at OI stations are supported via a v-weighted
 * balanced-fairness balance Phi^v(n)=(1/mu(n)) sum_r v_r Phi^v(n-e_r).
 */
public final class Solver_nc_oi {
    private Solver_nc_oi() {}

    /**
     * True when the model is a closed queueing network with at least one OI
     * station and every other station a BCMP product-form station: infinite
     * server (delay), PS, LCFS-PR, SIRO or class-independent-rate FCFS. The OI
     * requirement keeps pure-BCMP networks on the standard (faster) path.
     */
    public static boolean nc_is_oi_model(NetworkStruct sn) {
        for (int r = 0; r < sn.njobs.getNumElements(); r++) {
            if (Double.isInfinite(sn.njobs.get(r))) {
                return false;
            }
        }
        boolean hasOI = false;
        for (int ist = 0; ist < sn.nstations; ist++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(ist));
            if (s == SchedStrategy.INF) {
                continue;
            } else if (s == SchedStrategy.PAS || s == SchedStrategy.OI) {
                Node node = sn.nodes.get((int) sn.stationToNode.get(ist));
                NodeParam np = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
                if (!(np instanceof QueueNodeParam)) {
                    return false;
                }
                Matrix sg = ((QueueNodeParam) np).swapGraph;
                if (sg == null || hasNonZero(sg)) {
                    return false;  // genuine pass-and-swap: not order-independent
                }
                hasOI = true;
            } else if (isBcmpSched(s)) {
                if ((s == SchedStrategy.FCFS || s == SchedStrategy.SIRO) && !fcfsRateOk(sn, ist)) {
                    return false;  // class-dependent FCFS/SIRO: not product form
                }
            } else {
                return false;      // an unsupported (non-product-form) station
            }
        }
        return hasOI;
    }

    public static SolverNC.SolverNCReturn solver_nc_oi(NetworkStruct sn, SolverOptions options) {
        long tStart = System.nanoTime();
        int M = sn.nstations;
        int K = sn.nclasses;

        // ---- reject class switching (OI rank rates are per raw class) ------
        for (int c = 0; c < sn.nchains; c++) {
            if (sn.inchain.get(c).getNumElements() > 1) {
                throw new RuntimeException("solver_nc_oi requires one class per chain (no class switching).");
            }
        }
        for (int r = 0; r < sn.njobs.getNumElements(); r++) {
            if (Double.isInfinite(sn.njobs.get(r))) {
                throw new RuntimeException("solver_nc_oi requires a closed queueing network.");
            }
        }
        int[] N = new int[K];
        for (int r = 0; r < K; r++) {
            N[r] = (int) Math.round(sn.njobs.get(r));
        }

        // ---- classify stations ---------------------------------------------
        // OI: rank-rate balanced-fairness station. INF: aggregated into delay Z.
        // Q : ordinary BCMP product-form station (PS / LCFS-PR / FCFS / SIRO).
        boolean[] isOI = new boolean[M];
        boolean[] isINF = new boolean[M];
        boolean[] isQ = new boolean[M];
        List<Integer> oiList = new ArrayList<Integer>();
        List<SerializableFunction<Matrix, Double>> oiSvc = new ArrayList<SerializableFunction<Matrix, Double>>();
        for (int ist = 0; ist < M; ist++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(ist));
            if (s == SchedStrategy.INF) {
                isINF[ist] = true;
            } else if (s == SchedStrategy.PAS || s == SchedStrategy.OI) {
                Node node = sn.nodes.get((int) sn.stationToNode.get(ist));
                NodeParam np = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
                Matrix sg = (np instanceof QueueNodeParam) ? ((QueueNodeParam) np).swapGraph : null;
                if (sg == null || hasNonZero(sg)) {
                    throw new RuntimeException("solver_nc_oi supports OI stations only (PAS with a non-empty swap graph is not order-independent).");
                }
                isOI[ist] = true;
                SerializableFunction<Matrix, Double> fun = ((QueueNodeParam) np).svcRateFun;
                if (fun == null) {
                    throw new RuntimeException("OI station " + ist + " has no service rate function; set it via setServiceRateFunction.");
                }
                oiList.add(ist);
                oiSvc.add(fun);
            } else if (isBcmpSched(s)) {
                isQ[ist] = true;
                if ((s == SchedStrategy.FCFS || s == SchedStrategy.SIRO) && !fcfsRateOk(sn, ist)) {
                    throw new RuntimeException("Station " + ist + " has class-dependent FCFS/SIRO rates and is not product form; solver_nc_oi requires class-independent rates.");
                }
            } else {
                throw new RuntimeException("solver_nc_oi supports only INF (delay), OI, PS, LCFS-PR, SIRO and class-independent FCFS stations.");
            }
        }

        // ---- per-class visits (chain == class); normalize to ref station ---
        Matrix V = new Matrix(M, K);
        V.zero();
        for (int r = 0; r < K; r++) {
            int c = -1;
            for (int cc = 0; cc < sn.nchains; cc++) {
                if (sn.chains.get(cc, r) != 0) {
                    c = cc;
                    break;
                }
            }
            Matrix vis = sn.visits.get(c);               // (nstateful x nclasses)
            for (int ist = 0; ist < M; ist++) {
                int isf = (int) sn.stationToStateful.get(ist);
                V.set(ist, r, vis.get(isf, r));
            }
            double vref = V.get((int) sn.refstat.get(r), r);
            if (vref > 0) {
                for (int ist = 0; ist < M; ist++) {
                    V.set(ist, r, V.get(ist, r) / vref);
                }
            }
        }

        // ---- per-class demand and aggregated delay demand Z_r --------------
        Matrix ST = sn.rates.copy();
        for (int i = 0; i < ST.getNumRows(); i++) {
            for (int j = 0; j < ST.getNumCols(); j++) {
                double v = 1.0 / ST.get(i, j);
                ST.set(i, j, Double.isFinite(v) ? v : 0.0);
            }
        }
        double[] Z = new double[K];
        for (int ist = 0; ist < M; ist++) {
            if (!isINF[ist]) continue;
            for (int r = 0; r < K; r++) {
                Z[r] += V.get(ist, r) * ST.get(ist, r);
            }
        }
        double[][] D = new double[M][K];             // per-class demand at BCMP queues
        for (int ist = 0; ist < M; ist++) {
            if (!isQ[ist]) continue;
            for (int r = 0; r < K; r++) {
                D[ist][r] = V.get(ist, r) * ST.get(ist, r);
            }
        }

        // OI-station class visit ratios feed the v-weighted balance; general
        // (non-unit) visits are supported. oivis[m] is the 1xR vector of OI
        // station oiList.get(m).
        double[][] oivis = new double[oiList.size()][K];
        for (int m = 0; m < oiList.size(); m++) {
            int ist = oiList.get(m);
            for (int r = 0; r < K; r++) {
                oivis[m][r] = V.get(ist, r);
            }
        }

        // ---- OI rank-rate handles on a per-class count vector --------------
        // The ordered-list rate mu(c) is permutation-invariant for an OI
        // station, hence a function of the multiset (per-class counts) of n;
        // evaluate it on a canonical microstate with n_r copies of class r.
        List<ToDoubleFunction<int[]>> rates = new ArrayList<ToDoubleFunction<int[]>>();
        for (int m = 0; m < oiList.size(); m++) {
            rates.add(makeRankRate(oiSvc.get(m)));
        }

        // ---- population lattice --------------------------------------------
        int[] shp = new int[K];
        int[] stride = new int[K];
        int total = 1;
        for (int d = 0; d < K; d++) {
            shp[d] = N[d] + 1;
            total *= shp[d];
        }
        stride[0] = 1;
        for (int d = 1; d < K; d++) {
            stride[d] = stride[d - 1] * shp[d - 1];
        }

        // ---- core normalizing-constant table (OI stations + delay) ---------
        // One call: the balanced-fairness convolution is a lattice convolution,
        // so its internal table already holds G(n) for every 0 <= n <= N on the
        // same column-major stride used here. Re-calling it per population would
        // cost a needless factor total = prod_r (N_r + 1).
        double[] Gfull = Pfqn_ncoi.pfqn_ncoi(Z, N, rates, oivis).Gtab;

        // ---- fold the BCMP queueing stations by lattice convolution --------
        List<Integer> qList = new ArrayList<Integer>();
        for (int ist = 0; ist < M; ist++) {
            if (isQ[ist]) qList.add(ist);
        }
        for (int ist : qList) {
            double[] Wq = ldTable(D[ist], sn.nservers.get(ist), shp, total);
            Gfull = convLattice(Gfull, Wq, shp, stride, total);
        }

        double G = Gfull[total - 1];
        double lG = Math.log(G);

        // ---- per-class throughput X_r = G(N - e_r)/G(N) -------------------
        Matrix X = new Matrix(1, K);
        X.zero();
        for (int r = 0; r < K; r++) {
            if (N[r] > 0) {
                int idx = total - 1 - stride[r];
                X.set(r, Gfull[idx] / G);
            }
        }

        // ---- per-station per-class mean queue length via the FNC identity --
        Matrix Q = new Matrix(M, K);
        Q.zero();
        for (int m = 0; m < oiList.size(); m++) {
            int ist = oiList.get(m);
            Matrix Phi = oi_phi(rates.get(m), N, oivis[m]);
            for (int r = 0; r < K; r++) {
                if (N[r] > 0) {
                    Ret.pfqnOifnc fres = Pfqn_oi_fnc.pfqn_oi_fnc(Phi, N, makeTarget(r));
                    Q.set(ist, r, fncMean(fres.Psi, Gfull, shp, stride, total) / G - 1);
                }
            }
        }
        for (int ist : qList) {
            Matrix Wq = toMatrix(ldTable(D[ist], sn.nservers.get(ist), shp, total));
            for (int r = 0; r < K; r++) {
                if (N[r] > 0) {
                    Ret.pfqnOifnc fres = Pfqn_oi_fnc.pfqn_oi_fnc(Wq, N, makeTarget(r));
                    Q.set(ist, r, fncMean(fres.Psi, Gfull, shp, stride, total) / G - 1);
                }
            }
        }
        for (int ist = 0; ist < M; ist++) {
            if (!isINF[ist]) continue;
            for (int r = 0; r < K; r++) {
                Q.set(ist, r, X.get(r) * V.get(ist, r) * ST.get(ist, r));
            }
        }

        // ---- throughput, utilization, response time per station -----------
        Matrix T = new Matrix(M, K);
        Matrix U = new Matrix(M, K);
        Matrix R = new Matrix(M, K);
        T.zero();
        U.zero();
        R.zero();
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
                T.set(ist, r, X.get(r) * V.get(ist, r));
            }
        }
        for (int ist = 0; ist < M; ist++) {
            if (!isINF[ist]) continue;
            for (int r = 0; r < K; r++) {
                U.set(ist, r, Q.get(ist, r));   // INF utilization convention
            }
        }
        for (int ist : qList) {                 // BCMP queue: offered-load per server
            double c = sn.nservers.get(ist);
            if (!Double.isFinite(c) || c <= 0) {
                c = 1;
            }
            for (int r = 0; r < K; r++) {
                U.set(ist, r, X.get(r) * D[ist][r] / c);
            }
        }
        for (int m = 0; m < oiList.size(); m++) {
            // In-service utilization U_r = E[sir_r]/c, with sir_r the number of
            // class-r jobs receiving a strictly positive rank rate.
            // E[sir_r | n] = g(n,r) is a function of the count vector
            // (Pfqn_oi_insvc), so its mean is read off the functional-server
            // identity E[f(n)] = G^{+}/G - 1 of Pfqn_oi_fnc. This matches the
            // exact CTMC and LDES convention; it coincides with the offered-load
            // form T/mu(e_r)/c only when a job engages a single server.
            int ist = oiList.get(m);
            double S = sn.nservers.get(ist);
            if (!Double.isFinite(S) || S <= 0) {
                S = 1;
            }
            Matrix PhiU = oi_phi(rates.get(m), N, oivis[m]);
            Pfqn_oi_insvc.Result ins = Pfqn_oi_insvc.pfqn_oi_insvc(rates.get(m), N);
            for (int r = 0; r < K; r++) {
                if (N[r] > 0) {
                    Ret.pfqnOifnc fres = Pfqn_oi_fnc.pfqn_oi_fnc(PhiU, N, makeInsvcTarget(ins.g, stride, r));
                    U.set(ist, r, (fncMean(fres.Psi, Gfull, shp, stride, total) / G - 1) / S);
                }
            }
        }
        for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
                if (T.get(ist, r) > 0) {
                    R.set(ist, r, Q.get(ist, r) / T.get(ist, r));
                }
            }
        }

        Matrix STeff = new Matrix(M, K);
        STeff.zero();
        double runtime = (System.nanoTime() - tStart) / 1.0e9;
        return new SolverNC.SolverNCReturn(Q, U, R, T, sn.nchains, X, lG, STeff, 1, runtime, "oi");
    }

    // ======================================================================
    private static boolean isBcmpSched(SchedStrategy s) {
        return s == SchedStrategy.PS || s == SchedStrategy.LCFSPR
                || s == SchedStrategy.FCFS || s == SchedStrategy.SIRO;
    }

    private static boolean fcfsRateOk(NetworkStruct sn, int ist) {
        // False when station ist has class-dependent FCFS/SIRO rates.
        double lo = Double.POSITIVE_INFINITY, hi = 0;
        for (int r = 0; r < sn.nclasses; r++) {
            double rate = sn.rates.get(ist, r);
            if (Double.isFinite(rate) && sn.njobs.get(r) > 0) {
                lo = Math.min(lo, rate);
                hi = Math.max(hi, rate);
            }
        }
        if (hi == 0) {
            return true;
        }
        return (hi - lo) <= 1e-9 * hi;
    }

    private static ToDoubleFunction<int[]> makeRankRate(final SerializableFunction<Matrix, Double> fun) {
        // OI rank rate on a per-class count vector n. The service function
        // svcRateFun(c) takes an ordered microstate list; for an order-
        // independent station it is permutation-invariant, hence a function of
        // the multiset of present jobs (class multiplicities included, not
        // merely the support). Evaluate it on a canonical 0-based microstate
        // holding n_r copies of class r, passed as a (1 x sum n) row Matrix.
        return new ToDoubleFunction<int[]>() {
            @Override
            public double applyAsDouble(int[] n) {
                int len = 0;
                for (int v : n) {
                    if (v > 0) len += v;
                }
                Matrix c = new Matrix(1, len);
                int col = 0;
                for (int r = 0; r < n.length; r++) {
                    for (int t = 0; t < n[r]; t++) {
                        c.set(0, col++, r);
                    }
                }
                return fun.apply(c);
            }
        };
    }

    private static ToDoubleFunction<int[]> makeInsvcTarget(final double[][] g, final int[] stride, final int r) {
        // Target f(n) = E[sir_r | n] for the FNC balance construction.
        return new ToDoubleFunction<int[]>() {
            @Override
            public double applyAsDouble(int[] n) {
                int idx = 0;
                for (int d = 0; d < n.length; d++) {
                    idx += n[d] * stride[d];
                }
                return g[idx][r];
            }
        };
    }

    private static ToDoubleFunction<int[]> makeTarget(final int r) {
        // Target f(n) = n_r for the FNC balance construction.
        return new ToDoubleFunction<int[]>() {
            @Override
            public double applyAsDouble(int[] n) {
                return n[r];
            }
        };
    }

    private static int[] decodeSub(int i, int[] shp) {
        int R = shp.length;
        int[] n = new int[R];
        int li = i;
        for (int d = 0; d < R; d++) {
            n[d] = li % shp[d];
            li /= shp[d];
        }
        return n;
    }

    /**
     * BCMP load-dependent weight table over the lattice:
     *   W(n) = (sum n)!/prod(n_r!) * prod_r D_r^{n_r} / prod_{k=1}^{sum n} beta(k),
     * beta(k) = min(k, c) for a c-server queue (c=1 -> single server). W(0)=1.
     */
    private static double[] ldTable(double[] Dq, double c, int[] shp, int total) {
        int R = shp.length;
        int cc = (int) c;
        if (!Double.isFinite(c) || cc <= 0) {
            cc = 1;
        }
        double[] W = new double[total];
        for (int i = 0; i < total; i++) {
            int[] n = decodeSub(i, shp);
            int tot = 0;
            for (int r = 0; r < R; r++) {
                tot += n[r];
            }
            double logf = lgamma(tot + 1);
            boolean ok = true;
            for (int r = 0; r < R; r++) {
                if (n[r] > 0) {
                    if (Dq[r] <= 0) {
                        ok = false;
                        break;
                    }
                    logf += n[r] * Math.log(Dq[r]) - lgamma(n[r] + 1);
                }
            }
            if (!ok) {
                continue;
            }
            for (int k = 1; k <= tot; k++) {
                logf -= Math.log(Math.min(k, cc));
            }
            W[i] = Math.exp(logf);
        }
        return W;
    }

    /** Lattice convolution Cv(m) = sum_{0<=a<=m} Av(a) Bv(m-a) over 0..N. */
    private static double[] convLattice(double[] Av, double[] Bv, int[] shp, int[] stride, int total) {
        int R = shp.length;
        int[][] subs = new int[total][R];
        for (int i = 0; i < total; i++) {
            subs[i] = decodeSub(i, shp);
        }
        double[] Cv = new double[total];
        for (int i = 0; i < total; i++) {
            int[] m = subs[i];
            double acc = 0.0;
            for (int j = 0; j <= i; j++) {
                int[] a = subs[j];
                boolean le = true;
                int idx = 0;
                for (int d = 0; d < R; d++) {
                    if (a[d] > m[d]) {
                        le = false;
                        break;
                    }
                    idx += (m[d] - a[d]) * stride[d];
                }
                if (le) {
                    acc += Av[j] * Bv[idx];
                }
            }
            Cv[i] = acc;
        }
        return Cv;
    }

    /**
     * G^{+} = sum_{0<=b<=N} Psi(b) G(N-b): the FNC of the target station
     * convolved against the full-network normalizing-constant table at n = N.
     */
    private static double fncMean(Matrix Psi, double[] Gfull, int[] shp, int[] stride, int total) {
        int R = shp.length;
        double val = 0.0;
        for (int i = 0; i < total; i++) {
            double psi = Psi.get(i, 0);
            if (psi == 0) {
                continue;
            }
            int[] b = decodeSub(i, shp);
            int idx = 0;
            for (int d = 0; d < R; d++) {
                idx += ((shp[d] - 1) - b[d]) * stride[d];
            }
            val += psi * Gfull[idx];
        }
        return val;
    }

    private static Matrix toMatrix(double[] v) {
        Matrix out = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) {
            out.set(i, 0, v[i]);
        }
        return out;
    }

    /**
     * Forward balanced-fairness fill of the OI balance function over the
     * lattice: Phi(0)=1, Phi(n) = (1/mu(n)) sum_{r: n_r>0} Phi(n - e_r).
     * Returned as a flat column-major (total x 1) Matrix.
     */
    private static Matrix oi_phi(ToDoubleFunction<int[]> oirate, int[] N, double[] vis) {
        int R = N.length;
        int[] shp = new int[R];
        int total = 1;
        for (int d = 0; d < R; d++) {
            shp[d] = N[d] + 1;
            total *= shp[d];
        }
        int[] stride = new int[R];
        stride[0] = 1;
        for (int d = 1; d < R; d++) {
            stride[d] = stride[d - 1] * shp[d - 1];
        }
        double[] Phiv = new double[total];
        for (int i = 0; i < total; i++) {
            int li = i;
            int[] n = new int[R];
            int sum = 0;
            for (int d = 0; d < R; d++) {
                n[d] = li % shp[d];
                li /= shp[d];
                sum += n[d];
            }
            if (sum == 0) {
                Phiv[i] = 1.0;
                continue;
            }
            double s = 0.0;
            for (int r = 0; r < R; r++) {
                if (n[r] > 0) {
                    double vr = (vis == null) ? 1.0 : vis[r];
                    s += vr * Phiv[i - stride[r]];
                }
            }
            Phiv[i] = s / oirate.applyAsDouble(n);
        }
        return toMatrix(Phiv);
    }

    private static boolean hasNonZero(Matrix m) {
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                if (m.get(i, j) != 0) {
                    return true;
                }
            }
        }
        return false;
    }

    private static double lgamma(double x) {
        // Lanczos approximation of ln Gamma(x), x > 0.
        double[] g = {
            676.5203681218851, -1259.1392167224028, 771.32342877765313,
            -176.61502916214059, 12.507343278686905, -0.13857109526572012,
            9.9843695780195716e-6, 1.5056327351493116e-7
        };
        if (x < 0.5) {
            return Math.log(Math.PI / Math.sin(Math.PI * x)) - lgamma(1 - x);
        }
        x -= 1;
        double a = 0.99999999999980993;
        double t = x + 7.5;
        for (int i = 0; i < g.length; i++) {
            a += g[i] / (x + i + 1);
        }
        return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
    }
}
