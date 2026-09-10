/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.analyzers;

import java.util.Map;
import java.util.function.DoubleUnaryOperator;

import jline.api.mam.Map_cdf;
import jline.api.qsys.Qsys_ggingi_tga;
import jline.api.qsys.Qsys_ggisgi_fluid;
import jline.api.qsys.Qsys_gtmtst_fluid;
import jline.api.qsys.Qsys_mtginf;
import jline.api.qsys.Qsys_mtgs0_mol;
import jline.api.qsys.QsysFluidAbandonResult;
import jline.api.qsys.QsysMtginfResult;
import jline.api.qsys.QsysTvFluidResult;
import jline.api.sn.SnArrivalRateFun;
import jline.api.sn.SnPatienceHandles;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * The single-station fluid limits: a Source -&gt; Queue -&gt; Sink model with one
 * class, answered by a closed-form fluid or Gaussian limit rather than by
 * integrating the network drift.
 *
 * <p>WHY THESE ARE FLUID METHODS AND NOT MVA ONES. Each depends on the service
 * or patience law BEYOND ITS MEAN -- the stationary point of the Liu-Whitt model
 * is where the patience ccdf crosses 1/rho, the Mt/G/inf mean is a convolution
 * with the service ccdf -- and each is the limit of a sequence of systems, not
 * an approximation to a fixed one. That is the fluid solver's contract.
 *
 * <p>Methods: {@code ggisgi.fluid} (Liu and Whitt, Operations Research 60(5),
 * 2012), {@code ggingi.tga} (Liu, Whitt and Yu, Naval Research Logistics 63(3),
 * 2016), {@code tvms} at constant staffing (Liu and Whitt, INFORMS J. Computing
 * 26(1), 2014), {@code mtginf} (Eick, Massey and Whitt, Management Science
 * 39(2), 1993) and {@code mol} (Massey and Whitt, Annals of Applied Probability
 * 4(4), 1994).
 *
 * <p>Java twin of {@code matlab/src/solvers/FLD/solver_fluid_qsys_analyzer.m}.
 *
 * @since LINE 3.1.0
 */
public class QsysLimitAnalyzer implements FluidAnalyzer {

    private Matrix xvecIt = Matrix.zeros(1, 1);
    private final String method;

    /**
     * @param method the concrete single-station limit to run
     */
    public QsysLimitAnalyzer(String method) {
        this.method = method;
    }

    /** The method names this analyzer answers. */
    public static boolean handles(String method) {
        return "ggisgi.fluid".equals(method) || "fluid.ggisgi".equals(method)
                || "ggisgi".equals(method) || "tga".equals(method)
                || "ggingi.tga".equals(method) || "fluid.tga".equals(method)
                || "tvms".equals(method) || "fluid.tvms".equals(method)
                || "mtginf".equals(method) || "fluid.mtginf".equals(method)
                || "mol".equals(method) || "fluid.mol".equals(method);
    }

    /** The canonical name of an alias. */
    public static String canonical(String method) {
        // the short spellings too, as the C++ fluid_qsys_canonical maps them
        if ("fluid.ggisgi".equals(method) || "ggisgi".equals(method)) return "ggisgi.fluid";
        if ("fluid.tga".equals(method) || "tga".equals(method)) return "ggingi.tga";
        if ("fluid.tvms".equals(method)) return "tvms";
        if ("fluid.mtginf".equals(method)) return "mtginf";
        if ("fluid.mol".equals(method)) return "mol";
        return method;
    }

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {
        int M = sn.nstations;
        int K = sn.nclasses;
        result.QN = Matrix.zeros(M, K);
        result.UN = Matrix.zeros(M, K);
        result.RN = Matrix.zeros(M, K);
        result.TN = Matrix.zeros(M, K);
        result.AN = Matrix.zeros(M, K);
        result.CN = Matrix.zeros(1, K);
        result.XN = Matrix.zeros(1, K);
        result.WN = new Matrix(0, 0);

        int src = -1;
        int qi = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                src = (int) sn.nodeToStation.get(i);
            } else if (sn.nodetype.get(i) == NodeType.Queue || sn.nodetype.get(i) == NodeType.Delay) {
                qi = (int) sn.nodeToStation.get(i);
            }
        }
        // THE SHAPE THESE LIMITS ARE STATED FOR, refused by name rather than
        // answered on a model they do not describe: one open class through one
        // queueing station. The MVA qsys analyzer is reached by a structural
        // dispatch that guarantees it; these methods are selected by NAME, so
        // the check has to live here.
        if (src < 0 || qi < 0) {
            throw new RuntimeException("the single-station fluid limits need a Source and a queueing station");
        }
        if (K != 1 || sn.nclosedjobs > 0) {
            throw new RuntimeException("the '" + this.method + "' method is a single-station limit: "
                    + "it needs one open class through one Source and one queueing station");
        }

        int statefulIdx = (int) sn.stationToStateful.get(qi);
        double Vq = sn.visits.get(0).get(statefulIdx);
        double lambda = sn.rates.get(src, 0) * Vq;
        double mu = sn.rates.get(qi, 0);
        double nserv = sn.nservers.get(qi);
        double scvS = sn.scv.get(qi, 0);
        double ca = Math.sqrt(sn.scv.get(src, 0));
        double cs = Math.sqrt(scvS);
        SnPatienceHandles.Handles h = SnPatienceHandles.snPatienceHandles(sn, qi, 0);

        // The service ccdf, needed by the two Mt/G methods: they are exact in
        // the service DISTRIBUTION, not in its mean, which is the whole point
        // of the Eick-Massey-Whitt lag.
        final Matrix D0;
        final Matrix D1;
        Station qStation = sn.stations.get(qi);
        JobClass jobClass = sn.jobclasses.get(0);
        Map<JobClass, MatrixCell> procMap = sn.proc == null ? null : sn.proc.get(qStation);
        MatrixCell svc = procMap == null ? null : procMap.get(jobClass);
        if (svc != null && svc.size() >= 2 && svc.get(0) != null && svc.get(1) != null
                && svc.get(0).getNumRows() == svc.get(1).getNumRows()) {
            D0 = svc.get(0);
            D1 = svc.get(1);
        } else {
            D0 = null;
            D1 = null;
        }
        final double muF = mu;
        DoubleUnaryOperator serviceCcdf = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double x) {
                if (D0 == null) {
                    return Math.exp(-muF * x);
                }
                Matrix pt = new Matrix(1, 1);
                pt.set(0, 0, x);
                return 1.0 - Map_cdf.map_cdf(D0, D1, pt).get(0);
            }
        };
        double ES = 1.0 / mu;
        double ES2 = (1.0 + scvS) * ES * ES;

        String m = canonical(this.method);
        if ("ggisgi.fluid".equals(m)) {
            requirePatience(h, m);
            QsysFluidAbandonResult res = Qsys_ggisgi_fluid.qsys_ggisgi_fluid(
                    lambda, mu, (int) Math.round(nserv), h.ccdf);
            stationary(result, qi, src, Vq, lambda, res.meanNumber, res.throughput, res.utilization);
        } else if ("ggingi.tga".equals(m)) {
            requirePatience(h, m);
            requireFiniteServers(nserv, m);
            Map<String, Double> res = Qsys_ggingi_tga.qsys_ggingi_tga(
                    lambda, mu, (int) Math.round(nserv), ca, cs, h.ccdf, h.pdf, serviceCcdf);
            double Tq = lambda * (1.0 - res.get("probAbandon"));
            stationary(result, qi, src, Vq, lambda, res.get("meanNumber"), Tq,
                    Math.min(res.get("meanNumberInService") / nserv, 1.0));
        } else if ("tvms".equals(m)) {
            requirePatience(h, m);
            requireFiniteServers(nserv, m);
            SnArrivalRateFun.RateFun rf = SnArrivalRateFun.snArrivalRateFun(sn, src, 0);
            double[] window = horizon(options);
            final double sFixed = nserv;
            // CONSTANT STAFFING. Nothing in a Network declares a time-varying
            // server count, so s(t) is the station's own s; the time variation
            // the method is for enters through lambda(t) alone. A staffing
            // schedule would need a model feature that does not exist, and
            // inventing one here would make the solver answer a model the user
            // did not build.
            QsysTvFluidResult res = Qsys_gtmtst_fluid.qsys_gtmtst_fluid(
                    rf.lambda,
                    new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double t) { return sFixed; }
                    },
                    new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double t) { return muF; }
                    },
                    h.ccdf, window[1] - window[0], Double.NaN, 0.0, 0.0, null, h.pdf, null);
            double[] times = shift(res.times, window[0]);
            double[] served = scale(res.B, mu);
            transient_(result, sn, qi, src, Vq, times, res.X, res.utilization, served, res.arrivalRate);
        } else if ("mtginf".equals(m)) {
            SnArrivalRateFun.RateFun rf = SnArrivalRateFun.snArrivalRateFun(sn, src, 0);
            double[] window = horizon(options);
            double[] tvals = linspace(window[0], window[1], 200);
            QsysMtginfResult res = Qsys_mtginf.qsys_mtginf(rf.lambda, serviceCcdf, ES, tvals,
                    Double.NEGATIVE_INFINITY, ES2, null, Qsys_mtginf.DEFAULT_TOL,
                    Qsys_mtginf.DEFAULT_PANELS, 1e12);
            // An infinite-server station serves everything that arrives, so the
            // throughput is the arrival rate and the busy-server count is what
            // a utilization column can carry.
            transient_(result, sn, qi, src, Vq, res.times, res.meanNumber, res.meanNumber,
                    res.arrivalRate, res.arrivalRate);
        } else if ("mol".equals(m)) {
            requireFiniteServers(nserv, m);
            SnArrivalRateFun.RateFun rf = SnArrivalRateFun.snArrivalRateFun(sn, src, 0);
            double[] window = horizon(options);
            double[] tvals = linspace(window[0], window[1], 200);
            double cap = sn.cap == null ? Double.POSITIVE_INFINITY : sn.cap.get(qi);
            // An uncapped station carries Integer.MAX_VALUE here, not Inf as in
            // MATLAB and Python; a finite buffer beyond the servers is not part
            // of the loss model the approximation is for.
            boolean unbounded = !Double.isFinite(cap) || cap >= Integer.MAX_VALUE;
            boolean useDelay = !unbounded && cap > nserv;
            Map<String, double[]> res = Qsys_mtgs0_mol.qsys_mtgs0_mol(rf.lambda, serviceCcdf, ES,
                    (int) Math.round(nserv), tvals, Double.NEGATIVE_INFINITY, useDelay);
            double[] busy = res.get("meanBusyMOL");
            transient_(result, sn, qi, src, Vq, res.get("times"), busy, scale(busy, 1.0 / nserv),
                    scale(busy, mu), res.get("arrivalRate"));
        } else {
            throw new RuntimeException("the '" + this.method + "' method is not a single-station fluid limit");
        }
        result.method = m;
        xvecIt = Matrix.zeros(1, 1);
    }

    @Override
    public Matrix getXVecIt() {
        return xvecIt;
    }

    private static void requirePatience(SnPatienceHandles.Handles h, String method) {
        if (h == null) {
            throw new RuntimeException("the '" + method + "' method needs a reneging patience law "
                    + "on the queue (Queue.setPatience)");
        }
    }

    private static void requireFiniteServers(double nserv, String method) {
        if (!Double.isFinite(nserv) || nserv < 1) {
            throw new RuntimeException("the '" + method + "' method needs a finite number of servers");
        }
    }

    /**
     * The integration window of a time-varying single-station fluid limit,
     * asked as a predicate rather than thrown.
     *
     * <p>The time-varying limits ("mol", "mtginf", "tvms") report a TRAJECTORY,
     * so an infinite or absent upper end of options.timespan leaves them
     * nothing to report.</p>
     *
     * <p>A horizon is a solver OPTION and not a model feature, so the feature
     * registry has no name for it and {@code SolverFluid.supportsModelMethod}
     * has to ask this predicate directly. {@link #horizon} asks the same one on
     * the solve path, which is what keeps the report and the run from
     * disagreeing about whether a method can be asked for.</p>
     *
     * @param options solver options carrying the timespan
     * @return empty string when the window is a finite non-empty interval, else the refusal
     */
    public static String horizonReason(SolverOptions options) {
        double t0 = 0.0;
        double t1 = 1.0;
        if (options != null && options.timespan != null && options.timespan.length >= 2) {
            if (Double.isFinite(options.timespan[0])) {
                t0 = options.timespan[0];
            }
            t1 = options.timespan[1];
        }
        if (!Double.isFinite(t1) || t1 <= t0) {
            return "A time-varying fluid method needs a finite horizon: set "
                    + "options.timespan = [t0, t1].";
        }
        return "";
    }

    /** The integration window, refused by name when it is not a finite interval. */
    private static double[] horizon(SolverOptions options) {
        String reason = horizonReason(options);
        if (!reason.isEmpty()) {
            throw new RuntimeException(reason);
        }
        double t0 = 0.0;
        if (options != null && options.timespan != null && options.timespan.length >= 2
                && Double.isFinite(options.timespan[0])) {
            t0 = options.timespan[0];
        }
        double t1 = (options != null && options.timespan != null && options.timespan.length >= 2)
                ? options.timespan[1] : 1.0;
        return new double[]{t0, t1};
    }

    /** Little's law on the CARRIED rate, as every LINE solver reports a station that loses work. */
    private static void stationary(SolverResult result, int qi, int src, double Vq, double lambda,
                                   double Lsys, double Tq, double Uq) {
        double R = Tq > 0 ? Lsys / Tq : 0.0;
        result.RN.set(qi, 0, R);
        result.QN.set(qi, 0, Lsys);
        result.UN.set(qi, 0, Uq);
        result.TN.set(qi, 0, Tq);
        result.TN.set(src, 0, lambda / Vq);
        // The OFFERED rate, so that the loss table reads the abandonment as
        // ArvR - Tput.
        result.AN.set(qi, 0, lambda);
        result.XN.set(0, 0, Tq);
        result.CN.set(0, 0, R * Vq);
    }

    /**
     * The steady-state row of a time-varying model is the TIME AVERAGE over the
     * horizon, which is what a stationary reader of a periodic system measures;
     * the trajectory itself is returned beside it.
     */
    private static void transient_(SolverResult result, NetworkStruct sn, int qi, int src, double Vq,
                                   double[] t, double[] Lt, double[] Ut, double[] Tt,
                                   double[] arrival) {
        double span = t[t.length - 1] - t[0];
        double Lbar;
        double Ubar;
        double Tbar;
        double Abar;
        if (span <= 0) {
            Lbar = Lt[0];
            Ubar = Ut[0];
            Tbar = Tt[0];
            Abar = arrival[0];
        } else {
            Lbar = trapz(t, Lt) / span;
            Ubar = trapz(t, Ut) / span;
            Tbar = trapz(t, Tt) / span;
            Abar = trapz(t, arrival) / span;
        }
        result.QN.set(qi, 0, Lbar);
        result.UN.set(qi, 0, Ubar);
        result.TN.set(qi, 0, Tbar);
        result.TN.set(src, 0, Abar / Vq);
        result.RN.set(qi, 0, Tbar > 0 ? Lbar / Tbar : 0.0);
        result.AN.set(qi, 0, Abar);
        result.XN.set(0, 0, Tbar);
        result.CN.set(0, 0, result.RN.get(qi, 0) * Vq);

        int M = sn.nstations;
        int K = sn.nclasses;
        result.t = new Matrix(t.length, 1);
        for (int i = 0; i < t.length; i++) {
            result.t.set(i, 0, t[i]);
        }
        result.QNt = new Matrix[M][K];
        result.UNt = new Matrix[M][K];
        result.TNt = new Matrix[M][K];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                result.QNt[i][r] = Matrix.zeros(t.length, 1);
                result.UNt[i][r] = Matrix.zeros(t.length, 1);
                result.TNt[i][r] = Matrix.zeros(t.length, 1);
            }
        }
        for (int i = 0; i < t.length; i++) {
            result.QNt[qi][0].set(i, 0, Lt[i]);
            result.UNt[qi][0].set(i, 0, Ut[i]);
            result.TNt[qi][0].set(i, 0, Tt[i]);
            result.TNt[src][0].set(i, 0, arrival[i] / Vq);
        }
    }

    private static double trapz(double[] x, double[] y) {
        double s = 0;
        for (int i = 1; i < x.length; i++) {
            s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
        }
        return s;
    }

    private static double[] linspace(double a, double b, int n) {
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            out[i] = a + (b - a) * i / (n - 1);
        }
        return out;
    }

    private static double[] shift(double[] v, double by) {
        double[] out = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            out[i] = v[i] + by;
        }
        return out;
    }

    private static double[] scale(double[] v, double by) {
        double[] out = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            out[i] = v[i] * by;
        }
        return out;
    }
}
