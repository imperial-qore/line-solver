/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.analyzers;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.MAPt;
import jline.lang.processes.PHt;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Fluid and diffusion limits of the (MAP_t/Ph_t/inf)^N network of Y. M. Ko and J. Pender,
 * "Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett. 45 (2017)
 * 248-253.
 *
 * <p>The mean and the covariance of the limit are integrated jointly:
 *
 * <pre>
 *   dq/dt     = F(t,q) = A f(t,q)
 *   dSigma/dt = J Sigma + Sigma J' + G,   J = A df/dq,  G = A diag(f) A'
 * </pre>
 *
 * <p>with A the jump matrix whose column e is the jump vector of event e and f the event rate
 * vector. G is exactly dH dH' of Theorem 3.3, each independent Poisson term contributing
 * l_e l_e' f_e. Where f is affine in q -- infinite-server stations and the arrival phase
 * process -- J does not depend on q and both equations close exactly, so for the
 * (MAP_t/Ph_t/inf)^N case the mean and the covariance are exact rather than asymptotic.
 * Finite-server stations are admitted through the usual fluid min(x,c) term, where the
 * covariance degrades to a linear-noise approximation.
 *
 * <p>This analyzer does NOT reuse the closing ODE. That formulation routes a departure from the
 * source to the destination station and returns mass through the STATIONARY arrival-instant
 * vector pie, replacing the D1' operator by the rank-one map pie*(D1*e)', i.e. by the PH renewal
 * process with representation (pie, D0). Its stationary arrival rate is exact but its
 * autocorrelation is gone, and a non-renewal arrival stream is the entire point of a MAP.
 *
 * <p>State layout, station-major, arrival phases before service phases: one u-block per
 * (EXT station, class) holding the arrival MAP phase occupancy, which sums to 1, and one x-block
 * per (queueing station, class) holding the fluid count in each service phase.
 */
public class KoPenderAnalyzer implements FluidAnalyzer {

    /** kind codes of the five Ko-Pender event families. */
    private static final int A0 = 1;
    private static final int A1 = 2;
    private static final int SVC = 3;
    private static final int DEP = 4;
    private static final int ROUTE = 5;

    private Matrix xvec_it;

    /** One state block: a (station, class) pair with its offset and phase count. */
    private static final class Block {
        final int station;
        final int jobclass;
        final int offset;
        final int phases;

        Block(int station, int jobclass, int offset, int phases) {
            this.station = station;
            this.jobclass = jobclass;
            this.offset = offset;
            this.phases = phases;
        }
    }

    /** One event: its jump column and the descriptor from which its rate is assembled. */
    private static final class Event {
        final int kind;
        final int station;
        final int jobclass;
        final int fromPhase;
        final int toPhase;
        final int destStation;
        final int destClass;
        final int destPhase;
        final double[] jump;

        Event(int kind, int station, int jobclass, int fromPhase, int toPhase,
              int destStation, int destClass, int destPhase, double[] jump) {
            this.kind = kind;
            this.station = station;
            this.jobclass = jobclass;
            this.fromPhase = fromPhase;
            this.toPhase = toPhase;
            this.destStation = destStation;
            this.destClass = destClass;
            this.destPhase = destPhase;
            this.jump = jump;
        }
    }

    /** A MAPt or PHt attached to a (station, class), unpacked into per-segment MAP pairs. */
    private static final class Schedule {
        final int station;
        final int jobclass;
        final double[] breakpoints;
        final boolean cyclic;
        final List<Matrix> d0;
        final List<Matrix> d1;

        Schedule(int station, int jobclass, double[] breakpoints, boolean cyclic,
                 List<Matrix> d0, List<Matrix> d1) {
            this.station = station;
            this.jobclass = jobclass;
            this.breakpoints = breakpoints;
            this.cyclic = cyclic;
            this.d0 = d0;
            this.d1 = d1;
        }

        int segmentAt(double t) {
            double period = breakpoints[breakpoints.length - 1] - breakpoints[0];
            double offset = t - breakpoints[0];
            if (cyclic) {
                offset = offset % period;
                if (offset < 0) {
                    offset += period;
                }
            } else if (offset < 0.0 || offset >= period) {
                return -1;
            }
            double pos = breakpoints[0] + offset;
            for (int k = 0; k < d0.size(); k++) {
                if (pos < breakpoints[k + 1]) {
                    return k;
                }
            }
            return d0.size() - 1;
        }
    }

    private NetworkStruct sn;
    private List<Block> ublocks;
    private List<Block> xblocks;
    private List<Event> events;
    private List<Schedule> schedules;
    private Matrix[][] nomD0;
    private Matrix[][] nomD1;
    private Matrix[][] nomPie;
    private int dim;

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {
        this.sn = sn;
        int M = sn.nstations;
        int K = sn.nclasses;

        for (int r = 0; r < K; r++) {
            if (Double.isFinite(sn.njobs.get(0, r))) {
                throw new RuntimeException(
                        "the 'kp' method analyses the open (MAP_t/Ph_t/inf)^N network of Ko and "
                                + "Pender (2017); a closed class has no arrival process to "
                                + "modulate. Use 'closing' or 'matrix' for closed models.");
            }
        }

        buildNominals(M, K);
        buildBlocks(M, K);
        if (ublocks.isEmpty()) {
            throw new RuntimeException(
                    "the 'kp' method needs at least one Source with an arrival process");
        }
        buildSchedules(M, K);
        buildEvents(K);

        double t0 = options.timespan[0];
        double tend = options.timespan[1];
        boolean unbounded = !Double.isFinite(tend);
        if (!Double.isFinite(t0)) {
            t0 = 0.0;
        }
        double period = 0.0;
        for (int e = 0; e < schedules.size(); e++) {
            Schedule sc = schedules.get(e);
            if (sc.cyclic) {
                period = Math.max(period,
                        sc.breakpoints[sc.breakpoints.length - 1] - sc.breakpoints[0]);
            }
        }
        if (unbounded) {
            double slow = Double.POSITIVE_INFINITY;
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    double rate = sn.rates.get(i, r);
                    if (Double.isFinite(rate) && rate > 0 && rate < slow) {
                        slow = rate;
                    }
                }
            }
            if (!Double.isFinite(slow) || slow <= 0) {
                slow = 1.0;
            }
            tend = t0 + Math.max(10.0, 30.0 / slow);
            if (period > 0) {
                tend = Math.max(tend, t0 + 10.0 * period);
            }
        }

        double[] q0 = new double[dim];
        double[][] sigma0 = new double[dim][dim];
        for (int b = 0; b < ublocks.size(); b++) {
            Block blk = ublocks.get(b);
            double[] theta = stationaryPhase(pairAt(blk.station, blk.jobclass, t0)[0],
                    pairAt(blk.station, blk.jobclass, t0)[1]);
            for (int i = 0; i < blk.phases; i++) {
                q0[blk.offset + i] = theta[i];
            }
            // The arrival phase is drawn from the stationary vector rather than known, so its
            // indicator carries covariance diag(u0) - u0 u0'. Zero would assert a known initial
            // phase and understate the variance until the initial condition has washed out.
            for (int i = 0; i < blk.phases; i++) {
                for (int j = 0; j < blk.phases; j++) {
                    sigma0[blk.offset + i][blk.offset + j] =
                            (i == j ? theta[i] : 0.0) - theta[i] * theta[j];
                }
            }
        }

        // NOT options.init_sol: that one is laid out for the CLOSING state vector, and this
        // method has its own layout (u-blocks before x-blocks, no returning mass), so consuming
        // it would silently zero the source phase mass and the whole network.
        //
        // A WRONG-SIZED SEED IS REFUSED, not ignored. Dropping it would integrate from the
        // default initial condition under the caller's name and return a plausible trajectory
        // for a model the caller did not ask about.
        if (options.config != null && options.config.kp_init_sol != null
                && options.config.kp_init_sol.length > 0) {
            if (options.config.kp_init_sol.length != dim) {
                throw new RuntimeException(String.format(
                        "config.kp_init_sol has %d entries but the 'kp' state vector of this model "
                                + "has %d, laid out station-major over the (station, class) blocks. "
                                + "It is NOT laid out like options.init_sol.",
                        options.config.kp_init_sol.length, dim));
            }
            System.arraycopy(options.config.kp_init_sol, 0, q0, 0, dim);
        }
        // Companion seed for the covariance. A caller that carries a DISTRIBUTION across a
        // handoff supplies the second moment beside the mean, so the next stage does not restart
        // from a point mass it never had. Same layout as kp_init_sol.
        if (options.config != null && options.config.init_cov != null
                && !options.config.init_cov.isEmpty()) {
            Matrix cov = options.config.init_cov;
            if (cov.getNumRows() != dim || cov.getNumCols() != dim) {
                throw new RuntimeException(String.format(
                        "config.init_cov is %dx%d but the 'kp' state vector of this model has %d "
                                + "entries, so the covariance must be %dx%d.",
                        cov.getNumRows(), cov.getNumCols(), dim, dim, dim));
            }
            double asym = 0.0;
            double scale = 0.0;
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    double d = cov.get(i, j) - cov.get(j, i);
                    asym += d * d;
                    scale += cov.get(i, j) * cov.get(i, j);
                }
            }
            // Loose enough for the rounding of a covariance that was itself integrated, tight
            // enough to catch a matrix that is simply not one.
            if (Math.sqrt(asym) > 1e-6 * Math.max(1.0, Math.sqrt(scale))) {
                throw new RuntimeException("config.init_cov must be symmetric.");
            }
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    sigma0[i][j] = cov.get(i, j);
                }
            }
        }

        double[] z0 = new double[dim + dim * dim];
        System.arraycopy(q0, 0, z0, 0, dim);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                z0[dim + i * dim + j] = sigma0[i][j];
            }
        }

        double maxStep = (tend - t0) / 10.0;
        if (period > 0) {
            double narrowest = Double.POSITIVE_INFINITY;
            for (int e = 0; e < schedules.size(); e++) {
                double[] bp = schedules.get(e).breakpoints;
                for (int k = 0; k + 1 < bp.length; k++) {
                    narrowest = Math.min(narrowest, bp[k + 1] - bp[k]);
                }
            }
            if (Double.isFinite(narrowest)) {
                maxStep = Math.min(maxStep, narrowest / 4.0);
            }
        }

        final double horizonEnd = tend;
        final List<double[]> traj = new ArrayList<double[]>();
        final List<Double> times = new ArrayList<Double>();
        // The period average below is a QUADRATURE, and the integrator picks its grid for
        // accuracy of the SOLUTION: at a tight tolerance it takes long steps through the
        // smooth stretches and the trapezoid over them loses more than the integration
        // gained. Force output at a uniform mesh over the averaging window, and bracket each
        // schedule boundary so that no trapezoid interval straddles the jump in the arrival
        // rate. The interpolator is asked for those points explicitly because a step handler
        // only sees the steps the integrator chose.
        final double[] mesh = averagingMesh(t0, tend, period, unbounded);
        FirstOrderDifferentialEquations ode = new AugmentedODE();
        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                1e-10, maxStep, options.tol * 1e-3, options.tol);
        integrator.addStepHandler(new org.apache.commons.math3.ode.sampling.StepHandler() {
            @Override
            public void init(double t0i, double[] y0i, double tf) {
                times.add(t0i);
                traj.add(y0i.clone());
            }

            @Override
            public void handleStep(org.apache.commons.math3.ode.sampling.StepInterpolator interp,
                                   boolean isLast) {
                double prev = interp.getPreviousTime();
                double cur = interp.getCurrentTime();
                for (int i = 0; i < mesh.length; i++) {
                    double tm = mesh[i];
                    if (tm > prev && tm < cur) {
                        interp.setInterpolatedTime(tm);
                        times.add(tm);
                        traj.add(interp.getInterpolatedState().clone());
                    }
                }
                interp.setInterpolatedTime(cur);
                times.add(cur);
                traj.add(interp.getInterpolatedState().clone());
            }
        });
        double[] zEnd = new double[z0.length];
        integrator.integrate(ode, t0, z0, horizonEnd, zEnd);

        int nt = times.size();
        Matrix tMat = new Matrix(nt, 1);
        for (int n = 0; n < nt; n++) {
            tMat.set(n, 0, times.get(n));
        }

        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix RN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);
        Matrix[][] QNt = new Matrix[M][K];
        Matrix[][] UNt = new Matrix[M][K];
        Matrix[][] TNt = new Matrix[M][K];
        Matrix[][] QVart = new Matrix[M][K];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                QNt[i][r] = new Matrix(nt, 1);
                UNt[i][r] = new Matrix(nt, 1);
                TNt[i][r] = new Matrix(nt, 1);
                QVart[i][r] = new Matrix(nt, 1);
            }
        }

        // A cyclic schedule has no fixed point, so a steady-state request is answered by the
        // time average over the last full period; the value at tend would be an arbitrary point
        // of the cycle and source and station throughput would disagree.
        double windowLo = (unbounded && period > 0) ? Math.max(t0, tend - period) : Double.NaN;

        for (int b = 0; b < xblocks.size(); b++) {
            Block blk = xblocks.get(b);
            double[] qser = new double[nt];
            double[] vser = new double[nt];
            double[] tser = new double[nt];
            for (int n = 0; n < nt; n++) {
                double[] z = traj.get(n);
                double q = 0.0;
                for (int i = 0; i < blk.phases; i++) {
                    q += z[blk.offset + i];
                }
                double v = 0.0;
                for (int i = 0; i < blk.phases; i++) {
                    for (int j = 0; j < blk.phases; j++) {
                        v += z[dim + (blk.offset + i) * dim + (blk.offset + j)];
                    }
                }
                Matrix d1 = pairAt(blk.station, blk.jobclass, times.get(n))[1];
                double scale = capacityScale(z, blk.station);
                double tp = 0.0;
                for (int i = 0; i < blk.phases; i++) {
                    double rowsum = 0.0;
                    for (int j = 0; j < d1.getNumCols(); j++) {
                        rowsum += d1.get(i, j);
                    }
                    tp += rowsum * Math.max(z[blk.offset + i], 0.0);
                }
                qser[n] = q;
                vser[n] = v;
                tser[n] = tp * scale;
                QNt[blk.station][blk.jobclass].set(n, 0, q);
                QVart[blk.station][blk.jobclass].set(n, 0, v);
                TNt[blk.station][blk.jobclass].set(n, 0, tser[n]);
                double util;
                if (sn.sched.get(sn.stations.get(blk.station)) == SchedStrategy.INF
                        || !Double.isFinite(sn.nservers.get(blk.station, 0))) {
                    util = q;
                } else {
                    util = Math.min(q, sn.nservers.get(blk.station, 0))
                            / sn.nservers.get(blk.station, 0);
                }
                UNt[blk.station][blk.jobclass].set(n, 0, util);
            }
            QN.set(blk.station, blk.jobclass, summarise(qser, times, windowLo, tend));
            double[] userArr = new double[nt];
            for (int n = 0; n < nt; n++) {
                userArr[n] = UNt[blk.station][blk.jobclass].get(n, 0);
            }
            UN.set(blk.station, blk.jobclass, summarise(userArr, times, windowLo, tend));
            TN.set(blk.station, blk.jobclass, summarise(tser, times, windowLo, tend));
            double tput = TN.get(blk.station, blk.jobclass);
            RN.set(blk.station, blk.jobclass,
                    tput > 0 ? QN.get(blk.station, blk.jobclass) / tput : 0.0);
        }

        for (int b = 0; b < ublocks.size(); b++) {
            Block blk = ublocks.get(b);
            double[] aser = new double[nt];
            for (int n = 0; n < nt; n++) {
                double[] z = traj.get(n);
                Matrix d1 = pairAt(blk.station, blk.jobclass, times.get(n))[1];
                double arr = 0.0;
                for (int i = 0; i < blk.phases; i++) {
                    double rowsum = 0.0;
                    for (int j = 0; j < d1.getNumCols(); j++) {
                        rowsum += d1.get(i, j);
                    }
                    arr += rowsum * Math.max(z[blk.offset + i], 0.0);
                }
                aser[n] = arr;
                TNt[blk.station][blk.jobclass].set(n, 0, arr);
            }
            TN.set(blk.station, blk.jobclass, summarise(aser, times, windowLo, tend));
        }

        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.QNt = QNt;
        result.UNt = UNt;
        result.TNt = TNt;
        result.t = tMat;
        // The covariance is the point of this method, so it must reach the caller rather
        // than stay a local: FluidResult carries the per-block variance and the full Sigma(t).
        if (result instanceof jline.solvers.fluid.FluidResult) {
            jline.solvers.fluid.FluidResult fr = (jline.solvers.fluid.FluidResult) result;
            fr.QVart = QVart;
            Matrix[] sigmaSeries = new Matrix[nt];
            for (int n = 0; n < nt; n++) {
                Matrix sig = new Matrix(dim, dim);
                double[] z = traj.get(n);
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        sig.set(i, j, z[dim + i * dim + j]);
                    }
                }
                sigmaSeries[n] = sig;
            }
            fr.Sigmat = sigmaSeries;
        }

        this.xvec_it = new Matrix(1, dim);
        double[] zLast = traj.get(nt - 1);
        for (int j = 0; j < dim; j++) {
            this.xvec_it.set(0, j, zLast[j]);
        }
    }

    /**
     * Extra output times over the averaging window: a uniform mesh plus a bracket around each
     * schedule boundary. Empty when the run is not answering a steady-state request for a
     * cyclic schedule, in which case no period average is taken.
     */
    private double[] averagingMesh(double t0, double tend, double period, boolean unbounded) {
        if (!unbounded || !(period > 0)) {
            return new double[0];
        }
        double lo = Math.max(t0, tend - period);
        int uniform = 2001;
        java.util.TreeSet<Double> pts = new java.util.TreeSet<Double>();
        for (int i = 0; i < uniform; i++) {
            pts.add(lo + (tend - lo) * i / (uniform - 1.0));
        }
        double eps = Math.max(1e-9, 1e-7 * (tend - lo));
        for (int e = 0; e < schedules.size(); e++) {
            Schedule sc = schedules.get(e);
            double per = sc.breakpoints[sc.breakpoints.length - 1] - sc.breakpoints[0];
            int kmax = (per > 0) ? (int) Math.ceil((tend - lo) / per) + 2 : 0;
            for (int k = -1; k <= kmax; k++) {
                for (int b = 0; b < sc.breakpoints.length; b++) {
                    double bp = sc.breakpoints[b] + (sc.cyclic ? k * per : 0.0);
                    for (double cand : new double[]{bp - eps, bp, bp + eps}) {
                        if (cand > lo && cand < tend) {
                            pts.add(cand);
                        }
                    }
                    if (!sc.cyclic) {
                        break;
                    }
                }
            }
        }
        double[] out = new double[pts.size()];
        int i = 0;
        for (Double d : pts) {
            out[i++] = d;
        }
        return out;
    }

    @Override
    public Matrix getXVecIt() {
        return xvec_it;
    }

    // ------------------------------------------------------------------ structure

    private void buildNominals(int M, int K) {
        nomD0 = new Matrix[M][K];
        nomD1 = new Matrix[M][K];
        nomPie = new Matrix[M][K];
        for (int i = 0; i < M; i++) {
            Station st = sn.stations.get(i);
            for (int r = 0; r < K; r++) {
                JobClass jc = sn.jobclasses.get(r);
                MatrixCell slot = sn.proc.get(st).get(jc);
                ProcessType pt = sn.procid.get(st).get(jc);
                if (pt == ProcessType.MAPT || pt == ProcessType.PHT) {
                    MatrixCell nominal = nominalFromSchedule(slot, pt == ProcessType.MAPT);
                    nomD0[i][r] = nominal.get(0);
                    nomD1[i][r] = nominal.get(1);
                } else if (slot != null && slot.size() >= 2
                        && slot.get(0).getNumRows() == slot.get(1).getNumRows()
                        && slot.get(0).getNumCols() == slot.get(1).getNumCols()) {
                    nomD0[i][r] = slot.get(0);
                    nomD1[i][r] = slot.get(1);
                } else {
                    double lam = sn.rates.get(i, r);
                    Matrix d0 = new Matrix(1, 1, 1);
                    Matrix d1 = new Matrix(1, 1, 1);
                    d0.set(0, 0, -lam);
                    d1.set(0, 0, lam);
                    nomD0[i][r] = d0;
                    nomD1[i][r] = d1;
                }
                nomPie[i][r] = jline.api.mam.Map_pie.map_pie(nomD0[i][r], nomD1[i][r]);
            }
        }
    }

    /**
     * Width-weighted time-averaged (D0, D1) of a MAPt or PHt slot of sn.proc.
     *
     * <p>The slot is flat, [breakpoints, A_1..A_n, B_1..B_n, cyclic]; for a PHt the pairs are
     * (alpha, S) whose equivalent MAP pair is (S, s*alpha).
     */
    static MatrixCell nominalFromSchedule(MatrixCell slot, boolean isMapt) {
        int n = (slot.size() - 2) / 2;
        Matrix bp = slot.get(0);
        int h = slot.get(1 + n).getNumRows();
        Matrix d0bar = new Matrix(h, h);
        Matrix d1bar = new Matrix(h, h);
        double total = bp.get(0, n) - bp.get(0, 0);
        for (int k = 0; k < n; k++) {
            double w = (bp.get(0, k + 1) - bp.get(0, k)) / total;
            Matrix[] pair = segmentPair(slot, isMapt, k, n);
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < h; j++) {
                    d0bar.set(i, j, d0bar.get(i, j) + w * pair[0].get(i, j));
                    d1bar.set(i, j, d1bar.get(i, j) + w * pair[1].get(i, j));
                }
            }
        }
        MatrixCell out = new MatrixCell();
        out.set(0, d0bar);
        out.set(1, d1bar);
        return out;
    }

    static Matrix[] segmentPair(MatrixCell slot, boolean isMapt, int k, int n) {
        if (isMapt) {
            return new Matrix[]{slot.get(1 + k), slot.get(1 + n + k)};
        }
        Matrix alpha = slot.get(1 + k);
        Matrix s = slot.get(1 + n + k);
        int h = s.getNumRows();
        Matrix d1 = new Matrix(h, h);
        for (int i = 0; i < h; i++) {
            double rowsum = 0.0;
            for (int j = 0; j < h; j++) {
                rowsum += s.get(i, j);
            }
            for (int j = 0; j < h; j++) {
                d1.set(i, j, -rowsum * alpha.get(0, j));
            }
        }
        return new Matrix[]{s, d1};
    }

    private void buildBlocks(int M, int K) {
        ublocks = new ArrayList<Block>();
        xblocks = new ArrayList<Block>();
        int off = 0;
        for (int i = 0; i < M; i++) {
            boolean isExt = sn.sched.get(sn.stations.get(i)) == SchedStrategy.EXT;
            for (int r = 0; r < K; r++) {
                int h = (int) sn.phases.get(i, r);
                double rate = sn.rates.get(i, r);
                if (h <= 0 || !Double.isFinite(rate) || rate <= 0) {
                    continue; // disabled class at this station
                }
                Block blk = new Block(i, r, off, h);
                if (isExt) {
                    ublocks.add(blk);
                } else {
                    xblocks.add(blk);
                }
                off += h;
            }
        }
        dim = off;
    }

    private void buildSchedules(int M, int K) {
        schedules = new ArrayList<Schedule>();
        for (int i = 0; i < M; i++) {
            Station st = sn.stations.get(i);
            for (int r = 0; r < K; r++) {
                JobClass jc = sn.jobclasses.get(r);
                ProcessType pt = sn.procid.get(st).get(jc);
                if (pt != ProcessType.MAPT && pt != ProcessType.PHT) {
                    continue;
                }
                MatrixCell slot = sn.proc.get(st).get(jc);
                int n = (slot.size() - 2) / 2;
                Matrix bpMat = slot.get(0);
                double[] bp = new double[n + 1];
                for (int k = 0; k <= n; k++) {
                    bp[k] = bpMat.get(0, k);
                }
                boolean cyclic = slot.get(slot.size() - 1).get(0, 0) != 0.0;
                List<Matrix> d0 = new ArrayList<Matrix>();
                List<Matrix> d1 = new ArrayList<Matrix>();
                for (int k = 0; k < n; k++) {
                    Matrix[] pair = segmentPair(slot, pt == ProcessType.MAPT, k, n);
                    d0.add(pair[0]);
                    d1.add(pair[1]);
                }
                schedules.add(new Schedule(i, r, bp, cyclic, d0, d1));
            }
        }
    }

    /** (D0, D1) in force at t; the nominal when the process carries no schedule. */
    private Matrix[] pairAt(int station, int jobclass, double t) {
        for (int e = 0; e < schedules.size(); e++) {
            Schedule sc = schedules.get(e);
            if (sc.station == station && sc.jobclass == jobclass) {
                int idx = sc.segmentAt(t);
                if (idx < 0) {
                    int h = sc.d0.get(0).getNumRows();
                    return new Matrix[]{new Matrix(h, h), new Matrix(h, h)};
                }
                return new Matrix[]{sc.d0.get(idx), sc.d1.get(idx)};
            }
        }
        return new Matrix[]{nomD0[station][jobclass], nomD1[station][jobclass]};
    }

    private void buildEvents(int K) {
        events = new ArrayList<Event>();
        for (int b = 0; b < ublocks.size(); b++) {
            Block blk = ublocks.get(b);
            for (int k = 0; k < blk.phases; k++) {
                for (int j = 0; j < blk.phases; j++) {
                    if (k != j) {
                        double[] jump = new double[dim];
                        jump[blk.offset + k] -= 1.0;
                        jump[blk.offset + j] += 1.0;
                        events.add(new Event(A0, blk.station, blk.jobclass, k, j, -1, -1, -1, jump));
                    }
                }
            }
        }
        for (int b = 0; b < ublocks.size(); b++) {
            Block blk = ublocks.get(b);
            for (int d = 0; d < xblocks.size(); d++) {
                Block dst = xblocks.get(d);
                if (sn.rt.get(blk.station * K + blk.jobclass, dst.station * K + dst.jobclass) <= 0) {
                    continue;
                }
                for (int k = 0; k < blk.phases; k++) {
                    for (int j = 0; j < blk.phases; j++) {
                        for (int ip = 0; ip < dst.phases; ip++) {
                            double[] jump = new double[dim];
                            jump[blk.offset + k] -= 1.0;
                            jump[blk.offset + j] += 1.0;
                            jump[dst.offset + ip] += 1.0;
                            events.add(new Event(A1, blk.station, blk.jobclass, k, j,
                                    dst.station, dst.jobclass, ip, jump));
                        }
                    }
                }
            }
        }
        for (int b = 0; b < xblocks.size(); b++) {
            Block blk = xblocks.get(b);
            for (int p = 0; p < blk.phases; p++) {
                for (int q = 0; q < blk.phases; q++) {
                    if (p != q) {
                        double[] jump = new double[dim];
                        jump[blk.offset + p] -= 1.0;
                        jump[blk.offset + q] += 1.0;
                        events.add(new Event(SVC, blk.station, blk.jobclass, p, q, -1, -1, -1, jump));
                    }
                }
            }
        }
        for (int b = 0; b < xblocks.size(); b++) {
            Block blk = xblocks.get(b);
            if (exitProbability(blk, K) > 0) {
                for (int p = 0; p < blk.phases; p++) {
                    double[] jump = new double[dim];
                    jump[blk.offset + p] -= 1.0;
                    events.add(new Event(DEP, blk.station, blk.jobclass, p, -1, -1, -1, -1, jump));
                }
            }
            for (int d = 0; d < xblocks.size(); d++) {
                Block dst = xblocks.get(d);
                if (sn.rt.get(blk.station * K + blk.jobclass, dst.station * K + dst.jobclass) <= 0) {
                    continue;
                }
                for (int p = 0; p < blk.phases; p++) {
                    for (int ip = 0; ip < dst.phases; ip++) {
                        double[] jump = new double[dim];
                        jump[blk.offset + p] -= 1.0;
                        jump[dst.offset + ip] += 1.0;
                        events.add(new Event(ROUTE, blk.station, blk.jobclass, p, -1,
                                dst.station, dst.jobclass, ip, jump));
                    }
                }
            }
        }
    }

    /**
     * Probability that a completion leaves the network: sn.rt is closed through the Source, so
     * any destination that is not a service block is an exit.
     */
    private double exitProbability(Block blk, int K) {
        double pout = 0.0;
        for (int j = 0; j < sn.nstations; j++) {
            for (int l = 0; l < K; l++) {
                boolean isBlock = false;
                for (int d = 0; d < xblocks.size(); d++) {
                    if (xblocks.get(d).station == j && xblocks.get(d).jobclass == l) {
                        isBlock = true;
                        break;
                    }
                }
                if (!isBlock) {
                    pout += sn.rt.get(blk.station * K + blk.jobclass, j * K + l);
                }
            }
        }
        return pout;
    }

    // ------------------------------------------------------------------ dynamics

    private double capacityScale(double[] z, int station) {
        if (sn.sched.get(sn.stations.get(station)) == SchedStrategy.INF
                || !Double.isFinite(sn.nservers.get(station, 0))) {
            return 1.0;
        }
        double ni = 0.0;
        for (int b = 0; b < xblocks.size(); b++) {
            Block blk = xblocks.get(b);
            if (blk.station == station) {
                for (int i = 0; i < blk.phases; i++) {
                    ni += Math.max(z[blk.offset + i], 0.0);
                }
            }
        }
        double c = sn.nservers.get(station, 0);
        return ni <= c ? 1.0 : c / ni;
    }

    private Block blockOf(List<Block> list, int station, int jobclass) {
        for (int b = 0; b < list.size(); b++) {
            if (list.get(b).station == station && list.get(b).jobclass == jobclass) {
                return list.get(b);
            }
        }
        return null;
    }

    private double[] rates(double t, double[] q, int K) {
        double[] f = new double[events.size()];
        for (int e = 0; e < events.size(); e++) {
            Event ev = events.get(e);
            Matrix[] pair = pairAt(ev.station, ev.jobclass, t);
            if (ev.kind == A0) {
                Block blk = blockOf(ublocks, ev.station, ev.jobclass);
                f[e] = pair[0].get(ev.fromPhase, ev.toPhase)
                        * Math.max(q[blk.offset + ev.fromPhase], 0.0);
            } else if (ev.kind == A1) {
                Block blk = blockOf(ublocks, ev.station, ev.jobclass);
                double p = sn.rt.get(ev.station * K + ev.jobclass,
                        ev.destStation * K + ev.destClass);
                double beta = nomPie[ev.destStation][ev.destClass].get(0, ev.destPhase);
                f[e] = pair[1].get(ev.fromPhase, ev.toPhase) * p * beta
                        * Math.max(q[blk.offset + ev.fromPhase], 0.0);
            } else if (ev.kind == SVC) {
                Block blk = blockOf(xblocks, ev.station, ev.jobclass);
                f[e] = pair[0].get(ev.fromPhase, ev.toPhase)
                        * Math.max(q[blk.offset + ev.fromPhase], 0.0)
                        * capacityScale(q, ev.station);
            } else {
                Block blk = blockOf(xblocks, ev.station, ev.jobclass);
                double rowsum = 0.0;
                for (int j = 0; j < pair[1].getNumCols(); j++) {
                    rowsum += pair[1].get(ev.fromPhase, j);
                }
                double weight;
                if (ev.kind == DEP) {
                    weight = exitProbability(blk, K);
                } else {
                    double p = sn.rt.get(ev.station * K + ev.jobclass,
                            ev.destStation * K + ev.destClass);
                    weight = p * nomPie[ev.destStation][ev.destClass].get(0, ev.destPhase);
                }
                f[e] = rowsum * weight * Math.max(q[blk.offset + ev.fromPhase], 0.0)
                        * capacityScale(q, ev.station);
            }
        }
        return f;
    }

    private double[] drift(double t, double[] q, int K) {
        double[] f = rates(t, q, K);
        double[] dq = new double[dim];
        for (int e = 0; e < events.size(); e++) {
            double[] jump = events.get(e).jump;
            for (int i = 0; i < dim; i++) {
                if (jump[i] != 0.0) {
                    dq[i] += jump[i] * f[e];
                }
            }
        }
        return dq;
    }

    private static double[] stationaryPhase(Matrix d0, Matrix d1) {
        int h = d0.getNumRows();
        double[][] a = new double[h + 1][h];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < h; j++) {
                a[j][i] = d0.get(i, j) + d1.get(i, j); // transpose of the generator
            }
            a[h][i] = 1.0;
        }
        // least squares on [Q'; 1'] x = [0; 1] by normal equations, h being small
        double[][] ata = new double[h][h];
        double[] atb = new double[h];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < h; j++) {
                double s = 0.0;
                for (int k = 0; k <= h; k++) {
                    s += a[k][i] * a[k][j];
                }
                ata[i][j] = s;
            }
            atb[i] = a[h][i];
        }
        double[] theta = solveDense(ata, atb);
        double total = 0.0;
        for (int i = 0; i < h; i++) {
            theta[i] = Math.max(theta[i], 0.0);
            total += theta[i];
        }
        for (int i = 0; i < h; i++) {
            theta[i] = total > 0 ? theta[i] / total : 1.0 / h;
        }
        return theta;
    }

    private static double[] solveDense(double[][] a, double[] b) {
        int n = b.length;
        double[][] m = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(a[i], 0, m[i], 0, n);
            m[i][n] = b[i];
        }
        for (int c = 0; c < n; c++) {
            int piv = c;
            for (int r = c + 1; r < n; r++) {
                if (Math.abs(m[r][c]) > Math.abs(m[piv][c])) {
                    piv = r;
                }
            }
            double[] tmp = m[c];
            m[c] = m[piv];
            m[piv] = tmp;
            if (Math.abs(m[c][c]) < 1e-14) {
                continue;
            }
            for (int r = 0; r < n; r++) {
                if (r == c) {
                    continue;
                }
                double factor = m[r][c] / m[c][c];
                for (int k = c; k <= n; k++) {
                    m[r][k] -= factor * m[c][k];
                }
            }
        }
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = Math.abs(m[i][i]) < 1e-14 ? 0.0 : m[i][n] / m[i][i];
        }
        return x;
    }

    private static double summarise(double[] series, List<Double> times, double windowLo,
                                    double tend) {
        if (Double.isNaN(windowLo)) {
            return series[series.length - 1];
        }
        double area = 0.0;
        double span = 0.0;
        for (int n = 1; n < series.length; n++) {
            double ta = times.get(n - 1);
            double tb = times.get(n);
            if (ta < windowLo || tb > tend) {
                continue;
            }
            area += 0.5 * (series[n - 1] + series[n]) * (tb - ta);
            span += tb - ta;
        }
        return span > 0 ? area / span : series[series.length - 1];
    }

    /** The coupled (mean, covariance) system of Theorems 3.1 and 3.3. */
    private final class AugmentedODE implements FirstOrderDifferentialEquations {

        @Override
        public int getDimension() {
            return dim + dim * dim;
        }

        @Override
        public void computeDerivatives(double t, double[] z, double[] dz) {
            int K = sn.nclasses;
            double[] q = new double[dim];
            System.arraycopy(z, 0, q, 0, dim);
            double[] f = rates(t, q, K);
            double[] dq = drift(t, q, K);
            System.arraycopy(dq, 0, dz, 0, dim);

            // Jacobian by central differences on the assembled rates, so every capacity term is
            // differentiated consistently with the drift actually integrated.
            double scale = 1.0;
            for (int i = 0; i < dim; i++) {
                scale = Math.max(scale, Math.abs(q[i]));
            }
            double hstep = 1e-6 * scale;
            double[][] jac = new double[dim][dim];
            for (int m = 0; m < dim; m++) {
                double[] qp = q.clone();
                double[] qm = q.clone();
                qp[m] += hstep;
                qm[m] -= hstep;
                double[] fp = drift(t, qp, K);
                double[] fm = drift(t, qm, K);
                for (int i = 0; i < dim; i++) {
                    jac[i][m] = (fp[i] - fm[i]) / (2.0 * hstep);
                }
            }

            double[][] g = new double[dim][dim];
            for (int e = 0; e < events.size(); e++) {
                double[] jump = events.get(e).jump;
                double fe = f[e];
                if (fe == 0.0) {
                    continue;
                }
                for (int i = 0; i < dim; i++) {
                    if (jump[i] == 0.0) {
                        continue;
                    }
                    for (int j = 0; j < dim; j++) {
                        if (jump[j] != 0.0) {
                            g[i][j] += jump[i] * jump[j] * fe;
                        }
                    }
                }
            }

            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    double acc = g[i][j];
                    for (int k = 0; k < dim; k++) {
                        acc += jac[i][k] * z[dim + k * dim + j];
                        acc += z[dim + i * dim + k] * jac[j][k];
                    }
                    dz[dim + i * dim + j] = acc;
                }
            }
        }
    }
}
