/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.moments;

import jline.GlobalConstants;
import jline.api.sn.SnRtStations;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.FluidNhpp;
import jline.solvers.fluid.handlers.PassageTimeODE;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Event-based representation of the fluid population process, as required by the
 * moment-closure methods of SolverFluid. Java twin of the MATLAB
 * {@code fluid_moment_terms}.
 *
 * <p>The closing ODEs are a density-dependent Markov population process</p>
 *
 * <pre>
 *   dx/dt = F(x) = D * r(x),     r_e(x) = rateBase(e) * g_e(x)
 * </pre>
 *
 * <p>with D the jump matrix and g the rate-factor vector of
 * {@link FluidRateFactors}. The closing analyzer discards D and r once it has
 * composed the right-hand side, but the covariance equation of the linear noise
 * approximation needs them separately: the diffusion matrix is
 * {@code D*diag(r(x))*D'}, which cannot be recovered from F alone. This class
 * rebuilds that representation and exposes drift, rate and Jacobian evaluations
 * that all take the closure variance as an explicit argument.</p>
 *
 * @see FluidLyapunov
 * @see jline.solvers.fluid.analyzers.MinNormalAnalyzer
 */
public class FluidMomentTerms {

    public final int M;
    public final int K;
    public final int nstate;

    /** Jump matrix (nstate x nevents). */
    public final Matrix D;
    /** Constant rate factor of each event (nevents x 1). */
    public final Matrix rateBase;
    /** Source state coordinate of each event. */
    public final int[] eventIdx;

    public final Matrix qIndices;
    public final Matrix Kic;
    public final boolean[][] enabled;
    public final SchedStrategy[] sched;
    /** Effective server counts, with an infinite server replaced by the closed population. */
    public final Matrix S;
    public final Matrix lldscaling;

    /** State coordinates of each station, and of each (station, class). */
    public final int[][] stationBlock;
    public final int[][] classBlock;

    /** Coordinates carrying a real population, i.e. everything but the EXT source pools. */
    public final int[] covIdx;
    public final boolean[] isExt;
    /**
     * Stations whose occupancy cannot reach their server count, where min(n,c) is the
     * identity and the Gaussian closure must stay first order. See buildMinExact.
     */
    public final boolean[] minExact;

    /** Event classification, so throughputs can be read off the rate vector. */
    public final boolean[] evIsDeparture;
    public final int[] evStation;
    public final int[] evClass;
    /**
     * [nEventsReduced x nEventsOriginal] expected firings of each original event per firing of each
     * reduced one; the identity when no immediate coordinate was eliminated. Event attributes above
     * are indexed by ORIGINAL event, so a throughput is read as
     * {@code r' * (Emap * indicatorOverOriginalEvents)}.
     */
    public final Matrix Emap;
    /**
     * Projector taking an initial condition onto the surviving coordinates: the
     * identity on a coordinate the immediate reduction kept and the absorption
     * distribution on one it folded away, so a zero diagonal marks an eliminated
     * coordinate. Null when nothing was eliminated.
     */
    public final Matrix immediateAbsorb;

    private final FluidRateFactors factors;

    /**
     * @param sn      network structure, after the non-Markovian to phase-type conversion
     * @param options solver options
     */
    public FluidMomentTerms(NetworkStruct sn, SolverOptions options) {
        this.M = sn.nstations;
        this.K = sn.nclasses;
        int N = sn.nclosedjobs;

        // the moment methods need an autonomous drift
        if (options.config != null
                && (options.config.rate_traj_mmat != null || options.config.nhpp_sched != null
                    || options.config.rate_sched != null)) {
            line_error(mfilename(new Object() {
            }), "Moment-closure methods require an autonomous drift, but options.config.rate_traj/nhpp_sched "
                    + "make the rates time-varying. Use options.method=\"closing\" or \"matrix\".");
        }

        Map<Station, Map<JobClass, Matrix>> mu = new HashMap<Station, Map<JobClass, Matrix>>();
        Map<Station, Map<JobClass, Matrix>> phi = new HashMap<Station, Map<JobClass, Matrix>>();
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            Map<JobClass, Matrix> muCopy = new HashMap<JobClass, Matrix>();
            Map<JobClass, Matrix> phiCopy = new HashMap<JobClass, Matrix>();
            for (int k = 0; k < K; k++) {
                JobClass jobClass = sn.jobclasses.get(k);
                Matrix muik = sn.mu.get(station).get(jobClass);
                Matrix phiik = sn.phi.get(station).get(jobClass);
                if (muik == null || muik.hasNaN()) {
                    muCopy.put(jobClass, new Matrix(0, 0));
                    phiCopy.put(jobClass, new Matrix(0, 0));
                } else {
                    muCopy.put(jobClass, muik.copy());
                    phiCopy.put(jobClass, phiik.copy());
                }
            }
            mu.put(station, muCopy);
            phi.put(station, phiCopy);
        }

        this.S = sn.nservers.copy();
        for (int i = 0; i < M; i++) {
            if (Double.isInfinite(S.get(i, 0))) {
                // a pure open model has no closed population, so N = 0 would make the
                // utilization divisor vanish
                S.set(i, 0, Math.max(N, 1));
            }
        }

        this.sched = new SchedStrategy[M];
        this.isExt = new boolean[M];
        for (int i = 0; i < M; i++) {
            sched[i] = sn.sched.get(sn.stations.get(i));
            isExt[i] = sched[i] == SchedStrategy.EXT;
        }
        this.lldscaling = (sn.lldscaling == null || sn.lldscaling.isEmpty()) ? null : sn.lldscaling;

        Map<Station, Map<JobClass, MatrixCell>> proc = FluidNhpp.substituteNhppProc(sn, mu, phi);
        SolverOptions odeOptions = options.copy();
        odeOptions.method = "closing";
        int dim = 0;
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            for (int k = 0; k < K; k++) {
                Matrix muik = mu.get(station).get(sn.jobclasses.get(k));
                dim += (muik == null || muik.isEmpty()) ? 0 : muik.length();
            }
        }
        // station-major routing; sn.rt is indexed by stateful node (see SnRtStations)
        Matrix rtStations = SnRtStations.snRtStations(sn).getLeft();
        PassageTimeODE ode = new PassageTimeODE(sn, mu, phi, proc, rtStations, S, odeOptions, dim);

        this.D = ode.getAllJumps();
        this.rateBase = ode.getRateBase();
        this.qIndices = ode.getQIndices();
        this.Kic = ode.getKic();
        this.enabled = ode.getEnabled();
        this.factors = ode.getRateFactors();
        this.nstate = dim;

        Matrix eventIdxMat = ode.getEventIdx();
        int nevents = eventIdxMat.getNumRows();
        this.eventIdx = new int[nevents];
        for (int e = 0; e < nevents; e++) {
            eventIdx[e] = (int) eventIdxMat.get(e, 0);
        }

        // THE MOMENT CLOSURE READS THE SAME REDUCED EVENT SET AS EVERY OTHER ROUTE. It used to
        // refuse the reduction, on the grounds that it needs the untransformed event set; what it
        // actually needs is to be able to say which (station,class) each event is a completion of,
        // and Emap carries exactly that across the composition -- an event folded through an
        // immediate coordinate keeps a row with weight on every original event it stands for,
        // including the two completions a pass-through realises at once. The diffusion D*diag(r)*D'
        // is then the diffusion of the reduced process, which is the right one: the eliminated
        // coordinate holds O(1/InfRate) mass and contributes noise of the same order.
        Matrix emap = ode.getImmediateEmap();
        Matrix origEventIdxMat = ode.getOriginalEventIdx();
        if (origEventIdxMat == null) {
            origEventIdxMat = eventIdxMat;
        }
        int nevents0 = origEventIdxMat.getNumRows();
        this.Emap = (emap != null) ? emap : Matrix.eye(nevents);
        this.immediateAbsorb = ode.getImmediateAbsorb();
        int[] eventIdx0 = new int[nevents0];
        for (int e = 0; e < nevents0; e++) {
            eventIdx0[e] = (int) origEventIdxMat.get(e, 0);
        }

        this.stationBlock = new int[M][];
        this.classBlock = new int[M * K][];
        for (int i = 0; i < M; i++) {
            int lo = (int) qIndices.get(i, 0);
            int hi = (int) qIndices.get(i, K - 1) + (int) Kic.get(i, K - 1);
            stationBlock[i] = range(lo, hi);
            for (int k = 0; k < K; k++) {
                int klo = (int) qIndices.get(i, k);
                int khi = klo + (int) Kic.get(i, k);
                classBlock[i * K + k] = range(klo, khi);
            }
        }

        this.covIdx = buildCovIdx();
        this.minExact = buildMinExact(sn);
        // Indexed by ORIGINAL event, which is what Emap maps onto. Without a reduction Emap is the
        // identity and the two indexings coincide, exactly as before.
        this.evIsDeparture = buildDepartureFlags(rtStations, nevents0);

        int[] coordStation = new int[nstate];
        int[] coordClass = new int[nstate];
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                int[] blk = classBlock[i * K + k];
                for (int j = 0; j < blk.length; j++) {
                    coordStation[blk[j]] = i;
                    coordClass[blk[j]] = k;
                }
            }
        }
        this.evStation = new int[nevents0];
        this.evClass = new int[nevents0];
        for (int e = 0; e < nevents0; e++) {
            evStation[e] = coordStation[eventIdx0[e]];
            evClass[e] = coordClass[eventIdx0[e]];
        }
    }

    /**
     * Coordinates the covariance is solved on: everything but the EXT source
     * pools.
     *
     * <p>The closing representation models a Source as an EXT pseudo-station
     * holding unit mass, so its coordinate is a normalisation constant, not a job
     * count. Building {@code D*diag(r)*D'} over it would invent noise for a
     * direction that has no population. Projecting those coordinates away leaves
     * exactly the right open event set, because the closing form already emits
     * the correct events: with a single-phase source the EXT rate factor is
     * {@code 1 - sum(of nothing) = 1} identically, so an arrival is a
     * CONSTANT-rate event whose jump, once the source row is dropped, is a lone
     * +1 into the destination queue, the canonical exogenous Poisson arrival with
     * diffusion intensity lambda. The return leg that LINE routes Sink to Source
     * becomes a lone -1, the departure. The EXT row of the Jacobian is identically
     * zero for a single-phase source, so A(covIdx,covIdx) is exactly the Jacobian
     * of the projected drift, not an approximation of it.</p>
     *
     * <p>A MULTI-PHASE source is refused: those coordinates are the phase of ONE
     * arrival process, a single Markov chain rather than a population, so its
     * fluctuations are O(1) and the linear noise approximation does not apply to
     * them at any scale.</p>
     */
    private int[] buildCovIdx() {
        boolean[] covMask = new boolean[nstate];
        for (int j = 0; j < nstate; j++) {
            covMask[j] = true;
        }
        for (int i = 0; i < M; i++) {
            if (!isExt[i]) {
                continue;
            }
            for (int k = 0; k < K; k++) {
                int kic = (int) Kic.get(i, k);
                if (kic <= 0) {
                    continue;
                }
                if (kic > 1) {
                    line_error(mfilename(new Object() {
                    }), String.format("The moment-closure methods need a Poisson arrival stream, but the source of "
                            + "class %d is a %d-phase process. Those coordinates track the phase of a single arrival "
                            + "process rather than a population, so they carry no linear noise approximation. Use an "
                            + "exponential inter-arrival time, or options.method=\"matrix\".", k + 1, kic));
                }
                int klo = (int) qIndices.get(i, k);
                for (int j = klo; j < klo + kic; j++) {
                    covMask[j] = false;
                }
            }
        }
        List<Integer> keep = new ArrayList<Integer>();
        for (int j = 0; j < nstate; j++) {
            if (covMask[j]) {
                keep.add(j);
            }
        }
        int[] idx = new int[keep.size()];
        for (int j = 0; j < idx.length; j++) {
            idx[j] = keep.get(j);
        }
        return idx;
    }

    /**
     * A STATION THAT CANNOT FILL ITS SERVERS HAS NOTHING TO CLOSE. min(n,c) is the
     * identity on the whole support whenever the occupancy of a station is bounded
     * above by its server count, and there the Gaussian closure is not an improvement
     * on the first-order one, it is an ERROR: it spreads a normal marginal over
     * n &gt; c, mass the station can never hold, and returns E[min(n,c)] &lt; n. On a
     * closed model with one job per chain the exact answer is R = D at every queue (a
     * job cannot queue behind itself), which the first-order closure reproduces to
     * machine precision while the closure reads 0.4758 against 0.5 on the queue length
     * and 1.1017 against 1 on the response time. The bound is the total population of
     * every chain that VISITS the station -- a station may declare a service time for
     * every class while the routing never sends most of them there -- and an open
     * chain contributes an infinite population and never qualifies. MinNormalAnalyzer
     * holds the drift variance of these stations at zero, exactly as it does for the
     * delay stations, whose min() is likewise absent.
     */
    private boolean[] buildMinExact(NetworkStruct sn) {
        boolean[] exact = new boolean[M];
        if (sn.chains == null || sn.chains.isEmpty() || sn.njobs == null) {
            return exact;
        }
        int nchains = sn.chains.getNumRows();
        for (int i = 0; i < M; i++) {
            if (isExt[i] || sched[i] == SchedStrategy.INF
                    || Double.isInfinite(sn.nservers.get(i, 0))) {
                continue;
            }
            double bound = 0;
            boolean[] covered = new boolean[K];
            for (int ch = 0; ch < nchains; ch++) {
                boolean here = false;
                double pop = 0;
                Matrix vis = (sn.visits == null) ? null : sn.visits.get(ch);
                for (int k = 0; k < K; k++) {
                    if (sn.chains.get(ch, k) == 0) {
                        continue;
                    }
                    covered[k] = true;
                    pop += sn.njobs.get(k);
                    // sn.visits is indexed by STATEFUL node, not by station
                    // stationToStateful is a 1 x nstations ROW vector, so read it linearly
                    if (vis != null && !vis.isEmpty() && sn.stationToStateful != null) {
                        here = here || vis.get((int) sn.stationToStateful.get(i), k) > 0;
                    } else {
                        here = here || enabled[i][k];
                    }
                }
                if (here) {
                    bound += pop;
                }
            }
            boolean uncovered = false;
            for (int k = 0; k < K; k++) {
                if (enabled[i][k] && !covered[k]) {
                    uncovered = true;
                }
            }
            if (uncovered) {
                continue; // a class outside every chain carries no population bound
            }
            exact[i] = !Double.isInfinite(bound)
                    && bound <= sn.nservers.get(i, 0) + GlobalConstants.FineTol;
        }
        return exact;
    }

    /**
     * Which events are service completions. The rate-base builder emits every
     * service completion first, then every intra-PH phase change, so summing the
     * completion rates sourced at (i,c) gives the class-c throughput at station i
     * exactly: the routing probabilities and the PH entry vector each sum to one
     * over the destinations enumerated there.
     */
    private boolean[] buildDepartureFlags(Matrix rtStations, int nevents) {
        int nDeparture = 0;
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < K; c++) {
                if (!enabled[i][c]) {
                    continue;
                }
                for (int j = 0; j < M; j++) {
                    for (int l = 0; l < K; l++) {
                        if (rtStations.get(i * K + c, j * K + l) > 0) {
                            nDeparture += (int) Kic.get(i, c) * (int) Kic.get(j, l);
                        }
                    }
                }
            }
        }
        boolean[] flags = new boolean[nevents];
        for (int e = 0; e < Math.min(nDeparture, nevents); e++) {
            flags[e] = true;
        }
        return flags;
    }

    /** Per-coordinate service shares at the given closure variance. */
    public Matrix factors(double[] x, double[] sigma2, Matrix[] covblk) {
        return factors.factors(x, sigma2, covblk);
    }

    /** Event rates, i.e. the rate base scaled by the service shares. */
    public Matrix rates(double[] x, double[] sigma2, Matrix[] covblk) {
        Matrix g = factors.factors(x, sigma2, covblk);
        Matrix r = new Matrix(eventIdx.length, 1);
        for (int e = 0; e < eventIdx.length; e++) {
            r.set(e, 0, rateBase.get(e, 0) * g.get(eventIdx[e], 0));
        }
        return r;
    }

    /** Fluid drift F(x) = D*r(x). */
    public Matrix drift(double[] x, double[] sigma2, Matrix[] covblk) {
        return D.mult(rates(x, sigma2, covblk), null);
    }

    /**
     * First station sitting on the saturation kink of the first-order rate
     * factor, or -1 when none does. See
     * {@link FluidRateFactors#driftKinkStation(double[], double[])}: at such a
     * point {@link #jacobian} exists only one-sidedly, so a caller needing a
     * differentiable drift must check this rather than trust the returned slope.
     */
    public int driftKinkStation(double[] x, double[] sigma2) {
        return factors.driftKinkStation(x, sigma2);
    }

    /**
     * The state with every kinked station moved strictly onto one side of its
     * kink. See {@link FluidRateFactors#nudgedOffKink(double[], double[], double)}.
     */
    public double[] nudgedOffKink(double[] x, double[] sigma2, double rel) {
        return factors.nudgedOffKink(x, sigma2, rel);
    }

    /** Analytic Jacobian dF/dx = D*(rateBase .* G(eventIdx,:)). */
    public Matrix jacobian(double[] x, double[] sigma2, Matrix[] covblk) {
        Matrix G = factors.jacobian(x, sigma2, covblk);
        Matrix RG = new Matrix(eventIdx.length, nstate);
        for (int e = 0; e < eventIdx.length; e++) {
            double rb = rateBase.get(e, 0);
            if (rb == 0.0) {
                continue;
            }
            for (int j = 0; j < nstate; j++) {
                RG.set(e, j, rb * G.get(eventIdx[e], j));
            }
        }
        return D.mult(RG, null);
    }

    /** Refuses any station whose discipline has no drift branch. */
    public void checkSupported() {
        factors.checkSupported();
    }

    private static int[] range(int lo, int hi) {
        int[] r = new int[Math.max(0, hi - lo)];
        for (int j = 0; j < r.length; j++) {
            r[j] = lo + j;
        }
        return r;
    }
}
