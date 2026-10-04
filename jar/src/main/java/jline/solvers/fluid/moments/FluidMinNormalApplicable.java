/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.moments;

import jline.lang.JobClass;
import jline.lang.NodeParam;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Whether {@link jline.solvers.fluid.analyzers.MinNormalAnalyzer} can answer a
 * model, used by the {@code default} method of SolverFluid to prefer
 * {@code minnormal} over {@code matrix} when it applies. Java twin of the MATLAB
 * {@code fluid_minnormal_applicable}.
 *
 * <p>The test is STATIC: it inspects the model and the options, never the
 * solution. A method chosen by trial and rollback would make the reported method
 * depend on a failed run, and the fallback would silently absorb real defects in
 * the closure; anything the check cannot decide in advance is left to fail
 * loudly under an explicit {@code options.method="minnormal"}.</p>
 *
 * <p>Every condition below mirrors a guard that {@code MinNormalAnalyzer} or
 * {@link FluidMomentTerms} would otherwise raise, so this class and those guards
 * must move together.</p>
 */
public final class FluidMinNormalApplicable {

    private FluidMinNormalApplicable() {
    }

    /**
     * @param sn      network structure, after the non-Markovian to phase-type conversion
     * @param options solver options
     * @return null when {@code minnormal} can be selected, otherwise one line
     *         naming the blocking feature
     */
    public static String reasonToDecline(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;

        // Open and mixed models are supported: FluidMomentTerms projects the EXT source
        // pool out of the covariance. What it cannot take is a NON-POISSON arrival
        // stream, whose source coordinates track the phase of a single arrival process
        // rather than a population.
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            if (sn.sched.get(station) != SchedStrategy.EXT) {
                continue;
            }
            for (int r = 0; r < K; r++) {
                int phases = phaseCount(sn, station, sn.jobclasses.get(r), i, r);
                if (phases > 1) {
                    return String.format("class %d has a %d-phase (non-Poisson) arrival process", r + 1, phases);
                }
            }
        }
        // A cache model is answered through the decomposition analyzer, with the
        // closure in its network step, so cache nodes no longer decline the
        // method. What still declines is a replacement strategy with no
        // drift-based fluid model, mirroring the runtime guard in RMFAnalyzer:
        // LRU, HLRU, CLIMB and QLRU are answered by a characteristic-time fixed
        // point, which is not a fluid method and has no covariance.
        if (sn.nodetype != null) {
            for (int i = 0; i < sn.nodetype.size(); i++) {
                if (sn.nodetype.get(i) != NodeType.Cache) {
                    continue;
                }
                NodeParam param = sn.nodeparam == null ? null : sn.nodeparam.get(sn.nodes.get(i));
                if (!(param instanceof CacheNodeParam)) {
                    return "a cache uses a replacement strategy with no drift-based fluid model";
                }
                ReplacementStrategy strat = ((CacheNodeParam) param).replacestrat;
                if (strat != ReplacementStrategy.RR && strat != ReplacementStrategy.FIFO
                        && strat != ReplacementStrategy.SFIFO) {
                    return "a cache uses a replacement strategy with no drift-based fluid model";
                }
            }
        }

        for (int i = 0; i < M; i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            if (!FluidRateFactors.isSupported(s)) {
                return String.format("station %d uses %s, which has no fluid drift branch", i + 1, s);
            }
        }

        // the Lyapunov solve is cubic in the phase-resolved state, so the same cap the
        // analyzer enforces decides selection rather than being hit later
        int maxstate = (options.config == null) ? 200 : options.config.moment_maxstate;
        if (maxstate <= 0) {
            maxstate = 200;
        }
        int nstate = 0;
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            for (int r = 0; r < K; r++) {
                nstate += phaseCount(sn, station, sn.jobclasses.get(r), i, r);
            }
        }
        if (nstate > maxstate) {
            return String.format("the phase-resolved state has %d coordinates, above the moment_maxstate limit of %d",
                    nstate, maxstate);
        }

        // the moment methods need the untransformed event set and an autonomous drift
        if (options.config != null && options.config.rate_traj_mmat != null) {
            return "options.config.rate_traj makes the drift time-varying";
        }
        if (options.config != null && (options.config.nhpp_sched != null || options.config.rate_sched != null)) {
            return "options.config.nhpp_sched makes the drift time-varying";
        }
        return null;
    }

    private static int phaseCount(NetworkStruct sn, Station station, JobClass jobClass, int i, int r) {
        if (sn.rates != null && Double.isNaN(sn.rates.get(i, r))) {
            return 0;
        }
        Matrix muik = sn.mu.get(station).get(jobClass);
        if (muik == null || muik.isEmpty() || muik.hasNaN()) {
            return 0;
        }
        return muik.length();
    }
}
