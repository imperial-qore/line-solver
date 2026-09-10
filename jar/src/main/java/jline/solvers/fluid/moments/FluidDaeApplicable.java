/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.moments;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Whether {@link jline.solvers.fluid.analyzers.DaeAnalyzer} can answer a model,
 * used by SolverFluid to try {@code dae} before dropping a declined
 * {@code minnormal} to a first-order method. Java twin of the MATLAB
 * {@code fluid_dae_applicable}.
 *
 * <p>WHY THIS EXISTS SEPARATELY FROM {@link FluidMinNormalApplicable}. The two
 * methods state the SAME closure and differ only in how the coupled equations
 * are discharged, so a model {@code minnormal} accepts is almost always one
 * {@code dae} accepts too. Almost: {@code dae} carries a finite-difference
 * Jacobian over the whole unknown vector rather than one Lyapunov solve, so its
 * state cap is lower; it closes on the per-station variance only, so DPS and GPS
 * are out; and it has no decomposition route, so a cache model is out. Those
 * three are exactly the difference set, and naming them here keeps the fallback
 * ladder from entering a rung that would refuse the model a moment later.</p>
 *
 * <p>The test is STATIC, for the same reason the min-normal one is: a rung
 * chosen by trial and rollback would make the reported method depend on a failed
 * run. The one condition that cannot be static is the NON-HYPERBOLIC fixed point
 * the ladder exists to route around -- it exists only once the mean is solved --
 * and {@code dae} fails on it loudly, which is what moves the ladder to its last
 * rung.</p>
 *
 * <p>Every condition below mirrors a refusal {@code DaeAnalyzer} would otherwise
 * raise, so this class and those refusals must move together.</p>
 */
public final class FluidDaeApplicable {

    private FluidDaeApplicable() {
    }

    /**
     * @param sn      network structure, after the non-Markovian to phase-type conversion
     * @param options solver options
     * @return null when {@code dae} can be selected, otherwise one line naming
     *         the blocking feature
     */
    public static String reasonToDecline(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;
        int K = sn.nclasses;

        // A cache model is a DECOMPOSITION, and dae has no arm for it: minnormal
        // on a cache model routes through the rmf alternation with the closure
        // inside its network step, while DaeAnalyzer would read the cache nodes
        // as ordinary stations, because a decomposition has no single drift for
        // the algebraic constraint to attach to.
        if (sn.nodetype != null) {
            for (int i = 0; i < sn.nodetype.size(); i++) {
                if (sn.nodetype.get(i) == NodeType.Cache) {
                    return "a cache model is answered by the decomposition analyzer, which has no dae route";
                }
            }
        }

        // The DPS and GPS shares close on the covariance BETWEEN station
        // coordinates, not on the station total, so their closure state is a
        // matrix block rather than the scalar the Newton vector carries. Mirrors
        // the DaeAnalyzer refusal.
        for (int i = 0; i < M; i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            if (s == SchedStrategy.DPS || s == SchedStrategy.GPS) {
                return String.format("station %d uses %s, whose share closes on the covariance between its "
                        + "class coordinates rather than on the station variance", i + 1, s);
            }
        }

        // The simultaneous solve is quartic overall, against the cubic of one
        // Lyapunov solve, so its crossover is lower than the 200 minnormal
        // permits and it carries its own cap. Same count FluidMinNormalApplicable
        // forms, so the two limits are read on the same scale.
        int maxstate = (options.config == null) ? 0 : options.config.dae_maxstate;
        if (maxstate <= 0) {
            maxstate = 100;
        }
        int nstate = 0;
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            for (int r = 0; r < K; r++) {
                nstate += phaseCount(sn, station, sn.jobclasses.get(r), i, r);
            }
        }
        if (nstate > maxstate) {
            return String.format("the phase-resolved state has %d coordinates, above the dae_maxstate limit of %d",
                    nstate, maxstate);
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
