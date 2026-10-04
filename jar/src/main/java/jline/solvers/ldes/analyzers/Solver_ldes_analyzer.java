/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ldes.analyzers;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.lang.processes.Distribution;
import jline.lang.processes.Markovian;
import jline.solvers.SolverOptions;
import jline.solvers.ldes.LDESResult;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.ldes.handlers.Solver_ssj;
import jline.solvers.ldes.handlers.Solver_ssj;
import jline.streaming.Collector;

public final class Solver_ldes_analyzer {
    private Solver_ldes_analyzer() {}

    /**
     * LDES analyzer that selects and executes the appropriate analysis method.
     */
    public static LDESResult solver_ldes_analyzer(NetworkStruct sn, SolverOptions options, SolverLDES solverLDES) {
        long Tstart = System.nanoTime();

        LDESResult res = new LDESResult();
        String requestedMethod = options.method;
        String backendMethod = requestedMethod;

        boolean isValidJackson = validateJacksonNetwork(sn);
        if (!isValidJackson) {
            throw new RuntimeException(
                    "solver_ldes_analyzer: Currently only Jackson queueing networks (possibly multiclass) are supported. "
                            + "The model must have: Source, one or more Queues with FCFS and exponential service, Sink, "
                            + "and probabilistic routing.");
        }

        boolean isTransient = options.timespan != null
                && options.timespan.length >= 2
                && Double.isFinite(options.timespan[1]);

        if ("default".equals(backendMethod)) {
            backendMethod = "ssj";
        }

        Collector stream = solverLDES.getStream();

        // see _kb/09-ldes-and-cache.md (Ensemble transient section (explicit initial state honored))
        if (options.init_sol == null || options.init_sol.isEmpty()) {
            jline.util.matrix.Matrix derived = deriveInitSolFromState(sn);
            if (derived != null) {
                options.init_sol = derived;
            }
        }

        if ("ssj".equals(backendMethod)) {
            if (isTransient) {
                res = Solver_ssj.solver_ssj_transient(sn, options, stream);
            } else {
                res = Solver_ssj.solver_ssj(sn, options, stream);
            }
        } else {
            throw new RuntimeException("solver_ldes_analyzer:UnknownMethod - Unknown analysis method: "
                    + options.method);
        }

        res.method = requestedMethod;
        res.runtime = (System.nanoTime() - Tstart) / 1.0e9;
        return res;
    }

    /**
     * Build a station-major init_sol vector [station0_class0, ..., stationM-1_classK-1] from
     * the model's initialized marginal state (sn.state). Returns null when no usable state is
     * present or when the per-class totals do not conserve the closed-class populations, in
     * which case the SSJ backend falls back to reference-station placement.
     */
    private static jline.util.matrix.Matrix deriveInitSolFromState(NetworkStruct sn) {
        if (sn.state == null || sn.state.isEmpty()) {
            return null;
        }
        int M = sn.nstations;
        int R = sn.nclasses;
        jline.util.matrix.Matrix initSol = new jline.util.matrix.Matrix(1, M * R);
        initSol.zero();
        double[] classTotals = new double[R];
        for (int i = 0; i < M; i++) {
            Station station = sn.stations.get(i);
            if (!(station instanceof jline.lang.nodes.StatefulNode)) {
                return null;
            }
            jline.util.matrix.Matrix nodeState = sn.state.get((jline.lang.nodes.StatefulNode) station);
            if (nodeState == null || nodeState.isEmpty()) {
                return null;
            }
            int nodeIdx = (int) sn.stationToNode.get(i);
            jline.lang.state.State.StateMarginalStatistics marg =
                    jline.lang.state.ToMarginal.toMarginal(sn, nodeIdx, nodeState, null, null, null, null, null);
            jline.util.matrix.Matrix nir = marg.nir;
            for (int r = 0; r < R; r++) {
                double v = (nir != null && r < nir.length()) ? nir.get(r) : 0.0;
                initSol.set(0, i * R + r, v);
                classTotals[r] += v;
            }
        }
        // Only use the derived placement when it exactly conserves every closed-class population.
        for (int r = 0; r < R; r++) {
            double njobs = sn.njobs.get(r);
            if (!Double.isInfinite(njobs) && Math.abs(classTotals[r] - njobs) > 1e-6) {
                return null;
            }
        }
        return initSol;
    }

    private static boolean validateJacksonNetwork(NetworkStruct sn) {
        boolean hasSource = false;
        boolean hasSink = false;
        boolean hasServiceNode = false;
        boolean hasPlace = false;
        boolean hasTransition = false;
        boolean hasCache = false;

        for (NodeType nodeType : sn.nodetype) {
            if (nodeType == NodeType.Source) hasSource = true;
            else if (nodeType == NodeType.Sink) hasSink = true;
            else if (nodeType == NodeType.Queue) hasServiceNode = true;
            else if (nodeType == NodeType.Delay) hasServiceNode = true;
            else if (nodeType == NodeType.Place) hasPlace = true;
            else if (nodeType == NodeType.Transition) hasTransition = true;
            else if (nodeType == NodeType.Cache) hasCache = true;
        }

        boolean hasOpenClasses = false;
        for (int k = 0; k < sn.nclasses; k++) {
            if (Double.isInfinite(sn.njobs.get(k))) {
                hasOpenClasses = true;
                break;
            }
        }

        boolean hasPetriNet = hasPlace || hasTransition;
        if (!hasServiceNode && !hasPetriNet && !hasCache) {
            return false;
        }

        if (hasOpenClasses && (!hasSource || !hasSink)) {
            return false;
        }

        Set<SchedStrategy> supportedStrategies = new HashSet<SchedStrategy>(Arrays.asList(
                SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRIO,
                SchedStrategy.FCFSPR, SchedStrategy.FCFSPI,
                SchedStrategy.FCFSPRPRIO, SchedStrategy.FCFSPIPRIO,
                SchedStrategy.LCFS, SchedStrategy.LCFSPR, SchedStrategy.LCFSPI,
                SchedStrategy.LCFSPRIO, SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO,
                SchedStrategy.EXT,
                SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
                SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO,
                SchedStrategy.LPS,
                SchedStrategy.SIRO,
                SchedStrategy.SJF, SchedStrategy.LJF,
                SchedStrategy.LEPT, SchedStrategy.SEPT,
                SchedStrategy.SRPT, SchedStrategy.SRPTPRIO,
                SchedStrategy.PSJF, SchedStrategy.FB, SchedStrategy.LRPT,
                SchedStrategy.FSP,
                SchedStrategy.POLLING,
                SchedStrategy.PAS, SchedStrategy.OI,
                SchedStrategy.EDD, SchedStrategy.EDF, SchedStrategy.SETF
        ));

        for (Map.Entry<Station, SchedStrategy> entry : sn.sched.entrySet()) {
            Station station = entry.getKey();
            SchedStrategy sched = entry.getValue();
            int stationIdx = station.stationIdx;
            NodeType nodeType = sn.nodetype.get((int) sn.stationToNode.get(stationIdx));
            if (nodeType == NodeType.Queue) {
                if (!supportedStrategies.contains(sched)) {
                    return false;
                }
            } else if (nodeType == NodeType.Delay) {
                if (sched != SchedStrategy.INF && sched != SchedStrategy.EXT) {
                    return false;
                }
            }
        }

        Set<ProcessType> allowedProcTypes = new HashSet<ProcessType>(Arrays.asList(
                ProcessType.DISABLED,
                ProcessType.EXP,
                ProcessType.PH,
                ProcessType.APH,
                ProcessType.HYPEREXP,
                ProcessType.COXIAN,
                ProcessType.COX2,
                ProcessType.ERLANG,
                ProcessType.MAP,
                ProcessType.DMAP,
                ProcessType.GEOMETRIC,
                ProcessType.BERNOULLI,
                ProcessType.BINOMIAL,
                ProcessType.POISSON,
                ProcessType.MMPP2,
                ProcessType.BMAP,
                ProcessType.MMAP,
                ProcessType.ME,
                ProcessType.RAP,
                ProcessType.IMMEDIATE,
                ProcessType.REPLAYER,
                ProcessType.DET,
                ProcessType.UNIFORM,
                ProcessType.GAMMA,
                ProcessType.PARETO,
                ProcessType.WEIBULL,
                ProcessType.LOGNORMAL,
                ProcessType.NHPP
        ));
        for (Map.Entry<Station, Map<JobClass, ProcessType>> entry : sn.procid.entrySet()) {
            for (Map.Entry<JobClass, ProcessType> inner : entry.getValue().entrySet()) {
                ProcessType procType = inner.getValue();
                if (!allowedProcTypes.contains(procType)) {
                    return false;
                }
            }
        }

        Set<String> supportedDists = new HashSet<String>(Arrays.asList(
                "EXP", "EXPONENTIAL", "ERLANG", "HYPEREXP", "PH", "APH", "COXIAN", "COX2"));

        Set<SchedStrategy> psVariants = new HashSet<SchedStrategy>(Arrays.asList(
                SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
                SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO));

        for (Station station : sn.stations) {
            if (station instanceof Queue && ((Queue) station).isDelayOffEnabled()) {
                SchedStrategy sched = sn.sched.get(station);

                if (psVariants.contains(sched)) {
                    throw new RuntimeException(
                            "solver_ldes_analyzer: Setup/delayoff is incompatible with Processor Sharing (PS) scheduling. "
                                    + "Station '" + station.getName() + "' uses " + sched
                                    + " which does not support server on/off transitions.");
                }

                Queue q = (Queue) station;
                for (JobClass jobClass : sn.jobclasses) {
                    Distribution setupDist = q.getSetupTime(jobClass);
                    Distribution delayoffDist = q.getDelayOffTime(jobClass);

                    if (setupDist != null && !setupDist.isDisabled()) {
                        String distName = setupDist.getName().toUpperCase();
                        if (!supportedDists.contains(distName) && !(setupDist instanceof Markovian)) {
                            throw new RuntimeException(
                                    "solver_ldes_analyzer: Setup time distribution for station '" + station.getName()
                                            + "' class '" + jobClass.getName()
                                            + "' must be Exponential, Erlang, HyperExp, PH, APH, or Coxian. Distribution '"
                                            + setupDist.getName() + "' is not supported.");
                        }
                    }

                    if (delayoffDist != null && !delayoffDist.isDisabled()) {
                        String distName = delayoffDist.getName().toUpperCase();
                        if (!supportedDists.contains(distName) && !(delayoffDist instanceof Markovian)) {
                            throw new RuntimeException(
                                    "solver_ldes_analyzer: Delayoff time distribution for station '" + station.getName()
                                            + "' class '" + jobClass.getName()
                                            + "' must be Exponential, Erlang, HyperExp, PH, APH, or Coxian. Distribution '"
                                            + delayoffDist.getName() + "' is not supported.");
                        }
                    }
                }
            }
        }

        return true;
    }
}
