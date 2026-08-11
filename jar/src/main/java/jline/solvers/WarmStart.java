/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.util.matrix.Matrix;

/**
 * Warm-start placement computation shared by the solvers that can start from
 * an auxiliary solver's steady-state solution (LDES, SSA, JMT, Fluid).
 *
 * The placement is an integer (nstations x nclasses) job assignment decided by
 * the steady-state distribution of the auxiliary solver: with SolverCTMC it is
 * the mode of the exact stationary distribution over the aggregate state
 * space, with any other solver the rounded steady-state mean queue lengths,
 * conserving each closed-class population. Only service stations (Queue/Delay)
 * receive an initial population.
 */
public final class WarmStart {
    private WarmStart() {}

    /**
     * Integer job placement (nstations x nclasses) decided by the steady-state
     * solution of the auxiliary solver.
     *
     * @param initSolver auxiliary solver used to compute the steady-state distribution
     * @param sn struct of the model to be warm-started
     * @return placement matrix (nstations x nclasses)
     */
    public static Matrix warmStartPlacement(NetworkSolver initSolver, NetworkStruct sn) {
        if (initSolver instanceof jline.solvers.ctmc.SolverCTMC) {
            return placementFromCtmcSteadyState((jline.solvers.ctmc.SolverCTMC) initSolver, sn);
        }
        return placementFromMeanQLen(initSolver, sn);
    }

    /**
     * Placement from the exact CTMC stationary distribution: solve the CTMC,
     * aggregate the stationary probabilities over the aggregate (per-station,
     * per-class job count) state space, and return the aggregate state of
     * maximum stationary probability.
     */
    private static Matrix placementFromCtmcSteadyState(jline.solvers.ctmc.SolverCTMC ctmcSolver, NetworkStruct sn) {
        int M = sn.nstations;
        int K = sn.nclasses;
        NetworkStruct snCtmc = ctmcSolver.getModel().getStruct(true);
        jline.solvers.ctmc.ResultCTMC ctmcRes =
                jline.solvers.ctmc.handlers.Solver_ctmc.solver_ctmc(snCtmc, ctmcSolver.getOptions());
        Matrix Q = ctmcRes.getQ();
        Matrix SSq = ctmcRes.getStateSpaceAggr();
        snCtmc = ctmcRes.getSn();

        jline.util.Pair<Matrix, List<List<Integer>>> sol =
                jline.api.mc.Ctmc_solve_reducible.ctmc_solve_reducible(Q);
        Matrix pi = sol.getLeft();

        // Aggregate the stationary probability over identical aggregate states and
        // locate the mode of the aggregate distribution.
        Map<String, Double> aggrProb = new LinkedHashMap<String, Double>();
        Map<String, Integer> aggrRep = new LinkedHashMap<String, Integer>();
        int nStates = SSq.getNumRows();
        int nCols = SSq.getNumCols();
        for (int row = 0; row < nStates; row++) {
            StringBuilder keyBuilder = new StringBuilder();
            for (int col = 0; col < nCols; col++) {
                keyBuilder.append((long) Math.round(SSq.get(row, col))).append(',');
            }
            String key = keyBuilder.toString();
            Double acc = aggrProb.get(key);
            double p = pi.get(row);
            aggrProb.put(key, acc == null ? p : acc + p);
            if (!aggrRep.containsKey(key)) {
                aggrRep.put(key, row);
            }
        }
        String bestKey = null;
        double bestProb = -1.0;
        for (Map.Entry<String, Double> entry : aggrProb.entrySet()) {
            if (entry.getValue() > bestProb) {
                bestProb = entry.getValue();
                bestKey = entry.getKey();
            }
        }
        Matrix mode = SSq.getRow(aggrRep.get(bestKey));

        // Map the aggregate-state columns (nclasses per stateful station, in station
        // order) onto the placement matrix; only service stations (Queue/Delay)
        // receive an initial population.
        Matrix placement = new Matrix(M, K);
        placement.zero();
        int col = 0;
        for (int i = 0; i < snCtmc.nstations; i++) {
            int nodeIdx = (int) snCtmc.stationToNode.get(i);
            if (snCtmc.isstateful.get(nodeIdx) == 0.0) {
                continue;
            }
            NodeType nodeType = snCtmc.nodetype.get(nodeIdx);
            for (int r = 0; r < K; r++) {
                double v = mode.get(col++);
                if (i < M && (nodeType == NodeType.Queue || nodeType == NodeType.Delay)) {
                    placement.set(i, r, v);
                }
            }
        }
        return placement;
    }

    /**
     * Placement from the steady-state mean queue lengths of a generic network
     * solver: floor the per-station means and distribute the residual
     * closed-class jobs by largest remainder so each closed population is
     * conserved.
     */
    private static Matrix placementFromMeanQLen(NetworkSolver initSolver, NetworkStruct sn) {
        Matrix qlen = initSolver.getAvgQLen();
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix placement = new Matrix(M, K);
        placement.zero();
        for (int r = 0; r < K; r++) {
            double njobs = sn.njobs.get(r);
            boolean isClosed = !Double.isInfinite(njobs);
            int[] floors = new int[M];
            double[] fracs = new double[M];
            boolean[] eligible = new boolean[M];
            int placed = 0;
            for (int i = 0; i < M; i++) {
                int nodeIdx = (int) sn.stationToNode.get(i);
                NodeType nodeType = sn.nodetype.get(nodeIdx);
                eligible[i] = (nodeType == NodeType.Queue || nodeType == NodeType.Delay);
                if (!eligible[i]) {
                    continue;
                }
                double mean = Math.max(0.0, qlen.get(i, r));
                if (isClosed) {
                    floors[i] = (int) Math.floor(mean);
                    fracs[i] = mean - floors[i];
                    placed += floors[i];
                } else {
                    floors[i] = (int) Math.round(mean);
                }
                placement.set(i, r, floors[i]);
            }
            if (isClosed) {
                // Largest-remainder apportionment of the residual jobs.
                int residual = (int) Math.round(njobs) - placed;
                while (residual > 0) {
                    int bestStation = -1;
                    double bestFrac = -1.0;
                    for (int i = 0; i < M; i++) {
                        if (eligible[i] && fracs[i] > bestFrac) {
                            bestFrac = fracs[i];
                            bestStation = i;
                        }
                    }
                    if (bestStation < 0) {
                        // All remainders consumed: place the leftover jobs at the
                        // reference station so the population is conserved.
                        int refStation = (int) sn.refstat.get(r);
                        if (!eligible[refStation]) {
                            for (int i = 0; i < M; i++) {
                                if (eligible[i]) {
                                    refStation = i;
                                    break;
                                }
                            }
                        }
                        placement.set(refStation, r, placement.get(refStation, r) + residual);
                        residual = 0;
                        break;
                    }
                    placement.set(bestStation, r, placement.get(bestStation, r) + 1);
                    fracs[bestStation] = -1.0;
                    residual--;
                }
            }
        }
        return placement;
    }
}
