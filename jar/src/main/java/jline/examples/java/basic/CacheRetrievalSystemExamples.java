/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.nodes.Cache;
import jline.solvers.NetworkAvgNodeTable;
import jline.solvers.ctmc.CTMC;
import jline.solvers.ssa.SSA;
import jline.util.matrix.Matrix;

/**
 * Runnable examples for cache models with a retrieval system (delayed hits).
 *
 * <p>Each example solves a model built by {@link CacheRetrievalSystemModel} with
 * the CTMC and SSA solvers and prints the cache hit/miss ratios and the expected
 * latency reported by {@link Cache#getResidT()}.</p>
 *
 * @see CacheRetrievalSystemModel
 */
public class CacheRetrievalSystemExamples {

    /**
     * Single-queue retrieval system over 3 items, cache capacity 1.
     * Solved with both CTMC (exact) and SSA (simulation).
     */
    public static void retrieval_simple() {
        Network model = CacheRetrievalSystemModel.simple_retrieval_system_model(
                new double[] {0.6, 0.3, 0.1},
                new double[] {2.0, 2.0, 2.0},
                "[1]", SchedStrategy.INF);
        solveCtmc(model);
        solveSsa(model);
    }

    /**
     * PS variant of retrieval_simple, with the paper's motivating-example parameters.
     * Solved with SSA — the CTMC state space is too large for an exact solve; the
     * exact hit rate is available via the product-form recurrence (api/retrieval).
     */
    public static void retrieval_ps() {
        double[] accessProb = {49, 49, 49, 49, 7, 1, 1};
        double total = 0.0;
        for (double p : accessProb) total += p;
        for (int i = 0; i < accessProb.length; i++) accessProb[i] /= total;   // lambda = (49,49,49,49,7,1,1)/205
        Network model = CacheRetrievalSystemModel.simple_retrieval_system_model(
                accessProb,
                new double[] {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},             // mu_i = 1
                "[6]", SchedStrategy.PS, ReplacementStrategy.RR);
        solveSsa(model);
    }

    /**
     * Two-queue chained retrieval system over 3 items, cache capacity 1.
     * Solved with SSA — the CTMC state space is too large for an exact solve.
     */
    public static void retrieval_chain() {
        Matrix serviceRates = new Matrix(new double[][] {
                {2.0, 2.0, 2.0},  // Queue_1 per-item rates
                {3.0, 3.0, 3.0}   // Queue_2 per-item rates
        });
        Network model = CacheRetrievalSystemModel.chain_retrieval_system_model(
                new double[] {0.6, 0.3, 0.1}, serviceRates, 2, "[1]", SchedStrategy.FCFS);
        solveSsa(model);
    }

    /**
     * Probabilistic per-item routing through a 3-queue retrieval system, cache capacity 2.
     * Solved with SSA — the CTMC state space is too large for an exact solve.
     */
    public static void retrieval_routing() {
        Matrix serviceRates = new Matrix(new double[][] {
                {2.0, 2.0, 2.0},  // IS queue
                {3.0, 3.0, 3.0},  // Queue 1
                {3.0, 3.0, 3.0}   // Queue 2
        });
        Network model = CacheRetrievalSystemModel.retrieval_system_with_probabilistic_routing(
                new double[] {0.6, 0.3, 0.1}, serviceRates, "[2]", SchedStrategy.FCFS);
        solveSsa(model);
    }

    private static void solveCtmc(Network model) {
        try {
            CTMC ctmc = new CTMC(model, "cutoff", 1, "seed", 1);
            NetworkAvgNodeTable ctmcTable = ctmc.getAvgNodeTable();
            ctmcTable.print();
            Cache cache = (Cache) model.getNodeByName("Cache");
            System.out.println("CTMC  hit=" + cache.getHitRatio()
                    + " miss=" + cache.getMissRatio()
                    + " expectedLatency=" + cache.getResidT());
            model.reset();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void solveSsa(Network model) {
        try {
            SSA ssa = new SSA(model, "samples", 20000, "method", "serial", "seed", 1);
            NetworkAvgNodeTable ssaTable = ssa.getAvgNodeTable();
            ssaTable.print();
            Cache cache = (Cache) model.getNodeByName("Cache");
            System.out.println("SSA   hit=" + cache.getHitRatio()
                    + " miss=" + cache.getMissRatio()
                    + " expectedLatency=" + cache.getResidT());
            model.reset();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        retrieval_simple();
        retrieval_ps();
        retrieval_chain();
        retrieval_routing();
    }
}
