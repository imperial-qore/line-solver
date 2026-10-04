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
        solveAll(model);
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
        solveAll(model);
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
        solveAll(model);
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
        solveAll(model);
    }

    /**
     * Retrieval system inheriting per-item routing and service from the read
     * class, with item-level overrides. Solved with SSA and, analytically,
     * with MVA and NC.
     */
    public static void retrieval_default() {
        Network model = CacheRetrievalSystemModel.retrieval_default();
        solveAll(model);
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

    /**
     * The four solves every retrieval reference runs, on the CACHE table.
     *
     * <p>The seed and run length are the reference's own: the goldens for these
     * examples were recorded from that run, so a row simulated for a different
     * number of samples is a different measurement of the same model rather than
     * a disagreement about it.
     */
    private static void solveAll(Network model) {
        try {
            new SSA(model, "samples", 100000, "method", "serial", "seed", 1)
                    .getAvgCacheTable().print();
            model.reset();
        } catch (Exception e) {
            System.out.println("SSA failed: " + e.getMessage());
        }
        try {
            new jline.solvers.ldes.LDES(model, "samples", 1000000, "seed", 1)
                    .getAvgCacheTable().print();
            model.reset();
        } catch (Exception e) {
            System.out.println("LDES failed: " + e.getMessage());
        }
        try {
            new jline.solvers.mva.MVA(model).getAvgCacheTable().print();
            model.reset();
        } catch (Exception e) {
            System.out.println("MVA failed: " + e.getMessage());
        }
        try {
            new jline.solvers.nc.NC(model).getAvgCacheTable().print();
            model.reset();
        } catch (Exception e) {
            System.out.println("NC failed: " + e.getMessage());
        }
    }


    public static void main(String[] args) {
        retrieval_simple();
        retrieval_ps();
        retrieval_chain();
        retrieval_routing();
        retrieval_default();
    }
}
