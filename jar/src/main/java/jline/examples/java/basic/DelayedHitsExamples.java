/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.DiscreteSampler;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgNodeTable;
import jline.solvers.ctmc.CTMC;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

/**
 * Runnable examples for cache models with a retrieval system that exhibit
 * <em>delayed hits</em>: while an item is being fetched by the retrieval system
 * after a miss, subsequent requests for the same item neither hit nor miss
 * outright but are delayed until the in-flight retrieval completes.
 *
 * <p>Each example builds a {@link Network} with a {@link Cache} node whose miss
 * traffic is routed through a retrieval system, solves it with the CTMC solver,
 * prints the aggregate node table together with the generated state space and
 * infinitesimal generator, and reports the hit/miss ratios and expected latency
 * via {@link Cache#getHitRatio()}, {@link Cache#getMissRatio()} and
 * {@link Cache#getResidT()}.</p>
 *
 * @see CacheRetrievalSystemModel
 * @see CacheRetrievalSystemExamples
 */
public class DelayedHitsExamples {

    /**
     * Creates a caching model with a simple retrieval system, where there is only one queue in the retrieval system.
     */
    public static Network simple_retrieval_system_model() {
        Network model = new Network("simple model");

        int n = 2;
        int m = 1;

        Delay clientDelay = new Delay(model, "Client");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.RR);
        Queue queue = new Queue(model, "Queue", SchedStrategy.INF);

        // 2 jobs to allow hits/misses and delayed hits
        ClosedClass jobClass = new ClosedClass(model, "InitClass", 2, clientDelay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, clientDelay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, clientDelay, 0);

        clientDelay.setService(jobClass, new Exp(1.0));

        DiscreteSampler pAccess = new DiscreteSampler(new Matrix(new double[][]{{0.5, 0.5}}));

        cacheNode.setRead(jobClass, pAccess);
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        // Single-queue retrieval system: read-class service (rate 1.0) and cache <-> queue
        // routing for every item, wired via the per-item methods.
        queue.setService(jobClass, new Exp(1.0));
        cacheNode.setRetrievalSystem(jobClass, missClass, new Queue[] {queue});
        for (int item = 0; item < n; item++) {
            cacheNode.setItemRoutingProb(jobClass, item, cacheNode, queue, 1.0); // cache -> queue
            cacheNode.setItemRoutingProb(jobClass, item, queue, cacheNode, 1.0);  // queue -> cache
        }

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, clientDelay, cacheNode, 1.0);
        routingMatrix.set(hitClass, jobClass, cacheNode, clientDelay, 1.0);
        routingMatrix.set(missClass, jobClass, cacheNode, clientDelay, 1.0);

        model.link(routingMatrix);

        return model;
    }

    public static Network chain_retrieval_system_model() {
        Network model = new Network("model");

        int n = 3;
        int m = 1;

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.FIFO);
        Queue queue_1 = new Queue(model, "Queue Item 1", SchedStrategy.PS);
        Queue queue_2 = new Queue(model, "Queue Item 2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        Queue[] allQueues = new Queue[]{queue_1, queue_2};

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(1));

        DiscreteSampler pAccess = new DiscreteSampler(new Matrix(new double[][]{{0.6, 0.3, 0.1}}));

        cacheNode.setRead(jobClass, pAccess);
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        // Read-class service at each queue = default per-item retrieval service.
        for (Queue queue : allQueues) {
            queue.setService(jobClass, new Exp(1.0));
        }
        cacheNode.setRetrievalSystem(jobClass, missClass, allQueues);

        // Per-item service override (redundant here, kept to exercise the API).
        for (Queue queue : allQueues) {
            for (int i = 0; i < n; i++) {
                queue.setItemServiceRate(cacheNode, jobClass, i, 1.0);
            }
        }

        // All items enter at queue 1
        cacheNode.setItemRoutingProb(jobClass, 0, cacheNode, queue_1, 1);
        cacheNode.setItemRoutingProb(jobClass, 1, cacheNode, queue_1, 1);
        cacheNode.setItemRoutingProb(jobClass, 2, cacheNode, queue_1, 1);

        // After service at queue 1, route to queue 2
        cacheNode.setItemRoutingProb(jobClass, 0, queue_1, queue_2, 1.0);
        cacheNode.setItemRoutingProb(jobClass, 1, queue_1, queue_2, 1.0);
        cacheNode.setItemRoutingProb(jobClass, 2, queue_1, queue_2, 1.0);

        // After service at queue 2, always exit
        cacheNode.setItemRoutingProb(jobClass, 0, queue_2, cacheNode, 1.0);
        cacheNode.setItemRoutingProb(jobClass, 1, queue_2, cacheNode, 1.0);
        cacheNode.setItemRoutingProb(jobClass, 2, queue_2, cacheNode, 1.0);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);
        model.link(routingMatrix);

        return model;
    }

    public static void run_solvers(Network model) {
        System.out.println("=== Running Solvers ===");

        try {
            CTMC solver1 = new CTMC(model, "keep", false, "cutoff", 1, "seed", 1);
            NetworkAvgNodeTable avgTable1 = solver1.getAvgNodeTable();
            System.out.println("--- CTMC Solver ---");
            avgTable1.print();

            SolverCTMC.StateSpace stateSpace = solver1.getStateSpace();
            SolverCTMC.generatorResult infGen = solver1.getInfGen();

            System.out.println("stateSpace: \n" + stateSpace.stateSpace);

            SolverCTMC.printInfGen(infGen, stateSpace);

            Cache cacheNode = (Cache) model.getNodeByName("Cache");
            double hitRatio = cacheNode.getHitRatio().get(0);
            double missRatio = cacheNode.getMissRatio().get(0);
            double expectedLatency = cacheNode.getResidT().get(0);
            System.out.println("Hit Ratio: " + hitRatio);
            System.out.println("Miss Ratio: " + missRatio);
            System.out.println("Expected Latency: " + expectedLatency);

        } catch (Exception e) {
            System.out.println("Error running solvers: " + e.getMessage());
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        System.out.println("--- Simple Model ---");
        Network simpleModel = simple_retrieval_system_model();
        run_solvers(simpleModel);

        System.out.println("--- Chain Model ---");
        Network chainModel = chain_retrieval_system_model();
        run_solvers(chainModel);
    }
}
