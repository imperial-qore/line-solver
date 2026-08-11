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
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.*;
import jline.lang.nodes.*;
import jline.lang.processes.Disabled;
import jline.lang.processes.DiscreteSampler;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.processes.Zipf;
import jline.solvers.nc.NC;
import jline.solvers.ssa.SSA;
import jline.util.Maths;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Examples of caching models
 */
public class CacheModel {

    /**
     * Basic open cache model with Round Robin (RR) replacement strategy.
     * <p>
     * Features:
     * - Cache with 5 items, capacity 2, RR replacement
     * - Three classes: InitClass (requests), HitClass, MissClass
     * - Zipf access pattern with alpha=1.4 (skewed popularity)
     * - Simple Source → Cache → Sink topology
     * - Exponential arrival process with rate 2
     *
     * @return configured cache network model
     */
    public static Network cache_replc_rr() {
        Network model = new Network("model");

        int n = 5; // Number of items
        int m = 2;

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.RR);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(2));

        // Zipf-like item references
        Zipf pAccess = new Zipf(1.4, n);

        cacheNode.setRead(jobClass, pAccess);

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Closed cache model with feedback from hits and misses.
     * <p>
     * Features:
     * - Cache with 5 items, capacity 2, FIFO replacement
     * - Closed system with 1 job circulating
     * - Delay node with exponential service (rate 1.0)
     * - Both hits and misses return to delay node as same class
     * - Uniform access pattern across all items
     *
     * @return configured closed cache model
     */
    public static Network cache_replc_fifo() {
        Network model = new Network("model");

        int n = 5; // Number of items
        int m = 2; // Cache capacity

        Delay delay = new Delay(model, "Delay");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.FIFO);

        ClosedClass jobClass = new ClosedClass(model, "JobClass", 1, delay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, delay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, delay, 0);

        delay.setService(jobClass, new Exp(1.0));

        Matrix p = new Matrix(1, n).fill(1.0 / n);
        DiscreteSampler pAccess = new DiscreteSampler(p);
        cacheNode.setRead(jobClass, pAccess);

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, delay, cacheNode, 1.0);
        routingMatrix.set(hitClass, jobClass, cacheNode, delay, 1.0);
        routingMatrix.set(missClass, jobClass, cacheNode, delay, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Cache model with multiple delay nodes and random routing.
     * <p>
     * Features:
     * - Cache with 5 items, capacity 2, LRU replacement
     * - Router node directing hits/misses to two different delay nodes
     * - Different service rates for hits vs misses at each delay
     * - Random routing strategy for load balancing
     * - Demonstrates cache integration with complex topologies
     *
     * @return configured cache network with routing
     */
    public static Network cache_replc_routing() {
        Network model = new Network("model");

        int n = 5; // Number of items
        int m = 2;

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.FIFO);
        Router routerNode = new Router(model, "Router");
        Delay delay1 = new Delay(model, "Delay1");
        Delay delay2 = new Delay(model, "Delay2");
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(2));
        source.setArrival(hitClass, new Disabled());
        source.setArrival(missClass, new Disabled());

        delay1.setService(hitClass, new Exp(10));
        delay1.setService(missClass, new Exp(1.0));

        delay2.setService(hitClass, new Exp(20));
        delay2.setService(missClass, new Exp(2));

        Matrix p = new Matrix(1, n).fill(1.0 / n);
        DiscreteSampler pAccess = new DiscreteSampler(p);

        cacheNode.setRead(jobClass, pAccess);
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        model.addLink(source, cacheNode);
        model.addLink(cacheNode, routerNode);
        model.addLink(routerNode, delay1);
        model.addLink(routerNode, delay2);
        model.addLink(delay1, sink);
        model.addLink(delay2, sink);

        source.setProbRouting(jobClass, cacheNode, 1.0);

        cacheNode.setProbRouting(hitClass, routerNode, 1.0);
        cacheNode.setProbRouting(missClass, routerNode, 1.0);

        routerNode.setRouting(hitClass, RoutingStrategy.RAND);
        routerNode.setRouting(missClass, RoutingStrategy.RAND);

        delay1.setProbRouting(hitClass, sink, 1.0);
        delay1.setProbRouting(missClass, sink, 1.0);

        delay2.setProbRouting(hitClass, sink, 1.0);
        delay2.setProbRouting(missClass, sink, 1.0);

        return model;
    }

    /**
     * Cache model with Zipf access pattern and Round Robin replacement.
     * <p>
     * Features:
     * - Cache with 5 items, multi-level capacity [2,1], RR replacement
     * - Zipf distribution for realistic access patterns (skewed popularity)
     * - Alpha parameter controls skewness (1.0 = moderately skewed)
     * - Round Robin replacement strategy instead of FIFO
     * - Exponential arrivals with rate 1
     *
     * @return configured cache model with Zipf access
     */
    public static Network cache_compare_replc() {
        Network model = new Network("model");

        int n = 5; // Number of items
        Matrix m = new Matrix("[2,1]");
        double alpha = 1.0;

        Source source = new Source(model, "Source");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.RR, null);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(2));

        Zipf pAccess = new Zipf(alpha, n);

        cacheNode.setRead(jobClass, pAccess);
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Layered cache queueing model example 1.
     * <p>
     * Features:
     * - Layered network with client processor and cache processor
     * - Client task (T1) with reference scheduling and PS processor
     * - Cache task (C2) with 4 items, capacity 2, Round Robin replacement
     * - Item entry (I2) with uniform access pattern
     * - Activities: A1 (client), AC2 (cache access), AC2h (hit), AC2m (miss)
     * - Cache access precedence with hit and miss paths
     * - LN with MVA backend for solution
     *
     * @return configured layered cache queueing model
     */
    public static LayeredNetwork lcq_singlehost() {
        LayeredNetwork model = new LayeredNetwork("cacheInLayeredNetwork");

        // Client processor and task
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "E1").on(T1);

        // Cache processor and task
        int totalitems = 4;
        int cachecapacity = 2;
        Matrix pAccess = new Matrix(1, totalitems);
        pAccess.fill(1.0 / totalitems);
        DiscreteSampler discreteSampler = new DiscreteSampler(pAccess);
        
        Processor PC = new Processor(model, "PC", 1, SchedStrategy.PS);
        CacheTask C2 = new CacheTask(model, "C2");
        C2.on(PC);
        ItemEntry I2 = new ItemEntry(model, "I2", totalitems, discreteSampler).on(C2);

        // Activities
        Activity A1 = new Activity(model, "A1", Immediate.getInstance()).on(T1);
        A1.boundTo(E1);
        A1.synchCall(I2, 1);

        Activity AC2 = new Activity(model, "AC2", Immediate.getInstance()).on(C2);
        AC2.boundTo(I2);
        
        Activity AC2h = new Activity(model, "AC2h", new Exp(1.0)).on(C2);
        AC2h.repliesTo(I2);
        
        Activity AC2m = new Activity(model, "AC2m", new Exp(0.5)).on(C2);
        AC2m.repliesTo(I2);

        // Cache access precedence
        List<Activity> cacheActivities = new ArrayList<>();
        cacheActivities.add(AC2h);
        cacheActivities.add(AC2m);
        C2.addPrecedence(ActivityPrecedence.CacheAccess(AC2, cacheActivities));

        return model;
    }

    /**
     * Closed cache model with LRU replacement strategy.
     * <p>
     * Features:
     * - Cache with 5 items, capacity 2, LRU replacement
     * - Closed system with 1 job circulating
     * - Delay node with exponential service (rate 1.0)
     * - Both hits and misses return to delay node as same class
     * - Uniform access pattern across all items
     *
     * @return configured closed cache model with LRU
     */
    public static Network cache_replc_lru() {
        Network model = new Network("model");

        int n = 5; // Number of items
        int m = 2; // Cache capacity

        Delay delay = new Delay(model, "Delay");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.LRU);

        ClosedClass jobClass = new ClosedClass(model, "JobClass", 1, delay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, delay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, delay, 0);

        delay.setService(jobClass, new Exp(1.0));

        Matrix p = new Matrix(1, n).fill(1.0 / n);
        DiscreteSampler pAccess = new DiscreteSampler(p);
        cacheNode.setRead(jobClass, pAccess);

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, delay, cacheNode, 1.0);
        routingMatrix.set(hitClass, jobClass, cacheNode, delay, 1.0);
        routingMatrix.set(missClass, jobClass, cacheNode, delay, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Layered cache queueing model example 2 with multi-level cache and downstream service.
     * <p>
     * Features:
     * - Layered network with client, cache, and downstream service processors
     * - Client task (T1) with 1 user, reference scheduling
     * - Cache task (CT) with 4 items, multi-level capacity [1,1], Round Robin replacement
     * - Item entry (IE) with uniform access pattern
     * - Downstream service task (T2) with FCFS scheduling and exponential service
     * - Cache miss calls downstream service synchronously
     * - Activities: A1 (client), Ac (cache access), Ac_hit (hit), Ac_miss (miss with service call)
     * - Cache access precedence with hit and miss paths
     * - LN with NC and MVA backends for solution
     *
     * @return configured layered cache queueing model with downstream service
     */
    public static LayeredNetwork lcq_threehosts() {
        LayeredNetwork model = new LayeredNetwork("LQNwithCaching");

        int nusers = 1;
        int ntokens = 1;

        // Client processor and task
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", nusers, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "E1").on(T1);

        // Cache processor and task
        int totalitems = 4;
        Matrix cachecapacity = new Matrix("[1,1]");
        Matrix pAccess = new Matrix(1, totalitems);
        pAccess.fill(1.0 / totalitems);
        DiscreteSampler discreteSampler = new DiscreteSampler(pAccess);
        
        Processor PC = new Processor(model, "Pc", 1, SchedStrategy.PS);
        CacheTask CT = new CacheTask(model, "CT");
        CT.on(PC);
        ItemEntry IE = new ItemEntry(model, "IE", totalitems, discreteSampler).on(CT);

        // Downstream service processor and task
        Processor P3 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Task T3 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P3);
        Entry E3 = new Entry(model, "E2").on(T3);
        Activity A3 = new Activity(model, "A2", new Exp(5.0)).on(T3);
        A3.boundTo(E3);
        A3.repliesTo(E3);

        // Client activity
        Activity A1 = new Activity(model, "A1", Immediate.getInstance()).on(T1);
        A1.boundTo(E1);
        A1.synchCall(IE, 1);

        // Cache activities
        Activity AC2 = new Activity(model, "Ac", Immediate.getInstance()).on(CT);
        AC2.boundTo(IE);
        
        Activity AC2h = new Activity(model, "Ac_hit", new Exp(1.0)).on(CT);
        AC2h.repliesTo(IE);
        
        Activity AC2m = new Activity(model, "Ac_miss", new Exp(0.5)).on(CT);
        AC2m.synchCall(E3, 1);
        AC2m.repliesTo(IE);

        // Cache access precedence
        List<Activity> cacheActivities = new ArrayList<>();
        cacheActivities.add(AC2h);
        cacheActivities.add(AC2m);
        CT.addPrecedence(ActivityPrecedence.CacheAccess(AC2, cacheActivities));

        return model;
    }

    /**
     * Layered cache queueing model with asynchronous (non-blocking) cache access.
     * <p>
     * Features:
     * - Client makes async call to cache (fire-and-forget, non-blocking)
     * - Client continues immediately without waiting for cache response
     * - Cache still uses POST_CACHE precedence for hit/miss determination
     * - Demonstrates async cache access pattern for prefetching scenarios
     * - Based on lcq_singlehost() but with asynchCall instead of synchCall
     *
     * @return configured layered cache queueing model with async cache access
     */
    public static LayeredNetwork lcq_async_prefetch() {
        LayeredNetwork model = new LayeredNetwork("AsyncCachePrefetch");

        // Client processor and task
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "E1").on(T1);

        // Cache processor and task
        int totalitems = 4;
        int cachecapacity = 2;
        Matrix pAccess = new Matrix(1, totalitems);
        pAccess.fill(1.0 / totalitems);
        DiscreteSampler discreteSampler = new DiscreteSampler(pAccess);

        Processor PC = new Processor(model, "PC", 1, SchedStrategy.PS);
        CacheTask C2 = new CacheTask(model, "C2", totalitems, cachecapacity, ReplacementStrategy.LRU, 1);
        C2.on(PC);
        ItemEntry I2 = new ItemEntry(model, "I2", totalitems, discreteSampler).on(C2);

        // Client activity with ASYNC call to cache
        Activity A1 = new Activity(model, "A1", Immediate.getInstance()).on(T1);
        A1.boundTo(E1);
        A1.asynchCall(I2, 1); // ASYNC - fire and forget, non-blocking

        // Cache activities (unchanged from sync version)
        Activity AC2 = new Activity(model, "AC2", Immediate.getInstance()).on(C2);
        AC2.boundTo(I2);

        Activity AC2h = new Activity(model, "AC2h", new Exp(1.0)).on(C2);
        AC2h.repliesTo(I2);

        Activity AC2m = new Activity(model, "AC2m", new Exp(0.5)).on(C2);
        AC2m.repliesTo(I2);

        // Cache access precedence (unchanged)
        List<Activity> cacheActivities = new ArrayList<>();
        cacheActivities.add(AC2h);
        cacheActivities.add(AC2m);
        C2.addPrecedence(ActivityPrecedence.CacheAccess(AC2, cacheActivities));

        return model;
    }

    /**
     * Comparison of synchronous vs asynchronous cache access patterns.
     * <p>
     * Creates two models:
     * - Sync version: Client blocks waiting for cache response
     * - Async version: Client continues without waiting (fire-and-forget)
     * <p>
     * Use this to compare:
     * - Client response time (async should be lower)
     * - Client throughput (async should be higher)
     * - Cache hit/miss ratios (should be identical)
     *
     * @return configured layered cache queueing model for comparison
     */
    public static LayeredNetwork lcq_async_vs_sync_comparison() {
        // Return the sync version (lcq_singlehost) for comparison
        // Users can compare with lcq_async_prefetch() manually
        return lcq_singlehost();
    }

    /**
     * Main method for testing and demonstrating cache model examples.
     *
     * <p>Currently configured to:
     * - Set MATLAB-compatible random number generation
     * - Test both regular cache models and layered cache queueing models
     * - Solve using multiple solvers for comparison
     * - Measure and display execution time
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        Maths.setRandomNumbersMatlab(true);

        System.out.println("=== Testing cache_replc_rr ===");
        Network model1 = cache_replc_rr();
        
        try {
            System.out.println("--- NC Solver ---");
            new NC(model1).getAvgNodeTable().print();
        } catch (Exception e) {
            System.out.println("NC Solver failed: " + e.getMessage());
        }
        
        try {
            System.out.println("--- MVA Solver ---");
            new jline.solvers.mva.MVA(model1).getAvgNodeTable().print();
        } catch (Exception e) {
            System.out.println("MVA Solver failed: " + e.getMessage());
        }

        System.out.println("\n=== Testing cache_replc_fifo ===");
        Network model2 = cache_replc_fifo();
        
        try {
            System.out.println("--- NC Solver ---");
            new NC(model2).getAvgNodeTable().print();
        } catch (Exception e) {
            System.out.println("NC Solver failed: " + e.getMessage());
        }
        
        try {
            System.out.println("--- MVA Solver ---");
            new jline.solvers.mva.MVA(model2).getAvgNodeTable().print();
        } catch (Exception e) {
            System.out.println("MVA Solver failed: " + e.getMessage());
        }

        System.out.println("\n=== Testing lcq_singlehost ===");
        LayeredNetwork layeredModel1 = lcq_singlehost();
        
        try {
            System.out.println("--- Layered Network Model Created ---");
            System.out.println("Model: " + layeredModel1.getName());
            System.out.println("Layered cache model 1 created successfully");
            // Note: LN would be used here for solving layered networks
        } catch (Exception e) {
            System.out.println("Layered cache model 1 failed: " + e.getMessage());
        }

        System.out.println("\n=== Testing lcq_threehosts ===");
        LayeredNetwork layeredModel2 = lcq_threehosts();
        
        try {
            System.out.println("--- Layered Network Model Created ---");
            System.out.println("Model: " + layeredModel2.getName());
            System.out.println("Layered cache model 2 created successfully");
            // Note: LN would be used here for solving layered networks
        } catch (Exception e) {
            System.out.println("Layered cache model 2 failed: " + e.getMessage());
        }
        
        System.out.println("\n=== Testing cache_replc_routing ===");
        Network model3 = cache_replc_routing();
        
        System.out.println("\n=== Testing cache_replc_lru ===");
        Network model4 = cache_replc_lru();
        
        try {
            System.out.println("--- SSA Solver for routing model ---");
            SSA solver = new SSA(model3, "samples", 1000, "verbose", true, "method", "serial", "seed", 1);
            solver.getAvgNodeTable().print();
            System.out.println("SSA Solver succeeded");
        } catch (Exception e) {
            System.out.println("SSA Solver failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            System.out.println("--- NC Solver for LRU model ---");
            new NC(model4).getAvgNodeTable().print();
        } catch (Exception e) {
            System.out.println("NC Solver failed: " + e.getMessage());
        }
        
        try {
            System.out.println("--- MVA Solver for LRU model ---");
            new jline.solvers.mva.MVA(model4).getAvgNodeTable().print();
        } catch (Exception e) {
            System.out.println("MVA Solver failed: " + e.getMessage());
        }
    }
}
