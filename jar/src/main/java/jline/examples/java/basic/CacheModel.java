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
import jline.lib.rmf.CacheRMF;
import jline.solvers.SolverOptions;
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
     * Cache with CLIMB (transposition) replacement, closed model over 5 items.
     *
     * <p>On a hit an item moves up one position; on a miss it enters at the tail.
     * Exact in CTMC, simulated in SSA/LDES. Not product-form, so MVA/NC/FLD
     * reject it.</p>
     *
     * @return configured cache network model
     */
    public static Network cache_replc_climb() {
        Network model = new Network("model");

        int n = 5; // Number of items
        int m = 2; // Cache capacity

        Delay delay = new Delay(model, "Delay");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.CLIMB);

        ClosedClass jobClass = new ClosedClass(model, "JobClass", 1, delay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, delay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, delay, 0);

        delay.setService(jobClass, new Exp(1.0));
        cacheNode.setRead(jobClass, new Zipf(1.2, n));
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
     * Cache with q-LRU replacement, closed model over 5 items.
     *
     * <p>On a miss the item is admitted (LRU head insert) with probability q,
     * otherwise it passes through uncached. Admission filtering can raise the
     * hit ratio over plain LRU under skewed popularity. Exact in CTMC,
     * simulated in SSA/LDES. Not product-form, so MVA/NC/FLD reject it.</p>
     *
     * @return configured cache network model
     */
    public static Network cache_replc_qlru() {
        Network model = new Network("model");

        int n = 5; // Number of items
        int m = 2; // Cache capacity

        Delay delay = new Delay(model, "Delay");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.QLRU);
        cacheNode.setAdmissionProb(0.5); // admit a missed item with probability q=0.5

        ClosedClass jobClass = new ClosedClass(model, "JobClass", 1, delay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, delay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, delay, 0);

        delay.setService(jobClass, new Exp(1.0));
        cacheNode.setRead(jobClass, new Zipf(1.2, n));
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
     * Cache with h-LRU / LRU(m) replacement, open model over 6 items.
     *
     * <p>h LRU lists of capacities m[0..h-1]; a miss inserts the item at the head
     * of list 1, a hit in list l exchanges the item with the tail of list l+1.
     * Exact in CTMC, simulated in SSA/LDES; MVA uses the characteristic-time
     * (TTL) approximation of Gast and Van Houdt (SIGMETRICS 2015), which reduces
     * to the Che approximation for h=1.</p>
     *
     * @return configured cache network model
     */
    public static Network cache_replc_hlru() {
        Network model = new Network("model");

        int n = 6; // Number of items
        Matrix m = new Matrix(1, 2); // list 1 holds 2 items, list 2 holds 1
        m.set(0, 0, 2);
        m.set(0, 1, 1);

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.HLRU);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(1.0));
        cacheNode.setRead(jobClass, new Zipf(1.2, n));
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
     * Cache with per-item storage costs (sizes) and per-list cost caps.
     *
     * <p>Each item i carries a storage cost sigma_i and list j may hold items of
     * total cost at most k_j. A promotion that would breach a cap serves the
     * request without changing the cache state. SolverNC evaluates the
     * constrained normalizing constant E(m,k) of Casale-Gast (IEEE/ACM ToN
     * 29(2), 2021), Sec. IX; SolverLDES simulates the same rule directly.</p>
     *
     * @return configured cache network model
     */
    public static Network cache_itemsize_costcap() {
        Network model = new Network("model");

        int n = 6;                   // number of items
        Matrix m = new Matrix(1, 2); // two lists, one item each
        m.set(0, 0, 1);
        m.set(0, 1, 1);

        Matrix sizes = new Matrix(1, n); // small and large items
        double[] sizeValues = {1, 1, 1, 2, 2, 2};
        for (int i = 0; i < n; i++) sizes.set(0, i, sizeValues[i]);
        Matrix caps = new Matrix(1, 2); // list 2 admits small items only
        caps.set(0, 0, 2);
        caps.set(0, 1, 1);

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.RR);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(2.0));
        cacheNode.setRead(jobClass, new DiscreteSampler(new Matrix(1, n).fill(1.0 / n)));
        cacheNode.setItemSizes(sizes);
        cacheNode.setCostCaps(caps);
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
     * Refined mean field (RMF) transient and steady-state analysis of a
     * two-list RANDOM(m) cache, driven through {@link CacheRMF} directly.
     *
     * <p>Prints the steady-state hit and miss probabilities with the 1/N
     * correction, then the transient evolution X(t) + V(t)/N of the hit rate.
     * This example has no Network: RMF is a mean-field limit of the cache
     * chain, not a queueing model, so there is nothing to hand a solver.</p>
     */
    public static void cache_rmf_transient() {
        int n = 10;             // number of items
        int[] m = {3, 2};       // list capacities (2-list cache)
        double alpha = 0.8;     // Zipf exponent

        double[] p = new double[n];
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            p[i] = Math.pow(i + 1, -alpha);
            sum += p[i];
        }
        for (int i = 0; i < n; i++) p[i] /= sum;

        System.out.println("Cache parameters: n=" + n + ", m=[" + m[0] + ", " + m[1]
                + "], Zipf(" + alpha + ")");

        CacheRMF rmf = new CacheRMF(p, m);

        Object[] ss = rmf.meanFieldExpansionSteadyState(1);
        double[] pi = (double[]) ss[0];
        double[] V = (double[]) ss[1];
        double[] piRefined = new double[pi.length];
        for (int i = 0; i < pi.length; i++) piRefined[i] = pi[i] + V[i] / n;

        System.out.println("\nSteady-state results (refined mean field):");
        double totalHit = 0.0;
        for (int k = 1; k <= m.length; k++) {
            double hr = rmf.hitRate(piRefined, k);
            totalHit += hr;
            // The golden holds these under CACHE, one row per quantity; nothing a
            // result table carries can answer them, so the example records its own.
            jline.io.LineResultRecorder.scalar("CACHE", "HitRate_L" + k, "Steady", hr);
            System.out.printf("  Hit rate (list %d): %.6f%n", k, hr);
        }
        double missRate = rmf.hitRate(piRefined, 0);
        jline.io.LineResultRecorder.scalar("CACHE", "MissRate", "Steady", missRate);
        jline.io.LineResultRecorder.scalar("CACHE", "TotalHitProb", "Steady", totalHit);
        jline.io.LineResultRecorder.scalar("CACHE", "TotalMissProb", "Steady", missRate);
        System.out.printf("  Miss rate:         %.6f%n", missRate);
        System.out.printf("  Total hit prob:    %.6f%n", totalHit);
        System.out.printf("  Total miss prob:   %.6f%n", missRate);

        Object[] tr = rmf.meanFieldExpansionTransient(50.0, 200, 1);
        double[] T = (double[]) tr[0];
        double[][] X = (double[][]) tr[1];
        double[][] Vt = (double[][]) tr[2];

        System.out.println("\nTransient hit rates (refined, N=" + n + "):");
        int[] timeIndices = {0, 20, 50, 100, 199}; // t=0, ~5, ~12.5, ~25, 50
        for (int idx : timeIndices) {
            double[] xt = new double[X[idx].length];
            for (int i = 0; i < xt.length; i++) xt[i] = X[idx][i] + Vt[idx][i] / n;
            double hr = 0.0;
            for (int k = 1; k <= m.length; k++) hr += rmf.hitRate(xt, k);
            jline.io.LineResultRecorder.scalar("CACHE",
                    String.format("t=%.3f", T[idx]), "Transient", hr);
            System.out.printf("  t=%7.3f: hit_rate=%.6f%n", T[idx], hr);
        }
    }

    /**
     * MMAP-fed small RR cache with two correlated classes.
     *
     * <p>A marked MMPP2 arrival stream feeds a small Round-Robin cache. Its two
     * marks are bound to two open read classes that share the modulating chain,
     * so the classes are cross-correlated and autocorrelated in time, and each
     * reads the cache with a DIFFERENT item popularity. Phase 1 (bursty) emits
     * mostly class-1 references at a high rate, phase 2 (calm) mostly class-2 at
     * a low rate, so the shared chain couples "which class arrives" with "how
     * fast requests arrive".</p>
     *
     * @param lambda1 Read1 arrival rate, or a non-positive value to bind the MMAP
     * @param lambda2 Read2 arrival rate, ignored when the MMAP is bound
     * @param mmap    the marked arrival process, or null for phase-conditional Poisson
     * @return configured cache network model
     */
    public static Network cache_mmap_rr_env(double lambda1, double lambda2,
                                            jline.lang.processes.MarkedMAP mmap) {
        Network model = new Network("MMAPCache");

        int n = 4; // number of items
        int m = 2; // cache capacity

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.RR);
        Sink sink = new Sink(model, "Sink");

        OpenClass rd1 = new OpenClass(model, "Read1", 0);
        OpenClass rd2 = new OpenClass(model, "Read2", 0);
        OpenClass hit1 = new OpenClass(model, "Hit1", 0);
        OpenClass mis1 = new OpenClass(model, "Miss1", 0);
        OpenClass hit2 = new OpenClass(model, "Hit2", 0);
        OpenClass mis2 = new OpenClass(model, "Miss2", 0);

        // Per-class item-popularity distributions (deliberately different)
        Matrix p1 = new Matrix(1, n);   // class 1 favors low-index items
        Matrix p2 = new Matrix(1, n);   // class 2 favors high-index items
        double[] w1 = {8, 4, 2, 1};
        for (int i = 0; i < n; i++) {
            p1.set(0, i, w1[i] / 15.0);
            p2.set(0, i, w1[n - 1 - i] / 15.0);
        }
        cacheNode.setRead(rd1, new DiscreteSampler(p1));
        cacheNode.setRead(rd2, new DiscreteSampler(p2));
        cacheNode.setHitClass(rd1, hit1);
        cacheNode.setMissClass(rd1, mis1);
        cacheNode.setHitClass(rd2, hit2);
        cacheNode.setMissClass(rd2, mis2);

        if (mmap != null) {
            List<jline.lang.JobClass> marks = new ArrayList<jline.lang.JobClass>();
            marks.add(rd1);
            marks.add(rd2);
            source.setMarkedArrival(mmap, marks);
        } else {
            source.setArrival(rd1, new Exp(lambda1));
            source.setArrival(rd2, new Exp(lambda2));
        }

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(rd1, rd1, source, cacheNode, 1.0);
        routingMatrix.set(rd2, rd2, source, cacheNode, 1.0);
        routingMatrix.set(hit1, hit1, cacheNode, sink, 1.0);
        routingMatrix.set(mis1, mis1, cacheNode, sink, 1.0);
        routingMatrix.set(hit2, hit2, cacheNode, sink, 1.0);
        routingMatrix.set(mis2, mis2, cacheNode, sink, 1.0);
        model.link(routingMatrix);

        return model;
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

        // CLIMB and q-LRU are not product-form: CTMC is the reference solver.
        System.out.println("\n=== Testing cache_replc_climb ===");
        try {
            new jline.solvers.ctmc.CTMC(cache_replc_climb(), "keep", false).getAvgNodeTable().print();
        } catch (Exception e) {
            System.out.println("CTMC Solver failed: " + e.getMessage());
        }

        System.out.println("\n=== Testing cache_replc_qlru ===");
        try {
            new jline.solvers.ctmc.CTMC(cache_replc_qlru(), "keep", false).getAvgNodeTable().print();
        } catch (Exception e) {
            System.out.println("CTMC Solver failed: " + e.getMessage());
        }

        System.out.println("\n=== Testing cache_replc_hlru ===");
        Network model5 = cache_replc_hlru();
        try {
            System.out.println("--- CTMC Solver (exact) ---");
            new jline.solvers.ctmc.CTMC(model5, "keep", false, "cutoff", 1).getAvgNodeTable().print();
            model5.reset();
        } catch (Exception e) {
            System.out.println("CTMC Solver failed: " + e.getMessage());
        }
        try {
            System.out.println("--- MVA Solver (TTL approximation) ---");
            new jline.solvers.mva.MVA(model5).getAvgNodeTable().print();
            model5.reset();
        } catch (Exception e) {
            System.out.println("MVA Solver failed: " + e.getMessage());
        }
        try {
            System.out.println("--- SSA Solver ---");
            new SSA(model5, "samples", 10000, "seed", 23000).getAvgNodeTable().print();
        } catch (Exception e) {
            System.out.println("SSA Solver failed: " + e.getMessage());
        }

        System.out.println("\n=== Testing cache_itemsize_costcap ===");
        Network model6 = cache_itemsize_costcap();
        try {
            NC nc = new NC(model6, "exact");
            nc.getAvgNodeTable().print();
            nc.getAvgCacheTable().print();
            nc.getAvgItemTable().print();
            jline.lang.nodes.Cache cache6 = (jline.lang.nodes.Cache) model6.getNodeByName("Cache");
            System.out.println("Hit Ratio: " + cache6.getHitRatio());
            System.out.println("Mean per-list storage cost: " + cache6.getListCost());
        } catch (Exception e) {
            System.out.println("NC Solver failed: " + e.getMessage());
        }

        System.out.println("\n=== Testing cache_rmf_transient ===");
        try {
            cache_rmf_transient();
        } catch (Exception e) {
            System.out.println("CacheRMF failed: " + e.getMessage());
        }

        System.out.println("\n=== Testing cache_mmap_rr_env ===");
        // Marked MMPP2 (M3A layout D = {D0, D11, D12}); D1 = D11 + D12 diagonal.
        // Phase 1 bursty (rate 4, 90% class1); phase 2 calm (rate 1, 80% class2).
        // Off-diagonal of D0 are the phase-switch rates (both 0.5).
        Matrix D0 = new Matrix(new double[][] {{-4.5, 0.5}, {0.5, -1.5}});
        Matrix D11 = new Matrix(new double[][] {{3.6, 0.0}, {0.0, 0.2}});
        Matrix D12 = new Matrix(new double[][] {{0.4, 0.0}, {0.0, 0.8}});
        jline.util.matrix.MatrixCell cell = new jline.util.matrix.MatrixCell(3);
        cell.set(0, D0);
        cell.set(1, D11);
        cell.set(2, D12);
        jline.lang.processes.MarkedMAP mmap = new jline.lang.processes.MarkedMAP(cell, 2);

        Network trueModel = cache_mmap_rr_env(0, 0, mmap);
        jline.lang.nodes.Cache cacheTrue =
                (jline.lang.nodes.Cache) trueModel.getNodeByName("Cache");
        try {
            System.out.println("--- LDES (simulation of the true MMAP-fed cache) ---");
            new jline.solvers.ldes.LDES(trueModel, "samples", 200000, "seed", 23000)
                    .getAvgNodeTable().print();
            System.out.println("LDES hit ratio: " + cacheTrue.getHitRatio());
            trueModel.reset();
        } catch (Exception e) {
            System.out.println("LDES Solver failed: " + e.getMessage());
        }
        try {
            System.out.println("--- CTMC (exact solution of the true system) ---");
            new jline.solvers.ctmc.CTMC(trueModel, "exact", "keep", false, "cutoff", 1)
                    .getAvgNodeTable().print();
            System.out.println("CTMC hit ratio: " + cacheTrue.getHitRatio());
        } catch (Exception e) {
            System.out.println("CTMC Solver failed: " + e.getMessage());
        }

        // The MMPP2 phase as a random environment modulating phase-conditional
        // Poisson arrivals (D1 diagonal); the environment switches at the MMPP2
        // phase-transition rates, i.e. the -D0 diagonal minus that phase's total.
        jline.lang.Environment env = new jline.lang.Environment("MMPPphase", 2);
        env.addStage(0, "Phase1", "bursty",
                cache_mmap_rr_env(D11.get(0, 0), D12.get(0, 0), null));
        env.addStage(1, "Phase2", "calm",
                cache_mmap_rr_env(D11.get(1, 1), D12.get(1, 1), null));
        env.addTransition(0, 1, new Exp(-D0.get(0, 0) - (D11.get(0, 0) + D12.get(0, 0))));
        env.addTransition(1, 0, new Exp(-D0.get(1, 1) - (D11.get(1, 1) + D12.get(1, 1))));
        env.init();
        System.out.println(env.getStageTable());

        jline.solvers.ln.SolverFactory ctmcFactory = new jline.solvers.ln.SolverFactory() {
            public jline.solvers.NetworkSolver at(Network mdl) {
                return new jline.solvers.ctmc.CTMC(mdl, "exact", "keep", false, "cutoff", 1);
            }
        };
        String[] methods = {"avg", "dec", "blend"};
        for (int k = 0; k < methods.length; k++) {
            try {
                SolverOptions opt = new SolverOptions();
                opt.method = methods[k];
                opt.verbose = jline.VerboseLevel.SILENT;
                opt.iter_max = 100;
                opt.iter_tol = 1e-4;
                jline.solvers.env.ENV envSolver =
                        new jline.solvers.env.ENV(env, ctmcFactory, opt);
                envSolver.getAvg();
                jline.lang.nodes.Cache c0 = (jline.lang.nodes.Cache)
                        env.getEnsemble().get(0).getNodeByName("Cache");
                System.out.println("ENV (" + methods[k] + ") hit ratio: " + c0.getHitRatio());
            } catch (Exception e) {
                System.out.println("ENV (" + methods[k] + ") failed: " + e.getMessage());
            }
        }
    }
}
