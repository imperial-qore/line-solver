/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.DiscreteSampler;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

/**
 * Builders for cache models that use a retrieval system (delayed hits).
 *
 * <p>When a request misses the cache, the missed item is fetched through a
 * retrieval system — a user-defined sub-network of queues. While the retrieval
 * is pending, repeat requests for the same item become "delayed hits". See
 * {@link jline.lang.nodes.Cache#setRetrievalSystem}.</p>
 *
 * @see jline.lang.nodes.Cache
 * @see CacheRetrievalSystemExamples
 */
public class CacheRetrievalSystemModel {

    /**
     * Apply per-(queue,item) service to the retrieval system using the new API:
     * a read-class service is set at each queue (so setRetrievalSystem can inherit it)
     * and each item's rate is then applied via {@link Queue#setItemServiceRate}.
     */
    private static void applyRetrievalService(Cache cacheNode, OpenClass jobClass, OpenClass missClass,
                                              Queue[] queues, Matrix serviceRates) {
        int nQueues = queues.length;
        int n = serviceRates.getNumCols();
        for (int q = 0; q < nQueues; q++) {
            queues[q].setService(jobClass, new Exp(serviceRates.get(q, 0)));
        }
        cacheNode.setRetrievalSystem(jobClass, missClass, queues);
        for (int q = 0; q < nQueues; q++) {
            for (int item = 0; item < n; item++) {
                queues[q].setItemServiceRate(cacheNode, jobClass, item, serviceRates.get(q, item));
            }
        }
    }

    /**
     * Apply per-(queue,item) service and per-item routing matrices via the new API
     * (read-class service + {@link Queue#setItemServiceRate} + the setItem* routing methods).
     * Routing index nQueues (last row/col) is the cache; row = from, col = to.
     */
    private static void applyRetrieval(Cache cacheNode, OpenClass jobClass, OpenClass missClass,
                                       Queue[] queues, Matrix serviceRates, Matrix[] routingMatrices) {
        int nQueues = queues.length;
        int n = serviceRates.getNumCols();
        applyRetrievalService(cacheNode, jobClass, missClass, queues, serviceRates);
        for (int item = 0; item < n; item++) {
            Matrix R = routingMatrices[item];
            for (int a = 0; a < nQueues; a++) {
                cacheNode.setItemRoutingProb(jobClass, item, cacheNode, queues[a], R.get(nQueues, a));
                cacheNode.setItemRoutingProb(jobClass, item, queues[a], cacheNode, R.get(a, nQueues));
                for (int b = 0; b < nQueues; b++) {
                    cacheNode.setItemRoutingProb(jobClass, item, queues[a], queues[b], R.get(a, b));
                }
            }
        }
    }

    /**
     * Cache model with a single retrieval queue. The number of items {@code n} is
     * inferred from {@code serviceRates.length}.
     *
     * @param accessProb   access probabilities for each item
     * @param serviceRates per-item retrieval service rates
     * @param itemLevelCap MATLAB-style string defining per-level cache capacity (e.g. "[1]", "[1,2]")
     * @param sched        scheduling strategy for the retrieval queue
     * @return Network of the simple cache model
     */
    public static Network simple_retrieval_system_model(double[] accessProb, double[] serviceRates,
                                                        String itemLevelCap, SchedStrategy sched) {
        return simple_retrieval_system_model(accessProb, serviceRates, itemLevelCap, sched, ReplacementStrategy.FIFO);
    }

    public static Network simple_retrieval_system_model(double[] accessProb, double[] serviceRates,
                                                        String itemLevelCap, SchedStrategy sched,
                                                        ReplacementStrategy replStrat) {
        Network model = new Network("Simple Model");

        int n = serviceRates.length;
        Matrix capacity = new Matrix(itemLevelCap);

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, capacity, replStrat);
        Queue queue = new Queue(model, "Queue", sched);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(1));

        DiscreteSampler pAccess = new DiscreteSampler(new Matrix(accessProb));

        cacheNode.setRead(jobClass, pAccess);
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        // Read-class service (per-item rates applied as overrides below); routing is
        // cache <-> queue for every item.
        queue.setService(jobClass, new Exp(serviceRates[0]));
        cacheNode.setRetrievalSystem(jobClass, missClass, new Queue[] {queue});
        for (int item = 0; item < n; item++) {
            queue.setItemServiceRate(cacheNode, jobClass, item, serviceRates[item]);
            cacheNode.setItemRoutingProb(jobClass, item, cacheNode, queue, 1.0); // cache -> queue
            cacheNode.setItemRoutingProb(jobClass, item, queue, cacheNode, 1.0);  // queue -> cache
        }

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);
        model.link(routingMatrix);

        return model;
    }

    /**
     * Cache model with a chained retrieval system. On a miss, jobs traverse
     * Queue_1 -&gt; ... -&gt; Queue_nQueues before returning to the cache. The number of
     * items {@code n} is inferred from the columns of {@code serviceRates}.
     *
     * @param accessProb   access probabilities for each item
     * @param serviceRates per-queue per-item retrieval service rates [nQueues x nItems]
     * @param nQueues      number of queues in the retrieval chain
     * @param itemLevelCap MATLAB-style string defining per-level cache capacity (e.g. "[1]", "[1,2]")
     * @param sched        scheduling strategy for the retrieval queues
     * @return Network of the chain retrieval system model
     */
    public static Network chain_retrieval_system_model(double[] accessProb, Matrix serviceRates, int nQueues,
                                                       String itemLevelCap, SchedStrategy sched) {
        Network model = new Network("Chain Model");

        int n = serviceRates.getNumCols();
        Matrix capacity = new Matrix(itemLevelCap);

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, capacity, ReplacementStrategy.FIFO);
        Queue[] queues = new Queue[nQueues];
        for (int i = 0; i < nQueues; i++) {
            queues[i] = new Queue(model, "Queue_" + (i + 1), sched);
        }
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(1));

        DiscreteSampler pAccess = new DiscreteSampler(new Matrix(accessProb));

        cacheNode.setRead(jobClass, pAccess);
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        applyRetrievalService(cacheNode, jobClass, missClass, queues, serviceRates);

        for (int item = 0; item < n; item++) {
            cacheNode.setItemRoutingProb(jobClass, item, cacheNode, queues[0], 1.0);
            for (int queueIdx = 0; queueIdx < nQueues - 1; queueIdx++) {
                cacheNode.setItemRoutingProb(jobClass, item, queues[queueIdx], queues[queueIdx + 1], 1.0);
            }
            cacheNode.setItemRoutingProb(jobClass, item, queues[nQueues - 1], cacheNode, 1.0);
        }

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);
        model.link(routingMatrix);

        return model;
    }

    /**
     * Cache model with probabilistic per-item routing through a retrieval system of 3 queues
     * (one infinite-server, two scheduled). The number of items and queues is inferred from
     * the dimensions of {@code serviceRates}.
     *
     * @param accessProb      access probabilities for each item
     * @param serviceRates    per-queue per-item retrieval service rates [nQueues x nItems]
     * @param routingMatrices per-item routing matrices, one per item [(nQueues+1) x (nQueues+1)]
     * @param itemLevelCap    MATLAB-style string defining per-level cache capacity (e.g. "[2]", "[1,2]")
     * @param sched           scheduling strategy for the non-IS retrieval queues
     * @return Network of the probabilistic routing retrieval system model
     */
    public static Network retrieval_system_with_probabilistic_routing(double[] accessProb, Matrix serviceRates,
                                                                      Matrix[] routingMatrices, String itemLevelCap,
                                                                      SchedStrategy sched) {
        Network model = new Network("Probabilistic Routing");

        int n = serviceRates.getNumCols();
        Matrix capacity = new Matrix(itemLevelCap);

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, capacity, ReplacementStrategy.FIFO);
        Queue is_queue = new Queue(model, "IS Queue", SchedStrategy.INF);
        Queue queue_1 = new Queue(model, "Queue 1", sched);
        Queue queue_2 = new Queue(model, "Queue 2", sched);
        Queue[] queues = new Queue[] {is_queue, queue_1, queue_2};
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(1));

        DiscreteSampler pAccess = new DiscreteSampler(new Matrix(accessProb));

        cacheNode.setRead(jobClass, pAccess);
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        applyRetrieval(cacheNode, jobClass, missClass, queues, serviceRates, routingMatrices);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);
        model.link(routingMatrix);

        return model;
    }

    /** Overload using the default 3-item routing matrices defined for the probabilistic routing model. */
    public static Network retrieval_system_with_probabilistic_routing(double[] accessProb, Matrix serviceRates,
                                                                      String itemLevelCap, SchedStrategy sched) {
        Matrix[] routingMatrices = new Matrix[] {
            new Matrix(new double[][] {
                {0.00, 0.50, 0.00, 0.50}, // from IS(0)
                {0.00, 0.00, 0.70, 0.30}, // from Queue1(1)
                {0.00, 0.00, 0.00, 1.00}, // from Queue2(2)
                {0.70, 0.30, 0.00, 0.00}, // from Cache(3)
            }),
            new Matrix(new double[][] {
                {0.00, 0.30, 0.00, 0.70}, // from IS(0)
                {0.00, 0.00, 0.50, 0.50}, // from Queue1(1)
                {0.00, 0.00, 0.00, 1.00}, // from Queue2(2)
                {0.20, 0.80, 0.00, 0.00}, // from Cache(3)
            }),
            new Matrix(new double[][] {
                {0.00, 0.60, 0.00, 0.40}, // from IS(0)
                {0.00, 0.00, 0.40, 0.60}, // from Queue1(1)
                {0.00, 0.00, 0.00, 1.00}, // from Queue2(2)
                {0.50, 0.50, 0.00, 0.00}, // from Cache(3)
            }),
        };
        return retrieval_system_with_probabilistic_routing(accessProb, serviceRates, routingMatrices, itemLevelCap, sched);
    }

    /**
     * Cache model with probabilistic per-item routing including self-loops through a retrieval
     * system of 3 queues (one infinite-server, two scheduled). The number of items and queues is
     * inferred from the dimensions of {@code serviceRates}.
     *
     * @param accessProb      access probabilities for each item
     * @param serviceRates    per-queue per-item retrieval service rates [nQueues x nItems]
     * @param routingMatrices per-item routing matrices, one per item [(nQueues+1) x (nQueues+1)]
     * @param itemLevelCap    MATLAB-style string defining per-level cache capacity (e.g. "[2]", "[1,2]")
     * @param sched           scheduling strategy for the non-IS retrieval queues
     * @return Network of the probabilistic routing with self-loops retrieval system model
     */
    public static Network retrieval_system_with_probabilistic_routing_and_self_loops(double[] accessProb,
                                                                                    Matrix serviceRates,
                                                                                    Matrix[] routingMatrices,
                                                                                    String itemLevelCap,
                                                                                    SchedStrategy sched) {
        Network model = new Network("Probabilistic Routing Self Loops");

        int n = serviceRates.getNumCols();
        Matrix capacity = new Matrix(itemLevelCap);

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, capacity, ReplacementStrategy.FIFO);
        Queue is_queue = new Queue(model, "IS Queue", SchedStrategy.INF);
        Queue queue_1 = new Queue(model, "Queue 1", sched);
        Queue queue_2 = new Queue(model, "Queue 2", sched);
        Queue[] queues = new Queue[] {is_queue, queue_1, queue_2};
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(1));

        DiscreteSampler pAccess = new DiscreteSampler(new Matrix(accessProb));

        cacheNode.setRead(jobClass, pAccess);
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        applyRetrieval(cacheNode, jobClass, missClass, queues, serviceRates, routingMatrices);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);
        model.link(routingMatrix);

        return model;
    }

    /** Overload using the default 3-item routing matrices (with self-loops) defined for the probabilistic routing model. */
    public static Network retrieval_system_with_probabilistic_routing_and_self_loops(double[] accessProb,
                                                                                    Matrix serviceRates,
                                                                                    String itemLevelCap,
                                                                                    SchedStrategy sched) {
        Matrix[] routingMatrices = new Matrix[] {
            new Matrix(new double[][] {
                {0.30, 0.50, 0.00, 0.20}, // from IS(0):     loop=0.30, ->Queue1=0.50, ->Cache=0.20
                {0.00, 0.20, 0.60, 0.20}, // from Queue1(1): loop=0.20, ->Queue2=0.60, ->Cache=0.20
                {0.00, 0.00, 0.10, 0.90}, // from Queue2(2): loop=0.10, ->Cache=0.90
                {0.80, 0.20, 0.00, 0.00}, // from Cache(3):  ->IS=0.80, ->Queue1=0.20
            }),
            new Matrix(new double[][] {
                {0.10, 0.40, 0.00, 0.50}, // from IS(0):     loop=0.10, ->Queue1=0.40, ->Cache=0.50
                {0.00, 0.40, 0.30, 0.30}, // from Queue1(1): loop=0.40, ->Queue2=0.30, ->Cache=0.30
                {0.00, 0.00, 0.05, 0.95}, // from Queue2(2): loop=0.05, ->Cache=0.95
                {0.25, 0.75, 0.00, 0.00}, // from Cache(3):  ->IS=0.25, ->Queue1=0.75
            }),
            new Matrix(new double[][] {
                {0.50, 0.30, 0.00, 0.20}, // from IS(0):     loop=0.50, ->Queue1=0.30, ->Cache=0.20
                {0.00, 0.10, 0.50, 0.40}, // from Queue1(1): loop=0.10, ->Queue2=0.50, ->Cache=0.40
                {0.00, 0.00, 0.20, 0.80}, // from Queue2(2): loop=0.20, ->Cache=0.80
                {0.50, 0.50, 0.00, 0.00}, // from Cache(3):  ->IS=0.50, ->Queue1=0.50
            }),
        };
        return retrieval_system_with_probabilistic_routing_and_self_loops(accessProb, serviceRates, routingMatrices,
                itemLevelCap, sched);
    }

    /**
     * Retrieval system whose per-item miss routing AND service are taken, by
     * default, from the read class.
     *
     * <p>Routing comes from the read class's edges among the retrieval queues in
     * the top-level routing matrix P, and service from the read class's service
     * distribution at each queue. Per-item overrides
     * ({@link Cache#setItemRoutingProb} with the cache as source or destination,
     * and {@link Queue#setItemServiceRate}) then reconfigure one item at finer
     * granularity: here item 0 skips Queue_2 and is fetched faster at
     * Queue_1.</p>
     *
     * <p>The retrieval stations are PS: per-item (class-dependent) service rates
     * are admissible in the analytical retrieval algorithm, whereas FCFS/SIRO
     * would require identical rates.</p>
     *
     * @return configured delayed-hit cache network model
     */
    public static Network retrieval_default() {
        double[] accessProb = {0.6, 0.3, 0.1};   // per-item access probabilities

        Network model = new Network("DelayedHits");

        int n = accessProb.length;               // number of items
        Matrix capacity = new Matrix("[1]");     // per-level cache capacity

        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, capacity, ReplacementStrategy.FIFO);
        Queue queue1 = new Queue(model, "Queue_1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue_2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(1));

        // Read-class service at each retrieval queue = default per-item fetch service.
        queue1.setService(jobClass, new Exp(2.0));
        queue2.setService(jobClass, new Exp(3.0));

        cacheNode.setRead(jobClass, new DiscreteSampler(new Matrix(accessProb)));
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        // No service rates and no routing matrices: both inherited from the read class.
        cacheNode.setRetrievalSystem(jobClass, missClass, new Queue[] {queue1, queue2});

        // Item-level overrides for item 0: skip Queue_2 and fetch faster at Queue_1.
        cacheNode.setItemRoutingProb(jobClass, 0, queue1, queue2, 0.0);      // delete default edge
        cacheNode.setItemRoutingProb(jobClass, 0, queue1, cacheNode, 1.0);   // exit after Queue_1
        queue1.setItemServiceRate(cacheNode, jobClass, 0, 5.0);              // faster item-0 fetch

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        // Default retrieval topology, drawn once for the read class: cache -> Q1 -> Q2 -> cache
        routingMatrix.set(jobClass, jobClass, cacheNode, queue1, 1.0);   // entry into retrieval
        routingMatrix.set(jobClass, jobClass, queue1, queue2, 1.0);      // Queue_1 -> Queue_2
        routingMatrix.set(jobClass, jobClass, queue2, cacheNode, 1.0);   // exit back to cache
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);
        model.link(routingMatrix);

        return model;
    }
}
