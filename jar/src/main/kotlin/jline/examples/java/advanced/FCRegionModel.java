/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.advanced;

import jline.lang.Region;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;

import java.util.Arrays;
import java.util.Collections;

/**
 * Examples of models with Finite Capacity Regions (FCR).
 * FCRs allow defining capacity constraints across multiple nodes,
 * with either blocking (wait) or dropping behavior when full.
 */
public class FCRegionModel {

    /**
     * Multiclass open network with FCR blocking.
     * <p>
     * Features:
     * - Two open classes with different arrival/service rates
     * - Two queues covered by a single FCR
     * - Global max: 8 jobs, Class1 max: 5, Class2 max: 4
     * - Blocking behavior: jobs wait when region is full
     *
     * @return configured FCR blocking model
     */
    public static Network fcr_oqnwaitq() {
        Network model = new Network("FCR Blocking Example");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass class1 = new OpenClass(model, "Class1", 0);
        OpenClass class2 = new OpenClass(model, "Class2", 1);

        source.setArrival(class1, Exp.fitRate(0.4));
        source.setArrival(class2, Exp.fitRate(0.3));
        queue1.setService(class1, Exp.fitRate(1.0));
        queue1.setService(class2, Exp.fitRate(0.9));
        queue2.setService(class1, Exp.fitRate(1.1));
        queue2.setService(class2, Exp.fitRate(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        // Class1 routing
        P.set(class1, class1, source, queue1, 0.5);
        P.set(class1, class1, source, queue2, 0.5);
        P.set(class1, class1, queue1, queue2, 0.3);
        P.set(class1, class1, queue1, sink, 0.7);
        P.set(class1, class1, queue2, sink, 1.0);

        // Class2 routing
        P.set(class2, class2, source, queue1, 0.6);
        P.set(class2, class2, source, queue2, 0.4);
        P.set(class2, class2, queue1, queue2, 0.5);
        P.set(class2, class2, queue1, sink, 0.5);
        P.set(class2, class2, queue2, sink, 1.0);

        model.link(P);

        // Add FCR with blocking
        model.addRegion(Arrays.asList(queue1, queue2));
        Region fcr = model.getRegions().get(0);
        fcr.setGlobalMaxJobs(8);
        fcr.setClassMaxJobs(class1, 5);
        fcr.setClassMaxJobs(class2, 4);
        fcr.setDropRule(class1, false);  // false = block
        fcr.setDropRule(class2, false);

        return model;
    }

    /**
     * Multiclass open network with FCR dropping.
     * <p>
     * Features:
     * - Two open classes with different arrival/service rates
     * - Two queues covered by a single FCR
     * - Global max: 8 jobs, Class1 max: 5, Class2 max: 4
     * - Dropping behavior: jobs are lost when region is full
     *
     * @return configured FCR dropping model
     */
    public static Network fcr_oqndrop() {
        Network model = new Network("FCR Dropping Example");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass class1 = new OpenClass(model, "Class1", 0);
        OpenClass class2 = new OpenClass(model, "Class2", 1);

        source.setArrival(class1, Exp.fitRate(0.4));
        source.setArrival(class2, Exp.fitRate(0.3));
        queue1.setService(class1, Exp.fitRate(1.0));
        queue1.setService(class2, Exp.fitRate(0.9));
        queue2.setService(class1, Exp.fitRate(1.1));
        queue2.setService(class2, Exp.fitRate(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        // Class1 routing
        P.set(class1, class1, source, queue1, 0.5);
        P.set(class1, class1, source, queue2, 0.5);
        P.set(class1, class1, queue1, queue2, 0.3);
        P.set(class1, class1, queue1, sink, 0.7);
        P.set(class1, class1, queue2, sink, 1.0);

        // Class2 routing
        P.set(class2, class2, source, queue1, 0.6);
        P.set(class2, class2, source, queue2, 0.4);
        P.set(class2, class2, queue1, queue2, 0.5);
        P.set(class2, class2, queue1, sink, 0.5);
        P.set(class2, class2, queue2, sink, 1.0);

        model.link(P);

        // Add FCR with dropping
        model.addRegion(Arrays.asList(queue1, queue2));
        Region fcr = model.getRegions().get(0);
        fcr.setGlobalMaxJobs(8);
        fcr.setClassMaxJobs(class1, 5);
        fcr.setClassMaxJobs(class2, 4);
        fcr.setDropRule(class1, true);  // true = drop
        fcr.setDropRule(class2, true);

        return model;
    }

    /**
     * Simple M/M/1 with FCR blocking.
     * <p>
     * This model demonstrates that FCR with blocking around a single
     * queue behaves identically to a standard M/M/1 queue.
     *
     * @return configured M/M/1 with FCR blocking
     */
    public static Network fcr_mm1waitq() {
        Network model = new Network("FCR MM1 Blocking");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass = new OpenClass(model, "Class1", 0);

        source.setArrival(jobclass, Exp.fitRate(0.5));
        queue.setService(jobclass, Exp.fitRate(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, source, queue, 1.0);
        P.set(jobclass, jobclass, queue, sink, 1.0);
        model.link(P);

        // Add FCR with blocking
        model.addRegion(Arrays.asList(queue));
        Region fcr = model.getRegions().get(0);
        fcr.setGlobalMaxJobs(10);
        fcr.setDropRule(jobclass, false);  // false = block

        return model;
    }

    /**
     * Simple M/M/1 without FCR (for comparison).
     *
     * @return standard M/M/1 model
     */
    public static Network mm1() {
        Network model = new Network("MM1");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass = new OpenClass(model, "Class1", 0);

        source.setArrival(jobclass, Exp.fitRate(0.5));
        queue.setService(jobclass, Exp.fitRate(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, source, queue, 1.0);
        P.set(jobclass, jobclass, queue, sink, 1.0);
        model.link(P);

        return model;
    }

    /**
     * M/M/1/K with FCR dropping (K=2).
     * <p>
     * This model demonstrates that FCR with dropping around a single
     * queue behaves like an M/M/1/K queue where K is the FCR capacity.
     *
     * @return configured M/M/1/K with FCR dropping
     */
    public static Network fcr_mm1kdrop() {
        Network model = new Network("FCR MM1K Dropping");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass = new OpenClass(model, "Class1", 0);

        source.setArrival(jobclass, Exp.fitRate(0.8));
        queue.setService(jobclass, Exp.fitRate(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, source, queue, 1.0);
        P.set(jobclass, jobclass, queue, sink, 1.0);
        model.link(P);

        // Add FCR with dropping (K=2)
        model.addRegion(Arrays.asList(queue));
        Region fcr = model.getRegions().get(0);
        fcr.setGlobalMaxJobs(2);
        fcr.setDropRule(jobclass, true);  // true = drop

        return model;
    }

    /**
     * M/M/1/K using queue capacity (K=2, for comparison).
     *
     * @return standard M/M/1/K model
     */
    public static Network mm1k() {
        Network model = new Network("MM1K");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(1);
        queue.setCapacity(2);  // K=2
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass = new OpenClass(model, "Class1", 0);

        source.setArrival(jobclass, Exp.fitRate(0.8));
        queue.setService(jobclass, Exp.fitRate(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, source, queue, 1.0);
        P.set(jobclass, jobclass, queue, sink, 1.0);
        model.link(P);

        return model;
    }

    /**
     * FCR with multiple constraint types demonstration.
     * <p>
     * Features:
     * - Global max jobs: limits total jobs across all classes
     * - Per-class max jobs: limits jobs of a specific class
     * - Different drop rules per class: blocking vs dropping
     *
     * @return configured FCR model with multiple constraint types
     */
    public static Network fcr_constraints() {
        Network model = new Network("FCR Constraints Demo");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Two classes with different priorities
        OpenClass highPriority = new OpenClass(model, "HighPriority", 0);
        OpenClass lowPriority = new OpenClass(model, "LowPriority", 1);

        // Arrival and service rates
        source.setArrival(highPriority, Exp.fitRate(0.3));
        source.setArrival(lowPriority, Exp.fitRate(0.5));
        queue1.setService(highPriority, Exp.fitRate(1.0));
        queue1.setService(lowPriority, Exp.fitRate(0.8));
        queue2.setService(highPriority, Exp.fitRate(1.2));
        queue2.setService(lowPriority, Exp.fitRate(1.0));

        // Routing: both classes go through both queues
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(highPriority, highPriority, source, queue1, 1.0);
        P.set(highPriority, highPriority, queue1, queue2, 1.0);
        P.set(highPriority, highPriority, queue2, sink, 1.0);
        P.set(lowPriority, lowPriority, source, queue1, 1.0);
        P.set(lowPriority, lowPriority, queue1, queue2, 1.0);
        P.set(lowPriority, lowPriority, queue2, sink, 1.0);
        model.link(P);

        // Add FCR with multiple constraint types
        model.addRegion(Arrays.asList(queue1, queue2));
        Region fcr = model.getRegions().get(0);

        // Global constraint: max 6 jobs total in the region
        fcr.setGlobalMaxJobs(6);

        // Per-class constraints: high priority gets more space
        fcr.setClassMaxJobs(highPriority, 4);  // max 4 high priority jobs
        fcr.setClassMaxJobs(lowPriority, 3);   // max 3 low priority jobs

        // Drop rules: high priority jobs wait, low priority jobs are dropped
        fcr.setDropRule(highPriority, false);  // block (wait)
        fcr.setDropRule(lowPriority, true);    // drop

        return model;
    }

    /**
     * Loss network with FCR (for NC solver lossn method).
     * <p>
     * This model demonstrates the NC solver's ability to analyze open
     * loss networks using the Erlang fixed-point approximation.
     * Features:
     * - Single Delay node (infinite server) inside an FCR
     * - Multiple classes with different arrival/service rates
     * - DROP policy: jobs are lost when region is full
     * <p>
     * The NC solver can analytically solve this using the lossn method.
     *
     * @return configured loss network model
     */
    public static Network fcr_lossn() {
        Network model = new Network("FCR Loss Network");

        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Sink sink = new Sink(model, "Sink");

        OpenClass class1 = new OpenClass(model, "Class1", 0);
        OpenClass class2 = new OpenClass(model, "Class2", 1);

        // Arrival and service rates
        double lambda1 = 0.3, lambda2 = 0.2;  // arrival rates
        double mu1 = 1.0, mu2 = 0.8;          // service rates

        source.setArrival(class1, Exp.fitRate(lambda1));
        source.setArrival(class2, Exp.fitRate(lambda2));
        delay.setService(class1, Exp.fitRate(mu1));
        delay.setService(class2, Exp.fitRate(mu2));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, source, delay, 1.0);
        P.set(class1, class1, delay, sink, 1.0);
        P.set(class2, class2, source, delay, 1.0);
        P.set(class2, class2, delay, sink, 1.0);
        model.link(P);

        // Add FCR with dropping
        model.addRegion(Collections.singletonList(delay));
        Region fcr = model.getRegions().get(0);
        fcr.setGlobalMaxJobs(5);           // Global: max 5 jobs
        fcr.setClassMaxJobs(class1, 3);    // Class1: max 3 jobs
        fcr.setClassMaxJobs(class2, 3);    // Class2: max 3 jobs
        fcr.setDropRule(class1, true);     // true = drop
        fcr.setDropRule(class2, true);

        return model;
    }
}
