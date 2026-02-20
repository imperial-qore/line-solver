/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.solvers.SolverOptions;
import jline.solvers.jmt.JMTOptions;
import jline.solvers.ssa.SSA;
import org.jetbrains.annotations.NotNull;

/**
 * Examples of queueing models with priorities
 */
public class PrioModel {

    /**
     * Open network with multiple scheduling strategies and priority classes.
     * <p>
     * Features:
     * - Three open classes with different priority levels (Class2 has priority 1)
     * - Five stations with different scheduling: FCFS, SIRO, PS, HOL
     * - HOL (Head-of-Line) queue demonstrates priority scheduling
     * - Complex routing with 25% probability splits from central queue
     * - All classes have identical arrival rates but different service rates
     * - Demonstrates impact of scheduling strategies on class performance
     *
     * @return configured priority-based open network model
     */
    public static Network prio_hol_open() {
        Network model = new Network("MyNetwork");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "WebServer", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Storage1", SchedStrategy.SIRO);
        Queue queue3 = new Queue(model, "Storage2", SchedStrategy.PS);
        Queue queue4 = new Queue(model, "Storage3", SchedStrategy.HOL);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 1);
        OpenClass jobclass3 = new OpenClass(model, "Class3", 0);

        source.setArrival(jobclass1, Exp.fitMean(10.00)); // (Source,Class1)
        source.setArrival(jobclass2, Exp.fitMean(10.00)); // (Source,Class2)
        source.setArrival(jobclass3, Exp.fitMean(10.00)); // (Source,Class3)
        queue1.setService(jobclass1, Exp.fitMean(0.300)); // (Queue1,Class1)
        queue1.setService(jobclass2, Exp.fitMean(0.500)); // (Queue1,Class2)
        queue1.setService(jobclass3, Exp.fitMean(0.600)); // (Queue1,Class3)
        queue2.setService(jobclass1, Exp.fitMean(1.100)); // (Queue2,Class1)
        queue2.setService(jobclass2, Exp.fitMean(1.300)); // (Queue2,Class2)
        queue2.setService(jobclass3, Exp.fitMean(1.500)); // (Queue2,Class3)
        queue3.setService(jobclass1, Exp.fitMean(2.00)); // (Queue3,Class1)
        queue3.setService(jobclass2, Exp.fitMean(2.100)); // (Queue3,Class2)
        queue3.setService(jobclass3, Exp.fitMean(1.900)); // (Queue3,Class3)
        queue4.setService(jobclass1, Exp.fitMean(2.500)); // (Queue4,Class1)
        queue4.setService(jobclass2, Exp.fitMean(1.900)); // (Queue4,Class2)
        queue4.setService(jobclass3, Exp.fitMean(4.300)); // (Queue4,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, source, queue1, 1.00); // (Source,Class1) -> (Queue1,Class1)

        routingMatrix.set(jobclass1, jobclass1, queue1, queue2, 0.250000); // (Queue1,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, queue1, queue3, 0.250000); // (Queue1,Class1) -> (Queue3,Class1)
        routingMatrix.set(jobclass1, jobclass1, queue1, queue4, 0.250000); // (Queue1,Class1) -> (Queue4,Class1)
        routingMatrix.set(jobclass1, jobclass1, queue1, sink, 0.250000); // (Queue1,Class1) -> (Sink,Class1)

        routingMatrix.set(jobclass1, jobclass1, queue2, queue1, 1.00); // (Queue2,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, queue3, queue1, 1.00); // (Queue3,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, queue4, queue1, 1.00); // (Queue4,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass2, jobclass2, source, queue1, 1.00); // (Source,Class2) -> (Queue1,Class2)

        routingMatrix.set(jobclass2, jobclass2, queue1, queue2, 0.250000); // (Queue1,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, queue1, queue3, 0.250000); // (Queue1,Class2) -> (Queue3,Class2)
        routingMatrix.set(jobclass2, jobclass2, queue1, queue4, 0.250000); // (Queue1,Class2) -> (Queue4,Class2)
        routingMatrix.set(jobclass2, jobclass2, queue1, sink, 0.250000); // (Queue1,Class2) -> (Sink,Class2)

        routingMatrix.set(jobclass2, jobclass2, queue2, queue1, 1.00); // (Queue2,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, queue3, queue1, 1.00); // (Queue3,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, queue4, queue1, 1.00); // (Queue4,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass3, jobclass3, source, queue1, 1.00); // (Source,Class3) -> (Queue1,Class3)

        routingMatrix.set(jobclass3, jobclass3, queue1, queue2, 0.250000); // (Queue1,Class3) -> (Queue2,Class3)
        routingMatrix.set(jobclass3, jobclass3, queue1, queue3, 0.250000); // (Queue1,Class3) -> (Queue3,Class3)
        routingMatrix.set(jobclass3, jobclass3, queue1, queue4, 0.250000); // (Queue1,Class3) -> (Queue4,Class3)
        routingMatrix.set(jobclass3, jobclass3, queue1, sink, 0.250000); // (Queue1,Class3) -> (Sink,Class3)

        routingMatrix.set(jobclass3, jobclass3, queue2, queue1, 1.00); // (Queue2,Class3) -> (Queue1,Class3)
        routingMatrix.set(jobclass3, jobclass3, queue3, queue1, 1.00); // (Queue3,Class3) -> (Queue1,Class3)
        routingMatrix.set(jobclass3, jobclass3, queue4, queue1, 1.00); // (Queue4,Class3) -> (Queue1,Class3)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Closed network comparing scheduling strategies with priorities.
     * <p>
     * Features:
     * - Three closed classes (18 jobs each) with Class2 having priority 1
     * - Six stations: SlowDelay, FCFS, SIRO, PS, HOL queues, FastDelay
     * - Central FCFS queue distributes traffic (25% each direction)
     * - HOL queue provides priority service for Class2
     * - Cyclic topology demonstrating sustained priority effects
     * - Identical service rates at delay stations for fair comparison
     *
     * @return configured closed priority network model
     */
    public static Network prio_hol_closed() {
        Network model = new Network("MyNetwork");

        // Block 1: nodes
        Delay node1 = new Delay(model, "SlowDelay");
        Queue node2 = new Queue(model, "FCFSQueue", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "SIROQueue", SchedStrategy.SIRO);
        Queue node4 = new Queue(model, "PSQueue", SchedStrategy.PS);
        Queue node5 = new Queue(model, "HOLQueue", SchedStrategy.HOL);
        Delay node6 = new Delay(model, "FastDelay");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 18, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 18, node1, 1);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 18, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(10.00)); // (SlowDelay,Class1)
        node1.setService(jobclass2, Exp.fitMean(10.00)); // (SlowDelay,Class2)
        node1.setService(jobclass3, Exp.fitMean(10.00)); // (SlowDelay,Class3)

        node2.setService(jobclass1, Exp.fitMean(0.300)); // (FCFSQueue,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.500)); // (FCFSQueue,Class2)
        node2.setService(jobclass3, Exp.fitMean(0.600)); // (FCFSQueue,Class3)

        node3.setService(jobclass1, Exp.fitMean(1.100)); // (SIROQueue,Class1)
        node3.setService(jobclass2, Exp.fitMean(1.300)); // (SIROQueue,Class2)
        node3.setService(jobclass3, Exp.fitMean(1.500)); // (SIROQueue,Class3)

        node4.setService(jobclass1, Exp.fitMean(1.00)); // (PSQueue,Class1)
        node4.setService(jobclass2, Exp.fitMean(1.100)); // (PSQueue,Class2)
        node4.setService(jobclass3, Exp.fitMean(1.900)); // (PSQueue,Class3)

        node5.setService(jobclass1, Exp.fitMean(2.500)); // (HOLQueue,Class1)
        node5.setService(jobclass2, Exp.fitMean(1.900)); // (HOLQueue,Class2)
        node5.setService(jobclass3, Exp.fitMean(4.300)); // (HOLQueue,Class3)

        node6.setService(jobclass1, Exp.fitMean(1.00)); // (FastDelay,Class1)
        node6.setService(jobclass2, Exp.fitMean(1.00)); // (FastDelay,Class2)
        node6.setService(jobclass3, Exp.fitMean(1.00)); // (FastDelay,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (SlowDelay,Class1) -> (FCFSQueue,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 0.250000); // (FCFSQueue,Class1) -> (SIROQueue,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node4, 0.250000); // (FCFSQueue,Class1) -> (PSQueue,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node5, 0.250000); // (FCFSQueue,Class1) -> (HOLQueue,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node6, 0.250000); // (FCFSQueue,Class1) -> (FastDelay,Class1)

        routingMatrix.set(jobclass1, jobclass1, node3, node2, 1.00); // (SIROQueue,Class1) -> (FCFSQueue,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node2, 1.00); // (PSQueue,Class1) -> (FCFSQueue,Class1)
        routingMatrix.set(jobclass1, jobclass1, node5, node2, 1.00); // (HOLQueue,Class1) -> (FCFSQueue,Class1)

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (SlowDelay,Class2) -> (FCFSQueue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 0.250000); // (FCFSQueue,Class2) -> (SIROQueue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node4, 0.250000); // (FCFSQueue,Class2) -> (PSQueue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node5, 0.250000); // (FCFSQueue,Class2) -> (HOLQueue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node6, 0.250000); // (FCFSQueue,Class2) -> (FastDelay,Class2)

        routingMatrix.set(jobclass2, jobclass2, node3, node2, 1.00); // (SIROQueue,Class2) -> (FCFSQueue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node2, 1.00); // (PSQueue,Class2) -> (FCFSQueue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node5, node2, 1.00); // (HOLQueue,Class2) -> (FCFSQueue,Class2)

        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.00); // (SlowDelay,Class3) -> (FCFSQueue,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node3, 0.250000); // (FCFSQueue,Class3) -> (SIROQueue,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node4, 0.250000); // (FCFSQueue,Class3) -> (PSQueue,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node5, 0.250000); // (FCFSQueue,Class3) -> (HOLQueue,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node6, 0.250000); // (FCFSQueue,Class3) -> (FastDelay,Class3)

        routingMatrix.set(jobclass3, jobclass3, node3, node2, 1.00); // (SIROQueue,Class3) -> (FCFSQueue,Class3)
        routingMatrix.set(jobclass3, jobclass3, node4, node2, 1.00); // (PSQueue,Class3) -> (FCFSQueue,Class3)
        routingMatrix.set(jobclass3, jobclass3, node5, node2, 1.00); // (HOLQueue,Class3) -> (FCFSQueue,Class3)

        routingMatrix.set(jobclass1, jobclass1, node6, node1, 1.00); // (FastDelay,Class1) -> (SlowDelay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node6, node1, 1.00); // (FastDelay,Class2) -> (SlowDelay,Class2)
        routingMatrix.set(jobclass3, jobclass3, node6, node1, 1.00); // (FastDelay,Class3) -> (SlowDelay,Class3)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Simple priority network with PSPRIO scheduling.
     * <p>
     * Features:
     * - Two closed classes (2 jobs each) with Class2 having priority 1
     * - PSPRIO (Processor Sharing with Priorities) scheduling
     * - Mixed service distributions: Erlang for Class1, HyperExp for Class2
     * - Serial two-station topology: Delay → PSPRIO Queue → Delay
     * - Demonstrates preemptive priority in processor sharing
     *
     * @return configured PSPRIO network model
     */
    @NotNull
    public static Network prio_psprio() {
        Network model = new Network("MyNetwork");

        // Block 1: nodes
        Delay node1 = new Delay(model, "SlowDelay");
        Queue node2 = new Queue(model, "PSPRIOQueue", SchedStrategy.PSPRIO);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 1);

        node1.setService(jobclass1, new Erlang(3.0, 2)); // (SlowDelay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (SlowDelay,Class2)

        node2.setService(jobclass1, new HyperExp(0.1, 1.0, 10.0)); // (FCFSQueue,Class1)
        node2.setService(jobclass2, new Exp(1)); // (FCFSQueue,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass1, Network.serialRouting(node1, node2));
        routingMatrix.set(jobclass2, Network.serialRouting(node1, node2));
        model.link(routingMatrix);
        return model;
    }

    /**
     * Network with GPSPRIO scheduling and multiple priority classes.
     * <p>
     * Features:
     * - Four closed classes with different populations and priorities
     * - Two stations: Delay node and GPSPRIO queue with weights
     * - Class1: 6 jobs, priority 0; Class2,3: 4 jobs, priority 1; Class4: 1 job, priority 2
     * - Demonstrates GPS (Generalized Processor Sharing) with priorities
     * - Different service rates and weights for weighted fair sharing
     * - Simple serial routing: Delay → GPSPRIO Queue → Delay
     *
     * @return configured network with GPSPRIO scheduling
     */
    public static Network prio_identical() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.GPSPRIO);

        // Block 2: classes with different populations and priorities
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 6, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 4, node1, 1);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 4, node1, 1);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 1, node1, 2);

        // Service rates at delay station
        node1.setService(jobclass1, new Erlang(3, 2));
        node1.setService(jobclass2, new Exp(1));
        node1.setService(jobclass3, new Exp(1.0));
        node1.setService(jobclass4, new Exp(2.0));

        // GPSPRIO queue with weights and service rates
        int w1 = 12; node2.setService(jobclass1, new Exp(30), w1);
        int w2 = 3;  node2.setService(jobclass2, new Exp(2), w2);
        int w3 = 5;  node2.setService(jobclass3, new Exp(12), w3);
        int w4 = 1;  node2.setService(jobclass4, new Exp(1), w4);

        // Block 3: topology - simple serial routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        
        routingMatrix.set(jobclass1, Network.serialRouting(node1, node2));
        routingMatrix.set(jobclass2, Network.serialRouting(node1, node2));
        routingMatrix.set(jobclass3, Network.serialRouting(node1, node2));
        routingMatrix.set(jobclass4, Network.serialRouting(node1, node2));

        model.link(routingMatrix);
        return model;
    }

    /**
     * Main method for testing and demonstrating priority examples.
     *
     * <p>Currently configured to:
     * - Run prio_psprio() with PSPRIO scheduling
     * - Use JMT options with specified seed and sample count
     * - Solve using SSA solver and print average metrics
     * - Model visualization and JMT solving are commented out
     *
     * @param args command line arguments (not used)
     * @throws Exception if solver encounters an error
     */
    public static void main(String[] args) throws Exception {
        Network model = prio_psprio();
        SolverOptions options = new JMTOptions();
        options.seed = 23000;
        options.samples = 10000;
        //model.view();
        //new JMT(model, options).getAvgTable().print();
        new SSA(model, options).getAvgTable().print();
    }

}
