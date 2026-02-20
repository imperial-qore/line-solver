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
import jline.lang.nodes.*;
import jline.lang.processes.Disabled;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.solvers.jmt.JMT;
import jline.solvers.mva.MVA;

/**
 * Examples of fork-join queueing networks
 */
public class ForkJoinModel {

    /**
     * Basic open fork-join network with single class.
     * <p>
     * <p>
     * Features:
     * - Single open class with exponential arrivals (rate 0.05)
     * - Fork node splits jobs to two parallel FCFS queues
     * - Join node synchronizes completion from both queues
     * - Different service rates: Queue1 (rate 1.0), Queue2 (rate 2.0)
     * - Demonstrates basic fork-join synchronization
     *
     * @return configured fork-join network model
     */
    public static Network fj_basic_open() {
        Network model = new Network("model");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass1 = new OpenClass(model, "class1");

        source.setArrival(jobclass1, new Exp(0.05));
        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(2.0));

        RoutingMatrix P = model.initRoutingMatrix();

        P.set(jobclass1, jobclass1, source, fork, 1);
        P.set(jobclass1, jobclass1, fork, queue1, 1);
        P.set(jobclass1, jobclass1, fork, queue2, 1);
        P.set(jobclass1, jobclass1, queue1, join, 1);
        P.set(jobclass1, jobclass1, queue2, join, 1);
        P.set(jobclass1, jobclass1, join, sink, 1);

        model.link(P);

        return model;
    }

    /**
     * Asymmetric closed fork-join network with serial queues.
     * <p>
     * Features:
     * - Single closed class with 10 jobs
     * - Delay node followed by asymmetric fork-join
     * - Branch 1: Queue1 (FCFS), Branch 2: Queue2→Queue3 serial (FCFS)
     * - Different service rates creating asymmetric load
     * - Demonstrates closed system with asymmetric parallel branches
     *
     * @return configured asymmetric fork-join model
     */
    public static Network fj_asymm() {
        Network model = new Network("model");

        Delay delay = new Delay(model, "Delay1");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 10, delay, 0);

        delay.setService(jobclass1, new Exp(0.5));
        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(2.0));
        queue3.setService(jobclass1, new Exp(1.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, delay, fork, 1.0);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1.0);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1.0);
        routingMatrix.set(jobclass1, jobclass1, queue2, queue3, 1.0);
        routingMatrix.set(jobclass1, jobclass1, queue3, join, 1.0);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Closed fork-join network with two delays and two parallel queues.
     * <p>
     * Features:
     * - Single closed class with 10 jobs
     * - Delay1→Delay2→Fork→{Queue1,Queue2}→Join→Delay1 cycle
     * - FCFS scheduling for both parallel queues
     * - Different service rates for delays and queues
     * - Demonstrates multi-stage delays before fork-join
     *
     * @return configured closed fork-join model with delays
     */
    public static Network fj_delays() {
        Network model = new Network("model");

        Delay delay = new Delay(model, "Delay1");
        Delay delay2 = new Delay(model, "Delay2");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 10, delay, 0);

        delay.setService(jobclass1, new Exp(0.5));
        delay2.setService(jobclass1, new Exp(2.0));
        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(1.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, delay, delay2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, delay2, fork, 1.0);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1.0);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1.0);
        routingMatrix.set(jobclass1, jobclass1, queue2, join, 1.0);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Complex closed fork-join with serial processing within branches.
     * <p>
     * Features:
     * - Single closed class with 10 jobs
     * - Fork to two branches with series queues within each branch
     * - Branch 1: Queue1 → Queue4 → Queue5
     * - Branch 2: Queue2 → Queue3
     * - Different service rates creating load imbalance
     * - Demonstrates fork-join with complex internal structure
     *
     * @return configured complex fork-join model
     */
    public static Network fj_complex_serial() {
        Network model = new Network("model");

        Delay delay = new Delay(model, "Delay1");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Queue queue5 = new Queue(model, "Queue5", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 10, delay, 0);

        delay.setService(jobclass1, Exp.fitMean(2.00));
        queue1.setService(jobclass1, Exp.fitMean(1.00));
        queue2.setService(jobclass1, Exp.fitMean(0.500));
        queue3.setService(jobclass1, Exp.fitMean(1.00));
        queue4.setService(jobclass1, Exp.fitMean(0.333333));
        queue5.setService(jobclass1, Exp.fitMean(1.250000));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, delay, fork, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue1, queue4, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue2, queue3, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue3, join, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue4, queue5, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue5, join, 1.00);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1.00);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1.00);

        model.link(routingMatrix);
        return model;
    }

    /**
     * Multi-class closed fork-join with asymmetric routing.
     * <p>
     * Features:
     * - Two closed classes with 10 jobs each
     * - Fork to three PS queues with asymmetric connections
     * - Queue2 → Queue3 connection creates series within branch
     * - Different service rates for each class at each queue
     * - Demonstrates multi-class load balancing in fork-join
     *
     * @return configured multi-class asymmetric fork-join model
     */
    public static Network fj_threebranches() {
        Network model = new Network("model");

        Delay delay = new Delay(model, "Delay1");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.PS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 10, delay, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "class2", 10, delay, 0);

        delay.setService(jobclass1, new Exp(0.5));
        delay.setService(jobclass2, new Exp(0.8));
        queue1.setService(jobclass1, new Exp(1.5));
        queue1.setService(jobclass2, new Exp(2.8));
        queue2.setService(jobclass1, new Exp(1.1));
        queue2.setService(jobclass2, new Exp(3.0));
        queue3.setService(jobclass1, new Exp(2.5));
        queue3.setService(jobclass2, new Exp(1.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, delay, fork, 1.00);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1.00);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue2, queue3, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue3, join, 1.00);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1.00);

        routingMatrix.set(jobclass2, jobclass2, delay, fork, 1.00);
        routingMatrix.set(jobclass2, jobclass2, fork, queue1, 1.00);
        routingMatrix.set(jobclass2, jobclass2, fork, queue2, 1.00);
        routingMatrix.set(jobclass2, jobclass2, queue1, join, 1.00);
        routingMatrix.set(jobclass2, jobclass2, queue2, queue3, 1.00);
        routingMatrix.set(jobclass2, jobclass2, queue3, join, 1.00);
        routingMatrix.set(jobclass2, jobclass2, join, delay, 1.00);

        model.link(routingMatrix);
        return model;
    }

    /**
     * Multi-class fork-join network with different task multiplicity.
     * <p>
     * Features:
     * - Two open classes with different arrival rates (0.25 each)
     * - Fork node creates 2 tasks per link for each job
     * - PS scheduling for both parallel queues
     * - Class2 has immediate service at Queue1, exponential at Queue2
     * - Demonstrates multi-class fork-join with task multiplication
     *
     * @return configured multi-class fork-join model
     */
    public static Network fj_twoclasses_forked() {
        Network model = new Network("model");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Fork fork = new Fork(model, "Fork");
        fork.setTasksPerLink(2);
        Join join = new Join(model, "Join", fork);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass1 = new OpenClass(model, "class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "class2", 0);

        source.setArrival(jobclass1, new Exp(0.25));
        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(0.75));

        source.setArrival(jobclass2, new Exp(0.25));
        queue1.setService(jobclass2, Immediate.getInstance());
        queue2.setService(jobclass2, new Exp(2.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, source, fork, 1.00);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1.00);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue2, join, 1.00);
        routingMatrix.set(jobclass1, jobclass1, join, sink, 1.00);

        routingMatrix.set(jobclass2, jobclass2, source, fork, 1.00);
        routingMatrix.set(jobclass2, jobclass2, fork, queue1, 1.00);
        routingMatrix.set(jobclass2, jobclass2, fork, queue2, 1.00);
        routingMatrix.set(jobclass2, jobclass2, queue1, join, 1.00);
        routingMatrix.set(jobclass2, jobclass2, queue2, join, 1.00);
        routingMatrix.set(jobclass2, jobclass2, join, sink, 1.00);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Complex closed fork-join network with nested fork-join structures.
     * <p>
     * Features:
     * - Two closed classes with populations [5, 2]
     * - Nested fork-join: Class2 goes through Fork1_1 then Fork1
     * - Multiple join synchronization points
     * - PS queues with different service rates per class
     * - Demonstrates complex fork-join topologies in closed systems
     *
     * @return configured nested fork-join model
     */
    public static Network fj_basic_nesting() {
        Network model = new Network("model");

        Delay delay = new Delay(model, "Delay");
        Fork fork1 = new Fork(model, "Fork1");
        fork1.setTasksPerLink(1);
        Fork fork11 = new Fork(model, "Fork1_1");
        fork11.setTasksPerLink(2);
        Join join1 = new Join(model, "Join1", fork1);
        Join join11 = new Join(model, "Join1_1", fork11);
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 5, delay, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "class2", 2, delay, 0);

        delay.setService(jobclass1, new Exp(0.25));
        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(0.75));

        delay.setService(jobclass2, new Exp(0.25));
        queue1.setService(jobclass2, new Exp(2.0));
        queue2.setService(jobclass2, new Exp(2.0));

        RoutingMatrix P = model.initRoutingMatrix();

        P.set(jobclass1, jobclass1, delay, fork1, 1.00);
        P.set(jobclass1, jobclass1, fork1, queue1, 1.00);
        P.set(jobclass1, jobclass1, fork1, queue2, 1.00);
        P.set(jobclass1, jobclass1, queue1, join1, 1.00);
        P.set(jobclass1, jobclass1, queue2, join1, 1.00);
        P.set(jobclass1, jobclass1, join1, delay, 1.00);

        P.set(jobclass2, jobclass2, delay, fork11, 1.00);
        P.set(jobclass2, jobclass2, fork11, fork1, 1.00);
        P.set(jobclass2, jobclass2, fork1, queue1, 1.00);
        P.set(jobclass2, jobclass2, fork1, queue2, 1.00);
        P.set(jobclass2, jobclass2, queue1, join1, 1.00);
        P.set(jobclass2, jobclass2, queue2, join1, 1.00);
        P.set(jobclass2, jobclass2, join1, join11, 1.00);
        P.set(jobclass2, jobclass2, join11, delay, 1.00);

        model.link(P);

        return model;
    }

    /**
     * Open fork-join network without explicit join node.
     * <p>
     * Features:
     * - Single open class with exponential arrivals
     * - Fork node splits to three parallel PS queues
     * - No explicit Join node - jobs exit independently
     * - Different service rates for each parallel branch
     * - Demonstrates fork without synchronization requirement
     *
     * @return configured fork-only network model
     */
    public static Network fj_nojoin() {
        Network model = new Network("model");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue node4 = new Queue(model, "Queue3", SchedStrategy.PS);
        Fork node5 = new Fork(model, "Fork");
        Sink node6 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "class1", 0);

        node1.setArrival(jobclass1, Exp.fitMean(2.00)); // (Source,class1)
        node2.setService(jobclass1, Exp.fitMean(1.00)); // (Queue1,class1)
        node3.setService(jobclass1, Exp.fitMean(0.500)); // (Queue2,class1)
        node4.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue3,class1)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node5, 1.00); // (Source,class1) -> (Fork,class1)
        routingMatrix.set(jobclass1, jobclass1, node5, node2, 1.0); // (Fork,class1) -> (Queue1,class1)
        routingMatrix.set(jobclass1, jobclass1, node5, node3, 1.0); // (Fork,class1) -> (Queue2,class1)
        routingMatrix.set(jobclass1, jobclass1, node5, node4, 1.0); // (Fork,class1) -> (Queue3,class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node6, 1.00); // (Queue1,class1) -> (Sink,class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node6, 1.00); // (Queue2,class1) -> (Sink,class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node6, 1.00); // (Queue3,class1) -> (Sink,class1)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Simple closed fork-join network with symmetric service rates.
     * <p>
     * Features:
     * - Single closed class with 5 jobs
     * - Delay node and two PS queues with identical service rates (1.0)
     * - Basic closed fork-join cycle: Delay → Fork → Queue1&Queue2 → Join → Delay
     * - Demonstrates balanced load sharing in closed systems
     *
     * @return configured simple closed fork-join model
     */
    public static Network fj_basic_closed() {
        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 5, delay);

        delay.setService(jobclass1, new Exp(1.0));
        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(1.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, delay, fork, 1.0);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1.0);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1.0);
        routingMatrix.set(jobclass1, jobclass1, queue2, join, 1.0);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1.0);
        model.link(routingMatrix);

        return model;
    }

    /**
     * Open network with cascaded fork-join stages.
     * <p>
     * Features:
     * - Single open class with exponential arrivals (rate 2.5)
     * - Two fork-join stages in series
     * - First stage: Fork1 → Queue1&Queue2 → Join1
     * - Second stage: Fork2 → Queue3&Queue4 → Join2
     * - All queues FCFS with identical service rates
     * - Demonstrates staged parallel processing
     *
     * @return configured cascaded fork-join model
     */
    public static Network fj_serialfjs_open() {
        Network model = new Network("model");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Fork fork1 = new Fork(model, "Fork1");
        Join join1 = new Join(model, "Join1", fork1);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Fork fork2 = new Fork(model, "Fork2");
        Join join2 = new Join(model, "Join2", fork2);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass1 = new OpenClass(model, "class1", 0);

        source.setArrival(jobclass1, Exp.fitMean(2.500));
        queue1.setService(jobclass1, Exp.fitMean(1.00));
        queue2.setService(jobclass1, Exp.fitMean(1.00));
        queue3.setService(jobclass1, Exp.fitMean(1.00));
        queue4.setService(jobclass1, Exp.fitMean(1.00));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, source, fork1, 1.00);
        routingMatrix.set(jobclass1, jobclass1, fork1, queue1, 1.00);
        routingMatrix.set(jobclass1, jobclass1, fork1, queue2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue1, join1, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue2, join1, 1.00);
        routingMatrix.set(jobclass1, jobclass1, join1, fork2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, fork2, queue3, 1.00);
        routingMatrix.set(jobclass1, jobclass1, fork2, queue4, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue3, join2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, queue4, join2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, join2, sink, 1.00);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Closed fork-join network with class switching after fork (post-fork).
     * <p>
     * Based on matlab/examples/basic/forkJoin/fj_cs_postfork.m
     * <p>
     * Features:
     * - Two closed classes with 1 job each
     * - Class 1: normal fork-join flow without class switching
     * - Class 2: switches to class 1 after the fork (post-fork switching)
     * - PS queues with exponential service rates
     * - Demonstrates fork-join with class transformation after forking
     *
     * @return configured fork-join model with post-fork class switching
     */
    public static Network fj_cs_postfork() {
        Network model = new Network("model");

        // Create nodes exactly as in MATLAB
        Delay delay = new Delay(model, "Delay");
        Fork fork1 = new Fork(model, "Fork1");
        Join join1 = new Join(model, "Join1", fork1);
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Create classes with 1 job each, reference node is delay
        ClosedClass jobclass1 = new ClosedClass(model, "class1", 1, delay);
        ClosedClass jobclass2 = new ClosedClass(model, "class2", 1, delay);

        // Set service rates exactly as in MATLAB
        delay.setService(jobclass1, new Exp(0.25));  // Exp(0.25) in MATLAB
        queue1.setService(jobclass1, new Exp(2.0)); // Exp(2.0) in MATLAB
        queue2.setService(jobclass1, new Exp(2.0)); // Exp(2.0) in MATLAB

        delay.setService(jobclass2, new Exp(0.25));  // Exp(0.25) in MATLAB
        queue1.setService(jobclass2, new Exp(2.0)); // Exp(2.0) in MATLAB  
        queue2.setService(jobclass2, new Exp(2.0)); // Exp(2.0) in MATLAB

        // Initialize routing matrix
        RoutingMatrix P = model.initRoutingMatrix();

        // Class 1 routing - normal fork-join flow without class switching
        P.set(jobclass1, jobclass1, delay, fork1, 1.0);
        P.set(jobclass1, jobclass1, fork1, queue1, 1.0);
        P.set(jobclass1, jobclass1, fork1, queue2, 1.0);
        P.set(jobclass1, jobclass1, queue1, join1, 1.0);
        P.set(jobclass1, jobclass1, queue2, join1, 1.0);
        P.set(jobclass1, jobclass1, join1, delay, 1.0);

        // Class 2 routing - switches to class 1 after fork (post-fork switching)
        P.set(jobclass2, jobclass2, delay, fork1, 1.0);
        P.set(jobclass2, jobclass1, fork1, queue1, 1.0);  // Class switch here!
        P.set(jobclass2, jobclass1, fork1, queue2, 1.0);  // Class switch here!
        // Note: No routing from queue1/queue2 to join1 for class2 
        // because class2 jobs have switched to class1

        model.link(P);

        return model;
    }

    /**
     * Open fork-join with feedback loop and class switching.
     * <p>
     * Features:
     * - Two open classes: Class1 with arrivals, Class2 generated internally
     * - Fork-join structure with PS queues
     * - Router node creates feedback from Join back to Fork
     * - Class switching capabilities at router
     * - Demonstrates complex routing in fork-join networks
     *
     * @return configured fork-join model with feedback
     */
    public static Network fj_cs_multi_visits() {
        // fork-join with multiple visits within the same chain
        Network model = new Network("model");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass1 = new OpenClass(model, "class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "class2", 0);

        source.setArrival(jobclass1, new Exp(0.1));  // Exp(0.1) in MATLAB
        queue1.setService(jobclass1, Exp.fitMean(1.0));
        queue2.setService(jobclass1, Exp.fitMean(1.0));
        queue1.setService(jobclass2, Exp.fitMean(1.0));
        queue2.setService(jobclass2, Exp.fitMean(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        
        P.set(jobclass1, jobclass1, source, fork, 1.0);
        P.set(jobclass1, jobclass1, fork, queue1, 1.0);
        P.set(jobclass1, jobclass1, fork, queue2, 1.0);
        P.set(jobclass1, jobclass1, queue1, join, 1.0);
        P.set(jobclass1, jobclass1, queue2, join, 1.0);

        // now loop back in class 2
        P.set(jobclass1, jobclass2, join, fork, 1.0);
        P.set(jobclass2, jobclass2, fork, queue1, 1.0);
        P.set(jobclass2, jobclass2, fork, queue2, 1.0);
        P.set(jobclass2, jobclass2, queue1, join, 1.0);
        P.set(jobclass2, jobclass2, queue2, join, 1.0);
        P.set(jobclass2, jobclass2, join, sink, 1.0);

        model.link(P);

        return model;
    }

    /**
     * Fork-join with overlapping routes (fj_route_overlap.m).
     * <p>
     * Two-class closed fork-join model matching MATLAB implementation.
     * <p>
     * Features:
     * - Two closed classes with 10 jobs each
     * - Delay node with different service rates per class
     * - Fork splits to two parallel FCFS queues
     * - Queue service rates: Queue1(1.0), Queue2(2.0) for both classes
     * - Join synchronizes both branches
     * - Closed cycle: Delay → Fork → Queues → Join → Delay
     *
     * @return configured two-class fork-join model
     */
    public static Network fj_route_overlap() {
        // Based on matlab/examples/basic/forkJoin/fj_route_overlap.m
        Network model = new Network("model");

        // Create nodes matching MATLAB implementation
        Delay delay = new Delay(model, "Delay1");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);

        // Create two closed classes with 10 jobs each, reference node is delay
        ClosedClass jobclass1 = new ClosedClass(model, "class1", 10, delay, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "class2", 10, delay, 0);

        // Set service rates for class1 (matching MATLAB)
        delay.setService(jobclass1, new Exp(0.5));   // Exp(0.5) in MATLAB
        queue1.setService(jobclass1, new Exp(1.0));  // Exp(1.0) in MATLAB
        queue2.setService(jobclass1, new Exp(2.0));  // Exp(2.0) in MATLAB

        // Set service rates for class2 (matching MATLAB)
        delay.setService(jobclass2, new Exp(0.2));   // Exp(0.2) in MATLAB
        queue1.setService(jobclass2, new Exp(1.0));  // Exp(1.0) in MATLAB
        queue2.setService(jobclass2, new Exp(2.0));  // Exp(2.0) in MATLAB

        // Initialize routing matrix
        RoutingMatrix P = model.initRoutingMatrix();

        // Set routing for class1 (matching MATLAB P{jobclass1, jobclass1} syntax)
        P.set(jobclass1, jobclass1, delay, fork, 1.0);
        P.set(jobclass1, jobclass1, fork, queue1, 1.0);
        P.set(jobclass1, jobclass1, fork, queue2, 1.0);
        P.set(jobclass1, jobclass1, queue1, join, 1.0);
        P.set(jobclass1, jobclass1, queue2, join, 1.0);
        P.set(jobclass1, jobclass1, join, delay, 1.0);

        // Set routing for class2 (matching MATLAB P{jobclass2, jobclass2} syntax)
        P.set(jobclass2, jobclass2, delay, fork, 1.0);
        P.set(jobclass2, jobclass2, fork, queue1, 1.0);
        P.set(jobclass2, jobclass2, fork, queue2, 1.0);
        P.set(jobclass2, jobclass2, queue1, join, 1.0);
        P.set(jobclass2, jobclass2, queue2, join, 1.0);
        P.set(jobclass2, jobclass2, join, delay, 1.0);

        model.link(P);

        return model;
    }

    /**
     * Fork-join network with class switching demonstration.
     * <p>
     * Features:
     * - Two closed classes with potential class switching
     * - Multiple delay and queue nodes with PS scheduling
     * - Fork-join structure with symmetric branches
     * - Identical service rates for class switching analysis
     * - Demonstrates class transformation in fork-join context
     *
     * @return configured fork-join model with class switching
     */
    public static Network fj_cs_prefork() {
        Network model = new Network("model");

        Delay delay = new Delay(model, "Delay1");
        Delay delay2 = new Delay(model, "Delay2");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 10, delay, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "class2", 10, delay, 0);

        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(1.0));
        delay.setService(jobclass1, new Exp(0.5));
        delay.setService(jobclass2, new Exp(0.5));
        delay2.setService(jobclass2, new Exp(2.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass2, delay, delay2, 1.0);
        routingMatrix.set(jobclass2, jobclass1, delay2, fork, 1.0);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1.0);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1.0);
        routingMatrix.set(jobclass1, jobclass1, queue2, join, 1.0);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1.0);

        model.link(routingMatrix);
        return model;
    }

    public static Network fj_deep_nesting() {
        Network model = new Network("model");

        Delay delay = new Delay(model, "Delay1");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Fork fork2 = new Fork(model, "Fork2");
        Join join2 = new Join(model, "Join2", fork2);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 1, delay, 0);

        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(1.0));
        delay.setService(jobclass1, new Exp(0.5));
        queue3.setService(jobclass1, new Exp(2.0));
        queue4.setService(jobclass1, new Exp(2.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, delay, fork, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1);
        routingMatrix.set(jobclass1, jobclass1, queue1, fork2, 1);
        routingMatrix.set(jobclass1, jobclass1, fork2, queue3, 1);
        routingMatrix.set(jobclass1, jobclass1, fork2, queue4, 1);
        routingMatrix.set(jobclass1, jobclass1, queue3, join2, 1);
        routingMatrix.set(jobclass1, jobclass1, queue4, join2, 1);
        routingMatrix.set(jobclass1, jobclass1, join2, join, 1);
        routingMatrix.set(jobclass1, jobclass1, queue2, join, 1);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1);

        model.link(routingMatrix);

        return model;
    }

    public static Network fj_serialfjs_closed() {
        Network model = new Network("model");

        Delay node1 = new Delay(model, "Delay1");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        Fork node4 = new Fork(model, "Fork");
        Join node5 = new Join(model, "Join", node4);
        Queue node6 = new Queue(model, "Queue3", SchedStrategy.PS);
        Queue node7 = new Queue(model, "Queue4", SchedStrategy.PS);
        Fork node8 = new Fork(model, "Fork2");
        Join node9 = new Join(model, "Join2", node8);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 10, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(2.00));
        node2.setService(jobclass1, Exp.fitMean(1.00));
        node3.setService(jobclass1, Exp.fitMean(1.00));
        node6.setService(jobclass1, Exp.fitMean(1.00));
        node7.setService(jobclass1, Exp.fitMean(1.00));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node4, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node2, node5, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node3, node5, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node4, node2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node4, node3, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node5, node8, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node6, node9, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node7, node9, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node8, node6, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node8, node7, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node9, node1, 1.00);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Main method for testing and demonstrating fork-join model examples.
     *
     * <p>Currently configured to:
     * - Run fj_basic_open() basic fork-join example
     * - Solve using MVA solver and display results
     * - Open JSimwView visualization (if available)
     * - Solve using JMT solver for comparison
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        Network model = fj_basic_open();
        MVA solver = new MVA(model);
        solver.getAvgTable().print();
        model.jsimwView();
        new JMT(model).getAvgTable().print();
    }

    public static Network test_forkJoinCS_1() {
        Network model = new Network("model");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass1 = new OpenClass(model, "class1");
        OpenClass jobclass2 = new OpenClass(model, "class2");

        source.setArrival(jobclass1, new Exp(0.05));
        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass2, new Exp(2.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, source, fork, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1);
        routingMatrix.set(jobclass1, jobclass2, fork, queue2, 1);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1);
        routingMatrix.set(jobclass2, jobclass1, queue2, join, 1);
        routingMatrix.set(jobclass1, jobclass1, join, sink, 1);

        model.link(routingMatrix);

        return model;
    }
}
