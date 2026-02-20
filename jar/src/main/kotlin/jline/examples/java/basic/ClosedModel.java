/*
 * Copyright (c) 2012-2024, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.ClassSwitch;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.*;
import jline.solvers.jmt.JMT;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.lang.processes.APH;

import java.util.Arrays;

/**
 * Examples of closed queueing networks
 */
public class ClosedModel {

    /**
     * Simple closed network with a delay node and FCFS queue.
     * <p>
     * Features:
     * - Single closed class with 10 jobs
     * - Delay node with exponential service (mean 1.0)
     * - FCFS queue with exponential service (mean 1.5)
     * - Probabilistic routing: 70% self-loop at delay, 30% to queue
     *
     * @return configured closed network model
     */
    public static Network cqn_repairmen() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.00)); // (Delay,Class1)
        node2.setService(jobclass1, Exp.fitMean(1.50)); // (Queue1,Class1)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.70); // (Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.3); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.00); // (Queue1,Class1) -> (Delay,Class1)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Two-class closed network with class switching.
     * <p>
     * Features:
     * - Two closed classes with 2 jobs each
     * - Delay node with Erlang and HyperExp distributions
     * - PS (Processor Sharing) queue
     * - Class switching between Class1 and Class2
     *
     * @return configured closed network model
     */
    public static Network cqn_twoclass_hyperl() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,Class2)
        node2.setService(jobclass1, new HyperExp(0.1, 1.0, 10.0)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(1.0)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.3);
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.1);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 0.2);

        routingMatrix.set(jobclass1, jobclass2, node1, node1, 0.6);
        routingMatrix.set(jobclass1, jobclass2, node2, node1, 0.8);

        routingMatrix.set(jobclass2, jobclass1, node2, node1, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Three-class closed network with multi-server PS queue.
     * <p>
     * Features:
     * - Three closed classes with populations [2, 0, 1]
     * - Two-server PS queue
     * - Mixed service distributions (Erlang, HyperExp, Exp)
     * - Complex class switching patterns
     * - Uses circular routing matrix for Class3
     *
     * @return configured closed network model
     */
    public static Network cqn_threeclass_hyperl() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        node2.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 1, node1, 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,Class2)
        node1.setService(jobclass3, Exp.fitMean(1.00)); // (Delay,Class3)

        node2.setService(jobclass1, new HyperExp(0.1, 1.0, 10.0)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(2.0)); // (Queue1,Class2)
        node2.setService(jobclass3, new Exp(3.0)); // (Queue1,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.3);
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.1);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 0.2);

        routingMatrix.set(jobclass1, jobclass2, node1, node1, 0.6);
        routingMatrix.set(jobclass1, jobclass2, node2, node1, 0.8);

        routingMatrix.set(jobclass2, jobclass1, node2, node1, 1);

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1);

        Matrix c = Maths.circul(model.getNumberOfStations());

        routingMatrix.set(jobclass3, jobclass3, node1, node1, c.value());
        routingMatrix.set(jobclass3, jobclass3, node1, node2, c.get(0, 1));
        routingMatrix.set(jobclass3, jobclass3, node2, node1, c.get(1, 0));
        routingMatrix.set(jobclass3, jobclass3, node2, node2, c.get(1, 1));

        model.link(routingMatrix);

        return model;
    }

    /**
     * Four-class closed network with multiple FCFS queues.
     * <p>
     * Features:
     * - Four closed classes with populations [2, 2, 2, 1]
     * - Three nodes: Delay, Queue1 (3 servers), Queue2 (3 servers)
     * - Queue2 disabled for Class1 and Class2
     * - Complex class switching topology
     * - Mix of exponential and Erlang service distributions
     *
     * @return configured closed network model
     */
    public static Network cqn_multiserver() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        node2.setNumberOfServers(3);
        node3.setNumberOfServers(3);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 2, node1, 0);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 1, node1, 0);

        node1.setService(jobclass1, new Exp(1)); // (Delay,Class1)
        node1.setService(jobclass2, new Exp(1)); // (Delay,Class2)
        node1.setService(jobclass3, new Exp(10)); // (Delay,Class3)
        node1.setService(jobclass4, new Exp(1)); // (Delay,Class4)

        node2.setService(jobclass1, new Exp(1)); // (Queue1,Class1)
        node2.setService(jobclass2, new Erlang(1, 2)); // (Queue1,Class2)
        node2.setService(jobclass3, new Exp(10)); // (Queue1,Class3)
        node2.setService(jobclass4, new Exp(1)); // (Queue1,Class4)

        node3.setService(jobclass1, Disabled.getInstance()); // (Queue2,Class1)
        node3.setService(jobclass2, Disabled.getInstance()); // (Queue2,Class2)
        node3.setService(jobclass3, new Erlang(1, 2)); // (Queue2,Class3)
        node3.setService(jobclass4, new Exp(1)); // (Queue2,Class4)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.50);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.00);

        routingMatrix.set(jobclass1, jobclass2, node1, node2, 0.50);

        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.00);
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.00);

        routingMatrix.set(jobclass2, jobclass1, node1, node2, 1.00);

        routingMatrix.set(jobclass3, jobclass3, node1, node2, 0.25);
        routingMatrix.set(jobclass3, jobclass3, node1, node3, 0.25);
        routingMatrix.set(jobclass3, jobclass3, node2, node1, 1.00);
        routingMatrix.set(jobclass3, jobclass3, node3, node1, 1.00);

        routingMatrix.set(jobclass3, jobclass4, node1, node2, 0.50);

        routingMatrix.set(jobclass4, jobclass4, node2, node1, 1.00);
        routingMatrix.set(jobclass4, jobclass4, node3, node1, 1.00);

        routingMatrix.set(jobclass4, jobclass3, node1, node2, 1.00);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Two-class closed network with four-node cyclic topology.
     * <p>
     * Features:
     * - Two closed classes with populations [1, 2]
     * - Four nodes: two delays and two PS queues
     * - Cyclic routing: Delay1 → Delay2 → Queue1 → Queue2 → Delay1
     * - Different service rates for each class-node combination
     * - Uses Network.cyclicPsInf(N,D,Z) with:
     *   - D = [10,5; 5,9] (Queue service times)
     *   - Z = [91,92; 93,94] (Delay times)
     *   - N = [1,2] (Population)
     *
     * @return configured closed network model
     */
    public static Network cqn_oneline() {
        // D(i,r) - mean service time of class r at queue i
        Matrix D = new Matrix(2, 2);
        D.set(0, 0, 10.0); D.set(0, 1, 5.0);
        D.set(1, 0, 5.0);  D.set(1, 1, 9.0);

        // N(r) - number of jobs of class r
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 1.0); N.set(0, 1, 2.0);

        // Z(i,r) - mean think time of class r at delay station i
        Matrix Z = new Matrix(2, 2);
        Z.set(0, 0, 91.0); Z.set(0, 1, 92.0);
        Z.set(1, 0, 93.0); Z.set(1, 1, 94.0);

        return Network.cyclicPsInf(N, D, Z);
    }

    /**
     * Closed network with class switching using probabilistic routing.
     * <p>
     * Features:
     * - Two closed classes with populations [15, 5]
     * - ClassSwitch node with class switching matrix [0,1;1,0]
     * - Mixed scheduling: PS queues and infinite server delay
     * - Probabilistic routing strategies (RAND, RROBIN, WRROBIN)
     * - Demonstrates advanced routing configuration
     *
     * @return configured closed network model
     */
    public static Network cqn_twoclass_erl() {
        Network model = new Network("Model");

        Matrix csMatrix = new Matrix("[0,1;1,0]");

        ClassSwitch node1 = new ClassSwitch(model, "CS", csMatrix);
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue node4 = new Queue(model, "Delay", SchedStrategy.INF);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 15, node4, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 5, node4, 0);

        node2.setService(jobclass1, Exp.fitMean(1.50)); // (Queue1,Class1)
        node2.setService(jobclass2, Erlang.fitMeanAndOrder(1.5, 2)); // (Queue1,Class2)

        node3.setService(jobclass1, Erlang.fitMeanAndOrder(1.5, 2)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(1.50)); // (Queue2,Class2)
        node4.setService(jobclass1, Exp.fitMean(1.00)); // (Delay,Class1)
        node4.setService(jobclass2, Exp.fitMean(1.00)); // (Delay,Class2)

        // Block 3: topology
        model.addLink(node2, node1);
        model.addLink(node3, node1);
        model.addLink(node1, node4);
        model.addLink(node4, node2);
        model.addLink(node4, node3);

        node1.setRouting(jobclass1, RoutingStrategy.RAND);
        node2.setRouting(jobclass1, RoutingStrategy.RAND);
        node3.setRouting(jobclass1, RoutingStrategy.RAND);
        node4.setRouting(jobclass1, RoutingStrategy.RROBIN);

        node1.setRouting(jobclass2, RoutingStrategy.RAND);
        node2.setRouting(jobclass2, RoutingStrategy.RAND);
        node3.setRouting(jobclass2, RoutingStrategy.RAND);
        node4.setRouting(jobclass2, RoutingStrategy.WRROBIN, node2, 1);
        node4.setRouting(jobclass2, RoutingStrategy.WRROBIN, node3, 2);

        return model;
    }

    /**
     * Two-class FCFS network demonstrating First-Come-First-Served scheduling.
     * <p>
     * Features:
     * - Two closed classes with 2 jobs each
     * - Delay node with Erlang and HyperExp distributions
     * - FCFS queue with exponential service
     * - Simple serial routing between delay and queue
     *
     * @return configured closed network model
     */
    public static Network cqn_bcmp_theorem_fcfs() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(1.00)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.00)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.00); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.00); // (Queue1,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Two-class LCFSPR network demonstrating Last-Come-First-Served Preemptive Resume.
     * <p>
     * Features:
     * - Two closed classes with 2 jobs each
     * - Delay node with APH (Acyclic Phase-type) distributions
     * - LCFSPR queue with preemptive resume scheduling
     * - Different SCV values for delay services
     *
     * @return configured closed network model
     */
    public static Network cqn_bcmp_theorem_lcfspr() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.LCFSPR);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(1.00)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.00)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.00); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.00); // (Queue1,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Simple two-class PS network for processor sharing demonstration.
     * <p>
     * Features:
     * - Two closed classes with 2 jobs each
     * - Delay node with Erlang and HyperExp distributions
     * - PS queue with exponential service for both classes
     * - Simple serial routing between delay and queue
     *
     * @return configured closed network model
     */
    public static Network cqn_bcmp_theorem_ps() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,Class2)
        node2.setService(jobclass1, new Exp(1)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(1)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.00); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.00); // (Queue1,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Two-class closed network with multi-server queue.
     * <p>
     * Features:
     * - Two closed classes with populations [4, 2]
     * - Delay node with exponential service for both classes
     * - FCFS queue with 3 servers
     * - Different service rates (Class2 has 10x faster service)
     *
     * @return configured closed network model
     */
    public static Network cqn_repairmen_multi() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        node2.setNumberOfServers(3);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 4, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, new Exp(1)); // (Delay,Class1)
        node1.setService(jobclass2, new Exp(1));

        node2.setService(jobclass1, new Exp(1)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(10));


        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.00); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.00); // (Queue1,Class2) -> (Delay,Class2)
        model.link(routingMatrix);

        return model;
    }

    /**
     * Two-class three-node closed network with serial topology.
     * <p>
     * Features:
     * - Two closed classes with 10 jobs each
     * - Three nodes: Delay, Queue1, Queue2 in series
     * - All exponential service times with different rates
     * - Cyclic routing: Delay → Queue1 → Queue2 → Delay
     *
     * @return configured closed network model
     */
    public static Network cqn_twoqueues_multi() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 10, node1, 0);


        node1.setService(jobclass1, Exp.fitMean(1.00));
        node2.setService(jobclass1, Exp.fitMean(1.50));
        node3.setService(jobclass1, Exp.fitMean(3.00));
        node1.setService(jobclass2, Exp.fitMean(1.00));
        node2.setService(jobclass2, Exp.fitMean(1.50));
        node3.setService(jobclass2, Exp.fitMean(3.00));

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.00); // (Queue1,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.00); // (Queue2,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue1,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.00); // (Queue2,Class2) -> (Delay,Class2)
        model.link(routingMatrix);

        return model;
    }

    /**
     * Single-class closed network with reducible routing matrix.
     * <p>
     * Features:
     * - Single closed class with 1 job
     * - Three nodes: Delay, Queue1 (FCFS), Queue2 (FCFS)
     * - Reducible routing matrix where Queue1 and Queue2 have self-loops
     * - Initial routing from Delay splits between Queue1 (30%) and Queue2 (50%)
     * - Demonstrates handling of absorbing states in routing
     *
     * @return configured closed network model
     */
    public static Network cqn_twoqueues() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.00)); // (Delay,Class1)
        node2.setService(jobclass1, Exp.fitMean(1.50)); // (Queue1,Class1)
        node3.setService(jobclass1, Exp.fitMean(3.00)); // (Queue2,Class1)

        // Block 3: topology - reducible routing matrix
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.20); // (Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.30); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node3, 0.50); // (Delay,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node2, 1.00); // (Queue1,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node3, 1.00); // (Queue2,Class1) -> (Queue2,Class1)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Mixed scheduling strategies with probabilistic routing.
     * <p>
     * Features:
     * - Two closed classes: Class1 (2 jobs), Class2 (1 job)
     * - Three stations: Delay, PS Queue, DPS Queue
     * - Probabilistic routing with opposite preferences per class
     * - Class1: 30% to Queue1, 70% to Queue2
     * - Class2: 70% to Queue1, 30% to Queue2
     * - Demonstrates DPS (Discriminatory Processor Sharing) scheduling
     *
     * @return configured mixed scheduling network model
     */
    public static Network cqn_scheduling_dps() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.DPS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 1, node1, 0);

        node1.setService(jobclass1, new Exp(3)); // (Delay,Class1) - rate 3
        node1.setService(jobclass2, new Exp(0.5)); // (Delay,Class2) - rate 0.5

        node2.setService(jobclass1, new Exp(0.1), 5.0); // (Queue1,Class1) - rate 0.1, weight 5 (PS ignores weights)
        node2.setService(jobclass2, new Exp(1), 1.0); // (Queue1,Class2) - rate 1, weight 1 (PS ignores weights)
        node3.setService(jobclass1, new Exp(0.1), 1.0); // (Queue2,Class1) - rate 0.1, weight 1
        node3.setService(jobclass2, new Exp(1), 5.0); // (Queue2,Class2) - rate 1, weight 5

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.300); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node3, 0.700); // (Delay,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.00); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.00); // (Queue2,Class1) -> (Delay,Class1)
        
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 0.700); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node1, node3, 0.300); // (Delay,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.00); // (Queue1,Class2) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.00); // (Queue2,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    /**
     * APH service distributions with complex routing patterns.
     * <p>
     * Features:
     * - Two closed classes (1 job each) with APH service distributions
     * - All service times fitted using APH with specific mean and SCV values
     * - Class1: probabilistic routing (0% self-loop, 30% Queue1, 70% Queue2)
     * - Class2: random routing strategy
     * - High variability service times (SCV up to 5.038)
     * - Demonstrates APH fitting for realistic service modeling
     *
     * @return configured APH-based network model
     */
    public static Network cqn_mmpp2_service() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 1, node1, 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,Class2)
        node2.setService(jobclass1, new HyperExp(0.1, 1.0, 10.0)); // (Queue1,Class1)
        node2.setService(jobclass2, new MMPP2(1, 2, 3, 4)); // (Queue1,Class2)
        node3.setService(jobclass1, new HyperExp(0.1, 1.0, 10.0)); // (Queue2,Class1)
        node3.setService(jobclass2, new Erlang(1, 2)); // (Queue2,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        // Class 1 routing probabilities
        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.3);
        routingMatrix.set(jobclass1, jobclass1, node1, node3, 0.7);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.0);

        // Class 2 routing probabilities (RAND = uniform among connected neighbors)
        // node1 connects to node1, node2, node3 -> 1/3 each
        routingMatrix.set(jobclass2, jobclass2, node1, node1, 1.0 / 3.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0 / 3.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node3, 1.0 / 3.0);
        // node2 connects only to node1 -> 1.0
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.0);
        // node3 connects only to node1 -> 1.0
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Two-class closed network with LCFS and LCFSPR scheduling.
     * <p>
     * Features:
     * - Two closed classes with 1 job each
     * - Queue1 with LCFS (Last Come First Served) scheduling
     * - Queue2 with LCFSPR (LCFS with Preemptive Resume) scheduling
     * - Exponential service times at all stations
     * - Cyclic routing between the two queues
     *
     * @return configured LCFS/LCFSPR network model
     */
    public static Network cqn_lcfs_lcfspr() {
        Network model = new Network("model");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.LCFS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.LCFSPR);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 1, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(2.000000)); // (Queue1,Class1)
        node1.setService(jobclass2, Exp.fitMean(3.000000)); // (Queue1,Class2)
        node2.setService(jobclass1, Exp.fitMean(5.000000)); // (Queue2,Class1)
        node2.setService(jobclass2, Exp.fitMean(7.000000)); // (Queue2,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Queue1,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue2,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Queue1,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue2,Class2) -> (Queue1,Class2)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Creates a closed queueing network with LCFS and LCFSPR scheduling - 3 class variant.
     *
     * Model structure:
     * - Queue1: LCFS (Last Come First Served) scheduling
     * - Queue2: LCFSPR (Last Come First Served Preemptive Resume) scheduling
     * - Three closed classes with 1 job each
     * - Cyclic routing between the two queues
     *
     * Expected CTMC ground truth:
     * Queue1, Class1: QLen=0.34369, Util=0.16883, RespT=4.0714, Tput=0.084416
     * Queue1, Class2: QLen=0.28527, Util=0.14752, RespT=5.8013, Tput=0.049173
     * Queue1, Class3: QLen=0.17611, Util=0.14794, RespT=15.476, Tput=0.01138
     * Queue2, Class1: QLen=0.65631, Util=0.42208, RespT=7.7747, Tput=0.084416
     * Queue2, Class2: QLen=0.71473, Util=0.34421, RespT=14.535, Tput=0.049173
     * Queue2, Class3: QLen=0.82389, Util=0.12518, RespT=72.401, Tput=0.01138
     *
     * @return configured 3-class LCFS/LCFSPR network model
     */
    public static Network cqn_lcfs_lcfspr_3class() {
        Network model = new Network("model");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.LCFS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.LCFSPR);

        // Block 2: classes (3 closed classes, 1 job each)
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 1, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 1, node1, 0);

        // Service times (mean values)
        node1.setService(jobclass1, Exp.fitMean(2.0));   // (Queue1,Class1)
        node1.setService(jobclass2, Exp.fitMean(3.0));   // (Queue1,Class2)
        node1.setService(jobclass3, Exp.fitMean(13.0));  // (Queue1,Class3)
        node2.setService(jobclass1, Exp.fitMean(5.0));   // (Queue2,Class1)
        node2.setService(jobclass2, Exp.fitMean(7.0));   // (Queue2,Class2)
        node2.setService(jobclass3, Exp.fitMean(11.0));  // (Queue2,Class3)

        // Block 3: topology - cyclic routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.0);
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.0);
        routingMatrix.set(jobclass3, jobclass3, node2, node1, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Creates a closed queueing network with LCFS and LCFSPR scheduling - 4 class variant.
     *
     * Model structure:
     * - Queue1: LCFS (Last Come First Served) scheduling
     * - Queue2: LCFSPR (Last Come First Served Preemptive Resume) scheduling
     * - Four closed classes with 1 job each
     * - Cyclic routing between the two queues
     *
     * Expected CTMC ground truth:
     * Queue1, Class1: QLen=0.3969, Util=0.1551, RespT=5.1181, Tput=0.077549
     * Queue1, Class2: QLen=0.33858, Util=0.13422, RespT=7.5679, Tput=0.044739
     * Queue1, Class3: QLen=0.19303, Util=0.12099, RespT=20.741, Tput=0.0093067
     * Queue1, Class4: QLen=0.14419, Util=0.095919, RespT=25.555, Tput=0.0056423
     * Queue2, Class1: QLen=0.6031, Util=0.38774, RespT=7.777, Tput=0.077549
     * Queue2, Class2: QLen=0.66142, Util=0.31317, RespT=14.784, Tput=0.044739
     * Queue2, Class3: QLen=0.80697, Util=0.10237, RespT=86.709, Tput=0.0093067
     * Queue2, Class4: QLen=0.85581, Util=0.1072, RespT=151.68, Tput=0.0056423
     *
     * @return configured 4-class LCFS/LCFSPR network model
     */
    public static Network cqn_lcfs_lcfspr_4class() {
        Network model = new Network("model");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.LCFS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.LCFSPR);

        // Block 2: classes (4 closed classes, 1 job each)
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 1, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 1, node1, 0);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 1, node1, 0);

        // Service times (mean values)
        node1.setService(jobclass1, Exp.fitMean(2.0));   // (Queue1,Class1)
        node1.setService(jobclass2, Exp.fitMean(3.0));   // (Queue1,Class2)
        node1.setService(jobclass3, Exp.fitMean(13.0));  // (Queue1,Class3)
        node1.setService(jobclass4, Exp.fitMean(17.0));  // (Queue1,Class4)
        node2.setService(jobclass1, Exp.fitMean(5.0));   // (Queue2,Class1)
        node2.setService(jobclass2, Exp.fitMean(7.0));   // (Queue2,Class2)
        node2.setService(jobclass3, Exp.fitMean(11.0));  // (Queue2,Class3)
        node2.setService(jobclass4, Exp.fitMean(19.0));  // (Queue2,Class4)

        // Block 3: topology - cyclic routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.0);
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.0);
        routingMatrix.set(jobclass3, jobclass3, node2, node1, 1.0);
        routingMatrix.set(jobclass4, jobclass4, node1, node2, 1.0);
        routingMatrix.set(jobclass4, jobclass4, node2, node1, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Main method for testing and demonstrating closed model examples.
     *
     * <p>Currently configured to run cqn_threeclass_hyperl() and solve it
     * using the JMT solver to display system-level metrics.
     *
     * @param args command line arguments (not used)
     * @throws Exception if solver encounters an error
     */
    public static void main(String[] args) throws Exception {
        Network model = cqn_threeclass_hyperl();
//        jlqnGetPath();
//        new JMT(model).getAvgTable().print();
//        new MVA(model, new MVAOptions().verbose(VerboseLevel.DEBUG)).getAvgTable().print();
//        new MVA(model).getAvgChainTable().print();
//        new MVA(model).getAvgNodeTable().print();
//        new MVA(model).getAvgNodeChainTable().print();
//        NetworkStruct sn = model.getStruct();
        new JMT(model).getAvgSysTable().print();
    }

}
