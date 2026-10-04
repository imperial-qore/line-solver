/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.ClassSwitch;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.mva.MVA;

/**
 * Examples of models with class switching
 */
public class ClassSwitchingModel {

    /**
     * Class switching example with probabilistic transitions.
     * <p>
     * Features:
     * - Two open classes with different arrival rates (10.0 and 2.0)
     * - ClassSwitch node with probabilistic class transitions
     * - Class1 → Class1 (30%) or Class1 → Class2 (70%)
     * - Class2 → Class1 (100%)
     * - FCFS queue with identical service rates for both classes
     * - Demonstrates class switching in open networks
     *
     * @return configured class switching network model
     */
    public static Network cs_implicit() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source 1");
        Queue node2 = new Queue(model, "Queue 1", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink 1");
        ClassSwitch node4 = new ClassSwitch(model, "ClassSwitch 1"); // Dummy node, class switching is embedded in the routing matrix P

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, Exp.fitMean(10.00)); // (Source 1,Class1)
        node1.setArrival(jobclass2, Exp.fitMean(2.00));  // (Source 1,Class2)
        node2.setService(jobclass1, Exp.fitMean(1.00));  // (Queue 1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.00));  // (Queue 1,Class2)

        ClassSwitchMatrix csMatrix = node4.initClassSwitchMatrix();
        csMatrix.set(jobclass1, jobclass1, 0.3);
        csMatrix.set(jobclass1, jobclass2, 0.7);
        csMatrix.set(jobclass2, jobclass1, 1.0);
        csMatrix.set(jobclass2, jobclass2, 0.0);
        node4.setClassSwitchingMatrix(csMatrix);

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node4, 1.00); // (Source 1,Class1) -> (ClassSwitch 1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.00); // (Queue 1,Class1) -> (Sink 1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node2, 1.00); // (ClassSwitch 1,Class1) -> (Queue 1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node4, 1.00); // (Source 1,Class2) -> (ClassSwitch 1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue 1,Class2) -> (Sink 1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node2, 1.00); // (ClassSwitch 1,Class2) -> (Queue 1,Class2)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Multi-diamond class switching example.
     * <p>
     * Features:
     * - Three FCFS queues with different service rates
     * - Three open classes with only Class1 having external arrivals
     * - Class switching through routing matrix creating multiple paths
     * - Demonstrates implicit class switching without ClassSwitch node
     * 
     * @return configured class switching network model
     */
    public static Network cs_multi_diamond() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source 1");
        Queue node2 = new Queue(model, "Queue 0", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue 1", SchedStrategy.FCFS);
        Queue node4 = new Queue(model, "Queue 2", SchedStrategy.FCFS);
        Sink node5 = new Sink(model, "Sink 1");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);
        OpenClass jobclass3 = new OpenClass(model, "Class3", 0);

        node1.setArrival(jobclass1, new Exp(1.0)); // (Source 1,Class1)
        node2.setService(jobclass1, new Exp(10.0)); // (Queue 0,Class1)
        node3.setService(jobclass2, new Exp(20.0)); // (Queue 1,Class2)
        node4.setService(jobclass3, new Exp(30.0)); // (Queue 2,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0); // Source -> Queue 0
        routingMatrix.set(jobclass1, jobclass1, node2, node2, 0.2); // Queue 0 -> Queue 0 (self-loop)
        routingMatrix.set(jobclass1, jobclass2, node2, node3, 0.3); // Queue 0 -> Queue 1 (class switch)
        routingMatrix.set(jobclass1, jobclass3, node2, node4, 0.5); // Queue 0 -> Queue 2 (class switch)
        routingMatrix.set(jobclass2, jobclass2, node3, node5, 1.0); // Queue 1 -> Sink
        routingMatrix.set(jobclass3, jobclass3, node4, node5, 1.0); // Queue 2 -> Sink

        model.link(routingMatrix);

        return model;
    }

    /**
     * Single-diamond class switching example.
     * <p>
     * Features:
     * - Three Delay nodes (infinite servers) in a closed network
     * - Three closed classes with only Class1 having initial population
     * - Creates a diamond-shaped flow pattern through class switching
     * - Jobs can switch from Class1 to Class2/3, then return to Class1
     * 
     * @return configured class switching network model
     */
    public static Network cs_single_diamond() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Queue 0");
        Delay node2 = new Delay(model, "Queue 1");
        Delay node3 = new Delay(model, "Queue 2");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 0, node1);

        node1.setService(jobclass1, Exp.fitMean(1.0)); // (Queue 0,Class1)
        node2.setService(jobclass2, Exp.fitMean(2.0)); // (Queue 1,Class2)
        node3.setService(jobclass3, Exp.fitMean(3.0)); // (Queue 2,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        
        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.2); // Queue 0 -> Queue 0 (stay in Class1)
        routingMatrix.set(jobclass1, jobclass2, node1, node2, 0.3); // Queue 0 -> Queue 1 (switch to Class2)
        routingMatrix.set(jobclass1, jobclass3, node1, node3, 0.5); // Queue 0 -> Queue 2 (switch to Class3)
        routingMatrix.set(jobclass2, jobclass1, node2, node1, 1.0); // Queue 1 -> Queue 0 (return to Class1)
        routingMatrix.set(jobclass3, jobclass1, node3, node1, 1.0); // Queue 2 -> Queue 0 (return to Class1)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Transient class switching example.
     * <p>
     * Features:
     * - Three Delay nodes with Class1 served at all nodes
     * - Class2 and Class3 only served at specific nodes
     * - Creates an absorbing Markov chain where jobs eventually get stuck
     * - Once jobs switch to Class2 or Class3, they remain in those classes
     * 
     * @return configured class switching network model
     */
    public static Network cs_transient_class() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Queue 0");
        Delay node2 = new Delay(model, "Queue 1");
        Delay node3 = new Delay(model, "Queue 2");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 0, node1);

        node1.setService(jobclass1, Exp.fitMean(1.0)); // (Queue 0,Class1)
        node1.setService(jobclass2, Exp.fitMean(1.0)); // (Queue 0,Class2)
        node1.setService(jobclass3, Exp.fitMean(1.0)); // (Queue 0,Class3)
        node2.setService(jobclass2, Exp.fitMean(1.0)); // (Queue 1,Class2)
        node3.setService(jobclass3, Exp.fitMean(1.0)); // (Queue 2,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        
        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.2); // Queue 0 -> Queue 0 (stay in Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0); // Queue 0 -> Queue 1 (Class2)
        routingMatrix.set(jobclass3, jobclass3, node1, node3, 1.0); // Queue 0 -> Queue 2 (Class3)
        routingMatrix.set(jobclass1, jobclass2, node1, node2, 0.3); // Queue 0 -> Queue 1 (switch to Class2)
        routingMatrix.set(jobclass1, jobclass3, node1, node3, 0.5); // Queue 0 -> Queue 2 (switch to Class3)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.0); // Queue 1 -> Queue 0 (stay in Class2)
        routingMatrix.set(jobclass3, jobclass3, node3, node1, 1.0); // Queue 2 -> Queue 0 (stay in Class3)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Main method for testing and demonstrating class switching examples.
     *
     * <p>Currently configured to:
     * - Run cs_implicit() with class switching demonstration
     * - Print network structure and routing matrix details
     * - Solve using MVA solver and display various performance tables
     * - Show average metrics, system metrics, and chain metrics
     *
     * @param args command line arguments (not used)
     * @throws Exception if solver encounters an error
     */
    public static void main(String[] args) throws Exception {
        Network model = cs_implicit();

        NetworkStruct sn = model.getStruct(false);
        sn.rt.print();
        model.printRoutingMatrix();
        MVA solver = new MVA(model);
        solver.getAvgTable().print();
        solver.getAvgSysTable().print();
        solver.getAvgChainTable().print();
    }
}