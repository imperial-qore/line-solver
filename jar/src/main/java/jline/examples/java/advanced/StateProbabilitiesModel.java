/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

/**
 * Examples of state probability computations
 */
public class StateProbabilitiesModel {

    /**
     * Basic closed network for state probability analysis.
     * <p>
     * Features:
     * - Two closed classes: Class1 (2 jobs), Class2 (0 jobs)
     * - Three stations: Delay, PS Queue1, 2-server PS Queue2
     * - Serial routing: Delay → Queue1 → Queue2 → Delay
     * - Different service rates for each class at each station
     * - Simple topology for state space analysis
     *
     * @return configured network for state probability computation
     */
    public static Network statepr_aggr() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        node3.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);

        node1.setService(jobclass1, new Exp(1)); // (Delay,Class1)
        node1.setService(jobclass2, new Exp(1)); // (Delay,Class2)
        node2.setService(jobclass1, new Exp(3)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(4)); // (Queue1,Class2)
        node3.setService(jobclass1, new Exp(1)); // (Queue2,Class1)
        node3.setService(jobclass2, new Exp(3)); // (Queue2,Class2)

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
     * Four-class closed network with matrix-based routing.
     * <p>
     * Features:
     * - Four closed classes with populations [1,0,4,0]
     * - Complex class switching patterns defined by routing matrices
     * - Class1 → Class2 after Queue2, Class2 → Class1 after Queue2
     * - Class3 → Class4 after Queue2, Class4 → Class3 after Queue2
     * - Demonstrates state-dependent class transformations
     *
     * @return configured multi-class state probability model
     */
    public static Network statepr_aggr_large() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        node3.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 4, node1, 0);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 0, node1, 0);

        node1.setService(jobclass1, new Exp(1)); // (Delay,Class1)
        node1.setService(jobclass2, new Exp(2)); // (Delay,Class2)
        node1.setService(jobclass3, new Exp(1)); // (Delay,Class3)
        node1.setService(jobclass4, new Exp(1)); // (Delay,Class4)
        node2.setService(jobclass1, new Exp(3)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(4)); // (Queue1,Class2)
        node2.setService(jobclass3, new Exp(5)); // (Queue1,Class3)
        node2.setService(jobclass4, new Exp(1)); // (Queue1,Class4)
        node3.setService(jobclass1, new Exp(1)); // (Queue2,Class1)
        node3.setService(jobclass2, new Exp(3)); // (Queue2,Class2)
        node3.setService(jobclass3, new Exp(5)); // (Queue2,Class3)
        node3.setService(jobclass4, new Exp(2)); // (Queue2,Class4)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        routingMatrix.set(jobclass1, jobclass2, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        routingMatrix.set(jobclass1, jobclass3, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass1, jobclass4, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));

        routingMatrix.set(jobclass2, jobclass1, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        routingMatrix.set(jobclass2, jobclass2, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        routingMatrix.set(jobclass2, jobclass3, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass2, jobclass4, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));

        routingMatrix.set(jobclass3, jobclass1, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass3, jobclass2, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass3, jobclass3, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        routingMatrix.set(jobclass3, jobclass4, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));

        routingMatrix.set(jobclass4, jobclass1, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass4, jobclass2, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass4, jobclass3, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        routingMatrix.set(jobclass4, jobclass4, new Matrix("[0,0,1; 0,0,0; 0,0,0]"));

        model.link(routingMatrix);

        return model;
    }

    /**
     * Complex four-class network with class switching.
     * <p>
     * Features:
     * - Four closed classes with populations [1,0,3,0]
     * - Three nodes: Delay, Queue1 (PS), Queue2 (2-server PS)
     * - Class1 ↔ Class2 and Class3 ↔ Class4 transformations
     * - Matrix-based routing with class switching
     * - Demonstrates complex class switching patterns
     *
     * @return configured state probability model
     */
    public static Network statepr_sys_aggr() {
        Network model = new Network("model");

        // Block 1: nodes (3 nodes only, matching MATLAB)
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        node3.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 3, node1, 0);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 0, node1, 0);

        // Service rates from MATLAB example
        node1.setService(jobclass1, new Exp(1)); // (Delay,Class1)
        node1.setService(jobclass2, new Exp(2)); // (Delay,Class2)
        node1.setService(jobclass3, new Exp(1)); // (Delay,Class3)
        node1.setService(jobclass4, new Exp(1)); // (Delay,Class4)
        node2.setService(jobclass1, new Exp(3)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(4)); // (Queue1,Class2)
        node2.setService(jobclass3, new Exp(5)); // (Queue1,Class3)
        node2.setService(jobclass4, new Exp(1)); // (Queue1,Class4)
        node3.setService(jobclass1, new Exp(1)); // (Queue2,Class1)
        node3.setService(jobclass2, new Exp(3)); // (Queue2,Class2)
        node3.setService(jobclass3, new Exp(5)); // (Queue2,Class3)
        node3.setService(jobclass4, new Exp(2)); // (Queue2,Class4)

        // Block 3: topology using matrix-based routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        // P{1,1} = [0,1,0; 0,0,1; 0,0,0]
        routingMatrix.set(jobclass1, jobclass1, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        // P{1,2} = [0,0,0; 0,0,0; 1,0,0]
        routingMatrix.set(jobclass1, jobclass2, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        // P{1,3} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass1, jobclass3, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{1,4} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass1, jobclass4, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));

        // P{2,1} = [0,0,0; 0,0,0; 1,0,0]
        routingMatrix.set(jobclass2, jobclass1, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        // P{2,2} = [0,1,0; 0,0,1; 0,0,0]
        routingMatrix.set(jobclass2, jobclass2, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        // P{2,3} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass2, jobclass3, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{2,4} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass2, jobclass4, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));

        // P{3,1} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass3, jobclass1, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{3,2} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass3, jobclass2, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{3,3} = [0,1,0; 0,0,1; 0,0,0]
        routingMatrix.set(jobclass3, jobclass3, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        // P{3,4} = [0,0,0; 0,0,0; 1,0,0]
        routingMatrix.set(jobclass3, jobclass4, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));

        // P{4,1} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass4, jobclass1, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{4,2} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass4, jobclass2, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{4,3} = [0,0,0; 0,0,0; 1,0,0]
        routingMatrix.set(jobclass4, jobclass3, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        // P{4,4} = [0,0,1; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass4, jobclass4, new Matrix("[0,0,1; 0,0,0; 0,0,0]"));

        model.link(routingMatrix);

        return model;
    }

    /**
     * Three-queue network with symmetric class populations.
     * <p>
     * Features:
     * - Four closed classes with 1 job each (symmetric populations)
     * - Three PS queues: Queue1, Queue2, 3-server Queue3
     * - Class switching patterns defined by routing matrices
     * - Class1 ↔ Class2 and Class3 ↔ Class4 transformations
     * - Demonstrates balanced state space with equal populations
     *
     * @return configured symmetric state probability model
     */
    public static Network statepr_sys_aggr_large() {
        Network model = new Network("model");

        // Block 1: nodes (3 queues only, no Delay or Router)
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.PS);
        node3.setNumberOfServers(3);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 1, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 1, node1, 0);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 1, node1, 0);

        // Service rates from MATLAB example
        node1.setService(jobclass1, new Exp(1)); // (Queue1,Class1)
        node1.setService(jobclass2, new Exp(2)); // (Queue1,Class2)
        node1.setService(jobclass3, new Exp(1)); // (Queue1,Class3)
        node1.setService(jobclass4, new Exp(1)); // (Queue1,Class4)
        node2.setService(jobclass1, new Exp(3)); // (Queue2,Class1)
        node2.setService(jobclass2, new Exp(4)); // (Queue2,Class2)
        node2.setService(jobclass3, new Exp(5)); // (Queue2,Class3)
        node2.setService(jobclass4, new Exp(1)); // (Queue2,Class4)
        node3.setService(jobclass1, new Exp(1)); // (Queue3,Class1)
        node3.setService(jobclass2, new Exp(3)); // (Queue3,Class2)
        node3.setService(jobclass3, new Exp(5)); // (Queue3,Class3)
        node3.setService(jobclass4, new Exp(2)); // (Queue3,Class4)

        // Block 3: topology using matrix-based routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        // P{1,1} = [0,1,0; 0,0,1; 0,0,0]
        routingMatrix.set(jobclass1, jobclass1, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        // P{1,2} = [0,0,0; 0,0,0; 1,0,0]
        routingMatrix.set(jobclass1, jobclass2, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        // P{1,3} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass1, jobclass3, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{1,4} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass1, jobclass4, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));

        // P{2,1} = [0,0,0; 0,0,0; 1,0,0]
        routingMatrix.set(jobclass2, jobclass1, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        // P{2,2} = [0,1,0; 0,0,1; 0,0,0]
        routingMatrix.set(jobclass2, jobclass2, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        // P{2,3} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass2, jobclass3, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{2,4} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass2, jobclass4, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));

        // P{3,1} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass3, jobclass1, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{3,2} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass3, jobclass2, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{3,3} = [0,1,0; 0,0,1; 0,0,0]
        routingMatrix.set(jobclass3, jobclass3, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        // P{3,4} = [0,0,0; 0,0,0; 1,0,0]
        routingMatrix.set(jobclass3, jobclass4, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));

        // P{4,1} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass4, jobclass1, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{4,2} = [0,0,0; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass4, jobclass2, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        // P{4,3} = [0,0,0; 0,0,0; 1,0,0]
        routingMatrix.set(jobclass4, jobclass3, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        // P{4,4} = [0,0,1; 0,0,0; 0,0,0]
        routingMatrix.set(jobclass4, jobclass4, new Matrix("[0,0,1; 0,0,0; 0,0,0]"));

        model.link(routingMatrix);

        return model;
    }

    /**
     * Two-class network with bidirectional class switching.
     * <p>
     * Features:
     * - Two closed classes: Class1 (2 jobs), Class2 (0 jobs)
     * - Class switching occurs at Queue1: both classes can transform
     * - Class1 → Class2 from Queue1 to Queue2
     * - Class2 → Class1 from Delay to Queue1
     * - Simplified topology for analyzing bidirectional switching
     *
     * @return configured bidirectional switching model
     */
    public static Network statepr_allprobs_ps() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        node3.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);

        node1.setService(jobclass1, new Exp(1)); // (Delay,Class1)
        node1.setService(jobclass2, new Exp(1)); // (Delay,Class2)
        node2.setService(jobclass1, new Exp(3)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(4)); // (Queue1,Class2)
        node3.setService(jobclass1, new Exp(1)); // (Queue2,Class1)
        node3.setService(jobclass2, new Exp(3)); // (Queue2,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.00);
        routingMatrix.set(jobclass1, jobclass2, node2, node3, 1.00);
        routingMatrix.set(jobclass2, jobclass1, node1, node2, 1.00);
        routingMatrix.set(jobclass2, jobclass1, node3, node1, 1.00);
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Mixed scheduling network with class switching.
     * <p>
     * Features:
     * - Two closed classes: Class1 (2 jobs), Class2 (0 jobs)
     * - Mixed scheduling: PS Queue1, 2-server FCFS Queue2
     * - Same bidirectional class switching pattern as example 5
     * - Demonstrates impact of scheduling strategy on state probabilities
     * - FCFS vs PS comparison with identical routing
     *
     * @return configured mixed scheduling state probability model
     */
    public static Network statepr_allprobs_fcfs() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        node3.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);

        node1.setService(jobclass1, new Exp(1)); // (Delay,Class1)
        node1.setService(jobclass2, new Exp(1)); // (Delay,Class2)
        node2.setService(jobclass1, new Exp(3)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(4)); // (Queue1,Class2)
        node3.setService(jobclass1, new Exp(3)); // (Queue2,Class1)
        node3.setService(jobclass2, new Exp(3)); // (Queue2,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.00);
        routingMatrix.set(jobclass1, jobclass2, node2, node3, 1.00);
        routingMatrix.set(jobclass2, jobclass1, node1, node2, 1.00);
        routingMatrix.set(jobclass2, jobclass1, node3, node1, 1.00);
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Main method for testing and demonstrating state probability examples.
     *
     * <p>Currently configured to create statepr_aggr_large()
     * without running any solvers or analysis.
     *
     * @param args command line arguments (not used)
     * @throws Exception if model creation fails
     */
    public static void main(String[] args) throws Exception {
        Network model = statepr_aggr_large();
    }
}
