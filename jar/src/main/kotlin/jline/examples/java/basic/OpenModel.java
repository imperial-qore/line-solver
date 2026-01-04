/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.VerboseLevel;
import jline.lang.ClassSwitchMatrix;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Disabled;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.lang.processes.Replayer;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.CTMC;
import jline.solvers.fluid.FLD;
import jline.solvers.jmt.JMT;
import jline.solvers.mam.MAM;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.ssa.SSA;
import jline.util.matrix.Matrix;

/**
 * Examples of open queueing networks
 */
public class OpenModel {
    /**
     * Simple open network with delay and queue in series.
     * <p>
     * Features:
     * - Single open class with exponential arrivals (rate 0.1)
     * - Delay node with HyperExp service distribution
     * - FCFS queue with exponential service
     * - Serial routing: Source → Delay → Queue → Sink
     *
     * @return configured open network model
     */
    public static Network oqn_basic() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Source node3 = new Source(model, "Source");
        Sink node4 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);

        node1.setService(jobclass1, new HyperExp(0.5, 3.0, 10.0));
        node2.setService(jobclass1, new Exp(1.0));
        node3.setArrival(jobclass1, new Exp(0.1));

        // Block 3: topology
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass1, jobclass1, new Matrix("[0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0]"));
        model.link(P);

        return model;
    }

    /**
     * Tandem PS network with delay station using matrix notation.
     * <p>
     * Features:
     * - Uses Matrix-based constructor for PS queues with delay
     * - Lambda matrix defines arrival rates for two classes
     * - D matrix defines service demands at stations
     * - Z matrix defines service times at delay stations
     *
     * @return configured tandem PS network
     */
    public static Network oqn_oneline() {
        Matrix lambda = new Matrix("[1,2]");
        lambda.scaleEq(1.0 / 50); // arrival rate of class r
        Matrix D = new Matrix("[10,5;5,9]");  // mean service time of class r at station i
        Matrix Z = new Matrix("[91,92]"); // Z(r)  mean service time of class r at delay station i
        return Network.tandemPsInf(lambda, D, Z);
    }

    /**
     * Three-class open network with class switching.
     * <p>
     * Features:
     * - Three open classes (A, B, C) with Class C having no arrivals
     * - Class switching from A→C and B→C at ClassSwitch node
     * - Two PS queues with different service rates per class
     * - Demonstrates class transformation in open networks
     *
     * @return configured open network with class switching
     */
    public static Network oqn_cs_routing() {
        Network model = new Network("model");

        // Block 1: nodes
        Source node1 = new Source(model, "Source 1");
        Queue node2 = new Queue(model, "Queue 1", SchedStrategy.PS);
        ClassSwitch node3 = new ClassSwitch(model, "ClassSwitch 1");
        Sink node4 = new Sink(model, "Sink 1");
        Queue node5 = new Queue(model, "Queue 2", SchedStrategy.PS);

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class A", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class B", 0);
        OpenClass jobclass3 = new OpenClass(model, "Class C", 0);

        node1.setArrival(jobclass1, Exp.fitMean(0.500)); // (Source 1,Class A)
        node1.setArrival(jobclass2, Exp.fitMean(1.00)); // (Source 1,Class B)
        node1.setArrival(jobclass3, Disabled.getInstance()); // (Source 1,Class C)
        node2.setService(jobclass1, Exp.fitMean(0.200)); // (Queue 1,Class A)
        node2.setService(jobclass2, Exp.fitMean(0.300)); // (Queue 1,Class B)
        node2.setService(jobclass3, Exp.fitMean(0.333333)); // (Queue 1,Class C)
        node5.setService(jobclass1, Exp.fitMean(1.00)); // (Queue 2,Class A)
        node5.setService(jobclass2, Exp.fitMean(1.00)); // (Queue 2,Class B)
        node5.setService(jobclass3, Exp.fitMean(0.150000)); // (Queue 2,Class C)

        ClassSwitchMatrix C = node3.initClassSwitchMatrix();
        C = new ClassSwitchMatrix(Matrix.eye(C.getNumRows()));
        node3.setClassSwitchingMatrix(C);

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source 1,Class A) -> (Queue 1,Class A)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.00); // (Queue 1,Class A) -> (ClassSwitch 1,Class A)
        routingMatrix.set(jobclass1, jobclass1, node5, node4, 1.00); // (Queue 2,Class A) -> (Sink 1,Class A)
        routingMatrix.set(jobclass1, jobclass3, node3, node5, 1.00); // (CS_ClassSwitch 1_to_Queue 2,Class A) -> (Queue 2,Class C)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Source 1,Class B) -> (Queue 1,Class B)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue 1,Class B) -> (ClassSwitch 1,Class B)
        routingMatrix.set(jobclass2, jobclass2, node5, node4, 1.00); // (Queue 2,Class B) -> (Sink 1,Class B)
        routingMatrix.set(jobclass2, jobclass3, node3, node5, 1.00); // (CS_ClassSwitch 1_to_Queue 2,Class B) -> (Queue 2,Class C)
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.00); // (Source 1,Class C) -> (Queue 1,Class C)
        routingMatrix.set(jobclass3, jobclass3, node2, node3, 1.00); // (Queue 1,Class C) -> (ClassSwitch 1,Class C)
        routingMatrix.set(jobclass3, jobclass3, node3, node5, 1.00); // (ClassSwitch 1,Class C) -> (CS_ClassSwitch 1_to_Queue 2,Class C)
        routingMatrix.set(jobclass3, jobclass3, node5, node4, 1.00); // (Queue 2,Class C) -> (Sink 1,Class C)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Open network with trace-driven service using Replayer.
     * <p>
     * Features:
     * - Single open class with exponential arrivals
     * - Queue service driven by trace file (example_trace.txt)
     * - Demonstrates empirical service time distributions
     * - Simple Source → Queue → Sink topology
     *
     * @return configured trace-driven network model
     */
    public static Network oqn_trace_driven() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "OpenClass", 0);

        node1.setArrival(jobclass1, Exp.fitMean(1.00)); // (Source,OpenClass)
        // Load trace file from resources
        String traceFile = OpenModel.class.getResource("/example_trace.txt").getPath();
        node2.setService(jobclass1, new Replayer(traceFile)); // (Queue,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,OpenClass) -> (Queue,OpenClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.00); // (Queue,OpenClass) -> (Sink,OpenClass)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Open network with probabilistic routing to multiple virtual sinks.
     * <p>
     * Features:
     * - Two open classes with different routing patterns
     * - Class1: 60% to VSink1, 40% to VSink2
     * - Class2: 10% to VSink1, 90% to VSink2
     * - Router nodes as intermediate destinations
     * - Demonstrates probabilistic routing patterns
     *
     * @return configured network with probabilistic routing
     */
    public static Network oqn_vsinks() {
        Network model = new Network("model");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        Router vsink1 = new Router(model, "VSink1");
        Router vsink2 = new Router(model, "VSink2");

        // Block 2: classes
        OpenClass ocl1 = new OpenClass(model, "Class1", 0);
        OpenClass ocl2 = new OpenClass(model, "Class2", 0);

        source.setArrival(ocl1, new Exp(1.0));
        queue.setService(ocl1, new Exp(100.0));

        source.setArrival(ocl2, new Exp(1.0));
        queue.setService(ocl2, new Exp(100.0));

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(ocl1, ocl1, source, queue, 1.0);
        routingMatrix.set(ocl1, ocl1, queue, vsink1, 0.6);
        routingMatrix.set(ocl1, ocl1, queue, vsink2, 0.4);
        routingMatrix.set(ocl1, ocl1, vsink1, sink, 1.0);
        routingMatrix.set(ocl1, ocl1, vsink2, sink, 1.0);

        routingMatrix.set(ocl2, ocl2, source, queue, 1.0);
        routingMatrix.set(ocl2, ocl2, queue, vsink1, 0.1);
        routingMatrix.set(ocl2, ocl2, queue, vsink2, 0.9);
        routingMatrix.set(ocl2, ocl2, vsink1, sink, 1.0);
        routingMatrix.set(ocl2, ocl2, vsink2, sink, 1.0);

        model.link(routingMatrix);

        return model;
    }

    /**
     * Complex multi-class open network with four queues.
     * <p>
     * Features:
     * - Three open classes with different arrival rates
     * - Four queues: FCFS (Queue1,2,4) and PS (Queue3)
     * - Complex routing patterns with probabilistic splits
     * - Different service rates for each class at each queue
     *
     * @return configured multi-class open network
     */
    public static Network oqn_fourqueues() {
        Network model = new Network("MyNetwork");
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue node4 = new Queue(model, "Queue3", SchedStrategy.PS);
        Queue node5 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Sink node6 = new Sink(model, "Sink");

        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);
        OpenClass jobclass3 = new OpenClass(model, "Class3", 0);

        node1.setArrival(jobclass1, Exp.fitMean(5));
        node1.setArrival(jobclass2, Exp.fitMean(8));
        node1.setArrival(jobclass3, Exp.fitMean(7));
        node2.setService(jobclass1, Exp.fitMean(0.3));
        node2.setService(jobclass2, Exp.fitMean(0.5));
        node2.setService(jobclass3, Exp.fitMean(0.6));
        node3.setService(jobclass1, Exp.fitMean(1.1));
        node3.setService(jobclass2, Exp.fitMean(1.3));
        node3.setService(jobclass3, Exp.fitMean(1.5));
        node4.setService(jobclass1, Exp.fitMean(2.0));
        node4.setService(jobclass2, Exp.fitMean(2.1));
        node4.setService(jobclass3, Exp.fitMean(1.9));
        node5.setService(jobclass1, Exp.fitMean(1.5));
        node5.setService(jobclass2, Exp.fitMean(0.9));
        node5.setService(jobclass3, Exp.fitMean(2.3));

        RoutingMatrix p = new RoutingMatrix(model, model.getJobClasses(), model.getNodes());
        p.addConnection(jobclass1, jobclass1, node1, node2, 1);
        p.addConnection(jobclass1, jobclass1, node2, node3, 0.25);
        p.addConnection(jobclass1, jobclass1, node2, node4, 0.25);
        p.addConnection(jobclass1, jobclass1, node2, node5, 0.25);
        p.addConnection(jobclass1, jobclass1, node2, node6, 0.25);
        p.addConnection(jobclass1, jobclass1, node3, node2, 1);
        p.addConnection(jobclass1, jobclass1, node4, node2, 1);
        p.addConnection(jobclass1, jobclass1, node5, node2, 1);
        p.addConnection(jobclass2, jobclass2, node1, node2, 1);
        p.addConnection(jobclass2, jobclass2, node2, node3, 0.25);
        p.addConnection(jobclass2, jobclass2, node2, node4, 0.25);
        p.addConnection(jobclass2, jobclass2, node2, node5, 0.25);
        p.addConnection(jobclass2, jobclass2, node2, node6, 0.25);
        p.addConnection(jobclass2, jobclass2, node3, node2, 1);
        p.addConnection(jobclass2, jobclass2, node4, node2, 1);
        p.addConnection(jobclass2, jobclass2, node5, node2, 1);
        p.addConnection(jobclass3, jobclass3, node1, node2, 1);
        p.addConnection(jobclass3, jobclass3, node2, node3, 0.25);
        p.addConnection(jobclass3, jobclass3, node2, node4, 0.25);
        p.addConnection(jobclass3, jobclass3, node2, node5, 0.25);
        p.addConnection(jobclass3, jobclass3, node2, node6, 0.25);
        p.addConnection(jobclass3, jobclass3, node3, node2, 1);
        p.addConnection(jobclass3, jobclass3, node4, node2, 1);
        p.addConnection(jobclass3, jobclass3, node5, node2, 1);

        model.link(p);

        return model;
    }

    /**
     * Main method for testing and demonstrating open model examples.
     *
     * <p>Currently configured to run oqn_cs_routing() and solve it
     * using multiple solvers (CTMC, Fluid, MVA, MAM, NC, JMT, SSA) for comparison.
     * This matches the structure of the oqn_cs_routing.ipynb notebook.
     *
     * @param args command line arguments (not used)
     * @throws Exception if solver encounters an error
     */
    public static void main(String[] args) throws Exception {
        Network model = oqn_cs_routing();
        
        // Solver options configuration
        SolverOptions options = new SolverOptions();
        options.keep = true;
        options.verbose = VerboseLevel.STD;
        options.cutoff = new Matrix(new double[][]{{1,1,0},{3,3,0},{0,0,3}});
        options.seed = 23000;
        options.samples = 100000;
        
        // Create array of solvers to test
        NetworkSolver[] solvers = {
            new CTMC(model, options),
            new FLD(model, options),
            new MVA(model, options),
            new MAM(model, options),
            new NC(model, options),
            new JMT(model, options),
            new SSA(model, options)
        };
        
        // Execute all solvers and print results
        for (int i = 0; i < solvers.length; i++) {
            System.out.println("\n=== " + solvers[i].getClass().getSimpleName() + " Results ===");
            try {
                solvers[i].getAvgTable().print();
            } catch (Exception e) {
                System.out.println("Solver failed: " + e.getMessage());
            }
        }
    }

}
