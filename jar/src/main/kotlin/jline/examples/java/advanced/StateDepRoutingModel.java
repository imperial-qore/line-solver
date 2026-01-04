/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.solvers.jmt.JMT;
import jline.util.matrix.Matrix;

/**
 * Examples of models with state-dependent routing
 */
public class StateDepRoutingModel {

    /**
     * Single-class closed network with round-robin routing.
     * <p>
     * Features:
     * - Single closed class with 1 job
     * - Delay node with high-variability HyperExp service (SCV=25)
     * - Two PS queues with different service rates (1.0 and 2.0)
     * - Round-robin routing from delay to queues
     * - Demonstrates basic state-dependent routing behavior
     *
     * @return configured state-dependent routing model
     */
    public static Network sdroute_closed() {
        Network model = new Network("model");

        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);

        node1.setService(jobclass1, HyperExp.fitMeanAndSCV(1, 25)); // (Delay,Class1)
        node2.setService(jobclass1, new Exp(1.0)); // (Queue1,Class1)
        node3.setService(jobclass1, new Exp(2.0)); // (Queue2,Class1)

        model.addLink(node1, node1);
        model.addLink(node1, node2);
        model.addLink(node1, node3);
        model.addLink(node2, node1);
        model.addLink(node3, node1);

        node1.setRouting(jobclass1, RoutingStrategy.RROBIN);
        node2.setProbRouting(jobclass1, node1, 1.0);
        node3.setProbRouting(jobclass1, node1, 1.0);

        return model;
    }

    /**
     * Multi-class closed network with MAP and APH service processes.
     * <p>
     * Features:
     * - Two closed classes: Class1 (1 job), Class2 (2 jobs)
     * - Advanced service processes: APH, PH, MAP distributions
     * - Renewal (APH, PH) and non-renewal (MAP) service processes
     * - Mixed scheduling: PS Queue1, FCFS Queue2
     * - Round-robin routing for both classes from delay
     * - Demonstrates advanced stochastic processes in routing
     *
     * @return configured advanced state-dependent routing model
     */
    public static Network sdroute_twoclasses_closed() {
        Network model = new Network("model");

        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

//renewal
        APH map21 = new APH(new Matrix("[1, 0]"), new Matrix("[-2, 2; 0, -0.5]"));
        PH map22 = new PH(Matrix.singleton(1.0), Matrix.singleton(-1));
//non - renewal
        MAP map31 = new MAP(new Matrix("[-20, 0; 0, -1]"), new Matrix("[0, 20; 0.8, 0.2]"));
        MAP map32 = new MAP(new Matrix("[-4, 3; 4, -6]"), new Matrix("[1, 0;0, 2]"));
        //MMPP2 map32 = new MMPP2(1, 2, 3, 4);

        node1.setService(jobclass1, HyperExp.fitMeanAndSCV(1, 25));
        node2.setService(jobclass1, map21);
        node3.setService(jobclass1, map31);

        node1.setService(jobclass2, HyperExp.fitMeanAndSCV(1, 25));
        node2.setService(jobclass2, map22);
        node3.setService(jobclass2, map32);

        model.addLink(node1, node1);
        model.addLink(node1, node2);
        model.addLink(node1, node3);
        model.addLink(node2, node1);
        model.addLink(node3, node1);

        node1.setRouting(jobclass1, RoutingStrategy.RROBIN);
        node2.setProbRouting(jobclass1, node1, 1.0);
        node3.setProbRouting(jobclass1, node1, 1.0);

        node1.setRouting(jobclass2, RoutingStrategy.RROBIN);
        node2.setProbRouting(jobclass2, node1, 1.0);
        node3.setProbRouting(jobclass2, node1, 1.0);
        return model;
    }

    /**
     * Open network with load balancing router.
     * <p>
     * Features:
     * - Open network: Source → Router → Queues → Sink
     * - Router uses round-robin strategy for load balancing
     * - Two FCFS queues with identical service rates (0.5)
     * - Exponential arrivals (rate 1.0) distributed evenly
     * - Demonstrates state-dependent load balancing in open systems
     *
     * @return configured open load balancing model
     */
    public static Network sdroute_open() {
        Network model = new Network("model");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Router router = new Router(model, "Router");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass oclass = new OpenClass(model, "Class1", 0);

        source.setArrival(oclass, new Exp(1.00)); // (Source,Class1)
        queue1.setService(oclass, new Exp(2.00)); // (Queue1,Class1)
        queue2.setService(oclass, new Exp(2.00)); // (Queue2,Class1)

        // Block 3: topology
        model.addLink(source, router);
        model.addLink(router, queue1);
        model.addLink(router, queue2);
        model.addLink(queue1, sink);
        model.addLink(queue2, sink);

        router.setRouting(oclass, RoutingStrategy.RROBIN);

        return model;
    }

    /**
     * Main method for testing and demonstrating state-dependent routing examples.
     *
     * <p>Currently configured to:
     * - Run sdroute_twoclasses_closed() with advanced service processes
     * - Solve using JMT solver and print average performance metrics
     *
     * @param args command line arguments (not used)
     * @throws Exception if solver encounters an error
     */
    public static void main(String[] args) throws Exception {
        Network model = sdroute_twoclasses_closed();
        new JMT(model).getAvgTable().print();
    }

}