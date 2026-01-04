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
import jline.lang.processes.*;
import jline.solvers.nc.NC;

/**
 * Examples of mixed queueing networks
 */
public class MixedModel {

    /**
     * Simple mixed open/closed network with delay and PS queue.
     * <p>
     * Features:
     * - One closed class (2 jobs) and one open class
     * - Delay node with Erlang service for closed class, HyperExp for open
     * - PS queue shared by both classes with different service distributions
     * - Closed class circulates between delay and queue
     * - Open class flows: Source → Delay → Queue → Sink
     *
     * @return configured mixed network model
     */
    public static Network mqn_basic() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Source node3 = new Source(model, "Source");
        Sink node4 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 2, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,ClosedClass)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,OpenClass)
        node2.setService(jobclass1, new HyperExp(0.1, 1.0, 10.0)); // (Queue1,ClosedClass)
        node2.setService(jobclass2, new Exp(1)); // (Queue1,OpenClass)
        node3.setArrival(jobclass2, new Exp(0.1)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Delay,ClosedClass) -> (Queue1,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.00); // (Queue1,ClosedClass) -> (Delay,ClosedClass)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Delay,OpenClass) -> (Queue1,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node4, 1.00); // (Queue1,OpenClass) -> (Sink,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.00); // (Source,OpenClass) -> (Delay,OpenClass)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Complex mixed network with five PS queues and different server counts.
     * <p>
     * Features:
     * - One closed class (3 jobs) and one open class
     * - Five PS queues with 1, 2, 3, 4, 5 servers respectively
     * - Closed class cycles through Queue1→Queue2→Queue3→Queue4→Queue1
     * - Open class follows Queue1→Queue2→Queue3→Queue5→Sink
     * - Different service rates optimized for server counts
     *
     * @return configured mixed network model
     */
    public static Network mqn_multiserver_ps() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.PS);
        node1.setNumberOfServers(1);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.PS);
        node2.setNumberOfServers(2);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.PS);
        node3.setNumberOfServers(3);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.PS);
        node4.setNumberOfServers(4);
        Queue node5 = new Queue(model, "Queue5", SchedStrategy.PS);
        node5.setNumberOfServers(5);
        Source node6 = new Source(model, "Source");
        Sink node7 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 3, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, new Exp(1)); // (Queue1,ClosedClass)
        node1.setService(jobclass2, new Exp(1)); // (Queue1,OpenClass)
        node2.setService(jobclass1, new Exp(2)); // (Queue2,ClosedClass)
        node2.setService(jobclass2, new Exp(Math.sqrt(2))); // (Queue2,OpenClass)
        node3.setService(jobclass1, new Exp(3)); // (Queue3,ClosedClass)
        node3.setService(jobclass2, new Exp(Math.sqrt(3))); // (Queue3,OpenClass)
        node4.setService(jobclass1, new Exp(4)); // (Queue4,ClosedClass)
        node4.setService(jobclass2, new Exp(2)); // (Queue4,OpenClass)
        node5.setService(jobclass1, new Exp(5)); // (Queue5,ClosedClass)
        node5.setService(jobclass2, new Exp(Math.sqrt(5))); // (Queue5,OpenClass)

        node6.setArrival(jobclass2, new Exp(0.3)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Queue1,ClosedClass) -> (Queue2,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.00); // (Queue2,ClosedClass) -> (Queue3,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.00); // (Queue3,ClosedClass) -> (Queue4,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.00); // (Queue4,ClosedClass) -> (Queue1,ClosedClass)

        routingMatrix.set(jobclass2, jobclass2, node6, node1, 1.00); // (Source,OpenClass) -> (Queue1,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Queue1,OpenClass) -> (Queue2,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue2,OpenClass) -> (Queue3,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node5, 1.00); // (Queue3,OpenClass) -> (Queue5,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node5, node7, 1.00); // (Queue5,OpenClass) -> (Sink,OpenClass)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Mixed network similar to example 2 but with FCFS scheduling.
     * <p>
     * Features:
     * - Same topology as mqn_multiserver_ps() but FCFS queues
     * - One closed class (3 jobs) and one open class
     * - Five FCFS queues with varying server counts (1-5)
     * - Demonstrates scheduling strategy differences in mixed networks
     * - Identical service rates and routing patterns
     *
     * @return configured mixed network model
     */
    public static Network mqn_multiserver_fcfs() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        node2.setNumberOfServers(2);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        node3.setNumberOfServers(3);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        node4.setNumberOfServers(4);
        Queue node5 = new Queue(model, "Queue5", SchedStrategy.FCFS);
        node5.setNumberOfServers(5);
        Source node6 = new Source(model, "Source");
        Sink node7 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 3, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, new Exp(1)); // (Queue1,ClosedClass)
        node1.setService(jobclass2, new Exp(1)); // (Queue1,OpenClass)
        node2.setService(jobclass1, new Exp(2)); // (Queue2,ClosedClass)
        node2.setService(jobclass2, new Exp(Math.sqrt(2))); // (Queue2,OpenClass)
        node3.setService(jobclass1, new Exp(3)); // (Queue3,ClosedClass)
        node3.setService(jobclass2, new Exp(Math.sqrt(3))); // (Queue3,OpenClass)
        node4.setService(jobclass1, new Exp(4)); // (Queue4,ClosedClass)
        node4.setService(jobclass2, new Exp(2)); // (Queue4,OpenClass)
        node5.setService(jobclass1, new Exp(5)); // (Queue5,ClosedClass)
        node5.setService(jobclass2, new Exp(Math.sqrt(5))); // (Queue5,OpenClass)
        node6.setArrival(jobclass1, Disabled.getInstance()); // (Source,ClosedClass)
        node6.setArrival(jobclass2, new Exp(0.3)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Queue1,ClosedClass) -> (Queue2,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.00); // (Queue2,ClosedClass) -> (Queue3,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.00); // (Queue3,ClosedClass) -> (Queue4,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.00); // (Queue4,ClosedClass) -> (Queue1,ClosedClass)

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Queue1,OpenClass) -> (Queue2,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue2,OpenClass) -> (Queue3,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node5, 1.00); // (Queue3,OpenClass) -> (Queue5,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node5, node7, 1.00); // (Queue5,OpenClass) -> (Sink,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node6, node1, 1.00); // (Source,OpenClass) -> (Queue1,OpenClass)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Mixed network with large closed class population and APH arrivals.
     * <p>
     * Features:
     * - Large closed class with 100 jobs circulating through 4 FCFS queues
     * - Open class with APH (Acyclic Phase-type) arrival process
     * - High variability arrivals (SCV = 64.0) for bursty traffic
     * - Open class exits after Queue3, closed class cycles through all 4
     * - Demonstrates heavy traffic mixed system behavior
     *
     * @return configured mixed network model
     */
    public static Network mqn_singleserver_fcfs() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Source node5 = new Source(model, "Source");
        Sink node6 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 100, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, new Exp(1)); // (Queue1,ClosedClass)
        node1.setService(jobclass2, new Exp(1)); // (Queue1,OpenClass)
        node2.setService(jobclass1, new Exp(2)); // (Queue2,ClosedClass)
        node2.setService(jobclass2, new Exp(Math.sqrt(2))); // (Queue2,OpenClass)
        node3.setService(jobclass1, new Exp(3)); // (Queue3,ClosedClass)
        node3.setService(jobclass2, new Exp(Math.sqrt(3))); // (Queue3,OpenClass)
        node4.setService(jobclass1, new Exp(4)); // (Queue4,ClosedClass)
        node4.setService(jobclass2, new Exp(2)); // (Queue4,OpenClass)
        node5.setArrival(jobclass1, Disabled.getInstance()); // (Source,ClosedClass)
        node5.setArrival(jobclass2, APH.fitMeanAndSCV(3.000000, 64.00)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Queue1,ClosedClass) -> (Queue2,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.00); // (Queue2,ClosedClass) -> (Queue3,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.00); // (Queue3,ClosedClass) -> (Queue4,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.00); // (Queue4,ClosedClass) -> (Queue1,ClosedClass)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Queue1,OpenClass) -> (Queue2,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue2,OpenClass) -> (Queue3,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node6, 1.00); // (Queue3,OpenClass) -> (Sink,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node4, node1, 1.00); // (Queue4,OpenClass) -> (Queue1,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node5, node1, 1.00); // (Source,OpenClass) -> (Queue1,OpenClass)


        model.link(routingMatrix);

        return model;
    }

    /**
     * Mixed PS network with large closed population and simplified open routing.
     * <p>
     * Features:
     * - Large closed class (100 jobs) with 4 PS queues in cycle
     * - Open class with exponential arrivals (rate 0.3)
     * - Open class follows shortened path: Queue1→Queue2→Queue3→Sink
     * - PS scheduling for all queues
     * - Demonstrates processor sharing in mixed environments
     *
     * @return configured mixed network model
     */
    public static Network mqn_singleserver_ps() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.PS);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.PS);
        Source node5 = new Source(model, "Source");
        Sink node6 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 100, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, new Exp(1.0)); // (Queue1,ClosedClass) - lambda=1.0
        node1.setService(jobclass2, new Exp(1.0)); // (Queue1,OpenClass) - lambda=1.0
        node2.setService(jobclass1, new Exp(2.0)); // (Queue2,ClosedClass) - lambda=2.0
        node2.setService(jobclass2, new Exp(1.414213562373)); // (Queue2,OpenClass) - lambda=1.414213562373
        node3.setService(jobclass1, new Exp(3.0)); // (Queue3,ClosedClass) - lambda=3.0
        node3.setService(jobclass2, new Exp(1.732050807569)); // (Queue3,OpenClass) - lambda=1.732050807569
        node4.setService(jobclass1, new Exp(4.0)); // (Queue4,ClosedClass) - lambda=4.0
        node4.setService(jobclass2, new Exp(2.0)); // (Queue4,OpenClass) - lambda=2.0

        node5.setArrival(jobclass1, Disabled.getInstance()); // (Source,ClosedClass)
        node5.setArrival(jobclass2, new Exp(0.3)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Queue1,ClosedClass) -> (Queue2,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.00); // (Queue2,ClosedClass) -> (Queue3,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.00); // (Queue3,ClosedClass) -> (Queue4,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.00); // (Queue4,ClosedClass) -> (Queue1,ClosedClass)


        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Queue1,OpenClass) -> (Queue2,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue2,OpenClass) -> (Queue3,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node6, 1.00); // (Queue3,OpenClass) -> (Sink,OpenClass)

        routingMatrix.set(jobclass2, jobclass2, node5, node1, 1.00); // (Source,OpenClass) -> (Queue1,OpenClass)

        model.link(routingMatrix);

        return model;
    }

    /**
     * Main method for testing and demonstrating mixed model examples.
     *
     * <p>Currently configured to run mqn_singleserver_fcfs() and solve it
     * using the NC (Normalizing Constant) solver.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        Network model = mqn_singleserver_fcfs();
        new NC(model).getAvgTable().print();
    }

}
