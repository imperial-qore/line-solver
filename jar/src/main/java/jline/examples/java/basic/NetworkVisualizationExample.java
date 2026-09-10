/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.Network;
import jline.lang.ClosedClass;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Exp;

/**
 * Demonstrates JUNG visualization of Network models using the plot() method.
 *
 * <p>This example creates queueing network models and displays them
 * using the built-in plot() method. The visualization showcases the structure
 * of queueing networks.
 *
 * <p>The visualization uses shape coding:
 * <ul>
 * <li>White right-triangle: Source nodes</li>
 * <li>Black left-triangle: Sink nodes</li>
 * <li>White rectangle: Queue stations</li>
 * <li>White circle: Delay stations</li>
 * <li>White diamond: Fork/Join nodes</li>
 * <li>White hexagon: ClassSwitch nodes</li>
 * </ul>
 *
 * <p>All edges are black solid lines indicating routing paths.
 */
public class NetworkVisualizationExample {

    /**
     * Creates a simple open queueing network with Source, Queue, Delay, and Sink.
     *
     * @return configured open network model
     */
    public static Network createSimpleOpenNetwork() {
        Network model = new Network("SimpleOpenNetwork");

        // Create nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Delay delay = new Delay(model, "Delay");
        Sink sink = new Sink(model, "Sink");

        // Create job class
        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        // Set service processes
        source.setArrival(jobClass, new Exp(0.5));
        queue.setService(jobClass, new Exp(1.0));
        delay.setService(jobClass, new Exp(2.0));

        // Create routing matrix: Source -> Queue -> Delay -> Sink
        RoutingMatrix P = model.initRoutingMatrix();
        P.addConnection(source, queue, jobClass);
        P.addConnection(queue, delay, jobClass);
        P.addConnection(delay, sink, jobClass);
        model.link(P);

        return model;
    }

    /**
     * Creates a network with a Fork-Join structure.
     *
     * @return configured fork-join network model
     */
    public static Network createForkJoinNetwork() {
        Network model = new Network("ForkJoinNetwork");

        // Create nodes
        Source source = new Source(model, "Source");
        Fork fork = new Fork(model, "Fork");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Join join = new Join(model, "Join", fork);
        Sink sink = new Sink(model, "Sink");

        // Create job class
        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        // Set service processes
        source.setArrival(jobClass, new Exp(0.2));
        queue1.setService(jobClass, new Exp(1.0));
        queue2.setService(jobClass, new Exp(1.5));

        // Create routing matrix
        RoutingMatrix P = model.initRoutingMatrix();
        P.addConnection(source, fork, jobClass);
        P.addConnection(fork, queue1, jobClass);
        P.addConnection(fork, queue2, jobClass);
        P.addConnection(queue1, join, jobClass);
        P.addConnection(queue2, join, jobClass);
        P.addConnection(join, sink, jobClass);
        model.link(P);

        return model;
    }

    /**
     * Creates a more complex network with multiple paths.
     *
     * @return configured multi-path network model
     */
    public static Network createMultiPathNetwork() {
        Network model = new Network("MultiPathNetwork");

        // Create nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Delay delay = new Delay(model, "Delay");
        Sink sink = new Sink(model, "Sink");

        // Create job class
        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        // Set service processes
        source.setArrival(jobClass, new Exp(0.3));
        queue1.setService(jobClass, new Exp(1.0));
        queue2.setService(jobClass, new Exp(1.5));
        delay.setService(jobClass, new Exp(0.5));

        // Create routing matrix with probabilistic routing
        RoutingMatrix P = model.initRoutingMatrix();
        P.addConnection(source, queue1, jobClass, 0.5);
        P.addConnection(source, queue2, jobClass, 0.5);
        P.addConnection(queue1, delay, jobClass);
        P.addConnection(queue2, delay, jobClass);
        P.addConnection(delay, sink, jobClass);
        model.link(P);

        return model;
    }

    /**
     * Creates a larger mixed-class network with open and closed classes.
     *
     * @return configured mixed-class network model
     */
    public static Network createMixedClassNetwork() {
        Network model = new Network("MixedClassNetwork");

        // Create nodes: Source/Sink for open class, multiple queues and delays
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");

        // Application tier
        Queue webServer = new Queue(model, "WebServer", SchedStrategy.PS);
        Queue appServer1 = new Queue(model, "AppServer1", SchedStrategy.FCFS);
        Queue appServer2 = new Queue(model, "AppServer2", SchedStrategy.FCFS);

        // Database tier
        Queue dbPrimary = new Queue(model, "DB_Primary", SchedStrategy.FCFS);
        Queue dbReplica = new Queue(model, "DB_Replica", SchedStrategy.FCFS);

        // Cache layer
        Delay cache = new Delay(model, "Cache");

        // Background processing
        Queue batchQueue = new Queue(model, "BatchQueue", SchedStrategy.FCFS);
        Delay thinkTime = new Delay(model, "ThinkTime");

        // Create open class (external requests)
        OpenClass requests = new OpenClass(model, "Requests", 0);

        // Create closed class (background batch jobs)
        ClosedClass batchJobs = new ClosedClass(model, "BatchJobs", 5, thinkTime, 0);

        // Set arrival rate for open class
        source.setArrival(requests, new Exp(2.0));

        // Set service times for queues
        webServer.setService(requests, new Exp(10.0));
        appServer1.setService(requests, new Exp(5.0));
        appServer2.setService(requests, new Exp(5.0));
        dbPrimary.setService(requests, new Exp(20.0));
        dbReplica.setService(requests, new Exp(15.0));
        cache.setService(requests, new Exp(50.0));

        // Batch jobs service times
        batchQueue.setService(batchJobs, new Exp(2.0));
        dbPrimary.setService(batchJobs, new Exp(10.0));
        thinkTime.setService(batchJobs, new Exp(0.1));

        // Create routing matrix
        RoutingMatrix P = model.initRoutingMatrix();

        // Open class routing: Source -> WebServer -> (AppServer1 | AppServer2) -> Cache -> (DB_Primary | DB_Replica) -> Sink
        P.addConnection(source, webServer, requests);
        P.addConnection(webServer, appServer1, requests, 0.5);
        P.addConnection(webServer, appServer2, requests, 0.5);
        P.addConnection(appServer1, cache, requests);
        P.addConnection(appServer2, cache, requests);
        P.addConnection(cache, dbPrimary, requests, 0.3);
        P.addConnection(cache, dbReplica, requests, 0.2);
        P.addConnection(cache, sink, requests, 0.5);  // Cache hit - skip DB
        P.addConnection(dbPrimary, sink, requests);
        P.addConnection(dbReplica, sink, requests);

        // Closed class routing: ThinkTime -> BatchQueue -> DB_Primary -> ThinkTime
        P.addConnection(thinkTime, batchQueue, batchJobs);
        P.addConnection(batchQueue, dbPrimary, batchJobs);
        P.addConnection(dbPrimary, thinkTime, batchJobs);

        model.link(P);

        return model;
    }

    /**
     * Creates a central server model (closed network with CPU and multiple disks).
     *
     * @return configured central server network model
     */
    public static Network createCentralServerModel() {
        Network model = new Network("CentralServerModel");

        // Think time (terminals)
        Delay thinkTime = new Delay(model, "Terminals");

        // Central CPU
        Queue cpu = new Queue(model, "CPU", SchedStrategy.PS);

        // Multiple disk drives
        Queue disk1 = new Queue(model, "Disk1", SchedStrategy.FCFS);
        Queue disk2 = new Queue(model, "Disk2", SchedStrategy.FCFS);
        Queue disk3 = new Queue(model, "Disk3", SchedStrategy.FCFS);

        // Closed class: N users
        ClosedClass users = new ClosedClass(model, "Users", 20, thinkTime, 0);

        // Service times
        thinkTime.setService(users, new Exp(0.1));   // 10 sec think time
        cpu.setService(users, new Exp(20.0));         // 50 ms CPU time
        disk1.setService(users, new Exp(10.0));       // 100 ms disk time
        disk2.setService(users, new Exp(10.0));
        disk3.setService(users, new Exp(8.0));

        // Routing: Terminals -> CPU -> (Disks or back to Terminals) -> CPU
        RoutingMatrix P = model.initRoutingMatrix();

        // From terminals to CPU
        P.addConnection(thinkTime, cpu, users);

        // From CPU: 20% exit to terminals, 80% go to disks
        P.addConnection(cpu, thinkTime, users, 0.2);
        P.addConnection(cpu, disk1, users, 0.3);
        P.addConnection(cpu, disk2, users, 0.3);
        P.addConnection(cpu, disk3, users, 0.2);

        // From disks back to CPU
        P.addConnection(disk1, cpu, users);
        P.addConnection(disk2, cpu, users);
        P.addConnection(disk3, cpu, users);

        model.link(P);

        return model;
    }

    /**
     * Main method demonstrating the plot() method on Network models.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("=== Central Server Model Visualization ===\n");

        // Show central server model
        System.out.println("Visualizing central server model...");
        Network centralServerModel = createCentralServerModel();
        centralServerModel.plot("Central Server Model", 900, 700);

        System.out.println("\nVisualization window opened.");
        System.out.println("Node shapes: Queue=buffer+server, Delay=circle");
        System.out.println("Close window to exit.");
    }
}
