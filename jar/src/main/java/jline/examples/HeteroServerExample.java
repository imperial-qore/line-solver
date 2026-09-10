/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples;

import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.solvers.NetworkAvgTable;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.util.matrix.Matrix;

import java.util.HashMap;
import java.util.Map;

/**
 * Example demonstrating heterogeneous server support in LINE.
 *
 * Heterogeneous servers allow defining queues with multiple server types,
 * each with different service rates and class compatibilities.
 */
public class HeteroServerExample {

    public static void main(String[] args) {
        example1_basic();
        example2_compatibility();
        example3_policies();
        example4_parallelism();
    }

    /**
     * Example 1: Basic heterogeneous server queue
     *
     * A queue with two server types:
     * - Fast servers (2 servers, rate 2.0)
     * - Slow servers (3 servers, rate 1.0)
     */
    public static void example1_basic() {
        System.out.println("=== Example 1: Basic Heterogeneous Server Queue ===\n");

        Network model = new Network("HeteroBasic");

        // Create nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "HeteroQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Create job class
        OpenClass jobClass = new OpenClass(model, "Jobs");

        // Set arrival rate
        source.setArrival(jobClass, new Exp(3.0));  // 3 jobs/sec

        // Create server types
        ServerType fastServers = new ServerType("Fast", 2);  // 2 fast servers
        ServerType slowServers = new ServerType("Slow", 3);  // 3 slow servers

        // Add server types to queue
        queue.addServerType(fastServers);
        queue.addServerType(slowServers);

        // Set scheduling policy for heterogeneous servers
        queue.setHeteroSchedPolicy(HeteroSchedPolicy.FSF);  // Fastest Server First

        // The class-level service is what JMT actually runs every pool at: its
        // loader keys the ServiceStrategy array by class, so the per-pool rates
        // below reach SolverLDES and SolverCTMC but not the JMT engine.
        queue.setService(jobClass, new Exp(1.5));

        // Set service rates per server type
        queue.setService(jobClass, fastServers, new Exp(2.0));  // Fast: rate 2.0
        queue.setService(jobClass, slowServers, new Exp(1.0));  // Slow: rate 1.0

        // Routing
        model.link(Network.serialRouting(source, queue, sink));

        System.out.println("Model configuration:");
        System.out.println("  - Arrival rate: 3.0 jobs/sec");
        System.out.println("  - Fast servers: 2 x rate 2.0 = capacity 4.0");
        System.out.println("  - Slow servers: 3 x rate 1.0 = capacity 3.0");
        System.out.println("  - Total capacity: 7.0 jobs/sec");
        System.out.println("  - Scheduling policy: FSF (Fastest Server First)\n");

        // Solve using JMT
        try {
            SolverJMT solver = new SolverJMT(model);
            solver.options.seed = 23000;
            NetworkAvgTable results = solver.getAvgTable();
            System.out.println("Results:");
            System.out.println(results);
        } catch (Exception e) {
            System.out.println("Error: " + e.getMessage());
        }
        System.out.println();
    }

    /**
     * Example 2: Heterogeneous servers with class compatibility restrictions
     *
     * Two job classes with different server compatibility:
     * - ClassA: Can be served by both Fast and Slow servers
     * - ClassB: Can only be served by Slow servers
     */
    public static void example2_compatibility() {
        System.out.println("=== Example 2: Heterogeneous Servers with Compatibility ===\n");

        Network model = new Network("HeteroCompat");

        // Create nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "HeteroQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Create job classes
        OpenClass classA = new OpenClass(model, "ClassA");
        OpenClass classB = new OpenClass(model, "ClassB");

        // Set arrival rates
        source.setArrival(classA, new Exp(1.5));  // 1.5 jobs/sec
        source.setArrival(classB, new Exp(1.0));  // 1.0 jobs/sec

        // Create server types
        ServerType fastServers = new ServerType("Fast", 2);
        ServerType slowServers = new ServerType("Slow", 3);

        // Set compatibility - Fast only serves ClassA, Slow serves both
        fastServers.addCompatibleClass(classA);
        slowServers.addCompatibleClass(classA);
        slowServers.addCompatibleClass(classB);

        // Add server types to queue
        queue.addServerType(fastServers);
        queue.addServerType(slowServers);

        // Set scheduling policy
        queue.setHeteroSchedPolicy(HeteroSchedPolicy.FSF);

        // The class-level service is what the JMT engine runs every pool at
        queue.setService(classA, new Exp(1.5));
        queue.setService(classB, new Exp(1.5));

        // Set service rates per (class, server type)
        queue.setService(classA, fastServers, new Exp(2.0));
        queue.setService(classA, slowServers, new Exp(1.0));
        queue.setService(classB, slowServers, new Exp(1.5));

        // Routing
        model.link(Network.serialRouting(source, queue, sink));

        System.out.println("Model configuration:");
        System.out.println("  - ClassA arrival: 1.5 jobs/sec");
        System.out.println("  - ClassB arrival: 1.0 jobs/sec");
        System.out.println("  - Fast servers: 2 (ClassA only, rate 2.0)");
        System.out.println("  - Slow servers: 3 (ClassA rate 1.0, ClassB rate 1.5)\n");

        try {
            SolverJMT solver = new SolverJMT(model);
            solver.options.seed = 23000;
            NetworkAvgTable results = solver.getAvgTable();
            System.out.println("Results:");
            System.out.println(results);
            System.out.println("\nExpected: ClassB has higher response time due to limited server access");
        } catch (Exception e) {
            System.out.println("Error: " + e.getMessage());
        }
        System.out.println();
    }

    /**
     * Example 3: Compare different heterogeneous scheduling policies
     */
    public static void example3_policies() {
        System.out.println("=== Example 3: Comparing Scheduling Policies ===\n");

        HeteroSchedPolicy[] policies = {
            HeteroSchedPolicy.FSF,
            HeteroSchedPolicy.ALIS,
            HeteroSchedPolicy.FAIRNESS
        };

        for (HeteroSchedPolicy policy : policies) {
            System.out.println("Policy: " + policy.toText());

            Network model = new Network("HeteroPolicy_" + policy.toText());

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "HeteroQueue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Jobs");
            source.setArrival(jobClass, new Exp(2.5));

            ServerType fast = new ServerType("Fast", 2);
            ServerType slow = new ServerType("Slow", 2);

            queue.addServerType(fast);
            queue.addServerType(slow);
            queue.setHeteroSchedPolicy(policy);

            queue.setService(jobClass, new Exp(1.5));
            queue.setService(jobClass, fast, new Exp(2.0));
            queue.setService(jobClass, slow, new Exp(1.0));

            model.link(Network.serialRouting(source, queue, sink));

            try {
                SolverJMT solver = new SolverJMT(model);
                solver.options.seed = 23000;
                NetworkAvgTable results = solver.getAvgTable();

                // Extract queue metrics
                for (int i = 0; i < results.getStationNames().size(); i++) {
                    if ("HeteroQueue".equals(results.getStationNames().get(i))) {
                        System.out.printf("  Response Time: %.4f%n", results.getRespT().get(i));
                        System.out.printf("  Utilization: %.4f%n", results.getUtil().get(i));
                        break;
                    }
                }
            } catch (Exception e) {
                System.out.println("  Error: " + e.getMessage());
            }
            System.out.println();
        }

        System.out.println("Examples completed!");
    }

    /**
     * Example 4: Job parallelism, the servers a job seizes for its whole service.
     *
     * A saturated four-server station at rate 10 clears floor(4/n) jobs at a time,
     * so its throughput drops from 40 to 20 to 10 as n goes 1, 2, 4.
     */
    public static void example4_parallelism() {
        System.out.println("=== Example 4: Job Parallelism ===\n");

        for (int n : new int[]{1, 2, 4}) {
            Network model = new Network("Parallelism" + n);

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "ParQueue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Jobs");
            source.setArrival(jobClass, new Exp(1000.0));  // saturated
            queue.setNumberOfServers(4);
            queue.setService(jobClass, new Exp(10.0));
            queue.setServerParallelism(jobClass, n);

            model.link(Network.serialRouting(source, queue, sink));

            try {
                SolverJMT solver = new SolverJMT(model);
                solver.options.seed = 23000;
                NetworkAvgTable results = solver.getAvgTable();
                for (int i = 0; i < results.getStationNames().size(); i++) {
                    if ("ParQueue".equals(results.getStationNames().get(i))) {
                        System.out.printf("  servers seized = %d, throughput = %.2f (expected %d)%n",
                                n, results.getTput().get(i), (4 / n) * 10);
                        break;
                    }
                }
            } catch (Exception e) {
                System.out.println("  Error: " + e.getMessage());
            }
        }
        System.out.println();
    }
}
