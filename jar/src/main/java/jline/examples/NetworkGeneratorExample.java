/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples;

import jline.gen.NetworkGenerator;
import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodes.*;
import jline.solvers.NetworkAvgTable;
import jline.solvers.auto.AUTO;
import jline.util.matrix.Matrix;

/**
 * Example demonstrating the NetworkGenerator class
 */
public class NetworkGeneratorExample {
    
    public static void main(String[] args) {
        // Example 1: Generate a random network with default settings
        System.out.println("=== Example 1: Random Network Generation ===");
        NetworkGenerator generator = new NetworkGenerator();
        Network randomNetwork = generator.generate(3, 1, 1, 2);
        
        System.out.println("Generated network: " + randomNetwork.getName());
        System.out.println("Number of nodes: " + randomNetwork.getNumberOfNodes());
        System.out.println("Number of stations: " + randomNetwork.getNumberOfStations());
        System.out.println("Number of classes: " + randomNetwork.getNumberOfClasses());
        
        // Example 2: Generate a closed network with specific settings
        System.out.println("\n=== Example 2: Closed Network with FCFS ===");
        NetworkGenerator fcfsGen = new NetworkGenerator(
            "fcfs",              // All queues use FCFS
            "Probabilities",     // Probabilistic routing
            "Exp",               // Exponential service times
            "medium",            // Medium job load (11-20 jobs per class)
            false,               // No varying service rates
            false,               // No multi-server queues
            false,               // No random class switching
            false,               // No multi-chain class switching
            NetworkGenerator::cyclicGraph  // Cyclic topology
        );
        
        Network closedNetwork = fcfsGen.generate(4, 0, 0, 3);
        System.out.println("Generated closed network with " + 
                         closedNetwork.getNumberOfStations() + " stations and " + 
                         closedNetwork.getNumberOfClasses() + " classes");
        
        // Print station details
        for (Node node : closedNetwork.getNodes()) {
            if (node instanceof Queue) {
                Queue queue = (Queue) node;
                System.out.println("  Queue " + queue.getName() + 
                                 ": " + queue.getSchedStrategy() + 
                                 ", servers=" + queue.getNumberOfServers());
            }
        }
        
        // Example 3: Generate an open network
        System.out.println("\n=== Example 3: Open Network ===");
        NetworkGenerator openGen = new NetworkGenerator(
            "ps",                // Processor sharing
            "Random",            // Random routing
            "HyperExp",          // Hyper-exponential service times
            "low",               // Low job load
            true,                // Varying service rates
            true,                // Multi-server queues allowed
            true,                // Random class switching allowed
            false,               // Single chain per job type
            NetworkGenerator::randGraph  // Random topology
        );
        
        Network openNetwork = openGen.generate(2, 1, 2, 0);
        System.out.println("Generated open network");
        System.out.println("Has source: " + (openNetwork.getSource() != null));
        System.out.println("Has sink: " + (openNetwork.getSink() != null));
        
        // Example 4: Mixed network with class switching
        System.out.println("\n=== Example 4: Mixed Network with Class Switching ===");
        NetworkGenerator mixedGen = new NetworkGenerator();
        Network mixedNetwork = mixedGen.generate(3, 1, 1, 2);
        
        // Check for class switching nodes
        int csCount = 0;
        for (Node node : mixedNetwork.getNodes()) {
            if (node instanceof ClassSwitch) {
                csCount++;
            }
        }
        System.out.println("Class switch nodes created: " + csCount);
        
        // Example 5: Custom topology function
        System.out.println("\n=== Example 5: Custom Topology ===");
        NetworkGenerator customGen = new NetworkGenerator(
            "randomize", "randomize", "randomize", "randomize",
            true, true, false, false,
            NetworkGeneratorExample::starTopology
        );
        
        Network customNetwork = customGen.generate(4, 0, 0, 2);
        System.out.println("Generated network with star topology");
        
        // Try to solve one of the networks
        System.out.println("\n=== Solving a Generated Network ===");
        try {
            AUTO solver = new AUTO(closedNetwork);
            NetworkAvgTable results = solver.getAvgTable();
            if (results != null) {
                System.out.println("Successfully solved the generated network!");
            }
        } catch (Exception e) {
            System.out.println("Note: Solver not available or network too complex for automatic solving");
        }
    }
    
    /**
     * Custom topology function - creates a star topology
     * where node 0 is connected to all other nodes
     */
    public static Matrix starTopology(int n) {
        Matrix adj = new Matrix(n, n);
        for (int i = 1; i < n; i++) {
            adj.set(0, i, 1.0);  // From center to all others
            adj.set(i, 0, 1.0);  // From all others to center
        }
        return adj;
    }
}