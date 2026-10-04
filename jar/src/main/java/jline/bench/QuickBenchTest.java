/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench;

import jline.bench.cqn.BenchCQN_PS;

/**
 * Quick test of benchmark package functionality
 */
public class QuickBenchTest {
    
    public static void main(String[] args) {
        System.out.println("=== Quick Benchmark Test ===");
        
        // Test one PS benchmark iteration
        System.out.println("Testing CQN PS Light Load...");
        try {
            BenchCQN_PS.runLightLoad();
            System.out.println("✅ CQN PS benchmark completed successfully!");
        } catch (Exception e) {
            System.out.println("❌ CQN PS benchmark failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Quick Benchmark Test Complete ===");
    }
}