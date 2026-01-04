/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java;

/**
 * Runs all examples by invoking BasicExamples and AdvancedExamples.
 */
public class AllExamples {

    public static void main(String[] args) throws Exception {
        System.out.println("=== Running All Examples ===\n");

        BasicExamples.main(args);
        System.out.println();

        AdvancedExamples.main(args);

        System.out.println("=== All Examples Completed ===");
    }
}