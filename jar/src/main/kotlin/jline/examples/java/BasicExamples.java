/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java;

import jline.examples.java.basic.ClassSwitchExamples;
import jline.examples.java.basic.ClosedExamples;
import jline.examples.java.basic.ForkJoinExamples;
import jline.examples.java.basic.LayeredExamples;
import jline.examples.java.basic.MixedExamples;
import jline.examples.java.basic.OpenExamples;
import jline.examples.java.basic.PrioExamples;
import jline.examples.java.basic.StochPetriNetExamples;

import java.util.Scanner;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/**
 * Runs all basic examples by invoking their main methods.
 */
public class BasicExamples {
    
    private static final long TIMEOUT_SECONDS = 30; // 30 seconds timeout per example group
    private static final ExecutorService executor = Executors.newSingleThreadExecutor();
    private static final Scanner scanner = new Scanner(System.in);

    public static void main(String[] args) throws Exception {
        System.out.println("=== Running Basic Examples ===");
        System.out.println("Each example group has a " + TIMEOUT_SECONDS + " second timeout.\n");

        try {
            System.out.println("--- Running ClassSwitchExamples ---");
            runWithTimeout(() -> { ClassSwitchExamples.main(args); return null; }, "ClassSwitchExamples");
            pauseForUser();

            System.out.println("--- Running ClosedExamples ---");
            runWithTimeout(() -> { ClosedExamples.main(args); return null; }, "ClosedExamples");
            pauseForUser();

            System.out.println("--- Running ForkJoinExamples ---");
            runWithTimeout(() -> { ForkJoinExamples.main(args); return null; }, "ForkJoinExamples");
            pauseForUser();

            System.out.println("--- Running LayeredExamples ---");
            runWithTimeout(() -> { LayeredExamples.main(args); return null; }, "LayeredExamples");
            pauseForUser();

            System.out.println("--- Running MixedExamples ---");
            runWithTimeout(() -> { MixedExamples.main(args); return null; }, "MixedExamples");
            pauseForUser();

            System.out.println("--- Running OpenExamples ---");
            runWithTimeout(() -> { OpenExamples.main(args); return null; }, "OpenExamples");
            pauseForUser();

            System.out.println("--- Running PrioExamples ---");
            runWithTimeout(() -> { PrioExamples.main(args); return null; }, "PrioExamples");
            pauseForUser();

            System.out.println("--- Running StochPetriNetExamples ---");
            runWithTimeout(() -> { StochPetriNetExamples.main(args); return null; }, "StochPetriNetExamples");
            System.out.println();

            System.out.println("=== All Basic Examples Completed ===");
        } finally {
            executor.shutdown();
            try {
                if (!executor.awaitTermination(5, TimeUnit.SECONDS)) {
                    executor.shutdownNow();
                }
            } catch (InterruptedException e) {
                executor.shutdownNow();
            }
            scanner.close();
        }
    }
    
    private static void runWithTimeout(Callable<Void> task, String taskName) {
        Future<Void> future = executor.submit(task);
        try {
            future.get(TIMEOUT_SECONDS, TimeUnit.SECONDS);
        } catch (TimeoutException e) {
            future.cancel(true);
            System.err.println(taskName + " timed out after " + TIMEOUT_SECONDS + " seconds and was cancelled.");
        } catch (Exception e) {
            System.err.println(taskName + " failed with error: " + e.getMessage());
            e.printStackTrace();
        }
    }
    
    private static void pauseForUser() {
        // Skip pause if running in non-interactive mode (e.g., Maven exec)
        if (System.console() == null) {
            System.out.println("\n[Running in non-interactive mode, continuing...]");
            return;
        }
        System.out.println("\nPress Enter to continue...");
        try {
            scanner.nextLine();
        } catch (Exception e) {
            // Ignore scanner errors in case of pipe or redirection
        }
    }
}