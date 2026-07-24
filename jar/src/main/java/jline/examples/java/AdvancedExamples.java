/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java;

import jline.examples.java.advanced.CDFRespTExamples;
import jline.examples.java.advanced.CyclicPollingExamples;
import jline.examples.java.advanced.InitStateExamples;
import jline.examples.java.advanced.LayeredCQExamples;
import jline.examples.java.advanced.LoadDependentExamples;
import jline.examples.java.advanced.RandomEnvExamples;
import jline.examples.java.advanced.StateDepRoutingExamples;
import jline.examples.java.advanced.StateProbabilitiesExamples;
import jline.examples.java.advanced.SwitchoverTimesExamples;

import java.util.Scanner;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/**
 * Runs all advanced examples by invoking their main methods.
 */
public class AdvancedExamples {
    
    private static final long TIMEOUT_SECONDS = 30; // 30 seconds timeout per example group
    private static final ExecutorService executor = Executors.newSingleThreadExecutor();
    private static final Scanner scanner = new Scanner(System.in);

    public static void main(String[] args) throws Exception {
        System.out.println("=== Running Advanced Examples ===");
        System.out.println("Each example group has a " + TIMEOUT_SECONDS + " second timeout.\n");

        try {
            System.out.println("--- Running CDFRespTExamples ---");
            runWithTimeout(() -> { CDFRespTExamples.main(args); return null; }, "CDFRespTExamples");
            System.out.println();

            System.out.println("--- Running CyclicPollingExamples ---");
            runWithTimeout(() -> { CyclicPollingExamples.main(args); return null; }, "CyclicPollingExamples");
            System.out.println();

            System.out.println("--- Running InitStateExamples ---");
            runWithTimeout(() -> { InitStateExamples.main(args); return null; }, "InitStateExamples");
            System.out.println();

            System.out.println("--- Running LayeredCQExamples ---");
            runWithTimeout(() -> { LayeredCQExamples.main(args); return null; }, "LayeredCQExamples");
            System.out.println();

            System.out.println("--- Running LoadDependentExamples ---");
            runWithTimeout(() -> { LoadDependentExamples.main(args); return null; }, "LoadDependentExamples");
            System.out.println();

            System.out.println("--- Running RandomEnvExamples ---");
            runWithTimeout(() -> { RandomEnvExamples.main(args); return null; }, "RandomEnvExamples");
            System.out.println();

            System.out.println("--- Running StateDepRoutingExamples ---");
            runWithTimeout(() -> { StateDepRoutingExamples.main(args); return null; }, "StateDepRoutingExamples");
            System.out.println();

            System.out.println("--- Running StateProbabilitiesExamples ---");
            runWithTimeout(() -> { StateProbabilitiesExamples.main(args); return null; }, "StateProbabilitiesExamples");
            System.out.println();

            System.out.println("--- Running SwitchoverTimesExamples ---");
            runWithTimeout(() -> { SwitchoverTimesExamples.main(args); return null; }, "SwitchoverTimesExamples");
            System.out.println();

            System.out.println("=== All Advanced Examples Completed ===");
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