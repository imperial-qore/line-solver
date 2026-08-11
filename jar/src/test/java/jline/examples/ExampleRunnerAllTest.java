/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.List;

/**
 * Runs every example main in {@code jline.examples.java.basic} and
 * {@code jline.examples.java.advanced}, every example main outside those two
 * owned packages (so the partition of the example suite stays exhaustive as
 * new packages appear), and solves every static no-arg Network factory of the
 * model-factory classes (Gallery, ForkJoinModel), whose mains are GUI-bound
 * and hence not run by the example-main runners.
 */
public class ExampleRunnerAllTest extends ExampleRunnerSupport {

    public static List<String> basicExampleMains() throws Exception {
        return mainClassesInPackage(BASIC_PACKAGE);
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("basicExampleMains")
    public void runBasicExampleMain(String fqcn) {
        runExampleMainClass(fqcn);
    }

    public static List<String> advancedExampleMains() throws Exception {
        return mainClassesInPackage(ADVANCED_PACKAGE);
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("advancedExampleMains")
    public void runAdvancedExampleMain(String fqcn) {
        runExampleMainClass(fqcn);
    }

    public static List<String> miscExampleMains() throws Exception {
        return mainClassesOutsideOwnedPackages();
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("miscExampleMains")
    public void runMiscExampleMain(String fqcn) {
        runExampleMainClass(fqcn);
    }

    public static List<String> factorySpecs() throws Exception {
        return modelFactories();
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("factorySpecs")
    public void solveFactoryModel(String spec) {
        solveFactorySpec(spec);
    }
}
