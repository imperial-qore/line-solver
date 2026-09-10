/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.List;

/**
 * Runs every example main outside the packages owned by
 * {@link ExampleRunnerBasicTest} and {@link ExampleRunnerAdvancedTest}, so the
 * partition of the example suite stays exhaustive as new packages appear, and
 * solves every static no-arg Network factory of the model-factory classes
 * (Gallery, ForkJoinModel), whose mains are GUI-bound and hence not run by the
 * example-main runners.
 */
public class ExampleRunnerMiscTest extends ExampleRunnerSupport {

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
