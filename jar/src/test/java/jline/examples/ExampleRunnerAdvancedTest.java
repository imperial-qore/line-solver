/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.List;

/** Runs every example main in {@code jline.examples.java.advanced}. */
public class ExampleRunnerAdvancedTest extends ExampleRunnerSupport {

    public static List<String> advancedExampleMains() throws Exception {
        return mainClassesInPackage(ADVANCED_PACKAGE);
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("advancedExampleMains")
    public void runAdvancedExampleMain(String fqcn) {
        runExampleMainClass(fqcn);
    }
}
