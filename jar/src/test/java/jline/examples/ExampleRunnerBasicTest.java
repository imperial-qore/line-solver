/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.List;

/** Runs every example main in {@code jline.examples.java.basic}. */
public class ExampleRunnerBasicTest extends ExampleRunnerSupport {

    public static List<String> basicExampleMains() throws Exception {
        return mainClassesInPackage(BASIC_PACKAGE);
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("basicExampleMains")
    public void runBasicExampleMain(String fqcn) {
        runExampleMainClass(fqcn);
    }
}
