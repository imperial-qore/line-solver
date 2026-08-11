/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.testutil;

import jline.GlobalConstants;
import jline.VerboseLevel;
import org.junit.jupiter.api.extension.BeforeAllCallback;
import org.junit.jupiter.api.extension.ExtensionContext;

/**
 * Suite-wide silencer for the JAR test run.
 *
 * <p>Auto-registered as a JUnit Jupiter {@link org.junit.jupiter.api.extension.Extension}
 * through {@code META-INF/services/org.junit.jupiter.api.extension.Extension} plus
 * {@code junit.jupiter.extensions.autodetection.enabled=true} in
 * {@code junit-platform.properties}, so it applies to every test class without any
 * per-class annotation. It forces the global verbosity to
 * {@link VerboseLevel#SILENT}, which suppresses the informational solver output that
 * would otherwise clutter the test console: the per-solver
 * "analysis [method: ...] completed in ...s." profiling line (gated by
 * {@code options.verbose}, which every {@code SolverOptions} inherits from
 * {@link GlobalConstants#getVerbose()} at construction) and the {@code line_warning}
 * notices emitted by the solver runAnalyzer / refreshStruct / snNonmarkovToPh paths.</p>
 *
 * <p>This is test-scoped only; it does not alter the production default
 * ({@code VerboseLevel.STD}). Individual tests that need a different verbosity remain
 * free to set {@code GlobalConstants.Verbose} themselves.</p>
 */
public class LineVerboseSilencer implements BeforeAllCallback {

    @Override
    public void beforeAll(ExtensionContext context) {
        GlobalConstants.Verbose = VerboseLevel.SILENT;
    }
}
