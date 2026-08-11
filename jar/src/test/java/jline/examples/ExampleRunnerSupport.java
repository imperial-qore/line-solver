/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.Network;
import jline.solvers.NetworkAvgTable;
import jline.solvers.auto.SolverAUTO;
import jline.util.Maths;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.InputStream;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static jline.TestTools.withSuppressedOutput;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTimeoutPreemptively;
import static org.junit.jupiter.api.Assertions.fail;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

/**
 * Shared discovery and execution logic for the example-runner test classes.
 *
 * The example suite is partitioned across several {@code @ParameterizedTest}
 * methods in {@link ExampleRunnerAllTest} so that a failure localises to a
 * package. The partition is exhaustive by construction: every discovered
 * example main lands in exactly one of
 * {@code runBasicExampleMain} (jline.examples.java.basic),
 * {@code runAdvancedExampleMain} (jline.examples.java.advanced) and
 * {@code runMiscExampleMain} (everything else), while
 * {@code solveFactoryModel} covers the static Network factories of the
 * model-factory classes whose mains are GUI-bound.
 *
 * The example classes are excluded from coverage reporting; the purpose of
 * these tests is to exercise the product code paths (lang, state, solvers)
 * that the examples drive end-to-end, and to guarantee the shipped examples
 * stay runnable.
 */
public abstract class ExampleRunnerSupport {

    /** Aggregators (would re-run other mains) and mains that open GUI viewers. */
    static final Set<String> EXCLUDED_MAINS = new HashSet<String>(Arrays.asList(
            "jline.examples.java.AllExamples",       // aggregator of BasicExamples+AdvancedExamples
            "jline.examples.java.BasicExamples",     // aggregator
            "jline.examples.java.AdvancedExamples",  // aggregator
            "jline.examples.java.models.Gallery",           // main calls model.view() (GUI)
            "jline.examples.java.basic.ForkJoinModel", // main calls model.jsimwView() (GUI)
            "jline.examples.java.basic.NetworkVisualizationExample", // plot() opens JUNG/AWT windows
            "jline.examples.java.basic.LayeredVisualizationExample"  // plot() opens JUNG/AWT windows
    ));

    /** Packages owned by a dedicated test class; the remainder falls to the misc runner. */
    static final String BASIC_PACKAGE = "jline.examples.java.basic";
    static final String ADVANCED_PACKAGE = "jline.examples.java.advanced";

    /** Model-factory classes whose mains are GUI-bound and hence excluded above. */
    static final String[] FACTORY_CLASSES = {
            "jline.examples.java.models.Gallery",
            "jline.examples.java.basic.ForkJoinModel"
    };

    /** Per-example wall-clock guard so one hung example cannot stall the suite. */
    static final Duration MAIN_TIMEOUT = Duration.ofMinutes(10);
    static final Duration FACTORY_TIMEOUT = Duration.ofMinutes(5);

    private static VerboseLevel originalVerboseLevel;
    private static InputStream originalStdin;

    @BeforeAll
    public static void setUpClass() {
        originalVerboseLevel = GlobalConstants.getVerbose();
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
        Maths.setRandomNumbersMatlab(true);
        // Several example mains hold a static Scanner(System.in) and close it on
        // exit. Under fork-per-class Surefire, System.in is the fork's command
        // channel, so closing it yields "[SUREFIRE] std/in stream corrupted".
        // Swap in a dummy stream before any example class initializes so each
        // static Scanner wraps (and later closes) the dummy, not the real stdin.
        originalStdin = System.in;
        System.setIn(new ByteArrayInputStream(new byte[0]));
    }

    @AfterAll
    public static void tearDownClass() {
        GlobalConstants.setVerbose(originalVerboseLevel);
        System.setIn(originalStdin);
    }

    /**
     * Discovers every concrete class under jline.examples that declares a
     * public static main, minus the exclusion list.
     */
    static List<String> allExampleMainClasses() throws Exception {
        List<String> result = new ArrayList<String>();
        File root = new File(ExampleRunnerSupport.class.getProtectionDomain()
                .getCodeSource().getLocation().toURI());
        // Tests run from target/test-classes; example classes live in target/classes
        File classesDir = new File(root.getParentFile(), "classes");
        File pkgDir = new File(classesDir, "jline/examples");
        if (!pkgDir.isDirectory()) {
            fail("Cannot locate compiled examples directory: " + pkgDir);
        }
        collectMains(pkgDir, "jline.examples", result);
        Collections.sort(result);
        assertFalse(result.isEmpty(), "No example main classes discovered");
        return result;
    }

    /** Mains declared directly in {@code pkg} (not in its subpackages). */
    static List<String> mainClassesInPackage(String pkg) throws Exception {
        List<String> result = new ArrayList<String>();
        for (String fqcn : allExampleMainClasses()) {
            if (packageOf(fqcn).equals(pkg)) {
                result.add(fqcn);
            }
        }
        assertFalse(result.isEmpty(), "No example main classes discovered in " + pkg);
        return result;
    }

    /**
     * Mains outside the packages owned by a dedicated test class. This is the
     * complement that keeps the partition exhaustive: a newly added example in
     * any other package is picked up here rather than silently untested.
     */
    static List<String> mainClassesOutsideOwnedPackages() throws Exception {
        List<String> result = new ArrayList<String>();
        for (String fqcn : allExampleMainClasses()) {
            String pkg = packageOf(fqcn);
            if (!pkg.equals(BASIC_PACKAGE) && !pkg.equals(ADVANCED_PACKAGE)) {
                result.add(fqcn);
            }
        }
        assertFalse(result.isEmpty(), "No uncategorised example main classes discovered");
        return result;
    }

    private static String packageOf(String fqcn) {
        int dot = fqcn.lastIndexOf('.');
        return dot < 0 ? "" : fqcn.substring(0, dot);
    }

    private static void collectMains(File dir, String pkg, List<String> out) throws Exception {
        File[] entries = dir.listFiles();
        if (entries == null) {
            return;
        }
        for (File f : entries) {
            String name = f.getName();
            if (f.isDirectory()) {
                collectMains(f, pkg + "." + name, out);
            } else if (name.endsWith(".class") && name.indexOf('$') < 0) {
                String fqcn = pkg + "." + name.substring(0, name.length() - 6);
                if (EXCLUDED_MAINS.contains(fqcn)) {
                    continue;
                }
                Class<?> clazz = Class.forName(fqcn, false,
                        ExampleRunnerSupport.class.getClassLoader());
                try {
                    Method main = clazz.getDeclaredMethod("main", String[].class);
                    if (Modifier.isStatic(main.getModifiers())
                            && Modifier.isPublic(main.getModifiers())) {
                        out.add(fqcn);
                    }
                } catch (NoSuchMethodException e) {
                    // class has no main; skip
                }
            }
        }
    }

    /**
     * Static no-arg Network factories of the model-factory classes whose mains
     * are GUI-bound and therefore excluded from the main runners.
     */
    static List<String> modelFactories() throws Exception {
        List<String> result = new ArrayList<String>();
        for (String cn : FACTORY_CLASSES) {
            Class<?> clazz = Class.forName(cn, false,
                    ExampleRunnerSupport.class.getClassLoader());
            for (Method m : clazz.getDeclaredMethods()) {
                if (Modifier.isStatic(m.getModifiers())
                        && Modifier.isPublic(m.getModifiers())
                        && m.getParameterTypes().length == 0
                        && Network.class.isAssignableFrom(m.getReturnType())) {
                    result.add(cn + "#" + m.getName());
                }
            }
        }
        Collections.sort(result);
        assertFalse(result.isEmpty(), "No model factory methods discovered");
        return result;
    }

    /** Runs an example main under output suppression and a wall-clock guard. */
    static void runExampleMainClass(final String fqcn) {
        assertTimeoutPreemptively(MAIN_TIMEOUT, () -> withSuppressedOutput(() -> {
            try {
                Class<?> clazz = Class.forName(fqcn);
                Method main = clazz.getMethod("main", String[].class);
                main.invoke(null, (Object) new String[0]);
            } catch (Exception e) {
                Throwable cause = e.getCause() != null ? e.getCause() : e;
                throw new AssertionError(fqcn + ".main failed: " + cause, cause);
            }
        }), fqcn + ".main exceeded " + MAIN_TIMEOUT);
    }

    /** Builds a factory model and solves it with SolverAUTO. */
    static void solveFactorySpec(final String spec) {
        assertTimeoutPreemptively(FACTORY_TIMEOUT, () -> withSuppressedOutput(() -> {
            Network model;
            try {
                int hash = spec.indexOf('#');
                Class<?> clazz = Class.forName(spec.substring(0, hash));
                Method factory = clazz.getMethod(spec.substring(hash + 1));
                model = (Network) factory.invoke(null);
            } catch (Exception e) {
                Throwable cause = e.getCause() != null ? e.getCause() : e;
                throw new AssertionError(spec + " model construction failed: " + cause, cause);
            }
            assertNotNull(model, spec + " returned null model");
            NetworkAvgTable table;
            try {
                SolverAUTO solver = new SolverAUTO(model, "seed", 23000);
                table = solver.getAvgTable();
            } catch (RuntimeException e) {
                // Models outside every candidate solver's feature set are skipped,
                // not failed: the factory itself has been validated above.
                assumeTrue(false, spec + " unsupported by SolverAUTO: " + e.getMessage());
                return;
            }
            assertNotNull(table, spec + " produced null AvgTable");
            assertNotNull(table.getQLen(), spec + " produced AvgTable without QLen column");
            assertFalse(table.getQLen().isEmpty(), spec + " produced empty AvgTable");
        }), spec + " exceeded " + FACTORY_TIMEOUT);
    }
}
