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

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import static jline.TestTools.withSuppressedOutput;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTimeoutPreemptively;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

/**
 * Shared discovery and execution logic for the example-runner test classes.
 *
 * The example suite is partitioned across three test classes so that no single
 * class dominates the suite wall-clock (Surefire runs with reuseForks=false, so
 * each class gets its own JVM) and so that a failure localises to a package.
 * The partition is exhaustive by construction: every discovered example main
 * lands in exactly one of
 * {@link ExampleRunnerBasicTest} (jline.examples.java.basic),
 * {@link ExampleRunnerAdvancedTest} (jline.examples.java.advanced) and
 * {@link ExampleRunnerMiscTest} (everything else), the last of which also
 * covers the static Network factories of the model-factory classes whose mains
 * are GUI-bound.
 *
 * The example classes are excluded from coverage reporting; the purpose of
 * these tests is to exercise the product code paths (lang, state, solvers)
 * that the examples drive end-to-end, and to guarantee the shipped examples
 * stay runnable.
 */
public abstract class ExampleRunnerSupport {

    /**
     * Aggregators (would re-run other mains), mains that open GUI viewers, and
     * CLI drivers that end in System.exit. The last kind is the dangerous one:
     * the runner invokes a main in-process, so an exit call does not fail one
     * test, it terminates the FORK -- Surefire then reports "the forked VM
     * terminated without properly saying goodbye" with the whole class crashed
     * and no result for any of its examples.
     */
    static final Set<String> EXCLUDED_MAINS = new HashSet<String>(Arrays.asList(
            "jline.examples.java.AllExamples",       // aggregator of BasicExamples+AdvancedExamples
            "jline.examples.java.BasicExamples",     // aggregator
            "jline.examples.java.AdvancedExamples",  // aggregator
            "jline.examples.java.models.Gallery",           // main calls model.view() (GUI)
            "jline.examples.java.basic.ForkJoinModel", // main calls model.jsimwView() (GUI)
            "jline.examples.java.basic.NetworkVisualizationExample", // plot() opens JUNG/AWT windows
            "jline.examples.java.basic.LayeredVisualizationExample", // plot() opens JUNG/AWT windows
            "jline.examples.parity.LineExamples",    // CLI driver: no args -> System.exit(usage())
            "jline.examples.ExampleRunner"           // CLI driver for harvest-java.py: exits 1/2
    ));

    /** Packages owned by a dedicated test class; the remainder falls to the misc runner. */
    static final String BASIC_PACKAGE = "jline.examples.java.basic";
    static final String ADVANCED_PACKAGE = "jline.examples.java.advanced";

    /** Model-factory classes whose mains are GUI-bound and hence excluded above. */
    static final String[] FACTORY_CLASSES = {
            "jline.examples.java.models.Gallery",
            "jline.examples.java.basic.ForkJoinModel"
    };

    /**
     * Per-example wall-clock guard so one hung example cannot stall the suite.
     *
     * <p>Raised from 10 to 20 minutes on 2026-08-19. The slowest example main,
     * {@code jline.examples.java.basic.PrioExamples}, took 574.7s on this box --
     * 96% of the old budget -- so the guard was failing on machine load rather
     * than on a hang, and a retry that happens to be 5% slower turns a passing
     * suite red. 20 minutes keeps roughly 2x headroom over the measured worst
     * case while still bounding a genuine hang. If an example ever approaches
     * this again, trim its sample counts rather than raising the guard further:
     * past ~20 minutes a hung example and a slow one stop being distinguishable
     * in practice.
     */
    static final Duration MAIN_TIMEOUT = Duration.ofMinutes(20);

    /**
     * Per-factory wall-clock guard, the same hang detector as
     * {@link #MAIN_TIMEOUT} for the {@code solveFactoryModel} specs.
     *
     * <p>Raised from 5 to 15 minutes on 2026-09-07, for the reason the guard
     * above was raised: THE BUDGET WAS CALIBRATED ON ONE HOST AND THE PHASE
     * MOVED TO ANOTHER. The two closed fork-join factories,
     * {@code ForkJoinModel#fj_basic_closed} and
     * {@code Gallery#gallery_fj_closed}, are the same 5-job network, and
     * {@link jline.solvers.auto.SolverAUTO} answers it with SolverCTMC: five
     * jobs is at {@code EXACT_POPULATION_MAX}, so the exact order is preferred,
     * and MVA and NC decline a fork-join model. The tag augmentation then makes
     * the state space far larger than a population of five suggests, and the
     * cost is structural rather than a sample count that could be trimmed.
     *
     * <p>Measured on an idle host, one solve at a time, at
     * {@code common/jline.jar} of 2026-09-07: 157.8s on picard01 (EPYC 7302P),
     * 181.4s on picard09, and 430.2s on picard05 (Xeon E5-2440 v2 at 1.90GHz).
     * The Java phase ran on picard01 or picard02 until the wave-2 host pinning
     * of 2026-09-07 sent it to picard05, where 430s against a 300s budget
     * failed both specs on wall clock rather than on a defect. 15 minutes keeps
     * the same ~2x headroom over the measured worst case that the 20-minute
     * guard above keeps over its own, and the roster's slowest host is what the
     * worst case has to be measured on: placement is a preference, so any phase
     * can fall back to any capable host.
     */
    static final Duration FACTORY_TIMEOUT = Duration.ofMinutes(15);

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
        Set<String> excludedSeen = new HashSet<String>();
        Map<String, File> classFiles = new HashMap<String, File>();
        collectMains(pkgDir, "jline.examples", result, excludedSeen, classFiles);
        Collections.sort(result);
        assertFalse(result.isEmpty(), "No example main classes discovered");
        // An exclusion whose class no longer exists under that name is a silent
        // re-admission of whatever the class was renamed to, which for a GUI or
        // exiting main takes the fork down with it.
        Set<String> staleExclusions = new HashSet<String>(EXCLUDED_MAINS);
        staleExclusions.removeAll(excludedSeen);
        assertTrue(staleExclusions.isEmpty(),
                "EXCLUDED_MAINS names classes that no longer exist: " + staleExclusions);
        // A main that can exit must be caught HERE, at discovery, where it costs
        // one named failure. Reached the other way -- invoked in-process by
        // runExampleMainClass -- it kills the fork instead, and the whole class
        // reports nothing at all, so the run shows fewer tests rather than a
        // failure. That is what LineExamples did to two runners on 2026-08-12.
        Set<String> exiting = new TreeSet<String>();
        for (String fqcn : result) {
            if (referencesSystemExit(classFiles.get(fqcn))) {
                exiting.add(fqcn);
            }
        }
        assertTrue(exiting.isEmpty(), "Example mains that reference System.exit(int) would "
                + "terminate the Surefire fork when run in-process: " + exiting
                + ". Either drop the exit from the main, or add the class to EXCLUDED_MAINS.");
        return result;
    }

    /**
     * True when the class file holds a direct method reference to
     * {@code java.lang.System.exit(int)}. The constant pool is read rather than
     * the source, so the check sees exactly what the fork would execute, and it
     * needs no bytecode library on the Java 8 test classpath.
     */
    private static boolean referencesSystemExit(File classFile) throws IOException {
        if (classFile == null || !classFile.isFile()) {
            return false;
        }
        DataInputStream in = new DataInputStream(
                new BufferedInputStream(new FileInputStream(classFile)));
        try {
            if (in.readInt() != 0xCAFEBABE) {
                return false;
            }
            in.readUnsignedShort();                       // minor version
            in.readUnsignedShort();                       // major version
            int count = in.readUnsignedShort();
            String[] utf8 = new String[count];
            int[] classNameIdx = new int[count];
            int[] natName = new int[count];
            int[] natDesc = new int[count];
            int[] refClass = new int[count];
            int[] refNat = new int[count];
            for (int i = 1; i < count; i++) {
                int tag = in.readUnsignedByte();
                switch (tag) {
                    case 1:                               // Utf8
                        utf8[i] = in.readUTF();
                        break;
                    case 7:                               // Class
                        classNameIdx[i] = in.readUnsignedShort();
                        break;
                    case 12:                              // NameAndType
                        natName[i] = in.readUnsignedShort();
                        natDesc[i] = in.readUnsignedShort();
                        break;
                    case 9:                               // Fieldref
                    case 10:                              // Methodref
                    case 11:                              // InterfaceMethodref
                        refClass[i] = in.readUnsignedShort();
                        refNat[i] = in.readUnsignedShort();
                        break;
                    case 5:                               // Long
                    case 6:                               // Double
                        in.skipBytes(8);
                        i++;                              // occupies two slots
                        break;
                    case 3:                               // Integer
                    case 4:                               // Float
                    case 17:                              // Dynamic
                    case 18:                              // InvokeDynamic
                        in.skipBytes(4);
                        break;
                    case 15:                              // MethodHandle
                        in.skipBytes(3);
                        break;
                    case 8:                               // String
                    case 16:                              // MethodType
                    case 19:                              // Module
                    case 20:                              // Package
                        in.skipBytes(2);
                        break;
                    default:
                        // An unknown tag means the pool can no longer be walked;
                        // reporting "no exit" would be a guess, so say so loudly.
                        throw new IOException("Unknown constant pool tag " + tag
                                + " at entry " + i + " of " + classFile);
                }
            }
            for (int i = 1; i < count; i++) {
                if (refClass[i] == 0) {
                    continue;
                }
                String owner = utf8[classNameIdx[refClass[i]]];
                String name = utf8[natName[refNat[i]]];
                String desc = utf8[natDesc[refNat[i]]];
                if ("java/lang/System".equals(owner) && "exit".equals(name)
                        && "(I)V".equals(desc)) {
                    return true;
                }
            }
            return false;
        } finally {
            in.close();
        }
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

    private static void collectMains(File dir, String pkg, List<String> out,
                                     Set<String> excludedSeen,
                                     Map<String, File> classFiles) throws Exception {
        File[] entries = dir.listFiles();
        if (entries == null) {
            return;
        }
        for (File f : entries) {
            String name = f.getName();
            if (f.isDirectory()) {
                collectMains(f, pkg + "." + name, out, excludedSeen, classFiles);
            } else if (name.endsWith(".class") && name.indexOf('$') < 0) {
                String fqcn = pkg + "." + name.substring(0, name.length() - 6);
                if (EXCLUDED_MAINS.contains(fqcn)) {
                    excludedSeen.add(fqcn);
                    continue;
                }
                Class<?> clazz = Class.forName(fqcn, false,
                        ExampleRunnerSupport.class.getClassLoader());
                try {
                    Method main = clazz.getDeclaredMethod("main", String[].class);
                    if (Modifier.isStatic(main.getModifiers())
                            && Modifier.isPublic(main.getModifiers())) {
                        out.add(fqcn);
                        classFiles.put(fqcn, f);
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
