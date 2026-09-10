/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Environment;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.Mode;
import jline.lang.OpenClass;
import jline.lang.Region;
import jline.lang.RoutingMatrix;
import jline.lang.SelfLoopingClass;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.ImpatienceType;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.layered.Activity;
import jline.lang.layered.ActivityPrecedence;
import jline.lang.layered.CacheTask;
import jline.lang.layered.Entry;
import jline.lang.layered.ItemEntry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Place;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Transition;
import jline.lang.processes.APH;
import jline.lang.processes.BMAP;
import jline.lang.processes.Cox2;
import jline.lang.processes.DMAP;
import jline.lang.processes.Distribution;
import jline.lang.processes.EmpiricalCDF;
import jline.lang.processes.Exp;
import jline.lang.processes.Expolynomial;
import jline.lang.processes.HyperExp;
import jline.lang.processes.MMDP2;
import jline.lang.processes.ME;
import jline.lang.processes.MarkedMMPP;
import jline.lang.processes.RAP;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import javax.xml.parsers.DocumentBuilderFactory;
import java.io.File;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static jline.TestTools.withSuppressedOutput;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

/**
 * Round-trip persistence tests over the Gallery model corpus:
 * - LINE JSON save/load/save reaches a fixed point (idempotent canonical form)
 *   and preserves the model topology;
 * - JMT JSIM export produces well-formed XML for every JMT-supported model;
 * - LQN XML write/parse round-trips a layered network.
 */
public class ModelIORoundTripTest {

    private static VerboseLevel originalVerboseLevel;

    @TempDir
    static Path tempDir;

    @BeforeAll
    public static void setUpClass() {
        originalVerboseLevel = GlobalConstants.getVerbose();
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
        Maths.setRandomNumbersMatlab(true);
    }

    @AfterAll
    public static void tearDownClass() {
        GlobalConstants.setVerbose(originalVerboseLevel);
    }

    /** Gallery model factories: public static no-arg methods returning Network. */
    public static List<String> galleryFactories() throws Exception {
        Class<?> gallery = Class.forName("jline.examples.java.models.Gallery");
        List<String> result = new ArrayList<String>();
        for (Method m : gallery.getDeclaredMethods()) {
            if (Modifier.isStatic(m.getModifiers())
                    && Modifier.isPublic(m.getModifiers())
                    && m.getParameterTypes().length == 0
                    && Network.class.isAssignableFrom(m.getReturnType())) {
                result.add(m.getName());
            }
        }
        Collections.sort(result);
        assertFalse(result.isEmpty(), "No Gallery factories discovered");
        return result;
    }

    private static Network buildGalleryModel(String factoryName) {
        try {
            Class<?> gallery = Class.forName("jline.examples.java.models.Gallery");
            return (Network) gallery.getMethod(factoryName).invoke(null);
        } catch (Exception e) {
            Throwable cause = e.getCause() != null ? e.getCause() : e;
            // Some gallery models depend on external assets (e.g. trace files);
            // their absence is a fixture issue, not an IO defect.
            assumeTrue(false, factoryName + " not constructible here: " + cause);
            return null;
        }
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("galleryFactories")
    public void jsonRoundTripReachesFixedPoint(String factoryName) throws Exception {
        Network original = buildGalleryModel(factoryName);
        File f1 = tempDir.resolve(factoryName + "_1.json").toFile();
        File f2 = tempDir.resolve(factoryName + "_2.json").toFile();

        withSuppressedOutput(() -> {
            try {
                LineModelIO.save(original, f1.getAbsolutePath());
                Object loaded1 = LineModelIO.load(f1.getAbsolutePath());
                assertTrue(loaded1 instanceof Network,
                        factoryName + ": loaded object is not a Network");
                Network network1 = (Network) loaded1;
                assertEquals(original.getNumberOfNodes(), network1.getNumberOfNodes(),
                        factoryName + ": node count changed in round trip");
                assertEquals(original.getNumberOfClasses(), network1.getNumberOfClasses(),
                        factoryName + ": class count changed in round trip");

                LineModelIO.save(network1, f2.getAbsolutePath());
                Object loaded2 = LineModelIO.load(f2.getAbsolutePath());
                File f3 = tempDir.resolve(factoryName + "_3.json").toFile();
                LineModelIO.save((Network) loaded2, f3.getAbsolutePath());

                // The canonical JSON form must be a fixed point of save/load
                String json2 = new String(Files.readAllBytes(f2.toPath()), StandardCharsets.UTF_8);
                String json3 = new String(Files.readAllBytes(f3.toPath()), StandardCharsets.UTF_8);
                assertEquals(json2, json3,
                        factoryName + ": JSON save/load is not idempotent");
            } catch (Exception e) {
                throw new AssertionError(factoryName + ": round trip failed: " + e, e);
            }
        });
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("galleryFactories")
    public void jsimExportIsWellFormedXml(String factoryName) throws Exception {
        Network model = buildGalleryModel(factoryName);
        SolverJMT solver;
        try {
            solver = new SolverJMT(model);
        } catch (RuntimeException e) {
            assumeTrue(false, factoryName + ": SolverJMT rejects model: " + e.getMessage());
            return;
        }
        assumeTrue(solver.supports(model),
                factoryName + ": outside SolverJMT feature set");

        File out = tempDir.resolve(factoryName + ".jsimg").toFile();
        withSuppressedOutput(() -> {
            try {
                String written = solver.writeJSIM(model.getStruct(true), out.getAbsolutePath());
                assertNotNull(written, factoryName + ": writeJSIM returned null path");
                File file = new File(written);
                assertTrue(file.isFile() && file.length() > 0,
                        factoryName + ": JSIM file missing or empty");
                // Must parse as well-formed XML
                DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
                dbf.newDocumentBuilder().parse(file);
            } catch (Exception e) {
                throw new AssertionError(factoryName + ": JSIM export failed: " + e, e);
            }
        });
    }

    /** LQN factories from the basic LayeredModel example class. */
    public static List<String> lqnFactories() throws Exception {
        Class<?> holder = Class.forName("jline.examples.java.basic.LayeredModel");
        List<String> result = new ArrayList<String>();
        for (Method m : holder.getDeclaredMethods()) {
            if (Modifier.isStatic(m.getModifiers())
                    && Modifier.isPublic(m.getModifiers())
                    && m.getParameterTypes().length == 0
                    && LayeredNetwork.class.isAssignableFrom(m.getReturnType())) {
                result.add(m.getName());
            }
        }
        Collections.sort(result);
        assertFalse(result.isEmpty(), "No LQN factories discovered");
        return result;
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("lqnFactories")
    public void lqnXmlRoundTrip(String factoryName) throws Exception {
        LayeredNetwork original;
        try {
            Class<?> holder = Class.forName("jline.examples.java.basic.LayeredModel");
            original = (LayeredNetwork) holder.getMethod(factoryName).invoke(null);
        } catch (Exception e) {
            Throwable cause = e.getCause() != null ? e.getCause() : e;
            assumeTrue(false, factoryName + " not constructible here: " + cause);
            return;
        }

        File f1 = tempDir.resolve(factoryName + ".lqnx").toFile();
        withSuppressedOutput(() -> {
            original.writeXML(f1.getAbsolutePath());
            assertTrue(f1.isFile() && f1.length() > 0,
                    factoryName + ": LQN XML missing or empty");
            LayeredNetwork reloaded = LayeredNetwork.parseXML(f1.getAbsolutePath());
            assertNotNull(reloaded, factoryName + ": parseXML returned null");
            assertEquals(original.getStruct().nhosts, reloaded.getStruct().nhosts,
                    factoryName + ": host count changed in round trip");
            assertEquals(original.getStruct().ntasks, reloaded.getStruct().ntasks,
                    factoryName + ": task count changed in round trip");
            assertEquals(original.getStruct().nentries, reloaded.getStruct().nentries,
                    factoryName + ": entry count changed in round trip");
            assertLqnFeaturesSurvive(factoryName, original, reloaded);
        });
    }

    /**
     * Asserts that the LINE .lqnx dialect carried a model's cache, item entry and
     * setup/delay-off times through the round trip.
     *
     * <p>COUNTS ALONE ARE NOT ENOUGH, and that is not a hypothetical: this test
     * asserted only nhosts, ntasks and nentries, so lqn_setup passed it for as long
     * as the writer silently dropped the setup and delay-off times entirely -- a
     * setup task survives all three count assertions while losing everything that
     * makes it a setup task. A round-trip test that compares only counts cannot
     * detect the loss it round-trips through. Every assertion below is conditional
     * on the ORIGINAL carrying the feature, so a model without one asserts nothing
     * extra and no existing case is weakened.</p>
     */
    private static void assertLqnFeaturesSurvive(String factoryName,
                                                 LayeredNetwork original,
                                                 LayeredNetwork reloaded) {
        Map<String, Task> reloadedTasks = new HashMap<>();
        for (Task t : reloaded.getTasks().values()) {
            reloadedTasks.put(t.getName(), t);
        }
        for (Task origTask : original.getTasks().values()) {
            Task backTask = reloadedTasks.get(origTask.getName());
            assertNotNull(backTask, factoryName + ": task " + origTask.getName() + " lost in round trip");

            if (origTask instanceof CacheTask) {
                assertTrue(backTask instanceof CacheTask,
                        factoryName + ": task " + origTask.getName() + " stopped being a CacheTask");
                CacheTask origCache = (CacheTask) origTask;
                CacheTask backCache = (CacheTask) backTask;
                assertEquals(origCache.getItems(), backCache.getItems(),
                        factoryName + ": cache item count changed");
                assertArrayEquals(origCache.getItemLevelCap(), backCache.getItemLevelCap(),
                        factoryName + ": cache level capacities changed (itemLevelCap is an array)");
                assertEquals(origCache.getReplacestrategy(), backCache.getReplacestrategy(),
                        factoryName + ": cache replacement strategy changed");
                assertEquals(origCache.hasRetrieval(), backCache.hasRetrieval(),
                        factoryName + ": cache retrieval flag changed");
            }

            if (origTask.hasSetupDelayoff()) {
                assertEquals(origTask.getSetupTimeMean(), backTask.getSetupTimeMean(), 1e-9,
                        factoryName + ": setup time mean changed on task " + origTask.getName());
                assertEquals(origTask.getSetupTimeSCV(), backTask.getSetupTimeSCV(), 1e-9,
                        factoryName + ": setup time SCV changed on task " + origTask.getName());
                assertEquals(origTask.getDelayOffTimeMean(), backTask.getDelayOffTimeMean(), 1e-9,
                        factoryName + ": delay-off time mean changed on task " + origTask.getName());
                assertEquals(origTask.getDelayOffTimeSCV(), backTask.getDelayOffTimeSCV(), 1e-9,
                        factoryName + ": delay-off time SCV changed on task " + origTask.getName());
            }

            Map<String, Entry> backEntries = new HashMap<>();
            for (Entry e : backTask.getEntries()) {
                backEntries.put(e.getName(), e);
            }
            for (Entry origEntry : origTask.getEntries()) {
                if (!(origEntry instanceof ItemEntry)) {
                    continue;
                }
                Entry backEntry = backEntries.get(origEntry.getName());
                assertTrue(backEntry instanceof ItemEntry,
                        factoryName + ": entry " + origEntry.getName() + " stopped being an ItemEntry");
                ItemEntry origItem = (ItemEntry) origEntry;
                ItemEntry backItem = (ItemEntry) backEntry;
                assertEquals(origItem.getCardinality(), backItem.getCardinality(),
                        factoryName + ": item cardinality changed on " + origEntry.getName());
                Distribution origPop = origItem.getPopularity();
                if (origPop != null) {
                    Distribution backPop = backItem.getPopularity();
                    assertNotNull(backPop,
                            factoryName + ": access popularity lost on " + origEntry.getName());
                    assertEquals(origPop.getName(), backPop.getName(),
                            factoryName + ": access popularity class changed on " + origEntry.getName());
                    Matrix origP = (Matrix) origPop.getParam(1).getValue();
                    Matrix backP = (Matrix) backPop.getParam(1).getValue();
                    assertEquals(origP.length(), backP.length(),
                            factoryName + ": access popularity length changed on " + origEntry.getName());
                    for (int k = 0; k < origP.length(); k++) {
                        assertEquals(origP.get(k), backP.get(k), 1e-9,
                                factoryName + ": access popularity p[" + k + "] changed on "
                                + origEntry.getName());
                    }
                }
            }

            for (ActivityPrecedence origPrec : origTask.getPrecedences()) {
                if (!ActivityPrecedenceType.POST_CACHE.equals(origPrec.getPostType())) {
                    continue;
                }
                boolean found = false;
                for (ActivityPrecedence backPrec : backTask.getPrecedences()) {
                    if (ActivityPrecedenceType.POST_CACHE.equals(backPrec.getPostType())
                            && backPrec.getPreActs().equals(origPrec.getPreActs())) {
                        assertEquals(origPrec.getPostActs(), backPrec.getPostActs(),
                                factoryName + ": POST_CACHE hit/miss assignment changed");
                        found = true;
                        break;
                    }
                }
                assertTrue(found, factoryName + ": POST_CACHE precedence lost in round trip");
            }
        }
    }

    /** JSON round trip for layered networks. */
    @ParameterizedTest(name = "{0}")
    @MethodSource("lqnFactories")
    public void lqnJsonRoundTrip(String factoryName) throws Exception {
        LayeredNetwork original;
        try {
            Class<?> holder = Class.forName("jline.examples.java.basic.LayeredModel");
            original = (LayeredNetwork) holder.getMethod(factoryName).invoke(null);
        } catch (Exception e) {
            Throwable cause = e.getCause() != null ? e.getCause() : e;
            assumeTrue(false, factoryName + " not constructible here: " + cause);
            return;
        }

        File f1 = tempDir.resolve(factoryName + "_lqn.json").toFile();
        withSuppressedOutput(() -> {
            try {
                LineModelIO.save(original, f1.getAbsolutePath());
                Object loaded = LineModelIO.load(f1.getAbsolutePath());
                assertTrue(loaded instanceof LayeredNetwork,
                        factoryName + ": loaded object is not a LayeredNetwork");
                LayeredNetwork reloaded = (LayeredNetwork) loaded;
                assertEquals(original.getStruct().ntasks, reloaded.getStruct().ntasks,
                        factoryName + ": task count changed in JSON round trip");
            } catch (Exception e) {
                throw new AssertionError(factoryName + ": LQN JSON round trip failed: " + e, e);
            }
        });
    }

    // ====================================================================
    // Targeted wire-contract coverage: constructs the Gallery corpus does
    // not exercise, each of which used to be silently dropped or downgraded.
    // ====================================================================

    /** Saves a model and loads it back through the JSON bridge. */
    private Network jsonRoundTrip(Network model, String tag) throws Exception {
        File f = tempDir.resolve(tag + "_wire.json").toFile();
        LineModelIO.save(model, f.getAbsolutePath());
        Object loaded = LineModelIO.load(f.getAbsolutePath());
        assertTrue(loaded instanceof Network, tag + ": loaded object is not a Network");
        return (Network) loaded;
    }

    /**
     * Source -> Queue -> Sink with one open class. The distribution under test
     * is placed on the Source arrival rather than the queue service, because
     * Queue.setService requires a Markovian while Source.setArrival accepts any
     * Distribution; both sides go through the same serializeDistribution /
     * deserializeDistribution pair.
     */
    private static Network openModel(String name, Distribution arv) {
        Network model = new Network(name);
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, arv);
        queue.setService(oclass, new Exp(10.0));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    private static Source sourceOf(Network model) {
        for (Node n : model.getNodes()) {
            if (n instanceof Source) {
                return (Source) n;
            }
        }
        throw new AssertionError("no Source in model");
    }

    private static Queue queueOf(Network model, String nodeName) {
        for (Node n : model.getNodes()) {
            if (n.getName().equals(nodeName)) {
                return (Queue) n;
            }
        }
        throw new AssertionError("node " + nodeName + " not found");
    }

    private static Matrix rowVector(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    /**
     * Distributions that previously had no wire branch at all: each must
     * survive as its own family, not as a moment-matched substitute.
     */
    public static List<Distribution> exoticDistributions() {
        List<Distribution> out = new ArrayList<Distribution>();
        // ME: alpha/A
        out.add(new ME(rowVector(1.0, 0.0), new Matrix(new double[][]{{-2.0, 1.0}, {0.0, -3.0}})));
        // RAP: H0/H1
        out.add(new RAP(new Matrix(new double[][]{{-2.0, 1.0}, {1.0, -3.0}}),
                new Matrix(new double[][]{{0.5, 0.5}, {1.0, 1.0}})));
        // DMAP: extends MarkovModulated, never matched instanceof MAP
        out.add(new DMAP(new Matrix(new double[][]{{0.2, 0.1}, {0.1, 0.3}}),
                new Matrix(new double[][]{{0.4, 0.3}, {0.4, 0.2}})));
        // MMDP2: scalar parameterization
        out.add(new MMDP2(2.0, 5.0, 0.5, 0.25));
        // EmpiricalCDF
        Matrix cdf = new Matrix(3, 1);
        Matrix xs = new Matrix(3, 1);
        double[] fv = {0.25, 0.5, 1.0};
        double[] xv = {1.0, 2.0, 4.0};
        for (int i = 0; i < 3; i++) {
            cdf.set(i, 0, fv[i]);
            xs.set(i, 0, xv[i]);
        }
        out.add(new EmpiricalCDF(cdf, xs));
        return out;
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("exoticDistributions")
    public void exoticDistributionSurvivesRoundTrip(Distribution dist) throws Exception {
        String tag = dist.getClass().getSimpleName();
        Network model = openModel("exotic_" + tag, dist);
        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, tag);
            } catch (Exception e) {
                throw new AssertionError(tag + ": round trip failed: " + e, e);
            }
        });
        Distribution arv = sourceOf(out[0]).getArrivalDistribution(out[0].getClasses().get(0));
        assertNotNull(arv, tag + ": arrival lost in round trip");
        assertEquals(dist.getClass(), arv.getClass(),
                tag + ": distribution family changed in round trip");
        assertEquals(dist.getMean(), arv.getMean(), 1e-9 * Math.abs(dist.getMean()) + 1e-9,
                tag + ": mean changed in round trip");
    }

    /**
     * Cox2 is a schema-legal type name. The writer maps it to "Coxian" (as
     * linemodel_save.m does: case {'Coxian','Cox2'} -> 'Coxian'), but the reader
     * had no Cox2 case, so a file declaring one reached the terminal fallback
     * and became a zero-service Immediate.
     */
    @Test
    public void cox2TypeIsNotDecodedAsImmediate() throws Exception {
        Cox2 oracle = new Cox2(3.0, 6.0, 0.4);
        File f = tempDir.resolve("cox2_type.json").toFile();
        String json = "{\"format\":\"line-model\",\"version\":\"1.0\",\"model\":{"
                + "\"type\":\"Network\",\"name\":\"cox2\","
                + "\"nodes\":[{\"name\":\"Source\",\"type\":\"Source\",\"service\":{"
                + "\"Class1\":{\"type\":\"Cox2\",\"params\":{\"mu1\":3.0,\"mu2\":6.0,"
                + "\"phi1\":0.4}}}},"
                + "{\"name\":\"Queue\",\"type\":\"Queue\",\"scheduling\":\"FCFS\",\"service\":{"
                + "\"Class1\":{\"type\":\"Exp\",\"params\":{\"lambda\":10.0}}}},"
                + "{\"name\":\"Sink\",\"type\":\"Sink\"}],"
                + "\"classes\":[{\"name\":\"Class1\",\"type\":\"Open\"}],"
                + "\"routing\":{\"type\":\"matrix\",\"matrix\":{\"Class1,Class1\":{"
                + "\"Source\":{\"Queue\":1.0},\"Queue\":{\"Sink\":1.0}}}}}}";
        Files.write(f.toPath(), json.getBytes(StandardCharsets.UTF_8));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = (Network) LineModelIO.load(f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("cox2: load failed: " + e, e);
            }
        });
        Distribution arv = sourceOf(out[0]).getArrivalDistribution(out[0].getClasses().get(0));
        assertFalse(arv instanceof jline.lang.processes.Immediate,
                "Cox2 decoded as a zero-service Immediate");
        assertEquals(oracle.getMean(), arv.getMean(), 1e-9, "Cox2 mean lost on load");
        assertEquals(oracle.getSCV(), arv.getSCV(), 1e-9, "Cox2 SCV lost on load");
    }

    /**
     * BMAP is excluded from the MMAP branch and had no branch of its own, so
     * it fell through to the mean/SCV fallback.
     */
    @Test
    public void bmapSurvivesRoundTrip() throws Exception {
        Matrix D0 = new Matrix(new double[][]{{-3.0, 0.5}, {0.5, -4.0}});
        Matrix D1 = new Matrix(new double[][]{{1.5, 0.5}, {1.0, 1.0}});
        Matrix D2 = new Matrix(new double[][]{{0.3, 0.2}, {1.0, 0.5}});
        BMAP bmap = new BMAP(D0, D1, D2);
        Network model = openModel("bmap", bmap);
        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "bmap");
            } catch (Exception e) {
                throw new AssertionError("bmap: round trip failed: " + e, e);
            }
        });
        Distribution arv = sourceOf(out[0]).getArrivalDistribution(out[0].getClasses().get(0));
        assertTrue(arv instanceof BMAP, "BMAP decoded as " + arv.getClass().getSimpleName());
        BMAP back = (BMAP) arv;
        assertEquals(2, back.getMaxBatchSize(), "BMAP batch size changed in round trip");
        assertEquals(D2.get(1, 0), back.getBatchMatrix(2).get(1, 0), 1e-12,
                "BMAP batch-2 matrix changed in round trip");
    }

    /**
     * The BMAP block array must sit under "params", which is where
     * linemodel_save.m and linemodel_io.py put it and where linemodel_io.py
     * reads it. This class used to write AND read it at the TOP LEVEL, so the
     * JAR-to-JAR round trip above passed while every cross-codebase load failed:
     * a MATLAB or Python model.json carrying a BMAP service was rejected with
     * "BMAP requires D0 and at least one batch matrix", which made LDES's BMAP
     * service (BMSP) unreachable from either. The top-level form is still
     * accepted on load so files written before the fix keep working.
     */
    @Test
    public void bmapWireFormatNestsTheBlocksUnderParams() throws Exception {
        Matrix D0 = new Matrix(new double[][]{{-3.0, 0.5}, {0.5, -4.0}});
        Matrix D1 = new Matrix(new double[][]{{1.5, 0.5}, {1.0, 1.0}});
        Matrix D2 = new Matrix(new double[][]{{0.3, 0.2}, {1.0, 0.5}});
        Network model = openModel("bmapwire", new BMAP(D0, D1, D2));
        File f = tempDir.resolve("bmapwire_wire.json").toFile();
        withSuppressedOutput(() -> {
            try {
                LineModelIO.save(model, f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("bmapwire: save failed: " + e, e);
            }
        });
        String json = new String(java.nio.file.Files.readAllBytes(f.toPath()), "UTF-8");
        int typeAt = json.indexOf("\"BMAP\"");
        assertTrue(typeAt >= 0, "no BMAP distribution in the emitted JSON");
        // the "params" wrapper must appear in the same distribution object
        String around = json.substring(Math.max(0, typeAt - 400),
                Math.min(json.length(), typeAt + 400));
        assertTrue(around.contains("params"),
                "BMAP blocks are not nested under \"params\": " + around);

        // and a document in exactly that shape must load
        Object loaded = LineModelIO.load(f.getAbsolutePath());
        assertTrue(loaded instanceof Network, "params-shaped BMAP failed to load");
        Distribution arv = sourceOf((Network) loaded)
                .getArrivalDistribution(((Network) loaded).getClasses().get(0));
        assertTrue(arv instanceof BMAP, "params-shaped BMAP decoded as "
                + arv.getClass().getSimpleName());
        assertEquals(D2.get(1, 0), ((BMAP) arv).getBatchMatrix(2).get(1, 0), 1e-12,
                "params-shaped BMAP lost its batch-2 matrix");
    }

    /**
     * MarkedMMPP extends Marked, not MarkedMAP, so it matched neither the MMAP
     * nor the MAP branch.
     */
    @Test
    public void markedMmppSurvivesRoundTrip() throws Exception {
        Matrix D0 = new Matrix(new double[][]{{-3.0, 0.5}, {0.5, -4.0}});
        Matrix D1 = new Matrix(new double[][]{{2.5, 0.0}, {0.0, 3.5}});
        Matrix D11 = new Matrix(new double[][]{{1.0, 0.0}, {0.0, 1.5}});
        Matrix D12 = new Matrix(new double[][]{{1.5, 0.0}, {0.0, 2.0}});
        MatrixCell cell = new MatrixCell();
        cell.set(0, D0);
        cell.set(1, D1);
        cell.set(2, D11);
        cell.set(3, D12);
        MarkedMMPP mmpp = new MarkedMMPP(cell);
        Network model = openModel("markedmmpp", mmpp);
        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "markedmmpp");
            } catch (Exception e) {
                throw new AssertionError("markedmmpp: round trip failed: " + e, e);
            }
        });
        Distribution arv = sourceOf(out[0]).getArrivalDistribution(out[0].getClasses().get(0));
        assertTrue(arv instanceof MarkedMMPP,
                "MarkedMMPP decoded as " + arv.getClass().getSimpleName());
        MarkedMMPP back = (MarkedMMPP) arv;
        assertEquals(4, back.getProcess().size(), "MarkedMMPP process cell changed in round trip");
        assertEquals(D12.get(1, 1), back.getProcess().get(3).get(1, 1), 1e-12,
                "MarkedMMPP per-mark matrix changed in round trip");
    }

    /**
     * d4: the Delay branch used to shadow the Queue branch, so a Delay never
     * emitted any of the station-level fields declared on Station/Queue.
     */
    @Test
    public void delayCarriesStationFields() throws Exception {
        Network model = new Network("delayFields");
        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.1));
        delay.setService(oclass, new Exp(1.0));
        delay.setCapacity(7);
        delay.setClassCap(oclass, 4);
        delay.setDropRule(oclass, DropStrategy.Drop);
        // Queue.setLoadDependence rejects INF scheduling at the lang level, so a
        // load-dependent Delay is unreachable; the class-dependence handle is the
        // scaling a Delay can actually carry.
        // The peak is MANDATORY and, on an open class, genuinely not derivable:
        // beta(n) = 1 + n grows without bound, so the declaration is the only
        // thing that says what utilization is a fraction OF. 11 is beta at the
        // open-class wire cutoff of 10.
        delay.setLimitedClassDependence((Matrix n) -> rowVector(1.0 + n.get(0)), rowVector(11.0));
        delay.setPatience(oclass, ImpatienceType.RENEGING, new Exp(0.25));
        model.link(Network.serialRouting(source, delay, sink));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "delayFields");
            } catch (Exception e) {
                throw new AssertionError("delayFields: round trip failed: " + e, e);
            }
        });
        Queue back = queueOf(out[0], "Delay");
        JobClass jc = out[0].getClasses().get(0);
        assertTrue(back instanceof Delay, "Delay decoded as " + back.getClass().getSimpleName());
        assertEquals(7.0, back.getCap(), 1e-12, "Delay buffer lost in round trip");
        assertEquals(4.0, back.getClassCap(jc), 1e-12, "Delay classCap lost in round trip");
        assertEquals(DropStrategy.Drop, back.getDropRule(jc), "Delay dropRule lost in round trip");
        assertNotNull(back.getLimitedClassDependence(),
                "Delay classDependence lost in round trip");
        assertEquals(3.0, back.getLimitedClassDependence().apply(rowVector(2.0)).get(0), 1e-12,
                "Delay classDependence scaling changed in round trip");
        assertTrue(back.hasPatience(jc), "Delay patience lost in round trip");
        assertEquals(4.0, back.getPatience(jc).getMean(), 1e-9,
                "Delay patience distribution changed in round trip");
    }

    /** d5: SQ (formerly KCHOICES) used to reset to the OutputStrategy default d=2. */
    @Test
    public void sqParametersSurviveRoundTrip() throws Exception {
        Network model = new Network("sq");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.1));
        q1.setService(oclass, new Exp(1.0));
        q2.setService(oclass, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, source, q1, 0.5);
        P.set(oclass, oclass, source, q2, 0.5);
        P.set(oclass, oclass, q1, sink, 1.0);
        P.set(oclass, oclass, q2, sink, 1.0);
        model.link(P);
        source.setSQRouting(oclass, 5);

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "sq");
            } catch (Exception e) {
                throw new AssertionError("sq: round trip failed: " + e, e);
            }
        });
        Node backSource = out[0].getNodes().get(0);
        JobClass jc = out[0].getClasses().get(0);
        assertEquals(RoutingStrategy.SQ, backSource.getRoutingStrategy(jc),
                "SQ strategy lost in round trip");
        boolean found = false;
        for (jline.lang.OutputStrategy os : backSource.getOutput().getOutputStrategyByClass(jc)) {
            if (os.getRoutingStrategy() == RoutingStrategy.SQ) {
                assertEquals(5, os.getSqD(), "SQ d silently reset to the default");
                found = true;
            }
        }
        assertTrue(found, "no SQ output strategy after round trip");
    }

    /** d10: SelfLoopingClass subclasses ClosedClass and must be tested first. */
    @Test
    public void selfLoopingClassSurvivesRoundTrip() throws Exception {
        Network model = new Network("selfloop");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        ClosedClass cclass = new ClosedClass(model, "Closed1", 2, delay);
        SelfLoopingClass slc = new SelfLoopingClass(model, "SelfLoop1", 1, queue);
        delay.setService(cclass, new Exp(1.0));
        queue.setService(cclass, new Exp(2.0));
        queue.setService(slc, new Exp(3.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(cclass, cclass, delay, queue, 1.0);
        P.set(cclass, cclass, queue, delay, 1.0);
        P.set(slc, slc, queue, queue, 1.0);
        model.link(P);

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "selfloop");
            } catch (Exception e) {
                throw new AssertionError("selfloop: round trip failed: " + e, e);
            }
        });
        JobClass back = null;
        for (JobClass jc : out[0].getClasses()) {
            if ("SelfLoop1".equals(jc.getName())) {
                back = jc;
            }
        }
        assertNotNull(back, "SelfLooping class lost in round trip");
        assertTrue(back instanceof SelfLoopingClass,
                "SelfLoopingClass downgraded to " + back.getClass().getSimpleName());
    }

    /** FCR linear constraints A*n <= b: implemented in Region, never bridged. */
    @Test
    public void fcrLinearConstraintsSurviveRoundTrip() throws Exception {
        Network model = new Network("fcrlincon");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.1));
        q1.setService(oclass, new Exp(1.0));
        model.link(Network.serialRouting(source, q1, sink));
        Region region = model.addRegion(Arrays.asList((Node) q1));
        region.setGlobalMaxJobs(6);
        Matrix A = new Matrix(new double[][]{{2.0}});
        Matrix b = new Matrix(1, 1);
        b.set(0, 0, 5.0);
        region.setLinearConstraints(A, b);

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "fcrlincon");
            } catch (Exception e) {
                throw new AssertionError("fcrlincon: round trip failed: " + e, e);
            }
        });
        assertEquals(1, out[0].getRegions().size(), "region lost in round trip");
        Region backRegion = out[0].getRegions().get(0);
        assertTrue(backRegion.hasLinearConstraints(), "FCR linear constraints lost in round trip");
        assertEquals(2.0, backRegion.getLinearConstraints()[0].get(0, 0), 1e-12,
                "FCR constraint matrix A changed in round trip");
        assertEquals(5.0, backRegion.getLinearConstraints()[1].get(0, 0), 1e-12,
                "FCR constraint vector b changed in round trip");
    }

    /** Class-level state that reached no wire key at all. */
    @Test
    public void classLevelFieldsSurviveRoundTrip() throws Exception {
        Network model = openModel("classfields", new Exp(1.0));
        OpenClass oclass = (OpenClass) model.getClasses().get(0);
        oclass.setReferenceClass(true);
        oclass.setPatience(ImpatienceType.RENEGING, new Exp(0.2));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "classfields");
            } catch (Exception e) {
                throw new AssertionError("classfields: round trip failed: " + e, e);
            }
        });
        JobClass back = out[0].getClasses().get(0);
        assertTrue(back.isReferenceClass(), "isReferenceClass lost in round trip");
        assertTrue(back.hasPatience(), "class-level patience lost in round trip");
        assertEquals(5.0, back.getPatience().getMean(), 1e-9,
                "class-level patience distribution changed in round trip");
        assertEquals(ImpatienceType.RENEGING, back.getImpatienceType(),
                "class-level impatienceType lost in round trip");
    }

    /**
     * replySignalClass carries sn.syncreply; without it a REPLY signal is inert
     * after a round trip.
     */
    @Test
    public void replySignalClassSurvivesRoundTrip() throws Exception {
        Network model = new Network("replysignal");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass work = new OpenClass(model, "Work");
        OpenClass reply = new OpenClass(model, "Reply");
        source.setArrival(work, new Exp(0.1));
        source.setArrival(reply, new Exp(0.1));
        queue.setService(work, new Exp(1.0));
        queue.setService(reply, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(work, work, source, queue, 1.0);
        P.set(work, work, queue, sink, 1.0);
        P.set(reply, reply, source, queue, 1.0);
        P.set(reply, reply, queue, sink, 1.0);
        model.link(P);
        // 1-based, as JobClass.setReplySignalClassIndex documents
        work.setReplySignalClassIndex(model.getClasses().indexOf(reply) + 1);

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "replysignal");
            } catch (Exception e) {
                throw new AssertionError("replysignal: round trip failed: " + e, e);
            }
        });
        JobClass backWork = null;
        JobClass backReply = null;
        for (JobClass jc : out[0].getClasses()) {
            if ("Work".equals(jc.getName())) backWork = jc;
            if ("Reply".equals(jc.getName())) backReply = jc;
        }
        assertNotNull(backWork, "Work class lost in round trip");
        assertTrue(backWork.expectsReply(), "replySignalClass lost in round trip");
        assertEquals(out[0].getClasses().indexOf(backReply) + 1,
                backWork.getReplySignalClassIndex(),
                "replySignalClass points at the wrong class after a round trip");
    }

    /** Station-level orbit impatience and immediate feedback. */
    @Test
    public void orbitImpatienceAndImmediateFeedbackSurviveRoundTrip() throws Exception {
        Network model = openModel("orbitimm", new Exp(1.0));
        Queue queue = queueOf(model, "Queue");
        JobClass jc = model.getClasses().get(0);
        queue.setOrbitImpatience(jc, new Exp(0.5));
        queue.setImmediateFeedback(jc);

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "orbitimm");
            } catch (Exception e) {
                throw new AssertionError("orbitimm: round trip failed: " + e, e);
            }
        });
        Queue back = queueOf(out[0], "Queue");
        JobClass backJc = out[0].getClasses().get(0);
        assertTrue(back.hasOrbitImpatience(backJc), "orbitImpatience lost in round trip");
        assertEquals(2.0, back.getOrbitImpatience(backJc).getMean(), 1e-9,
                "orbitImpatience distribution changed in round trip");
        assertTrue(back.hasImmediateFeedback(backJc.getIndex()),
                "immediateFeedback lost in round trip");
    }

    /**
     * A declared state is a TRIO -- state, state space and prior -- and the wire
     * carries only the first when the prior is the trivial [1] over one row, so
     * the READER owes the other two back. A node holding a state over an empty
     * space is not one any solver can start from: SolverFluid indexes the space
     * directly, and MATLAB's reader leaving it empty returned an all-zero table
     * for every stage of a reloaded random environment. The model is initialized
     * FIRST because a state on a strict subset of the stateful nodes does not
     * travel at all -- see serializeNetworkNodes and linemodel_save.m.
     */
    @Test
    public void aDeclaredStateRestoresSpaceAndPrior() throws Exception {
        Network model = openModel("statetrio", new Exp(1.0));
        model.initDefault();
        Queue queue = queueOf(model, "Queue");
        queue.setState(rowVector(2.0, 1.0));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "statetrio");
            } catch (Exception e) {
                throw new AssertionError("statetrio: round trip failed: " + e, e);
            }
        });
        Queue back = queueOf(out[0], "Queue");
        assertEquals(2.0, back.getState().get(0, 0), 1e-12, "state lost in round trip");
        assertEquals(1.0, back.getState().get(0, 1), 1e-12, "state lost in round trip");
        assertEquals(1, back.getStateSpace().getNumRows(),
                "a declared state must come back as a one-row state space");
        assertEquals(2.0, back.getStateSpace().get(0, 0), 1e-12,
                "the restored state space is not the declared row");
        assertEquals(1, back.getStatePrior().getNumRows(),
                "a declared state must come back under a prior of one");
        assertEquals(1.0, back.getStatePrior().get(0, 0), 1e-12,
                "the restored prior is not [1]");
    }

    /**
     * A NON-trivial prior does ride the wire, paired with the state space whose
     * rows it indexes, and must not be overwritten by the trivial pair above.
     */
    @Test
    public void anExplicitStatePriorSurvivesRoundTrip() throws Exception {
        Network model = openModel("statepr", new Exp(1.0));
        model.initDefault();
        Queue queue = queueOf(model, "Queue");
        Matrix space = new Matrix(2, 2);
        space.set(0, 0, 1.0); space.set(0, 1, 1.0);
        space.set(1, 0, 0.0); space.set(1, 1, 1.0);
        queue.setStateSpace(space);
        Matrix prior = new Matrix(2, 1);
        prior.set(0, 0, 0.25);
        prior.set(1, 0, 0.75);
        queue.setStatePrior(prior);
        queue.setState(rowVector(1.0, 1.0));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "statepr");
            } catch (Exception e) {
                throw new AssertionError("statepr: round trip failed: " + e, e);
            }
        });
        Queue back = queueOf(out[0], "Queue");
        assertEquals(2, back.getStateSpace().getNumRows(), "state space rows lost in round trip");
        assertEquals(2, back.getStatePrior().getNumRows(), "prior rows lost in round trip");
        assertEquals(0.25, back.getStatePrior().get(0, 0), 1e-12, "prior changed in round trip");
        assertEquals(0.75, back.getStatePrior().get(1, 0), 1e-12, "prior changed in round trip");
    }

    /**
     * The stage TYPE of a random environment. Every reader looked for it and no
     * writer emitted it, so an Environment came back from JSON with its stage
     * types blanked -- which getStageTable prints and every UP/DOWN consumer
     * keys off.
     */
    @Test
    public void environmentStageTypeSurvivesRoundTrip() throws Exception {
        Environment env = new Environment("stagetype", 2);
        env.addStage(0, "Fast", "operational", openModel("fast", new Exp(1.0)));
        env.addStage(1, "Slow", "degraded", openModel("slow", new Exp(1.0)));
        env.addTransition(0, 1, new Exp(0.5));
        env.addTransition(1, 0, new Exp(1.0));

        File f = tempDir.resolve("stagetype_wire.json").toFile();
        final Object[] out = new Object[1];
        withSuppressedOutput(() -> {
            try {
                LineModelIO.save(env, f.getAbsolutePath());
                out[0] = LineModelIO.load(f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("stagetype: round trip failed: " + e, e);
            }
        });
        assertTrue(out[0] instanceof Environment, "loaded object is not an Environment");
        Environment back = (Environment) out[0];
        assertEquals("operational", back.getStageType(0), "stage type lost in round trip");
        assertEquals("degraded", back.getStageType(1), "stage type lost in round trip");
    }

    /**
     * Server breakdown/repair, with the degraded down-server service. Without
     * the wire block the LDES clients, which serialize the model, would simulate
     * an always-up server and return a result indistinguishable from a correct one.
     */
    @Test
    public void breakdownSurvivesRoundTrip() throws Exception {
        Network model = openModel("breakdown", new Exp(1.0));
        Queue queue = queueOf(model, "Queue");
        JobClass jc = model.getClasses().get(0);
        queue.setBreakdown(new Exp(0.2), new Exp(1.0));
        queue.setDownService(jc, new Exp(0.5));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "breakdown");
            } catch (Exception e) {
                throw new AssertionError("breakdown: round trip failed: " + e, e);
            }
        });
        Queue back = queueOf(out[0], "Queue");
        JobClass backJc = out[0].getClasses().get(0);
        assertTrue(back.hasBreakdown(), "breakdown lost in round trip");
        assertEquals(5.0, back.getBreakdownFailure().getMean(), 1e-9,
                "time to failure changed in round trip");
        assertEquals(1.0, back.getBreakdownRepair().getMean(), 1e-9,
                "repair time changed in round trip");
        assertNotNull(back.getDownService(backJc), "downService lost in round trip");
        assertEquals(2.0, back.getDownService(backJc).getMean(), 1e-9,
                "downService distribution changed in round trip");
    }

    /**
     * An n-phase hyper-exponential on the wire must not be truncated to
     * (p[0], lambda[0], lambda[1]), and must load back as a HyperExp of n phases.
     */
    @Test
    public void nPhaseHyperExpIsNotTruncated() throws Exception {
        File f = tempDir.resolve("nphase_hyperexp.json").toFile();
        String json = "{\"format\":\"line-model\",\"version\":\"1.0\",\"model\":{"
                + "\"type\":\"Network\",\"name\":\"nphase\","
                + "\"nodes\":[{\"name\":\"Source\",\"type\":\"Source\",\"service\":{"
                + "\"Class1\":{\"type\":\"Exp\",\"params\":{\"lambda\":0.1}}}},"
                + "{\"name\":\"Queue\",\"type\":\"Queue\",\"scheduling\":\"FCFS\",\"service\":{"
                + "\"Class1\":{\"type\":\"HyperExp\",\"params\":{"
                + "\"p\":[0.2,0.3,0.5],\"lambda\":[1.0,2.0,8.0]}}}},"
                + "{\"name\":\"Sink\",\"type\":\"Sink\"}],"
                + "\"classes\":[{\"name\":\"Class1\",\"type\":\"Open\"}],"
                + "\"routing\":{\"type\":\"matrix\",\"matrix\":{\"Class1,Class1\":{"
                + "\"Source\":{\"Queue\":1.0},\"Queue\":{\"Sink\":1.0}}}}}}";
        Files.write(f.toPath(), json.getBytes(StandardCharsets.UTF_8));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = (Network) LineModelIO.load(f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("nphase: load failed: " + e, e);
            }
        });
        Distribution svc = queueOf(out[0], "Queue").getService(out[0].getClasses().get(0));
        // Exact mean of the declared 3-phase hyper-exponential
        double expectedMean = 0.2 / 1.0 + 0.3 / 2.0 + 0.5 / 8.0;
        assertEquals(expectedMean, svc.getMean(), 1e-9,
                "n-phase HyperExp truncated to its first two phases");
        double m2 = 2 * (0.2 / 1.0 + 0.3 / 4.0 + 0.5 / 64.0);
        double expectedScv = (m2 - expectedMean * expectedMean) / (expectedMean * expectedMean);
        assertEquals(expectedScv, svc.getSCV(), 1e-8,
                "n-phase HyperExp second moment changed on load");
        assertTrue(svc instanceof HyperExp,
                "n-phase HyperExp loaded back as " + svc.getClass().getSimpleName()
                        + " rather than HyperExp");
        assertEquals("HyperExp", svc.getName(), "n-phase HyperExp lost its type name on load");
        assertEquals(3L, ((HyperExp) svc).getNumberOfPhases(),
                "n-phase HyperExp lost phases on load");
        assertArrayEquals(new double[]{0.2, 0.3, 0.5}, ((HyperExp) svc).getP(), 1e-12,
                "n-phase HyperExp branch probabilities changed on load");
        assertArrayEquals(new double[]{1.0, 2.0, 8.0}, ((HyperExp) svc).getLambda(), 1e-12,
                "n-phase HyperExp rates changed on load");
    }

    /**
     * A 3-phase HyperExp built in Java must survive save/load unchanged: same
     * type, same phase count, same moments. The wire form carries p and lambda
     * both of length n, never of length 2n.
     */
    @Test
    public void nPhaseHyperExpSurvivesJsonRoundTrip() throws Exception {
        double[] p = new double[]{0.3, 0.4, 0.3};
        double[] lambda = new double[]{3.0, 1.0, 0.2};
        HyperExp he = new HyperExp(p, lambda);
        // Exact n-phase moments: E[X] = sum p_i/mu_i, E[X^2] = sum 2 p_i/mu_i^2
        double expectedMean = 0.3 / 3.0 + 0.4 / 1.0 + 0.3 / 0.2;
        double m2 = 2 * (0.3 / 9.0 + 0.4 / 1.0 + 0.3 / 0.04);
        double expectedScv = m2 / (expectedMean * expectedMean) - 1;
        assertEquals(expectedMean, he.getMean(), 1e-12, "3-phase HyperExp mean is wrong");
        assertEquals(expectedScv, he.getSCV(), 1e-12, "3-phase HyperExp SCV is wrong");
        assertEquals(3L, he.getNumberOfPhases(), "3-phase HyperExp reports the wrong phase count");

        Network model = new Network("nphase_rt");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.1));
        queue.setService(oclass, he);
        model.link(Network.serialRouting(source, queue, sink));

        File f = tempDir.resolve("nphase_rt_wire.json").toFile();
        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                LineModelIO.save(model, f.getAbsolutePath());
                out[0] = (Network) LineModelIO.load(f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("nphase_rt: round trip failed: " + e, e);
            }
        });

        // The written vectors are of length n, not 2n
        String wire = new String(Files.readAllBytes(f.toPath()), StandardCharsets.UTF_8);
        assertTrue(wire.contains("\"type\": \"HyperExp\"") || wire.contains("\"type\":\"HyperExp\""),
                "3-phase HyperExp not written as type HyperExp");

        Distribution svc = queueOf(out[0], "Queue").getService(out[0].getClasses().get(0));
        assertTrue(svc instanceof HyperExp,
                "3-phase HyperExp loaded back as " + svc.getClass().getSimpleName());
        HyperExp back = (HyperExp) svc;
        assertEquals(3L, back.getNumberOfPhases(), "3-phase HyperExp lost phases in round trip");
        assertArrayEquals(p, back.getP(), 1e-12, "branch probabilities changed in round trip");
        assertArrayEquals(lambda, back.getLambda(), 1e-12, "rates changed in round trip");
        assertEquals(expectedMean, back.getMean(), 1e-12, "3-phase mean changed in round trip");
        assertEquals(expectedScv, back.getSCV(), 1e-12, "3-phase SCV changed in round trip");
    }

    /**
     * Regression: a 2-phase HyperExp must be numerically and structurally
     * unchanged, whichever of the three constructors built it.
     */
    @Test
    public void twoPhaseHyperExpIsUnchanged() throws Exception {
        HyperExp legacy = new HyperExp(0.25, 2.0, 5.0);
        assertEquals(2L, legacy.getNumberOfPhases(), "2-phase HyperExp phase count changed");
        double expectedMean = 0.25 / 2.0 + 0.75 / 5.0;
        double m2 = 2 * (0.25 / 4.0 + 0.75 / 25.0);
        double expectedScv = m2 / (expectedMean * expectedMean) - 1;
        assertEquals(expectedMean, legacy.getMean(), 1e-12, "2-phase HyperExp mean changed");
        assertEquals(expectedScv, legacy.getSCV(), 1e-12, "2-phase HyperExp SCV changed");

        // The n-phase constructor must reproduce the 2-phase one exactly
        HyperExp viaVector = new HyperExp(new double[]{0.25, 0.75}, new double[]{2.0, 5.0});
        assertEquals(legacy.getNumberOfPhases(), viaVector.getNumberOfPhases(),
                "n-phase constructor disagrees on 2-phase phase count");
        assertEquals(legacy.getMean(), viaVector.getMean(), 0.0,
                "n-phase constructor disagrees on 2-phase mean");
        assertEquals(legacy.getSCV(), viaVector.getSCV(), 0.0,
                "n-phase constructor disagrees on 2-phase SCV");
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(legacy.D(0).get(i, j), viaVector.D(0).get(i, j), 0.0,
                        "n-phase constructor disagrees on 2-phase D0(" + i + "," + j + ")");
                assertEquals(legacy.D(1).get(i, j), viaVector.D(1).get(i, j), 0.0,
                        "n-phase constructor disagrees on 2-phase D1(" + i + "," + j + ")");
            }
        }

        // The single-rate constructor is a degenerate exponential
        HyperExp single = new HyperExp(0.4, 3.0);
        assertEquals(2L, single.getNumberOfPhases(), "HyperExp(p,lambda) phase count changed");
        assertEquals(1.0 / 3.0, single.getMean(), 1e-12, "HyperExp(p,lambda) mean changed");

        // And it round-trips through JSON unchanged
        Network model = openModel("hyperexp2", legacy);
        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "hyperexp2");
            } catch (Exception e) {
                throw new AssertionError("hyperexp2: round trip failed: " + e, e);
            }
        });
        Distribution arv = sourceOf(out[0]).getArrivalDistribution(out[0].getClasses().get(0));
        assertTrue(arv instanceof HyperExp,
                "2-phase HyperExp loaded back as " + arv.getClass().getSimpleName());
        assertEquals(2L, ((HyperExp) arv).getNumberOfPhases(),
                "2-phase HyperExp phase count changed in round trip");
        assertEquals(expectedMean, arv.getMean(), 1e-12,
                "2-phase HyperExp mean changed in round trip");
        assertEquals(expectedScv, arv.getSCV(), 1e-12,
                "2-phase HyperExp SCV changed in round trip");
    }

    /**
     * JMT's jmt.engine.random.HyperExpPar is 2-phase only, so an n-phase
     * HyperExp service must be exported as a PhaseTypeDistr/PhaseTypePar built
     * from the same PH. A 2-phase one is unaffected.
     */
    @Test
    public void nPhaseHyperExpExportsAsPhaseTypeInJsim() throws Exception {
        String threePhase = jsimServiceXml("hexp3",
                new HyperExp(new double[]{0.3, 0.4, 0.3}, new double[]{3.0, 1.0, 0.2}));
        assertTrue(threePhase.contains("jmt.engine.random.PhaseTypePar"),
                "3-phase HyperExp not exported as PhaseTypePar");
        assertFalse(threePhase.contains("jmt.engine.random.HyperExpPar"),
                "3-phase HyperExp exported as 2-phase HyperExpPar, truncating phases");

        String twoPhase = jsimServiceXml("hexp2", new HyperExp(0.25, 2.0, 5.0));
        assertTrue(twoPhase.contains("jmt.engine.random.HyperExpPar"),
                "2-phase HyperExp no longer exported as HyperExpPar");
        assertFalse(twoPhase.contains("jmt.engine.random.PhaseTypePar"),
                "2-phase HyperExp export changed to PhaseTypePar");
    }

    /** Writes a Source -> Queue -> Sink model with the given queue service to JSIM XML. */
    private String jsimServiceXml(String tag, Distribution service) throws Exception {
        Network model = new Network(tag);
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.1));
        queue.setService(oclass, service);
        model.link(Network.serialRouting(source, queue, sink));

        File out = tempDir.resolve(tag + ".jsimg").toFile();
        final String[] written = new String[1];
        withSuppressedOutput(() -> {
            try {
                SolverJMT solver = new SolverJMT(model);
                written[0] = solver.writeJSIM(model.getStruct(true), out.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError(tag + ": JSIM export failed: " + e, e);
            }
        });
        assertNotNull(written[0], tag + ": writeJSIM returned null path");
        File file = new File(written[0]);
        // Must remain well-formed XML
        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        dbf.newDocumentBuilder().parse(file);
        return new String(Files.readAllBytes(file.toPath()), StandardCharsets.UTF_8);
    }

    /** The manual's published LQN example uses Exp "rate" rather than "lambda". */
    @Test
    public void expRateAliasIsAccepted() throws Exception {
        File f = tempDir.resolve("exp_rate_alias.json").toFile();
        String json = "{\"format\":\"line-model\",\"version\":\"1.0\",\"model\":{"
                + "\"type\":\"Network\",\"name\":\"ratealias\","
                + "\"nodes\":[{\"name\":\"Source\",\"type\":\"Source\",\"service\":{"
                + "\"Class1\":{\"type\":\"Exp\",\"params\":{\"rate\":0.1}}}},"
                + "{\"name\":\"Queue\",\"type\":\"Queue\",\"scheduling\":\"FCFS\",\"service\":{"
                + "\"Class1\":{\"type\":\"Exp\",\"params\":{\"rate\":4.0}}}},"
                + "{\"name\":\"Sink\",\"type\":\"Sink\"}],"
                + "\"classes\":[{\"name\":\"Class1\",\"type\":\"Open\"}],"
                + "\"routing\":{\"type\":\"matrix\",\"matrix\":{\"Class1,Class1\":{"
                + "\"Source\":{\"Queue\":1.0},\"Queue\":{\"Sink\":1.0}}}}}}";
        Files.write(f.toPath(), json.getBytes(StandardCharsets.UTF_8));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = (Network) LineModelIO.load(f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("ratealias: load failed: " + e, e);
            }
        });
        Distribution svc = queueOf(out[0], "Queue").getService(out[0].getClasses().get(0));
        assertEquals(0.25, svc.getMean(), 1e-12, "Exp 'rate' alias not honoured");
    }

    /**
     * A spec-conformant cache described by the nested "cache" object used to
     * load with default parameters, because only the flat keys were read.
     */
    @Test
    public void nestedCacheObjectIsRead() throws Exception {
        File f = tempDir.resolve("nested_cache.json").toFile();
        String json = "{\"format\":\"line-model\",\"version\":\"1.0\",\"model\":{"
                + "\"type\":\"Network\",\"name\":\"nestedcache\","
                + "\"nodes\":[{\"name\":\"Source\",\"type\":\"Source\",\"service\":{"
                + "\"Class1\":{\"type\":\"Exp\",\"params\":{\"lambda\":0.1}}}},"
                + "{\"name\":\"Cache\",\"type\":\"Cache\",\"cache\":{"
                + "\"items\":9,\"capacity\":[3],\"replacement\":\"FIFO\",\"admissionProb\":0.75}},"
                + "{\"name\":\"Sink\",\"type\":\"Sink\"}],"
                + "\"classes\":[{\"name\":\"Class1\",\"type\":\"Open\"}],"
                + "\"routing\":{\"type\":\"matrix\",\"matrix\":{\"Class1,Class1\":{"
                + "\"Source\":{\"Cache\":1.0},\"Cache\":{\"Sink\":1.0}}}}}}";
        Files.write(f.toPath(), json.getBytes(StandardCharsets.UTF_8));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = (Network) LineModelIO.load(f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("nestedcache: load failed: " + e, e);
            }
        });
        jline.lang.nodes.Cache cache = null;
        for (Node n : out[0].getNodes()) {
            if (n instanceof jline.lang.nodes.Cache) {
                cache = (jline.lang.nodes.Cache) n;
            }
        }
        assertNotNull(cache, "nested cache object did not produce a Cache node");
        assertEquals(9, cache.getNumberOfItems(), "cache.items ignored");
        assertEquals(3.0, cache.getItemLevelCap().get(0), 1e-12, "cache.capacity ignored");
        assertEquals(jline.lang.constant.ReplacementStrategy.FIFO, cache.getReplacementStrategy(),
                "cache.replacement ignored");
        assertEquals(0.75, cache.getAdmissionProb(), 1e-12, "cache.admissionProb ignored");
    }

    /**
     * A CLIMB cache must survive the wire verbatim: policy CLIMB, one list of
     * capacity C. The writer used to export it remapped as FIFO on C unit-capacity
     * lists while emitting the original (nLevels+1)-square accessProb alongside,
     * so the reloaded cache carried C lists against a 2x2 access matrix and blew
     * up on the first hit. The CLIMB-to-FIFO rewrite belongs to refreshLocalVars,
     * which remaps itemcap and accost together.
     */
    @Test
    public void climbCacheSurvivesRoundTripVerbatim() throws Exception {
        Network model = new Network("climbcache");
        Delay delay = new Delay(model, "Delay");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(
                model, "Cache", 5, 2, jline.lang.constant.ReplacementStrategy.CLIMB);
        ClosedClass jobClass = new ClosedClass(model, "JobClass", 1, delay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, delay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, delay, 0);
        delay.setService(jobClass, new Exp(1.0));
        cacheNode.setRead(jobClass, new jline.lang.processes.Zipf(1.2, 5));
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, delay, cacheNode, 1.0);
        P.set(hitClass, jobClass, cacheNode, delay, 1.0);
        P.set(missClass, jobClass, cacheNode, delay, 1.0);
        model.link(P);

        // Populate the lazily-built accessProb before saving: this is the state in
        // which the remap used to emit a capacity and an accessProb of different
        // geometries.
        model.getStruct(false);

        Network reloaded = jsonRoundTrip(model, "climbcache");
        jline.lang.nodes.Cache rc = null;
        for (Node n : reloaded.getNodes()) {
            if (n instanceof jline.lang.nodes.Cache) {
                rc = (jline.lang.nodes.Cache) n;
            }
        }
        assertNotNull(rc, "CLIMB cache did not survive the round trip");
        assertEquals(jline.lang.constant.ReplacementStrategy.CLIMB, rc.getReplacementStrategy(),
                "CLIMB was downgraded to another policy on the wire");
        assertEquals(1, rc.getItemLevelCap().getNumCols(),
                "CLIMB capacity was remapped to unit-capacity lists on the wire");
        assertEquals(2.0, rc.getItemLevelCap().get(0), 1e-12, "CLIMB capacity was not preserved");
        assertEquals(5, rc.getNumberOfItems(), "CLIMB cache item count was not preserved");

        // The reloaded model must be solvable: refreshLocalVars rewrites CLIMB into
        // the FIFO unit-list form, remapping itemcap and accost consistently.
        reloaded.getStruct(false);
    }

    /**
     * An unrecognized type carrying mean+SCV must rebuild an APH matching BOTH
     * moments. It used to become an Exp (SCV discarded), or, absent params, an
     * Immediate: zero service in place of the real distribution.
     */
    @Test
    public void unknownTypeFallbackMatchesMeanAndScv() throws Exception {
        File f = tempDir.resolve("unknown_type.json").toFile();
        String json = "{\"format\":\"line-model\",\"version\":\"1.0\",\"model\":{"
                + "\"type\":\"Network\",\"name\":\"unknowntype\","
                + "\"nodes\":[{\"name\":\"Source\",\"type\":\"Source\",\"service\":{"
                + "\"Class1\":{\"type\":\"Exp\",\"params\":{\"lambda\":0.1}}}},"
                + "{\"name\":\"Queue\",\"type\":\"Queue\",\"scheduling\":\"FCFS\",\"service\":{"
                + "\"Class1\":{\"type\":\"NotARealFamily\",\"params\":{\"mean\":2.5,\"scv\":4.0}}}},"
                + "{\"name\":\"Sink\",\"type\":\"Sink\"}],"
                + "\"classes\":[{\"name\":\"Class1\",\"type\":\"Open\"}],"
                + "\"routing\":{\"type\":\"matrix\",\"matrix\":{\"Class1,Class1\":{"
                + "\"Source\":{\"Queue\":1.0},\"Queue\":{\"Sink\":1.0}}}}}}";
        Files.write(f.toPath(), json.getBytes(StandardCharsets.UTF_8));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = (Network) LineModelIO.load(f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("unknowntype: load failed: " + e, e);
            }
        });
        Distribution svc = queueOf(out[0], "Queue").getService(out[0].getClasses().get(0));
        assertTrue(svc instanceof APH,
                "unknown type rebuilt as " + svc.getClass().getSimpleName() + ", not an APH");
        assertEquals(2.5, svc.getMean(), 1e-8, "unknown-type fallback lost the mean");
        assertEquals(4.0, svc.getSCV(), 1e-8, "unknown-type fallback discarded the SCV");
    }

    /** fitCentral is documented and is a real APH API; it used to fit the mean alone. */
    @Test
    public void fitCentralIsHonoured() throws Exception {
        File f = tempDir.resolve("fit_central.json").toFile();
        String json = "{\"format\":\"line-model\",\"version\":\"1.0\",\"model\":{"
                + "\"type\":\"Network\",\"name\":\"fitcentral\","
                + "\"nodes\":[{\"name\":\"Source\",\"type\":\"Source\",\"service\":{"
                + "\"Class1\":{\"type\":\"Exp\",\"params\":{\"lambda\":0.1}}}},"
                + "{\"name\":\"Queue\",\"type\":\"Queue\",\"scheduling\":\"FCFS\",\"service\":{"
                + "\"Class1\":{\"type\":\"APH\",\"fit\":{\"method\":\"fitCentral\","
                + "\"moments\":[1.0,0.99,1.999]}}}},"
                + "{\"name\":\"Sink\",\"type\":\"Sink\"}],"
                + "\"classes\":[{\"name\":\"Class1\",\"type\":\"Open\"}],"
                + "\"routing\":{\"type\":\"matrix\",\"matrix\":{\"Class1,Class1\":{"
                + "\"Source\":{\"Queue\":1.0},\"Queue\":{\"Sink\":1.0}}}}}}";
        Files.write(f.toPath(), json.getBytes(StandardCharsets.UTF_8));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = (Network) LineModelIO.load(f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("fitcentral: load failed: " + e, e);
            }
        });
        Distribution svc = queueOf(out[0], "Queue").getService(out[0].getClasses().get(0));
        APH oracle = APH.fitCentral(1.0, 0.99, 1.999);
        assertEquals(oracle.getMean(), svc.getMean(), 1e-9,
                "fitCentral mean not honoured");
        // The variance is the moment fitMean would have discarded
        assertEquals(oracle.getSCV(), svc.getSCV(), 1e-9,
                "fitCentral degraded to a mean-only fit");
    }

    /**
     * Expolynomial had no class in the JAR at all, so a MATLAB- or Python-authored
     * model carrying one decoded to a zero-service Immediate. The density is an
     * expolynomial expression string and the type carries no moments, so nothing about it
     * can be moment-matched: it must survive verbatim as its own family.
     */
    @Test
    public void expolynomialSurvivesRoundTrip() throws Exception {
        Expolynomial oracle = new Expolynomial("2 * x^1 * Exp[-3 x]", 0.5, 4.0);
        Network model = openModel("expoly", oracle);
        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "expoly");
            } catch (Exception e) {
                throw new AssertionError("expoly: round trip failed: " + e, e);
            }
        });

        // The wire shape must stay identical to the Python writer's, which nests
        // the density under an "expolynomial" object rather than under "params".
        String wire = new String(
                Files.readAllBytes(tempDir.resolve("expoly_wire.json")), StandardCharsets.UTF_8);
        assertTrue(wire.contains("\"expolynomial\""),
                "Expolynomial not written under the nested \"expolynomial\" key");
        assertTrue(wire.contains("\"density\"") && wire.contains("\"eft\"")
                        && wire.contains("\"lft\""),
                "Expolynomial wire object lost one of density/eft/lft");

        Distribution arv = sourceOf(out[0]).getArrivalDistribution(out[0].getClasses().get(0));
        assertFalse(arv instanceof jline.lang.processes.Immediate,
                "Expolynomial decoded as a zero-service Immediate");
        assertTrue(arv instanceof Expolynomial,
                "Expolynomial decoded as " + arv.getClass().getSimpleName());
        Expolynomial back = (Expolynomial) arv;
        assertEquals(oracle.getDensity(), back.getDensity(),
                "Expolynomial density expression changed in round trip");
        assertEquals(oracle.getEft(), back.getEft(), 1e-12,
                "Expolynomial eft changed in round trip");
        assertEquals(oracle.getLft(), back.getLft(), 1e-12,
                "Expolynomial lft changed in round trip");
        // The moments are NaN by construction in all three codebases: a reader
        // that produced a finite mean here would have substituted another family.
        assertTrue(Double.isNaN(back.getMean()), "Expolynomial mean is no longer NaN");
        assertTrue(Double.isNaN(back.getSCV()), "Expolynomial SCV is no longer NaN");
        assertEquals(oracle.getEft(), back.getProcess().get(0).get(0, 0), 1e-12,
                "Expolynomial process representation lost eft");
        assertEquals(oracle.getLft(), back.getProcess().get(1).get(0, 0), 1e-12,
                "Expolynomial process representation lost lft");
    }

    /**
     * An unbounded latest firing time has no JSON number literal, so both the
     * MATLAB and the Python writer emit the string "Inf". The JAR must emit and
     * accept the same token, otherwise Inf silently becomes a finite bound.
     */
    @Test
    public void expolynomialUnboundedLftUsesInfToken() throws Exception {
        Expolynomial oracle = new Expolynomial("1 * x^0 * Exp[-1 x]", 0.0,
                Double.POSITIVE_INFINITY);
        Network model = openModel("expolyinf", oracle);
        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "expolyinf");
            } catch (Exception e) {
                throw new AssertionError("expolyinf: round trip failed: " + e, e);
            }
        });
        String wire = new String(
                Files.readAllBytes(tempDir.resolve("expolyinf_wire.json")), StandardCharsets.UTF_8);
        assertTrue(wire.contains("\"lft\": \"Inf\"") || wire.contains("\"lft\":\"Inf\""),
                "an unbounded lft was not written as the \"Inf\" token: " + wire);

        Distribution arv = sourceOf(out[0]).getArrivalDistribution(out[0].getClasses().get(0));
        assertTrue(arv instanceof Expolynomial,
                "Expolynomial decoded as " + arv.getClass().getSimpleName());
        assertTrue(Double.isInfinite(((Expolynomial) arv).getLft()),
                "an unbounded lft became finite on load");
    }

    /**
     * A file authored by MATLAB or Python is read with the same nested shape;
     * this pins the reader against the wire form rather than against our writer.
     */
    @Test
    public void expolynomialForeignFileIsRead() throws Exception {
        File f = tempDir.resolve("expoly_foreign.json").toFile();
        String json = "{\"format\":\"line-model\",\"version\":\"1.0\",\"model\":{"
                + "\"type\":\"Network\",\"name\":\"expolyforeign\","
                + "\"nodes\":[{\"name\":\"Source\",\"type\":\"Source\",\"service\":{"
                + "\"Class1\":{\"type\":\"Expolynomial\",\"expolynomial\":{"
                + "\"density\":\"0.5 * x^2 * Exp[-2 x]\",\"eft\":1.0,\"lft\":\"Inf\"}}}},"
                + "{\"name\":\"Queue\",\"type\":\"Queue\",\"scheduling\":\"FCFS\",\"service\":{"
                + "\"Class1\":{\"type\":\"Exp\",\"params\":{\"lambda\":10.0}}}},"
                + "{\"name\":\"Sink\",\"type\":\"Sink\"}],"
                + "\"classes\":[{\"name\":\"Class1\",\"type\":\"Open\"}],"
                + "\"routing\":{\"type\":\"matrix\",\"matrix\":{\"Class1,Class1\":{"
                + "\"Source\":{\"Queue\":1.0},\"Queue\":{\"Sink\":1.0}}}}}}";
        Files.write(f.toPath(), json.getBytes(StandardCharsets.UTF_8));

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = (Network) LineModelIO.load(f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("expolyforeign: load failed: " + e, e);
            }
        });
        Distribution arv = sourceOf(out[0]).getArrivalDistribution(out[0].getClasses().get(0));
        assertTrue(arv instanceof Expolynomial,
                "foreign Expolynomial decoded as " + arv.getClass().getSimpleName());
        Expolynomial back = (Expolynomial) arv;
        assertEquals("0.5 * x^2 * Exp[-2 x]", back.getDensity(),
                "foreign Expolynomial density expression not read");
        assertEquals(1.0, back.getEft(), 1e-12, "foreign Expolynomial eft not read");
        assertTrue(Double.isInfinite(back.getLft()),
                "foreign Expolynomial \"Inf\" lft not read as infinite");
    }

    /**
     * An explicit firing priority of 0 is a legal JMT value and must survive the
     * round trip. The writer omits the key only when it equals the builder
     * default of 1, so 0 travels on the wire and 1 does not. Guarding on
     * {@code > 0} instead dropped 0 and preserved the redundant 1, and every
     * reader then restored the default: see BUGS.md BUG-90.
     */
    @Test
    public void firingPriorityZeroSurvivesTheRoundTrip() throws Exception {
        double[] priorities = {0.0, 1.0, 3.0};
        for (int p = 0; p < priorities.length; p++) {
            final double prio = priorities[p];
            Network model = new Network("spnprio");
            Source source = new Source(model, "Source");
            Sink sink = new Sink(model, "Sink");
            Place place = new Place(model, "P1");
            Transition trans = new Transition(model, "T1");
            OpenClass jobclass = new OpenClass(model, "Class1", 0);
            source.setArrival(jobclass, Exp.fitMean(1.0));
            Mode mode = trans.addMode("Mode1");
            trans.setNumberOfServers(mode, Integer.MAX_VALUE);
            trans.setDistribution(mode, new Exp(4));
            trans.setEnablingConditions(mode, jobclass, place, 1);
            trans.setFiringOutcome(mode, jobclass, sink, 1);
            trans.setFiringPriorities(mode, (int) prio);
            model.link(Network.serialRouting(source, place, trans, sink));

            File f = tempDir.resolve("spnprio_" + p + ".json").toFile();
            final Network[] out = new Network[1];
            withSuppressedOutput(() -> {
                try {
                    LineModelIO.save(model, f.getAbsolutePath());
                    out[0] = (Network) LineModelIO.load(f.getAbsolutePath());
                } catch (Exception e) {
                    throw new AssertionError("firingPriority " + prio + ": round trip failed: " + e, e);
                }
            });

            String wire = new String(Files.readAllBytes(f.toPath()), StandardCharsets.UTF_8);
            assertEquals(prio != 1.0, wire.contains("\"firingPriority\""),
                    "firingPriority " + prio + ": the key is written iff it differs from the default 1");

            Transition back = null;
            for (Node n : out[0].getNodes()) {
                if (n instanceof Transition) {
                    back = (Transition) n;
                }
            }
            assertNotNull(back, "firingPriority " + prio + ": no Transition after the round trip");
            assertEquals(prio, back.firingPriorities.get(0), 1e-12,
                    "firingPriority " + prio + " did not survive save/load");
        }
    }

    /**
     * `fanIn` is a MAP keyed by the SOURCE task, exactly as `fanOut` is keyed
     * by the dest. This used to be written and read as
     * `{"source":..,"value":..}` by the JAR ALONE, so every MATLAB- or
     * Python-written layered document with a fan-in made loadLayeredNetwork
     * throw NullPointerException before a solver ran -- the lqn_sockshop [P2J]
     * parity row. The wire text is asserted and not just the round trip,
     * because a JAR-only shape round-trips through the JAR perfectly.
     */
    @Test
    public void fanInIsAMapKeyedBySourceTask() throws Exception {
        LayeredNetwork model = new LayeredNetwork("fanin_wire");
        Processor p1 = new Processor(model, "P1", 1, SchedStrategy.INF);
        Processor p2 = new Processor(model, "P2", 1, SchedStrategy.FCFS);
        Task caller = new Task(model, "T1", 1, SchedStrategy.REF).on(p1);
        Task callee = new Task(model, "T2", 1, SchedStrategy.FCFS).on(p2);
        callee.setFanIn("T1", 3);
        caller.setFanOut("T2", 2);
        Entry e1 = new Entry(model, "E1").on(caller);
        Entry e2 = new Entry(model, "E2").on(callee);
        Activity a1 = new Activity(model, "A1", new Exp(1.0)).on(caller).boundTo(e1).synchCall(e2, 1);
        Activity a2 = new Activity(model, "A2", new Exp(1.0)).on(callee).boundTo(e2).repliesTo(e2);

        File f = tempDir.resolve("fanin_wire.json").toFile();
        final Object[] out = new Object[1];
        withSuppressedOutput(() -> {
            try {
                LineModelIO.save(model, f.getAbsolutePath());
                out[0] = LineModelIO.load(f.getAbsolutePath());
            } catch (Exception e) {
                throw new AssertionError("fanIn round trip failed: " + e, e);
            }
        });

        String wire = new String(Files.readAllBytes(f.toPath()), StandardCharsets.UTF_8);
        assertTrue(wire.contains("\"fanIn\""), "fanIn was not written at all");
        assertFalse(wire.contains("\"source\""),
                "fanIn is keyed by the source task name, not by a \"source\" property");

        assertTrue(out[0] instanceof LayeredNetwork, "loaded object is not a LayeredNetwork");
        LayeredNetwork back = (LayeredNetwork) out[0];
        Task backCallee = null;
        Task backCaller = null;
        for (Task t : back.getTasks().values()) {
            if ("T2".equals(t.getName())) {
                backCallee = t;
            }
            if ("T1".equals(t.getName())) {
                backCaller = t;
            }
        }
        assertNotNull(backCallee, "T2 is missing after the round trip");
        assertNotNull(backCaller, "T1 is missing after the round trip");
        assertEquals("T1", backCallee.getFanInSource(), "the fan-in source did not survive");
        assertEquals(3, backCallee.getFanInValue(), "the fan-in value did not survive");
        assertEquals(Integer.valueOf(2), backCaller.getFanOutMap().get("T2"),
                "the fan-out did not survive");
    }
}
