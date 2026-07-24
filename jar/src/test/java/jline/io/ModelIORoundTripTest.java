/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.Region;
import jline.lang.RoutingMatrix;
import jline.lang.SelfLoopingClass;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.ImpatienceType;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.LayeredNetwork;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
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
import java.util.List;

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
        });
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
        delay.setLimitedClassDependence((Matrix n) -> rowVector(1.0 + n.get(0)));
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

    /**
     * d5: RL routing carried only the bare strategy name, so sub_rl found no
     * value function after a round trip and silently degraded to its JSQ
     * fallback.
     */
    @Test
    public void rlRoutingParametersSurviveRoundTrip() throws Exception {
        Network model = new Network("rlrouting");
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
        Matrix vf = rowVector(0.5, 1.5, 2.5, 3.5);
        int[] vfShape = new int[]{2, 2};
        int[] actionNodes = new int[]{source.getNodeIndex()};
        source.setRLRouting(oclass, vf, vfShape, actionNodes, 0);

        final Network[] out = new Network[1];
        withSuppressedOutput(() -> {
            try {
                out[0] = jsonRoundTrip(model, "rlrouting");
            } catch (Exception e) {
                throw new AssertionError("rlrouting: round trip failed: " + e, e);
            }
        });
        Node backSource = out[0].getNodes().get(0);
        JobClass jc = out[0].getClasses().get(0);
        assertEquals(RoutingStrategy.RL, backSource.getRoutingStrategy(jc),
                "RL strategy lost in round trip");
        jline.lang.NodeParam np = out[0].getStruct().nodeparam.get(backSource);
        assertNotNull(np, "RL node parameters lost in round trip");
        assertNotNull(np.rlValueFunction.get(jc), "RL value function lost in round trip");
        assertEquals(3.5, np.rlValueFunction.get(jc).get(3), 1e-12,
                "RL value function changed in round trip");
        assertArrayEquals(vfShape, np.rlValueFunctionShape.get(jc),
                "RL value-function shape lost in round trip");
        assertEquals(0, np.rlStateSize.get(jc).intValue(), "RL stateSize lost in round trip");
        assertArrayEquals(actionNodes, np.rlNodesNeedAction.get(jc),
                "RL action nodes lost in round trip");
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
     * model carrying one decoded to a zero-service Immediate. The density is a
     * Sirio expression string and the type carries no moments, so nothing about it
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
}
