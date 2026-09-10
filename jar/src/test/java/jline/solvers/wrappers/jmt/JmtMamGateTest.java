package jline.solvers.wrappers.jmt;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.processes.Zipf;
import jline.solvers.SolverOptions;
import jline.solvers.mam.SolverMAM;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The per-method support gates of SolverJMT and SolverMAM.
 *
 * <p>findSolver / model.help() reports one row per (solver, method) pair, and
 * the row comes from supportsModelMethod -- the same gate
 * SolverAUTO.chooseSolverRanked uses before delegating. A gate weaker than the
 * analyzer offers a pair that then RAISES, or worse returns a table of ZEROS
 * under a method labelled exact. Three such gaps are pinned here:
 *
 * <ul>
 *   <li>SolverJMT declared ONE envelope for TWO engines. The 'jmva.*' names
 *       drive the JMVA ANALYTICAL engine, whose document carries a station
 *       type, a per-chain demand, a per-chain visit count, the populations and
 *       a reference station -- no cache, no fork, no region, no discipline. On
 *       a three-class LRU cache model all eight closed-form jmva methods
 *       returned an entirely zero table with no error.</li>
 *   <li>'replication' integrates a transient mean over [0,T], which is an
 *       OPTION rather than a model feature and so invisible to a feature set.</li>
 *   <li>SolverMAM's 'dec.mmap' is an open-network departure-process fixed
 *       point and 'retrial' needs an impatience configuration to analyze -- a
 *       MUST BE PRESENT rule a feature set cannot state at all.</li>
 * </ul>
 *
 * <p>No JMT subprocess is started: the assertions are about the gate, and the
 * numeric behaviour of the engines is SolverJMTJmvaTest's subject.
 */
public class JmtMamGateTest {

    private static final List<String> JMVA_CLOSED_ONLY = Arrays.asList(
            "jmva.amva", "jmva.recal", "jmva.comom", "jmva.chow",
            "jmva.bs", "jmva.aql", "jmva.lin", "jmva.dmlin");
    private static final List<String> JMVA_EXACT = Arrays.asList("jmva", "jmva.mva");

    // ------------------------------------------------------------------
    // models
    // ------------------------------------------------------------------

    /** Source -&gt; FCFS Queue -&gt; Sink: open, single server, product form. */
    private static Network mm1() {
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass c = new OpenClass(model, "C1");
        source.setArrival(c, new Exp(1.0));
        queue.setService(c, new Exp(2.0));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /** Delay -&gt; FCFS Queue, N = 3: the shape every jmva algorithm serves. */
    private static Network closedSingleServer() {
        Network model = new Network("repairmen");
        Delay think = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C1", 3, think);
        think.setService(c, new Exp(1.0));
        queue.setService(c, new Exp(2.0));
        model.link(Network.serialRouting(think, queue));
        return model;
    }

    /** The same with c = 3, which the closed-form jmva algorithms refuse. */
    private static Network closedMultiserver() {
        Network model = new Network("multiserver");
        Delay think = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(3);
        ClosedClass c = new ClosedClass(model, "C1", 4, think);
        think.setService(c, new Exp(1.0));
        queue.setService(c, new Exp(2.0));
        model.link(Network.serialRouting(think, queue));
        return model;
    }

    /** Two open classes at a head-of-line priority queue. */
    private static Network openHol() {
        Network model = new Network("prio");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.HOL);
        Sink sink = new Sink(model, "Sink");
        OpenClass hi = new OpenClass(model, "Hi", 0);
        OpenClass lo = new OpenClass(model, "Lo", 1);
        source.setArrival(hi, new Exp(0.4));
        source.setArrival(lo, new Exp(0.4));
        queue.setService(hi, new Exp(2.0));
        queue.setService(lo, new Exp(2.0));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /**
     * Delay -&gt; FCFS Queue with cap 2 and N = 4: the buffer BINDS.
     *
     * <p>LINE blocks a closed job that finds no room; no JMT drop strategy
     * reproduces that, and the JMVA document has no capacity element at all.
     */
    private static Network closedBindingBuffer() {
        Network model = new Network("blk");
        Delay think = new Delay(model, "D");
        Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
        queue.setCapacity(2);
        ClosedClass c = new ClosedClass(model, "C", 4, think);
        think.setService(c, new Exp(1.0));
        queue.setService(c, new Exp(2.0));
        model.link(Network.serialRouting(think, queue));
        return model;
    }

    /**
     * M/M/1/2: an OPEN binding buffer, which JSIM DOES simulate.
     *
     * <p>JMT's queue section carries the drop rule directly, so a refused
     * arrival is lost in JMT exactly as it is in LINE. The over-tightening guard.
     */
    private static Network openBindingBuffer() {
        Network model = new Network("mm1k");
        Source source = new Source(model, "S");
        Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
        queue.setCapacity(2);
        Sink sink = new Sink(model, "K");
        OpenClass c = new OpenClass(model, "C");
        source.setArrival(c, new Exp(1.0));
        queue.setService(c, new Exp(2.0));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /**
     * Source -&gt; Delay -&gt; PS Queue -&gt; Sink with ONE OPEN and ONE CLOSED class:
     * jline.examples.java.basic.MixedModel.mqn_basic.
     *
     * <p>The shape the binding-buffer rule got wrong. Nobody caps anything here,
     * but refreshCapacity DERIVES a station capacity from the classes served,
     * and in THIS port an unbounded contribution is Integer.MAX_VALUE rather
     * than Inf -- so the Delay and Queue1 came out at 2147483647 + 2 =
     * 2147483649 and were read as binding buffers, refusing five
     * MixedExamplesTest models that have no buffer at all.
     */
    private static Network mixedOpenClosed() {
        Network model = new Network("mqn_basic");
        Delay think = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");
        ClosedClass cc = new ClosedClass(model, "ClosedClass", 2, think, 0);
        OpenClass oc = new OpenClass(model, "OpenClass", 0);
        think.setService(cc, new Exp(1.0));
        think.setService(oc, new Exp(3.0));
        queue.setService(cc, new Exp(2.0));
        queue.setService(oc, new Exp(1.0));
        source.setArrival(oc, new Exp(0.1));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(cc, cc, think, queue, 1.00);
        P.set(cc, cc, queue, think, 1.00);
        P.set(oc, oc, source, think, 1.00);
        P.set(oc, oc, think, queue, 1.00);
        P.set(oc, oc, queue, sink, 1.00);
        model.link(P);
        return model;
    }

    /**
     * Two closed classes at one PS queue: sn.cap is DERIVED as 2 x N here.
     *
     * <p>The multi-class shape is the one that makes "finite cap" and "binding
     * cap" different questions, so it is the model the buffer rule must leave
     * alone.
     */
    private static Network cqnTwoClass() {
        Network model = new Network("cqn2");
        Delay think = new Delay(model, "D");
        Queue queue = new Queue(model, "Q", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, think);
        ClosedClass c2 = new ClosedClass(model, "C2", 1, think);
        think.setService(c1, new Exp(1.0));
        think.setService(c2, new Exp(2.0));
        queue.setService(c1, new Exp(2.0));
        queue.setService(c2, new Exp(3.0));
        model.link(Network.serialRouting(think, queue));
        return model;
    }

    /**
     * Source -&gt; Fork -&gt; two FCFS queues -&gt; Join -&gt; Sink, one open class.
     *
     * <p>dec.mmap's sweep uses the PLAIN traffic step, which has no
     * synchronization; the topology router sends this shape to Solver_mam_fj
     * from 'default' and 'dec.source' and never to dec.mmap.
     */
    private static Network openForkJoin() {
        Network model = new Network("fj");
        Source source = new Source(model, "S");
        Fork fork = new Fork(model, "F");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Join join = new Join(model, "J", fork);
        Sink sink = new Sink(model, "K");
        OpenClass c = new OpenClass(model, "C");
        source.setArrival(c, new Exp(0.5));
        q1.setService(c, new Exp(2.0));
        q2.setService(c, new Exp(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c, c, source, fork, 1.0);
        P.set(c, c, fork, q1, 1.0);
        P.set(c, c, fork, q2, 1.0);
        P.set(c, c, q1, join, 1.0);
        P.set(c, c, q2, join, 1.0);
        P.set(c, c, join, sink, 1.0);
        model.link(P);
        return model;
    }

    /** Client -&gt; LRU Cache -&gt; CacheDelay: the all-zero case. */
    private static Network cacheModel() {
        Network model = new Network("cache");
        Delay client = new Delay(model, "Client");
        Cache cache = new Cache(model, "Cache", 4, 2, ReplacementStrategy.LRU);
        Delay hitmiss = new Delay(model, "CacheDelay");
        ClosedClass cc = new ClosedClass(model, "ClientClass", 1, client, 0);
        ClosedClass hc = new ClosedClass(model, "HitClass", 0, client, 0);
        ClosedClass mc = new ClosedClass(model, "MissClass", 0, client, 0);
        client.setService(cc, new Immediate());
        hitmiss.setService(hc, Exp.fitMean(0.2));
        hitmiss.setService(mc, Exp.fitMean(1.0));
        cache.setRead(cc, new Zipf(1.4, 4));
        cache.setHitClass(cc, hc);
        cache.setMissClass(cc, mc);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(cc, cc, client, cache, 1.0);
        P.set(hc, hc, cache, hitmiss, 1.0);
        P.set(mc, mc, cache, hitmiss, 1.0);
        P.set(hc, cc, hitmiss, client, 1.0);
        P.set(mc, cc, hitmiss, client, 1.0);
        model.link(P);
        return model;
    }

    private static SolverJMT jmt(Network model, String method) {
        SolverOptions options = new SolverOptions(SolverType.JMT);
        options.method = method;
        options.verbose = VerboseLevel.SILENT;
        return new SolverJMT(model, options);
    }

    private static SolverMAM mam(Network model, String method) {
        SolverOptions options = new SolverOptions(SolverType.MAM);
        options.method = method;
        options.verbose = VerboseLevel.SILENT;
        return new SolverMAM(model, options);
    }

    // ------------------------------------------------------------------
    // SolverJMT: the two engines have two envelopes
    // ------------------------------------------------------------------

    @Test
    @DisplayName("the JMVA envelope drops what write_jmva cannot serialize")
    public void jmvaEnvelopeIsNarrowerThanJsim() {
        FeatureSet jsim = SolverJMT.getFeatureSet();
        FeatureSet jmva = SolverJMT.getJMVAFeatureSet();
        // The constructs the JMVA document has no element for at all.
        String[] absent = new String[]{"Cache", "CacheClassSwitcher", "ReplacementStrategy_LRU",
                "Fork", "Join", "Region", "Place", "Transition", "Reneging", "Balking",
                "HeteroServers", "SetupDelayOff", "ServerParallelism"};
        for (int i = 0; i < absent.length; i++) {
            assertTrue(jsim.inspectFeature(absent[i]), absent[i] + " left the JSIM envelope");
            assertFalse(jmva.inspectFeature(absent[i]), absent[i] + " is still declared for JMVA");
        }
        // The writer emits NO discipline, so only the BCMP station types survive.
        String[] noDiscipline = new String[]{"SchedStrategy_HOL", "SchedStrategy_DPS",
                "SchedStrategy_GPS", "SchedStrategy_POLLING", "SchedStrategy_SRPT"};
        for (int i = 0; i < noDiscipline.length; i++) {
            assertFalse(jmva.inspectFeature(noDiscipline[i]),
                    noDiscipline[i] + " is still declared for JMVA");
        }
        String[] bcmp = new String[]{"SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_FCFS", "SchedStrategy_LCFSPR"};
        for (int i = 0; i < bcmp.length; i++) {
            assertTrue(jmva.inspectFeature(bcmp[i]), bcmp[i] + " was lost from JMVA");
        }
        // A mean demand is all JMVA reads, so the laws stay.
        String[] laws = new String[]{"Exp", "Erlang", "HyperExp", "Pareto", "Replayer"};
        for (int i = 0; i < laws.length; i++) {
            assertTrue(jmva.inspectFeature(laws[i]), laws[i] + " was lost from JMVA");
        }
    }

    @Test
    @DisplayName("the closed-only JMVA algorithms drop OpenClass and LoadDependence")
    public void closedOnlyAlgorithmsDropTwoMoreNames() {
        SolverJMT probe = jmt(closedSingleServer(), "jmva");
        for (String method : JMVA_EXACT) {
            assertFalse(SolverJMT.jmvaIsClosedOnly(method), method + " was called closed-only");
            assertTrue(probe.getMethodFeatureSet(method).inspectFeature("OpenClass"));
            assertTrue(probe.getMethodFeatureSet(method).inspectFeature("LoadDependence"));
        }
        for (String method : JMVA_CLOSED_ONLY) {
            assertTrue(SolverJMT.jmvaIsClosedOnly(method), method + " was not called closed-only");
            FeatureSet narrowed = probe.getMethodFeatureSet(method);
            assertFalse(narrowed.inspectFeature("OpenClass"), method + " still declares OpenClass");
            assertFalse(narrowed.inspectFeature("LoadDependence"),
                    method + " still declares LoadDependence");
            assertTrue(narrowed.inspectFeature("ClosedClass"));
        }
        // The simulation methods keep the whole JSIM envelope.
        for (String method : new String[]{"default", "jsim", "replication"}) {
            assertTrue(probe.getMethodFeatureSet(method).inspectFeature("Cache"),
                    method + " lost the JSIM cache");
        }
    }

    @Test
    @DisplayName("no jmva method is offered on a cache model, and the reason names Cache")
    public void jmvaRefusesACacheModel() {
        Network model = cacheModel();
        for (String method : JMVA_EXACT) {
            String reason = jmt(model, method).supportsModelMethod(method);
            assertNotEquals("", reason, method + " was offered on a cache model");
            assertTrue(reason.contains("Cache"), "the reason does not name Cache: " + reason);
        }
        for (String method : JMVA_CLOSED_ONLY) {
            String reason = jmt(model, method).supportsModelMethod(method);
            assertNotEquals("", reason, method + " was offered on a cache model");
            assertTrue(reason.contains("Cache"), "the reason does not name Cache: " + reason);
        }
        // The SIMULATOR does serialize a cache (SaveHandlers), so it stays.
        assertEquals("", jmt(model, "jsim").supportsModelMethod("jsim"));
    }

    @Test
    @DisplayName("the closed-only algorithms refuse an open model by feature name")
    public void closedOnlyAlgorithmsRefuseAnOpenModel() {
        Network model = mm1();
        for (String method : JMVA_CLOSED_ONLY) {
            String reason = jmt(model, method).supportsModelMethod(method);
            assertNotEquals("", reason, method + " was offered on an open model");
            assertTrue(reason.contains("OpenClass"), "the reason does not name OpenClass: " + reason);
        }
        for (String method : JMVA_EXACT) {
            assertEquals("", jmt(model, method).supportsModelMethod(method),
                    method + " was lost on an open product-form model");
        }
    }

    @Test
    @DisplayName("the closed-only algorithms refuse a multi-server station, gate and run alike")
    public void closedOnlyAlgorithmsRefuseAMultiserverModel() {
        Network model = closedMultiserver();
        for (String method : JMVA_CLOSED_ONLY) {
            String reason = jmt(model, method).supportsModelMethod(method);
            assertEquals(method + " does not support multi-server stations.", reason);
        }
        for (String method : JMVA_EXACT) {
            assertEquals("", jmt(model, method).supportsModelMethod(method),
                    method + " was lost on a multi-server model");
        }
        // ONE PREDICATE, TWO CALLERS: the writer raises what the gate says.
        RuntimeException raised = assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                jmt(closedMultiserver(), "jmva.amva").getAvgTable();
            }
        });
        assertTrue(String.valueOf(raised.getMessage()).contains("multi-server"),
                "the run did not repeat the gate's reason: " + raised.getMessage());
    }

    @Test
    @DisplayName("a priority model is not offered the analytical engine")
    public void priorityModelIsNotOfferedJmva() {
        Network model = openHol();
        for (String method : JMVA_EXACT) {
            String reason = jmt(model, method).supportsModelMethod(method);
            assertNotEquals("", reason, method + " was offered on a HOL model");
            assertTrue(reason.contains("SchedStrategy_HOL"), reason);
        }
        assertEquals("", jmt(model, "jsim").supportsModelMethod("jsim"));
    }

    @Test
    @DisplayName("replication needs a finite timespan, which no feature name states")
    public void replicationNeedsAFiniteTimespan() {
        Network model = mm1();
        SolverOptions unbounded = new SolverOptions(SolverType.JMT);
        unbounded.method = "replication";
        unbounded.verbose = VerboseLevel.SILENT;
        String reason = new SolverJMT(model, unbounded).supportsModelMethod("replication");
        assertNotEquals("", reason);
        assertTrue(reason.contains("finite timespan"), reason);

        SolverOptions bounded = new SolverOptions(SolverType.JMT);
        bounded.method = "replication";
        bounded.verbose = VerboseLevel.SILENT;
        bounded.timespan = new double[]{0.0, 10.0};
        assertEquals("", new SolverJMT(model, bounded).supportsModelMethod("replication"));
    }

    @Test
    @DisplayName("a closed single-server product-form network keeps every jmva method")
    public void jmvaIsNotOverTightened() {
        Network model = closedSingleServer();
        for (String method : JMVA_EXACT) {
            assertEquals("", jmt(model, method).supportsModelMethod(method), method);
        }
        for (String method : JMVA_CLOSED_ONLY) {
            assertEquals("", jmt(model, method).supportsModelMethod(method), method);
        }
        assertEquals("", jmt(model, "jsim").supportsModelMethod("jsim"));
    }

    @Test
    @DisplayName("a closed binding buffer is refused by BOTH engines, for opposite reasons")
    public void closedBindingBufferIsRefusedByBothEngines() {
        Network model = closedBindingBuffer();
        for (String method : new String[]{"default", "jsim"}) {
            String reason = jmt(model, method).supportsModelMethod(method);
            assertNotEquals("", reason, "jmt." + method + " was offered on a binding buffer");
            assertTrue(reason.contains("binds for the closed class"), reason);
        }
        for (String method : JMVA_EXACT) {
            String reason = jmt(model, method).supportsModelMethod(method);
            assertNotEquals("", reason, "jmt." + method + " was offered on a binding buffer");
            assertTrue(reason.contains("no capacity element"), reason);
        }
        for (String method : JMVA_CLOSED_ONLY) {
            assertNotEquals("", jmt(model, method).supportsModelMethod(method), method);
        }
    }

    @Test
    @DisplayName("the gate repeats the JSIM writer's own buffer sentence")
    public void theGateAndTheWriterGiveTheSameBufferSentence() {
        // ONE PREDICATE, TWO CALLERS: saveBufferCapacity raises what the gate
        // returns, because both ask SaveHandlers.jmtStationCapRefusal.
        Network model = closedBindingBuffer();
        String reason = jmt(model, "jsim").supportsModelMethod("jsim");
        assertEquals(jline.solvers.wrappers.jmt.handlers.SaveHandlers.jmtStationCapRefusal(
                model.getStruct(), 1), reason);
    }

    @Test
    @DisplayName("an open loss buffer keeps the simulator and still loses the analytical engine")
    public void openLossBufferKeepsTheSimulator() {
        // The over-tightening guard: a refused OPEN arrival is LOST, which JMT's
        // queue section expresses directly, so jsim genuinely simulates an
        // M/M/1/K and the row must survive.
        Network model = openBindingBuffer();
        assertEquals("", jmt(model, "jsim").supportsModelMethod("jsim"));
        assertEquals("", jmt(model, "default").supportsModelMethod("default"));
        for (String method : JMVA_EXACT) {
            String reason = jmt(model, method).supportsModelMethod(method);
            assertNotEquals("", reason, method);
            assertTrue(reason.contains("no capacity element"), reason);
        }
    }

    @Test
    @DisplayName("a derived capacity carrying an unbounded sentinel is not a buffer")
    public void mixedOpenClosedHasNoBindingBuffer() {
        // THE SENTINEL IS THE POINT. Station.cap is an int, so it cannot hold
        // Inf: its "no bound" value is Integer.MAX_VALUE, and refreshCapacity
        // SUMS that across the classes served, giving 2147483649 at a station
        // with one open and one closed class. Reading it as finite refused five
        // MixedExamplesTest models that declare no capacity at all.
        Network model = mixedOpenClosed();
        NetworkStruct sn = model.getStruct();
        boolean sawSentinel = false;
        for (int ist = 0; ist < sn.nstations; ist++) {
            if (sn.cap.get(ist) >= Integer.MAX_VALUE) {
                sawSentinel = true;
                assertTrue(jline.solvers.wrappers.jmt.handlers.SaveHandlers.jmtCapIsUnbounded(sn, ist),
                        "station " + ist + " cap " + sn.cap.get(ist) + " read as bounded");
            }
        }
        assertTrue(sawSentinel, "the mixed model no longer carries the derived sentinel");

        for (String method : new String[]{"default", "jsim"}) {
            assertEquals("", jmt(model, method).supportsModelMethod(method),
                    "jmt." + method + " was refused on a model with no buffer");
        }
        for (String method : JMVA_EXACT) {
            assertEquals("", jmt(model, method).supportsModelMethod(method),
                    "jmt." + method + " was refused on a model with no buffer");
        }
        // the closed-only algorithms go for the OPEN class, not for a buffer
        for (String method : JMVA_CLOSED_ONLY) {
            String reason = jmt(model, method).supportsModelMethod(method);
            assertNotEquals("", reason, method);
            assertTrue(reason.contains("OpenClass"), reason);
        }
    }

    @Test
    @DisplayName("a Source is never reported as carrying a buffer")
    public void aSourceNeverCarriesABuffer() {
        // A Source has no buffer that can bind: it IS the external world.
        // refreshCapacity still writes it a capacity row, and here that row is a
        // sum of unbounded sentinels, which is how a refusal came to name one.
        Network[] models = new Network[]{mixedOpenClosed(), openBindingBuffer(), mm1()};
        for (int i = 0; i < models.length; i++) {
            String[] methods = new String[]{"jsim", "jmva", "jmva.mva"};
            for (int j = 0; j < methods.length; j++) {
                String reason = jmt(models[i], methods[j]).supportsModelMethod(methods[j]);
                assertFalse(reason.contains("Source"), reason);
            }
        }
    }

    @Test
    @DisplayName("an uncapped model is untouched by the buffer rule")
    public void uncappedModelIsUntouchedByTheBufferRule() {
        // refreshCapacity DERIVES a finite sn.cap for a station nobody capped,
        // so the rule tests the cap against the population that can REACH the
        // station; reading "finite" as "binding" would refuse every closed model.
        Network[] models = new Network[]{closedSingleServer(), cqnTwoClass(), mixedOpenClosed()};
        for (int i = 0; i < models.length; i++) {
            assertEquals("", jmt(models[i], "jsim").supportsModelMethod("jsim"));
            for (String method : JMVA_EXACT) {
                assertEquals("", jmt(models[i], method).supportsModelMethod(method), method);
            }
        }
    }

    // ------------------------------------------------------------------
    // SolverMAM: dec.mmap and retrial
    // ------------------------------------------------------------------

    @Test
    @DisplayName("dec.mmap is refused on a closed model and on an INF station")
    public void decMmapIsAnOpenNetworkMethod() {
        String reason = mam(closedSingleServer(), "dec.mmap").supportsModelMethod("dec.mmap");
        assertNotEquals("", reason, "dec.mmap was offered on a closed model");
        assertTrue(reason.contains("ClosedClass"), reason);
        // A Delay is an INF station, which the analyzer's opening loop has no
        // branch for, so the reason has to name that too.
        assertTrue(reason.contains("SchedStrategy_INF"), reason);
    }

    @Test
    @DisplayName("dec.mmap is refused on a fork-join topology")
    public void decMmapIsRefusedOnForkJoin() {
        // The sweep has no synchronization, so a Fork is not a smaller model:
        // the router sends this shape to Solver_mam_fj from 'default' and
        // 'dec.source', and dec.mmap has no such route.
        Network model = openForkJoin();
        String reason = mam(model, "dec.mmap").supportsModelMethod("dec.mmap");
        assertNotEquals("", reason, "dec.mmap was offered on a fork-join model");
        assertTrue(reason.contains("Fork"), reason);
        // ... and the routes that ARE written for it stay
        assertEquals("", mam(model, "default").supportsModelMethod("default"));
        assertEquals("", mam(model, "dec.source").supportsModelMethod("dec.source"));
    }

    @Test
    @DisplayName("dec.mmap still serves the open models it is written for")
    public void decMmapIsNotOverTightened() {
        assertEquals("", mam(mm1(), "dec.mmap").supportsModelMethod("dec.mmap"));
        assertEquals("", mam(openHol(), "dec.mmap").supportsModelMethod("dec.mmap"));
    }

    @Test
    @DisplayName("dec.mmap raises rather than answering a table of zeros")
    public void decMmapRaisesOnAClosedModel() {
        // The analyzer used to warn and return the ZERO matrices it had
        // initialised, which SolverMAM reported as if it were the answer.
        RuntimeException raised = assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                mam(closedSingleServer(), "dec.mmap").getAvgTable();
            }
        });
        assertNotEquals("", String.valueOf(raised.getMessage()));
    }

    @Test
    @DisplayName("retrial is refused without an impatience configuration, gate and run alike")
    public void retrialNeedsARetrialTopology() {
        Network model = mm1();
        String reason = mam(model, "retrial").supportsModelMethod("retrial");
        assertNotEquals("", reason, "retrial was offered on a model with no orbit");
        assertTrue(reason.contains("retrial configuration"), reason);
        // ONE PREDICATE, TWO CALLERS: the analyzer repeats the gate's sentence.
        RuntimeException raised = assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                mam(mm1(), "retrial").getAvgTable();
            }
        });
        assertTrue(String.valueOf(raised.getMessage()).contains("retrial configuration"),
                "the run did not repeat the gate's reason: " + raised.getMessage());
        // and a closed model is refused for the reason qsys_is_retrial reports
        String closedReason = mam(closedSingleServer(), "retrial").supportsModelMethod("retrial");
        assertTrue(closedReason.contains("open queueing model"), closedReason);
    }

    @Test
    @DisplayName("the MAM methods a model can genuinely run are still offered")
    public void mamIsNotOverTightened() {
        Network model = closedSingleServer();
        for (String method : new String[]{"default", "dec.source", "dec.poisson", "ldqbd"}) {
            assertEquals("", mam(model, method).supportsModelMethod(method), method);
        }
    }
}
