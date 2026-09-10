package jline.solvers.ctmc;

import jline.lang.ClosedClass;
import jline.lang.FeatureSet;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.ssa.SolverSSA;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The per-method support gates of SolverCTMC, SolverFluid and SolverSSA.
 *
 * <p>WHAT THIS PINS. {@code Network.help()} / {@code findSolver()} reports one
 * row per (solver, method) pair and asks {@code supportsModelMethod(method)} --
 * the same gate {@code SolverAUTO} consults before delegating. Until the rules
 * below reached that gate it was much weaker than what the ANALYZERS enforce at
 * run time, so the report offered pairs that then threw:</p>
 *
 * <pre>
 *   ctmc cftp / cftp.approx   single-class and closed-only
 *   ctmc mdd                  single-class and closed-only, and no fork-join
 *   ctmc default/exact/gpu    the state space has to fit memory
 *   ctmc.* and ssa.*          the fork-join model class fjValidate admits
 *   fluid diffusion           closed-only, and no fork-join
 *   fluid kp                  no fork-join
 *   fluid refined             closed-only
 *   fluid tbi                 closed-only, and no cache
 *   fluid dae                 no OPEN fork-join model
 *   fluid mol / mtginf        a finite options.timespan
 * </pre>
 *
 * <p>Each test asserts the refusal AND its converse: a model the method IS
 * derived for must keep it. Over-tightening a gate hides a method the user
 * could have run, which is the same defect with the sign flipped -- and two of
 * the cases here ARE that flipped defect, found and fixed: native python's
 * SolverSSA declared no Fork/Join at all, and the C++ fluid envelope withheld
 * them from five methods that do integrate a fork-join model.</p>
 *
 * <p>The MATLAB twin is the same rules in {@code @SolverCTMC/supportsModelMethod.m},
 * {@code @SolverFLD/supportsModelMethod.m} and
 * {@code @SolverSSA/supportsModelMethod.m}; the Python twin is
 * python/tests/test_gate_ctmc_fluid.py and the C++ twin
 * cpp/tests/test_gate_ctmc_fluid.cpp.</p>
 */
public class CtmcFluidGateTest {

    /** Source -&gt; Queue -&gt; Sink, one open class: the smallest open model. */
    private static Network mm1() {
        Network model = new Network("mm1");
        Source source = new Source(model, "S");
        Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "K");
        OpenClass cl = new OpenClass(model, "C");
        source.setArrival(cl, new Exp(1.0));
        queue.setService(cl, new Exp(2.0));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    /** Delay -&gt; Queue -&gt; Delay, ONE closed class: what cftp and mdd are for. */
    private static Network repairmen() {
        Network model = new Network("rep");
        Delay delay = new Delay(model, "D");
        Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
        ClosedClass cl = new ClosedClass(model, "C", 3, delay, 0);
        delay.setService(cl, new Exp(1.0));
        queue.setService(cl, new Exp(2.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** The same shape with TWO closed classes: closed, but not single-class. */
    private static Network cqn2() {
        Network model = new Network("cqn2");
        Delay delay = new Delay(model, "D");
        Queue queue = new Queue(model, "Q", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, delay, 0);
        ClosedClass c2 = new ClosedClass(model, "C2", 1, delay, 0);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(2.0));
        queue.setService(c1, new Exp(2.0));
        queue.setService(c2, new Exp(3.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** Two open HOL classes: the state space the memory gate has to refuse. */
    private static Network prio() {
        Network model = new Network("prio");
        Source source = new Source(model, "S");
        Queue queue = new Queue(model, "Q", SchedStrategy.HOL);
        Sink sink = new Sink(model, "K");
        OpenClass hi = new OpenClass(model, "Hi", 0);
        OpenClass lo = new OpenClass(model, "Lo", 1);
        source.setArrival(hi, new Exp(0.4));
        source.setArrival(lo, new Exp(0.4));
        queue.setService(hi, new Exp(2.0));
        queue.setService(lo, new Exp(2.0));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    /**
     * Fork -&gt; two FCFS queues -&gt; Join.
     *
     * <p>{@code paired} names the Fork on the Join constructor, which is what
     * DECLARES the fork-join pairing: sn.fj is read off that declaration and not
     * derived from the routing, because a nested model (fj_basic_nesting) has
     * two forks and two joins the routing alone does not pair. A Join built
     * without it leaves sn.fj empty, which the exact construction refuses.</p>
     */
    private static Network forkjoin(boolean paired, boolean closed) {
        Network model = new Network("fj");
        Fork fork = new Fork(model, "F");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Join join = paired ? new Join(model, "J", fork) : new Join(model, "J");
        Node entry;
        Node exit;
        JobClass cl;
        if (closed) {
            Delay delay = new Delay(model, "D");
            ClosedClass cc = new ClosedClass(model, "C", 2, delay, 0);
            delay.setService(cc, new Exp(1.0));
            entry = delay;
            exit = delay;
            cl = cc;
        } else {
            Source source = new Source(model, "S");
            Sink sink = new Sink(model, "K");
            OpenClass oc = new OpenClass(model, "C");
            source.setArrival(oc, new Exp(0.5));
            entry = source;
            exit = sink;
            cl = oc;
        }
        q1.setService(cl, new Exp(2.0));
        q2.setService(cl, new Exp(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(cl, cl, entry, fork, 1.0);
        P.set(cl, cl, fork, q1, 1.0);
        P.set(cl, cl, fork, q2, 1.0);
        P.set(cl, cl, q1, join, 1.0);
        P.set(cl, cl, q2, join, 1.0);
        P.set(cl, cl, join, exit, 1.0);
        model.link(P);
        return model;
    }

    private static SolverCTMC ctmc(Network model, String method) {
        SolverOptions options = new SolverOptions(SolverType.CTMC);
        options.method = method;
        options.samples = 200;
        return new SolverCTMC(model, options);
    }

    private static SolverSSA ssa(Network model) {
        SolverOptions options = new SolverOptions(SolverType.SSA);
        options.method = "default";
        options.samples = 2000;
        options.seed = 23000;
        return new SolverSSA(model, options);
    }

    private static SolverFluid fluid(Network model, String method) {
        SolverOptions options = new SolverOptions(SolverType.FLUID);
        options.method = method;
        return new SolverFluid(model, options);
    }

    // -----------------------------------------------------------------------
    // SolverCTMC: cftp
    // -----------------------------------------------------------------------

    @Test
    public void cftpRefusesAnOpenModelInTheSamplersOwnWords() {
        String reason = ctmc(mm1(), "cftp").supportsModelMethod("cftp");
        assertFalse(reason.isEmpty());
        assertTrue(reason.contains("closed models only"), reason);
    }

    @Test
    public void cftpRefusesAMulticlassModelAndNamesTheClassCount() {
        String[] methods = new String[]{"cftp", "cftp.approx"};
        for (String method : methods) {
            String reason = ctmc(cqn2(), method).supportsModelMethod(method);
            assertFalse(reason.isEmpty(), method);
            assertTrue(reason.contains("single-class models only"), reason);
            assertTrue(reason.contains("2 classes"), reason);
        }
    }

    @Test
    public void cftpKeepsASingleClassClosedModel() {
        // The converse: refusing this one would hide a method that runs.
        String[] methods = new String[]{"cftp", "cftp.approx"};
        for (String method : methods) {
            assertEquals("", ctmc(repairmen(), method).supportsModelMethod(method), method);
        }
    }

    @Test
    public void theCftpGateAndTheAnalyzerAreOnePredicate() {
        Network model = cqn2();
        String reason = ctmc(model, "cftp").supportsModelMethod("cftp");
        RuntimeException err = assertThrows(RuntimeException.class,
                () -> ctmc(model, "cftp").getAvgQLen());
        assertTrue(err.getMessage().contains(reason), err.getMessage());
    }

    @Test
    public void theCftpEnvelopeDropsWhatTheRegistryCanName() {
        FeatureSet feats = ctmc(repairmen(), "cftp").getMethodFeatureSet("cftp");
        String[] dropped = new String[]{"OpenClass", "Source", "Sink", "Cache",
                "SchedStrategy_HOL", "Region", "LoadDependence", "RoutingStrategy_JSQ"};
        for (String name : dropped) {
            assertFalse(feats.inspectFeature(name), name);
        }
        // ... and keeps the product-form core it does serve.
        String[] kept = new String[]{"ClosedClass", "Queue", "Delay", "Exp",
                "SchedStrategy_FCFS", "SchedStrategy_PS", "SchedStrategy_INF",
                "SchedStrategy_SIRO", "SchedStrategy_LCFSPR"};
        for (String name : kept) {
            assertTrue(feats.inspectFeature(name), name);
            assertTrue(SolverCTMC.getFeatureSet().inspectFeature(name), name);
        }
    }

    // -----------------------------------------------------------------------
    // SolverCTMC: mdd
    // -----------------------------------------------------------------------

    @Test
    public void mddRefusesAnOpenModel() {
        String reason = ctmc(mm1(), "mdd").supportsModelMethod("mdd");
        assertFalse(reason.isEmpty());
        assertTrue(reason.contains("CLOSED networks") || reason.contains("OpenClass"), reason);
    }

    @Test
    public void mddRefusesAMulticlassModelAndNamesTheClassCount() {
        String reason = ctmc(cqn2(), "mdd").supportsModelMethod("mdd");
        assertFalse(reason.isEmpty());
        assertTrue(reason.contains("single-class networks"), reason);
        assertTrue(reason.contains("2 classes"), reason);
    }

    @Test
    public void mddKeepsASingleClassClosedModel() {
        assertEquals("", ctmc(repairmen(), "mdd").supportsModelMethod("mdd"));
    }

    @Test
    public void theMddGateAndTheAnalyzerAreOnePredicate() {
        Network model = cqn2();
        String reason = ctmc(model, "mdd").supportsModelMethod("mdd");
        RuntimeException err = assertThrows(RuntimeException.class,
                () -> ctmc(model, "mdd").getAvgQLen());
        assertTrue(err.getMessage().contains(reason), err.getMessage());
    }

    // -----------------------------------------------------------------------
    // SolverCTMC: the state-space methods
    // -----------------------------------------------------------------------

    @Test
    public void anIntractableChainIsRefusedByTheStateSpaceMethods() {
        Network model = prio();
        SolverOptions probe = new SolverOptions(SolverType.CTMC);
        assertFalse(SolverCTMC.isStateSpaceTractable(model, probe).ok,
                "the fixture must be intractable for this to test anything");
        String[] methods = new String[]{"default", "exact", "gpu"};
        for (String method : methods) {
            String reason = ctmc(model, method).supportsModelMethod(method);
            assertFalse(reason.isEmpty(), method);
            assertTrue(reason.toLowerCase().contains("memory"), reason);
        }
    }

    @Test
    public void aSmallChainKeepsTheStateSpaceMethods() {
        String[] methods = new String[]{"default", "exact", "gpu"};
        for (String method : methods) {
            assertEquals("", ctmc(repairmen(), method).supportsModelMethod(method), method);
        }
    }

    @Test
    public void theSamplerIsNotGatedOnAStateSpaceItNeverBuilds() {
        // cftp draws from the balance function, so the memory estimate that
        // stops "default" says nothing about it. Its own refusal here is the
        // class count, not the state space.
        String reason = ctmc(prio(), "cftp").supportsModelMethod("cftp");
        assertFalse(reason.isEmpty());
        assertFalse(reason.toLowerCase().contains("memory"), reason);
    }

    // -----------------------------------------------------------------------
    // SolverCTMC: the fork-join model class, which EVERY method has to clear
    // -----------------------------------------------------------------------

    @Test
    public void anUndeclaredPairingIsRefusedForEveryMethod() {
        Network model = forkjoin(false, true);
        String[] methods = new String[]{"default", "exact", "gpu", "mdd", "cftp", "cftp.approx"};
        for (String method : methods) {
            assertFalse(ctmc(model, method).supportsModelMethod(method).isEmpty(), method);
        }
    }

    @Test
    public void thePairingRefusalIsTheValidatorsOwnSentence() {
        String reason = ctmc(forkjoin(false, true), "default").supportsModelMethod("default");
        assertFalse(reason.isEmpty());
        assertTrue(reason.contains("without a matched Join")
                || reason.contains("no matched Join"), reason);
    }

    @Test
    public void anOpenClassThroughAForkIsRefused() {
        // The pairing is declared here, so this is the SECOND rule of the same
        // validator: the exact construction is stated for closed chains.
        String reason = ctmc(forkjoin(true, false), "default").supportsModelMethod("default");
        assertFalse(reason.isEmpty());
        assertTrue(reason.contains("Open classes routed through a Fork"), reason);
    }

    @Test
    public void aDeclaredClosedForkJoinKeepsTheStateSpaceMethods() {
        // The converse: the exact construction IS derived for this model, and
        // refusing it would hide the only exact answer a fork-join model has.
        Network model = forkjoin(true, true);
        String[] methods = new String[]{"default", "exact", "gpu"};
        for (String method : methods) {
            assertEquals("", ctmc(model, method).supportsModelMethod(method), method);
        }
    }

    @Test
    public void mddRefusesAForkJoinModelItCanNeverBeSingleClassFor() {
        // The tag augmentation adds one auxiliary class per branch, so the
        // struct that reaches the analyzer is never single-class however the
        // model was written. Stated as a feature so the reason names it.
        Network model = forkjoin(true, true);
        FeatureSet feats = ctmc(model, "mdd").getMethodFeatureSet("mdd");
        assertFalse(feats.inspectFeature("Fork"));
        assertFalse(feats.inspectFeature("Join"));
        assertTrue(SolverCTMC.getFeatureSet().inspectFeature("Fork"));
        String reason = ctmc(model, "mdd").supportsModelMethod("mdd");
        assertFalse(reason.isEmpty());
        assertTrue(reason.contains("Fork"), reason);
    }

    @Test
    public void theForkJoinGateAndTheAnalyzerAreOnePredicate() {
        Network model = forkjoin(false, true);
        String reason = ctmc(model, "default").supportsModelMethod("default");
        RuntimeException err = assertThrows(RuntimeException.class,
                () -> ctmc(model, "default").getAvgQLen());
        assertTrue(err.getMessage().contains(reason), err.getMessage());
    }

    // -----------------------------------------------------------------------
    // SolverSSA: the same construction, so the same gate
    // -----------------------------------------------------------------------

    @Test
    public void ssaDeclaresTheForkJoinNamesAndGatesTheWiring() {
        // The JAR has always DECLARED the four (native python alone omitted
        // them). Declaring without the structural gate would move the defect
        // rather than fix it: the wiring rules the names cannot state must
        // still refuse.
        String[] names = new String[]{"Fork", "Join", "Forker", "Joiner"};
        for (String name : names) {
            assertTrue(SolverSSA.getFeatureSet().inspectFeature(name), name);
        }
        assertEquals("", ssa(forkjoin(true, true)).supportsModelMethod("default"));

        String unpaired = ssa(forkjoin(false, true)).supportsModelMethod("default");
        assertFalse(unpaired.isEmpty());
        assertTrue(unpaired.contains("without a matched Join")
                || unpaired.contains("no matched Join"), unpaired);

        String open = ssa(forkjoin(true, false)).supportsModelMethod("default");
        assertFalse(open.isEmpty());
        assertTrue(open.contains("Open classes routed through a Fork"), open);

        // A model with no fork at all is never asked about.
        assertEquals("", ssa(repairmen()).supportsModelMethod("default"));
        assertEquals("", ssa(mm1()).supportsModelMethod("default"));
    }

    // -----------------------------------------------------------------------
    // SolverFluid
    // -----------------------------------------------------------------------

    @Test
    public void refinedIsClosedOnlyWhereverTheModelIsOpen() {
        // The reference's runAnalyzer has always refused this by name; the
        // featset never said so, and on an open fork-join model the restriction
        // surfaced as a failure inside the transform instead of a refusal.
        Network[] models = new Network[]{mm1(), forkjoin(true, false)};
        for (Network model : models) {
            String reason = fluid(model, "refined").supportsModelMethod("refined");
            assertFalse(reason.isEmpty());
            assertTrue(reason.contains("OpenClass"), reason);
        }
        assertFalse(fluid(mm1(), "refined").getMethodFeatureSet("refined").inspectFeature("OpenClass"));
        assertTrue(SolverFluid.getFeatureSet().inspectFeature("OpenClass"));
    }

    @Test
    public void refinedKeepsEveryClosedModelForkJoinIncluded() {
        Network[] models = new Network[]{repairmen(), cqn2(), forkjoin(true, true)};
        for (Network model : models) {
            assertEquals("", fluid(model, "refined").supportsModelMethod("refined"));
        }
    }

    @Test
    public void daeRefusesAnOpenForkJoinModelOnly() {
        // The conjunction the feature set cannot state: Fork and OpenClass are
        // both declared names, and it is having BOTH that the DAE form cannot
        // take -- the transform's auxiliary open classes have no unknown in it.
        String reason = fluid(forkjoin(true, false), "dae").supportsModelMethod("dae");
        assertFalse(reason.isEmpty());
        assertTrue(reason.contains("fork-join fixed point"), reason);
        assertTrue(reason.contains("minnormal"), reason);
        // ... and each half on its own stays runnable.
        Network[] models = new Network[]{mm1(), forkjoin(true, true)};
        for (Network model : models) {
            assertEquals("", fluid(model, "dae").supportsModelMethod("dae"));
        }
    }

    @Test
    public void diffusionAndKpDoNotIntegrateAForkJoinModelAtAll() {
        // Measured, not assumed. On the SYMMETRIC closed fork-join whose exact
        // chain is Q1 = Q2 = 0.664, "diffusion" put the whole population on ONE
        // station and zero elsewhere (a different station on a rerun), and "kp"
        // returned an all-zero table on a symmetric OPEN fork-join fed at rate
        // 0.5. Both answered instead of refusing, which is the reason the names
        // come off the envelope: a refusal is recoverable, a silent wrong answer
        // is not. C++ has always withheld them.
        String[] methods = new String[]{"diffusion", "kp"};
        String[] names = new String[]{"Fork", "Join", "Forker", "Joiner"};
        for (String method : methods) {
            FeatureSet feats = fluid(mm1(), method).getMethodFeatureSet(method);
            for (String name : names) {
                assertFalse(feats.inspectFeature(name), method + "/" + name);
            }
        }
        assertTrue(SolverFluid.getFeatureSet().inspectFeature("Fork"));
        String reason = fluid(forkjoin(true, true), "diffusion").supportsModelMethod("diffusion");
        assertFalse(reason.isEmpty());
        assertTrue(reason.contains("Fork"), reason);
    }

    @Test
    public void theFiveThatDoIntegrateItKeepTheFork() {
        // The converse, and the half that had to come BACK in C++: these five
        // answer the symmetric closed fork-join symmetrically, which is the one
        // property no approximation of a symmetric model may lose.
        Network model = forkjoin(true, true);
        String[] methods = new String[]{"statedep", "refined", "tbi", "mfq", "rmf"};
        for (String method : methods) {
            assertTrue(fluid(model, method).getMethodFeatureSet(method).inspectFeature("Fork"),
                    method);
            assertEquals("", fluid(model, method).supportsModelMethod(method), method);
        }
    }

    @Test
    public void theFluidForkJoinGateAndTheAnalyzerAreOnePredicate() {
        Network model = forkjoin(true, false);
        String reason = fluid(model, "dae").supportsModelMethod("dae");
        RuntimeException err = assertThrows(RuntimeException.class,
                () -> fluid(model, "dae").getAvgQLen());
        assertTrue(err.getMessage().contains("fork-join fixed point"), err.getMessage());
        assertTrue(reason.contains("fork-join fixed point"), reason);
    }

    @Test
    public void anOpenModelLosesDiffusionAndTbi() {
        String[] methods = new String[]{"diffusion", "tbi"};
        for (String method : methods) {
            String reason = fluid(mm1(), method).supportsModelMethod(method);
            assertFalse(reason.isEmpty(), method);
            assertTrue(reason.contains("OpenClass"), reason);
        }
    }

    @Test
    public void aClosedModelKeepsDiffusionAndTbi() {
        Network[] models = new Network[]{repairmen(), cqn2()};
        String[] methods = new String[]{"diffusion", "tbi"};
        for (Network model : models) {
            for (String method : methods) {
                assertEquals("", fluid(model, method).supportsModelMethod(method), method);
            }
        }
    }

    @Test
    public void tbiRefusesACacheStation() {
        assertFalse(fluid(repairmen(), "tbi").getMethodFeatureSet("tbi").inspectFeature("Cache"));
        assertTrue(SolverFluid.getFeatureSet().inspectFeature("Cache"));
    }

    @Test
    public void aTimeVaryingLimitNeedsAFiniteTimespan() {
        String[] methods = new String[]{"mol", "mtginf"};
        for (String method : methods) {
            SolverOptions options = new SolverOptions(SolverType.FLUID);
            options.method = method;
            options.timespan = new double[]{0.0, Double.POSITIVE_INFINITY};
            String reason = new SolverFluid(mm1(), options).supportsModelMethod(method);
            assertFalse(reason.isEmpty(), method);
            assertTrue(reason.contains("finite horizon"), reason);
            assertTrue(reason.contains("timespan"), reason);
        }
    }

    @Test
    public void aFiniteTimespanRestoresThem() {
        // The converse, and the whole point of the rule being on the OPTIONS:
        // the model never changed, only the horizon did.
        String[] methods = new String[]{"mol", "mtginf"};
        for (String method : methods) {
            SolverOptions options = new SolverOptions(SolverType.FLUID);
            options.method = method;
            options.timespan = new double[]{0.0, 10.0};
            assertEquals("", new SolverFluid(mm1(), options).supportsModelMethod(method), method);
        }
    }

    @Test
    public void aStationaryLimitIsNotGatedOnAHorizon() {
        // "ggisgi" and "tga" report a stationary point, so the horizon rule must
        // not touch them whatever options.timespan says.
        String[] methods = new String[]{"ggisgi.fluid", "ggingi.tga"};
        for (String method : methods) {
            SolverOptions options = new SolverOptions(SolverType.FLUID);
            options.method = method;
            options.timespan = new double[]{0.0, Double.POSITIVE_INFINITY};
            String reason = new SolverFluid(mm1(), options).supportsModelMethod(method);
            assertFalse(reason.contains("finite horizon"), reason);
        }
    }
}
