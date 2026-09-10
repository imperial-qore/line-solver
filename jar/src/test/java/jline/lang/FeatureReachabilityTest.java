/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Cox2;
import jline.lang.processes.Coxian;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Replayer;
import jline.lang.processes.Trace;
import jline.solvers.ldes.SolverLDES;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Registry entries that no model could ever set are dead surface: SUPPORTS
 * iterates the registry, so a name nothing emits gates nothing, and a solver
 * declaring it says nothing about what it accepts.
 *
 * Cox2, Trace and SchedStrategy_REF were three such names (see
 * _kb/06-solver-catalog.md): the first two were shadowed by the parent
 * distribution's getName(), the third by an empty case in
 * LayeredNetwork.getUsedLangFeatures. This test pins that they are now emitted,
 * AND that emitting them widened nothing -- a solver declaring only the general
 * name still accepts the model, through FeatureSet.generalizationOf.
 */
public class FeatureReachabilityTest {

    private static Network openModelWith(Distribution service) {
        Network model = new Network("featreach");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Exp(0.2));
        queue.setService(jobClass, service);
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    private static double[] samples() {
        return new double[]{0.5, 1.0, 1.5, 2.0};
    }

    @Test
    public void testTwoPhaseCoxianMarksCox2() {
        // MATLAB has no Cox2 OBJECT -- Cox2.fitMeanAndSCV returns a Coxian -- so
        // the phase count is the only reading of the entry the four codebases
        // share, and the one under which it is reachable at all.
        FeatureSet used = openModelWith(Cox2.fitMeanAndSCV(1.0, 2.0)).getUsedLangFeatures();
        assertTrue(used.inspectFeature("Cox2"), "a Cox2 must mark the Cox2 entry");
        assertFalse(used.inspectFeature("Coxian"), "the most specific name is the one marked");

        // A generic Coxian holding two phases is the same distribution and marks
        // the same entry; the class the user happened to call is not the test.
        Coxian twoPhase = new Coxian(Arrays.asList(2.0, 3.0), Arrays.asList(0.4, 1.0));
        assertEquals("Cox2", twoPhase.getFeatureName());
        assertTrue(openModelWith(twoPhase).getUsedLangFeatures().inspectFeature("Cox2"));
    }

    @Test
    public void testCoxianBeyondTwoPhasesStillMarksCoxian() {
        Coxian threePhase = new Coxian(Arrays.asList(2.0, 3.0, 4.0), Arrays.asList(0.3, 0.5, 1.0));
        assertEquals("Coxian", threePhase.getFeatureName());
        FeatureSet used = openModelWith(threePhase).getUsedLangFeatures();
        assertTrue(used.inspectFeature("Coxian"));
        assertFalse(used.inspectFeature("Cox2"));
    }

    @Test
    public void testTraceMarksTraceAndReplayerMarksReplayer() {
        FeatureSet usedTrace = openModelWith(new Trace(samples())).getUsedLangFeatures();
        assertTrue(usedTrace.inspectFeature("Trace"), "a Trace must mark the Trace entry");
        assertFalse(usedTrace.inspectFeature("Replayer"));

        FeatureSet usedReplayer = openModelWith(new Replayer(samples())).getUsedLangFeatures();
        assertTrue(usedReplayer.inspectFeature("Replayer"));
        assertFalse(usedReplayer.inspectFeature("Trace"));
    }

    @Test
    public void testGetNameIsUnchangedSoTheWireTypeIsUnchanged() {
        // getFeatureName exists precisely so that getName can stay put: it also
        // resolves the ProcessType and the JSON wire type, and moving it would
        // change what a saved model reloads as.
        assertEquals("Coxian", Cox2.fitMeanAndSCV(1.0, 2.0).getName());
        assertEquals("Replayer", new Trace(samples()).getName());
    }

    @Test
    public void testSpecializationFallsBackToItsGeneralization() {
        assertEquals("Coxian", FeatureSet.generalizationOf("Cox2"));
        assertEquals("Replayer", FeatureSet.generalizationOf("Trace"));
        // Nothing else may acquire a fallback silently: a fallback WIDENS what a
        // declared set accepts, so each one is a deliberate decision.
        assertNull(FeatureSet.generalizationOf("Coxian"));
        assertNull(FeatureSet.generalizationOf("Erlang"));
        assertNull(FeatureSet.generalizationOf("SchedStrategy_HOL"));

        FeatureSet used = new FeatureSet();
        used.setTrue("Cox2");
        FeatureSet general = new FeatureSet();
        general.setTrue("Coxian");
        assertTrue(FeatureSet.unsupportedFeatures(general, used).isEmpty(),
                "a solver declaring Coxian must still accept a two-phase Coxian");

        // One way only: declaring the SPECIAL case does not buy the general one.
        FeatureSet specific = new FeatureSet();
        specific.setTrue("Cox2");
        FeatureSet usedGeneral = new FeatureSet();
        usedGeneral.setTrue("Coxian");
        assertEquals(Arrays.asList("Coxian"), FeatureSet.unsupportedFeatures(specific, usedGeneral));

        // Neither declared is still a refusal, naming the specific feature.
        assertEquals(Arrays.asList("Cox2"), FeatureSet.unsupportedFeatures(new FeatureSet(), used));
    }

    @Test
    public void testEmittingTheSpecificNameRefusesNoSolverThatAcceptedBefore() {
        // The regression the fallback exists to prevent, stated on the featsets
        // that ship: every solver declaring Coxian (Replayer) must accept a
        // two-phase Coxian (a Trace), or the split silently narrowed it.
        List<FeatureSet> featsets = Arrays.asList(
                jline.solvers.mva.SolverMVA.getFeatureSet(),
                jline.solvers.nc.SolverNC.getFeatureSet(),
                jline.solvers.ctmc.SolverCTMC.getFeatureSet(),
                jline.solvers.ssa.SolverSSA.getFeatureSet(),
                jline.solvers.fluid.SolverFluid.getFeatureSet(),
                jline.solvers.mam.SolverMAM.getFeatureSet(),
                SolverLDES.getFeatureSet());
        FeatureSet cox2 = new FeatureSet();
        cox2.setTrue("Cox2");
        FeatureSet trace = new FeatureSet();
        trace.setTrue("Trace");
        for (FeatureSet featset : featsets) {
            if (featset.inspectFeature("Coxian")) {
                assertTrue(FeatureSet.unsupportedFeatures(featset, cox2).isEmpty(),
                        "a featset declaring Coxian must accept Cox2");
            }
            if (featset.inspectFeature("Replayer")) {
                assertTrue(FeatureSet.unsupportedFeatures(featset, trace).isEmpty(),
                        "a featset declaring Replayer must accept Trace");
            }
        }
    }

    @Test
    public void testReferenceTaskMarksSchedStrategyREF() {
        LayeredNetwork model = new LayeredNetwork("lqnfeat");
        Processor p = new Processor(model, "p1", 1, SchedStrategy.PS);
        Task t = new Task(model, "t1", 1, SchedStrategy.REF);
        t.on(p);
        Entry e = new Entry(model, "e1");
        e.on(t);
        Activity a = new Activity(model, "a1", new Exp(1.0));
        a.on(t);
        a.boundTo(e);

        FeatureSet used = model.getUsedLangFeatures();
        assertTrue(used.inspectFeature("SchedStrategy_REF"),
                "a reference task IS a discipline, and every LQN carries one");
        // The sole consumer of the LQN set already declared it, so marking it
        // rejects nothing that solved before.
        assertTrue(FeatureSet.unsupportedFeatures(SolverLDES.getLNFeatureSet(), used).isEmpty(),
                "SolverLDES declares SchedStrategy_REF, so the LQN must still pass");
    }

    @Test
    public void testTheLqnScanMarksTheConstructsItWalks() {
        // Host, Processor, Task, Entry, Activity, SyncCall and AsyncCall are all
        // declared by getLNFeatureSet and were emitted by nothing: the scan
        // walked each collection and marked only the scheduling discipline.
        LayeredNetwork model = new LayeredNetwork("lqnconstructs");
        Processor p = new Processor(model, "p1", 1, SchedStrategy.PS);
        Task client = new Task(model, "client", 1, SchedStrategy.REF);
        Task server = new Task(model, "server", 1, SchedStrategy.FCFS);
        client.on(p);
        server.on(p);
        Entry ce = new Entry(model, "ce");
        ce.on(client);
        Entry se = new Entry(model, "se");
        se.on(server);
        Activity ca = new Activity(model, "ca", new Exp(1.0));
        ca.on(client);
        ca.boundTo(ce);
        ca.synchCall(se, 1.0);
        Activity sa = new Activity(model, "sa", new Exp(2.0));
        sa.on(server);
        sa.boundTo(se);
        sa.repliesTo(se);

        FeatureSet used = model.getUsedLangFeatures();
        for (String feature : Arrays.asList("Host", "Processor", "Task", "Entry", "Activity",
                "SyncCall", "SchedStrategy_REF", "SchedStrategy_PS", "SchedStrategy_FCFS", "Exp")) {
            assertTrue(used.inspectFeature(feature), feature + " must be marked");
        }
        assertFalse(used.inspectFeature("AsyncCall"), "this model makes no asynchronous call");
        // Everything the scan can emit is declared by the only solver that reads
        // the set, so completing the scan refuses no LQN that solved before.
        assertTrue(FeatureSet.unsupportedFeatures(SolverLDES.getLNFeatureSet(), used).isEmpty(),
                "SolverLDES must still accept a plain two-task LQN");
    }

    @Test
    public void testActivityPrecedenceNamesAreEmitted() {
        LayeredNetwork model = new LayeredNetwork("lqnprec");
        Processor p = new Processor(model, "p1", 1, SchedStrategy.PS);
        Task t = new Task(model, "t1", 1, SchedStrategy.REF);
        t.on(p);
        Entry e = new Entry(model, "e1");
        e.on(t);
        Activity a0 = new Activity(model, "a0", new Exp(1.0));
        Activity b1 = new Activity(model, "b1", new Exp(1.0));
        Activity b2 = new Activity(model, "b2", new Exp(1.0));
        Activity c1 = new Activity(model, "c1", new Exp(1.0));
        a0.on(t);
        b1.on(t);
        b2.on(t);
        c1.on(t);
        a0.boundTo(e);
        t.addPrecedence(jline.lang.layered.ActivityPrecedence.AndFork(a0, Arrays.asList(b1, b2)));
        t.addPrecedence(jline.lang.layered.ActivityPrecedence.AndJoin(Arrays.asList(b1, b2), c1));

        FeatureSet used = model.getUsedLangFeatures();
        assertTrue(used.inspectFeature("ActivityPrecedence_POST_AND"), "the fork is a POST_AND");
        assertTrue(used.inspectFeature("ActivityPrecedence_PRE_AND"), "the join is a PRE_AND");
        assertTrue(FeatureSet.unsupportedFeatures(SolverLDES.getLNFeatureSet(), used).isEmpty(),
                "SolverLDES declares the AND precedences, so the LQN must still pass");
    }

    @Test
    public void testTheLqnFeatureSetIsBuildableAtAll() {
        // getLNFeatureSet declared CacheTask, ItemEntry and
        // ActivityPrecedence_POST_CACHE, none of which were registered, and
        // setTrue line_errors on an unknown name -- so building the set threw
        // and SolverLDES.supports(LayeredNetwork) could never run. That is also
        // why SchedStrategy_REF going unrecorded had gone unnoticed.
        FeatureSet ln = SolverLDES.getLNFeatureSet();
        assertTrue(ln.inspectFeature("CacheTask"));
        assertTrue(ln.inspectFeature("ItemEntry"));
        assertTrue(ln.inspectFeature("ActivityPrecedence_POST_CACHE"));
    }

    @Test
    public void testPrecedenceFeatureNamesAreRegistryNames() {
        // The helper returned "ActivityPrecedenceType_" names, which match no
        // registry entry, and took an argument type it could never be passed, so
        // every branch fell through to the throw.
        assertEquals("ActivityPrecedence_PRE_SEQ",
                jline.lang.constant.ActivityPrecedenceType.toFeature(
                        jline.lang.constant.ActivityPrecedenceType.PRE_SEQ));
        assertEquals("ActivityPrecedence_POST_CACHE",
                jline.lang.constant.ActivityPrecedenceType.toFeature(
                        jline.lang.constant.ActivityPrecedenceType.POST_CACHE));
        // POST_LOOP has no registry entry: a loop is modelled by its pseudo-task.
        assertEquals("", jline.lang.constant.ActivityPrecedenceType.toFeature(
                jline.lang.constant.ActivityPrecedenceType.POST_LOOP));
        FeatureSet registry = new FeatureSet();
        for (String type : Arrays.asList(
                jline.lang.constant.ActivityPrecedenceType.PRE_SEQ,
                jline.lang.constant.ActivityPrecedenceType.PRE_AND,
                jline.lang.constant.ActivityPrecedenceType.PRE_OR,
                jline.lang.constant.ActivityPrecedenceType.POST_SEQ,
                jline.lang.constant.ActivityPrecedenceType.POST_AND,
                jline.lang.constant.ActivityPrecedenceType.POST_OR,
                jline.lang.constant.ActivityPrecedenceType.POST_CACHE)) {
            // setTrue throws on an unregistered name, which is the assertion.
            registry.setTrue(jline.lang.constant.ActivityPrecedenceType.toFeature(type));
        }
    }

    @Test
    public void testHeterogeneousServerPoolsAreMarkedAndGated() {
        // A pooled station is NOT a station of the same total size: the pools
        // carry their own class compatibilities and their own rates. Every
        // solver that reads only sn.nservers therefore has to refuse rather
        // than flatten them, which it can only do once the pools are marked.
        Network model = new Network("hetfeat");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Exp(0.2));
        queue.setService(jobClass, new Exp(1.0));
        model.link(model.serialRouting(source, queue, sink));
        assertFalse(model.getUsedLangFeatures().inspectFeature("HeteroServers"),
                "a plain queue must not mark HeteroServers");

        queue.addServerType(new jline.lang.constant.ServerType("fast", 1));
        queue.addServerType(new jline.lang.constant.ServerType("slow", 2));
        FeatureSet used = model.getUsedLangFeatures();
        assertTrue(used.inspectFeature("HeteroServers"),
                "server pools must mark HeteroServers");

        // JMT serialises the pools and the LDES engine simulates them; the
        // analytical solvers read nservers and must refuse.
        assertTrue(jline.solvers.wrappers.jmt.SolverJMT.getFeatureSet()
                .inspectFeature("HeteroServers"));
        assertTrue(SolverLDES.getFeatureSet().inspectFeature("HeteroServers"));
        assertFalse(jline.solvers.mva.SolverMVA.getFeatureSet()
                .inspectFeature("HeteroServers"));
        assertFalse(jline.solvers.nc.SolverNC.getFeatureSet()
                .inspectFeature("HeteroServers"));
        assertFalse(jline.solvers.ctmc.SolverCTMC.getFeatureSet()
                .inspectFeature("HeteroServers"));

        assertTrue(FeatureSet.unsupportedFeatures(
                jline.solvers.wrappers.jmt.SolverJMT.getFeatureSet(), used).isEmpty());
        assertEquals(Arrays.asList("HeteroServers"), FeatureSet.unsupportedFeatures(
                jline.solvers.mva.SolverMVA.getFeatureSet(), used));
    }

    @Test
    public void testNonNormalDepartureDisciplineIsRefusedByEverySolver() {
        // A FIFO depository releases a served token only after the tokens that
        // entered service before it, so which output transitions are enabled
        // depends on the arrival order and not only on the marking. NO engine
        // in ANY codebase implements it, so no feature set declares it and the
        // model is refused rather than served as if it were Normal.
        Network model = new Network("depfeat");
        jline.lang.nodes.Place place = new jline.lang.nodes.Place(model, "P1");
        jline.lang.nodes.Transition transition =
                new jline.lang.nodes.Transition(model, "T1");
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, place);
        place.setService(jobClass, new Exp(1.0));
        assertFalse(model.getUsedLangFeatures().inspectFeature("DepartureDiscipline"),
                "a Normal depository must not mark DepartureDiscipline");

        place.setDepartureDiscipline(jobClass,
                jline.lang.constant.DepartureDiscipline.FIFO);
        FeatureSet used = model.getUsedLangFeatures();
        assertTrue(used.inspectFeature("DepartureDiscipline"),
                "a FIFO depository must mark DepartureDiscipline");

        List<FeatureSet> featsets = Arrays.asList(
                jline.solvers.wrappers.jmt.SolverJMT.getFeatureSet(),
                SolverLDES.getFeatureSet(),
                jline.solvers.ctmc.SolverCTMC.getFeatureSet(),
                jline.solvers.ssa.SolverSSA.getFeatureSet(),
                jline.solvers.mva.SolverMVA.getFeatureSet());
        for (FeatureSet featset : featsets) {
            assertFalse(featset.inspectFeature("DepartureDiscipline"),
                    "no solver implements a FIFO depository, so none may declare it");
            assertTrue(FeatureSet.unsupportedFeatures(featset, used)
                            .contains("DepartureDiscipline"),
                    "the refusal must name the offending construct");
        }
        // transition is part of the model shape the Place needs; referenced so
        // the local is not dead in a reader's eye
        assertEquals("T1", transition.getName());
    }
}
