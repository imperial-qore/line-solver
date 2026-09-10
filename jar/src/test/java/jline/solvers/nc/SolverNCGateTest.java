/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

/**
 * The per-method support gate of SolverNC: what the report offers must be what runs.
 *
 * <p>{@code model.help()} / {@code findSolver} reports one row per (solver, method) pair,
 * and it builds the nc rows by asking {@link SolverNC#supportsModelMethod} -- the same
 * gate SolverAUTO applies before delegating and the same one {@code listValidMethods}
 * projects. Until the rules below reached it the gate was far weaker than what the
 * ANALYZER enforces, so the report offered pairs that then threw: measured on a plain
 * M/M/1, 15 of the 42 nc.* rows called runnable raised, and on a single-class closed
 * network 10 did.</p>
 *
 * <p>TWO ROUTES REACH THE REPORT, and both are pinned here. What the feature registry
 * CAN name rides in {@link SolverNC#methodFeatureSet}: a closed population for 'is' and
 * the six load-dependent evaluators (drop OpenClass), no think time for 'divdiff' (drop
 * SchedStrategy_INF). What it CANNOT -- requires a cache, requires state-dependent
 * routing, requires a loss network, requires exactly two stations, requires normal
 * usage -- is {@link SolverNC#ncMethodRefusal}, ONE predicate that runAnalyzer throws on
 * and that the gate asks on the way in, so the two cannot drift apart.</p>
 *
 * <p>Each case asserts the refusal AND its converse: a model the method IS derived for
 * must keep it. Over-tightening hides an answer the caller could have had, which is the
 * same defect with the sign flipped.</p>
 *
 * <p>The MATLAB twin is matlab/src/solvers/NC/nc_method_refusal.m, the python twin
 * python/tests/test_gate_nc.py and the C++ twin cpp/tests/test_gate_nc.cpp.</p>
 */
public class SolverNCGateTest {

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    // ---- model builders ----

    /** Source -> FCFS Queue -> Sink, one OPEN class. */
    private static Network mm1() {
        Network model = new Network("ncGateMM1");
        Source source = new Source(model, "S");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "K");
        OpenClass c = new OpenClass(model, "C");
        source.setArrival(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(model.serialRouting(source, q, sink));
        return model;
    }

    /**
     * Delay -> FCFS Queue, one closed class, N = 3. Saturated: the queue's offered load
     * 1.5 exceeds its unit saturation rate, so the model is NOT in normal usage and
     * PANACEA's asymptotic expansion does not apply to it.
     */
    private static Network repairmen() {
        Network model = new Network("ncGateRepairmen");
        Delay d = new Delay(model, "D");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", 3, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(model.serialRouting(d, q));
        return model;
    }

    /**
     * The same closed shape with a three-server queue, N = 4: mu(n) = min(n,3) is 3 at
     * saturation against an offered load of 2, so this one IS in normal usage and
     * 'panald' must survive the gate.
     */
    private static Network multiserver() {
        Network model = new Network("ncGateMultiserver");
        Delay d = new Delay(model, "D");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        q.setNumberOfServers(3);
        ClosedClass c = new ClosedClass(model, "C", 4, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(model.serialRouting(d, q));
        return model;
    }

    /**
     * The repairmen shape with a rate LATTICE on the queue, N = 4. It matters because
     * Pfqn_ncld evaluates "pana" and "panald" with the same Pfqn_panaceald, so on a
     * load-dependent model the load-INDEPENDENT name reaches the load-dependent expansion.
     * Saturation rate 1.9 against an offered load of 2, so this one is NOT in normal usage
     * and both names must be refused; every other nc method on it runs.
     */
    private static Network loaddep() {
        Network model = new Network("ncGateLoadDep");
        Delay d = new Delay(model, "D");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", 4, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        Matrix alpha = new Matrix(1, 4);
        alpha.set(0, 0, 1.0);
        alpha.set(0, 1, 1.6);
        alpha.set(0, 2, 1.8);
        alpha.set(0, 3, 1.9);
        q.setLoadDependence(alpha);
        model.link(model.serialRouting(d, q));
        return model;
    }

    /**
     * Delay -> Q1 -> Q2, one closed class, N = 4: TWO queueing stations, which is what the
     * mmint2/gleint/comomld recursions are not stated for.
     */
    private static Network cqn3() {
        Network model = new Network("ncGateCqn3");
        Delay d = new Delay(model, "D");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "C", 4, d, 0);
        d.setService(c, new Exp(1.0));
        q1.setService(c, new Exp(2.0));
        q2.setService(c, new Exp(3.0));
        model.link(model.serialRouting(d, q1, q2));
        return model;
    }

    /**
     * Source -> Q1 -> Q2 -> Sink: two queueing stations but NO closed population, so
     * Pfqn_nc answers with the exact open formulas before its method switch and the
     * single-station rule must stay inactive.
     */
    private static Network openTandem() {
        Network model = new Network("ncGateOpenTandem");
        Source source = new Source(model, "S");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "K");
        OpenClass c = new OpenClass(model, "C");
        source.setArrival(c, new Exp(0.5));
        q1.setService(c, new Exp(2.0));
        q2.setService(c, new Exp(3.0));
        model.link(model.serialRouting(source, q1, q2, sink));
        return model;
    }

    /**
     * Queue1 -> Queue2, one closed class, NO delay station: zero think time, which is the
     * shape the divided-difference closed form 'divdiff' is derived for (Casale,
     * SIGMETRICS 2017, Eqs. 15-16).
     */
    private static Network cyclicDelayFree() {
        Network model = new Network("ncGateCyclic");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "C", 2, q1, 0);
        q1.setService(c, new Exp(1.0));
        q2.setService(c, new Exp(2.0));
        model.link(model.serialRouting(q1, q2));
        return model;
    }

    /** The gate the report and the ranking both ask: "" when the pair may run. */
    private static String gate(Network model, String method) {
        return new SolverNC(model, method).supportsModelMethod(method);
    }

    /** The methods the report offers: listValidMethods filtered by the gate. */
    private static List<String> offered(Network model) {
        SolverNC probe = new SolverNC(model);
        List<String> out = new ArrayList<String>();
        String[] declared = probe.listValidMethods();
        for (int i = 0; i < declared.length; i++) {
            if (gate(model, declared[i]).isEmpty()) {
                out.add(declared[i]);
            }
        }
        return out;
    }

    /**
     * The predicate asked the RUN's question: "" when the reference performs the method
     * by name, even where the report declines to offer it.
     */
    private static String runQuestion(Network model, String method) {
        SolverNC probe = new SolverNC(model, method);
        return SolverNC.ncMethodRefusal(model.getStruct(false), method, probe.getOptions(), false);
    }

    /** The message the ANALYZER refused with, or "" when it did not refuse. */
    private static String analyzerRefusal(Network model, String method) {
        try {
            new SolverNC(model, method).getAvgTable();
            return "";
        } catch (Exception e) {
            String msg = e.getMessage();
            return msg == null ? e.getClass().getName() : msg;
        }
    }

    // ---- what the registry CANNOT name: ncMethodRefusal ----

    @Test
    public void theCacheTokensNameTheMissingCache() {
        String[] methods = {"rayint", "spm"};
        for (int i = 0; i < methods.length; i++) {
            String reason = gate(mm1(), methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " was offered on a model with no Cache");
            assertTrue(reason.contains("Cache node"), reason);
            assertTrue(reason.contains("SPM saddle point"), reason);
        }
    }

    @Test
    public void theLossNetworkTokensNameTheMissingRegion() {
        String[] methods = {"ms", "erlangfp"};
        for (int i = 0; i < methods.length; i++) {
            String reason = gate(repairmen(), methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " was offered off a loss network");
            assertTrue(reason.contains("loss network"), reason);
        }
    }

    @Test
    public void recNamesBothRoutesItHas() {
        String reason = gate(repairmen(), "rec");
        assertFalse(reason.isEmpty(), "rec was offered on a plain queueing network");
        assertTrue(reason.contains("MDD-rec"), reason);
        assertTrue(reason.contains("Petri net"), reason);
    }

    @Test
    public void sdrNamesTheRoutingTheModelDoesNotDeclare() {
        String[] methods = {"sdr", "sdr.mva"};
        for (int i = 0; i < methods.length; i++) {
            String reason = gate(mm1(), methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " was offered without sdr");
            assertTrue(reason.contains("state-dependent routing"), reason);
        }
    }

    @Test
    public void morrisonNamesTheShapeItIsDerivedFor() {
        String reason = gate(repairmen(), "morrison");
        assertFalse(reason.isEmpty(), "morrison was offered on a model with no DPS station");
        assertTrue(reason.contains("DPS"), reason);
        assertTrue(reason.contains("two stations"), reason);
    }

    @Test
    public void panaceaLdNamesNormalUsageOnASaturatedClosedModel() {
        String reason = gate(repairmen(), "panald");
        assertFalse(reason.isEmpty(), "panald was offered outside normal usage");
        assertTrue(reason.contains("normal usage"), reason);
    }

    @Test
    public void panaceaIsRefusedOnALoadDependentModelOutsideNormalUsage() {
        // Pfqn_ncld's case label is {"pana", "panald"}: on a rate lattice the
        // load-independent NAME is evaluated by the load-dependent kernel, so it throws the
        // panald refusal. The gate has to know that aliasing.
        String reason = gate(loaddep(), "pana");
        assertFalse(reason.isEmpty(), "pana was offered on a load-dependent model "
                + "outside normal usage");
        assertTrue(reason.contains("normal usage"), reason);
        assertTrue(reason.contains("evaluates it as 'panald'"), reason);
        assertFalse(analyzerRefusal(loaddep(), "pana").isEmpty(),
                "the analyzer accepted a pair the gate refused");
    }

    @Test
    public void theLoadDependentModelKeepsEveryOtherMethod() {
        // The converse, and the whole point of not gating "pana" off the lattice.
        String[] kept = {"default", "exact", "ca", "clw", "comom", "comomld", "le", "ble",
                "mmint2", "gleint", "kt", "bkt", "lekt", "cub", "rd", "nrl", "nrp", "nre", "propfair",
                "ger", "rgf"};
        Network model = loaddep();
        for (int i = 0; i < kept.length; i++) {
            assertTrue(gate(model, kept[i]).isEmpty(), kept[i] + ": " + gate(model, kept[i]));
        }
    }

    @Test
    public void panaceaIsLeftAloneOffTheLoadDependentRoute() {
        // Without a rate lattice "pana" takes its own Pfqn_nc arm, which warns and
        // returns an empty constant rather than throwing, so the gate must not refuse it
        // there -- that would be over-tightening on a path this change is not about.
        assertTrue(gate(repairmen(), "pana").isEmpty(), gate(repairmen(), "pana"));
        assertTrue(gate(multiserver(), "pana").isEmpty(), gate(multiserver(), "pana"));
    }

    @Test
    public void theSingleStationRecursionsNameTheStationCount() {
        // Pfqn_nc states "a model with a delay and a single queueing station" for
        // mmint2/gleint, and Pfqn_comomrm_ld raises "accepts at most a single queueing
        // station". Neither is a feature name: it is a COUNT.
        Network model = cqn3();
        String[] methods = {"mmint2", "gleint", "comomld"};
        for (int i = 0; i < methods.length; i++) {
            String reason = gate(model, methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " was offered on two queueing stations");
            assertTrue(reason.contains("single queueing station"), reason);
            assertTrue(reason.contains("has 2"), reason);
        }
    }

    /**
     * THE RULING (2026-07-25, reaffirmed when this gate was added): the report answers
     * "should this be offered" and the run answers "what does the reference do".
     * Pfqn_nc answers mmint2/gleint outside their shape with an empty constant and a
     * ZERO TABLE (pfqn_nc.m case {'mmint2','gleint'}: lG = [] and return,
     * unconditionally), so a caller who names the method keeps that answer while
     * model.help() stops offering it.
     *
     * <p>comomld is NOT in that bucket: Pfqn_comomrm_ld raises natively, so it is
     * refused on both paths.</p>
     */
    @Test
    public void mmint2AndGleintAreGatedForTheReportOnly() {
        Network model = cqn3();
        String[] reportOnly = {"mmint2", "gleint"};
        for (int i = 0; i < reportOnly.length; i++) {
            assertFalse(gate(model, reportOnly[i]).isEmpty(),
                    reportOnly[i] + " must not be offered by the report");
            assertTrue(runQuestion(model, reportOnly[i]).isEmpty(),
                    reportOnly[i] + " must stay runnable when asked by name");
        }
        assertFalse(gate(model, "comomld").isEmpty());
        assertFalse(runQuestion(model, "comomld").isEmpty(),
                "comomld raises natively and must be refused on both paths");
    }

    @Test
    public void theSingleStationRecursionsSurviveOnOneQueueingStation() {
        String[] methods = {"mmint2", "gleint", "comomld"};
        for (int i = 0; i < methods.length; i++) {
            assertTrue(gate(repairmen(), methods[i]).isEmpty(),
                    methods[i] + ": " + gate(repairmen(), methods[i]));
        }
    }

    @Test
    public void theStationCountRuleIsInactiveWithoutAClosedPopulation() {
        // An open network never reaches Pfqn_nc's method switch, so a two-queue OPEN tandem
        // runs these names correctly and must keep them.
        String[] methods = {"mmint2", "gleint"};
        for (int i = 0; i < methods.length; i++) {
            assertTrue(gate(openTandem(), methods[i]).isEmpty(),
                    methods[i] + ": " + gate(openTandem(), methods[i]));
        }
    }

    @Test
    public void panaceaLdSurvivesWhereTheExpansionDoesApply() {
        // The converse: a saturation rate of 3 against an offered load of 2 IS normal
        // usage, so the row must stay and the run must go through.
        assertTrue(gate(multiserver(), "panald").isEmpty(),
                gate(multiserver(), "panald"));
        assertTrue(analyzerRefusal(multiserver(), "panald").isEmpty(),
                analyzerRefusal(multiserver(), "panald"));
    }

    // ---- what the registry CAN name: methodFeatureSet ----

    @Test
    public void explicitRefusesAThinkTimeByFeatureName() {
        String reason = gate(repairmen(), "divdiff");
        assertFalse(reason.isEmpty(), "divdiff was offered on a model with a think time");
        assertTrue(reason.contains("SchedStrategy_INF"), reason);
    }

    @Test
    public void explicitSurvivesOnADelayFreeClosedModel() {
        Network model = cyclicDelayFree();
        assertTrue(gate(model, "divdiff").isEmpty(), gate(model, "divdiff"));
        assertTrue(analyzerRefusal(cyclicDelayFree(), "divdiff").isEmpty(),
                analyzerRefusal(cyclicDelayFree(), "divdiff"));
        FeatureSet fs = SolverNC.methodFeatureSet("divdiff");
        assertFalse(fs.inspectFeature("SchedStrategy_INF"),
                "divdiff must not declare an infinite server");
        assertTrue(SolverNC.methodFeatureSet("default").inspectFeature("SchedStrategy_INF"),
                "the base envelope keeps the infinite server");
    }

    @Test
    public void theClosedPopulationMethodsRefuseAnOpenChainByFeatureName() {
        String[] methods = {"is", "rd", "nrp", "nrl", "nre", "comomld", "panald"};
        for (int i = 0; i < methods.length; i++) {
            String reason = gate(mm1(), methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " was offered on an open model");
            assertTrue(reason.contains("OpenClass"), methods[i] + ": " + reason);
        }
    }

    @Test
    public void theClosedPopulationMethodsSurviveOnAClosedModel() {
        String[] methods = {"is", "rd", "nrp", "nrl", "nre", "comomld"};
        for (int i = 0; i < methods.length; i++) {
            String reason = gate(repairmen(), methods[i]);
            assertTrue(reason.isEmpty(), methods[i] + ": " + reason);
        }
    }

    // ---- one predicate, two callers ----

    @Test
    public void everyRefusedPairIsRefusedByTheAnalyzerToo() {
        Network[] models = {mm1(), repairmen(), cyclicDelayFree(), loaddep(), cqn3(),
                openTandem()};
        for (int mi = 0; mi < models.length; mi++) {
            String[] declared = new SolverNC(models[mi]).listValidMethods();
            for (int i = 0; i < declared.length; i++) {
                if (!gate(models[mi], declared[i]).isEmpty()) {
                    if (runQuestion(models[mi], declared[i]).isEmpty()) {
                        // A REPORT-ONLY refusal, and the ruling says so: the reference
                        // performs this one by name (a warning and a zero table), which is
                        // exactly why the report declines to offer it.
                        continue;
                    }
                    String msg = analyzerRefusal(models[mi], declared[i]);
                    assertFalse(msg.isEmpty(), models[mi].getName() + "/" + declared[i]
                            + ": the gate refused but the run did not");
                }
            }
        }
    }

    @Test
    public void theModelKeepsTheMethodsItCanGenuinelyRun() {
        // A closed product-form network must not lose its normalizing-constant methods
        // to this change; these are the ones it has always answered with.
        List<String> ok = offered(repairmen());
        String[] kept = {"default", "exact", "ca", "comom", "le", "ble", "mmint2",
                "gleint", "pana", "propfair", "cub", "kt", "bkt", "lekt", "clw"};
        for (int i = 0; i < kept.length; i++) {
            assertTrue(ok.contains(kept[i]), kept[i] + " was dropped from a closed "
                    + "product-form network");
        }
    }
}
