package jline.solvers.auto;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.Pareto;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests of findSolver: the report of which solvers and solver methods can
 * analyze a model, why the others cannot, what kind of answer each returns and
 * which measures it can report.
 *
 * <p>WHAT IS ASSERTED, and why none of it is a value read back out of the
 * implementation.
 *
 * <p>1. THE PROJECTION IDENTITY. listValidMethods is defined as the runnable
 * rows of findSolver plus the method names that name no single method. That is the
 * point of the refactor -- one gate, two views -- so it is asserted directly
 * rather than trusted.
 *
 * <p>2. STRUCTURAL FACTS ABOUT THE MODELS. An M/M/1 has a product-form solution
 * and is one queueing station fed by a Source, so exact MVA and the QBD are
 * exact ON IT; a three-station closed network is still product form but is no
 * longer one queue, so the QBD is not. These are properties of the models.
 *
 * <p>3. THE REFUSAL REASONS NAME THE FEATURE. A report whose whole content is
 * the explanation is worthless if the explanation is "unsupported".
 */
public class FindSolverTest {

    /** Source -> FCFS Queue -> Sink, exponential throughout. */
    private static Network mm1() {
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oc = new OpenClass(model, "C1");
        source.setArrival(oc, new Exp(1.0));
        queue.setService(oc, new Exp(2.0));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /** The same shape with a Pareto service law: SolverNC has no Pareto. */
    private static Network mpareto() {
        Network model = new Network("mpareto");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oc = new OpenClass(model, "C1");
        source.setArrival(oc, new Exp(0.5));
        queue.setService(oc, new Pareto(2.5, 1.0));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /** Delay -> Queue1 -> Queue2, N = 4: a closed product-form network. */
    private static Network cqn() {
        Network model = new Network("cqn");
        Delay delay = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        ClosedClass cc = new ClosedClass(model, "C1", 4, delay);
        delay.setService(cc, new Exp(1.0));
        q1.setService(cc, new Exp(2.0));
        q2.setService(cc, new Exp(3.0));
        model.link(Network.serialRouting(delay, q1, q2));
        return model;
    }

    private static SolverCandidate row(List<SolverCandidate> rows, String method) {
        for (SolverCandidate r : rows) {
            if (r.method.equals(method)) {
                return r;
            }
        }
        return null;
    }

    @Test
    public void listValidMethodsIsTheRunnableRowsOfFindSolver() {
        Network model = mm1();
        SolverAUTO auto = new SolverAUTO(model);
        List<SolverCandidate> all = auto.findSolver("", true);
        Set<String> valid = new HashSet<String>(Arrays.asList(auto.listValidMethods()));
        assertFalse(all.isEmpty());

        int runnable = 0;
        for (SolverCandidate r : all) {
            if (r.runnable) {
                runnable++;
                // Every runnable pair is a method name a caller may ask AUTO for.
                assertTrue(valid.contains(r.method), r.method + " is runnable but not offered");
                assertTrue(valid.contains(r.solver));
                assertEquals("", r.reason);
            } else {
                // and a refused one is not offered.
                assertFalse(valid.contains(r.method), r.method + " is refused but still offered");
                assertFalse(r.reason.isEmpty(), r.method + " was refused without a reason");
            }
        }
        assertTrue(runnable > 0);

        // The default report is exactly the runnable half.
        List<SolverCandidate> runnableOnly = auto.findSolver();
        assertEquals(runnable, runnableOnly.size());
        for (SolverCandidate r : runnableOnly) {
            assertTrue(r.runnable);
        }

        // The selection intents name a ranking rather than an algorithm and are
        // listed whatever the model is.
        assertTrue(valid.contains("default"));
        assertTrue(valid.contains("exact"));
        assertTrue(valid.contains("bound"));
    }

    @Test
    public void aSelfQualifiedSpellingIsNotDoubled() {
        // SolverFluid declares both "dae" and "fluid.dae" so that its own gate
        // takes either; prefixing the family again would yield
        // "fluid.fluid.dae", a method name that resolves but names the same method
        // twice and would double every fluid row.
        List<SolverCandidate> rows = mm1().findSolver("", true);
        assertNotNull(row(rows, "fluid.dae"));
        assertNull(row(rows, "fluid.fluid.dae"));
        for (SolverCandidate r : rows) {
            assertFalse(r.method.contains(r.solver + "." + r.solver + "."));
        }
    }

    @Test
    public void modelEntryPointsAndTheirAliasesAgree() {
        Network model = mm1();
        List<SolverCandidate> t = model.findSolver();
        assertFalse(t.isEmpty());
        // findMethod and help are the same question asked in other words.
        assertEquals(t.size(), model.findMethod().size());
        assertEquals(t.size(), model.help().size());
        assertEquals(t.get(0).method, model.findMethod().get(0).method);
        assertEquals(t.get(0).method, model.help().get(0).method);
    }

    @Test
    public void exactnessIsClaimedOfTheModelNotTheAlgorithm() {
        Network open = mm1();
        assertTrue(open.hasProductFormSolution());
        List<SolverCandidate> a = open.findSolver("", true);
        // Exact MVA on a product-form model; the QBD on one queueing station
        // fed by a Source.
        assertEquals(SolverCandidate.CLASS_EXACT, row(a, "mva.exact").methodClass);
        assertEquals(SolverCandidate.CLASS_EXACT, row(a, "mam.default").methodClass);
        // A decomposition of a network into such queues approximates it.
        assertEquals(SolverCandidate.CLASS_APPROX, row(a, "mam.dec.source").methodClass);

        // Three stations: still product form, so MVA stays exact, but no longer
        // one queueing station, so the QBD is not.
        Network closed = cqn();
        List<SolverCandidate> b = closed.findSolver("", true);
        assertEquals(SolverCandidate.CLASS_EXACT, row(b, "mva.exact").methodClass);
        SolverCandidate mam = row(b, "mam.default");
        if (mam != null) {
            assertEquals(SolverCandidate.CLASS_APPROX, mam.methodClass);
        }
        // Every AMVA arm is an approximation whatever the model.
        assertEquals(SolverCandidate.CLASS_APPROX, row(b, "mva.amva").methodClass);
        // Bounds are what SolverBA is for, and a simulator is a simulator.
        for (SolverCandidate r : b) {
            if ("ba".equals(r.solver)) {
                assertEquals(SolverCandidate.CLASS_BOUND, r.methodClass);
            }
            if ("ssa".equals(r.solver) || "ldes".equals(r.solver)) {
                assertEquals(SolverCandidate.CLASS_SIMULATION, r.methodClass);
            }
        }
        // The CTMC generator is solved as written; "cftp.approx" says otherwise
        // in its own name.
        assertEquals(SolverCandidate.CLASS_EXACT, row(b, "ctmc.exact").methodClass);
        SolverCandidate cftp = row(b, "ctmc.cftp.approx");
        if (cftp != null) {
            assertEquals(SolverCandidate.CLASS_APPROX, cftp.methodClass);
        }
    }

    @Test
    public void aRefusalNamesTheFeatureThatCausedIt() {
        List<SolverCandidate> rows = mpareto().findSolver("", true);
        int refused = 0;
        for (SolverCandidate r : rows) {
            if (!r.runnable) {
                refused++;
                assertFalse(r.reason.isEmpty());
            }
        }
        assertTrue(refused > 0);
        // SolverNC has no Pareto in its feature set, so it must refuse and say so.
        SolverCandidate nc = row(rows, "nc.default");
        assertNotNull(nc);
        assertFalse(nc.runnable);
        assertTrue(nc.reason.contains("Pareto"), "the reason must name the feature: " + nc.reason);
    }

    @Test
    public void aMetricNarrowsTheReportToTheFamiliesThatAnswerIt() {
        Network model = mm1();
        List<SolverCandidate> byGroup = model.findSolver("cdf", false);
        List<SolverCandidate> byAccessor = model.findSolver("getCdfRespT", false);
        assertFalse(byGroup.isEmpty());
        // A group name and the accessor that returns it ask the same question.
        assertEquals(byGroup.size(), byAccessor.size());

        // MVA computes no passage-time law and must not appear; the simulators
        // and the transform solvers do.
        boolean sawMva = false;
        boolean sawLdes = false;
        for (SolverCandidate r : byGroup) {
            if ("mva".equals(r.solver)) {
                sawMva = true;
            }
            if ("ldes".equals(r.solver)) {
                sawLdes = true;
            }
            assertTrue(r.metrics.contains("cdf"));
        }
        assertFalse(sawMva);
        assertTrue(sawLdes);

        // Every family answers the mean measures, so "avg" narrows nothing away.
        assertEquals(model.findSolver().size(), model.findSolver("avg", false).size());

        // A name that is neither a group nor an accessor is a caller error, not
        // an empty answer that would read as "nothing can do this".
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                mm1().findSolver("nosuchmeasure", false);
            }
        });
    }

    @Test
    public void theMetricRegistryIsSelfConsistent() {
        Set<String> groups = new HashSet<String>(Arrays.asList(SolverAUTO.metricGroups()));
        // A group name maps to itself.
        for (String g : SolverAUTO.metricGroups()) {
            assertEquals(g, SolverAUTO.metricGroupOf(g));
        }
        // Every group a family declares is a registered one; a typo here would
        // silently hide the family from a caller asking for that measure.
        for (String fam : SolverAUTO.familyNames()) {
            String[] declared = SolverAUTO.familyMetrics(fam);
            for (String g : declared) {
                assertTrue(groups.contains(g), fam + " declares an unregistered measure " + g);
            }
            // Every family answers the mean measures, which is what a solver is for.
            assertTrue(Arrays.asList(declared).contains("avg"));
        }
        assertEquals("", SolverAUTO.metricGroupOf(""));
        assertEquals("", SolverAUTO.metricGroupOf("any"));
        assertEquals("avg", SolverAUTO.metricGroupOf("getAvgTable"));
        assertEquals("tranprob", SolverAUTO.metricGroupOf("getTranProbAggr"));
    }

    @Test
    public void theTableRendersEveryRow() {
        List<SolverCandidate> rows = mm1().findSolver();
        String table = SolverCandidate.toTable(rows);
        // One header line plus one line per row.
        int lines = 0;
        for (int i = 0; i < table.length(); i++) {
            if (table.charAt(i) == '\n') {
                lines++;
            }
        }
        assertEquals(rows.size() + 1, lines);
        assertTrue(table.startsWith("Solver"));
        for (SolverCandidate r : rows) {
            assertTrue(table.contains(r.method));
        }
        // No row of a runnable-only report ends in the padding of an empty reason.
        assertFalse(table.contains(" \n"));
        assertEquals("No solver method can analyze this model.\n",
                SolverCandidate.toTable(new java.util.ArrayList<SolverCandidate>()));
    }
}
