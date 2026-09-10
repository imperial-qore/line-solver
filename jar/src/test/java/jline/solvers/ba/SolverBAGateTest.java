package jline.solvers.ba;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.lang.ClosedClass;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Place;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Transition;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.Mode;
import jline.util.matrix.Matrix;

/**
 * The SolverBA method gate: what the report offers must be what runs.
 *
 * <p>{@code model.help()} / {@code findSolver} reports one row per (solver,
 * method) pair a model can run, and it builds that report from
 * {@code listValidMethods} filtered by {@code supportsModelMethod}. Both were
 * much weaker than the rules the ANALYZER enforces, so the report offered pairs
 * that raised on contact: measured on a two-class closed network, 30 of the 36
 * ba.* rows called runnable raised, and on a closed multiserver model 31 of 38
 * did.
 *
 * <p>The fix is one predicate, {@link SolverBA#methodRefusal}, asked by
 * runAnalyzer, by the gate and by the list, plus the per-method feature-set
 * deltas for the premises a feature NAME can state. These tests pin both
 * directions: a refused pair names the offending thing, and a model the bounds
 * ARE derived for keeps every method.
 */
public class SolverBAGateTest {

    // ---- model builders ----

    /** Single-class closed, no delay, one server each: the shape every family here is derived for. */
    private static Network cyclicDelayFree() {
        Network model = new Network("baGateCyclic");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "C", 3, q1);
        q1.setService(c, new Exp(1.0));
        q2.setService(c, new Exp(2.0));
        model.link(model.serialRouting(q1, q2));
        return model;
    }

    /** Two closed classes, DELAY-FREE so the class premise is the only one violated. */
    private static Network cyclicTwoClass() {
        Network model = new Network("baGate2Class");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, q1);
        ClosedClass c2 = new ClosedClass(model, "C2", 1, q1);
        q1.setService(c1, new Exp(1.0));
        q1.setService(c2, new Exp(2.0));
        q2.setService(c1, new Exp(2.0));
        q2.setService(c2, new Exp(3.0));
        model.link(model.serialRouting(q1, q2));
        return model;
    }

    /** Single-class closed with a three-server station: ssd, ldbcmp and auto survive. */
    private static Network multiserver() {
        Network model = new Network("baGateMultiserver");
        Delay d = new Delay(model, "D");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        q.setNumberOfServers(3);
        ClosedClass c = new ClosedClass(model, "C", 4, d);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(model.serialRouting(d, q));
        return model;
    }

    /** Single-class closed WITH a delay station: the think-time premise is the only one violated. */
    private static Network withDelay() {
        Network model = new Network("baGateDelay");
        Delay d = new Delay(model, "D");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", 3, d);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(model.serialRouting(d, q));
        return model;
    }

    /**
     * Two delay stations and nothing else: no queueing station at all, so the
     * ldbcmp bottleneck the open-network occupancy is built on does not exist.
     * The general form of the shape a Petri net presents, every Place being an
     * INF station -- and the shape that made methodDegenerate throw from inside
     * listValidMethods before it was guarded.
     */
    private static Network allDelayClosed() {
        Network model = new Network("baGateAllDelay");
        Delay d1 = new Delay(model, "D1");
        Delay d2 = new Delay(model, "D2");
        ClosedClass c = new ClosedClass(model, "C", 2, d1);
        d1.setService(c, new Exp(1.0));
        d2.setService(c, new Exp(2.0));
        model.link(model.serialRouting(d1, d2));
        return model;
    }

    /**
     * The fork-join stochastic Petri net of {@code SpnLpbndTest}: four Places
     * and three Transitions.
     *
     * <p>THE SHAPE THAT SEPARATES THE TWO INDEX SPACES. A Place extends Station
     * and so is a station AND a stateful node; a Transition extends ServiceNode
     * and so is stateful and NOT a station. This model therefore has 4 stations
     * against 7 stateful nodes, and sn.visits -- which SnRefreshVisits builds
     * (nstateful x nclasses) -- cannot be walked against sn.sched, sn.rates or
     * sn.stations, which are station-indexed. Reading it wrong threw
     * IndexOutOfBounds from inside listValidMethods.
     */
    private static Network forkJoinSpn(int ntokens) {
        Network model = new Network("baGateSpn");
        Place[] pl = new Place[4];
        for (int i = 0; i < 4; i++) {
            pl[i] = new Place(model, "P" + i);
        }
        Transition tf = new Transition(model, "Tf");
        Transition tj = new Transition(model, "Tj");
        Transition tb = new Transition(model, "Tb");
        ClosedClass jc = new ClosedClass(model, "C", ntokens, pl[0]);
        Mode mf = tf.addMode("f");
        tf.setDistribution(mf, new Exp(1.3));
        tf.setNumberOfServers(mf, 1);
        tf.setEnablingConditions(mf, jc, pl[0], 1);
        tf.setFiringOutcome(mf, jc, pl[1], 1);
        tf.setFiringOutcome(mf, jc, pl[2], 1);
        Mode mj = tj.addMode("j");
        tj.setDistribution(mj, new Exp(0.7));
        tj.setNumberOfServers(mj, 1);
        tj.setEnablingConditions(mj, jc, pl[1], 1);
        tj.setEnablingConditions(mj, jc, pl[2], 1);
        tj.setFiringOutcome(mj, jc, pl[3], 1);
        Mode mb = tb.addMode("b");
        tb.setDistribution(mb, new Exp(1.9));
        tb.setNumberOfServers(mb, 1);
        tb.setEnablingConditions(mb, jc, pl[3], 1);
        tb.setFiringOutcome(mb, jc, pl[0], 1);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jc, jc, pl[0], tf, 1.0);
        P.set(jc, jc, tf, pl[1], 1.0);
        P.set(jc, jc, tf, pl[2], 1.0);
        P.set(jc, jc, pl[1], tj, 1.0);
        P.set(jc, jc, pl[2], tj, 1.0);
        P.set(jc, jc, tj, pl[3], 1.0);
        P.set(jc, jc, pl[3], tb, 1.0);
        P.set(jc, jc, tb, pl[0], 1.0);
        model.link(P);
        for (int i = 0; i < 4; i++) {
            pl[i].setState(Matrix.singleton(i == 0 ? ntokens : 0));
        }
        return model;
    }

    private static Network mm1Open() {
        return openMM1(new Exp(1.0), new Exp(2.0), "baGateMM1");
    }

    /** Erlang SERVICE: what all three open families refuse at a queueing station. */
    private static Network mm1ErlangService() {
        return openMM1(new Exp(1.0), Erlang.fitMeanAndOrder(0.5, 3), "baGateErlSvc");
    }

    /**
     * Erlang SOURCE with exponential service. THE CONVERSE MODEL: 'snc' consumes
     * the arrival law and answers this one, so a gate that dropped the Erlang
     * feature outright would hide a bound the user could have had. 'bpt' and
     * 'bgt' must still be refused -- they read the mean alone and would bound the
     * Poisson system instead.
     */
    private static Network mm1ErlangSource() {
        return openMM1(Erlang.fitMeanAndOrder(1.0, 3), new Exp(2.0), "baGateErlSrc");
    }

    private static Network openMM1(jline.lang.processes.Distribution arrival,
                                   jline.lang.processes.Distribution service, String name) {
        Network model = new Network(name);
        Source source = new Source(model, "Source");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass c = new OpenClass(model, "C", 0);
        source.setArrival(c, arrival);
        q.setService(c, service);
        model.link(model.serialRouting(source, q, sink));
        return model;
    }

    /**
     * Delay + two PS queues, one class, N = 4: exactly the ldbcmp regime boundary
     * N == Qhat, where the bound degenerates to the trivial X &gt;= 0.
     */
    private static Network ldbcmpBoundary() {
        Network model = new Network("baGateLdbcmp");
        Delay d = new Delay(model, "D");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "C", 4, d);
        d.setService(c, new Exp(1.0));
        q1.setService(c, new Exp(2.0));
        q2.setService(c, new Exp(3.0));
        model.link(model.serialRouting(d, q1, q2));
        return model;
    }

    /**
     * The methods the REPORT offers: listValidMethods filtered by the gate,
     * which is exactly what findSolver builds its rows from. The two halves of
     * the gate live in different places on purpose -- the structural premises
     * narrow the list, and the ones a feature name states are applied by
     * supportsModelMethod, so a showAll report can still print the offending
     * FEATURE rather than dropping the row silently.
     */
    private static List<String> offered(Network model) {
        SolverBA probe = new SolverBA(model);
        List<String> out = new ArrayList<String>();
        String[] declared = probe.listValidMethods();
        for (int i = 0; i < declared.length; i++) {
            if (probe.supportsModelMethod(declared[i]).isEmpty()) {
                out.add(declared[i]);
            }
        }
        return out;
    }

    // ---- the gate names the offending thing ----

    @Test
    public void multiclassModelRefusesTheSingleClassFamiliesByName() {
        SolverBA probe = new SolverBA(cyclicTwoClass());
        String[] methods = {"aba.upper", "gb.lower", "pbh.upper", "harel.lower",
                "ssd.upper", "ldbcmp.lower", "auto.upper", "default"};
        for (int i = 0; i < methods.length; i++) {
            String reason = probe.supportsModelMethod(methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " was offered on a two-class model");
            assertTrue(reason.contains("single-class closed networks only"), reason);
        }
        // The reason names the method that will RUN, not the alias asked for.
        assertTrue(probe.supportsModelMethod("default").contains("'gb.upper'"));
    }

    @Test
    public void multiserverModelRefusesTheSingleServerFamiliesByName() {
        SolverBA probe = new SolverBA(multiserver());
        String[] methods = {"aba.upper", "bjb.lower", "gb.upper", "pbh.upper",
                "bjbk.lower", "default"};
        for (int i = 0; i < methods.length; i++) {
            String reason = probe.supportsModelMethod(methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " was offered on a multiserver model");
            assertTrue(reason.contains("multi-server stations"), reason);
            // 'ssd' IS the multiserver bound, so the refusal points at it.
            assertTrue(reason.contains("use 'ssd'"), reason);
        }
    }

    @Test
    public void multiclassChainFamiliesRefuseAMultiserverModelWithoutPointingAtSsd() {
        // mwba/cub/mbjb/looping survive a multiclass model but not a multiserver
        // one, and 'ssd' is single-class, so it is no alternative for them.
        SolverBA probe = new SolverBA(multiserver());
        String[] methods = {"mwba.upper", "cub.upper", "mbjb.lower", "looping.lower"};
        for (int i = 0; i < methods.length; i++) {
            String reason = probe.supportsModelMethod(methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " was offered on a multiserver model");
            assertTrue(reason.contains("multi-server stations"), reason);
            assertFalse(reason.contains("use 'ssd'"), reason);
        }
    }

    @Test
    public void openModelRefusesTheClosedFamiliesByName() {
        SolverBA probe = new SolverBA(mm1Open());
        String[] methods = {"aba.upper", "gb.upper", "mwba.upper", "cub.upper"};
        for (int i = 0; i < methods.length; i++) {
            String reason = probe.supportsModelMethod(methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " was offered on an open model");
            assertTrue(reason.contains("closed networks only"), reason);
        }
    }

    @Test
    public void aDelayStationRefusesTheThinkTimeFamiliesByFeatureName() {
        // This half of the gate is the FEATURE SET, not the structural
        // predicate: "does not accept a delay station" is expressible as
        // dropping SchedStrategy_INF, and a feature set can refuse a model for
        // HAVING a construct. The reason therefore names the feature.
        SolverBA probe = new SolverBA(withDelay());
        String[] refused = {"harel.upper", "harel.lower", "sb.upper", "sb.lower",
                "scb.upper", "sib.lower", "lr.upper", "lr"};
        for (int i = 0; i < refused.length; i++) {
            String reason = probe.supportsModelMethod(refused[i]);
            assertFalse(reason.isEmpty(), refused[i] + " was offered on a model with a delay");
            assertTrue(reason.contains("SchedStrategy_INF"), reason);
        }
        // and the ones that DO carry a think time are untouched
        String[] kept = {"gb.upper", "aba.lower", "pbh.upper", "default"};
        for (int i = 0; i < kept.length; i++) {
            assertEquals("", probe.supportsModelMethod(kept[i]),
                    kept[i] + " was dropped by the delay premise");
        }
    }

    @Test
    public void theOpenFamiliesDeclareTheMirrorPremise() {
        SolverBA probe = new SolverBA(mm1Open());
        assertFalse(probe.getMethodFeatureSet("bpt.lower").inspectFeature("ClosedClass"));
        assertFalse(probe.getMethodFeatureSet("bpt.lower").inspectFeature("SchedStrategy_INF"));
        assertTrue(probe.getMethodFeatureSet("bpt.lower").inspectFeature("OpenClass"));
    }

    // ---- one predicate, two callers ----

    @Test
    public void theRunRaisesTheSentenceTheGateReports() {
        Network[] models = {cyclicTwoClass(), multiserver(), multiserver(), mm1Open()};
        String[] methods = {"aba.upper", "gb.lower", "cub.upper", "mwba.upper"};
        for (int i = 0; i < models.length; i++) {
            SolverBA probe = new SolverBA(models[i]);
            String reason = probe.supportsModelMethod(methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " passed the gate");
            // Same predicate, same string: a caller gets one answer whichever
            // gate it meets first.
            assertEquals(reason,
                    SolverBA.methodRefusal(models[i].getStruct(false), methods[i]));
            final Network model = models[i];
            final String method = methods[i];
            assertThrows(Exception.class, new org.junit.jupiter.api.function.Executable() {
                public void execute() throws Throwable {
                    new SolverBA(model, method).getAvgTable();
                }
            }, method + " answered an inapplicable model instead of refusing");
        }
    }

    @Test
    public void theListIsAProjectionOfTheSamePredicate() {
        Network[] models = {cyclicDelayFree(), cyclicTwoClass(), multiserver(),
                withDelay(), mm1Open()};
        for (int i = 0; i < models.length; i++) {
            SolverBA probe = new SolverBA(models[i]);
            NetworkStruct sn = models[i].getStruct(false);
            String[] declared = probe.listValidMethods();
            for (int j = 0; j < declared.length; j++) {
                assertEquals("", SolverBA.methodRefusal(sn, declared[j]),
                        declared[j] + " is listed and structurally refused");
            }
        }
    }

    // ---- the bounds still answer the models they are derived for ----

    @Test
    public void aDelayFreeSingleClassClosedNetworkKeepsEveryMethod() {
        // Over-tightening is as bad as the leak. The only names this model may
        // lose are the three OPEN families and the four spnlp ones, both of
        // which were already gated before this change.
        List<String> listed = offered(cyclicDelayFree());
        List<String> missing = new ArrayList<String>(Arrays.asList(SolverBA.listAllMethods()));
        missing.removeAll(listed);
        List<String> expected = Arrays.asList("bpt.lower", "bgt.upper", "snc.upper",
                "spnlp.upper", "spnlp.lower", "spnlp.op.upper", "spnlp.op.lower");
        for (int i = 0; i < missing.size(); i++) {
            assertTrue(expected.contains(missing.get(i)),
                    missing.get(i) + " was dropped from a model the bounds are derived for");
        }
        for (int i = 0; i < expected.size(); i++) {
            assertFalse(listed.contains(expected.get(i)),
                    expected.get(i) + " is not a closed-network bound");
        }
    }

    @Test
    public void theMultiserverModelKeepsTheBoundsStatedForIt() {
        List<String> listed = offered(multiserver());
        String[] kept = {"ssd.upper", "ssd.lower", "ldbcmp.lower", "auto.upper", "auto.lower"};
        for (int i = 0; i < kept.length; i++) {
            assertTrue(listed.contains(kept[i]), kept[i] + " was dropped on a multiserver model");
        }
    }

    @Test
    public void theTwoClassModelKeepsTheMulticlassBounds() {
        List<String> listed = offered(cyclicTwoClass());
        String[] kept = {"mwba.upper", "mwba.lower", "cub.upper", "mbjb.lower",
                "looping.upper", "looping.lower"};
        for (int i = 0; i < kept.length; i++) {
            assertTrue(listed.contains(kept[i]), kept[i] + " was dropped on a two-class model");
        }
    }

    @Test
    public void everyOfferedMethodActuallyRuns() {
        // The audit condition, in the small: no offered ba.* pair may raise.
        Network[] models = {cyclicDelayFree(), cyclicTwoClass(), multiserver(),
                withDelay(), mm1Open()};
        for (int i = 0; i < models.length; i++) {
            List<String> list = offered(models[i]);
            assertFalse(list.isEmpty(), "no bound method offered for " + models[i].getName());
            for (int j = 0; j < list.size(); j++) {
                new SolverBA(models[i], list.get(j)).getAvgTable();
            }
        }
    }

    // ---- the exponential-service premise of the three open families ----

    @Test
    public void erlangServiceRefusesAllThreeOpenFamilies() {
        SolverBA probe = new SolverBA(mm1ErlangService());
        String[] methods = {"bpt.lower", "bgt.upper", "snc.upper"};
        for (int i = 0; i < methods.length; i++) {
            String reason = probe.supportsModelMethod(methods[i]);
            assertFalse(reason.isEmpty(), methods[i] + " was offered on Erlang service");
            assertTrue(reason.contains("Erlang") || reason.contains("exponential service"), reason);
        }
        assertTrue(offered(mm1ErlangService()).isEmpty());
    }

    @Test
    public void bptAndBgtDropEveryLawButExp() {
        // Registry-expressible: both read the mean alone, so a non-exponential
        // law ANYWHERE -- source included -- is silently bounded as if it were
        // Poisson rather than refused. Measured: swapping the Exp(1) source of an
        // M/M/1 for an Erlang of the same mean leaves bgt.upper at QLen 32.6667
        // and bpt.lower at 1.0, digit for digit.
        SolverBA probe = new SolverBA(mm1Open());
        String[] laws = {"APH", "Coxian", "Cox2", "Erlang", "HyperExp", "PH",
                "Det", "Lognormal", "Pareto", "Uniform", "Weibull"};
        String[] methods = {"bpt.lower", "bgt.upper"};
        for (int i = 0; i < methods.length; i++) {
            FeatureSet f = probe.getMethodFeatureSet(methods[i]);
            assertTrue(f.inspectFeature("Exp"), methods[i]);
            for (int j = 0; j < laws.length; j++) {
                assertFalse(f.inspectFeature(laws[j]), methods[i] + " still declares " + laws[j]);
            }
        }
    }

    @Test
    public void sncKeepsTheLawsBecauseItConsumesTheArrivalOne() {
        FeatureSet f = new SolverBA(mm1Open()).getMethodFeatureSet("snc.upper");
        String[] laws = {"Erlang", "Coxian", "APH", "PH", "HyperExp"};
        for (int i = 0; i < laws.length; i++) {
            assertTrue(f.inspectFeature(laws[i]), laws[i]);
        }
    }

    @Test
    public void anErlangSourceKeepsSncAndDropsBptAndBgt() {
        // The converse of the delta above, and the reason snc's rule is
        // structural: no feature name can say "Erlang at a Queue but not at a
        // Source", so dropping the law would refuse a model snc answers.
        Network model = mm1ErlangSource();
        SolverBA probe = new SolverBA(model);
        assertEquals("", probe.supportsModelMethod("snc.upper"),
                "snc was over-tightened by the Erlang source");
        assertEquals(Arrays.asList("snc.upper"), offered(model));
        String[] methods = {"bpt.lower", "bgt.upper"};
        for (int i = 0; i < methods.length; i++) {
            String reason = probe.supportsModelMethod(methods[i]);
            assertFalse(reason.isEmpty(), methods[i]);
            assertTrue(reason.contains("Erlang"), reason);
        }
        // ... and it really does run, which is what makes the refusal of the
        // other two a judgement about correctness rather than about coverage.
        new SolverBA(model, "snc.upper").getAvgTable();
    }

    @Test
    public void theSncRuleIsStructuralAndSkipsTheSource() {
        NetworkStruct svc = mm1ErlangService().getStruct(false);
        NetworkStruct src = mm1ErlangSource().getStruct(false);
        assertTrue(SolverBA.methodRefusal(svc, "snc.upper").contains("requires exponential service"),
                SolverBA.methodRefusal(svc, "snc.upper"));
        assertEquals("", SolverBA.methodRefusal(src, "snc.upper"));
        // bpt and bgt carry no structural rule at all: theirs is the feature set.
        assertEquals("", SolverBA.methodRefusal(svc, "bpt.lower"));
        assertEquals("", SolverBA.methodRefusal(svc, "bgt.upper"));
    }

    // ---- a degenerate bound is not offered, but is still answered ----

    @Test
    public void ldbcmpIsWithheldAtItsRegimeBoundary() {
        Network model = ldbcmpBoundary();
        NetworkStruct sn = model.getStruct(false);
        String why = SolverBA.methodDegenerate(sn, "ldbcmp.lower");
        assertTrue(why.contains("Qhat=4.0000"), why);
        assertTrue(why.contains("N=4"), why);
        assertTrue(why.contains("trivial bound X >= 0"), why);
        // The APPLICABILITY predicate stays silent: the model is inside the
        // method's domain, which is exactly why this is a separate question.
        assertEquals("", SolverBA.methodRefusal(sn, "ldbcmp.lower"));
        assertEquals(why, new SolverBA(model).supportsModelMethod("ldbcmp.lower"));
        assertFalse(offered(model).contains("ldbcmp.lower"));
    }

    @Test
    public void theRunStillAnswersADegenerateBoundWhenNamed() {
        // Not a refusal: X >= 0 IS a lower bound, so the run publishes it. What
        // changed is that nothing offers it.
        new SolverBA(ldbcmpBoundary(), "ldbcmp.lower").getAvgTable();
    }

    @Test
    public void ldbcmpSurvivesWhereItsBoundSaysSomething() {
        Network[] models = {withDelay(), multiserver()};
        for (int i = 0; i < models.length; i++) {
            assertEquals("", SolverBA.methodDegenerate(models[i].getStruct(false), "ldbcmp.lower"));
            assertTrue(offered(models[i]).contains("ldbcmp.lower"));
        }
    }

    @Test
    public void aModelWithNoQueueingStationIsAnsweredNotThrown() {
        // A PREDICATE MUST NOT THROW. methodDegenerate is asked once per name by
        // listValidMethods, which runs it before any later sieve has dropped
        // anything, so a model whose every station is an INF server reached
        // Pfqn_ldbcmp with an empty demand vector and took the whole listing
        // down with it.
        Network model = allDelayClosed();
        String why = SolverBA.methodDegenerate(model.getStruct(false), "ldbcmp.lower");
        assertTrue(why.contains("no queueing station"), why);
        List<String> listed = Arrays.asList(new SolverBA(model).listValidMethods());
        assertFalse(listed.isEmpty());
        assertFalse(listed.contains("ldbcmp.lower"));
    }

    @Test
    public void aPetriNetIsAnsweredNotThrown() {
        // THE TWO INDEX SPACES. sn.visits is stateful-indexed while sn.sched,
        // sn.rates and sn.stations are station-indexed, and on this net that is
        // 7 rows against 4 stations. Walking the visit rows against the station
        // space threw "Index 4 out of bounds for length 4" from inside
        // listValidMethods, which is what SpnLpbndTest caught.
        Network model = forkJoinSpn(3);
        NetworkStruct sn = model.getStruct(false);
        assertTrue(sn.nstations != sn.nstateful, "the net must separate the two spaces");
        String why = SolverBA.methodDegenerate(sn, "ldbcmp.lower");
        assertTrue(why.contains("no queueing station"), why);
        // ... and the listing completes, offering exactly the marking-indexed
        // family and nothing else.
        List<String> listed = Arrays.asList(new SolverBA(model).listValidMethods());
        assertEquals(4, listed.size(), "got " + listed);
        assertTrue(listed.containsAll(Arrays.asList("spnlp.upper", "spnlp.lower",
                "spnlp.op.upper", "spnlp.op.lower")), "got " + listed);
    }

    @Test
    public void noOtherMethodHasADegenerateRegime() {
        Network[] models = {cyclicDelayFree(), cyclicTwoClass(), multiserver(), withDelay(),
                mm1Open(), ldbcmpBoundary(), allDelayClosed(), forkJoinSpn(3)};
        String[] all = SolverBA.listAllMethods();
        for (int i = 0; i < models.length; i++) {
            NetworkStruct sn = models[i].getStruct(false);
            for (int j = 0; j < all.length; j++) {
                if ("ldbcmp.lower".equals(all[j])) {
                    continue;
                }
                assertEquals("", SolverBA.methodDegenerate(sn, all[j]), all[j]);
            }
        }
    }
}
