package jline.solvers.ln;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

/**
 * Method 'flat.ph': the squashed layering with the composed phase-type encoding.
 *
 * <p>ONE submodel holds a station for every processor and every called task, and
 * a caller task is one closed class that visits each server it uses once per
 * invocation, carrying there the composed law of the demand it places on that
 * server. It is the same composition 'srvn.ph' performs; what changes is that
 * the servers contend inside one network instead of meeting each other through
 * surrogate delays, so the client delay keeps only the think times.</p>
 *
 * <p>The goldens are MATLAB's, the reference implementation. Against lqsim they
 * read +7.16 per cent (two tasks), -0.32 (three tasks) and +0.87 (a three-deep
 * chain), which is the accuracy band of 'flat.cs' and 'srvn.ph' on the same
 * models and better than lqns on two of the three. Twin of the MATLAB
 * @SolverLN/buildLayersPH.m flat branch.</p>
 */
class SolverLNFlatPHTest {

    /**
     * The layer solve is an AMVA fixed point and the outer LN loop is another, so
     * a golden taken from MATLAB is met to the band the two codebases' iterations
     * already differ by, not to machine precision. Measured here: 3e-6 on the two
     * task models and 2.3e-4 on the three-deep chain, against the 7.7e-4 the
     * pre-existing 'srvn.cs' encoding shows on the same kind of corpus.
     */
    private static final double TOL = 1e-3;

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    private static double rowOf(LayeredNetworkAvgTable t, String nodeName, List<Double> col) {
        List<String> names = t.getNodeNames();
        int i = names.indexOf(nodeName);
        assertTrue(i >= 0, "row " + nodeName + " not in " + names);
        return col.get(i);
    }

    private static SolverLN flatPH(LayeredNetwork model) {
        SolverOptions o = new SolverOptions();
        o.method = "flat.ph";
        o.verbose = VerboseLevel.SILENT;
        return new SolverLN(model, o);
    }

    /** T1 -> T2, a single call level with a three-thread callee. */
    private static LayeredNetwork twoTasks() {
        LayeredNetwork m = new LayeredNetwork("twotasks");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.PS);
        Processor p2 = new Processor(m, "P2", 1, SchedStrategy.PS);
        Task t1 = new Task(m, "T1", 5, SchedStrategy.REF).on(p1).setThinkTime(new Exp(1));
        Task t2 = new Task(m, "T2", 3, SchedStrategy.FCFS).on(p2);
        Entry e1 = new Entry(m, "E1").on(t1);
        Entry e2 = new Entry(m, "E2").on(t2);
        new Activity(m, "A1", new Exp(2)).on(t1).boundTo(e1).synchCall(e2, 2);
        new Activity(m, "A2", new Exp(3)).on(t2).boundTo(e2).repliesTo(e2);
        return m;
    }

    /** T1 -> T2 -> T3, so the middle task is at once a server and a caller. */
    private static LayeredNetwork serialChain() {
        LayeredNetwork m = new LayeredNetwork("serialchain");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.PS);
        Processor p2 = new Processor(m, "P2", 1, SchedStrategy.PS);
        Processor p3 = new Processor(m, "P3", 1, SchedStrategy.PS);
        Task t1 = new Task(m, "T1", 4, SchedStrategy.REF).on(p1).setThinkTime(new Exp(1));
        Task t2 = new Task(m, "T2", 2, SchedStrategy.FCFS).on(p2);
        Task t3 = new Task(m, "T3", 1, SchedStrategy.FCFS).on(p3);
        Entry e1 = new Entry(m, "E1").on(t1);
        Entry e2 = new Entry(m, "E2").on(t2);
        Entry e3 = new Entry(m, "E3").on(t3);
        new Activity(m, "A1", new Exp(5)).on(t1).boundTo(e1).synchCall(e2, 1);
        new Activity(m, "A2", new Exp(5)).on(t2).boundTo(e2).synchCall(e3, 1).repliesTo(e2);
        new Activity(m, "A3", new Exp(5)).on(t3).boundTo(e3).repliesTo(e3);
        return m;
    }

    @Test
    void flatPHBuildsOneLayerHoldingEveryServer() {
        SolverLN s = flatPH(serialChain());
        s.getAvgTable();
        // the label reports what was BUILT, and the alias 'flat' does not reach it
        assertEquals("flat.ph", s.lnmethod);
        // ONE submodel, not one per server
        assertEquals(1, s.nlayers);
        // three processors and two called tasks, all stations of that one model
        assertEquals(3, s.getEnsemble().get(0).getAttribute().getHostStations().size());
        assertEquals(2, s.getEnsemble().get(0).getAttribute().getTaskStations().size());
        // one closed class per caller task, not one per entry, activity and call
        assertEquals(3, s.getEnsemble().get(0).getNumberOfClasses());
    }

    @Test
    void twoTasksMatchesMatlab() {
        LayeredNetworkAvgTable t = (LayeredNetworkAvgTable) flatPH(twoTasks()).getAvgTable();
        assertEquals(1.391660, rowOf(t, "T1", t.getTput()), TOL);
        assertEquals(2.783319, rowOf(t, "T2", t.getTput()), TOL);
        assertEquals(0.695830, rowOf(t, "P1", t.getUtil()), TOL);
        assertEquals(0.927773, rowOf(t, "P2", t.getUtil()), TOL);
    }

    @Test
    void serialChainMatchesMatlab() {
        LayeredNetworkAvgTable t = (LayeredNetworkAvgTable) flatPH(serialChain()).getAvgTable();
        assertEquals(2.126440, rowOf(t, "T1", t.getTput()), TOL);
        assertEquals(2.126440, rowOf(t, "T2", t.getTput()), TOL);
        assertEquals(2.126440, rowOf(t, "T3", t.getTput()), TOL);
    }

    /**
     * The composed law folds every call of an invocation into ONE visit, so the
     * dispatch order of a routed call group has nowhere to be expressed.
     * Squashing does not recover it, which is why the 'flat' alias stays on
     * 'flat.cs' rather than probing this method.
     */
    @Test
    void flatAliasStaysOnTheRoutingEncoding() {
        SolverOptions o = new SolverOptions();
        o.method = "flat";
        o.verbose = VerboseLevel.SILENT;
        SolverLN s = new SolverLN(serialChain(), o);
        s.getAvgTable();
        assertEquals("flat.cs", s.lnmethod);
    }

    /** A replicated element needs a submodel per replica, so squashing refuses it. */
    @Test
    void flatPHRefusesReplication() {
        LayeredNetwork m = serialChain();
        m.getTasks().get(1).setReplication(2);
        assertThrows(RuntimeException.class, () -> flatPH(m).getAvgTable());
    }
}
