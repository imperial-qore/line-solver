package jline.solvers.ln;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.*;
import jline.lang.processes.Exp;
import jline.solvers.LayeredNetworkAvgTable;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

/**
 * Admission constraints on an LN layer station.
 *
 * <p>A Task or Host may declare A*n &lt;= b on the station that represents it in
 * its layer (addConstraint / setConstraint), which is how a passive resource
 * (semaphore, buffer pool, admission method name) is expressed in an LQN.
 * SolverLN.buildLayersRecursive emits the rows as a Region on the layer's server
 * station, expanding an entry column to the CALL classes targeting that entry,
 * and DefaultSolverFactory escalates the layer solver to one whose feature set
 * covers Region.</p>
 *
 * <p>There is no external oracle: JMT's FCR XML cannot express a general
 * A*n &lt;= b and lqns has no admission-constraint concept, so these tests assert
 * structural and conservation properties rather than golden numbers. Twin of
 * line-test.git test_fcr_lincon_ln.m and python
 * tests/test_ln_admission_constraint.py. The population is deliberately small:
 * a constrained layer is solved by SolverCTMC, whose state-space generation is
 * quadratic in the JAR.</p>
 */
class SolverLNAdmissionConstraintTest {

    private static final double TOL_EQ = 1e-6;
    private static final double TOL_BAL = 1e-3;

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /** Two-tier LQN: T1 (reference, N=2) calls entries E2 and E3 of T2 in sequence. */
    private static LayeredNetwork mkLqn() {
        LayeredNetwork model = new LayeredNetwork("fcrlqn");
        Processor p1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor p2 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Task t1 = new Task(model, "T1", 2, SchedStrategy.REF).on(p1);
        Task t2 = new Task(model, "T2", 2, SchedStrategy.FCFS).on(p2);
        Entry e1 = new Entry(model, "E1").on(t1);
        Entry e2 = new Entry(model, "E2").on(t2);
        Entry e3 = new Entry(model, "E3").on(t2);
        Activity a1 = new Activity(model, "A1", new Exp(1.0)).on(t1).boundTo(e1).synchCall(e2, 1.0);
        Activity a1b = new Activity(model, "A1b", new Exp(1.0)).on(t1).synchCall(e3, 1.0);
        t1.addPrecedence(ActivityPrecedence.Serial(a1, a1b));
        new Activity(model, "A2", new Exp(2.0)).on(t2).boundTo(e2).repliesTo(e2);
        new Activity(model, "A3", new Exp(2.0)).on(t2).boundTo(e3).repliesTo(e3);
        return model;
    }

    private static Task taskOf(LayeredNetwork model, String name) {
        for (Task t : model.getTasks().values()) {
            if (t.getName().equals(name)) {
                return t;
            }
        }
        throw new IllegalArgumentException("no task " + name);
    }

    private static Entry entryOf(LayeredNetwork model, String name) {
        for (Task t : model.getTasks().values()) {
            for (Entry e : t.getEntries()) {
                if (e.getName().equals(name)) {
                    return e;
                }
            }
        }
        throw new IllegalArgumentException("no entry " + name);
    }

    /** Model whose T2 admits at most cap calls across its two entries. */
    private static LayeredNetwork withCap(double cap) {
        LayeredNetwork model = mkLqn();
        taskOf(model, "T2").addConstraint(
                Arrays.asList(entryOf(model, "E2"), entryOf(model, "E3")),
                new double[]{1, 1}, cap);
        return model;
    }

    private static Matrix[] linconOfT2(LayeredNetwork model) {
        LayeredNetworkStruct lqn = model.getStruct();
        // T2 is the second task and the element space is 0-based: tshift+0 is T1.
        return lqn.lincon.get(lqn.tshift + 1);
    }

    /** Value of column col on the row whose node name is nodeName. */
    private static double rowOf(LayeredNetworkAvgTable t, String nodeName, String col) {
        List<String> names = t.getNodeNames();
        int i = names.indexOf(nodeName);
        assertTrue(i >= 0, "row " + nodeName + " not in " + names);
        if ("QLen".equals(col)) {
            return t.getQLen().get(i);
        } else if ("Tput".equals(col)) {
            return t.getTput().get(i);
        } else if ("RespT".equals(col)) {
            return t.getRespT().get(i);
        }
        throw new IllegalArgumentException("no column " + col);
    }

    private static LayeredNetworkAvgTable avgTable(LayeredNetwork model) {
        return (LayeredNetworkAvgTable) new SolverLN(model).getAvgTable();
    }

    @Test
    void declarationFormsAgree() {
        LayeredNetwork m1 = mkLqn();
        Matrix a = new Matrix(2, 2);
        a.set(0, 0, 1);
        a.set(0, 1, 1);
        a.set(1, 0, 0);
        a.set(1, 1, 1);
        Matrix b = new Matrix(2, 1);
        b.set(0, 0, 2);
        b.set(1, 0, 1);
        taskOf(m1, "T2").setConstraint(a, b);
        Matrix[] ab1 = linconOfT2(m1);

        LayeredNetwork m2 = mkLqn();
        Task t2 = taskOf(m2, "T2");
        t2.addConstraint(Arrays.asList(entryOf(m2, "E2"), entryOf(m2, "E3")), new double[]{1, 1}, 2);
        t2.addConstraint(Arrays.asList(entryOf(m2, "E3")), new double[]{1}, 1);
        Matrix[] ab2 = linconOfT2(m2);

        LayeredNetwork m3 = mkLqn();
        Task t3 = taskOf(m3, "T2");
        // null coefficients default to all ones
        t3.addConstraintByName(Arrays.asList("E2", "E3"), null, 2);
        t3.addConstraintByName(Arrays.asList("E3"), new double[]{1}, 1);
        Matrix[] ab3 = linconOfT2(m3);

        assertEquals(2, ab1[0].getNumRows());
        assertEquals(2, ab1[0].getNumCols());
        for (int r = 0; r < 2; r++) {
            assertEquals(ab1[1].get(r, 0), ab2[1].get(r, 0), TOL_EQ);
            assertEquals(ab1[1].get(r, 0), ab3[1].get(r, 0), TOL_EQ);
            for (int c = 0; c < 2; c++) {
                assertEquals(ab1[0].get(r, c), ab2[0].get(r, c), TOL_EQ);
                assertEquals(ab1[0].get(r, c), ab3[0].get(r, c), TOL_EQ);
            }
        }
    }

    @Test
    void foreignOperandRejected() {
        // E1 belongs to T1, so it cannot appear in a constraint on T2
        LayeredNetwork model = mkLqn();
        taskOf(model, "T2").addConstraintByName(Arrays.asList("E1"), new double[]{1}, 1);
        IllegalArgumentException ex = assertThrows(IllegalArgumentException.class, model::getStruct);
        assertTrue(ex.getMessage().contains("not one of the entries"), ex.getMessage());
    }

    @Test
    void duplicateOperandRejected() {
        LayeredNetwork model = mkLqn();
        IllegalArgumentException ex = assertThrows(IllegalArgumentException.class,
                () -> taskOf(model, "T2").addConstraintByName(
                        Arrays.asList("E2", "E2"), new double[]{1, 1}, 2));
        assertTrue(ex.getMessage().contains("same operand more than once"), ex.getMessage());
    }

    @Test
    void wrongColumnCountRejected() {
        // the positional form checks the column count, which is its only guard
        LayeredNetwork model = mkLqn();
        Matrix a = new Matrix(1, 3);
        a.set(0, 0, 1);
        a.set(0, 1, 1);
        a.set(0, 2, 1);
        Matrix b = new Matrix(1, 1);
        b.set(0, 0, 2);
        taskOf(model, "T2").setConstraint(a, b);
        IllegalArgumentException ex = assertThrows(IllegalArgumentException.class, model::getStruct);
        assertTrue(ex.getMessage().contains("columns but there are"), ex.getMessage());
    }

    /**
     * Exactly one layer carries the Region, and the solver records it so the
     * chain recovery runs there. The escalation off MVA is asserted indirectly
     * by constraintIsEnforcedAndSatisfied: MVA has no Region support, so an
     * unescalated layer returns the unconstrained numbers.
     */
    @Test
    void oneLayerCarriesTheRegion() {
        SolverLN solver = new SolverLN(withCap(1));
        solver.getAvgTable();
        List<Integer> constrained = new ArrayList<Integer>();
        List<jline.lang.Network> ensemble = solver.getEnsemble();
        for (int e = 0; e < ensemble.size(); e++) {
            if (ensemble.get(e).getRegions() != null && !ensemble.get(e).getRegions().isEmpty()) {
                constrained.add(e);
            }
        }
        assertEquals(1, constrained.size());
        assertNotNull(solver.layerHasRegion);
        assertTrue(solver.layerHasRegion[constrained.get(0)]);
        assertNotNull(solver.layerChains[constrained.get(0)], "chain matrix not cached");
    }

    @Test
    void constraintIsEnforcedAndSatisfied() {
        LayeredNetworkAvgTable free = avgTable(mkLqn());
        LayeredNetworkAvgTable con = avgTable(withCap(1));
        double qFree = rowOf(free, "E2", "QLen");
        double qE2 = rowOf(con, "E2", "QLen");
        double qE3 = rowOf(con, "E3", "QLen");
        assertTrue(Math.abs(qE2 - qFree) > TOL_EQ, "constraint is inert");
        assertTrue(qE2 + qE3 <= 1 + TOL_BAL, "n(E2)+n(E3) = " + (qE2 + qE3));
    }

    /**
     * A job blocked at the region is counted at no station, so without the
     * chain-level Little recovery in SolverLN.regionWait the fixed point has the
     * calling activity and the called entry disagreeing on throughput.
     */
    @Test
    void flowBalanceAcrossLayerBoundary() {
        LayeredNetworkAvgTable con = avgTable(withCap(1));
        double xA1 = rowOf(con, "A1", "Tput");
        double xE2 = rowOf(con, "E2", "Tput");
        assertEquals(xA1, xE2, TOL_BAL * Math.abs(xA1));
    }

    /**
     * Relaxing the cap cannot reduce throughput. Compared only within the
     * constrained family: a constrained layer is solved by SolverCTMC and an
     * unconstrained one by SolverMVA, so the constrained-vs-free gap carries a
     * solver discrepancy on top of the constraint effect.
     */
    @Test
    void throughputMonotoneInCap() {
        double x1 = rowOf(avgTable(withCap(1)), "E2", "Tput");
        double x2 = rowOf(avgTable(withCap(2)), "E2", "Tput");
        double x3 = rowOf(avgTable(withCap(3)), "E2", "Tput");
        assertTrue(x1 <= x2 + TOL_EQ, x1 + " > " + x2);
        assertTrue(x2 <= x3 + TOL_EQ, x2 + " > " + x3);
    }

    /**
     * A cap that cannot bind (3 method names for 2 jobs) reproduces the unconstrained
     * throughput up to the CTMC-vs-MVA layer-solver discrepancy.
     */
    @Test
    void slackCapMatchesUnconstrained() {
        double xSlack = rowOf(avgTable(withCap(3)), "E2", "Tput");
        double xFree = rowOf(avgTable(mkLqn()), "E2", "Tput");
        assertEquals(xFree, xSlack, 0.01 * Math.abs(xFree));
    }

    @Test
    void jsonRoundTripPreservesConstraint() throws Exception {
        // the wire form is named, so it is order-independent by construction
        LayeredNetwork model = mkLqn();
        Matrix a = new Matrix(2, 2);
        a.set(0, 0, 1);
        a.set(0, 1, 1);
        a.set(1, 0, 0);
        a.set(1, 1, 1);
        Matrix b = new Matrix(2, 1);
        b.set(0, 0, 2);
        b.set(1, 0, 1);
        taskOf(model, "T2").setConstraint(a, b);

        java.io.File f = java.io.File.createTempFile("ln_lincon_rt", ".json");
        f.deleteOnExit();
        jline.io.LineModelIO.save(model, f.getAbsolutePath());
        LayeredNetwork back = (LayeredNetwork) jline.io.LineModelIO.load(f.getAbsolutePath());

        Matrix[] before = linconOfT2(model);
        Matrix[] after = linconOfT2(back);
        assertNotNull(after, "admission constraint lost on the wire");
        assertEquals(before[0].getNumRows(), after[0].getNumRows());
        for (int r = 0; r < before[0].getNumRows(); r++) {
            assertEquals(before[1].get(r, 0), after[1].get(r, 0), TOL_EQ);
            for (int c = 0; c < before[0].getNumCols(); c++) {
                assertEquals(before[0].get(r, c), after[0].get(r, c), TOL_EQ);
            }
        }
    }
}
