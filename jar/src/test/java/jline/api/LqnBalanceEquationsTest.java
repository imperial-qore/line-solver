/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.lqn.LqnBalanceEquations;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkStruct;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.processes.Exp;
import jline.solvers.ln.SolverLN;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the LQN conservation-law enumerator ({@link LqnBalanceEquations}).
 *
 * <p>Two things are checked, and they are different in kind. STRUCTURE, with no
 * solution supplied: the right relations are emitted with the right term sets,
 * branches and constants, which depends on the model only. CONSISTENCY, on a
 * converged SolverLN: every relation the enumerator emits is a law the solution must
 * satisfy, so each residual must vanish to the solver's own tolerance. The second is
 * the real content -- it turns the enumerator into a conservation check on SolverLN
 * itself.</p>
 *
 * <p>The model is the three-tier LQN of {@code lqn_basic}: T1 (reference, 50 threads,
 * think 2) calls T2 (50 threads) once, which calls T3 (25 threads) five times, on two
 * processors of multiplicity 2 and 3.</p>
 */
public class LqnBalanceEquationsTest {

    private static final double TOL = 1e-6;

    private static LayeredNetwork build() {
        LayeredNetwork model = new LayeredNetwork("cs");
        Processor P1 = new Processor(model, "P1", 2, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 3, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(1.0 / 2));
        Task T2 = new Task(model, "T2", 50, SchedStrategy.FCFS).on(P1);
        T2.setThinkTime(new Exp(1.0 / 3));
        Task T3 = new Task(model, "T3", 25, SchedStrategy.FCFS).on(P2);
        T3.setThinkTime(new Exp(1.0 / 4));
        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T3);
        new Activity(model, "AS1", new Exp(10)).on(T1).boundTo(E1).synchCall(E2, 1);
        new Activity(model, "AS2", new Exp(20)).on(T2).boundTo(E2).synchCall(E3, 5).repliesTo(E2);
        new Activity(model, "AS3", new Exp(50)).on(T3).boundTo(E3).repliesTo(E3);
        return model;
    }

    private static LqnBalanceEquations.Relation by(LqnBalanceEquations.Result out,
                                                   String kind, String name) {
        for (LqnBalanceEquations.Relation r : out.eqs) {
            if (!kind.equals(r.kind)) {
                continue;
            }
            // the JAR prefixes a hashname by kind (P:/R:/T:/E:/A:), so match on the suffix
            if (r.targetname.equals(name) || r.targetname.endsWith(name)) {
                return r;
            }
        }
        throw new AssertionError("no " + kind + " relation on " + name);
    }

    private static int count(LqnBalanceEquations.Result out, String kind) {
        int n = 0;
        for (LqnBalanceEquations.Relation r : out.eqs) {
            if (kind.equals(r.kind)) {
                n++;
            }
        }
        return n;
    }

    @Test
    public void testFamiliesAndCounts() {
        LqnBalanceEquations.Result out = LqnBalanceEquations.compute(build().getStruct());
        assertEquals(3, count(out, "little"));      // one per task
        assertEquals(2, count(out, "callflow"));    // one per call
        assertEquals(2, count(out, "entryflow"));   // E2, E3; E1 rides the reference cycle
        assertEquals(3, count(out, "actflow"));     // one per activity
        assertEquals(2, count(out, "hostutil"));    // one per processor
    }

    @Test
    public void testLittleBranchesAndConstants() {
        LqnBalanceEquations.Result out = LqnBalanceEquations.compute(build().getStruct());
        LqnBalanceEquations.Relation r1 = by(out, "little", "T1");
        assertEquals("ref", r1.branch);
        assertEquals(50.0, r1.rhsconst, 0.0);
        assertEquals(1, r1.terms.length);
        assertTrue(r1.termisentry[0]);              // its own entry drives its cycle
        LqnBalanceEquations.Relation r2 = by(out, "little", "T2");
        assertEquals("queueing", r2.branch);
        assertEquals(50.0, r2.rhsconst, 0.0);
        assertFalse(r2.termisentry[0]);             // one call class
        assertTrue(r2.scaled);                      // U is normalized to [0,1] here
        assertEquals(25.0, by(out, "little", "T3").rhsconst, 0.0);
    }

    @Test
    public void testHostUtilCoefficientsAreTheHostDemands() {
        LqnBalanceEquations.Result out = LqnBalanceEquations.compute(build().getStruct());
        LqnBalanceEquations.Relation p1 = by(out, "hostutil", "P1");
        assertEquals(2.0, p1.rhsconst, 0.0);        // the declared multiplicity of P1
        assertEquals(2, p1.terms.length);
        double sum = p1.coeff[0] + p1.coeff[1];
        assertEquals(0.15, sum, 1e-12);             // 1/10 on AS1 plus 1/20 on AS2
        LqnBalanceEquations.Relation p2 = by(out, "hostutil", "P2");
        assertEquals(3.0, p2.rhsconst, 0.0);
        assertEquals(0.02, p2.coeff[0], 1e-12);
    }

    @Test
    public void testCallflowCarriesTheMeanCallCount() {
        LqnBalanceEquations.Result out = LqnBalanceEquations.compute(build().getStruct());
        // the JAR prefixes BOTH sides of a call hashname: A:AS1=>E:E2
        assertEquals(1.0, by(out, "callflow", "A:AS1=>E:E2").coeff[0], 1e-12);
        assertEquals(5.0, by(out, "callflow", "A:AS2=>E:E3").coeff[0], 1e-12);
    }

    @Test
    public void testVisitsAreOneOnASingleActivityEntry() {
        LqnBalanceEquations.Result out = LqnBalanceEquations.compute(build().getStruct());
        for (LqnBalanceEquations.Relation r : out.eqs) {
            if ("actflow".equals(r.kind)) {
                assertEquals(1.0, r.coeff[0], 1e-12);
            }
        }
    }

    @Test
    public void testIncidenceMatricesCarryTheCallClassesOnly() {
        LayeredNetworkStruct lqn = build().getStruct();
        LqnBalanceEquations.Result out = LqnBalanceEquations.compute(lqn);
        assertEquals(1.0, out.A_little.get(2, 1), 0.0);   // T2 is a caller class of call 1
        assertEquals(1.0, out.A_little.get(3, 2), 0.0);   // T3 of call 2
        double t1row = 0;
        for (int c = 0; c <= lqn.ncalls; c++) {
            t1row += out.A_little.get(1, c);
        }
        assertEquals(0.0, t1row, 0.0);                    // T1 is entry-driven
    }

    @Test
    public void testSymbolicReportIsPrintableWithoutASolution() {
        LqnBalanceEquations.Result out = LqnBalanceEquations.compute(build().getStruct());
        String txt = out.toString();
        assertTrue(txt.contains("thread-pool Little"));
        assertTrue(txt.contains("host utilization law"));
        assertTrue(Double.isNaN(out.maxresidual));
    }

    @Test
    public void testEveryRelationVanishesOnAConvergedSolve() {
        LayeredNetwork model = build();
        LayeredNetworkStruct lqn = model.getStruct();
        SolverLN solver = new SolverLN(model, m -> new SolverMVA(m));
        Matrix un = reportedUtil(solver, lqn);      // this is what runs the fixed point
        LqnBalanceEquations.Solution sol = LqnBalanceEquations.Solution.fromLayeredIterates(
                solver.tput, solver.util, solver.thinkt, solver.servt, solver.residt);
        sol.un = un;
        LqnBalanceEquations.Result out = LqnBalanceEquations.compute(lqn, sol);
        StringBuilder bad = new StringBuilder();
        for (LqnBalanceEquations.Relation r : out.eqs) {
            if (r.degenerate || Double.isNaN(r.residual)) {
                continue;
            }
            if (Math.abs(r.residual) > TOL) {
                bad.append(String.format("%s %s: residual %.3e (lhs %.6g, rhs %.6g)%n",
                        r.kind, r.targetname, r.residual, r.lhs, r.rhs));
            }
        }
        assertTrue(bad.length() == 0, "conservation violated:\n" + bad);
        assertTrue(out.maxresidual <= TOL, "max residual " + out.maxresidual);
    }

    @Test
    public void testPerClassUtilizationReproducesTheSolverIterate() {
        // The per-call decomposition of the busy threads must add up to the utilization
        // SolverLN carries for that task, which is what makes the emitted per-class U a
        // refinement of the solver's aggregate rather than a new quantity. The two reach
        // it by different routes through the same iterate, so they agree to the fixed
        // point's tolerance and not to machine precision.
        LayeredNetwork model = build();
        LayeredNetworkStruct lqn = model.getStruct();
        SolverLN solver = new SolverLN(model, m -> new SolverMVA(m));
        solver.getEnsembleAvg();
        LqnBalanceEquations.Solution sol = LqnBalanceEquations.Solution.fromLayeredIterates(
                solver.tput, solver.util, solver.thinkt, solver.servt, solver.residt);
        LqnBalanceEquations.Result out = LqnBalanceEquations.compute(lqn, sol);
        for (String name : new String[]{"T2", "T3"}) {
            LqnBalanceEquations.Relation r = by(out, "little", name);
            double sum = 0;
            for (double u : r.perclassutil) {
                sum += u;
            }
            double want = solver.util.get(r.target - 1);   // 0-based over elements
            assertEquals(want, sum, Math.max(1e-9, 1e-6 * Math.abs(want)));
        }
    }

    @Test
    public void testHostUtilIsNotInstantiatedWithoutTheReportedUtilization() {
        // The host law needs the REPORTED utilization, which the `util` iterate does not
        // carry. Without it the record must stay symbolic rather than report a residual
        // against a zero it never meant.
        LayeredNetwork model = build();
        LayeredNetworkStruct lqn = model.getStruct();
        SolverLN solver = new SolverLN(model, m -> new SolverMVA(m));
        solver.getEnsembleAvg();
        LqnBalanceEquations.Solution sol = LqnBalanceEquations.Solution.fromLayeredIterates(
                solver.tput, solver.util, solver.thinkt, solver.servt, solver.residt);
        LqnBalanceEquations.Result out = LqnBalanceEquations.compute(lqn, sol);
        for (LqnBalanceEquations.Relation r : out.eqs) {
            if ("hostutil".equals(r.kind)) {
                assertTrue(Double.isNaN(r.residual));
            } else if (!r.degenerate) {
                assertFalse(Double.isNaN(r.residual), r.kind + " " + r.targetname);
            }
        }
    }

    /**
     * The reported utilization per element, read off the ensemble average table. The
     * `util` iterate is a different quantity -- a task's utilization as a server in its
     * own task layer -- and is left at zero on a host, so the host law needs this one.
     */
    private static Matrix reportedUtil(SolverLN solver, LayeredNetworkStruct lqn) {
        Matrix un = new Matrix(1, lqn.nidx + 1);
        for (int i = 0; i <= lqn.nidx; i++) {
            un.set(0, i, Double.NaN);
        }
        java.util.List<Double> util = solver.getEnsembleAvg().get(1);   // the Util column
        for (int i = 1; i <= lqn.nidx && i - 1 < util.size(); i++) {
            un.set(0, i, util.get(i - 1));
        }
        return un;
    }
}
