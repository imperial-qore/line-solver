/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.auto.SolverCandidate;
import jline.util.matrix.Matrix;

/**
 * SolverMVA's per-method support gate: every (solver, method) row findSolver
 * offers for the MVA family must actually run, and must produce the model's
 * answer rather than a table of zeros.
 *
 * <p>The gate is the same predicate SolverAUTO's ranking consults before it
 * delegates and that listValidMethods projects, so a gate weaker than the
 * analyzer is not a cosmetic defect in a report: it hands a caller a method name
 * that then answers with zeros. The closed-population AMVA family is where that
 * bit -- bs, aql, qsa, sqni, tay, scat, lcp, chow, pamb, pami, pamt, clust,
 * dmlin, ab, schmidt and schmidt-ext each recur on a CLOSED population vector N
 * and are handed (L, N, Z) alone, so on an open model the recursion runs over an
 * empty set of chains and falls out with every metric at zero, silently, with
 * the row still marked runnable and "approx".
 *
 * <p>The numbers asserted are properties of the models, not values read back out
 * of the implementation: an M/M/1 with lambda = 1 and mu = 2 has rho = 1/2 and
 * E[Q] = 1 exactly, and a closed network of N jobs holds N of them somewhere
 * whatever approximation is used.
 */
public class SolverMVAGateTest {

    /** Source -> FCFS Queue -> Sink with lambda = 1, mu = 2, so E[Q] = 1. */
    private static Network mm1() {
        Network model = new Network("mm1");
        Source src = new Source(model, "Source");
        Queue q = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink snk = new Sink(model, "Sink");
        OpenClass c = new OpenClass(model, "C1");
        src.setArrival(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(Network.serialRouting(src, q, snk));
        return model;
    }

    /** Delay -> FCFS Queue, N = 3: a closed product-form network. */
    private static Network repairmen() {
        Network model = new Network("repairmen");
        Delay d = new Delay(model, "Delay");
        Queue q = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C1", 3, d);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(Network.serialRouting(d, q));
        return model;
    }

    /** Delay -> Delay, N = 2: a closed network with NO queueing station. */
    private static Network twoDelays() {
        Network model = new Network("twodelays");
        Delay d1 = new Delay(model, "Delay1");
        Delay d2 = new Delay(model, "Delay2");
        ClosedClass c = new ClosedClass(model, "C1", 2, d1);
        d1.setService(c, new Exp(1.0));
        d2.setService(c, new Exp(2.0));
        model.link(Network.serialRouting(d1, d2));
        return model;
    }

    /** Source -> HOL Queue -> Sink, two open classes at different priorities. */
    private static Network holOpen() {
        Network model = new Network("hol");
        Source src = new Source(model, "Source");
        Queue q = new Queue(model, "Queue", SchedStrategy.HOL);
        Sink snk = new Sink(model, "Sink");
        OpenClass hi = new OpenClass(model, "Hi", 0);
        OpenClass lo = new OpenClass(model, "Lo", 1);
        src.setArrival(hi, new Exp(0.4));
        src.setArrival(lo, new Exp(0.4));
        q.setService(hi, new Exp(2.0));
        q.setService(lo, new Exp(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(hi, hi, Network.serialRouting(src, q, snk));
        P.set(lo, lo, Network.serialRouting(src, q, snk));
        model.link(P);
        return model;
    }

    /** Delay -> FCFS Queue with 3 servers, N = 4. */
    private static Network multiserver() {
        Network model = new Network("ms");
        Delay d = new Delay(model, "Delay");
        Queue q = new Queue(model, "Queue", SchedStrategy.FCFS);
        q.setNumberOfServers(3);
        ClosedClass c = new ClosedClass(model, "C1", 4, d);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(Network.serialRouting(d, q));
        return model;
    }

    /**
     * Delay -> PS Queue -> Delay with the class relabelled on each hop.
     *
     * <p>C2 is reached only by switching, so its own population is 0 while the
     * CHAIN holds 2. It is the shape that exposed the extended Schmidt leak.</p>
     */
    private static Network classSwitching() {
        return switchingModel("cs", SchedStrategy.PS);
    }

    /**
     * The same shape with an FCFS queue: the station the -ext correction is
     * formed at, and the empty class it has no customer of to tag.
     */
    private static Network classSwitchingFcfs() {
        return switchingModel("csfcfs", SchedStrategy.FCFS);
    }

    private static Network switchingModel(String name, SchedStrategy sched) {
        Network model = new Network(name);
        Delay d = new Delay(model, "D");
        Queue q = new Queue(model, "Q", sched);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, d);
        ClosedClass c2 = new ClosedClass(model, "C2", 0, d);
        d.setService(c1, new Exp(1.0));
        d.setService(c2, new Exp(1.0));
        q.setService(c1, new Exp(2.0));
        q.setService(c2, new Exp(3.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c2, d, q, 1.0);
        P.set(c2, c1, q, d, 1.0);
        model.link(P);
        return model;
    }

    /** Source -> Fork -> two FCFS queues -> Join -> Sink. */
    private static Network forkJoin() {
        Network model = new Network("fj");
        Source src = new Source(model, "Source");
        Fork f = new Fork(model, "Fork");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Join j = new Join(model, "Join", f);
        Sink snk = new Sink(model, "Sink");
        OpenClass c = new OpenClass(model, "C1");
        src.setArrival(c, new Exp(0.5));
        q1.setService(c, new Exp(2.0));
        q2.setService(c, new Exp(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c, c, src, f, 1.0);
        P.set(c, c, f, q1, 1.0);
        P.set(c, c, f, q2, 1.0);
        P.set(c, c, q1, j, 1.0);
        P.set(c, c, q2, j, 1.0);
        P.set(c, c, j, snk, 1.0);
        model.link(P);
        return model;
    }

    /** The MVA rows of the model's report, as bare method names. */
    private static List<String> mvaMethods(Network model) {
        List<String> out = new ArrayList<String>();
        for (SolverCandidate row : model.findSolver()) {
            if (!"mva".equals(row.solver) || !row.runnable) {
                continue;
            }
            out.add(row.method.substring(row.method.indexOf('.') + 1));
        }
        return out;
    }

    private static Matrix qlen(Network model, String method) {
        SolverMVA s = new SolverMVA(model, "method", method, "verbose", VerboseLevel.SILENT);
        s.getAvg();
        return s.result.QN;
    }

    private static boolean allZero(Matrix m) {
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                if (Math.abs(m.get(i, j)) > 1e-12) {
                    return false;
                }
            }
        }
        return true;
    }

    private static double maxAbs(Matrix m) {
        double best = 0.0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                best = Math.max(best, Math.abs(m.get(i, j)));
            }
        }
        return best;
    }

    private static double sumAll(Matrix m) {
        double s = 0.0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                s += m.get(i, j);
            }
        }
        return s;
    }

    // ---- the open direction ----

    @Test
    public void everyOfferedMethodReturnsTheMm1QueueLength() {
        // rho = 1/2 gives E[Q] = rho/(1-rho) = 1 at the queue. Every method the
        // report offers is offered as a solution of THIS model, so it has to land
        // on that number to its own accuracy; a table of zeros is the failure this
        // gate exists to make impossible.
        List<String> offered = mvaMethods(mm1());
        assertTrue(offered.size() > 10, "the report offers " + offered.size() + " mva rows");
        for (String name : offered) {
            Matrix q = qlen(mm1(), name);
            assertFalse(allZero(q), "method '" + name + "' returned an all-zero queue-length table");
            assertEquals(1.0, maxAbs(q), 0.8, "method '" + name + "' did not put ~1 job at the queue");
        }
    }

    @Test
    public void theClosedPopulationFamilyIsNotOfferedOnAnOpenModel() {
        List<String> offered = mvaMethods(mm1());
        for (String name : SolverMVA.CLOSED_POPULATION_METHODS) {
            assertFalse(offered.contains(name), "'" + name + "' is closed-population only");
            assertFalse(offered.contains("amva." + name), "'amva." + name + "' is closed-population only");
        }
    }

    @Test
    public void askingForAClosedPopulationMethodByNameRaises() {
        // The gate withholds the row and the analyzer refuses the run. Silence, or
        // zeros, would be worse than either.
        for (String name : SolverMVA.CLOSED_POPULATION_METHODS) {
            final String m = name;
            assertThrows(RuntimeException.class, () -> qlen(mm1(), m),
                    "method '" + name + "' ran on an open model");
        }
    }

    @Test
    public void qnaIsWithheldWhereItsStationUpdateHasNoArm() {
        // Solver_qna decomposes each station as an INF, PS or FCFS centre and has
        // no arm for a priority discipline, so it used to leave the HOL station's
        // row of Q, U, R and T at zero and return the table.
        assertFalse(mvaMethods(holOpen()).contains("qna"));
        assertTrue(mvaMethods(mm1()).contains("qna"));
    }

    @Test
    public void theRobustAnalyzersAreWithheldOnAMulticlassModel() {
        // RQNA and RQT build one uncertainty set per flow from the two moments of
        // a single stream.
        List<String> offered = mvaMethods(holOpen());
        assertFalse(offered.contains("rqna"));
        assertFalse(offered.contains("rqt"));
        assertTrue(mvaMethods(mm1()).contains("rqna"));
    }

    @Test
    public void theSummationMethodIsWithheldOnAPriorityStation() {
        // sum/esum pass every station to the summation kernel as an INF, PS,
        // LCFS-PR, FCFS or SIRO centre and refuse the rest by name.
        List<String> offered = mvaMethods(holOpen());
        assertFalse(offered.contains("sum"));
        assertFalse(offered.contains("esum"));
        assertTrue(mvaMethods(mm1()).contains("sum"));
    }

    // ---- the closed direction: nothing that genuinely applies may be lost ----

    @Test
    public void aClosedProductFormNetworkKeepsTheWholeFamily() {
        List<String> offered = mvaMethods(repairmen());
        for (String name : SolverMVA.CLOSED_POPULATION_METHODS) {
            assertTrue(offered.contains(name),
                    "'" + name + "' is defined for this closed product-form model");
        }
    }

    @Test
    public void theFamilyRunsAndConservesThePopulation() {
        for (String name : SolverMVA.CLOSED_POPULATION_METHODS) {
            assertEquals(3.0, sumAll(qlen(repairmen(), name)), 1e-3,
                    "method '" + name + "' lost the closed population");
        }
    }

    @Test
    public void aPureDelayNetworkKeepsTheFamilyAndSolvesItExactly() {
        // With no queueing station there is no arrival-instant correction to make,
        // so every one of these algorithms coincides with the exact delay solution.
        // Refusing them there would be an over-tightening, and handing them a
        // zero-row demand matrix is what made them raise or report zeros.
        List<String> offered = mvaMethods(twoDelays());
        Matrix exact = qlen(twoDelays(), "default");
        for (String name : SolverMVA.CLOSED_POPULATION_METHODS) {
            if ("sqni".equals(name)) {
                // pfqn_sqni is a closed form for ONE queueing station with a delay,
                // so listValidMethods withholds it on any other shape; that is a
                // shape rule of its own, not the closed-chain rule.
                continue;
            }
            assertTrue(offered.contains(name), "'" + name + "' on a pure-delay network");
            Matrix q = qlen(twoDelays(), name);
            for (int i = 0; i < exact.getNumRows(); i++) {
                for (int j = 0; j < exact.getNumCols(); j++) {
                    // "Exactly" means the same ANSWER, to the fixed point's own
                    // stopping rule. 'default' reaches this model through the
                    // exact recursion, while every AMVA name goes through the
                    // load-dependent sweep, which halts at iter_tol (1e-6, the
                    // reference's default) and leaves a 2.6e-7 residual --
                    // identical across all of them, since with no queueing
                    // station they are the same iteration. The C++ twin
                    // (test_gate_mva.cpp) carries the same bound.
                    assertEquals(exact.get(i, j), q.get(i, j), 1e-6,
                            "method '" + name + "' on a pure-delay network");
                }
            }
        }
    }

    @Test
    public void mvacIsWithheldWhereItHasNoRecursion() {
        // Pfqn_mvac recurs over single-server fixed-rate queueing centres; it
        // refused a multiserver station by name, and is handed a zero-row demand
        // matrix on a delay-only model, while the report went on offering it.
        // MATLAB and native python answer the delay-only case from the delay
        // closed form instead; the JAR and the C++ port refuse it, and the gate
        // now says so rather than offering a row the analyzer throws on.
        assertFalse(mvaMethods(multiserver()).contains("mvac"));
        assertFalse(mvaMethods(twoDelays()).contains("mvac"));
        assertTrue(mvaMethods(repairmen()).contains("mvac"));
    }

    // ---- shapes whose population vector or topology is not what a kernel assumes ----

    @Test
    public void schmidtExtIsWithheldWhereItHasNoCustomerToTag() {
        // Schmidt's EXTENSION corrects an FCFS station from the network with one
        // class-r customer TAGGED, i.e. at population N - 1_r. A class reached only
        // by switching holds no customer of its own, so N_r - 1 is negative and the
        // state lattice prod(N+1) collapses to zero. Plain "schmidt" forms no such
        // sub-problem and must keep running.
        //
        // THIS ARM RECURS ON THE CHAIN POPULATIONS under class switching, which is
        // what Pfqn_schmidt_ext is handed, and schmidtExtReasonForStruct asks about
        // that same vector -- one chain holding 2 here, so the rule is inert and
        // the name stays offered. That is what "asked about the numbers this arm
        // passes" means, and the C++ twin (test_gate_mva.cpp) carries the identical
        // case. What must hold either way is that nothing OFFERED then throws.
        List<String> offered = mvaMethods(classSwitchingFcfs());
        assertTrue(offered.contains("schmidt"));
        for (String name : offered) {
            qlen(classSwitchingFcfs(), name);
        }
    }

    @Test
    public void aClassSwitchingModelOffersNoRowThatRaises() {
        // The reported leak: findSolver offered "schmidt-ext" on this model and
        // running it raised. Every row the report offers must run.
        List<String> offered = mvaMethods(classSwitching());
        assertTrue(offered.contains("schmidt-ext"));
        for (String name : offered) {
            Matrix q = qlen(classSwitching(), name);
            assertFalse(allZero(q), "method '" + name + "' answered with an all-zero table");
        }
    }

    @Test
    public void theRobustAnalyzersAreWithheldOnAForkJoinModel() {
        // A Join is a synchronisation node, not a queue: it carries no service
        // process, so the index-of-dispersion curve RQNA and RQT read off every
        // station does not exist for it. QNA keeps Fork/Join: its station loop has
        // an explicit Join arm.
        List<String> offered = mvaMethods(forkJoin());
        assertFalse(offered.contains("rqna"));
        assertFalse(offered.contains("rqt"));
        assertTrue(offered.contains("qna"));
        assertThrows(RuntimeException.class, () -> qlen(forkJoin(), "rqna"));
        assertThrows(RuntimeException.class, () -> qlen(forkJoin(), "rqt"));
    }

    // ---- the gate and the analyzer must be one predicate, not two copies ----

    /** A fresh instance of shape `i`; a solve may refresh a model's struct. */
    private static Network shape(int i) {
        switch (i) {
            case 0: return mm1();
            case 1: return repairmen();
            case 2: return twoDelays();
            case 3: return holOpen();
            case 4: return classSwitching();
            case 5: return classSwitchingFcfs();
            case 6: return forkJoin();
            default: return multiserver();
        }
    }

    @Test
    public void everyOfferedRowRunsOnEveryShape() {
        for (int i = 0; i < 8; i++) {
            String shapeName = shape(i).getName();
            for (String name : mvaMethods(shape(i))) {
                Matrix q;
                try {
                    q = qlen(shape(i), name);
                } catch (RuntimeException exc) {
                    throw new AssertionError(shapeName + " offers '" + name
                            + "' but the run raised: " + exc.getMessage(), exc);
                }
                assertFalse(allZero(q), shapeName + " offers '" + name
                        + "' and it answered with an all-zero table");
            }
        }
    }
}
