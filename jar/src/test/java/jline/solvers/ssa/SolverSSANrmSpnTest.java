package jline.solvers.ssa;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.Mode;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the NRM stochastic-Petri-net path (Solver_ssa_nrm.nrm_spn) against
 * the exact SPN CTMC. Each test asserts (a) the solver actually ran the NRM
 * method (result.method contains "nrm", never a silent serial fallback) and
 * (b) the simulated marking means / transition throughputs agree with the CTMC.
 *
 * Two net classes are covered:
 * <ul>
 *   <li>a closed net with an inhibitor arc (spn_inhibiting), validated directly
 *       against the SPN CTMC that handles it;</li>
 *   <li>a closed net with an IMMEDIATE transition (vanishing markings). The SPN
 *       CTMC path has a known open gap on immediate transitions, so ground truth
 *       is an EQUIVALENT all-timed net: the immediate transition merely relays
 *       a token from a vanishing place, so the reduced net (with that place and
 *       its immediate transition removed) has identical steady-state marking on
 *       every tangible place. The test asserts the vanishing place holds exactly
 *       zero tokens and the two timed places match the reduced-net CTMC.</li>
 * </ul>
 *
 * The single-seed tolerance is 3% relative, several standard deviations above
 * the ~0.3% noise floor at this sample count, so the fixed-seed assertions are
 * not flaky yet still fail on a systematic firing-accounting defect.
 */
public class SolverSSANrmSpnTest {

    private static final int SAMPLES = 300000;
    private static final int SEED = 23000;
    private static final double RTOL = 0.03;

    // ------------------------------------------------------------------
    // Fixtures
    // ------------------------------------------------------------------

    /** Closed 3-place net, 4 tokens, two firing modes and one inhibitor arc. */
    private static Network inhibiting() {
        Network model = new Network("spn_inhibiting");
        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Place P3 = new Place(model, "P3");
        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        Transition T3 = new Transition(model, "T3");
        ClosedClass jobclass = new ClosedClass(model, "Class1", 4, P1, 0);
        Mode m1 = T1.addMode("Mode1");
        T1.setDistribution(m1, new Exp(2));
        T1.setEnablingConditions(m1, jobclass, P1, 2);
        T1.setFiringOutcome(m1, jobclass, P2, 2);
        Mode m2 = T1.addMode("Mode2");
        T1.setDistribution(m2, new Exp(1));
        T1.setEnablingConditions(m2, jobclass, P1, 1);
        T1.setFiringOutcome(m2, jobclass, P3, 1);
        Mode m3 = T2.addMode("Mode3");
        T2.setDistribution(m3, new Exp(4));
        T2.setEnablingConditions(m3, jobclass, P2, 1);
        T2.setFiringOutcome(m3, jobclass, P1, 1);
        Mode m4 = T3.addMode("Mode4");
        T3.setDistribution(m4, new Exp(1));
        T3.setEnablingConditions(m4, jobclass, P3, 3);
        T3.setInhibitingConditions(m4, jobclass, P2, 1);
        T3.setFiringOutcome(m4, jobclass, P1, 3);
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(jobclass, jobclass, P1, T1, 1.0);
        rm.set(jobclass, jobclass, P2, T2, 1.0);
        rm.set(jobclass, jobclass, P2, T3, 1.0);
        rm.set(jobclass, jobclass, P3, T3, 1.0);
        rm.set(jobclass, jobclass, T1, P2, 1.0);
        rm.set(jobclass, jobclass, T1, P3, 1.0);
        rm.set(jobclass, jobclass, T2, P1, 1.0);
        rm.set(jobclass, jobclass, T3, P1, 1.0);
        model.link(rm);
        P1.setState(Matrix.singleton(jobclass.getPopulation()));
        P2.setState(Matrix.singleton(0));
        P3.setState(Matrix.singleton(0));
        return model;
    }

    /**
     * Closed net P1 --T1(Exp3)--> P2 --T2(immediate)--> P3 --T3(Exp2)--> P1,
     * 2 tokens. P2 is a vanishing place: every token entering it is relayed to
     * P3 in zero time by the immediate transition T2.
     */
    private static Network immediate() {
        Network model = new Network("spn_immediate");
        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Place P3 = new Place(model, "P3");
        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        Transition T3 = new Transition(model, "T3");
        ClosedClass jobclass = new ClosedClass(model, "Class1", 2, P1, 0);
        Mode m1 = T1.addMode("M1");
        T1.setDistribution(m1, new Exp(3));
        T1.setEnablingConditions(m1, jobclass, P1, 1);
        T1.setFiringOutcome(m1, jobclass, P2, 1);
        Mode m2 = T2.addMode("M2");
        T2.setDistribution(m2, new Immediate());
        T2.setTimingStrategy(m2, TimingStrategy.IMMEDIATE);
        T2.setFiringPriorities(m2, 1);
        T2.setFiringWeights(m2, 1.0);
        T2.setEnablingConditions(m2, jobclass, P2, 1);
        T2.setFiringOutcome(m2, jobclass, P3, 1);
        Mode m3 = T3.addMode("M3");
        T3.setDistribution(m3, new Exp(2));
        T3.setEnablingConditions(m3, jobclass, P3, 1);
        T3.setFiringOutcome(m3, jobclass, P1, 1);
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(jobclass, jobclass, P1, T1, 1.0);
        rm.set(jobclass, jobclass, P2, T2, 1.0);
        rm.set(jobclass, jobclass, P3, T3, 1.0);
        rm.set(jobclass, jobclass, T1, P2, 1.0);
        rm.set(jobclass, jobclass, T2, P3, 1.0);
        rm.set(jobclass, jobclass, T3, P1, 1.0);
        model.link(rm);
        P1.setState(Matrix.singleton(jobclass.getPopulation()));
        P2.setState(Matrix.singleton(0));
        P3.setState(Matrix.singleton(0));
        return model;
    }

    /** All-timed net equivalent to {@link #immediate()} with P2/T2 removed. */
    private static Network immediateReduced() {
        Network model = new Network("spn_immediate_reduced");
        Place P1 = new Place(model, "P1");
        Place P3 = new Place(model, "P3");
        Transition T1 = new Transition(model, "T1");
        Transition T3 = new Transition(model, "T3");
        ClosedClass jobclass = new ClosedClass(model, "Class1", 2, P1, 0);
        Mode m1 = T1.addMode("M1");
        T1.setDistribution(m1, new Exp(3));
        T1.setEnablingConditions(m1, jobclass, P1, 1);
        T1.setFiringOutcome(m1, jobclass, P3, 1);
        Mode m3 = T3.addMode("M3");
        T3.setDistribution(m3, new Exp(2));
        T3.setEnablingConditions(m3, jobclass, P3, 1);
        T3.setFiringOutcome(m3, jobclass, P1, 1);
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(jobclass, jobclass, P1, T1, 1.0);
        rm.set(jobclass, jobclass, P3, T3, 1.0);
        rm.set(jobclass, jobclass, T1, P3, 1.0);
        rm.set(jobclass, jobclass, T3, P1, 1.0);
        model.link(rm);
        P1.setState(Matrix.singleton(jobclass.getPopulation()));
        P3.setState(Matrix.singleton(0));
        return model;
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    private static SolverSSA nrm(Network model) {
        SolverSSA solver = new SolverSSA(model);
        solver.options.method = "nrm";
        solver.options.samples = SAMPLES;
        solver.options.seed = SEED;
        return solver;
    }

    private static void assertRanNrm(SolverSSA solver) {
        assertTrue(solver.result != null && solver.result.method != null
                        && solver.result.method.contains("nrm"),
                "SPN model must run the NRM method, got: "
                        + (solver.result == null ? "<no result>" : solver.result.method));
    }

    private static void assertClose(double expected, double got, String what) {
        double denom = Math.max(Math.abs(expected), 1e-9);
        double err = Math.abs(got - expected) / denom;
        assertTrue(err < RTOL, what + ": NRM " + got + " deviates from CTMC "
                + expected + " by " + (100.0 * err) + "%");
    }

    // ------------------------------------------------------------------
    // Tests
    // ------------------------------------------------------------------

    @Test
    public void testInhibitingVsCtmc() {
        NetworkAvgTable ctmc = new SolverCTMC(inhibiting(), "cutoff", 10, "seed", 1).getAvgTable();
        SolverSSA solver = nrm(inhibiting());
        NetworkAvgTable ssa = solver.getAvgTable();
        assertRanNrm(solver);
        List<Double> qCtmc = ctmc.getQLen();
        List<Double> tCtmc = ctmc.getTput();
        List<Double> qSsa = ssa.getQLen();
        List<Double> tSsa = ssa.getTput();
        for (int i = 0; i < qCtmc.size(); i++) {
            assertClose(qCtmc.get(i), qSsa.get(i), "QLen[" + i + "]");
            assertClose(tCtmc.get(i), tSsa.get(i), "Tput[" + i + "]");
        }
    }

    @Test
    public void testImmediateVanishingMarking() {
        // Ground truth: the reduced all-timed net (P2/T2 collapsed).
        NetworkAvgTable red = new SolverCTMC(immediateReduced(), "cutoff", 10, "seed", 1).getAvgTable();
        double qP1exact = red.getQLen().get(0);   // reduced station 0 = P1
        double qP3exact = red.getQLen().get(1);   // reduced station 1 = P3
        double tExact = red.getTput().get(0);

        SolverSSA solver = nrm(immediate());
        NetworkAvgTable ssa = solver.getAvgTable();
        assertRanNrm(solver);
        List<Double> q = ssa.getQLen();   // stations: P1(0), P2(1), P3(2)
        List<Double> t = ssa.getTput();

        // The vanishing place holds exactly zero tokens (never observed with a
        // token because the immediate transition drains it in zero time).
        assertEquals(0.0, q.get(1), 1e-12, "vanishing place P2 mean tokens");
        assertClose(qP1exact, q.get(0), "QLen P1");
        assertClose(qP3exact, q.get(2), "QLen P3");
        // T1 and T3 throughputs equal the reduced-net cycle throughput.
        assertClose(tExact, t.get(0), "Tput T1");
        assertClose(tExact, t.get(2), "Tput T3");
    }
}
