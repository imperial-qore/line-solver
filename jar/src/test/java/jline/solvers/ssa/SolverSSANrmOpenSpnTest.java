package jline.solvers.ssa;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.Mode;
import jline.lang.nodes.Place;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Transition;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.wrappers.jmt.SolverJMT;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the NRM stochastic-Petri-net path (Solver_ssa_nrm.nrm_spn) on OPEN
 * nets: a Source feeds a Place whose tokens drain through a Transition to a Sink.
 * Before the Source-arrival reaction was added the fed Place stayed empty and the
 * run threw "Deadlock: no transition is enabled".
 *
 * Each test asserts (a) the solver actually ran the NRM method (result.method
 * contains "nrm", never a silent serial fallback) and (b) the simulated marking
 * means / throughputs match the analytic M/M/1 result. A Source Exp(lambda)
 * feeding a single-server Transition Exp(mu) is an M/M/1 queue at the Place:
 * mean tokens = rho/(1-rho), throughput = lambda (rho = lambda/mu &lt; 1). The
 * canonical net is additionally cross-checked against SolverJMT.
 *
 * Single-seed tolerance is 4% relative at 3e5 samples, several sigma above the
 * ~0.3% noise floor, so the fixed-seed assertions are not flaky yet still fail on
 * a systematic arrival- or firing-accounting defect.
 */
public class SolverSSANrmOpenSpnTest {

    private static final int SAMPLES = 300000;
    private static final int SEED = 23000;
    private static final double RTOL = 0.04;

    /** Open M/M/1-as-SPN: Source Exp(lambda) -> P1 -> T1 Exp(mu) -> Sink. */
    private static Network mm1(double lambda, double mu) {
        Network model = new Network("mm1spn");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");
        Place P1 = new Place(model, "P1");
        Transition T1 = new Transition(model, "T1");
        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, new Exp(lambda));
        Mode m1 = T1.addMode("Mode1");
        T1.setDistribution(m1, new Exp(mu));
        T1.setEnablingConditions(m1, jobclass, P1, 1);
        T1.setFiringOutcome(m1, jobclass, sink, 1);
        model.link(Network.serialRouting(source, P1, T1, sink));
        return model;
    }

    /**
     * Open two-place tandem: Source Exp(lambda) -> P1 -> T1 Exp(mu1) -> P2 ->
     * T2 Exp(mu2) -> Sink. Both places are independent M/M/1 queues fed at
     * lambda: mean tokens P_i = rho_i/(1-rho_i), rho_i = lambda/mu_i.
     */
    private static Network tandem(double lambda, double mu1, double mu2) {
        Network model = new Network("tandemspn");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");
        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, new Exp(lambda));
        Mode m1 = T1.addMode("Mode1");
        T1.setDistribution(m1, new Exp(mu1));
        T1.setEnablingConditions(m1, jobclass, P1, 1);
        T1.setFiringOutcome(m1, jobclass, P2, 1);
        Mode m2 = T2.addMode("Mode1");
        T2.setDistribution(m2, new Exp(mu2));
        T2.setEnablingConditions(m2, jobclass, P2, 1);
        T2.setFiringOutcome(m2, jobclass, sink, 1);
        model.link(Network.serialRouting(source, P1, T1, P2, T2, sink));
        return model;
    }

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
                "open SPN must run the NRM method, got: "
                        + (solver.result == null ? "<no result>" : solver.result.method));
    }

    private static void assertClose(double expected, double got, double rtol, String what) {
        double denom = Math.max(Math.abs(expected), 1e-9);
        double err = Math.abs(got - expected) / denom;
        assertTrue(err < rtol, what + ": NRM " + got + " deviates from " + expected
                + " by " + (100.0 * err) + "%");
    }

    @Test
    public void testMM1OpenSpnRho05() {
        double lambda = 0.5, mu = 1.0, rho = lambda / mu;
        double qExact = rho / (1.0 - rho);   // = 1.0

        SolverSSA solver = nrm(mm1(lambda, mu));
        NetworkAvgTable ssa = solver.getAvgTable();
        assertRanNrm(solver);
        List<Double> q = ssa.getQLen();      // stations: Source(0), P1(1)
        List<Double> t = ssa.getTput();

        assertClose(qExact, q.get(1), RTOL, "P1 mean tokens");
        assertClose(lambda, t.get(0), RTOL, "Source throughput");
        assertClose(lambda, t.get(1), RTOL, "P1 (T1) throughput");

        // Cross-check the mean tokens against JMT's simulation of the same net.
        NetworkAvgTable jmt = new SolverJMT(mm1(lambda, mu), "seed", SEED,
                "samples", SAMPLES).getAvgTable();
        assertClose(jmt.getQLen().get(1), q.get(1), RTOL, "P1 mean tokens vs JMT");
    }

    @Test
    public void testMM1OpenSpnRho08() {
        double lambda = 0.8, mu = 1.0, rho = lambda / mu;
        double qExact = rho / (1.0 - rho);   // = 4.0

        SolverSSA solver = nrm(mm1(lambda, mu));
        NetworkAvgTable ssa = solver.getAvgTable();
        assertRanNrm(solver);
        assertClose(qExact, ssa.getQLen().get(1), RTOL, "P1 mean tokens");
        assertClose(lambda, ssa.getTput().get(0), RTOL, "Source throughput");
        assertClose(lambda, ssa.getTput().get(1), RTOL, "P1 (T1) throughput");
    }

    @Test
    public void testOpenTandemTwoPlaces() {
        double lambda = 0.5, mu1 = 1.0, mu2 = 2.0;
        double q1 = (lambda / mu1) / (1.0 - lambda / mu1);   // = 1.0
        double q2 = (lambda / mu2) / (1.0 - lambda / mu2);   // = 1/3

        SolverSSA solver = nrm(tandem(lambda, mu1, mu2));
        NetworkAvgTable ssa = solver.getAvgTable();
        assertRanNrm(solver);
        List<Double> q = ssa.getQLen();      // stations: Source(0), P1(1), P2(2)
        List<Double> t = ssa.getTput();

        assertClose(q1, q.get(1), RTOL, "P1 mean tokens");
        assertClose(q2, q.get(2), RTOL, "P2 mean tokens");
        assertClose(lambda, t.get(1), RTOL, "P1 (T1) throughput");
        assertClose(lambda, t.get(2), RTOL, "P2 (T2) throughput");
    }
}
