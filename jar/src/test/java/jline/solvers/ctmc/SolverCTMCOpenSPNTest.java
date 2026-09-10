package jline.solvers.ctmc;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.processes.Exp;
import jline.lang.nodes.Place;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Transition;
import jline.lang.Mode;
import jline.lang.OpenClass;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Bounded OPEN stochastic Petri net solved by SolverCTMC.
 *
 * Source Exp(0.5) -> Place P1 -> Transition T1 Exp(1.0) -> Sink is an M/M/1/K
 * queue (rho = lambda/mu = 0.5) once the P1 marking is bounded, by a finite
 * per-class Place capacity (setClassCapacity) or by the solver cutoff. The
 * enabling arc of T1 consumes one token from P1 on firing (routed to the Sink);
 * NO firing outcome is set on P1 (that would double-consume). The exact
 * stationary mean tokens at P1 and the throughput are
 *   E[P1] = sum_{n=0..K} n rho^n / sum_{n=0..K} rho^n,
 *   Tput  = lambda (1 - rho^K / sum_{n=0..K} rho^n).
 *
 * Before the fix, three defects hid here: (1) the Source->Place arrival threw
 * "Index -1 out of bounds" (the FIRE enabling-degree loop indexed outglspace
 * with the Sink's nodeToStateful = -1); (2) once past that, the ARV event wrote
 * the arriving token into a server-phase slot the ordinary place lacks, so the
 * marking never grew and the generator became a pure-death chain; (3) a firing
 * consumed enDegree*weight tokens instead of the arc weight. Closed SPN stays
 * exact (see {@link SpnMarkingInvariantTputTest}).
 */
public class SolverCTMCOpenSPNTest {

    private static final double TOL = 1e-6;
    private static final double LAMBDA = 0.5;
    private static final double MU = 1.0;
    private static final double RHO = LAMBDA / MU;

    private static double mm1kMean(int K) {
        double den = 0.0, num = 0.0;
        for (int n = 0; n <= K; n++) { den += Math.pow(RHO, n); }
        for (int n = 0; n <= K; n++) { num += n * Math.pow(RHO, n); }
        return num / den;
    }

    private static double mm1kTput(int K) {
        double den = 0.0;
        for (int n = 0; n <= K; n++) { den += Math.pow(RHO, n); }
        double pBlock = Math.pow(RHO, K) / den;
        return LAMBDA * (1.0 - pBlock);
    }

    /** Build Source Exp(0.5) -> P1 -> T1 Exp(1.0) -> Sink; cap<0 means no capacity. */
    private static Network buildOpenSpn(int cap) {
        Network model = new Network("ospn");
        Source source = new Source(model, "source");
        Sink sink = new Sink(model, "sink");
        Place p1 = new Place(model, "P1");
        Transition t1 = new Transition(model, "T1");
        OpenClass jc = new OpenClass(model, "jobs");

        source.setArrival(jc, Exp.fitMean(1.0 / LAMBDA));
        if (cap >= 0) {
            p1.setClassCapacity(jc, cap);
        }
        Mode fire = t1.addMode("fire");
        t1.setDistribution(fire, Exp.fitMean(1.0 / MU));
        t1.setEnablingConditions(fire, jc, p1, 1);

        RoutingMatrix R = model.initRoutingMatrix();
        R.set(jc, jc, source, p1, 1.0);
        R.set(jc, jc, p1, t1, 1.0);
        R.set(jc, jc, t1, sink, 1.0);
        model.link(R);
        return model;
    }

    private static double[] solveP1(Network model, int cutoff) {
        SolverCTMC solver = new SolverCTMC(model, "cutoff", cutoff, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable table = solver.getAvgTable();
        List<String> names = table.getStationNames();
        int pidx = names.indexOf("P1");
        return new double[] { table.getQLen().get(pidx), table.getTput().get(pidx) };
    }

    private void assertMM1K(double[] qx, int K, String how) {
        assertEquals(mm1kMean(K), qx[0], TOL, "E[P1] (" + how + ", K=" + K + ")");
        assertEquals(mm1kTput(K), qx[1], TOL, "Tput (" + how + ", K=" + K + ")");
    }

    @Test
    public void testCapacityBoundedK8() {
        assertMM1K(solveP1(buildOpenSpn(8), 12), 8, "capacity");
    }

    @Test
    public void testCapacityBoundedK12() {
        assertMM1K(solveP1(buildOpenSpn(12), 16), 12, "capacity");
    }

    @Test
    public void testCutoffBoundedK8() {
        assertMM1K(solveP1(buildOpenSpn(-1), 8), 8, "cutoff");
    }

    @Test
    public void testCutoffBoundedK12() {
        assertMM1K(solveP1(buildOpenSpn(-1), 12), 12, "cutoff");
    }
}
