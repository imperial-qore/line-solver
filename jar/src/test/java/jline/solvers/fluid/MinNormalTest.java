package jline.solvers.fluid;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

/**
 * Regression tests for the second-order moment closure of SolverFluid
 * ({@code options.method='minnormal'}).
 *
 * <p>Every expected value is the answer of the MATLAB reference implementation
 * (SolverFLD with the same method on the same model), which is ground truth for
 * this solver. The JAR and MATLAB integrate the ODE with different solvers, so
 * they agree to integration tolerance rather than exactly; the observed spread
 * on these models is below 5e-5, and the assertions use the LINE coarse
 * tolerance of 1e-3.</p>
 *
 * <p>Mirrors the MATLAB and Python twins of the same closure.</p>
 */
public class MinNormalTest {

    private static final double TOL = 1e-3;

    private static SolverFluid solver(Network model, String method) {
        SolverOptions o = SolverFluid.defaultOptions();
        o.method = method;
        o.verbose = VerboseLevel.SILENT;
        return new SolverFluid(model, o);
    }

    /** Delay(1) -> Queue(PS, mu=1, c=2), closed, population n. */
    private static Network closedPs2(int n) {
        Network model = new Network("cqn_ps2");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        queue.setNumberOfServers(2);
        ClosedClass cclass = new ClosedClass(model, "Class1", n, delay);
        delay.setService(cclass, new Exp(1.0));
        queue.setService(cclass, new Exp(1.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** Delay(z) -> Queue(PS, mu=1, c=1), closed, N=6. */
    private static Network closedThink(double z) {
        Network model = new Network("cqn_think");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        ClosedClass cclass = new ClosedClass(model, "Class1", 6, delay);
        delay.setService(cclass, new Exp(z));
        queue.setService(cclass, new Exp(1.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** Source(lambda) -> Queue(FCFS, mu=1) -> Sink. */
    private static Network openMm1(double lambda) {
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(lambda));
        queue.setService(oclass, new Exp(1.0));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    /** Delay(1) -> Queue(sched, mu=1, c=1), two closed classes of 2 jobs, weights [1 w2]. */
    private static Network twoClassShare(SchedStrategy sched, double w2) {
        Network model = new Network("share");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", sched);
        ClosedClass c1 = new ClosedClass(model, "Class1", 2, delay);
        ClosedClass c2 = new ClosedClass(model, "Class2", 2, delay);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(1.0));
        queue.setService(c1, new Exp(1.0));
        queue.setService(c2, new Exp(1.0));
        queue.setSchedStrategyPar(c1, 1.0);
        queue.setSchedStrategyPar(c2, w2);
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** Delay -> Queue(PS, c=2): the min-normal closure corrects the mean near saturation. */
    @Test
    public void testClosedMultiserver() {
        Matrix q4 = solver(closedPs2(4), "minnormal").getAvgQLen();
        assertEquals(1.695278, q4.get(0, 0), TOL);
        assertEquals(2.304722, q4.get(1, 0), TOL);
        assertEquals(4.0, q4.get(0, 0) + q4.get(1, 0), TOL); // population conserved

        Matrix q6 = solver(closedPs2(6), "minnormal").getAvgQLen();
        assertEquals(1.960633, q6.get(0, 0), TOL);
        assertEquals(4.039367, q6.get(1, 0), TOL);
        assertEquals(6.0, q6.get(0, 0) + q6.get(1, 0), TOL);
    }

    /** Think-rate sweep across the load range where the first-order closure is worst. */
    @Test
    public void testClosedThinkSweep() {
        assertEquals(2.482226, solver(closedThink(0.25), "minnormal").getAvgQLen().get(1, 0), TOL);
        assertEquals(1.352681, solver(closedThink(0.15), "minnormal").getAvgQLen().get(1, 0), TOL);
    }

    /**
     * Open M/M/1. The EXT source pool is projected out of the covariance, so the
     * source reports no queue and flow balance makes the throughput EXACT at any
     * load, all the closure error landing on the mean.
     *
     * <p>THE rho = 0.9 VALUE WAS RE-PINNED (2026-08-31), and the number it replaces
     * was not wrong when it was written. The alternation and its inner mean solve
     * converge to {@code mom_tol = 1e-6} rather than to {@code CoarseTol}, which
     * MOVES the converged answer, and only the SATURATED cell moves by more than
     * TOL -- 0.2, 0.5 and 0.7 still hold their pre-tightening values. The old
     * 7.021525 is what the LOOSE closure settled on. 7.026368 is this solver's
     * converged answer and native python returns 7.026387 for the same model, so
     * the two still agree to 2e-5 and the pin is a stale record rather than a
     * divergence. See _kb/11-conventions-and-gotchas.md.</p>
     */
    @Test
    public void testOpenMm1() {
        double[] lambda = new double[]{0.2, 0.5, 0.7, 0.9};
        double[] expectedQ = new double[]{0.207600, 0.750274, 1.738811, 7.026368};
        for (int i = 0; i < lambda.length; i++) {
            SolverFluid s = solver(openMm1(lambda[i]), "minnormal");
            Matrix q = s.getAvgQLen();
            Matrix t = s.getAvgTput();
            assertEquals(0.0, q.get(0, 0), TOL);
            assertEquals(expectedQ[i], q.get(1, 0), TOL);
            assertEquals(lambda[i], t.get(1, 0), TOL); // exact by flow balance
        }
    }

    /**
     * DPS capacity share, closed at second order by the delta method. The share is
     * a RATIO, so evaluating it at the mean is a separate closure from the min().
     */
    @Test
    public void testDpsShare() {
        double[] w2 = new double[]{1, 2, 4, 8};
        double[] expected1 = new double[]{1.503799, 1.611981, 1.712525, 1.798943};
        double[] expected2 = new double[]{1.503799, 1.395617, 1.295069, 1.208651};
        for (int i = 0; i < w2.length; i++) {
            Matrix q = solver(twoClassShare(SchedStrategy.DPS, w2[i]), "minnormal").getAvgQLen();
            assertEquals(expected1[i], q.get(1, 0), TOL);
            assertEquals(expected2[i], q.get(1, 1), TOL);
            assertEquals(4.0, q.get(0, 0) + q.get(0, 1) + q.get(1, 0) + q.get(1, 1), TOL);
        }
    }

    /**
     * GPS splits the server by weight among the BACKLOGGED classes, so its share is
     * a function of the backlog indicator. A first-order closure cannot express it
     * at all: with continuous x_k &gt; 0 every class is always backlogged and the
     * share collapses to the constant w_k/sum_j w_j regardless of load. These
     * answers therefore depend on the second moment for their entire content.
     */
    @Test
    public void testGpsShare() {
        double[] w2 = new double[]{1, 2, 4, 8};
        double[] expected1 = new double[]{1.502339, 1.609835, 1.693551, 1.747619};
        double[] expected2 = new double[]{1.502339, 1.393981, 1.308763, 1.253704};
        for (int i = 0; i < w2.length; i++) {
            Matrix q = solver(twoClassShare(SchedStrategy.GPS, w2[i]), "minnormal").getAvgQLen();
            assertEquals(expected1[i], q.get(1, 0), TOL);
            assertEquals(expected2[i], q.get(1, 1), TOL);
            assertEquals(4.0, q.get(0, 0) + q.get(0, 1) + q.get(1, 0) + q.get(1, 1), TOL);
        }
    }

    /**
     * The DPS drift takes the class-k coordinates to w_k*x/ni of the station
     * capacity psi(xi), so equal weights must reduce DPS to PS exactly. A
     * regulariser in the denominator, or the full server count in place of
     * psi(xi), breaks this and depresses utilization.
     */
    @Test
    public void testEqualWeightDpsIsPs() {
        Matrix q = solver(twoClassShare(SchedStrategy.DPS, 1.0), "closing").getAvgQLen();
        assertEquals(q.get(1, 0), q.get(1, 1), 1e-9);
        Matrix u = solver(twoClassShare(SchedStrategy.DPS, 1.0), "closing").getAvgUtil();
        assertEquals(1.0, u.get(1, 0) + u.get(1, 1), TOL); // the shares sum to one
    }

    /** The DPS utilization split against the MATLAB reference. */
    @Test
    public void testClosingDpsUtilization() {
        Matrix u = solver(twoClassShare(SchedStrategy.DPS, 4.0), "closing").getAvgUtil();
        assertEquals(0.257334, u.get(1, 0), TOL);
        assertEquals(0.742666, u.get(1, 1), TOL);
        assertEquals(1.0, u.get(1, 0) + u.get(1, 1), TOL);
    }

    /**
     * The closure produces a genuine second moment, unlike every first-order
     * method, which computes none at all. On Delay -> Queue(PS,c=2), N=6 the
     * MATLAB reference reports QVar = 1.838762 at the queue (QStd = 1.356009).
     */
    @Test
    public void testMomentReport() {
        SolverFluid s = solver(closedPs2(6), "minnormal");
        s.getAvgQLen();
        jline.solvers.fluid.FluidResult fr = (jline.solvers.fluid.FluidResult) s.result;
        assertTrue(fr.momentQVar != null, "minnormal must report a queue-length variance");
        assertEquals(1.838762, fr.momentQVar.get(1, 0), TOL);
        assertEquals(1.356009, Math.sqrt(fr.momentQVar.get(1, 0)), TOL);
        assertTrue(fr.momentOuterIters >= 1);
    }

    /**
     * getProbAggr answers from the covariance under the closure: the probability
     * the multivariate normal puts on the unit cell around the state, folded at
     * the boundaries of the state space so the cells tile it. On
     * Delay(1) -> Queue(PS, mean 0.8), N=4 the MATLAB reference returns
     * 0.016942 / 0.100103 / 0.281035 / 0.351526 / 0.250395 for n = 0..4 at the
     * queue, and they sum to one.
     */
    @Test
    public void testProbAggrGaussianCell() {
        double[] expected = {0.016942, 0.100103, 0.281035, 0.351526, 0.250395};
        double total = 0;
        for (int n = 0; n <= 4; n++) {
            Network model = new Network("probaggr");
            Delay delay = new Delay(model, "Delay");
            Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
            ClosedClass cclass = new ClosedClass(model, "Class1", 4, delay);
            delay.setService(cclass, Exp.fitMean(1.0));
            queue.setService(cclass, Exp.fitMean(0.8));
            model.link(model.serialRouting(delay, queue));
            Matrix state = new Matrix(2, 1);
            state.set(0, 0, 4 - n);
            state.set(1, 0, n);
            model.initFromMarginal(state);

            double p = solver(model, "minnormal").getProbAggr(1).probability.get(0, 0);
            assertEquals(expected[n], p, TOL, "cell probability at n=" + n);
            total += p;
        }
        assertEquals(1.0, total, TOL, "the folded cells must tile the state space");
    }

    /**
     * An OPEN class keeps the first-order answer, which is not an independence
     * heuristic but the EXACT product form: M/M/1 at rho = 0.5 must report 0.5
     * for the empty queue, where the Gaussian cell would return 0.391.
     */
    @Test
    public void testProbAggrKeepsOpenProductForm() {
        SolverFluid s = solver(openMm1(0.5), "minnormal");
        assertEquals(0.5, s.getProbAggr(1).probability.get(0, 0), TOL);
    }

    /** With no cache node and an applicable model, 'default' now resolves to minnormal. */
    @Test
    public void testDefaultResolvesToMinnormal() {
        SolverFluid s = solver(closedPs2(4), "default");
        Matrix q = s.getAvgQLen();
        assertEquals(2.304722, q.get(1, 0), TOL);
        assertTrue(s.result.method.contains("minnormal"),
                "default should resolve to minnormal, got " + s.result.method);
    }

    /**
     * Disciplines with no branch in the closing drift must be REFUSED, not
     * integrated as an infinite server. On this model the fall-through returns
     * 2.0000 against the exact 3.0154, and the closing metric reader accepts SIRO
     * as FCFS so the error was silent. The matrix method builds a PS drift for
     * every queueing station and keeps answering.
     */
    @Test
    public void testUnsupportedDisciplinesRefused() {
        SchedStrategy[] unbranched = new SchedStrategy[]{
                SchedStrategy.SIRO, SchedStrategy.LCFS, SchedStrategy.LCFSPR};
        for (int i = 0; i < unbranched.length; i++) {
            final SchedStrategy sched = unbranched[i];
            assertThrows(RuntimeException.class,
                    () -> solver(unbranchedModel(sched), "minnormal").getAvgQLen(),
                    "minnormal must refuse " + sched);
            assertThrows(RuntimeException.class,
                    () -> solver(unbranchedModel(sched), "closing").getAvgQLen(),
                    "closing must refuse " + sched);
            assertEquals(3.0, solver(unbranchedModel(sched), "matrix").getAvgQLen().get(1, 0), TOL,
                    "matrix must still answer " + sched);
        }
    }

    private static Network unbranchedModel(SchedStrategy sched) {
        Network model = new Network("guard");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", sched);
        ClosedClass cclass = new ClosedClass(model, "Class1", 4, delay);
        delay.setService(cclass, new Exp(1.0));
        queue.setService(cclass, new Exp(1.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }
}
