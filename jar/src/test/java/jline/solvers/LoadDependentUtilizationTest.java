package jline.solvers;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.mam.MAMOptions;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * One load-dependent fixture, every solver, ONE utilization convention.
 *
 * <p>Utilization has two readings that coincide for a load-independent single
 * server and part company the moment alpha(n) != 1:
 *
 * <pre>
 *   work-based   U = X * E[S] / max(c, max alpha)   the fraction of the
 *                                                   station's PEAK capacity
 *                                                   actually being delivered
 *   time-based   U = P(at least one server busy)    the fraction of time the
 *                                                   station is occupied
 * </pre>
 *
 * <p>A server running alpha(n) times faster does the same work in less time, so
 * the time-based reading calls it no busier than one at its nominal rate. On the
 * fixture below that is 0.9587 against 0.6612, a 45% spread on the same model.
 *
 * <p>LINE reports the WORK-BASED number everywhere. It used to be split: MAM's
 * LD-QBD, the NRM SSA engine and the C++ LDES engine measured busy time while
 * CTMC, MVA, NC, serial SSA and the Java LDES engine measured work, so SolverSSA
 * disagreed with itself depending on which engine ran.
 */
public class LoadDependentUtilizationTest {

    private static final int QI = 1;        // 0 = Delay, 1 = Queue
    private static final int N = 4;
    private static final double[] ALPHA = {1.0, 1.5, 2.0, 2.5};

    /**
     * Closed form for the fixture: the queue is a birth-death chain on n = 0..4
     * with birth (N-n)*1.0 and death alpha(n)*1.0.
     */
    private static final double X_EXACT = 1.6528925619834711;
    private static final double U_EXACT = X_EXACT / 2.5;    // 0.6611570247933884
    private static final double Q_EXACT = 2.3471074380165293;
    private static final double TOL = 1e-9;

    private static Network model() {
        Network model = new Network("lld_util");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass jobclass = new ClosedClass(model, "Class1", N, delay);
        delay.setService(jobclass, new Exp(1.0));
        queue.setService(jobclass, new Exp(1.0));
        Matrix alpha = new Matrix(1, ALPHA.length);
        for (int i = 0; i < ALPHA.length; i++) {
            alpha.set(0, i, ALPHA[i]);
        }
        queue.setLoadDependence(alpha);
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static void assertQueueRow(NetworkSolver solver, String what, double tol) {
        assertEquals(Q_EXACT, solver.getAvgQLen().get(QI, 0), tol, what + ": QLen");
        assertEquals(U_EXACT, solver.getAvgUtil().get(QI, 0), tol, what + ": Util");
        assertEquals(X_EXACT, solver.getAvgTput().get(QI, 0), tol, what + ": Tput");
    }

    @Test
    public void testTheClosedFormIsWhatCtmcComputes() {
        // anchors the other assertions to arithmetic rather than to a solver
        SolverCTMC solver = new SolverCTMC(model());
        assertQueueRow(solver, "CTMC", TOL);
        // and it is emphatically NOT the time-based reading
        assertTrue(Math.abs(solver.getAvgUtil().get(QI, 0) - 0.9586776860) > 0.29,
                "CTMC must not report P(busy) here");
    }

    @Test
    public void testMvaAgrees() {
        assertQueueRow(new SolverMVA(model()), "MVA", TOL);
    }

    @Test
    public void testNcAgrees() {
        assertQueueRow(new SolverNC(model()), "NC", TOL);
    }

    @Test
    public void testMamLdqbdAgrees() {
        // an exact chain like CTMC's, so the bar is not a percentage. This is
        // the regression: it used to report 1 - pi(0) = 0.9587 here.
        MAMOptions options = new MAMOptions();
        options.method = "ldqbd";
        assertQueueRow(new SolverMAM(model(), options), "MAM ldqbd", TOL);
    }

    @Test
    public void testSsaAgrees() {
        // the NRM engine integrates busy time, which is the other convention
        // until the load-dependent stations are overridden with T*S/peak
        SolverSSA solver = new SolverSSA(model(), "samples", 20000, "seed", 23000);
        assertEquals(U_EXACT, solver.getAvgUtil().get(QI, 0), 3e-2, "SSA: Util");
        assertEquals(Q_EXACT, solver.getAvgQLen().get(QI, 0), 8e-2, "SSA: QLen");
    }

    @Test
    public void testSsaParallelAgrees() {
        // the replica analyzer divided by the SERVER COUNT alone and never read
        // lldscaling, so it reported 1.6529 here -- a utilization above one
        SolverSSA solver = new SolverSSA(model(), "method", "para", "samples", 20000, "seed", 23000);
        double u = solver.getAvgUtil().get(QI, 0);
        assertTrue(u <= 1.0, "utilization above one: " + u);
        assertEquals(U_EXACT, u, 3e-2, "SSA para: Util");
    }

    @Test
    public void testUtilizationStaysBelowOneAtSaturation() {
        // the time-based reading saturates at 1 while the station still has
        // capacity left; the work-based one reaches 1 only at alpha's peak
        Network m = new Network("lld_sat");
        Delay delay = new Delay(m, "Delay");
        Queue queue = new Queue(m, "Queue", SchedStrategy.FCFS);
        ClosedClass jobclass = new ClosedClass(m, "Class1", 6, delay);
        delay.setService(jobclass, new Exp(100.0));    // a near-zero think time
        queue.setService(jobclass, new Exp(1.0));
        Matrix alpha = new Matrix(1, 6);
        for (int i = 0; i < 6; i++) {
            alpha.set(0, i, i + 1.0);
        }
        queue.setLoadDependence(alpha);
        m.link(m.serialRouting(delay, queue));

        double uCtmc = new SolverCTMC(m).getAvgUtil().get(QI, 0);
        MAMOptions options = new MAMOptions();
        options.method = "ldqbd";
        double uMam = new SolverMAM(m, options).getAvgUtil().get(QI, 0);
        assertEquals(uCtmc, uMam, TOL, "saturated load-dependent Util");
        assertTrue(uMam > 0.0 && uMam < 1.0, "Util out of range: " + uMam);
    }
}
