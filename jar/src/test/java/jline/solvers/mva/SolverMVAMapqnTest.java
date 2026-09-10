package jline.solvers.mva;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.api.mapqn.Mapqn_amva;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.lang.processes.MAP;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

/**
 * SolverMVA method 'amva.mapqn': the horizontal-cut mean value analysis (Casale-Smirni DSN
 * 2009 balances closed by the arrival theorem) for a closed model of one exponential delay
 * and one FCFS single-server queue with a MAP service per class. Reference values are the
 * prototype's, reproduced by every codebase to solver precision; SolverCTMC bounds the
 * accuracy at the level the method has.
 */
public class SolverMVAMapqnTest {

    private static Matrix m2(double a, double b, double c, double d) {
        Matrix M = new Matrix(2, 2);
        M.set(0, 0, a); M.set(0, 1, b); M.set(1, 0, c); M.set(1, 1, d);
        return M;
    }

    private static MAP mapA() { return new MAP(m2(-3.0, 0.5, 0.2, -0.4), m2(2.5, 0.0, 0.2, 0.0)); }
    private static MAP mapB() { return new MAP(m2(-1.8, 0.3, 0.6, -0.9), m2(1.5, 0.0, 0.0, 0.3)); }

    /** delay first, queue second; every class cycles delay -> queue -> delay */
    private static Network model(double[] think, Object[] services, int[] njobs, SchedStrategy sched) {
        Network model = new Network("mapqn");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", sched);
        for (int r = 0; r < njobs.length; r++) {
            ClosedClass c = new ClosedClass(model, "C" + (r + 1), njobs[r], delay);
            delay.setService(c, new Exp(1.0 / think[r]));
            if (services[r] instanceof MAP) {
                queue.setService(c, (MAP) services[r]);
            } else {
                queue.setService(c, (Exp) services[r]);
            }
        }
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static NetworkAvgTable run(Network model, String method) {
        return new SolverMVA(model, method).getAvgTable();
    }

    @Test
    public void reproducesTheReferenceValuesOnTwoClasses() {
        Network model = model(new double[]{2.0, 8.0}, new Object[]{mapA(), mapB()}, new int[]{2, 2}, SchedStrategy.FCFS);
        NetworkAvgTable t = run(model, "amva.mapqn");
        List<Double> tput = t.getTput();
        List<Double> qlen = t.getQLen();
        // rows are station-major: delay class 1, delay class 2, queue class 1, queue class 2
        assertEquals(0.565971964, tput.get(0), 1e-7);
        assertEquals(0.1934693276, tput.get(1), 1e-7);
        assertEquals(0.8680560719, qlen.get(2), 1e-7);
        assertEquals(0.4522453793, qlen.get(3), 1e-7);
        assertEquals(1.1319439281, qlen.get(0), 1e-7);
        assertEquals(1.5477546207, qlen.get(1), 1e-7);
    }

    @Test
    public void reproducesTheReferenceValuesOnOneClass() {
        Network model = model(new double[]{2.0}, new Object[]{mapA()}, new int[]{4}, SchedStrategy.FCFS);
        NetworkAvgTable t = run(model, "amva.mapqn");
        assertEquals(0.9717108597, t.getTput().get(0), 1e-7);
        assertEquals(2.0565782806, t.getQLen().get(1), 1e-7);
    }

    @Test
    public void directApiMatches() {
        Mapqn_amva.Result r = Mapqn_amva.solve(new double[]{0.5, 0.125},
                new Matrix[]{m2(-3.0, 0.5, 0.2, -0.4), m2(-1.8, 0.3, 0.6, -0.9)},
                new Matrix[]{m2(2.5, 0.0, 0.2, 0.0), m2(1.5, 0.0, 0.0, 0.3)}, new int[]{2, 2});
        assertEquals(0.565971964, r.X[0], 1e-7);
        assertEquals(0.1934693276, r.X[1], 1e-7);
        assertEquals(0.8680560719, r.Qq[0], 1e-7);
        assertEquals(0.4522453793, r.Qq[1], 1e-7);
    }

    @Test
    public void staysWithinTheMethodAccuracyOfTheExactChain() {
        Network model = model(new double[]{2.0, 8.0}, new Object[]{mapA(), mapB()}, new int[]{2, 2}, SchedStrategy.FCFS);
        NetworkAvgTable a = run(model, "amva.mapqn");
        NetworkAvgTable e = new SolverCTMC(model).getAvgTable();
        for (int i = 0; i < 4; i++) {
            assertEquals(e.getTput().get(i), a.getTput().get(i), 0.10 * e.getTput().get(i), "Tput row " + i);
            assertEquals(e.getQLen().get(i), a.getQLen().get(i), 0.10 * e.getQLen().get(i), "QLen row " + i);
            assertEquals(e.getRespT().get(i), a.getRespT().get(i), 0.10 * e.getRespT().get(i), "RespT row " + i);
        }
    }

    @Test
    public void exponentialServiceIsExactMva() {
        Network model = model(new double[]{2.0, 8.0}, new Object[]{new Exp(1.2), new Exp(1.2)}, new int[]{3, 3}, SchedStrategy.FCFS);
        NetworkAvgTable a = run(model, "amva.mapqn");
        NetworkAvgTable e = new SolverMVA(model, "exact").getAvgTable();
        for (int i = 0; i < 4; i++) {
            assertEquals(e.getTput().get(i), a.getTput().get(i), 1e-9, "Tput row " + i);
            assertEquals(e.getQLen().get(i), a.getQLen().get(i), 1e-9, "QLen row " + i);
            assertEquals(e.getUtil().get(i), a.getUtil().get(i), 1e-9, "Util row " + i);
        }
    }

    @Test
    public void offeredAndSupportedOnItsShapeOnly() {
        Network model = model(new double[]{2.0, 8.0}, new Object[]{mapA(), mapB()}, new int[]{2, 2}, SchedStrategy.FCFS);
        SolverMVA solver = new SolverMVA(model, "amva.mapqn");
        assertTrue(Arrays.asList(solver.listValidMethods()).contains("amva.mapqn"));
        assertEquals("", solver.supportsModelMethod("amva.mapqn"));
        Network ps = model(new double[]{2.0}, new Object[]{mapA()}, new int[]{2}, SchedStrategy.PS);
        SolverMVA psSolver = new SolverMVA(ps, "amva.mapqn");
        assertFalse(Arrays.asList(psSolver.listValidMethods()).contains("amva.mapqn"));
        assertFalse(psSolver.supportsModelMethod("amva.mapqn").isEmpty());
        assertThrows(Exception.class, () -> psSolver.getAvgTable());
    }
}
