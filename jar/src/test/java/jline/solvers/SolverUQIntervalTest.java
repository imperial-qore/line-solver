package jline.solvers;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import jline.api.pfqn.mva.Pfqn_mva_interval;
import jline.lang.ClosedClass;
import jline.GlobalConstants;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.VerboseLevel;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Prior;
import jline.solvers.mva.SolverMVA;
import jline.solvers.uq.SolverUQ;
import jline.util.matrix.Matrix;

/**
 * The support-only interval reported by SolverUQ must agree with the api-level hull of
 * Luthi and Haring (1998), and must fall back to the sampled range when the monotonicity
 * theorems do not apply.
 */
public class SolverUQIntervalTest {

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    /** Terminals Z = 15, D_cpu = [12,16] as a Prior, D_disk = 10, n = 10. */
    private static Network cpuDiskModel() {
        Network model = new Network("IntervalUQ");
        Delay terminals = new Delay(model, "Terminals");
        Queue cpu = new Queue(model, "CPU", SchedStrategy.PS);
        Queue disk = new Queue(model, "Disk", SchedStrategy.PS);
        ClosedClass jobs = new ClosedClass(model, "Jobs", 10, terminals);
        terminals.setService(jobs, new Exp(1.0 / 15));
        List<Distribution> alts = Arrays.asList((Distribution) new Exp(1.0 / 12), new Exp(1.0 / 16));
        cpu.setService(jobs, new Prior(alts, new double[]{0.5, 0.5}));
        disk.setService(jobs, new Exp(1.0 / 10));
        model.link(model.serialRouting(terminals, cpu, disk));
        return model;
    }

    @Test
    public void intervalMatchesTheApiHull() {
        SolverUQ solver = new SolverUQ(cpuDiskModel(), new SolverUQ.SolverFactory() {
            public NetworkSolver create(Network m) {
                return new SolverMVA(m);
            }
        });
        assertNull(solver.qualifiesForIntervalMVA());

        SolverUQ.Interval ival = solver.getInterval();
        assertTrue(ival.exact);
        assertEquals("mvainterval", ival.method);

        Matrix L = new Matrix(2, 2);
        L.set(0, 0, 12);
        L.set(0, 1, 16);
        L.set(1, 0, 10);
        L.set(1, 1, 10);
        Pfqn_mva_interval.Result ref = Pfqn_mva_interval.pfqn_mva_interval(L, row(10), row(15));

        assertEquals(ref.X.get(0, 0), ival.X.get(0, 0), 1e-12);
        assertEquals(ref.X.get(0, 1), ival.X.get(0, 1), 1e-12);
        // station 0 is the delay, stations 1 and 2 are the CPU and the disk
        assertEquals(ref.Q.get(0, 0), ival.Qlo.get(1, 0), 1e-12);
        assertEquals(ref.Q.get(0, 1), ival.Qup.get(1, 0), 1e-12);
        assertEquals(ref.Q.get(1, 0), ival.Qlo.get(2, 0), 1e-12);
        assertEquals(ref.Q.get(1, 1), ival.Qup.get(2, 0), 1e-12);
        assertEquals(15.0, ival.Rlo.get(0, 0), 1e-12);
        assertEquals(15.0, ival.Rup.get(0, 0), 1e-12);
        for (int i = 0; i < 3; i++) {
            assertTrue(ival.Tlo.get(i, 0) <= ival.Tup.get(i, 0) + 1e-12);
        }
    }

    /** Every solved alternative must sit inside the exact hull. */
    @Test
    public void everyAlternativeLiesInsideTheHull() {
        SolverUQ solver = new SolverUQ(cpuDiskModel(), new SolverUQ.SolverFactory() {
            public NetworkSolver create(Network m) {
                return new SolverMVA(m);
            }
        });
        SolverUQ.Interval exact = solver.getInterval();
        SolverUQ.Interval sampled = solver.intervalBySampling("forced");
        assertFalse(sampled.exact);
        for (int i = 0; i < 3; i++) {
            assertTrue(sampled.Qlo.get(i, 0) >= exact.Qlo.get(i, 0) - 1e-9);
            assertTrue(sampled.Qup.get(i, 0) <= exact.Qup.get(i, 0) + 1e-9);
        }
    }

    /** A multiclass model is refused by the gate and falls back to the sampled range. */
    @Test
    public void multiclassFallsBackToTheSampledRange() {
        Network model = new Network("TwoClass");
        Delay terminals = new Delay(model, "Terminals");
        Queue cpu = new Queue(model, "CPU", SchedStrategy.PS);
        ClosedClass a = new ClosedClass(model, "A", 4, terminals);
        ClosedClass b = new ClosedClass(model, "B", 3, terminals);
        terminals.setService(a, new Exp(1.0 / 5));
        terminals.setService(b, new Exp(1.0 / 5));
        List<Distribution> alts = Arrays.asList((Distribution) new Exp(1.0 / 2), new Exp(1.0 / 3));
        cpu.setService(a, new Prior(alts, new double[]{0.5, 0.5}));
        cpu.setService(b, new Exp(1.0 / 2));
        model.link(model.serialRouting(terminals, cpu));

        SolverUQ solver = new SolverUQ(model, new SolverUQ.SolverFactory() {
            public NetworkSolver create(Network m) {
                return new SolverMVA(m);
            }
        });
        assertEquals("the theorems are proved for a single class only", solver.qualifiesForIntervalMVA());
        SolverUQ.Interval ival = solver.getInterval();
        assertFalse(ival.exact);
        assertEquals("sampled", ival.method);
        assertNotNull(ival.reason);
        assertEquals(2, ival.Qlo.getNumRows());
        assertEquals(2, ival.Qlo.getNumCols());
    }
}
