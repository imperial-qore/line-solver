package jline.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.pfqn.mva.Pfqn_bs;
import jline.api.pfqn.mva.Pfqn_linearizer;
import jline.api.pfqn.mva.Pfqn_mva;
import jline.api.pfqn.mva.Pfqn_scat;
import jline.io.Ret;
import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Node;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;

/**
 * Neuse-Chandy SCAT (Pfqn_scat).
 *
 * SCAT is Linearizer with ONE Delta refresh instead of three, so what pins it is
 * not a closed form but its POSITION: it must sit strictly between
 * Bard-Schweitzer (no refresh) and Linearizer (three), and equal neither. The
 * absolute figures asserted here are MATLAB's own, printed by pfqn_scat.m and
 * by SolverMVA with options.method = 'scat' on the same two models the C++
 * method matrix uses, so a drift in the refresh count moves them at once.
 */
public class PfqnScatTest {

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    private static Matrix mat(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                m.set(i, j, a[i][j]);
            }
        }
        return m;
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    /** Three PS stations, two closed classes, a delay. MATLAB pfqn_scat.m values. */
    @Test
    public void apiAgreesWithMatlabAndSitsBetweenBsAndLinearizer() {
        Matrix L = mat(new double[][] {{1.0, 0.5}, {0.4, 1.2}, {0.8, 0.3}});
        Matrix N = row(4.0, 3.0);
        Matrix Z = row(1.0, 2.0);
        SchedStrategy[] type = {SchedStrategy.PS, SchedStrategy.PS, SchedStrategy.PS};

        Ret.pfqnAMVA sc = Pfqn_scat.pfqn_scat(L, N, Z, type, 1e-8, 1000);
        assertEquals(0.617697394, sc.X.get(0), 1e-9);
        assertEquals(0.431056566, sc.X.get(1), 1e-9);

        Ret.pfqnAMVA lin = Pfqn_linearizer.pfqn_linearizer(L, N, Z, type, 1e-8, 1000);
        Ret.pfqnAMVA bs = Pfqn_bs.pfqn_bs(L, N, Z);
        Ret.pfqnMVA ex = Pfqn_mva.pfqn_mva(L, N, Z, Matrix.ones(1, L.getNumRows()));

        // one refresh beats none and loses to three, on every class
        for (int r = 0; r < 2; r++) {
            double eScat = Math.abs(sc.X.get(r) - ex.X.get(r));
            double eLin = Math.abs(lin.X.get(r) - ex.X.get(r));
            double eBs = Math.abs(bs.X.get(r) - ex.X.get(r));
            assertTrue(eLin < eScat, "linearizer must beat scat on class " + r);
            assertTrue(eScat < eBs, "scat must beat bard-schweitzer on class " + r);
        }
    }

    /** Delay(Z=1) + FCFS Queue(D=0.5), one closed class of 3. */
    private static Network modelA() {
        Network m = new Network("A");
        Delay d = new Delay(m, "Delay");
        Queue q = new Queue(m, "Queue", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(m, "C1", 3, d, 0);
        d.setService(c, Exp.fitRate(1.0));
        q.setService(c, Exp.fitRate(2.0));
        List<JobClass> classes = new ArrayList<JobClass>();
        classes.add(c);
        List<Node> stations = new ArrayList<Node>();
        stations.add(d);
        stations.add(q);
        RoutingMatrix P = new RoutingMatrix(m, classes, stations);
        P.addConnection(c, c, d, q, 1.0);
        P.addConnection(c, c, q, d, 1.0);
        m.link(P);
        return m;
    }

    /** The same with two servers at the queue and a population of 5. */
    private static Network modelB() {
        Network m = new Network("B");
        Delay d = new Delay(m, "Delay");
        Queue q = new Queue(m, "Queue", SchedStrategy.FCFS);
        q.setNumberOfServers(2);
        ClosedClass c = new ClosedClass(m, "C1", 5, d, 0);
        d.setService(c, Exp.fitRate(1.0));
        q.setService(c, Exp.fitRate(2.0));
        List<JobClass> classes = new ArrayList<JobClass>();
        classes.add(c);
        List<Node> stations = new ArrayList<Node>();
        stations.add(d);
        stations.add(q);
        RoutingMatrix P = new RoutingMatrix(m, classes, stations);
        P.addConnection(c, c, d, q, 1.0);
        P.addConnection(c, c, q, d, 1.0);
        m.link(P);
        return m;
    }

    private static double[] solve(Network m, String method) {
        SolverMVA s = new SolverMVA(m);
        s.options.method = method;
        NetworkAvgTable t = s.getAvgTable();
        List<Double> q = t.getQLen();
        List<Double> u = t.getUtil();
        List<Double> r = t.getRespT();
        return new double[] {q.get(0), q.get(1), u.get(1), r.get(1)};
    }

    private static void check(double[] got, double qd, double qq, double uq, double rq) {
        assertEquals(qd, got[0], 1e-8);
        assertEquals(qq, got[1], 1e-8);
        assertEquals(uq, got[2], 1e-8);
        assertEquals(rq, got[3], 1e-8);
    }

    @Test
    public void solverMethodAgreesWithMatlabOnASingleServerModel() {
        Network m = modelA();
        check(solve(m, "scat"), 1.593017294, 1.406982706, 0.796508647, 0.883218727);
        // the advertised amva.* spelling must resolve to the same algorithm, not
        // survive the alias switch and fall through to a qd-equivalent path
        check(solve(m, "amva.scat"), 1.593017294, 1.406982706, 0.796508647, 0.883218727);
        assertTrue(Arrays.asList(new SolverMVA(m).listValidMethods()).contains("scat"));
        assertTrue(Arrays.asList(new SolverMVA(m).listValidMethods()).contains("amva.scat"));
    }

    @Test
    public void solverMethodTakesAMultiserverModelThroughSeidmann() {
        // Unlike aql/qsa/tay, SCAT is not refused on a multiserver model: it
        // reaches the algorithm Seidmann-scaled, exactly as bs does.
        Network m = modelB();
        check(solve(m, "scat"), 2.883835883, 2.116164117, 0.720958971, 0.733801854);
        assertTrue(Arrays.asList(new SolverMVA(m).listValidMethods()).contains("scat"));
    }
}
