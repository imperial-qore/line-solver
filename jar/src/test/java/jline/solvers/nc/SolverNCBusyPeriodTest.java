package jline.solvers.nc;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.pfqn.Pfqn_busyp;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

/**
 * {@link SolverNC#getAvgBusyPeriod} against the analysis it wraps.
 *
 * The solver method's job is not the mathematics -- {@link Pfqn_busyp} owns
 * that, and {@code SolverLDESBusyPeriodTest} already checks it against a sample
 * path. What is tested here is the LAYER between the model and the paper's
 * inputs: the chain demands, the visit-weighted station-to-station routing, the
 * load-dependent rate function, and the removal of the Source in the open case
 * (which is not a node of the Jackson network -- its outflow is gamma). Each of
 * those is a place the port could silently produce a well-formed number for the
 * wrong network, so every case builds the paper's inputs independently and
 * requires the two to agree to machine precision.
 */
public class SolverNCBusyPeriodTest {

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /** Closed three-station network: every station, and a two-station subnetwork. */
    @Test
    public void testClosedBusyPeriodMatchesPfqnBusyp() {
        int N = 5;
        double[] rate = {1.5, 0.9, 2.0};
        double[][] Pc = {{0, 0.6, 0.4}, {0.7, 0, 0.3}, {0.5, 0.5, 0}};

        Network model = new Network("ncbpclosed");
        Queue[] q = new Queue[3];
        for (int i = 0; i < 3; i++) {
            q[i] = new Queue(model, "Q" + (i + 1), SchedStrategy.FCFS);
        }
        ClosedClass cls = new ClosedClass(model, "C1", N, q[0], 0);
        for (int i = 0; i < 3; i++) {
            q[i].setService(cls, new Exp(rate[i]));
        }
        RoutingMatrix rm = model.initRoutingMatrix();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (Pc[i][j] > 0) {
                    rm.set(cls, cls, q[i], q[j], Pc[i][j]);
                }
            }
        }
        model.link(rm);

        Matrix alpha = visits(Pc);
        Matrix mu = new Matrix(3, N);
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < N; k++) {
                mu.set(i, k, rate[i]);
            }
        }
        Matrix P = matrix(Pc);

        SolverNC solver = new SolverNC(model);
        int[][] subnets = {{0}, {1}, {2}, {0, 1}};
        int[] orders = {1, 2, 3};
        for (int[] subnet : subnets) {
            double[] got = solver.getAvgBusyPeriod(subnet, orders);
            for (int t = 0; t < orders.length; t++) {
                double exact = Pfqn_busyp.pfqn_busyp(alpha, mu, P, N, subnet, orders[t], null);
                assertEquals(exact, got[t], 1e-9 * Math.abs(exact),
                        "busy period order " + orders[t] + " of subnetwork "
                        + java.util.Arrays.toString(subnet));
            }
        }
    }

    /**
     * Open tandem. The M/M/1 busy period of ANY order lasts 1/(mu-lambda), which
     * pins the first station without reference to the implementation at all.
     */
    @Test
    public void testOpenBusyPeriodMatchesPfqnBusyp() {
        double lambda = 0.5;
        double[] rate = {1.0, 1.2};

        Network model = new Network("ncbpopen");
        Source src = new Source(model, "Source");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Sink snk = new Sink(model, "Sink");
        OpenClass cls = new OpenClass(model, "C1");
        src.setArrival(cls, new Exp(lambda));
        q1.setService(cls, new Exp(rate[0]));
        q2.setService(cls, new Exp(rate[1]));
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(cls, cls, src, q1, 1.0);
        rm.set(cls, cls, q1, q2, 1.0);
        rm.set(cls, cls, q2, snk, 1.0);
        model.link(rm);

        // the paper's inputs, with the Source absent and its stream carried by gamma
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, lambda);
        alpha.set(0, 1, lambda);
        Matrix mu = new Matrix(2, 1);
        mu.set(0, 0, rate[0]);
        mu.set(1, 0, rate[1]);
        Matrix P = new Matrix(2, 2);
        P.set(0, 1, 1.0);
        Matrix gamma = new Matrix(1, 2);
        gamma.set(0, 0, lambda);

        // station 0 of the model is the Source, so Q1 and Q2 are stations 1 and 2
        SolverNC solver = new SolverNC(model);
        int[] orders = {1, 2, 3};
        double[] single = solver.getAvgBusyPeriod(new int[]{1}, orders);
        double[] joint = solver.getAvgBusyPeriod(new int[]{1, 2}, orders);
        for (int t = 0; t < orders.length; t++) {
            assertEquals(1.0 / (rate[0] - lambda), single[t], 1e-6,
                    "M/M/1 busy period order " + orders[t]);
            double exactSingle = Pfqn_busyp.pfqn_busyp(alpha, mu, P,
                    Double.POSITIVE_INFINITY, new int[]{0}, orders[t], gamma);
            assertEquals(exactSingle, single[t], 1e-9 * Math.abs(exactSingle),
                    "Q1 busy period order " + orders[t]);
            double exactJoint = Pfqn_busyp.pfqn_busyp(alpha, mu, P,
                    Double.POSITIVE_INFINITY, new int[]{0, 1}, orders[t], gamma);
            assertEquals(exactJoint, joint[t], 1e-9 * Math.abs(exactJoint),
                    "tandem busy period order " + orders[t]);
        }
    }

    /**
     * A DELAY STATION IN THE SUBNETWORK is what the rate-function form of the
     * API exists for: its rate is k/S and grows without bound, so the dense mu
     * table the Matrix overloads read -- which reuses its last column past its
     * width -- cannot express it. The reference value here is built from a table
     * wide enough that the clamp never fires, so the two must agree exactly; a
     * port that passed a short table would disagree at the higher orders.
     */
    @Test
    public void testClosedBusyPeriodWithDelayStation() {
        int N = 6;
        double delayRate = 0.8;
        double queueRate = 1.4;
        double[][] Pc = {{0, 1.0}, {1.0, 0}};

        Network model = new Network("ncbpdelay");
        Delay think = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.FCFS);
        ClosedClass cls = new ClosedClass(model, "C1", N, think, 0);
        think.setService(cls, new Exp(delayRate));
        q.setService(cls, new Exp(queueRate));
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(cls, cls, think, q, 1.0);
        rm.set(cls, cls, q, think, 1.0);
        model.link(rm);

        Matrix alpha = visits(Pc);
        Matrix mu = new Matrix(2, N);
        for (int k = 0; k < N; k++) {
            mu.set(0, k, (k + 1) * delayRate);   // infinite server: k jobs, rate k/S
            mu.set(1, k, queueRate);
        }
        Matrix P = matrix(Pc);

        SolverNC solver = new SolverNC(model);
        int[] orders = {1, 2, 3, 4};
        double[] got = solver.getAvgBusyPeriod(new int[]{0}, orders);
        for (int t = 0; t < orders.length; t++) {
            double exact = Pfqn_busyp.pfqn_busyp(alpha, mu, P, N, new int[]{0}, orders[t], null);
            assertEquals(exact, got[t], 1e-9 * Math.abs(exact),
                    "delay-station busy period order " + orders[t]);
        }
    }

    /**
     * A closed subnetwork must be a PROPER subset: with every station in it no
     * job can enter from outside, so no busy period ever starts and the answer
     * is undefined rather than infinite.
     */
    @Test
    public void testWholeNetworkSubnetIsRefused() {
        int N = 3;
        Network model = new Network("ncbpwhole");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass cls = new ClosedClass(model, "C1", N, q1, 0);
        q1.setService(cls, new Exp(1.0));
        q2.setService(cls, new Exp(1.0));
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(cls, cls, q1, q2, 1.0);
        rm.set(cls, cls, q2, q1, 1.0);
        model.link(rm);

        SolverNC solver = new SolverNC(model);
        assertThrows(IllegalArgumentException.class,
                () -> solver.getAvgBusyPeriod(new int[]{0, 1}, new int[]{1}));
    }

    /** Stochastic solution of x = x*P, the paper's relative arrival rates. */
    private static Matrix visits(double[][] Pc) {
        int n = Pc.length;
        double[] x = new double[n];
        java.util.Arrays.fill(x, 1.0 / n);
        for (int it = 0; it < 100000; it++) {
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    y[j] += x[i] * Pc[i][j];
                }
            }
            double s = 0;
            for (int i = 0; i < n; i++) {
                s += y[i];
            }
            for (int i = 0; i < n; i++) {
                y[i] /= s;
            }
            x = y;
        }
        Matrix m = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            m.set(0, i, x[i]);
        }
        return m;
    }

    private static Matrix matrix(double[][] v) {
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                m.set(i, j, v[i][j]);
            }
        }
        return m;
    }
}
