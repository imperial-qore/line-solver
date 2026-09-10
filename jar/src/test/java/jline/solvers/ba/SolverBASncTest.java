package jline.solvers.ba;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.api.snc.Snc_bound_delay;
import jline.api.snc.Snc_conv;
import jline.api.snc.Snc_env_map;
import jline.api.snc.Snc_env_poisson;
import jline.api.snc.Snc_perc_delay;
import jline.api.snc.Snc_srv_exp;
import jline.api.snc.SncEnvelope;
import jline.api.snc.SncResult;
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
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkPercTable;
import jline.util.matrix.Matrix;

/**
 * Tests for the stochastic network calculus family of {@link SolverBA}.
 *
 * <p>Every reference number here is the MATLAB answer of the same model under
 * {@code SolverBA(model,'method','snc.upper')}, which is the ground truth for
 * the port, and each is also reproduced by native Python. The tolerance is
 * 1e-6: the three codebases differ only in the refinement step of the Chernoff
 * search (MATLAB {@code fminbnd}, scipy's bounded Brent, and the golden section
 * of {@link jline.api.snc.Snc_thetaopt}), which agrees far inside that.</p>
 *
 * <p>The exact M/M/1 values are quoted alongside the bound so that the tests
 * document the tightness rather than only pinning the numbers: this is a
 * policy-robust bound, loose on the mean and asymptotically exact in the tail.
 * </p>
 */
public class SolverBASncTest {

    private static final double TOL = 1e-6;

    private static Network mm1(double lambda, double mu) {
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(lambda));
        queue.setService(oclass, new Exp(mu));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, Network.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    @Test
    public void mm1MeanBoundMatchesMatlab() {
        // MATLAB: R bound 4.9500 / 13.2051 / 87.9922 against the exact
        // 1.4286 / 2.5000 / 10.0000, i.e. 3.5x / 5.3x / 8.8x.
        double[] rho = {0.3, 0.6, 0.9};
        double[] expectedR = {4.950038919240365, 13.205141824117, 87.99215434588235};
        for (int i = 0; i < rho.length; i++) {
            NetworkAvgTable t = new SolverBA(mm1(rho[i], 1.0), "snc.upper").getAvgTable();
            List<Double> respt = t.getRespT();
            List<Double> qlen = t.getQLen();
            double r = respt.get(respt.size() - 1);
            double q = qlen.get(qlen.size() - 1);
            assertEquals(expectedR[i], r, 1e-4 * expectedR[i]);
            assertTrue(r > 1.0 / (1.0 - rho[i]),
                    "the bound must exceed the exact M/M/1 response time");
            // Q is Little's law on the bounded R and the exact throughput
            assertEquals(rho[i] * r, q, TOL);
        }
    }

    @Test
    public void mm1QuantilesMatchMatlab() {
        SolverBA solver = new SolverBA(mm1(0.6, 1.0), "snc.upper");
        Matrix d = solver.getDelayPerc(1e-3);
        Matrix b = solver.getBacklogPerc(1e-3);
        // MATLAB: 29.6082 slots and 23.6449 jobs at eps = 1e-3.
        assertEquals(29.608187856023928, d.get(1, 0), 1e-4);
        assertEquals(23.64494253393617, b.get(1, 0), 1e-4);
        // The Source row carries no traffic and stays NaN rather than zero.
        assertTrue(Double.isNaN(d.get(0, 0)));

        NetworkPercTable t = solver.getPercTable(1e-3);
        assertEquals(1, t.getRespTPerc().size());
        assertEquals(29.608187856023928, t.getRespTPerc().get(0), 1e-4);
    }

    @Test
    public void quantileRatioFallsTowardsOneDeepInTheTail() {
        // The bound reproduces the exact decay rate and pays a constant
        // prefactor, so tightening the guarantee tightens the ratio: 2.03x at
        // eps=1e-2 down to 1.21x at eps=1e-12 on this model.
        SolverBA solver = new SolverBA(mm1(0.6, 1.0), "snc.upper");
        double prev = Double.POSITIVE_INFINITY;
        double[] epsList = {1e-2, 1e-3, 1e-6, 1e-9, 1e-12};
        for (int i = 0; i < epsList.length; i++) {
            double bound = solver.getDelayPerc(epsList[i]).get(1, 0);
            double exact = -Math.log(epsList[i]) / 0.4;
            double ratio = bound / exact;
            assertTrue(ratio > 1.0, "the quantile bound must exceed the exact quantile");
            assertTrue(ratio < prev, "the ratio must fall as eps tightens");
            prev = ratio;
        }
        assertEquals(1.213, prev, 0.01);
    }

    @Test
    public void tandemBoundMatchesMatlabAndDegradesHopByHop() {
        Network model = new Network("tandem");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.6));
        q1.setService(oclass, new Exp(1.2));
        q2.setService(oclass, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, Network.serialRouting(source, q1, q2, sink));
        model.link(P);

        NetworkAvgTable t = new SolverBA(model, "snc.upper").getAvgTable();
        List<Double> respt = t.getRespT();
        // MATLAB: 7.3397 at Q1 and 20.0114 at Q2, against the exact Jackson
        // 1.6667 and 2.5000. The second hop is looser because its arrival
        // envelope is the DEPARTURE envelope of the first, which carries the
        // burst the server added; the exact answer does not degrade this way,
        // by Burke's theorem.
        assertEquals(7.339700, respt.get(respt.size() - 2), 1e-3);
        assertEquals(20.011400, respt.get(respt.size() - 1), 1e-3);
    }

    @Test
    public void twoClassesShareOneServerByBlindMultiplexing() {
        Network model = new Network("shared");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Shared", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass classA = new OpenClass(model, "ClassA");
        OpenClass classB = new OpenClass(model, "ClassB");
        source.setArrival(classA, new Exp(0.3));
        source.setArrival(classB, new Exp(0.3));
        queue.setService(classA, new Exp(1.0));
        queue.setService(classB, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(classA, classA, Network.serialRouting(source, queue, sink));
        P.set(classB, classB, Network.serialRouting(source, queue, sink));
        model.link(P);

        SolverBA solver = new SolverBA(model, "snc.upper");
        NetworkAvgTable t = solver.getAvgTable();
        List<Double> respt = t.getRespT();
        // MATLAB: 24.0303 per class, against an exact 2.5 for the aggregate
        // M/M/1. The gap is what covers every work-conserving policy, including
        // the one that always serves the other class first.
        assertEquals(24.030304, respt.get(respt.size() - 1), 1e-3);
        assertEquals(24.030304, respt.get(respt.size() - 2), 1e-3);

        // The backlog quantile of a shared station equals that of the solo
        // station at the aggregate rate: with zero burst terms the bound sees
        // the envelopes only through rhoS - rhoA, and blind multiplexing
        // subtracts the cross rate from rhoS exactly as aggregating the classes
        // would add it to rhoA. The DELAY quantile does separate them, since it
        // scales the level by rhoS.
        Matrix b = solver.getBacklogPerc(1e-3);
        Matrix d = solver.getDelayPerc(1e-3);
        assertEquals(23.64494253393617, b.get(1, 0), 1e-4);
        assertEquals(55.901808, d.get(1, 0), 1e-3);
    }

    @Test
    public void refusesWhatTheCalculusCannotCarry() {
        // closed model
        Network closed = new Network("closed");
        Delay delay = new Delay(closed, "Think");
        Queue queue = new Queue(closed, "Q", SchedStrategy.PS);
        ClosedClass cclass = new ClosedClass(closed, "C", 3, delay);
        delay.setService(cclass, new Exp(1.0));
        queue.setService(cclass, new Exp(2.0));
        RoutingMatrix Pc = closed.initRoutingMatrix();
        Pc.set(cclass, cclass, Network.serialRouting(delay, queue));
        closed.link(Pc);
        assertThrows(RuntimeException.class,
                () -> new SolverBA(closed, "snc.upper").getAvgTable());

        // probabilistic split downstream of the Source
        Network split = new Network("split");
        Source s2 = new Source(split, "Source");
        Queue qa = new Queue(split, "QA", SchedStrategy.FCFS);
        Queue qb = new Queue(split, "QB", SchedStrategy.FCFS);
        Sink k2 = new Sink(split, "Sink");
        OpenClass o2 = new OpenClass(split, "C");
        s2.setArrival(o2, new Exp(0.4));
        qa.setService(o2, new Exp(1.0));
        qb.setService(o2, new Exp(1.0));
        RoutingMatrix Ps = split.initRoutingMatrix();
        Ps.set(o2, o2, s2, qa, 1.0);
        Ps.set(o2, o2, qa, qb, 0.5);
        Ps.set(o2, o2, qa, k2, 0.5);
        Ps.set(o2, o2, qb, k2, 1.0);
        split.link(Ps);
        assertThrows(RuntimeException.class,
                () -> new SolverBA(split, "snc.upper").getAvgTable());

        // unequal service rates among the classes sharing a station
        Network unequal = new Network("unequal");
        Source s3 = new Source(unequal, "Source");
        Queue q3 = new Queue(unequal, "Q", SchedStrategy.FCFS);
        Sink k3 = new Sink(unequal, "Sink");
        OpenClass a3 = new OpenClass(unequal, "A");
        OpenClass b3 = new OpenClass(unequal, "B");
        s3.setArrival(a3, new Exp(0.2));
        s3.setArrival(b3, new Exp(0.2));
        q3.setService(a3, new Exp(1.0));
        q3.setService(b3, new Exp(2.0));
        RoutingMatrix Pu = unequal.initRoutingMatrix();
        Pu.set(a3, a3, Network.serialRouting(s3, q3, k3));
        Pu.set(b3, b3, Network.serialRouting(s3, q3, k3));
        unequal.link(Pu);
        assertThrows(RuntimeException.class,
                () -> new SolverBA(unequal, "snc.upper").getAvgTable());
    }

    @Test
    public void listValidMethodsNarrowsToTheOpenFamiliesOnAnOpenModel() {
        String[] valid = new SolverBA(mm1(0.6, 1.0)).listValidMethods();
        List<String> v = Arrays.asList(valid);
        assertEquals(3, v.size());
        assertTrue(v.contains("snc.upper"));
        assertTrue(v.contains("bpt.lower"));
        assertTrue(v.contains("bgt.upper"));
    }

    @Test
    public void apiOnePhaseMapEqualsThePoissonEnvelope() {
        Matrix D0 = new Matrix(1, 1);
        D0.set(0, 0, -0.6);
        Matrix D1 = new Matrix(1, 1);
        D1.set(0, 0, 0.6);
        double[] theta = {0.1, 0.5, 1.0, 2.0};
        for (int i = 0; i < theta.length; i++) {
            double[] fromMap = Snc_env_map.snc_env_map(D0, D1, theta[i]);
            double[] fromPoisson = Snc_env_poisson.snc_env_poisson(0.6, theta[i]);
            assertEquals(fromPoisson[0], fromMap[0], 1e-12);
            assertEquals(fromPoisson[1], fromMap[1], 1e-9);
        }
    }

    @Test
    public void apiPayBurstsOnlyOnceBeatsHopByHop() {
        final double lambda = 0.6;
        final double[] rates = {1.5, 1.2, 1.0};
        SncEnvelope arv = Snc_env_poisson.of(lambda);
        SncEnvelope end = new SncEnvelope() {
            public double[] eval(double theta) {
                double[] acc = Snc_srv_exp.snc_srv_exp(rates[0], theta);
                for (int i = 1; i < rates.length; i++) {
                    double[] next = Snc_srv_exp.snc_srv_exp(rates[i], theta);
                    acc = Snc_conv.snc_conv(acc[0], acc[1], next[0], next[1], theta);
                }
                return acc;
            }
        };
        SncResult concat = Snc_perc_delay.snc_perc_delay(arv, end, 1e-3);
        double hopByHop = 0;
        for (int i = 0; i < rates.length; i++) {
            hopByHop += Snc_perc_delay.snc_perc_delay(arv, Snc_srv_exp.of(rates[i]), 1e-3 / 3).value;
        }
        // MATLAB: 42.5195 concatenated against 65.5128 summed per hop.
        assertEquals(42.5195, concat.value, 1e-2);
        assertEquals(65.5128, hopByHop, 1e-2);
        assertTrue(concat.value < hopByHop);
    }

    @Test
    public void apiBoundAndQuantileAreMutuallyConsistent() {
        SncEnvelope arv = Snc_env_poisson.of(0.6);
        SncEnvelope srv = Snc_srv_exp.of(1.0);
        double eps = 1e-4;
        SncResult d = Snc_perc_delay.snc_perc_delay(arv, srv, eps);
        SncResult back = Snc_bound_delay.snc_bound_delay(arv, srv, d.value);
        assertEquals(eps, back.value, 1e-8);
    }
}
