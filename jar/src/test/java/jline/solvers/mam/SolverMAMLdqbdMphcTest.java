package jline.solvers.mam;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Distribution;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The three things the JAR LD-QBD path could not do before 2026-08-18.
 *
 * <p>1. EXACT MULTISERVER PH. The c parallel PH servers used to be collapsed
 * into one PH process run at min(n,c) times its speed, which keeps the aggregate
 * service rate right but forgets which phase each busy server is in (~1e-2
 * relative against SolverCTMC). The level now carries the MULTISET of the phases
 * the min(n,c) busy servers sit in, so the chain is exact and the bar here is
 * machine precision, not a percentage.
 *
 * <p>2. LOAD DEPENDENCE. {@code SolverMAM.methodFeatureSet} declares
 * LoadDependence for 'default'/'ldqbd', but the JAR's own copy of the block
 * construction never read sn.lldscaling: it accepted the model and returned the
 * NON-load-dependent answer (measured 0.98462 against the true 1.65290 on the
 * fixture below, a 68% throughput error, with no warning). Both paths now share
 * one builder.
 *
 * <p>3. THE OPEN REGIME. MATLAB, Python and C++ all serve a Source+Queue model
 * by truncating the level space; the JAR refused it by name.
 *
 * <p>The oracle is SolverCTMC throughout, exact for phase-type service on these
 * shapes, and all four metrics are compared including under load dependence.
 * MAM used to report P(busy) = 1 - pi(0) there, which reads a server running
 * alpha(n) times faster as no busier than one at its nominal rate: 0.9587
 * against CTMC's 0.6612 on the fixture below. It now reports the work-based
 * sum_n pi(n)*sf(n)/max(c, max(alpha)), CTMC's own convention.
 */
public class SolverMAMLdqbdMphcTest {

    private static final double TOL = 1e-9;   // exact vs exact
    private static final int QI = 1;          // 0 = Delay/Source, 1 = Queue

    private static Network closedModel(int N, double lamDelay, Distribution service,
                                       int servers, double[] lld) {
        Network model = new Network("ldqbd_mphc");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(servers);
        ClosedClass jobclass = new ClosedClass(model, "Class1", N, delay);
        delay.setService(jobclass, new Exp(lamDelay));
        queue.setService(jobclass, service);
        if (lld != null) {
            Matrix alpha = new Matrix(1, lld.length);
            for (int i = 0; i < lld.length; i++) {
                alpha.set(0, i, lld[i]);
            }
            queue.setLoadDependence(alpha);
        }
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** Queue-row (QLen, Util, RespT, Tput) from the LDQBD method. */
    private static double[] mam(Network model) {
        MAMOptions options = new MAMOptions();
        options.method = "ldqbd";
        SolverMAM solver = new SolverMAM(model, options);
        return new double[]{solver.getAvgQLen().get(QI, 0), solver.getAvgUtil().get(QI, 0),
                            solver.getAvgRespT().get(QI, 0), solver.getAvgTput().get(QI, 0)};
    }

    private static double[] ctmc(Network model) {
        SolverCTMC solver = new SolverCTMC(model);
        return new double[]{solver.getAvgQLen().get(QI, 0), solver.getAvgUtil().get(QI, 0),
                            solver.getAvgRespT().get(QI, 0), solver.getAvgTput().get(QI, 0)};
    }

    private static void assertMatchesCtmc(Network model, String what) {
        double[] got = mam(model);
        double[] want = ctmc(model);
        String[] name = {"QLen", "Util", "RespT", "Tput"};
        for (int i = 0; i < got.length; i++) {
            assertEquals(want[i], got[i], TOL, what + ": " + name[i]);
        }
    }

    @Test
    public void testErlangMultiserverIsExact() {
        // Erlang-2, mean 1 (two phases at rate 2): the collapsed chain was ~1e-2 off
        for (int c = 1; c <= 3; c++) {
            assertMatchesCtmc(closedModel(5, 0.8, new Erlang(2.0, 2), c, null),
                    "Erlang-2 with c=" + c);
        }
    }

    @Test
    public void testErlang3MultiserverIsExact() {
        for (int c = 1; c <= 3; c++) {
            assertMatchesCtmc(closedModel(4, 0.9, new Erlang(3.0, 3), c, null),
                    "Erlang-3 with c=" + c);
        }
    }

    @Test
    public void testHyperExpMultiserverIsExact() {
        for (int c = 1; c <= 3; c++) {
            assertMatchesCtmc(closedModel(5, 0.7, new HyperExp(0.6, 2.0, 0.5), c, null),
                    "HyperExp with c=" + c);
        }
    }

    @Test
    public void testExponentialMultiserverUnchanged() {
        // the exponential path never enters the multiset builder and must keep
        // matching the M/M/c boundary exactly
        for (int c = 1; c <= 4; c++) {
            assertMatchesCtmc(closedModel(6, 0.5, new Exp(1.0), c, null),
                    "Exp with c=" + c);
        }
    }

    @Test
    public void testLoadDependenceIsHonoured() {
        // the regression: the JAR used to drop sn.lldscaling on the floor here
        double[] lld = {1.0, 1.5, 2.0, 2.5};
        Network model = closedModel(4, 1.0, new Exp(1.0), 1, lld);
        double[] got = mam(model);
        double[] want = ctmc(model);
        assertEquals(want[0], got[0], TOL, "load-dependent QLen");
        assertEquals(want[1], got[1], TOL, "load-dependent Util");
        assertEquals(want[2], got[2], TOL, "load-dependent RespT");
        assertEquals(want[3], got[3], TOL, "load-dependent Tput");
        // and it is emphatically NOT the load-independent answer any more
        double[] noLld = mam(closedModel(4, 1.0, new Exp(1.0), 1, null));
        assertTrue(Math.abs(got[3] - noLld[3]) > 0.5,
                "load-dependent throughput must differ from the unscaled one, got "
                        + got[3] + " against " + noLld[3]);
    }

    @Test
    public void testLoadDependenceWithPhService() {
        double[] lld = {1.0, 1.5, 2.0, 2.5};
        Network model = closedModel(4, 1.0, new Erlang(2.0, 2), 1, lld);
        double[] got = mam(model);
        double[] want = ctmc(model);
        assertEquals(want[0], got[0], TOL, "load-dependent PH QLen");
        assertEquals(want[1], got[1], TOL, "load-dependent PH Util");
        assertEquals(want[3], got[3], TOL, "load-dependent PH Tput");
    }

    @Test
    public void testOpenRegimeIsServed() {
        // the JAR refused this shape outright before; the level space is truncated
        // at options.cutoff, so the CTMC oracle gets the same bound
        Network model = new Network("ldqbd_open");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1", 0);
        source.setArrival(oclass, new Exp(0.9));
        queue.setService(oclass, new Erlang(2.0, 2));
        model.link(model.serialRouting(source, queue, sink));

        MAMOptions options = new MAMOptions();
        options.method = "ldqbd";
        options.cutoff = new Matrix(1, 1);
        options.cutoff.set(0, 0, 60);
        SolverMAM solver = new SolverMAM(model, options);
        double qlen = solver.getAvgQLen().get(QI, 0);
        double util = solver.getAvgUtil().get(QI, 0);
        double tput = solver.getAvgTput().get(QI, 0);

        // lambda = 0.9, mean service 1.0 over two servers: rho = 0.45 per server,
        // and the truncation at 60 leaves nothing measurable in the tail
        assertEquals(0.9, tput, 1e-9, "open Tput is the offered rate");
        assertEquals(0.45, util, 1e-9, "open Util is rho per server");
        assertTrue(qlen > 0.9 && qlen < 1.3, "open QLen out of range: " + qlen);
    }
}
