package jline.solvers.ag;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/**
 * RCAT with phase-type processes: each component is a QBD over (queue length,
 * phase) rather than a scalar birth-death chain.
 *
 * WHAT MAKES THESE ORACLES. An isolated M/PH/1 has no synchronizing action, so
 * whatever the reversed-rate iterate does, its marginal is the M/G/1 one and its
 * mean is the Pollaczek-Khinchine value rho + rho^2 (1+scv) / (2 (1-rho)),
 * written out here from the arrival rate and the SCV alone. The same holds for
 * the first station of a tandem. So these assert the phase construction itself,
 * not the fixed point -- and they would ALL have failed before it existed,
 * because the analyzer then returned the M/M/1 answer for every SCV.
 *
 * The downstream station of the tandem has no closed form (the departure stream
 * of an M/PH/1 is not Poisson), so what is asserted there is flow conservation,
 * which every RCAT estimator must satisfy and which the phase expansion must not
 * break.
 */
public class RcatPhaseTypeTest {

    private static final double TOL = 1e-6;

    /** M/G/1 mean number in system. */
    private static double pk(double rho, double scv) {
        return rho + rho * rho * (1.0 + scv) / (2.0 * (1.0 - rho));
    }

    private Network single(Distribution arr, Distribution svc) {
        Network model = new Network("mg1");
        Source source = new Source(model, "S");
        Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "K");
        OpenClass oc = new OpenClass(model, "C", 0);
        source.setArrival(oc, arr);
        queue.setService(oc, svc);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    private Network tandem(Distribution arr, Distribution s1, Distribution s2) {
        Network model = new Network("t");
        Source source = new Source(model, "S");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "K");
        OpenClass oc = new OpenClass(model, "C", 0);
        source.setArrival(oc, arr);
        q1.setService(oc, s1);
        q2.setService(oc, s2);
        model.link(Network.serialRouting(source, q1, q2, sink));
        return model;
    }

    /** Erlang(k) service of mean 1, i.e. SCV = 1/k. */
    private static Erlang erlangMeanOne(int k) {
        return new Erlang((double) k, k);
    }

    @Test
    public void testMPh1MarginalIsTheExactPollaczekKhinchineMean() {
        double[] rhos = {0.5, 0.8};
        for (double rho : rhos) {
            for (int k : new int[]{2, 4}) {
                Network model = single(new Exp(rho), erlangMeanOne(k));
                for (String m : new String[]{"inap", "inapinf"}) {
                    SolverAG s = new SolverAG(model, m);
                    double q = s.getAvgQLen().get(1, 0);
                    assertEquals(pk(rho, 1.0 / k), q, TOL,
                        m + ": M/Er" + k + "/1 at rho=" + rho + " must be the exact P-K mean");
                    assertEquals(rho, s.getAvgUtil().get(1, 0), TOL, m + ": utilization is rho");
                    assertEquals(rho, s.getAvgTput().get(1, 0), TOL, m + ": throughput is lambda");
                }
            }
        }
    }

    /**
     * The exponential answer must be reachable through the phase machinery: an
     * Erlang(1) is an exponential, and a component with one phase per level is
     * the scalar birth-death chain the analyzer built before.
     */
    @Test
    public void testExponentialIsUnchangedByThePhaseConstruction() {
        double rho = 0.5;
        Network erl1 = single(new Exp(rho), erlangMeanOne(1));
        Network exp = single(new Exp(rho), new Exp(1.0));
        for (String m : new String[]{"inap", "inapplus", "inapinf"}) {
            assertEquals(new SolverAG(exp, m).getAvgQLen().get(1, 0),
                    new SolverAG(erl1, m).getAvgQLen().get(1, 0), 1e-9,
                    m + ": Erlang(1) service must reproduce Exp service");
        }
    }

    /** A non-Poisson Source needs the arrival-phase dimension of the component. */
    @Test
    public void testPhArrivalsAreNotAnsweredAsPoisson() {
        // Er2/M/1 with mean interarrival 2 and mean service 1: the M/M/1 answer
        // at the same rho would be 1.0, so a Poisson reading is visible here.
        Network model = single(new Erlang(1.0, 2), new Exp(1.0));
        for (String m : new String[]{"inap", "inapinf"}) {
            double q = new SolverAG(model, m).getAvgQLen().get(1, 0);
            assertTrue(q < 0.95, m + ": Er2/M/1 must be below the M/M/1 value 1.0, got " + q);
            assertTrue(q > 0.6, m + ": Er2/M/1 must be a sane queue length, got " + q);
        }
    }

    /**
     * Flow conservation through a tandem whose first station is phase-type. The
     * first component is an isolated M/PH/1 whatever the fixed point does, so its
     * mean is the P-K value; the reversed rate must then be the arrival rate, so
     * both stations must carry it.
     */
    @Test
    public void testTandemWithPhServiceConservesFlow() {
        double lambda = 0.5;
        Network model = tandem(new Exp(lambda), erlangMeanOne(2), new Exp(2.0));
        for (String m : new String[]{"inap", "inapplus", "inapinf"}) {
            SolverAG s = new SolverAG(model, m);
            Matrix q = s.getAvgQLen();
            Matrix t = s.getAvgTput();
            assertEquals(pk(lambda, 0.5), q.get(1, 0), TOL,
                m + ": the upstream station is the exact M/Er2/1");
            assertEquals(lambda, t.get(1, 0), TOL, m + ": upstream throughput is lambda");
            assertEquals(lambda, t.get(2, 0), TOL, m + ": downstream throughput is lambda");
        }
    }

    /**
     * 'inapinf' drops the truncation, so on a heavy-tailed service law it lands
     * on the exact P-K mean where the 100-level truncation of 'inap' cannot.
     */
    @Test
    public void testInapinfRemovesTheTruncationOnAPhaseComponent() {
        double lambda = 0.8;
        // HyperExp of mean 1 and SCV 4: the geometric tail is slow enough that
        // 100 levels lose more than a percent of the first moment.
        HyperExp svc = HyperExp.fitMeanAndSCV(1.0, 4.0);
        Network model = tandem(new Exp(lambda), svc, new Exp(4.0));
        double exact = pk(lambda, 4.0);
        double inap = new SolverAG(model, "inap").getAvgQLen().get(1, 0);
        double inapinf = new SolverAG(model, "inapinf").getAvgQLen().get(1, 0);
        assertEquals(exact, inapinf, 1e-6,
            "inapinf must be the untruncated M/PH/1 mean, got " + inapinf);
        assertTrue(Math.abs(inap - exact) > Math.abs(inapinf - exact),
            "the truncated inap must be further from the exact mean than inapinf");
    }
}
