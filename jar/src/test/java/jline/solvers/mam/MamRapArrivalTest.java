package jline.solvers.mam;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Distribution;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.lang.processes.ME;
import jline.lang.processes.RAP;
import jline.solvers.NetworkAvgTable;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * SolverMAM must not answer a correlated arrival with a closed form that reads
 * only the arrival's renewal marginal.
 *
 * <p>qsys_phmc solves PH/M/c from the pair (pie, D0), i.e. from the interarrival
 * distribution alone. It used to be reached whenever the source process was
 * merely non-exponential, which let a correlated MAP, RAP or ME arrival onto it
 * and discarded the autocorrelation: for the RAP below the answer came out at
 * 1.98523360058856, close to what Kingman gives from the arrival SCV alone,
 * against a true value near 6.844. The gate is now RENEWAL, so a correlated
 * arrival falls through to MMAPPH1FCFS, which is handed the arrival (D0, D1)
 * and therefore carries its correlation. A renewal arrival, phase-type or
 * matrix-exponential, still takes the fast path.</p>
 */
public class MamRapArrivalTest {

    private static final double TOL = 1e-8;

    /**
     * Rate-1 RAP with a negative off-diagonal entry in H0, so it is a genuine
     * rational arrival process and not a MAP. Mean 1, SCV 4.5448028674,
     * lag-1 autocorrelation 0.3432018450.
     */
    private static RAP correlatedRap() {
        double s = 30.0 / 119.0;
        Matrix H0 = new Matrix(new double[][]{{-9.9 * s, -0.2 * s}, {0.1 * s, -1.0 * s}});
        Matrix H1 = new Matrix(new double[][]{{9.7 * s, 0.4 * s}, {0.0, 0.9 * s}});
        return new RAP(H0, H1);
    }

    /**
     * The non-phase-type ME used across the ME anchors: alpha carries a negative
     * entry and the density has an interior zero, so no phase-type
     * representation of any order exists. It is nonetheless a RENEWAL process.
     */
    private static ME nonPhaseTypeME() {
        double w = 2 * Math.PI;
        Matrix alpha = new Matrix(new double[]{
                0.984694494294579, -0.040430911430916, 0.0557364171363366});
        Matrix A = new Matrix(new double[][]{
                {-0.5, 0, 0},
                {0, -1, w},
                {0, -w, -1}
        });
        return new ME(alpha, A);
    }

    private static double queueLength(Distribution arrival, Distribution service, int nservers) {
        Network model = new Network("mam_rap_arrival");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, arrival);
        queue.setService(oclass, service);
        queue.setNumberOfServers(nservers);
        model.link(model.serialRouting(source, queue, sink));
        NetworkAvgTable table = new SolverMAM(model).getAvgTable();
        return table.getQLen().get(1);
    }

    /**
     * The correlated RAP arrival. 6.84449125861694 is MMAPPH1FCFS applied to the
     * true (H0, H1); an independent SolverLDES run of 60 replications of 4e6
     * samples brackets it at [6.83540045, 6.86599427], so the analytic value is
     * pinned against simulation and not merely against itself.
     */
    @Test
    public void testCorrelatedRapArrivalCarriesAutocorrelation() {
        double qlen = queueLength(correlatedRap(), new Exp(2.0), 1);
        assertEquals(6.84449125861694, qlen, TOL);
        assertTrue(qlen > 6.83540045 && qlen < 6.86599427,
                "QLen " + qlen + " outside the SolverLDES interval [6.83540045, 6.86599427]");
    }

    /**
     * The discrimination that makes this class of bug visible: two arrival
     * processes with the SAME mean and SAME SCV but different lag-1
     * autocorrelation must give materially different queue lengths. A solver
     * that reads only the marginal cannot tell them apart and returns nearly the
     * same number for both, which is exactly the failure mode this test exists
     * to catch.
     */
    @Test
    public void testMomentIdenticalArrivalsWithDifferentAcfDiffer() {
        RAP correlated = correlatedRap();
        double mean = correlated.getMean();
        double scv = correlated.getSCV();
        // Renewal companion matched on the first two moments: acf1 = 0 by
        // construction, since a hyperexponential is a renewal process.
        HyperExp renewal = HyperExp.fitMeanAndSCV(mean, scv);
        assertEquals(mean, renewal.getMean(), 1e-9);
        assertEquals(scv, renewal.getSCV(), 1e-9);

        double qCorrelated = queueLength(correlated, new Exp(2.0), 1);
        double qRenewal = queueLength(renewal, new Exp(2.0), 1);

        assertTrue(qCorrelated > 3.0 * qRenewal,
                "moment-identical arrivals differing only in autocorrelation gave "
                        + qCorrelated + " and " + qRenewal
                        + "; the correlation is not reaching the queue-length computation");
    }

    /**
     * A matrix-exponential arrival is RENEWAL, so it must keep taking the
     * marginal-only fast path and must not move. Pinned to ten digits.
     */
    @Test
    public void testRenewalMEArrivalUnchangedOnFastPath() {
        assertEquals(1.0370455525, queueLength(nonPhaseTypeME(), new Exp(1 / 0.977420), 1), 1e-9);
    }

    /**
     * Renewal phase-type arrivals, single and multi-server, keep using
     * qsys_phmc. PH/M/2 at 1.1867726850 is inside the SolverLDES interval
     * [1.18614125, 1.18750740], so the fast path retained here is the accurate
     * one.
     */
    @Test
    public void testRenewalPhArrivalsUnchangedOnFastPath() {
        assertEquals(0.8090169944, queueLength(Erlang.fitMeanAndOrder(1.0, 2), new Exp(2.0), 1), 1e-9);
        assertEquals(1.1867726850, queueLength(Erlang.fitMeanAndOrder(0.5, 2), new Exp(2.0), 2), 1e-9);
        assertEquals(1.0, queueLength(new Exp(1.0), new Exp(2.0), 1), 1e-9);
        assertEquals(3.4285714286, queueLength(new Exp(3.0), new Exp(2.0), 2), 1e-9);
    }

    /**
     * ME SERVICE anchors, which this dispatch change must not touch: the gate
     * added here is on the ARRIVAL process only.
     */
    @Test
    public void testMEServiceAnchorsUnchanged() {
        assertEquals(0.8749999987, queueLength(new Exp(0.5), ME.fromErlang(2, 2.0), 1), 1e-9);
        assertEquals(1.5222222174, queueLength(new Exp(0.5),
                ME.fromHyperExp(new double[]{0.6, 0.4}, new double[]{2.0, 0.5}), 1), 1e-9);
        assertEquals(1.0152146751, queueLength(new Exp(0.255775446238906), nonPhaseTypeME(), 1), 1e-9);
    }

    /**
     * Phase-type SERVICE anchors, single and multi-server: unchanged.
     */
    @Test
    public void testPhServiceAnchorsUnchanged() {
        assertEquals(0.3125000000, queueLength(new Exp(0.5), Erlang.fitMeanAndOrder(0.5, 2), 1), 1e-9);
        // M/E2/2. Was 0.6964285714, the single-fast-server surrogate; the exact
        // MAP/PH/c multiset QBD (2026-08-16) answers it at 0.6462542940, and a
        // truncated CTMC written out state by state, sharing no code with that QBD,
        // gives 0.6462542936 at N=300 with 2.9e-16 of mass left in the tail. The
        // old value was 7.8% high.
        assertEquals(0.6462542940, queueLength(new Exp(1.2), Erlang.fitMeanAndOrder(0.5, 2), 2), 1e-9);
    }
}
