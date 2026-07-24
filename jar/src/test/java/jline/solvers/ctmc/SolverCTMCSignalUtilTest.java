package jline.solvers.ctmc;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.Signal;
import jline.lang.constant.RemovalPolicy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.DiscreteDistribution;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.Geometric;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Utilization convention on G-networks solved by SolverCTMC.
 *
 * Utilization at a station exposed to G-network signals is the CARRIED load,
 * i.e. the mean busy-server fraction. With exponential service the completion
 * rate is mu times the mean number of busy servers, so E[busy]/c = T*E[S]/c
 * holds exactly even though signals remove jobs without a completion; the
 * arrival-based estimator lambda*E[S]/c instead measures the OFFERED load and
 * is invalid there. {@link jline.solvers.ctmc.handlers.CtmcSignalLossy} marks
 * the classes a signal can annihilate at a station so that only the
 * departure-based estimator is used for them.
 *
 * Cases:
 * <ol>
 *   <li>single-server FCFS, untargeted negative customers (exact Gelenbe)</li>
 *   <li>single-server PS, untargeted negative customers (exact, same rho)</li>
 *   <li>single-server FCFS, catastrophes (exact, quadratic root)</li>
 *   <li>PS, two positive classes, signal targeting one (victim class only)</li>
 *   <li>multiserver c=2 (MATLAB CTMC golden, LDES-validated)</li>
 *   <li>Geometric batch removal (MATLAB CTMC golden, LDES-validated)</li>
 *   <li>Erlang-2 service (busy fraction strictly above T*E[S])</li>
 * </ol>
 *
 * Cases 5 to 7 have no closed form; their goldens come from the MATLAB CTMC,
 * itself cross-checked against the LDES sample path in
 * line-test.git/test/testsAdvFeatures/des/test_gnetwork_ctmc_util.m.
 */
public class SolverCTMCSignalUtilTest {

    private static final double MU = 1.0;
    private static final double TOL = 1e-6;

    /** Single queue with one positive class and one signal class. */
    private static Network gnetwork1(SchedStrategy sched, int nservers, double lambdaPos,
                                     double lambdaNeg, SignalType signalType,
                                     DiscreteDistribution remDist, RemovalPolicy remPolicy) {
        Network model = new Network("GNetworkUtil");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", sched);
        queue.setNumberOfServers(nservers);
        Sink sink = new Sink(model, "Sink");

        OpenClass pos = new OpenClass(model, "Positive");
        source.setArrival(pos, new Exp(lambdaPos));
        queue.setService(pos, new Exp(MU));

        Signal neg = (remDist == null)
                ? new Signal(model, "Negative", signalType)
                : new Signal(model, "Negative", signalType, 0, remDist, remPolicy);
        source.setArrival(neg, new Exp(lambdaNeg));
        queue.setService(neg, new Exp(MU));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(pos, pos, source, queue, 1.0);
        P.set(pos, pos, queue, sink, 1.0);
        P.set(neg, neg, source, queue, 1.0);
        P.set(neg, neg, queue, sink, 1.0);
        model.link(P);
        return model;
    }

    /** Row index of (station, class) in the flattened avg table. */
    private static int row(NetworkAvgTable t, String station, String jobClass) {
        List<String> stations = t.getStationNames();
        List<String> classes = t.getClassNames();
        for (int i = 0; i < stations.size(); i++) {
            if (stations.get(i).equals(station) && classes.get(i).equals(jobClass)) {
                return i;
            }
        }
        throw new IllegalStateException("row " + station + "/" + jobClass + " missing");
    }

    private static NetworkAvgTable solve(Network model, int cutoff) {
        return new SolverCTMC(model, "cutoff", cutoff, "verbose", VerboseLevel.SILENT).getAvgTable();
    }

    /**
     * rho = lambda+/(mu + lambda-) by flow balance: every arrival either
     * completes or is removed while the server is busy.
     */
    @Test
    public void fcfsNegativeCustomersMatchGelenbe() {
        double lambdaPos = 0.5, lambdaNeg = 0.4;
        NetworkAvgTable t = solve(gnetwork1(SchedStrategy.FCFS, 1, lambdaPos, lambdaNeg,
                SignalType.NEGATIVE, null, null), 30);
        int i = row(t, "Queue", "Positive");
        double rho = lambdaPos / (MU + lambdaNeg);
        assertEquals(rho, t.getUtil().get(i), TOL);
        assertEquals(rho / (1.0 - rho), t.getQLen().get(i), TOL);
        assertEquals(MU * rho, t.getTput().get(i), TOL);
        // The offered load lambda+/mu = 0.5 is what the arrival estimator gave.
        assertTrue(Math.abs(t.getUtil().get(i) - lambdaPos / MU) > 0.1,
                "Util fell back to the offered load");
    }

    /**
     * The Gelenbe product form is discipline-invariant here, so PS must return
     * the same rho and mean queue length as FCFS.
     */
    @Test
    public void psNegativeCustomersMatchGelenbe() {
        double lambdaPos = 0.5, lambdaNeg = 0.4;
        NetworkAvgTable t = solve(gnetwork1(SchedStrategy.PS, 1, lambdaPos, lambdaNeg,
                SignalType.NEGATIVE, null, null), 30);
        int i = row(t, "Queue", "Positive");
        double rho = lambdaPos / (MU + lambdaNeg);
        assertEquals(rho, t.getUtil().get(i), TOL);
        assertEquals(rho / (1.0 - rho), t.getQLen().get(i), TOL);
    }

    /**
     * Catastrophes flush the station: p_n = (1-r) r^n with
     * mu r^2 - (lambda + mu + delta) r + lambda = 0, so Util = P(busy) = r.
     */
    @Test
    public void catastropheMatchesQuadraticRoot() {
        double lambdaPos = 0.5, delta = 0.4;
        NetworkAvgTable t = solve(gnetwork1(SchedStrategy.FCFS, 1, lambdaPos, delta,
                SignalType.CATASTROPHE, null, null), 40);
        int i = row(t, "Queue", "Positive");
        double b = lambdaPos + MU + delta;
        double r = (b - Math.sqrt(b * b - 4 * lambdaPos * MU)) / (2 * MU);
        assertEquals(r, t.getUtil().get(i), TOL);
        assertEquals(r / (1.0 - r), t.getQLen().get(i), TOL);
        assertEquals(MU * r, t.getTput().get(i), TOL);
    }

    /**
     * A targeted signal downgrades ONLY its victim class: C1 is never removed,
     * so all its arrivals complete and its utilization stays at lambda1/mu.
     */
    @Test
    public void targetedSignalDowngradesOnlyItsVictim() {
        double lambda1 = 0.2, lambda2 = 0.2, lambdaNeg = 0.5;
        Network model = new Network("GNetworkTargeted");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass c1 = new OpenClass(model, "C1");
        OpenClass c2 = new OpenClass(model, "C2");
        source.setArrival(c1, new Exp(lambda1));
        queue.setService(c1, new Exp(MU));
        source.setArrival(c2, new Exp(lambda2));
        queue.setService(c2, new Exp(MU));
        Signal neg = new Signal(model, "Negative", SignalType.NEGATIVE).forJobClass(c2);
        source.setArrival(neg, new Exp(lambdaNeg));
        queue.setService(neg, new Exp(MU));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, source, queue, 1.0);
        P.set(c1, c1, queue, sink, 1.0);
        P.set(c2, c2, source, queue, 1.0);
        P.set(c2, c2, queue, sink, 1.0);
        P.set(neg, neg, source, queue, 1.0);
        P.set(neg, neg, queue, sink, 1.0);
        model.link(P);

        NetworkAvgTable t = solve(model, 12);
        int i1 = row(t, "Queue", "C1");
        int i2 = row(t, "Queue", "C2");
        assertEquals(lambda1 / MU, t.getUtil().get(i1), 1e-4);
        assertTrue(t.getUtil().get(i2) < lambda2 / MU - 0.01,
                "victim Util not below the offered load");
        assertEquals(t.getTput().get(i2) / MU, t.getUtil().get(i2), 1e-4);
        assertEquals(0.12482827, t.getUtil().get(i2), 1e-4);
    }

    /** No closed form; golden is the MATLAB CTMC value (LDES-validated). */
    @Test
    public void multiserverUtilIsBusyServerFraction() {
        double lambdaPos = 1.2, lambdaNeg = 0.3;
        int c = 2;
        NetworkAvgTable t = solve(gnetwork1(SchedStrategy.FCFS, c, lambdaPos, lambdaNeg,
                SignalType.NEGATIVE, null, null), 25);
        int i = row(t, "Queue", "Positive");
        assertEquals(t.getTput().get(i) / (c * MU), t.getUtil().get(i), 1e-8);
        assertTrue(t.getUtil().get(i) < lambdaPos / (c * MU) - 0.05,
                "multiserver Util not below the offered load");
        assertEquals(0.50119329, t.getUtil().get(i), 1e-4);
    }

    /**
     * Phase-type service. A job destroyed mid-service leaves busy time behind
     * with no completion, so T*E[S] (0.34941) under-counts; the exact
     * busy-server occupancy is 0.37648, which is what the LDES sample path
     * measures (0.37696 at 5e5 samples).
     *
     * This case also pins the ToMarginal phase-slicing repair: the event layer
     * threads the full (nstations x nclasses) sn.phasessz through, which was
     * indexed as a single station row, so a class-0 job in phase 2 was counted
     * as a class-1 job and the signal handler saw no eligible victim (the JAR
     * generator was missing the removal transition out of that state entirely,
     * giving T = 0.39003 against MATLAB's 0.34941).
     */
    @Test
    public void phaseTypeServiceBusyFraction() {
        double lambdaPos = 0.5, lambdaNeg = 0.4;
        Network model = new Network("GNetworkPH");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass pos = new OpenClass(model, "Positive");
        source.setArrival(pos, new Exp(lambdaPos));
        queue.setService(pos, Erlang.fitMeanAndOrder(1.0 / MU, 2));
        Signal neg = new Signal(model, "Negative", SignalType.NEGATIVE);
        source.setArrival(neg, new Exp(lambdaNeg));
        queue.setService(neg, new Exp(MU));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(pos, pos, source, queue, 1.0);
        P.set(pos, pos, queue, sink, 1.0);
        P.set(neg, neg, source, queue, 1.0);
        P.set(neg, neg, queue, sink, 1.0);
        model.link(P);

        NetworkAvgTable t = solve(model, 20);
        int i = row(t, "Queue", "Positive");
        assertEquals(0.34940753, t.getTput().get(i), TOL);
        assertEquals(0.55459104, t.getQLen().get(i), TOL);
        assertEquals(0.37648116, t.getUtil().get(i), TOL);
        // Strictly above the carried load: partial service is not a completion.
        assertTrue(t.getUtil().get(i) > t.getTput().get(i) + 0.02,
                "Util collapsed onto T*E[S]");
        // Flow balance: arrivals = completions + removals = T + lambda- * Util.
        assertEquals(lambdaPos, t.getTput().get(i) + lambdaNeg * t.getUtil().get(i), TOL);
    }

    /**
     * Geometric batch size clipped at the population; golden is the MATLAB CTMC
     * value (LDES-validated).
     */
    @Test
    public void batchRemovalUtilIsCarriedLoad() {
        double lambdaPos = 0.8, lambdaNeg = 0.3;
        NetworkAvgTable t = solve(gnetwork1(SchedStrategy.FCFS, 1, lambdaPos, lambdaNeg,
                SignalType.NEGATIVE, new Geometric(0.5), RemovalPolicy.RANDOM), 35);
        int i = row(t, "Queue", "Positive");
        assertEquals(t.getTput().get(i) / MU, t.getUtil().get(i), 1e-8);
        assertTrue(t.getUtil().get(i) < lambdaPos / MU - 0.05,
                "batch Util not below the offered load");
        assertEquals(0.56421833, t.getUtil().get(i), 1e-4);
    }
}
