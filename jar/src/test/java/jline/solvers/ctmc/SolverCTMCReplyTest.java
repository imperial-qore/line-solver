package jline.solvers.ctmc;

import java.util.List;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.ClosedSignal;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Synchronous calls (REPLY signals) under SolverCTMC.
 *
 * <p>A class whose jobs make a synchronous call (sn.syncreply >= 0) leaves the caller
 * for the callee but KEEPS its server there; the server is released only when the
 * matching REPLY signal class arrives back. LDES keys that hold by job id
 * (Solver_ssj.pendingReplyMap); a CTMC has no job identity, so the held servers are
 * counted per class in a local-state block ({@link jline.lang.state.ReplyBlock}) and the
 * caller is selected structurally, as the station the REPLY class is routed into.
 *
 * <p>Model (the canonical LQN shape RefTask -&gt; ClientTask --synchCall--&gt; ServerTask):
 *
 * <pre>
 *   ThinkDelay -&gt; ClientQueue -&gt; ServerQueue -&gt; [switch to Reply]
 *        ^             ^                              |
 *        |             |______ Reply signal __________|
 *        |____________ [switch back to Job] ___________|
 * </pre>
 *
 * <p>Conventions pinned here, all taken from LDES:
 * <ul>
 *   <li>QLen and Util at the caller COUNT the held server and the job blocked out at the
 *       callee (simultaneous resource possession, the point of the feature); RespT does
 *       NOT, since it measures time spent at the station.</li>
 *   <li>The REPLY does not queue at the caller: it takes the server it released and is
 *       routed on, so its queue length there is zero.</li>
 *   <li>Queue lengths therefore sum to N PLUS the mean number of outstanding calls (a
 *       blocked job is counted at the callee and at the caller alike); the token count is
 *       still N, since the reply IS the calling job.</li>
 * </ul>
 *
 * <p>Goldens are the MATLAB CTMC values, themselves validated against the LDES sample
 * path at 3e5 samples in line-test.git/test/testsAdvFeatures/des/test_ctmc_reply.m.
 */
public class SolverCTMCReplyTest {

    private static final double TOL = 1e-5;

    /** N, muClient, muServer, thinkTime, nservers(client). */
    private static final double[][] CFG = {
            {2, 4, 6, 1.0, 1},
            {3, 4, 6, 1.0, 1},
            {2, 2, 3, 0.5, 1},
            {4, 5, 5, 1.0, 1},
            {3, 4, 6, 1.0, 2},
    };

    /** MATLAB CTMC goldens per configuration: X, clientQLen, clientUtil, clientRespT, serverQLen. */
    private static final double[][] GOLDEN = {
            {1.31661, 0.68339, 0.54859, 0.35238, 0.21944},
            {1.79528, 1.20472, 0.74803, 0.50439, 0.29921},
            {1.10092, 1.44954, 0.91743, 0.98333, 0.36697},
            {2.17923, 1.82077, 0.87169, 0.63551, 0.43585},
            {2.03818, 0.96182, 0.45949, 0.27102, 0.40944},
    };

    private static Network replyModel(int N, double muC, double muS, double Z, int c) {
        Network model = new Network("ClosedReplySignal");
        Delay delay = new Delay(model, "ThinkDelay");
        Queue queue1 = new Queue(model, "ClientQueue", SchedStrategy.FCFS);
        queue1.setNumberOfServers(c);
        Queue queue2 = new Queue(model, "ServerQueue", SchedStrategy.FCFS);

        ClosedClass jobClass = new ClosedClass(model, "Job", N, delay);
        delay.setService(jobClass, new Exp(1 / Z));
        queue1.setService(jobClass, new Exp(muC));
        queue2.setService(jobClass, new Exp(muS));

        // The reply is the calling job, class-switched at the callee, so the population
        // is conserved; ClosedSignal shares the caller's reference station.
        ClosedSignal replyClass = new ClosedSignal(model, "Reply", SignalType.REPLY, delay);
        replyClass.forJobClass(jobClass);
        delay.setService(replyClass, new Immediate());
        queue1.setService(replyClass, new Immediate());
        queue2.setService(replyClass, new Immediate());

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, delay, queue1, 1.0);
        P.set(jobClass, jobClass, queue1, queue2, 1.0);
        P.set(jobClass, replyClass, queue2, queue1, 1.0); // callee replies
        P.set(replyClass, jobClass, queue1, delay, 1.0);  // caller unblocks, job resumes
        model.link(P);
        return model;
    }

    private static int row(NetworkAvgTable t, String station, String jobClass) {
        List<String> stations = t.getStationNames();
        List<String> classes = t.getClassNames();
        for (int i = 0; i < stations.size(); i++) {
            if (stations.get(i).equals(station) && classes.get(i).equals(jobClass)) {
                return i;
            }
        }
        throw new IllegalArgumentException("no row for " + station + "/" + jobClass);
    }

    private static NetworkAvgTable solve(int N, double muC, double muS, double Z, int c) {
        Network model = replyModel(N, muC, muS, Z, c);
        SolverCTMC solver = new SolverCTMC(model);
        solver.options.verbose = VerboseLevel.SILENT;
        return solver.getAvgTable();
    }

    /**
     * Reference configuration N=2, muC=4, muS=6, Z=1, c=1 against the MATLAB CTMC
     * goldens, to five decimals.
     */
    @Test
    public void testReplyBaseConfiguration() {
        NetworkAvgTable t = solve(2, 4, 6, 1.0, 1);
        int rDelay = row(t, "ThinkDelay", "Job");
        int rClient = row(t, "ClientQueue", "Job");
        int rClientReply = row(t, "ClientQueue", "Reply");
        int rServer = row(t, "ServerQueue", "Job");

        assertEquals(1.31661, t.getTput().get(rDelay), TOL);
        assertEquals(1.31661, t.getTput().get(rClient), TOL);
        assertEquals(1.31661, t.getTput().get(rServer), TOL);
        assertEquals(0.68339, t.getQLen().get(rClient), TOL);
        assertEquals(0.54859, t.getUtil().get(rClient), TOL);
        assertEquals(0.35238, t.getRespT().get(rClient), TOL);
        assertEquals(0.0, t.getQLen().get(rClientReply), 1e-8);
        assertEquals(0.21944, t.getQLen().get(rServer), TOL);
    }

    /**
     * Structural invariants that hold for every configuration: uniform throughput
     * around the cycle, no reply resident at the caller, the queue-length excess over
     * N equal to the number of outstanding calls, and a utilization strictly above the
     * carried load because the held server is counted.
     */
    @Test
    public void testReplyInvariants() {
        for (int i = 0; i < CFG.length; i++) {
            int N = (int) CFG[i][0];
            double muC = CFG[i][1];
            double muS = CFG[i][2];
            double Z = CFG[i][3];
            int c = (int) CFG[i][4];
            NetworkAvgTable t = solve(N, muC, muS, Z, c);
            int rDelay = row(t, "ThinkDelay", "Job");
            int rClient = row(t, "ClientQueue", "Job");
            int rClientReply = row(t, "ClientQueue", "Reply");
            int rServer = row(t, "ServerQueue", "Job");
            String msg = "config " + i;

            // MATLAB CTMC goldens, all five configurations, X / clientQLen / clientUtil
            // / clientRespT / serverQLen.
            assertEquals(GOLDEN[i][0], t.getTput().get(rDelay), TOL, msg + ": Tput");
            assertEquals(GOLDEN[i][1], t.getQLen().get(rClient), TOL, msg + ": client QLen");
            assertEquals(GOLDEN[i][2], t.getUtil().get(rClient), TOL, msg + ": client Util");
            assertEquals(GOLDEN[i][3], t.getRespT().get(rClient), TOL, msg + ": client RespT");
            assertEquals(GOLDEN[i][4], t.getQLen().get(rServer), TOL, msg + ": server QLen");

            double X = t.getTput().get(rDelay);
            assertTrue(X > 0, msg + ": zero throughput");
            // Throughput is uniform around the cycle (one token per call).
            assertEquals(X, t.getTput().get(rClient), 1e-6, msg);
            assertEquals(X, t.getTput().get(rServer), 1e-6, msg);

            // The reply never resides at the caller: it takes the server it released
            // and is routed on.
            assertEquals(0.0, t.getQLen().get(rClientReply), 1e-8, msg);

            // Simultaneous resource possession double-counts by construction: a job
            // blocked out at the callee is counted BOTH there and at the caller that
            // holds a server for it. So the station queue lengths sum to N plus the
            // mean number of outstanding calls, which here is exactly the callee's
            // queue length (every call in progress sits at the ServerQueue). The
            // underlying token count is still N -- the reply IS the calling job.
            double sumQ = 0;
            for (int r = 0; r < t.getQLen().size(); r++) {
                sumQ += t.getQLen().get(r);
            }
            assertEquals(t.getQLen().get(rServer), sumQ - N, 1e-6, msg + ": QLen excess");

            // Utilization counts the held server: it exceeds the carried load T*E[S].
            assertTrue(t.getUtil().get(rClient) > t.getTput().get(rClient) / muC / c + 0.05,
                    msg + ": client Util collapsed onto the carried load");

            // Response time excludes the time blocked out at the callee, so it is
            // strictly below QLen/Tput whenever calls are outstanding.
            assertTrue(t.getRespT().get(rClient) < t.getQLen().get(rClient) / t.getTput().get(rClient),
                    msg + ": RespT still counts the blocked-out job");
        }
    }

    /**
     * Holding a server across a call has no representation in the state of the
     * non-FCFS disciplines, so a caller declared PS must be REJECTED rather than
     * silently solved with the hold dropped -- which would report the carried
     * load and understate the caller by the whole call duration.
     */
    @Test
    public void nonFcfsCallerIsRejected() {
        Network model = new Network("ClosedReplyPS");
        Delay delay = new Delay(model, "ThinkDelay");
        Queue queue1 = new Queue(model, "ClientQueue", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "ServerQueue", SchedStrategy.FCFS);
        ClosedClass jobClass = new ClosedClass(model, "Job", 2, delay);
        delay.setService(jobClass, new Exp(1.0));
        queue1.setService(jobClass, new Exp(4.0));
        queue2.setService(jobClass, new Exp(6.0));
        ClosedSignal replyClass = (ClosedSignal) new ClosedSignal(model, "Reply", SignalType.REPLY, delay)
                .forJobClass(jobClass);
        delay.setService(replyClass, new Immediate());
        queue1.setService(replyClass, new Immediate());
        queue2.setService(replyClass, new Immediate());
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, delay, queue1, 1.0);
        P.set(jobClass, jobClass, queue1, queue2, 1.0);
        P.set(jobClass, replyClass, queue2, queue1, 1.0);
        P.set(replyClass, jobClass, queue1, delay, 1.0);
        model.link(P);

        boolean rejected = false;
        try {
            new SolverCTMC(model, "verbose", VerboseLevel.SILENT).getAvgTable();
        } catch (RuntimeException e) {
            rejected = true;
        }
        assertTrue(rejected, "a PS caller was accepted; the held server cannot be represented there");
    }
}
