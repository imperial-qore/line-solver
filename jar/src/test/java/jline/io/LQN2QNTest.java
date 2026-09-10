/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.ActivityPrecedence;
import jline.lang.layered.CacheTask;
import jline.lang.layered.Entry;
import jline.lang.layered.SetupTask;
import jline.lang.layered.ItemEntry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkStruct;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.lang.processes.Distribution;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.processes.DiscreteSampler;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.sections.Joiner;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * LQN features that LQN2QN must carry into the converted queueing network.
 *
 * <p>The AND-join quorum, the activity think time and the replication factor used to be read
 * into the LayeredNetworkStruct and then silently dropped by the conversion: a k-of-n join
 * became a wait-for-all Join, a think time vanished instead of appearing as a
 * processor-releasing delay in series with the host demand, and r replicas collapsed into a
 * single unscaled station.</p>
 */
public class LQN2QNTest {

    private static Matrix quorumOf(int k) {
        Matrix m = new Matrix(1, 1, 1);
        m.set(0, 0, k);
        return m;
    }

    /**
     * Client -> Server, the server entry forking into three branches joined with the given
     * quorum, optionally with a think time on the first branch.
     */
    private static LayeredNetwork forkJoinModel(int quorum, boolean withThink) {
        LayeredNetwork m = new LayeredNetwork("fj");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.INF);
        Processor p2 = new Processor(m, "P2", 1, SchedStrategy.PS);
        Task client = new Task(m, "Client", 5, SchedStrategy.REF).on(p1);
        client.setThinkTime(new Exp(1.0));
        Entry ce = new Entry(m, "CE").on(client);
        Activity ca = new Activity(m, "ca", new Exp(10.0)).on(client).boundTo(ce);
        ca.synchCall("SE", 1.0);

        Task server = new Task(m, "Server", 3, SchedStrategy.FCFS).on(p2);
        Entry se = new Entry(m, "SE").on(server);
        Activity a0 = new Activity(m, "a0", new Exp(1.0)).on(server).boundTo(se);
        List<Activity> branches = new ArrayList<Activity>();
        for (int i = 1; i <= 3; i++) {
            branches.add(new Activity(m, "b" + i, new Exp(2.0)).on(server));
        }
        Activity aj = new Activity(m, "aj", new Exp(3.0)).on(server);
        if (withThink) {
            branches.get(0).setThinkTime(new Exp(4.0));
        }
        server.addPrecedence(ActivityPrecedence.AndFork(a0, branches));
        server.addPrecedence(ActivityPrecedence.AndJoin(branches, aj, quorumOf(quorum)));
        aj.repliesTo(se);
        return m;
    }

    private static Join theJoin(Network qn) {
        for (Node n : qn.getNodes()) {
            if (n instanceof Join) {
                return (Join) n;
            }
        }
        return null;
    }

    private static Delay actThinkStation(Network qn) {
        for (Node n : qn.getNodes()) {
            if (n instanceof Delay && "ActivityThink".equals(n.getName())) {
                return (Delay) n;
            }
        }
        return null;
    }

    /** A genuine quorum k &lt; n reaches the Join node as a Quorum strategy requiring k. */
    @Test
    public void quorumJoinIsCarriedToTheJoinNode() {
        Network qn = LQN2QN.convert(forkJoinModel(2, false));
        Join join = theJoin(qn);
        assertNotNull(join, "the AND-fork must produce a Join node");
        Joiner joiner = (Joiner) join.getInput();
        int nset = 0;
        for (JobClass c : qn.getClasses()) {
            Double required = joiner.joinRequired.get(c);
            if (required != null && required > 0) {
                assertEquals(2.0, required, 1e-12);
                assertEquals(JoinStrategy.Quorum, joiner.joinStrategy.get(c));
                nset++;
            }
        }
        assertEquals(1, nset, "the quorum belongs to the single class that enters the fork");
    }

    /** A quorum equal to the branch count is an ordinary wait-for-all join, left as STD. */
    @Test
    public void fullQuorumLeavesTheJoinStandard() {
        Network qn = LQN2QN.convert(forkJoinModel(3, false));
        Joiner joiner = (Joiner) theJoin(qn).getInput();
        for (JobClass c : qn.getClasses()) {
            Double required = joiner.joinRequired.get(c);
            assertTrue(required == null || required <= 0,
                    "no quorum must be set when k equals the branch count");
            JoinStrategy js = joiner.joinStrategy.get(c);
            assertTrue(js == null || js == JoinStrategy.STD);
        }
    }

    /** An activity think time becomes a step on a shared INF station, so the processor is released. */
    @Test
    public void activityThinkTimeBecomesADelayStep() {
        Network qn = LQN2QN.convert(forkJoinModel(2, true));
        Delay think = actThinkStation(qn);
        assertNotNull(think, "an activity think time must produce the ActivityThink delay");
        int nserved = 0;
        for (JobClass c : qn.getClasses()) {
            if (think.getServiceProcess(c) != null && think.getServiceProcess(c).getMean() > 1e-8) {
                assertEquals(0.25, think.getServiceProcess(c).getMean(), 1e-12);
                nserved++;
            }
        }
        assertEquals(1, nserved, "exactly the one think-bearing activity has a think step");
    }

    /** Without any activity think time the shared delay is not created at all. */
    @Test
    public void noThinkTimeCreatesNoDelayStation() {
        assertNull(actThinkStation(LQN2QN.convert(forkJoinModel(2, false))));
    }

    // ------------------------------------------------------------- replication

    /**
     * Client -&gt; Server, both replicated rep-fold with one task replica per processor replica.
     * A positive fan-out makes each client replica call every server replica instead of its own.
     */
    private static LayeredNetwork replicatedModel(int rep, int fanOut, int servers, int taskMult) {
        LayeredNetwork m = new LayeredNetwork("repl");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.INF);
        Processor p2 = new Processor(m, "P2", servers, SchedStrategy.FCFS);
        p1.setReplication(rep);
        p2.setReplication(rep);
        Task client = new Task(m, "Client", 3, SchedStrategy.REF).on(p1);
        client.setReplication(rep);
        client.setThinkTime(new Exp(1.0));
        Entry ce = new Entry(m, "CE").on(client);
        Activity ca = new Activity(m, "ca", new Exp(10.0)).on(client).boundTo(ce);
        ca.synchCall("SE", 1.0);
        Task server = new Task(m, "Server", taskMult, SchedStrategy.FCFS).on(p2);
        server.setReplication(rep);
        if (fanOut > 0) {
            client.setFanOut("Server", fanOut);
        }
        Entry se = new Entry(m, "SE").on(server);
        Activity a0 = new Activity(m, "a0", new Exp(2.0)).on(server).boundTo(se);
        a0.repliesTo(se);
        return m;
    }

    private static Node nodeNamed(Network qn, String name) {
        for (Node n : qn.getNodes()) {
            if (name.equals(n.getName())) {
                return n;
            }
        }
        return null;
    }

    private static int classesServedAt(Network qn, Station st) {
        int n = 0;
        for (JobClass c : qn.getClasses()) {
            if (st.getServiceProcess(c) != null && st.getServiceProcess(c).getMean() > 1e-8) {
                n++;
            }
        }
        return n;
    }

    private static List<String> classNames(Network qn) {
        List<String> names = new ArrayList<String>();
        for (JobClass c : qn.getClasses()) {
            names.add(c.getName());
        }
        return names;
    }

    /** Three replicas give three processor stations and three closed chains. */
    @Test
    public void materialisedReplicasAreSeparateStationsAndChains() {
        Network qn = LQN2QN.convert(replicatedModel(3, 0, 1, 2), "materialize");
        assertNotNull(nodeNamed(qn, "P2"));
        assertNotNull(nodeNamed(qn, "P2_r2"));
        assertNotNull(nodeNamed(qn, "P2_r3"));
        int nthink = 0;
        for (JobClass c : qn.getClasses()) {
            if (c.getName().startsWith("Client_Think")) {
                assertEquals(3.0, ((ClosedClass) c).getPopulation(), 1e-12);
                nthink++;
            }
        }
        assertEquals(3, nthink, "one closed chain per reference task replica");
    }

    /** At fan-out 1 replica i calls replica i, so the chains stay disjoint. */
    @Test
    public void pairingReplicatesTheSubsystemWithoutCouplingIt() {
        Network qn = LQN2QN.convert(replicatedModel(2, 1, 1, 2), "materialize");
        assertEquals(1, classesServedAt(qn, (Station) nodeNamed(qn, "P2")));
        assertEquals(1, classesServedAt(qn, (Station) nodeNamed(qn, "P2_r2")));
    }

    /** At fan-out r each caller replica reaches all r callee replicas. */
    @Test
    public void broadcastFanOutCouplesEveryCallerReplicaToEveryReplica() {
        Network qn = LQN2QN.convert(replicatedModel(2, 2, 1, 2), "materialize");
        assertEquals(2, classesServedAt(qn, (Station) nodeNamed(qn, "P2")));
        assertEquals(2, classesServedAt(qn, (Station) nodeNamed(qn, "P2_r2")));
    }

    /** Pooling keeps one station of r times the capacity and one r-fold chain. */
    @Test
    public void pooledReplicasScaleServersPopulationAndAdmission() {
        Network qn = LQN2QN.convert(replicatedModel(3, 0, 2, 2), "pool");
        assertNull(nodeNamed(qn, "P2_r2"), "the replicas collapse into one station");
        assertEquals(6, ((Queue) nodeNamed(qn, "P2")).getNumberOfServers());
        int nthink = 0;
        for (JobClass c : qn.getClasses()) {
            if (c.getName().startsWith("Client_Think")) {
                assertEquals(9.0, ((ClosedClass) c).getPopulation(), 1e-12);
                nthink++;
            }
        }
        assertEquals(1, nthink);
        assertEquals(6.0, qn.getRegions().get(0).getLinearConstraints()[1].get(0), 1e-12);
    }

    /** Each task replica owns its own threads, so it owns its own admission row. */
    @Test
    public void materialisedThreadPoolsAreOneAdmissionRowPerReplica() {
        Network qn = LQN2QN.convert(replicatedModel(3, 0, 2, 2), "materialize");
        Matrix[] Ab = qn.getRegions().get(0).getLinearConstraints();
        Matrix A = Ab[0];
        Matrix b = Ab[1];
        assertEquals(3, A.getNumRows());
        for (int r = 0; r < 3; r++) {
            assertEquals(2.0, b.get(r), 1e-12);
        }
        for (int c = 0; c < A.getNumCols(); c++) {
            double col = 0.0;
            for (int r = 0; r < A.getNumRows(); r++) {
                col += A.get(r, c);
            }
            assertTrue(col <= 1.0, "no class is admitted by two replica rows");
        }
    }

    // --------------------------------------------------------------- retrieval

    /**
     * Client -&gt; CacheTask, the read forking into a hit and a miss activity, the miss being
     * the fetch. With retrieval the fetch becomes a retrieval system.
     */
    private static LayeredNetwork cacheModel(boolean retrieval) {
        LayeredNetwork m = new LayeredNetwork("lcq");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.INF);
        Task t1 = new Task(m, "T1", 4, SchedStrategy.REF).on(p1);
        t1.setThinkTime(new Exp(1.0));
        Entry e1 = new Entry(m, "E1").on(t1);
        Processor pc = new Processor(m, "PC", 1, SchedStrategy.PS);
        CacheTask c2 = new CacheTask(m, "C2", 4, 2, ReplacementStrategy.RR, 1);
        c2.on(pc);
        if (retrieval) {
            c2.setRetrieval(true);
        }
        Matrix pAccess = new Matrix(1, 4);
        for (int i = 0; i < 4; i++) {
            pAccess.set(0, i, 0.25);
        }
        ItemEntry i2 = new ItemEntry(m, "I2", 4, new DiscreteSampler(pAccess)).on(c2);
        new Activity(m, "A1", new Immediate()).on(t1).boundTo(e1).synchCall("I2", 1);
        Activity ac2 = new Activity(m, "AC2", new Immediate()).on(c2).boundTo(i2);
        Activity hit = new Activity(m, "AC2h", new Exp(1.0)).on(c2).repliesTo(i2);
        Activity miss = new Activity(m, "AC2m", new Exp(0.5)).on(c2).repliesTo(i2);
        c2.addPrecedence(ActivityPrecedence.CacheAccess(ac2, java.util.Arrays.asList(hit, miss)));
        return m;
    }

    private static JobClass classNamed(Network qn, String name) {
        for (JobClass c : qn.getClasses()) {
            if (name.equals(c.getName())) {
                return c;
            }
        }
        return null;
    }

    /** A delayed-hit CacheTask gets one PS fetch station holding the miss demand. */
    @Test
    public void retrievalCreatesAFetchStation() {
        Network qn = LQN2QN.convert(cacheModel(true));
        Node fetch = nodeNamed(qn, "C2_Cache_Fetch");
        assertNotNull(fetch, "the retrieval system must be a station of the QN");
        JobClass read = classNamed(qn, "AC2");
        assertEquals(2.0, ((Station) fetch).getServiceProcess(read).getMean(), 1e-12);
    }

    /** The fetch happens in the retrieval system, so it is not charged twice. */
    @Test
    public void retrievalMovesTheMissDemandOffTheProcessor() {
        Network qn = LQN2QN.convert(cacheModel(true));
        Station pc = (Station) nodeNamed(qn, "PC");
        assertEquals(1.0, pc.getServiceProcess(classNamed(qn, "AC2h")).getMean(), 1e-12);
        Distribution missSvc = pc.getServiceProcess(classNamed(qn, "AC2m"));
        assertTrue(missSvc == null || missSvc.getMean() < 1e-8,
                "the miss demand moved to the fetch station");
    }

    /** Without the retrieval flag the miss branch keeps its demand at the processor. */
    @Test
    public void noRetrievalCreatesNoFetchStation() {
        Network qn = LQN2QN.convert(cacheModel(false));
        assertNull(nodeNamed(qn, "C2_Cache_Fetch"));
        Station pc = (Station) nodeNamed(qn, "PC");
        assertEquals(2.0, pc.getServiceProcess(classNamed(qn, "AC2m")).getMean(), 1e-12);
    }

    // ----------------------------------------------------------- setup task

    /**
     * Client -&gt; SetupTask: with a setup and a delay-off time the function's server
     * shuts down when idle and pays a cold start on the next arrival.
     */
    private static LayeredNetwork functionModel(boolean setup, SchedStrategy hostSched) {
        LayeredNetwork m = new LayeredNetwork("faas");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.INF);
        Processor pf = new Processor(m, "PF", 1, hostSched);
        Task c = new Task(m, "Client", 2, SchedStrategy.REF).on(p1);
        c.setThinkTime(new Exp(0.5));
        Entry ce = new Entry(m, "CE").on(c);
        SetupTask f = new SetupTask(m, "F", 1, SchedStrategy.FCFS);
        f.on(pf);
        if (setup) {
            f.setSetupTime(new Exp(2.0));
            f.setDelayOffTime(new Exp(1.0));
        }
        Entry fe = new Entry(m, "FE").on(f);
        new Activity(m, "ca", new Exp(10.0)).on(c).boundTo(ce).synchCall("FE", 1.0);
        new Activity(m, "a0", new Exp(2.0)).on(f).boundTo(fe).repliesTo(fe);
        return m;
    }

    private static int armedClasses(Network qn, String stationName) {
        Queue st = (Queue) nodeNamed(qn, stationName);
        int n = 0;
        for (JobClass c : qn.getClasses()) {
            if (st.getSetupTime(c) != null) {
                assertEquals(0.5, st.getSetupTime(c).getMean(), 1e-12);
                assertEquals(1.0, st.getDelayOffTime(c).getMean(), 1e-12);
                n++;
            }
        }
        return n;
    }

    /** The setup/delay-off pair is set per step class of the setup task. */
    @Test
    public void functionTaskSetupReachesTheHostStation() {
        Network qn = LQN2QN.convert(functionModel(true, SchedStrategy.FCFS));
        assertEquals(1, armedClasses(qn, "PF"));
    }

    /** An ordinary task never arms the setup/delay-off pair. */
    @Test
    public void noSetupLeavesTheStationAlwaysOn() {
        Network qn = LQN2QN.convert(functionModel(false, SchedStrategy.FCFS));
        assertEquals(0, armedClasses(qn, "PF"));
    }

    /** An infinite-server processor never shuts down, so it never sets up. */
    @Test
    public void setupOnAnInfiniteServerIsDropped() {
        Network qn = LQN2QN.convert(functionModel(true, SchedStrategy.INF));
        assertEquals(0, armedClasses(qn, "PF"));
    }

    /** A model without replication converts identically in all three modes. */
    @Test
    public void noReplicationLeavesTheConversionUntouched() {
        List<String> auto = classNames(LQN2QN.convert(forkJoinModel(2, true), "auto"));
        List<String> mat = classNames(LQN2QN.convert(forkJoinModel(2, true), "materialize"));
        List<String> pool = classNames(LQN2QN.convert(forkJoinModel(2, true), "pool"));
        assertEquals(auto, mat);
        assertEquals(auto, pool);
    }

    // ------------------------------------------------- call count and think SCV

    /** Reference task with the given think time, calling a server callMean times. */
    private static LayeredNetwork momentModel(Distribution think, double callMean) {
        LayeredNetwork m = new LayeredNetwork("mom");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.INF);
        Processor p2 = new Processor(m, "P2", 1, SchedStrategy.PS);
        Task c = new Task(m, "C", 2, SchedStrategy.REF).on(p1);
        c.setThinkTime(think);
        Entry ce = new Entry(m, "CE").on(c);
        Task srv = new Task(m, "S", 1, SchedStrategy.FCFS).on(p2);
        Entry se = new Entry(m, "SE").on(srv);
        new Activity(m, "ca", new Exp(1.0)).on(c).boundTo(ce).synchCall("SE", callMean);
        new Activity(m, "a0", new Exp(2.0)).on(srv).boundTo(se).repliesTo(se);
        return m;
    }

    /** The think delay carries the task distribution, not its mean fitted to an Exp. */
    @Test
    public void thinkTimeKeepsItsScvThroughTheConversion() {
        Network qn = LQN2QN.convert(momentModel(Erlang.fitMeanAndSCV(2.0, 0.25), 1.0));
        Delay delay = (Delay) nodeNamed(qn, "C_Think");
        int seen = 0;
        for (JobClass c : qn.getClasses()) {
            Distribution d = delay.getService(c);
            if (d != null && d.getMean() > 1e-9) {
                assertEquals(2.0, d.getMean(), 1e-9);
                assertEquals(0.25, d.getSCV(), 1e-9);
                seen++;
            }
        }
        assertEquals(1, seen);
    }

    /**
     * A call that happens with probability p has mean p and SCV (1-p)/p, which
     * Geometric(1/p) cannot represent: its parameter would exceed 1.
     */
    @Test
    public void callCountBelowOneIsBernoulli() {
        LayeredNetworkStruct lsn = momentModel(new Exp(1.0), 0.4).getStruct();
        assertEquals(1, lsn.ncalls);
        Distribution d = lsn.callproc.get(0);   // call indices are 0-based
        assertNotNull(d);
        assertEquals("Bernoulli", d.getClass().getSimpleName());
        assertEquals(0.4, d.getMean(), 1e-9);
        assertEquals(0.6 / 0.4, d.getSCV(), 1e-9);
    }

    /** At or above one call the count is geometric with the declared mean. */
    @Test
    public void callCountAboveOneIsGeometric() {
        LayeredNetworkStruct lsn = momentModel(new Exp(1.0), 3.0).getStruct();
        assertEquals(1, lsn.ncalls);
        Distribution d = lsn.callproc.get(0);   // call indices are 0-based
        assertNotNull(d);
        assertEquals("Geometric", d.getClass().getSimpleName());
        assertEquals(3.0, d.getMean(), 1e-9);
        assertEquals(1.0 - 1.0 / 3.0, d.getSCV(), 1e-9);
    }

    /**
     * A call declared with mean 0 is not a call: the count is the zero-mean
     * placeholder, the same class MATLAB and Python store.
     */
    @Test
    public void callCountOfZeroIsADegeneratePlaceholder() {
        LayeredNetworkStruct lsn = momentModel(new Exp(1.0), 0.0).getStruct();
        assertEquals(1, lsn.ncalls);
        Distribution d = lsn.callproc.get(0);   // call indices are 0-based
        assertNotNull(d);
        assertEquals("Immediate", d.getClass().getSimpleName());
        assertEquals(0.0, d.getMean(), 1e-12);
    }
}
