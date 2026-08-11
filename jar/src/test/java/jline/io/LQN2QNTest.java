/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.ActivityPrecedence;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.processes.Exp;
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
 * <p>Both the AND-join quorum and the activity think time used to be read into the
 * LayeredNetworkStruct and then silently dropped by the conversion: a k-of-n join became a
 * wait-for-all Join, and a think time vanished instead of appearing as a processor-releasing
 * delay in series with the host demand.</p>
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
}
