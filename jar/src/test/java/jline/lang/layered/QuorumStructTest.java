/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.io.File;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The quorum count of an AND-join must survive LQNX parsing, the struct conversion and the
 * LQNX round-trip.
 *
 * <p>Before this was fixed, getStruct folded the quorum into the activity graph edge weight as
 * k/n, which both lost k and made an AND-join edge indistinguishable from a routing probability.
 * The fixtures are the 88-quorum and 89-quorum models of the LQNS regression corpus, plus
 * variants with k = 1 and k = n.</p>
 */
public class QuorumStructTest {

    private static final String DIR = "src/test/resources/lqn/quorum/";

    private static LayeredNetwork read(String name) {
        File f = new File(DIR + name);
        assertTrue(f.exists(), "missing fixture " + f.getPath());
        return LayeredNetwork.parseXML(f.getPath(), false);
    }

    /** Index of the single activity that is an AND-join target, or -1. */
    private static int joinTarget(LayeredNetworkStruct lsn) {
        for (int i = 1; i < lsn.actquorum.getNumCols(); i++) {
            if (lsn.actquorum.get(0, i) > 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * A 2-of-3 quorum must reach the struct as k = 2 on the join target.
     */
    @Test
    public void testQuorumReachesStruct() {
        LayeredNetworkStruct lsn = read("test_LQN_quorum_2of3.lqnx").getStruct();
        int aidx = joinTarget(lsn);
        assertTrue(aidx > 0, "no join target carried a quorum");
        assertEquals(2.0, lsn.actquorum.get(0, aidx), 1e-12);
    }

    /**
     * The struct value must track the model, so the three variants must differ.
     */
    @Test
    public void testQuorumVariantsAreDistinct() {
        double[] expected = {1.0, 2.0, 3.0};
        String[] files = {"test_LQN_quorum_1of3.lqnx", "test_LQN_quorum_2of3.lqnx",
                          "test_LQN_quorum_3of3.lqnx"};
        for (int i = 0; i < files.length; i++) {
            LayeredNetworkStruct lsn = read(files[i]).getStruct();
            int aidx = joinTarget(lsn);
            assertTrue(aidx > 0, files[i] + ": no join target carried a quorum");
            assertEquals(expected[i], lsn.actquorum.get(0, aidx), 1e-12, files[i]);
        }
    }

    /**
     * AND-join edges must carry weight 1, not the k/n fraction that was previously written
     * there, since the graph weight is read elsewhere as a routing probability.
     */
    @Test
    public void testAndJoinEdgeWeightsAreUnity() {
        LayeredNetworkStruct lsn = read("test_LQN_quorum_2of3.lqnx").getStruct();
        int aidx = joinTarget(lsn);
        assertTrue(aidx > 0);
        int branches = 0;
        for (int p = 1; p < lsn.graph.getNumRows(); p++) {
            double w = lsn.graph.get(p, aidx);
            if (w != 0) {
                assertEquals(1.0, w, 1e-12, "edge " + p + " -> " + aidx + " is not unity");
                branches++;
            }
        }
        assertEquals(3, branches, "expected three branches into the join");
    }

    /**
     * getQuorumCount must treat the several "no quorum" encodings alike and clamp nonsense.
     */
    @Test
    public void testGetQuorumCountSentinels() {
        assertEquals(3, ActivityPrecedence.getQuorumCount(null, 3));
        assertEquals(3, ActivityPrecedence.getQuorumCount(new Matrix(0, 0), 3));
        // A row vector of ones is not a quorum count; it was the old "all required" sentinel.
        assertEquals(3, ActivityPrecedence.getQuorumCount(Matrix.ones(1, 3), 3));
        assertEquals(2, ActivityPrecedence.getQuorumCount(Matrix.singleton(2), 3));
        // Out of range values fall back to "wait for all".
        assertEquals(3, ActivityPrecedence.getQuorumCount(Matrix.singleton(0), 3));
        assertEquals(3, ActivityPrecedence.getQuorumCount(Matrix.singleton(9), 3));
    }

    /**
     * A plain AND-join built from Activity objects must not acquire a spurious quorum, which
     * previously happened because the factory seeded preParams with ones(1, n) and the writer
     * then read its first element as k = 1.
     */
    @Test
    public void testPlainAndJoinHasNoQuorum() {
        LayeredNetwork model = new LayeredNetwork("plain");
        Processor p = new Processor(model, "p1", 1, jline.lang.constant.SchedStrategy.PS);
        Task t = new Task(model, "t1", 1, jline.lang.constant.SchedStrategy.REF);
        t.on(p);
        Entry e = new Entry(model, "e1");
        e.on(t);
        Activity a0 = new Activity(model, "a0", new jline.lang.processes.Exp(1.0));
        Activity b1 = new Activity(model, "b1", new jline.lang.processes.Exp(1.0));
        Activity b2 = new Activity(model, "b2", new jline.lang.processes.Exp(1.0));
        Activity c1 = new Activity(model, "c1", new jline.lang.processes.Exp(1.0));
        a0.on(t);
        b1.on(t);
        b2.on(t);
        c1.on(t);
        // A reference task drives the workload and never replies, so a0 is simply bound to e.
        a0.boundTo(e);
        t.addPrecedence(ActivityPrecedence.AndFork(a0, java.util.Arrays.asList(b1, b2)));
        t.addPrecedence(ActivityPrecedence.AndJoin(java.util.Arrays.asList(b1, b2), c1));

        LayeredNetworkStruct lsn = model.getStruct();
        int aidx = joinTarget(lsn);
        assertTrue(aidx > 0, "join target should record its branch count");
        // Two branches, no quorum: k equals n, so the join waits for all of them.
        assertEquals(2.0, lsn.actquorum.get(0, aidx), 1e-12);
        assertNotNull(lsn.actquorum);
    }
}
