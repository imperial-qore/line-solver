/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import org.junit.jupiter.api.Test;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Finite-capacity single-server FCFS with two classes: the exact MMAP[K]/G/1/K
 * branch must fire and give each class its OWN loss ratio.
 *
 * <p>The truncate-and-renormalize fallback returns a single blocking
 * probability and sets T_k = lambda_k (1 - p), which makes the loss ratio
 * identical across classes by construction. References are MATLAB, which has
 * carried this branch since qsys_mmapg1k landed.
 */
public class MamFiniteCapPerClassLossTest {

    private static Network build() {
        Network model = new Network("mamcap");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass c1 = new OpenClass(model, "C1");
        OpenClass c2 = new OpenClass(model, "C2");
        source.setArrival(c1, new Exp(0.4));
        source.setArrival(c2, Erlang.fitMeanAndSCV(1.0 / 0.5, 0.25));
        queue.setService(c1, new Exp(2.0));
        queue.setService(c2, HyperExp.fitMeanAndSCV(0.4, 3.0));
        queue.setNumberOfServers(1);
        queue.setCapacity(6);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(source, queue, sink));
        P.set(c2, c2, Network.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    @Test
    public void perClassLossMatchesMatlab() {
        SolverMAM s = new SolverMAM(build(), "verbose", VerboseLevel.SILENT);
        s.getAvg();
        // MATLAB SolverMAM, default method, on the same model. Rebased 2026-08-16,
        // when mmap_super_safe stopped letting its SCV sort decide the mark order:
        // the Exp chain (SCV 1) and the Erlang one (SCV 0.25) had their arrival
        // streams paired with each other's service laws. The oracle is the
        // relabelling invariant -- the PRE-FIX solver on the model with the two
        // classes exchanged already returned these numbers.
        assertEquals(0.30366872864559902, s.result.QN.get(1, 0), 1e-6);
        assertEquals(0.33037104539618389, s.result.QN.get(1, 1), 1e-6);
        assertEquals(0.39671466500670860, s.result.TN.get(1, 0), 1e-8);
        assertEquals(0.49645605237493995, s.result.TN.get(1, 1), 1e-8);
        assertEquals(0.76545879301049757, s.result.RN.get(1, 0), 1e-6);
        assertEquals(0.66545879301049760, s.result.RN.get(1, 1), 1e-6);

        // The two classes must NOT share a loss ratio: that is what the exact
        // branch buys over truncate-and-renormalize.
        double loss1 = 1.0 - s.result.TN.get(1, 0) / 0.4;
        double loss2 = 1.0 - s.result.TN.get(1, 1) / 0.5;
        assertTrue(Math.abs(loss1 - loss2) > 1e-4,
                "per-class loss ratios collapsed: " + loss1 + " vs " + loss2);
    }
}
