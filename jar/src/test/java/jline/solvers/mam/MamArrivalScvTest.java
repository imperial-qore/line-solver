/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.APH;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Regression test: the dec.source fixed point must preserve the arrival
 * variability of open chains (a defect previously overwrote the source
 * MMAP with a Poisson process each iteration, making PH/PH/1 results
 * arrival-scv-independent). References are the exact CTMC values.
 */
public class MamArrivalScvTest {

    private static Network build(double ca) {
        Network model = new Network("g");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "C");
        source.setArrival(oclass, APH.fitMeanAndSCV(2.0, ca));
        queue.setService(oclass, APH.fitMeanAndSCV(1.0, 4.0));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    @Test
    public void testPhPh1ArrivalScvDependence() {
        SolverMAM s2 = new SolverMAM(build(2.0), "verbose", VerboseLevel.SILENT);
        s2.getAvg();
        double L2 = s2.result.QN.get(1, 0);
        SolverMAM s16 = new SolverMAM(build(16.0), "verbose", VerboseLevel.SILENT);
        s16.getAvg();
        double L16 = s16.result.QN.get(1, 0);
        assertEquals(2.08423, L2, 5e-2);  // CTMC exact
        assertEquals(4.62730, L16, 2e-1); // CTMC exact
    }
}
