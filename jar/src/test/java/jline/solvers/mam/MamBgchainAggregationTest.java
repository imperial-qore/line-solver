/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The bgchain method's aggregation control, {@code config.bgaggr}, and the
 * per-class station support it stresses.
 *
 * <p>The model has three closed chains with DIFFERENT ROUTES -- one confined to
 * Queue1, one to Queue2, one crossing both -- plus an open class. That is the
 * shape where aggregating a group really distorts, because the aggregate carries
 * a flow-weighted MEAN of its members' routing matrices, and it is also the
 * shape that requires each background class to be enumerated over ITS OWN
 * stations: over the union, a chain confined to Queue1 would be given
 * configurations at Queue2 that it cannot reach and cannot leave (a zero routing
 * row row-normalizes to a self-loop), so the generator turns reducible and half
 * the chain's mass leaks into states it never visits.</p>
 *
 * <p>The oracle is SolverCTMC at cutoff 8 on the same model: Delay OnQ1 0.53043,
 * Delay Both 0.48639.</p>
 */
public class MamBgchainAggregationTest {

    private static final double TOL = 1e-6;

    private static Network build() {
        Network model = new Network("bgRoutes");
        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Q2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        OpenClass open = new OpenClass(model, "Open");
        ClosedClass onQ1 = new ClosedClass(model, "OnQ1", 1, delay);
        ClosedClass onQ2 = new ClosedClass(model, "OnQ2", 1, delay);
        ClosedClass both = new ClosedClass(model, "Both", 1, delay);

        source.setArrival(open, new Exp(0.4));
        queue1.setService(open, new Exp(4.0));
        queue2.setService(open, new Exp(4.0));
        delay.setService(onQ1, new Exp(1.0));
        queue1.setService(onQ1, new Exp(1.5));
        delay.setService(onQ2, new Exp(1.0));
        queue2.setService(onQ2, new Exp(1.5));
        delay.setService(both, new Exp(1.0));
        queue1.setService(both, new Exp(3.0));
        queue2.setService(both, new Exp(3.0));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(open, open, source, queue1, 1.0);
        P.set(open, open, queue1, queue2, 1.0);
        P.set(open, open, queue2, sink, 1.0);
        P.set(onQ1, onQ1, delay, queue1, 1.0);
        P.set(onQ1, onQ1, queue1, delay, 1.0);
        P.set(onQ2, onQ2, delay, queue2, 1.0);
        P.set(onQ2, onQ2, queue2, delay, 1.0);
        P.set(both, both, delay, queue1, 1.0);
        P.set(both, both, queue1, queue2, 1.0);
        P.set(both, both, queue2, delay, 1.0);
        model.link(P);
        return model;
    }

    private static Matrix solve(int bgaggr) {
        Network model = build();
        SolverOptions options = SolverMAM.defaultOptions();
        options.method = "bgchain";
        options.verbose = VerboseLevel.SILENT;
        options.config.put("bgaggr", Integer.valueOf(bgaggr));
        SolverMAM solver = new SolverMAM(model, options);
        return solver.getAvgQLen();
    }

    /** Stations are Source(0), Delay(1), Q1(2), Q2(3); classes Open(0), OnQ1(1), OnQ2(2), Both(3). */
    @Test
    public void testEachClosedChainConservesItsPopulation() {
        for (int g = 1; g <= 3; g++) {
            Matrix QN = solve(g);
            for (int c = 1; c <= 3; c++) {
                double held = QN.get(1, c) + QN.get(2, c) + QN.get(3, c);
                assertEquals(1.0, held, TOL,
                        "closed chain " + c + " must hold its single job at bgaggr=" + g);
            }
            // A chain confined to Q1 holds nothing at Q2, and conversely.
            assertEquals(0.0, QN.get(3, 1), TOL, "OnQ1 must hold nothing at Q2");
            assertEquals(0.0, QN.get(2, 2), TOL, "OnQ2 must hold nothing at Q1");
        }
    }

    @Test
    public void testRaisingBgaggrSeparatesMisrepresentedChains() {
        final double refDelayBoth = 0.48639;   // SolverCTMC, cutoff 8
        Matrix q1 = solve(1);
        Matrix q2 = solve(2);
        double err1 = Math.abs(q1.get(1, 3) - refDelayBoth) / refDelayBoth;
        double err2 = Math.abs(q2.get(1, 3) - refDelayBoth) / refDelayBoth;
        // One aggregate cannot represent two different routes; two carry them exactly.
        assertTrue(err1 > 5e-3, "bgaggr=1 should visibly misrepresent the routes, got " + err1);
        assertTrue(err2 < 5e-4, "bgaggr=2 should carry every chain exactly, got " + err2);
        assertTrue(err2 < err1, "raising bgaggr must not make the answer worse");
    }

    @Test
    public void testBgaggrClampsAboveTheChainCount() {
        // R = 3, so bgaggr = 3 is already "no aggregation" and must answer as 2
        // does: asking for no aggregation needs no magic value.
        Matrix q2 = solve(2);
        Matrix q3 = solve(3);
        for (int i = 0; i < 4; i++) {
            for (int c = 0; c < 4; c++) {
                assertEquals(q2.get(i, c), q3.get(i, c), 1e-12,
                        "bgaggr above R-1 must clamp, station " + i + " class " + c);
            }
        }
    }
}
