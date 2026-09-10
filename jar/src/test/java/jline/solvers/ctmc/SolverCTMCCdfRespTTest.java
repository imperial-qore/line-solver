/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ctmc;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

import jline.io.Ret;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.util.matrix.Matrix;

/**
 * Cross-codebase parity of the CTMC response-time distribution.
 *
 * <p>THE GOLDEN IS MATLAB AND C++, WHICH AGREE TO TEN DECIMALS. On the model
 * below (Delay Exp(2) + FCFS Queue Exp(1), one closed class of 2 jobs) both
 * report, on the same 0.999-truncated grid,
 *
 * <pre>
 *     Think : n = 3455 points, F(end) = 0.9990002447, E[T] = 0.4995002888
 *     Q     : n = 8840 points, F(end) = 0.9990007864, E[T] = 1.6655708356
 * </pre>
 *
 * <p>E[T] is the trapezoidal integral of the survival function over the grid,
 * so it carries the truncation and is NOT the exact mean response time; that is
 * deliberate, since it is what makes the three codebases comparable point for
 * point rather than only in the limit.
 */
public class SolverCTMCCdfRespTTest {

    private static Network buildModel() {
        Network model = new Network("cqn2");
        Delay think = new Delay(model, "Think");
        Queue q = new Queue(model, "Q", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", 2, think);
        think.setService(c, new Exp(2.0));
        q.setService(c, new Exp(1.0));
        model.link(model.serialRouting(think, q));
        return model;
    }

    /** Trapezoidal integral of the survival function: E[T] over the grid. */
    private static double meanFromCdf(Matrix curve) {
        double acc = 0.0;
        for (int i = 1; i < curve.getNumRows(); i++) {
            double t0 = curve.get(i - 1, 1);
            double t1 = curve.get(i, 1);
            double s0 = 1.0 - curve.get(i - 1, 0);
            double s1 = 1.0 - curve.get(i, 0);
            acc += 0.5 * (t1 - t0) * (s0 + s1);
        }
        return acc;
    }

    @Test
    public void getCdfRespTReturnsCurvesNotASummary() {
        Network model = buildModel();
        SolverCTMC solver = new SolverCTMC(model);

        Ret.DistributionResult d = solver.getCdfRespT(solver.getAvgRespTHandles());
        assertNotNull(d);
        assertNotNull(d.cdfData);
        assertEquals(2, d.numStations);

        // A DISTRIBUTION, not one number per pair: the defect this override
        // fixes was keeping only the last point of each curve.
        for (int ist = 0; ist < 2; ist++) {
            Matrix curve = d.cdfData.get(ist).get(0);
            assertNotNull(curve, "station " + ist + " has no curve");
            assertTrue(curve.getNumRows() > 100,
                    "station " + ist + " returned " + curve.getNumRows() + " points, not a curve");
            assertEquals(2, curve.getNumCols());
            // Column order is [F(t) t], as MATLAB and C++ store it.
            assertTrue(curve.get(0, 1) <= curve.get(curve.getNumRows() - 1, 1),
                    "column 1 must be the time and must increase");
            for (int i = 1; i < curve.getNumRows(); i++) {
                assertTrue(curve.get(i, 0) >= curve.get(i - 1, 0) - 1e-12,
                        "the CDF must be monotone");
            }
        }
    }

    @Test
    public void theCurvesMatchTheMatlabAndCppGolden() {
        Network model = buildModel();
        SolverCTMC solver = new SolverCTMC(model);
        Ret.DistributionResult d = solver.getCdfRespT(solver.getAvgRespTHandles());

        Matrix think = d.cdfData.get(0).get(0);
        Matrix queue = d.cdfData.get(1).get(0);

        assertEquals(0.9990002447, think.get(think.getNumRows() - 1, 0), 1e-6);
        assertEquals(0.9990007864, queue.get(queue.getNumRows() - 1, 0), 1e-6);
        assertEquals(0.4995002888, meanFromCdf(think), 1e-4);
        assertEquals(1.6655708356, meanFromCdf(queue), 1e-3);
    }

    @Test
    public void getCdfSysRespTIsTheCycleTimeLawPerChain() {
        // The system response time is the CYCLE TIME: one full trip of a tagged
        // job round the network, from one arrival at its reference station to
        // the next. For this model that is Z + R = 0.5 + 5/3 = 13/6.
        //
        // C++ reports, on the same 10000-interval grid truncated at 1-FineTol:
        //     sys chain 0 : n = 2183, F(end) = 0.9999999901, E[T] = 2.1666666563
        Network model = buildModel();
        SolverCTMC solver = new SolverCTMC(model);

        List<Matrix> rd = solver.getCdfSysRespT();
        assertEquals(1, rd.size(), "one law per chain");
        Matrix cycle = rd.get(0);
        assertNotNull(cycle);
        assertTrue(cycle.getNumRows() > 100,
                "a distribution, not a scalar: got " + cycle.getNumRows() + " points");
        assertEquals(2, cycle.getNumCols());
        for (int i = 1; i < cycle.getNumRows(); i++) {
            assertTrue(cycle.get(i, 0) >= cycle.get(i - 1, 0) - 1e-12, "the CDF must be monotone");
            assertTrue(cycle.get(i, 1) > cycle.get(i - 1, 1), "the grid must increase");
        }

        assertEquals(0.9999999901, cycle.get(cycle.getNumRows() - 1, 0), 1e-7);
        assertEquals(2.1666666563, meanFromCdf(cycle), 1e-4);
        // and it is the exact cycle time 13/6, since the truncation is at
        // 1 - FineTol rather than the coarser per-station one
        assertEquals(13.0 / 6.0, meanFromCdf(cycle), 1e-4);
    }
}
