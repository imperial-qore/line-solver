/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Regression test: the dec.source fixed point must conserve the closed
 * population at a multiserver station.
 *
 * <p>The post-loop second pass rescaled QN to the chain population, then floored
 * the response time at one full service time and restated QN = RN*TN from an
 * untouched TN, discarding the rescale. With c &gt; 1 the floor roughly doubled
 * RN, so this model (10 servers, N = 2) returned sum(QN) = 3.52 instead of 2.
 * Reference is the exact CTMC: with 10 servers and 2 jobs no job ever queues,
 * so R = S and MAM must reproduce CTMC outright.</p>
 */
public class MamClosedMultiserverPopulationTest {

    private static final double[][] RATES = {{96.0, 60.0}, {70.0, 84.0}, {96.0, 1.0}};

    private static Network build(int nservers) {
        Network model = new Network("cqn");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        ClosedClass classA = new ClosedClass(model, "ClassA", 1, delay);
        ClosedClass classB = new ClosedClass(model, "ClassB", 1, delay);
        Station[] st = {delay, queue1, queue2};
        JobClass[] jc = {classA, classB};
        for (int i = 0; i < 3; i++) {
            for (int r = 0; r < 2; r++) {
                st[i].setService(jc[r], new Exp(RATES[i][r]));
            }
            if (i > 0) {
                ((Queue) st[i]).setNumberOfServers(nservers);
            }
        }
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(classA, classA, delay, queue1, 0.6);
        P.set(classA, classA, delay, queue2, 0.4);
        P.set(classA, classA, queue1, delay, 1.0);
        P.set(classA, classA, queue2, delay, 1.0);
        P.set(classB, classB, delay, queue1, 1.0);
        P.set(classB, classB, queue1, delay, 1.0);
        P.set(classB, classB, queue2, delay, 1.0);
        model.link(P);
        return model;
    }

    private static double columnSum(Matrix m, int col) {
        double s = 0.0;
        for (int i = 0; i < m.getNumRows(); i++) {
            s += m.get(i, col);
        }
        return s;
    }

    @Test
    public void testClosedPopulationIsConserved() {
        for (int nservers : new int[]{1, 2, 10}) {
            SolverMAM solver = new SolverMAM(build(nservers), "verbose", VerboseLevel.SILENT);
            solver.getAvg();
            // one job per class, and each class is its own chain
            assertEquals(1.0, columnSum(solver.result.QN, 0), 1e-6,
                    "class A population, c=" + nservers);
            assertEquals(1.0, columnSum(solver.result.QN, 1), 1e-6,
                    "class B population, c=" + nservers);
        }
    }

    @Test
    public void testMultiserverMatchesCtmcExactly() {
        SolverMAM mam = new SolverMAM(build(10), "verbose", VerboseLevel.SILENT);
        mam.getAvg();
        SolverCTMC ctmc = new SolverCTMC(build(10), "verbose", VerboseLevel.SILENT);
        ctmc.getAvg();
        for (int i = 0; i < mam.result.QN.getNumRows(); i++) {
            for (int r = 0; r < mam.result.QN.getNumCols(); r++) {
                assertEquals(ctmc.result.QN.get(i, r), mam.result.QN.get(i, r), 1e-9);
                assertEquals(ctmc.result.RN.get(i, r), mam.result.RN.get(i, r), 1e-9);
                assertEquals(ctmc.result.TN.get(i, r), mam.result.TN.get(i, r), 1e-6);
            }
        }
    }
}
