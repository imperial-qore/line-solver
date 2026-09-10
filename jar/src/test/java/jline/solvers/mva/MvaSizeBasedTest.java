/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva;

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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * SolverMVA on the size-based M/G/1 disciplines (Wierman and Harchol-Balter,
 * SIGMETRICS 2003), against MATLAB SolverMVA on the same model.
 *
 * <p>Two open classes at a single Queue: C1 is Exp(1) at rate 0.3, C2 is
 * Erlang(mean 2, SCV 0.5) at rate 0.2, so rho = 0.7 and the classes differ in
 * both mean size and variability -- which is what the size-based orderings act
 * on. The utilizations and throughputs are discipline-independent and the
 * response times are not, which is the check that the right analyzer ran.
 */
public class MvaSizeBasedTest {

    private static Network build(SchedStrategy sched) {
        Network model = new Network("sb");
        Source src = new Source(model, "Source");
        Queue q = new Queue(model, "Q", sched);
        Sink snk = new Sink(model, "Sink");
        OpenClass c1 = new OpenClass(model, "C1");
        OpenClass c2 = new OpenClass(model, "C2");
        src.setArrival(c1, new Exp(0.3));
        src.setArrival(c2, new Exp(0.2));
        q.setService(c1, new Exp(1.0));
        q.setService(c2, Erlang.fitMeanAndSCV(2.0, 0.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(src, q, snk));
        P.set(c2, c2, Network.serialRouting(src, q, snk));
        model.link(P);
        return model;
    }

    private static void check(SchedStrategy sched, double q1, double q2, double r1, double r2) {
        SolverMVA s = new SolverMVA(build(sched), "verbose", VerboseLevel.SILENT);
        s.getAvg();
        assertEquals(q1, s.result.QN.get(1, 0), 1e-9, sched + " QLen C1");
        assertEquals(q2, s.result.QN.get(1, 1), 1e-9, sched + " QLen C2");
        assertEquals(r1, s.result.RN.get(1, 0), 1e-9, sched + " RespT C1");
        assertEquals(r2, s.result.RN.get(1, 1), 1e-9, sched + " RespT C2");
        // discipline-independent: rho_k = lambda_k / mu_k, and the flow is lossless
        assertEquals(0.3, s.result.UN.get(1, 0), 1e-12, sched + " Util C1");
        assertEquals(0.4, s.result.UN.get(1, 1), 1e-12, sched + " Util C2");
        assertEquals(0.3, s.result.TN.get(1, 0), 1e-12, sched + " Tput C1");
        assertEquals(0.2, s.result.TN.get(1, 1), 1e-12, sched + " Tput C2");
        // Little's law at the station
        assertEquals(s.result.QN.get(1, 0), s.result.TN.get(1, 0) * s.result.RN.get(1, 0), 1e-9);
        assertEquals(s.result.QN.get(1, 1), s.result.TN.get(1, 1) * s.result.RN.get(1, 1), 1e-9);
    }

    @Test
    public void srptMatchesMatlab() {
        check(SchedStrategy.SRPT, 0.4985567765823074, 0.8212546141301067,
                1.6618559219410247, 4.1062730706505333);
    }

    @Test
    public void psjfMatchesMatlab() {
        check(SchedStrategy.PSJF, 0.61224489795918369, 3.333333333333333,
                2.0408163265306123, 16.666666666666664);
    }

    @Test
    public void fbMatchesMatlab() {
        check(SchedStrategy.FB, 0.63587346812077949, 2.1712143660569709,
                2.1195782270692649, 10.856071830284854);
    }

    @Test
    public void lrptMatchesMatlab() {
        check(SchedStrategy.LRPT, 2.1333333333333333, 0.66666666666666674,
                7.1111111111111107, 3.3333333333333335);
    }

    @Test
    public void setfMatchesMatlab() {
        check(SchedStrategy.SETF, 1.2256856111717471, 2.8758520284392772,
                4.0856187039058236, 14.379260142196385);
    }

    @Test
    public void theDisciplinesAreNotInterchangeable() {
        // SRPT favours the short class and LRPT the long one: if the analyzer
        // were size-blind the two would agree.
        SolverMVA srpt = new SolverMVA(build(SchedStrategy.SRPT), "verbose", VerboseLevel.SILENT);
        srpt.getAvg();
        SolverMVA lrpt = new SolverMVA(build(SchedStrategy.LRPT), "verbose", VerboseLevel.SILENT);
        lrpt.getAvg();
        assertTrue(srpt.result.RN.get(1, 0) < lrpt.result.RN.get(1, 0),
                "SRPT must serve the short class faster than LRPT");
        assertTrue(srpt.result.RN.get(1, 1) > lrpt.result.RN.get(1, 1),
                "LRPT must serve the long class faster than SRPT");
    }
}
