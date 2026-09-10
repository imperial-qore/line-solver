/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.tr;

import jline.VerboseLevel;
import jline.api.pfqn.nc.Pfqn_bk;
import jline.io.Ret;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Load concealment as a model TRANSFORMATION
 * ({@code options.config.transform='lc'}).
 *
 * <p>Birman-Kogan Algorithm 2 exists twice in the tree and the two are not the
 * same computation. {@code pfqn_bklc} is the KERNEL: it sweeps chains on a
 * demand matrix with MVA as the inner single-chain solve. The transformation is
 * the same sweep, but each single-chain subproblem is a real single-class
 * Network solved by whichever solver was called.
 *
 * <p>THE ORACLE IS THE KERNEL'S OWN FIXED POINT. On a model whose per-chain
 * subproblem is single-class, PS and exponential, the chain aggregation is
 * exact, so the two must reach the SAME fixed point in the SAME number of
 * sweeps. Agreement on the fixed point alone is not enough: a Jacobi sweep
 * reaches the same point at a different sweep count, and the sweep count is
 * what the four codebases are pinned on.
 *
 * <p>The MATLAB twin is line-test.git test_tr_lc.m, the python twin is
 * python/tests/test_tr_lc.py and the C++ twin is cpp/tests/test_tr_lc.cpp.
 */
public class LcStrategyTest {

    /** Delay -&gt; Q1 -&gt; Q2 -&gt; Q3 -&gt; Delay, two closed classes, one chain each. */
    private static Network build() {
        Network m = new Network("lc");
        Delay d = new Delay(m, "Think");
        Queue q1 = new Queue(m, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(m, "Q2", SchedStrategy.PS);
        Queue q3 = new Queue(m, "Q3", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(m, "C1", 4, d);
        ClosedClass c2 = new ClosedClass(m, "C2", 3, d);
        d.setService(c1, new Exp(1.0));
        d.setService(c2, new Exp(2.0));
        q1.setService(c1, new Exp(2.0));
        q1.setService(c2, new Exp(3.0));
        q2.setService(c1, new Exp(3.0));
        q2.setService(c2, new Exp(1.5));
        q3.setService(c1, new Exp(1.8));
        q3.setService(c2, new Exp(2.5));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c1, c1, m.serialRouting(d, q1, q2, q3));
        P.set(c2, c2, m.serialRouting(d, q1, q2, q3));
        m.link(P);
        return m;
    }

    /** The demand matrix the model reduces to, as the kernel is given it. */
    private static Matrix demands() {
        Matrix L = new Matrix(4, 2);
        L.set(1, 0, 0.5);
        L.set(1, 1, 1.0 / 3);
        L.set(2, 0, 1.0 / 3);
        L.set(2, 1, 2.0 / 3);
        L.set(3, 0, 5.0 / 9);
        L.set(3, 1, 0.4);
        return L;
    }

    private static SolverOptions elevating() {
        SolverOptions o = new SolverOptions(SolverType.CTMC);
        o.verbose = VerboseLevel.SILENT;
        o.config.put("transform", "lc");
        o.iter_max = 1000;
        return o;
    }

    private static Ret.pfqnBkLc kernel() {
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 4);
        N.set(0, 1, 3);
        Matrix Z = new Matrix(1, 2);
        Z.set(0, 0, 1.0);
        Z.set(0, 1, 0.5);
        return Pfqn_bk.pfqn_bklc(demands(), N, Z, "mva", 1e-10, 1000);
    }

    @Test
    public void transformationReachesTheKernelFixedPoint() throws Exception {
        Ret.pfqnBkLc k = kernel();
        SolverCTMC s = new SolverCTMC(build(), elevating());
        Matrix X = s.getAvgSysTput();
        assertEquals(k.X.length(), X.length());
        for (int r = 0; r < X.length(); r++) {
            assertEquals(k.X.get(r), X.get(r), 1e-8, "chain " + r);
        }
    }

    @Test
    public void theSweepIsGaussSeidelNotJacobi() throws Exception {
        // Same fixed point in the same number of sweeps. Jacobi coupling reaches
        // the same point at a DIFFERENT sweep count, which is what this catches.
        SolverCTMC s = new SolverCTMC(build(), elevating());
        s.getAvgSysTput();
        assertEquals(kernel().it, ((jline.solvers.ctmc.CTMCResult) s.result).iter);
    }

    @Test
    public void theTransformedSolveNamesItself() throws Exception {
        SolverCTMC s = new SolverCTMC(build(), elevating());
        s.getAvgSysTput();
        assertEquals("default/lc", ((jline.solvers.ctmc.CTMCResult) s.result).method);
        assertEquals(TransformMethod.LC, TransformMethod.canonical("lc"));
    }

    @Test
    public void anUnknownTokenIsRefusedByName() {
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> TransformMethod.canonical("bogus"));
        assertEquals(true, e.getMessage().contains("not a known model transformation"));
    }
}
