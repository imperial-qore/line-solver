/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fj;

import jline.examples.java.models.Gallery;
import jline.lang.Network;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * The fork-join fixed point on a QUORUM join, against MATLAB.
 *
 * <p>{@code gallery_fj_quorum} is a closed 3-branch fork-join whose join fires on the SECOND
 * sibling. The reference values are MATLAB's on the same model
 * (X = 1.5854962587 under MVA, 1.64179603792 under NC); native python agrees with MATLAB to 1e-9
 * on both. They are NOT the exact throughput: the MMT floors the synchronisation delay at zero
 * here and under-states X, which SolverLDES puts at 1.847 -- see _kb/05-solvers-overview.md. What
 * this test pins is that the three analytical codebases charge the SAME order statistic.</p>
 */
public class FJQuorumDriverTest {

    /** Relative tolerance: MATLAB, the JAR and python differ by ~1e-7 on the NC route. */
    private static final double RELTOL = 1e-4;

    private static double throughputOfClass(Matrix TN) {
        return TN.get(0, 0);
    }

    @Test
    public void testQuorumUnderMVAMatchesMatlab() {
        Network model = Gallery.gallery_fj_quorum();
        SolverMVA solver = new SolverMVA(model);
        double X = throughputOfClass(solver.getAvgTput());
        assertEquals(1.5854962587, X, RELTOL * 1.5854962587);
    }

    @Test
    public void testQuorumUnderNCMatchesMatlab() {
        Network model = Gallery.gallery_fj_quorum();
        SolverNC solver = new SolverNC(model);
        double X = throughputOfClass(solver.getAvgTput());
        assertEquals(1.64179603792, X, RELTOL * 1.64179603792);
    }

    /**
     * The same model with the quorum removed must fall back on the ordinary AND-join, i.e. a
     * strictly LONGER cycle and so a lower throughput.
     */
    @Test
    public void testQuorumRaisesThroughputAboveTheFullJoin() {
        Network quorum = Gallery.gallery_fj_quorum();
        double Xq = throughputOfClass(new SolverMVA(quorum).getAvgTput());
        Network full = Gallery.gallery_fj_quorum();
        jline.lang.nodes.Join join = (jline.lang.nodes.Join) full.getNodeByName("Join");
        join.setStrategy(full.getClasses().get(0), jline.lang.constant.JoinStrategy.STD);
        join.setRequired(full.getClasses().get(0), -1);
        double Xf = throughputOfClass(new SolverMVA(full).getAvgTput());
        org.junit.jupiter.api.Assertions.assertTrue(Xq > Xf,
                "a 2-of-3 quorum must complete faster than a 3-of-3 join, got " + Xq + " vs " + Xf);
    }
}
