/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.VerboseLevel;
import jline.examples.java.advanced.FCRegionModel;
import jline.lang.Network;
import jline.solvers.NetworkAvgTable;
import jline.solvers.nc.NCResult;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * End-to-end test of loss-network method selection in SolverNC: the
 * Erlang fixed point ("erlangfp", default) and the Monte Carlo summation
 * ("mci", Ross-Wang 1992) must agree on carried throughput, and "mci" must
 * additionally report a finite log normalization constant.
 */
public class LossnMciSolverTest {

    private static double tputAt(NetworkAvgTable t, int row) {
        return t.getTput().get(row);
    }

    private static int nrows(NetworkAvgTable t) {
        return t.getTput().size();
    }

    @Test
    public void testErlangfpVsMciAgree() {
        Network m1 = FCRegionModel.fcr_lossn();
        SolverNC fp = new SolverNC(m1, "method", "erlangfp", "verbose", VerboseLevel.SILENT);
        NetworkAvgTable tfp = fp.getAvgTable();

        Network m2 = FCRegionModel.fcr_lossn();
        SolverNC mc = new SolverNC(m2, "method", "mci", "verbose", VerboseLevel.SILENT);
        mc.options.samples = 200000;
        mc.options.seed = 5;
        NetworkAvgTable tmc = mc.getAvgTable();

        // Carried throughput per class must agree between the two methods.
        for (int r = 0; r < nrows(tfp); r++) {
            double a = tputAt(tfp, r);
            double b = tputAt(tmc, r);
            assertEquals(a, b, 5e-3, "throughput row " + r + ": erlangfp=" + a + " mci=" + b);
        }

        // mci reports a finite normalization constant; erlangfp does not.
        double lG = ((NCResult) mc.result).prob.logNormConstAggr;
        assertTrue(Double.isFinite(lG), "mci lG must be finite, got " + lG);
    }
}
