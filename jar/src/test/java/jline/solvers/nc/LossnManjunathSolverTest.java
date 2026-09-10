/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import java.util.List;

import jline.VerboseLevel;
import jline.examples.java.advanced.FCRegionModel;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.Region;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * End-to-end test of the exact loss-network path in SolverNC.
 *
 * The expected values are the MATLAB reference implementation's output for the
 * same model (FCRegionModel.fcr_lossn, a Source feeding a two-class Delay inside
 * a DROP region with a global cap of 5 and a per-class cap of 3), obtained from
 * SolverNC(model) with the default method, which resolves to lossn.exact. They
 * are reproduced here to 16 digits because the transform is exact: there is no
 * tolerance to hide behind, and any drift is a defect rather than noise.
 *
 * Two things this pins that the erlangfp/mci tests cannot:
 *
 * 1. THE OFFERED LOAD IS SCALED. nu_r is the arrival rate times the mean holding
 *    time in the region, lambda_r V_r / mu_r, not the bare arrival rate. Class 2
 *    here has mu = 0.8, so the bare rate would give nu_2 = 0.2 instead of 0.25
 *    and report the blocking of a different network. The approximate methods share
 *    this input, so their own tests would have concealed it just as well.
 *
 * 2. QLen IS THE CARRIED LOAD, NOT A THROUGHPUT. lossn_manjunath returns E[n_r]
 *    directly, so Q is QLen and the throughput is lambda_r (1 - Loss_r); dividing
 *    QLen by mu again would double-apply Little's law.
 */
public class LossnManjunathSolverTest {

    /** Exact, so the assertion tolerance is numerical noise only. */
    private static final double TOL = 1e-9;

    // MATLAB SolverNC(fcr_lossn) with method 'default' -> 'default/lossn.exact'
    private static final double Q1 = 0.2989813602605746;
    private static final double Q2 = 0.2494742965193094;
    private static final double X1 = 0.2989813602605746;
    private static final double X2 = 0.1995794372154475;
    private static final double LG = 0.5495940110776684;
    private static final double R1 = 1.0;
    private static final double R2 = 1.25;

    /** Row order of the avg table is Source then Delay, two classes each. */
    private static final int DELAY_C1 = 2;
    private static final int DELAY_C2 = 3;

    @Test
    public void testDefaultResolvesToExactAndMatchesMatlab() {
        Network model = FCRegionModel.fcr_lossn();
        SolverNC solver = new SolverNC(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable t = solver.getAvgTable();

        assertEquals(Q1, t.getQLen().get(DELAY_C1), TOL, "carried load of class 1 at the delay");
        assertEquals(Q2, t.getQLen().get(DELAY_C2), TOL, "carried load of class 2 at the delay");
        // An infinite server never queues, so the "utilization" is the mean number
        // of busy servers, which is the population itself.
        assertEquals(Q1, t.getUtil().get(DELAY_C1), TOL, "busy servers of class 1");
        assertEquals(Q2, t.getUtil().get(DELAY_C2), TOL, "busy servers of class 2");
        assertEquals(R1, t.getRespT().get(DELAY_C1), TOL, "response time is the service time");
        assertEquals(R2, t.getRespT().get(DELAY_C2), TOL, "response time is the service time");
        assertEquals(X1, t.getTput().get(DELAY_C1), TOL, "carried throughput of class 1");
        assertEquals(X2, t.getTput().get(DELAY_C2), TOL, "carried throughput of class 2");
    }

    @Test
    public void testExactReportsTheExactNormalizingConstant() {
        Network model = FCRegionModel.fcr_lossn();
        SolverNC solver = new SolverNC(model, "verbose", VerboseLevel.SILENT);
        solver.getAvgTable();
        // erlangfp has no normalizing constant to report and leaves this NaN; the
        // transform computes g(C) itself, so a finite value here is what
        // distinguishes the exact path from the approximation.
        assertTrue(Double.isFinite(((NCResult) solver.result).prob.logNormConstAggr),
                "the exact transform must report a finite log g(C)");
        assertEquals(LG, ((NCResult) solver.result).prob.logNormConstAggr, TOL, "log of the exact normalizing constant");
        assertEquals(1, solver.result.iter, "the transform is direct, so one iteration");
    }

    @Test
    public void testExplicitMsTokenTakesTheExactPath() {
        // 'ms' names the transform explicitly; it must not silently fall back to
        // the Erlang approximation, which would report an approximation under the
        // name of an exact method.
        Network model = FCRegionModel.fcr_lossn();
        SolverNC solver = new SolverNC(model, "method", "ms", "verbose", VerboseLevel.SILENT);
        NetworkAvgTable t = solver.getAvgTable();
        assertEquals(Q1, t.getQLen().get(DELAY_C1), TOL, "'ms' must reach the exact transform");
        assertEquals(Q2, t.getQLen().get(DELAY_C2), TOL, "'ms' must reach the exact transform");
        assertEquals(LG, ((NCResult) solver.result).prob.logNormConstAggr, TOL, "and report the exact normalizing constant");
    }

    @Test
    public void testErlangfpIsCloseButNotExact() {
        // The reduced-load approximation treats the links as independent, so it
        // must land near the truth without reproducing it. Asserting BOTH bounds
        // is what keeps the test honest: an 'erlangfp' that silently dispatched to
        // the transform would pass a one-sided closeness check.
        Network model = FCRegionModel.fcr_lossn();
        SolverNC solver = new SolverNC(model, "method", "erlangfp", "verbose", VerboseLevel.SILENT);
        NetworkAvgTable t = solver.getAvgTable();
        double q1 = t.getQLen().get(DELAY_C1);
        double q2 = t.getQLen().get(DELAY_C2);
        assertEquals(Q1, q1, 5e-3, "erlangfp must be close on class 1");
        assertEquals(Q2, q2, 5e-3, "erlangfp must be close on class 2");
        assertTrue(Math.abs(q1 - Q1) > 1e-12 || Math.abs(q2 - Q2) > 1e-12,
                "erlangfp is an approximation and must not reproduce the exact answer");
    }

    /**
     * fcr_lossn with the region's per-class drop rules overridden.
     *
     * false means WaitingQueue, which holds the arrival back instead of
     * discarding it.
     */
    private static Network fcrLossnWithDropRules(boolean drop1, boolean drop2) {
        Network model = FCRegionModel.fcr_lossn();
        Region fcr = model.getRegions().get(0);
        List<JobClass> classes = model.getClasses();
        fcr.setDropRule(classes.get(0), drop1);
        fcr.setDropRule(classes.get(1), drop2);
        return model;
    }

    /**
     * The loss-network path requires EVERY class to be dropped, as the reference
     * tests with all(regionrule(1,:) == DROP).
     *
     * THE FIRST CASE IS THE REGRESSION. The dispatch used to read
     * sn.regionrule.get(0), i.e. class 0 alone, so a region that discards class 1
     * and holds class 2 back passed the test and was silently solved as a loss
     * network: a wrong number, not an error. Asserting BOTH orderings is what
     * separates the fix from the bug, since the reversed case was refused
     * correctly even before it.
     *
     * A mixed region is not a loss network at all. Its blocked class occupies the
     * region while it waits, so the per-class loss probabilities Kelly's
     * truncation implies are not the ones the model implies; there is no
     * truncated product form to evaluate.
     *
     * Verified against MATLAB SolverNC on the same four configurations: only
     * [1 1] solves, and [1 -1], [-1 1], [-1 -1] all raise this message.
     */
    @Test
    public void testMixedDropRulesAreRefused() {
        boolean[][] mixed = {{true, false}, {false, true}, {false, false}};
        String[] label = {"class 1 DROP, class 2 WAITQ", "class 1 WAITQ, class 2 DROP",
                          "both WAITQ"};
        for (int i = 0; i < mixed.length; i++) {
            final Network model = fcrLossnWithDropRules(mixed[i][0], mixed[i][1]);
            RuntimeException e = assertThrows(RuntimeException.class,
                    new org.junit.jupiter.api.function.Executable() {
                        public void execute() {
                            new SolverNC(model, "verbose", VerboseLevel.SILENT).getAvgTable();
                        }
                    }, label[i] + " must not reach the loss-network analyzer");
            assertTrue(e.getMessage().contains("WAITQ"),
                    label[i] + " must be refused by naming the blocking policy, got: "
                            + e.getMessage());
        }
    }

    @Test
    public void testAllDropStillSolvesAfterRewritingTheRules() {
        // The companion to the refusal above: setting both rules explicitly to
        // DROP must still take the exact path, so the all-classes test is not
        // simply rejecting everything.
        Network model = fcrLossnWithDropRules(true, true);
        SolverNC solver = new SolverNC(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable t = solver.getAvgTable();
        assertEquals(Q1, t.getQLen().get(DELAY_C1), TOL, "an all-DROP region still solves");
        assertEquals(Q2, t.getQLen().get(DELAY_C2), TOL, "an all-DROP region still solves");
    }

    @Test
    public void testMciBracketsTheExactAnswer() {
        // The Monte Carlo estimator is unbiased, so with a large sample it must
        // agree with the transform on both the metrics and log g(C).
        Network model = FCRegionModel.fcr_lossn();
        SolverNC solver = new SolverNC(model, "method", "mci", "verbose", VerboseLevel.SILENT);
        solver.options.samples = 200000;
        solver.options.seed = 42;
        NetworkAvgTable t = solver.getAvgTable();
        assertEquals(Q1, t.getQLen().get(DELAY_C1), 5e-3, "mci must bracket class 1");
        assertEquals(Q2, t.getQLen().get(DELAY_C2), 5e-3, "mci must bracket class 2");
        assertEquals(LG, ((NCResult) solver.result).prob.logNormConstAggr, 5e-3, "mci must bracket log g(C)");
    }
}
