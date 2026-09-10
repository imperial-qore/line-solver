package jline.solvers.fluid;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.VerboseLevel;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.NHPP;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

/**
 * Regression tests for the SolverFluid non-homogeneous Poisson (NHPP) transient.
 *
 * <p>A cyclic NHPP source intensity lambda(t) is injected into the closing fluid
 * ODE as a per-event rate multiplier. With a fast server the queue is never the
 * bottleneck, so its throughput must track lambda(t) segment by segment. The
 * steady-state getAvg is deliberately unaffected and reports the time-average
 * rate.
 *
 * <p>The SOURCE throughput is a fluid open-model artifact and is not asserted
 * on; observe the QUEUE. Mirrors line-test.git/test/testsFLD/test_solver_fld_nhpp.m and
 * python/line_solver/tests/test_solver_fld_nhpp.py.
 */
public class SolverFluidNHPPTest {

    private static final double[] BREAKPOINTS = new double[]{0, 3, 4, 6};
    private static final double[] RATES = new double[]{2, 8, 4};
    private static final double TEND = 12.0; // two full periods

    private static Network buildModel(jline.lang.processes.Distribution arrival) {
        Network model = new Network("fld_nhpp");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "Class1");
        source.setArrival(jobclass, arrival);
        queue.setService(jobclass, new Exp(50)); // fast server: queue tracks input
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    /** Linear interpolation of a (time, value) trajectory at tt. */
    private static double interp(Matrix traj, Matrix times, double tt) {
        int n = traj.getNumRows();
        double prevT = times.get(0, 0);
        for (int i = 1; i < n; i++) {
            double curT = times.get(i, 0);
            if (tt <= curT) {
                double dt = curT - prevT;
                double w = (dt > 0) ? (tt - prevT) / dt : 0.0;
                return (1 - w) * traj.get(i - 1, 0) + w * traj.get(i, 0);
            }
            prevT = curT;
        }
        return traj.get(n - 1, 0);
    }

    @Test
    public void testTransientTracksIntensity() {
        NHPP nhpp = new NHPP(BREAKPOINTS, RATES, true);
        Network model = buildModel(nhpp);
        SolverFluid solver = new SolverFluid(model);
        solver.options.timespan = new double[]{0, TEND};
        solver.options.verbose = VerboseLevel.SILENT;
        solver.getTranAvg();
        SolverResult result = solver.result;

        Matrix times = result.t;
        Matrix queueTput = result.TNt[1][0];

        double maxRel = 0;
        for (int k = 0; k < 2; k++) {
            for (int j = 0; j < RATES.length; j++) {
                double tt = k * nhpp.getPeriod() + (BREAKPOINTS[j] + BREAKPOINTS[j + 1]) / 2.0;
                if (tt <= 0.5 || tt >= TEND) {
                    continue; // skip the initial empty-system warmup
                }
                double expected = nhpp.getRateAt(tt);
                double actual = interp(queueTput, times, tt);
                maxRel = Math.max(maxRel, Math.abs(actual - expected) / expected);
            }
        }
        assertTrue(maxRel < 0.05,
                "queue throughput deviates by " + maxRel + " from lambda(t)");

        // The transient must actually vary: a solver that silently used the
        // time-average rate would give a flat trajectory.
        double lo = Double.POSITIVE_INFINITY;
        double hi = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < queueTput.getNumRows(); i++) {
            if (times.get(i, 0) > 0.5) {
                lo = Math.min(lo, queueTput.get(i, 0));
                hi = Math.max(hi, queueTput.get(i, 0));
            }
        }
        double spread = (hi - lo) / nhpp.getTimeAverageRate();
        assertTrue(spread > 0.5, "trajectory is flat (spread " + spread + "); lambda(t) was ignored");
    }

    @Test
    public void testSteadyStateIsTimeAverage() {
        NHPP nhpp = new NHPP(BREAKPOINTS, RATES, true);
        Network model = buildModel(nhpp);
        SolverFluid solver = new SolverFluid(model);
        solver.options.verbose = VerboseLevel.SILENT;
        solver.getAvg();
        double lamAvg = nhpp.getTimeAverageRate();
        assertEquals(1.0, solver.result.TN.get(1, 0) / lamAvg, 0.05);
    }

    @Test
    public void testHomogeneousModelUnperturbed() {
        double lamAvg = new NHPP(BREAKPOINTS, RATES, true).getTimeAverageRate();
        Network model = buildModel(new Exp(lamAvg));
        SolverFluid solver = new SolverFluid(model);
        solver.options.verbose = VerboseLevel.SILENT;
        solver.getAvg();
        assertEquals(lamAvg, solver.result.TN.get(1, 0), 1e-6);
        // The default method is now the second-order closure "minnormal", so the
        // queue length is its answer rather than the first-order rho = lambda/50
        // (0.0733333). MATLAB reference for this model: 0.0733548838, which the
        // first-order "closing" method still returns as exactly rho. Throughput is
        // unaffected: flow balance makes it exact under either closure.
        assertEquals(0.0733548838, solver.result.QN.get(1, 0), 1e-6);
    }

    @Test
    public void testFeatureSetDeclaresNhpp() {
        assertTrue(SolverFluid.getFeatureSet().inspectFeature("NHPP"));
    }
}
