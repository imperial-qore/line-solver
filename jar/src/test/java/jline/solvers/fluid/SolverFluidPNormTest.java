package jline.solvers.fluid;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.SolverResult;

/**
 * Wiring of the p-norm smoothing in the SolverFluid {@code matrix} method.
 *
 * <p>Ruuskanen et al., PEVA 151 (2021) smooth the PS min() of eq. (12) with the
 * inverse p-norm of eq. (26). Two consequences a mean-only comparison hides:</p>
 *
 * <ol>
 * <li>An INF station has k = infinity, so min(k, sum x) = sum x identically and
 * there is nothing to smooth. Sa carries the total population at a delay, so
 * smoothing it would run the delay as a k = N queue.</li>
 * <li>Eq. (23) reads the utilization off the SAME share the drift integrated,
 * k rho / E[sum X] = ghat, so the metrics may not revert to the hard min().</li>
 * </ol>
 *
 * <p>The M/M/1 case is the paper's own Example 1: at p = 1 the smoothed model IS
 * the Tipper/PSFFA model, whose fixed point is the exact mean rho/(1-rho), where
 * the unsmoothed mean-field model returns lambda. Mirrors
 * python/tests/test_fld_pnorm_wiring.py, cpp/tests/test_fluid_pnorm.cpp and
 * line-test.git/test_solver_fld_pnorm.m.
 */
public class SolverFluidPNormTest {

    /** Delay(rate 1) + PS Queue(rate 2), one closed class of 10. */
    private static Network cqn() {
        Network model = new Network("pnorm_cqn");
        Delay delay = new Delay(model, "D");
        Queue queue = new Queue(model, "Q", SchedStrategy.PS);
        ClosedClass jobclass = new ClosedClass(model, "C", 10, delay);
        delay.setService(jobclass, new Exp(1.0));
        queue.setService(jobclass, new Exp(2.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static Network mm1(double rho) {
        Network model = new Network("pnorm_mm1");
        Source source = new Source(model, "S");
        Queue queue = new Queue(model, "Q", SchedStrategy.PS);
        Sink sink = new Sink(model, "K");
        OpenClass jobclass = new OpenClass(model, "C");
        source.setArrival(jobclass, new Exp(rho));
        queue.setService(jobclass, new Exp(1.0));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    private static SolverResult solve(Network model, double pstar) {
        SolverFluid solver = new SolverFluid(model);
        solver.options.method = "matrix";
        solver.options.verbose = VerboseLevel.SILENT;
        if (pstar > 0) {
            List<Double> ps = new ArrayList<Double>();
            for (int i = 0; i < model.getNumberOfStations(); i++) {
                ps.add(pstar);
            }
            solver.options.config.pstar = ps;
        }
        solver.getAvg();
        return solver.result;
    }

    @Test
    public void testSmoothingLeavesAnInfStationAlone() {
        // every job at a delay is in service, so its departure rate is Q * mu
        // with no share to apply; smoothing it as a k = N queue returns Q * ghat
        Network model = cqn();
        double[] exponents = new double[]{1.0, 4.0, 20.0};
        for (int i = 0; i < exponents.length; i++) {
            SolverResult r = solve(cqn(), exponents[i]);
            assertEquals(r.QN.get(0, 0), r.TN.get(0, 0), 1e-9,
                    "delay smoothed at pstar=" + exponents[i]);
        }
        // a large exponent recovers the hard min(), i.e. the unsmoothed answer
        SolverResult hard = solve(model, 0.0);
        SolverResult soft = solve(cqn(), 20.0);
        assertEquals(hard.QN.get(0, 0), soft.QN.get(0, 0), 1e-6);
    }

    @Test
    public void testMetricsReadTheSmoothedShare() {
        // T is read off theta, so a metric taken from the hard min() while x
        // came from the smoothed drift shows as a flow imbalance in the cycle
        double[] exponents = new double[]{1.0, 4.0};
        for (int i = 0; i < exponents.length; i++) {
            SolverResult r = solve(cqn(), exponents[i]);
            assertEquals(r.TN.get(0, 0), r.TN.get(1, 0), 1e-6,
                    "flow imbalance at pstar=" + exponents[i]);
        }
    }

    @Test
    public void testUnitExponentRecoversTheExactMM1Mean() {
        double[] rhos = new double[]{0.3, 0.5, 0.7};
        for (int i = 0; i < rhos.length; i++) {
            double rho = rhos[i];
            SolverResult soft = solve(mm1(rho), 1.0);
            assertEquals(rho / (1.0 - rho), soft.QN.get(1, 0), 1e-6);
            assertEquals(rho, soft.UN.get(1, 0), 1e-6); // ODE tol; without the fix U would be 1.0
            // the unsmoothed mean-field model returns lambda, the failure the
            // smoothing exists to repair
            SolverResult hard = solve(mm1(rho), 0.0);
            assertEquals(rho, hard.QN.get(1, 0), 1e-6);
        }
    }
}
