package jline.solvers.ssa;

import jline.lang.Region;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import jline.lang.nodes.Node;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the NRM's DROP-rule finite capacity regions against the exact CTMC
 * solution, with the region on an INTERIOR station.
 * <p>
 * The interior placement is the whole point of this fixture. Under DROP the
 * refused job is DESTROYED: the departure fires at its full rate and the job
 * never reaches the destination. A previous implementation instead CENSORED the
 * refused transition by scaling the propensity with the admitted share, which
 * holds the job at its SOURCE. The two are observationally identical when the
 * region is fed by a Source -- a Source's population is fictitious, so "destroy
 * the arriving job" and "censor the arrival" cannot be told apart -- and every
 * pre-existing FCR fixture placed the region exactly there. They diverge
 * qualitatively as soon as the blocked job would be leaving a real queue: the
 * upstream queue then grows without bound and nothing is ever lost.
 * </p>
 * <p>
 * Fixture: Source(1.0) -&gt; Q1/PS(3.0) -&gt; Q2/PS(1.5) -&gt; Sink, region over
 * {Q2}, globalMaxJobs = 1, DropStrategy.Drop. The exact CTMC (cutoff-invariant
 * from 16 upward) gives Q1 = 0.5 and T_Q2 = 0.6, i.e. 40% of jobs are lost at
 * the boundary and Q1 is stable. The censoring implementation returned a
 * non-stationary Q1 (108.7 / 112.1 / 303.6 at 100k / 400k / 1.6M samples) and
 * T_Q2 -&gt; 0.995.
 * </p>
 * <p>
 * The exact values are hard-coded rather than obtained from a CTMC solve so the
 * assertion cannot drift with the reference solver, and because the Q1 = 0.5
 * figure is analytic: with the boundary loss, Q1 is an M/M/1-PS at
 * utilisation 1/3.
 * </p>
 */
public class SolverSSANrmFcrDropTest {

    private static final int SAMPLES = 400000;
    private static final int[] SEEDS = {23000, 24000, 25000, 26000};
    /** Relative tolerance on the simulated means, averaged over SEEDS. */
    private static final double RTOL = 0.02;

    private static Network model(boolean interior, int cap) {
        Network model = new Network("fcr_drop");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass c1 = new OpenClass(model, "C1");
        source.setArrival(c1, new Exp(1.0));
        q1.setService(c1, new Exp(3.0));
        q2.setService(c1, new Exp(1.5));
        model.link(model.serialRouting(source, q1, q2, sink));
        Region region = model.addRegion(
                Arrays.<Node>asList(interior ? q2 : q1));
        region.setGlobalMaxJobs(cap);
        region.setDropRule(c1, DropStrategy.Drop);
        return model;
    }

    /** Mean over SEEDS of (QLen at Q1, Tput at Q2). Station rows: 1 = Q1, 2 = Q2. */
    private static double[] meanQ1AndT2(boolean interior, int cap) {
        double q1 = 0.0;
        double t2 = 0.0;
        for (int seed : SEEDS) {
            SolverSSA solver = new SolverSSA(model(interior, cap));
            solver.options.method = "nrm";
            solver.options.samples = SAMPLES;
            solver.options.seed = seed;
            Matrix q = solver.getAvgQLen();
            Matrix t = solver.getAvgTput();
            // The NRM must actually have run. Solver_ssa_analyzer downgrades an
            // ineligible model to the serial engine, and the serial engine is
            // CORRECT here -- so a silent fallback would make this test pass
            // while never exercising the code it exists to cover.
            assertTrue(solver.result != null && solver.result.method != null
                            && solver.result.method.contains("nrm"),
                    "NRM did not run: method was "
                            + (solver.result == null ? "<no result>" : solver.result.method));
            q1 += q.get(1, 0);
            t2 += t.get(2, 0);
        }
        return new double[]{q1 / SEEDS.length, t2 / SEEDS.length};
    }

    private static void assertClose(double got, double exact, String what) {
        double err = Math.abs(got - exact) / exact;
        assertTrue(err < RTOL, what + ": NRM returned " + got + ", exact is " + exact
                + " (" + (100.0 * err) + "% off)");
    }

    /**
     * Region on the interior station Q2: the refused job must be destroyed, so
     * Q1 stays stable and 40% of the throughput is lost. Censoring the refused
     * departure instead makes Q1 diverge.
     */
    @Test
    public void testFcrDropInteriorRegionLosesJobs() {
        double[] r = meanQ1AndT2(true, 1);
        assertClose(r[0], 0.5, "interior region Q1");
        assertClose(r[1], 0.6, "interior region T_Q2");
    }

    /**
     * Region on the Source-fed station Q1: the pre-existing configuration, which
     * must not regress. Exact CTMC gives Q1 = 0.38462, T_Q2 = 0.92307.
     */
    @Test
    public void testFcrDropSourceFedRegionUnchanged() {
        double[] r = meanQ1AndT2(false, 2);
        assertClose(r[0], 0.38462, "source-fed region Q1");
        assertClose(r[1], 0.92307, "source-fed region T_Q2");
    }
}
