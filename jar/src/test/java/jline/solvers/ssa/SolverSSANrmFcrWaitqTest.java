package jline.solvers.ssa;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.Region;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the NRM's WAITQ-rule finite capacity regions against the exact CTMC,
 * with the region on an INTERIOR station.
 * <p>
 * The WAITQ rule PARKS a refused job in a per-region FIFO and admits it
 * head-of-line as capacity frees, in contrast to DROP, which destroys it. The
 * two are observationally identical when the region is fed by a Source (a
 * Source's population is fictitious), so a genuine WAITQ test must place the
 * region on a station fed by a real upstream queue. The interior fixture below
 * is a closed cycle Think(INF) -&gt; Q1/FCFS -&gt; Q2/FCFS -&gt; Think with the
 * region over {Q2}; the parked jobs are conserved (they re-enter Q2 later), so
 * the sum of station queue lengths is strictly below the closed population N,
 * the difference being the mean FIFO occupancy. A DROP rule there would destroy
 * jobs and drain the closed network -- a qualitatively different answer.
 * </p>
 * <p>
 * Ground truth is the exact CTMC solved in-process. The test also asserts the
 * result method is "nrm": Solver_ssa_analyzer downgrades an ineligible model to
 * the serial engine, and the serial engine is correct here, so a silent fallback
 * would pass this test while never exercising the WAITQ FIFO it exists to cover.
 * </p>
 */
public class SolverSSANrmFcrWaitqTest {

    private static final int SAMPLES = 300000;
    private static final int[] SEEDS = {23000, 24000, 25000, 26000};
    /** Absolute tolerance on the simulated means, averaged over SEEDS. */
    private static final double ATOL = 0.03;

    /** Closed interior-region WAITQ model: Think -> Q1 -> Q2{region cap} -> Think. */
    private static Network interiorModel() {
        Network model = new Network("wq_interior");
        Delay think = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c1 = new ClosedClass(model, "C1", 4, think);
        think.setService(c1, new Exp(1.0));
        q1.setService(c1, new Exp(1.5));
        q2.setService(c1, new Exp(1.2));
        model.link(model.serialRouting(think, q1, q2));
        Region region = model.addRegion(Arrays.<Node>asList((Node) q2));
        region.setGlobalMaxJobs(2);
        region.setDropRule(c1, DropStrategy.WaitingQueue);
        return model;
    }

    /** Open source-fed WAITQ model (non-regression): Source -> Q1{region cap} -> Sink. */
    private static Network sourceFedModel() {
        Network model = new Network("wq_srcfed");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass c1 = new OpenClass(model, "C1");
        source.setArrival(c1, new Exp(0.6));
        q1.setService(c1, new Exp(1.0));
        model.link(model.serialRouting(source, q1, sink));
        Region region = model.addRegion(Arrays.<Node>asList((Node) q1));
        region.setGlobalMaxJobs(3);
        region.setDropRule(c1, DropStrategy.WaitingQueue);
        return model;
    }

    private static Matrix ctmcQLen(Network model, int cutoff) {
        SolverOptions opt = new SolverOptions();
        opt.cutoff(cutoff).keep(false);
        return new SolverCTMC(model, opt).getAvgQLen();
    }

    private static Matrix ctmcTput(Network model, int cutoff) {
        SolverOptions opt = new SolverOptions();
        opt.cutoff(cutoff).keep(false);
        return new SolverCTMC(model, opt).getAvgTput();
    }

    /** Mean over SEEDS of an NRM per-station average, asserting method == nrm.
     *  metric==0 selects queue length, metric==1 selects throughput. */
    private static Matrix nrmMean(Network model, int metric) {
        Matrix acc = null;
        for (int seed : SEEDS) {
            SolverSSA solver = new SolverSSA(model);
            solver.options.method = "nrm";
            solver.options.samples = SAMPLES;
            solver.options.seed = seed;
            Matrix q = (metric == 0) ? solver.getAvgQLen() : solver.getAvgTput();
            assertTrue(solver.result != null && solver.result.method != null
                            && solver.result.method.contains("nrm"),
                    "NRM did not run: method was "
                            + (solver.result == null ? "<no result>" : solver.result.method));
            if (acc == null) {
                acc = q.copy();
            } else {
                acc = acc.add(1.0, q);
            }
        }
        return acc.scale(1.0 / SEEDS.length);
    }

    private static Matrix nrmQLenMean(Network model) {
        return nrmMean(model, 0);
    }

    private static void assertCloseAbs(double got, double exact, String what) {
        double err = Math.abs(got - exact);
        assertTrue(err < ATOL, what + ": NRM returned " + got + ", CTMC is " + exact
                + " (abs err " + err + ")");
    }

    /**
     * Interior region {Q2} in a closed cycle: refused jobs are parked in the
     * region FIFO, so every station queue length must match the exact CTMC and
     * the parked mass (N minus the sum of station queue lengths) must be strictly
     * positive -- the observable WAITQ signature that a DROP rule cannot produce.
     */
    @Test
    public void testFcrWaitqInteriorRegionMatchesCtmc() {
        Network m = interiorModel();
        Matrix ctmc = ctmcQLen(m, 8);
        Matrix nrm = nrmQLenMean(interiorModel());
        double sumCtmc = 0.0;
        for (int i = 0; i < ctmc.getNumRows(); i++) {
            assertCloseAbs(nrm.get(i, 0), ctmc.get(i, 0), "interior station " + i + " QLen");
            sumCtmc += ctmc.get(i, 0);
        }
        // The parked jobs live in the region FIFO, not at any station, so the
        // station queue lengths must sum to strictly below the closed population.
        assertTrue(sumCtmc < 3.9,
                "WAITQ FIFO must hold jobs off-station: sum of station QLen was "
                        + sumCtmc + " (population is 4)");
        // Throughput at every station must also match the exact CTMC.
        Matrix ctmcX = ctmcTput(m, 8);
        Matrix nrmX = nrmMean(interiorModel(), 1);
        for (int i = 0; i < ctmcX.getNumRows(); i++) {
            assertCloseAbs(nrmX.get(i, 0), ctmcX.get(i, 0), "interior station " + i + " Tput");
        }
    }

    /**
     * Source-fed region {Q1}: the pre-existing configuration, kept as a
     * non-regression. A source-fed WAITQ region blocks the source, so the model
     * is bounded and its station queue length matches the exact CTMC.
     */
    @Test
    public void testFcrWaitqSourceFedRegionMatchesCtmc() {
        Network m = sourceFedModel();
        Matrix ctmc = ctmcQLen(m, 20);
        Matrix nrm = nrmQLenMean(sourceFedModel());
        // Station 1 is Q1 (Source is station 0).
        assertCloseAbs(nrm.get(1, 0), ctmc.get(1, 0), "source-fed Q1 QLen");
        // Throughput at Q1 must also match the exact CTMC.
        Matrix ctmcX = ctmcTput(m, 20);
        Matrix nrmX = nrmMean(sourceFedModel(), 1);
        assertCloseAbs(nrmX.get(1, 0), ctmcX.get(1, 0), "source-fed Q1 Tput");
    }
}
