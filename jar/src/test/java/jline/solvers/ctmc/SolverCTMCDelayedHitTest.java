package jline.solvers.ctmc;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.examples.java.basic.CacheRetrievalSystemModel;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Cache;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * SolverCTMC reports the delayed-hit queue length of a retrieval cache exactly.
 *
 * <p>Block A of the cache local state marks the items being fetched and block B counts
 * the secondary requests merged onto each fetch, so per item i
 * phi_i = P(block A bit i set), d1_i = E[block B slots of i] and dfull_i = d1_i + phi_i
 * are state rewards of the stationary distribution. The delayed-hit RATE is a transition
 * reward over the generator (the transitions that clear block A bit i), which splits the
 * hit-class rate into true hits and delayed hits.</p>
 *
 * <p>Model: open, 3 items, cache capacity 1, FIFO, Exp(1) arrivals, Exp(2) fetch. The
 * cutoff truncates block B at maxPending = cutoff - 1, so every finite cutoff gives a
 * LOWER bound on the delayed-hit mass; the values below are exact for their cutoff and
 * approach the NC exact reference from below.</p>
 */
public class SolverCTMCDelayedHitTest {

    private static final double TOL = 1e-4;

    @BeforeAll
    public static void setUp() {
        GlobalConstants.getInstance().setVerbose(VerboseLevel.SILENT);
    }

    private static Network model() {
        return CacheRetrievalSystemModel.simple_retrieval_system_model(
                new double[]{0.6, 0.3, 0.1}, new double[]{2.0, 2.0, 2.0},
                "[1]", SchedStrategy.INF);
    }

    private static Cache cacheOf(Network model) {
        for (int i = 0; i < model.getNodes().size(); i++) {
            if (model.getNodes().get(i) instanceof Cache) {
                return (Cache) model.getNodes().get(i);
            }
        }
        throw new IllegalStateException("no cache node in model");
    }

    private static Cache solve(int cutoff) {
        Network model = model();
        SolverOptions o = new SolverOptions(SolverType.CTMC);
        Matrix co = new Matrix(1, 1);
        co.set(0, 0, cutoff);
        o.cutoff = co;
        new SolverCTMC(model, o).getAvgTable();
        return cacheOf(model);
    }

    @Test
    public void testDelayedHitQLenCutoff2() {
        Cache c = solve(2);
        Matrix d1 = c.getDelayedHitQLen();
        Matrix dfull = c.getDelayedHitQLenFull();
        assertEquals(3, d1.getNumCols());
        assertEquals(0.037531, d1.get(0, 0), TOL);
        assertEquals(0.020324, d1.get(0, 1), TOL);
        assertEquals(0.003728, d1.get(0, 2), TOL);
        assertEquals(0.138703, dfull.get(0, 0), TOL);
        assertEquals(0.110729, dfull.get(0, 1), TOL);
        assertEquals(0.047859, dfull.get(0, 2), TOL);
    }

    @Test
    public void testDfullMinusD1IsFetchProbability() {
        // dfull_i - d1_i = phi_i = P(a fetch of item i is in flight), by construction:
        // dfull counts the triggering request, d1 does not.
        Cache c = solve(2);
        Matrix d1 = c.getDelayedHitQLen();
        Matrix dfull = c.getDelayedHitQLenFull();
        for (int i = 0; i < d1.getNumCols(); i++) {
            double phi = dfull.get(0, i) - d1.get(0, i);
            assertTrue(phi > 0 && phi < 1,
                    "phi must be a probability, got " + phi + " for item " + (i + 1));
        }
    }

    @Test
    public void testHitDelayedMissSplitSumsToOne() {
        // Delayed hits depart in the hit class; the exact transition reward splits the
        // hit-class rate so that the three reported fractions partition the arrivals.
        Cache c = solve(2);
        double h = c.getHitRatio().get(0);
        double d = c.getDelayedHitRatio().get(0);
        double m = c.getMissRatio().get(0);
        assertEquals(0.467000, h, TOL);
        assertEquals(0.061583, d, TOL);
        assertEquals(0.471417, m, TOL);
        assertEquals(1.0, h + d + m, 1e-9);
    }

    @Test
    public void testDelayedHitMassIncreasesWithCutoff() {
        // maxPending = cutoff - 1 truncates block B, so a coarser cutoff can only
        // UNDERSTATE the delayed-hit mass and the delayed-hit queue length.
        Cache c2 = solve(2);
        Cache c3 = solve(3);
        double d2 = c2.getDelayedHitRatio().get(0);
        double d3 = c3.getDelayedHitRatio().get(0);
        assertTrue(d3 > d2, "delayed-hit mass must grow with the cutoff: " + d2 + " -> " + d3);
        assertTrue(c3.getDelayedHitQLen().get(0, 0) > c2.getDelayedHitQLen().get(0, 0),
                "d1 must grow with the cutoff");
    }

    @Test
    public void testConvergesTowardExactNormalizingConstantSolver() {
        // SolverNC solves this retrieval cache exactly (product form), so the truncated
        // CTMC must approach it from below as the cutoff grows.
        Network ncModel = model();
        new SolverNC(ncModel).getAvgTable();
        Cache ncCache = cacheOf(ncModel);
        double ncHitPlusDelayed = ncCache.getHitRatio().get(0) + ncCache.getDelayedHitRatio().get(0);

        double s2 = solve(2).getHitRatio().get(0) + solve(2).getDelayedHitRatio().get(0);
        double s3 = solve(3).getHitRatio().get(0) + solve(3).getDelayedHitRatio().get(0);
        assertTrue(s2 < s3, "cutoff 2 must undercount relative to cutoff 3");
        assertTrue(s3 < ncHitPlusDelayed + TOL,
                "truncated CTMC must not exceed the exact value: " + s3 + " vs " + ncHitPlusDelayed);
        assertTrue(ncHitPlusDelayed - s3 < ncHitPlusDelayed - s2,
                "cutoff 3 must be closer to the exact value than cutoff 2");
    }
}
