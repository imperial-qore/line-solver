package jline.solvers.fluid;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Disabled;
import jline.lang.processes.Exp;
import jline.lang.processes.Zipf;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

/**
 * The second-order moment closure on a cache-queueing model.
 *
 * <p>{@code minnormal} used to be refused on any model with a cache node. It is
 * now answered through the same decomposition as {@code rmf}, with the closure
 * in the network step instead of the matrix method, and it reports the second
 * moment of both layers: the queueing covariance and the linear noise
 * covariance of the cache item occupancy.</p>
 *
 * <p>Expected values are the MATLAB reference on the same model.</p>
 */
public class CacheMinNormalTest {

    private static final double TOL = 1e-2;
    private static final int NITEMS = 10;
    private static final int CAP = 3;
    private static final int N = 4;

    /** Think -> Cache -> Q1(PS) -> Think, Zipf(0.8) reads over 10 items, RANDOM(3). */
    private static Network cacheModel() {
        Network model = new Network("cacheqn");
        Delay think = new Delay(model, "Think");
        Cache cache = new Cache(model, "Cache", NITEMS, CAP, ReplacementStrategy.RR);
        Queue queue = new Queue(model, "Q1", SchedStrategy.PS);

        ClosedClass req = new ClosedClass(model, "Req", N, think);
        ClosedClass hit = new ClosedClass(model, "Hit", 0, think);
        ClosedClass miss = new ClosedClass(model, "Miss", 0, think);

        think.setService(req, Exp.fitMean(1.0));
        think.setService(hit, Exp.fitMean(0.1));
        think.setService(miss, Exp.fitMean(0.1));
        queue.setService(req, new Disabled());
        queue.setService(hit, Exp.fitMean(0.2));
        queue.setService(miss, Exp.fitMean(1.0));

        cache.setRead(req, new Zipf(0.8, NITEMS));
        cache.setHitClass(req, hit);
        cache.setMissClass(req, miss);

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(req, req, think, cache, 1.0);
        P.set(hit, hit, cache, queue, 1.0);
        P.set(miss, miss, cache, queue, 1.0);
        P.set(hit, req, queue, think, 1.0);
        P.set(miss, req, queue, think, 1.0);
        model.link(P);
        return model;
    }

    private static SolverFluid solver(Network model, String method) {
        SolverOptions o = SolverFluid.defaultOptions();
        o.method = method;
        o.verbose = VerboseLevel.SILENT;
        return new SolverFluid(model, o);
    }

    /**
     * The closure answers the cache model and conserves the population. MATLAB:
     * Think 1.3873, Q1 total 2.6127 under minnormal, against 1.4413 / 2.5587
     * under rmf.
     */
    @Test
    public void testClosureAnswersCacheModel() {
        Matrix q = solver(cacheModel(), "minnormal").getAvgQLen();
        double think = q.get(0, 0) + q.get(0, 1) + q.get(0, 2);
        double queue = q.get(1, 0) + q.get(1, 1) + q.get(1, 2);
        assertEquals(1.3873, think, TOL);
        assertEquals(2.6127, queue, TOL);
        assertEquals(N, think + queue, TOL, "the decomposition must conserve the population");
    }

    /** Only the queueing layer changes: the cache layer is the refined mean field either way. */
    @Test
    public void testCacheLayerIsShared() {
        Network rmfModel = cacheModel();
        SolverFluid rmf = solver(rmfModel, "rmf");
        rmf.getAvgQLen();
        Matrix hitRmf = ((Cache) rmfModel.getNodeByName("Cache")).getHitRatio();

        Network minModel = cacheModel();
        SolverFluid min = solver(minModel, "minnormal");
        min.getAvgQLen();
        Matrix hitMin = ((Cache) minModel.getNodeByName("Cache")).getHitRatio();

        assertEquals(hitRmf.get(0, 0), hitMin.get(0, 0), 1e-6,
                "the cache layer must be identical under both routes");
        assertEquals(0.3827, hitMin.get(0, 0), TOL);
    }

    /**
     * The moment report carries both layers, and the cache covariance respects
     * the capacity constraint: the number of uncached items is deterministic, so
     * the covariance of the miss indicators sums to zero.
     */
    @Test
    public void testCacheMomentReport() {
        Network model = cacheModel();
        SolverFluid s = solver(model, "minnormal");
        s.getAvgQLen();
        FluidResult fr = (FluidResult) s.result;

        assertNotNull(fr.momentQVar, "no queueing covariance reported");
        assertNotNull(fr.momentCacheSigma, "no cache covariance reported");
        assertEquals(1, fr.momentCacheSigma.length);

        Matrix sigma = fr.momentCacheSigma[0];
        Matrix pi0 = fr.momentCachePi0[0];
        assertEquals(NITEMS, pi0.getNumRows());

        double total = 0;
        for (int i = 0; i < NITEMS; i++) {
            for (int j = 0; j < NITEMS; j++) {
                total += sigma.get(i, j);
            }
        }
        assertEquals(0.0, total, 1e-8, "the occupancy covariance violates the capacity constraint");

        // MATLAB: pi0(1) = 0.4203, its variance 0.2135, miss-probability variance 0.01055
        assertEquals(0.4203, pi0.get(0, 0), TOL);
        assertEquals(0.2135, sigma.get(0, 0), TOL);
        assertEquals(0.01055, fr.momentCacheMissProbVar.get(0, 0), 1e-3);

        // the miss indicator is Bernoulli up to the O(1/n) error of the closure
        for (int i = 0; i < NITEMS; i++) {
            double bern = pi0.get(i, 0) * (1 - pi0.get(i, 0));
            assertTrue(Math.abs(sigma.get(i, i) - bern) / bern < 0.35,
                    "item " + i + ": variance " + sigma.get(i, i) + " against Bernoulli " + bern);
        }
    }
}
