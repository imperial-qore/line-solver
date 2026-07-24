package jline.solvers.ssa;

import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.DiscreteSampler;
import jline.lang.processes.Exp;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Source;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Queue;
import jline.lang.JobClass;
import jline.lang.ClosedClass;
import jline.lang.OpenClass;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * NRM cache support: the cache access (read -> hit/miss class switch with a
 * replacement-policy update) is simulated by the next-reaction engine at parity
 * with State.afterEventCache. The hit ratio must match the exact CTMC across
 * every replacement policy. Under a Zipf reference the policies DIVERGE, which is
 * the binding control -- an exp-collapsed or serial-fallback implementation would
 * return a policy-independent or wrong ratio.
 */
public class SolverSSANrmCacheTest {

    private static final int SAMPLES = 300000;
    private static final int SEED = 23000;
    private static final double RTOL = 0.03;

    private static Network build(int n, int m, ReplacementStrategy rep, double[] pop) {
        Network model = new Network("cache");
        Delay delay = new Delay(model, "Delay");
        Cache cache = new Cache(model, "Cache", n, m, rep);
        ClosedClass jobClass = new ClosedClass(model, "JobClass", 1, delay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, delay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, delay, 0);
        delay.setService(jobClass, new Exp(1.0));
        cache.setRead(jobClass, new DiscreteSampler(new Matrix(new double[][]{pop})));
        cache.setHitClass(jobClass, hitClass);
        cache.setMissClass(jobClass, missClass);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, delay, cache, 1.0);
        P.set(hitClass, jobClass, cache, delay, 1.0);
        P.set(missClass, jobClass, cache, delay, 1.0);
        model.link(P);
        return model;
    }

    private static double[] uniform(int n) {
        double[] p = new double[n];
        for (int i = 0; i < n; i++) p[i] = 1.0 / n;
        return p;
    }

    private static double[] zipf(int n, double a) {
        double[] p = new double[n];
        double s = 0.0;
        for (int i = 0; i < n; i++) { p[i] = Math.pow(i + 1, -a); s += p[i]; }
        for (int i = 0; i < n; i++) p[i] /= s;
        return p;
    }

    private static double nrmHit(Network model, Cache cache) {
        SolverSSA solver = new SolverSSA(model);
        solver.options.method = "nrm";
        solver.options.samples = SAMPLES;
        solver.options.seed = SEED;
        solver.getAvgNodeTable();
        return cache.getHitRatio().get(0);
    }

    private static double ctmcHit(Network model, Cache cache) {
        SolverCTMC solver = new SolverCTMC(model);
        solver.getAvgNodeTable();
        return cache.getHitRatio().get(0);
    }

    private void checkPolicy(ReplacementStrategy rep, int n, int m, double[] pop) {
        Network mn = build(n, m, rep, pop);
        Cache cn = (Cache) mn.getNodeByName("Cache");
        double hn = nrmHit(mn, cn);
        Network mc = build(n, m, rep, pop);
        Cache cc = (Cache) mc.getNodeByName("Cache");
        double hc = ctmcHit(mc, cc);
        assertEquals(hc, hn, RTOL * Math.max(hc, 1e-9),
                rep + " hit ratio NRM=" + hn + " CTMC=" + hc);
    }

    @Test public void uniformLRU()   { checkPolicy(ReplacementStrategy.LRU, 5, 2, uniform(5)); }
    @Test public void uniformFIFO()  { checkPolicy(ReplacementStrategy.FIFO, 5, 2, uniform(5)); }
    @Test public void uniformRR()    { checkPolicy(ReplacementStrategy.RR, 5, 2, uniform(5)); }

    @Test public void zipfLRU()   { checkPolicy(ReplacementStrategy.LRU, 6, 2, zipf(6, 1.2)); }
    @Test public void zipfFIFO()  { checkPolicy(ReplacementStrategy.FIFO, 6, 2, zipf(6, 1.2)); }
    @Test public void zipfRR()    { checkPolicy(ReplacementStrategy.RR, 6, 2, zipf(6, 1.2)); }
    @Test public void zipfSFIFO() { checkPolicy(ReplacementStrategy.SFIFO, 6, 2, zipf(6, 1.2)); }
    @Test public void zipfHLRU()  { checkPolicy(ReplacementStrategy.HLRU, 6, 2, zipf(6, 1.2)); }
    @Test public void zipfQLRU()  { checkPolicy(ReplacementStrategy.QLRU, 6, 2, zipf(6, 1.2)); }

    /** Binding control: under Zipf, LRU must cache the popular items better than FIFO. */
    @Test public void zipfPolicyDivergence() {
        Network ml = build(6, 2, ReplacementStrategy.LRU, zipf(6, 1.2));
        double hl = nrmHit(ml, (Cache) ml.getNodeByName("Cache"));
        Network mf = build(6, 2, ReplacementStrategy.FIFO, zipf(6, 1.2));
        double hf = nrmHit(mf, (Cache) mf.getNodeByName("Cache"));
        assertTrue(hl > hf + 0.005, "LRU (" + hl + ") should exceed FIFO (" + hf + ") under Zipf");
    }

    // ------------------------------------------------------------------
    // Retrieval (delayed-hit) system. An open Source -> Cache -> Queue ->
    // Cache -> Sink loop: a miss for an item not yet being fetched begins a
    // retrieval (the job goes to the fetch queue and returns to complete the
    // miss), and a concurrent request for an item already being fetched is
    // absorbed as a delayed hit. The retrieval state space OOMs the CTMC, so the
    // oracle is the serial SSA engine (there is a known serial-vs-LDES divergence
    // on delayed-hit counting; the NRM targets the serial engine).
    // ------------------------------------------------------------------
    private static Network buildRetrieval(ReplacementStrategy rep) {
        Network model = new Network("DelayedHits");
        double[] acc = {0.6, 0.3, 0.1};
        Source source = new Source(model, "Source");
        Cache cache = new Cache(model, "Cache", 3, 1, rep);
        Queue queue = new Queue(model, "Queue", SchedStrategy.INF);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);
        source.setArrival(jobClass, new Exp(1.0));
        queue.setService(jobClass, new Exp(2.0));
        cache.setRead(jobClass, new DiscreteSampler(new Matrix(new double[][]{acc})));
        cache.setHitClass(jobClass, hitClass);
        cache.setMissClass(jobClass, missClass);
        cache.setRetrievalSystem(jobClass, missClass, new Queue[]{queue});
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, cache, 1.0);
        P.set(jobClass, jobClass, cache, queue, 1.0);
        P.set(jobClass, jobClass, queue, cache, 1.0);
        P.set(hitClass, hitClass, cache, sink, 1.0);
        P.set(missClass, missClass, cache, sink, 1.0);
        model.link(P);
        return model;
    }

    private static double serialHit(Network model, Cache cache) {
        SolverSSA solver = new SolverSSA(model);
        solver.options.method = "serial";
        solver.options.samples = SAMPLES;
        solver.options.seed = SEED;
        solver.getAvgNodeTable();
        return cache.getHitRatio().get(0);
    }

    @Test public void retrievalFifoMatchesSerial() {
        Network mn = buildRetrieval(ReplacementStrategy.FIFO);
        double hn = nrmHit(mn, (Cache) mn.getNodeByName("Cache"));
        Network ms = buildRetrieval(ReplacementStrategy.FIFO);
        double hs = serialHit(ms, (Cache) ms.getNodeByName("Cache"));
        assertEquals(hs, hn, 0.04 * Math.max(hs, 1e-9),
                "retrieval FIFO hit NRM=" + hn + " serial=" + hs);
    }

    @Test public void retrievalLruMatchesSerial() {
        Network mn = buildRetrieval(ReplacementStrategy.LRU);
        double hn = nrmHit(mn, (Cache) mn.getNodeByName("Cache"));
        Network ms = buildRetrieval(ReplacementStrategy.LRU);
        double hs = serialHit(ms, (Cache) ms.getNodeByName("Cache"));
        assertEquals(hs, hn, 0.04 * Math.max(hs, 1e-9),
                "retrieval LRU hit NRM=" + hn + " serial=" + hs);
    }
}
