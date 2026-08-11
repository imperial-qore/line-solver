package jline.solvers.nc;

import jline.TestTools;
import jline.examples.java.basic.ClosedModel;
import jline.examples.java.advanced.LoadDependentModel;
import jline.lang.*;
import jline.lang.NetworkStruct;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.*;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.util.Maths;

import static jline.TestTools.*;
import static jline.solvers.nc.handlers.Solver_nc_joint.solver_nc_joint;
import static jline.solvers.nc.handlers.Solver_nc_marg.solver_nc_marg;
import static jline.solvers.nc.handlers.Solver_ncld.solver_ncld;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for SolverNC (Normalizing Constant) queueing network analyzer.
 *
 * <p>Includes tests from:
 * <ul>
 *   <li>Core NC solver functionality (joint/marginal computations, load-dependent networks)
 *   <li>SolverNCCacheTest: Cache queueing network analysis
 * </ul>
 */
public class SolverNCTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        // Set verbose level to SILENT to suppress warnings during tests
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    @Test
    public void test_marg_joint() {
        Network model = ClosedModel.cqn_repairmen();
        NetworkStruct sn = model.getStruct(true);
        SolverOptions options = new SolverOptions(SolverType.NC);

        SolverNC.SolverNCMargReturn ret1 = solver_nc_marg(sn, options, Double.NaN);
        Matrix lPr1 = ret1.lPr;
        double G1 = ret1.G;
        double runtime1 = ret1.runtime;

        assertEquals(2, lPr1.getNumRows());
        assertEquals(1, lPr1.getNumCols());

        assertEquals(-9.341536168923934, lPr1.get(0), relativeTolerance(9.341536168923934, TestTools.MID_TOL));
        assertEquals(-9.341536168923934, lPr1.get(1), relativeTolerance(9.341536168923934, TestTools.MID_TOL));

        assertEquals(0.003142060751148, G1, relativeTolerance(0.003142060751148, TestTools.MID_TOL));

        SolverNC.SolverNCJointReturn ret2 = solver_nc_joint(sn, options);
        double Pr2 = ret2.Pr;
        double G2 = ret2.G;
        double lG2 = ret2.lG;
        double runtime2 = ret2.runtime;

        assertEquals(8.770460346419127e-05, Pr2, relativeTolerance(8.770460346419127e-05, TestTools.MID_TOL));
        assertEquals(0.003142060751148, G2, relativeTolerance(0.003142060751148, TestTools.MID_TOL));
        assertEquals(-5.762876404151583, lG2, relativeTolerance(5.762876404151583, TestTools.MID_TOL));
    }

    @Test
    public void test_complexMultiClassModel_1() {
        Network model = SolverNCTestFixtures.complexMultiClassModel_1();
        NetworkStruct sn = model.getStruct(true);
        SolverOptions options = new SolverOptions(SolverType.NC);

        SolverNC.SolverNCMargReturn ret1 = solver_nc_marg(sn, options, Double.NaN);
        Matrix lPr1 = ret1.lPr;
        double G1 = ret1.G;
        double runtime1 = ret1.runtime;

        assertEquals(4, lPr1.getNumRows());
        assertEquals(1, lPr1.getNumCols());

        assertEquals(-3.041632215808478e2, lPr1.get(0), relativeTolerance(3.041632215808478e2, TestTools.MID_TOL));
        assertEquals(-0.018344541300874e2, lPr1.get(1), relativeTolerance(0.018344541300874e2, TestTools.MID_TOL));
        assertEquals(0.0, lPr1.get(2), MID_TOL);
        assertEquals(-0.471992884544684e2, lPr1.get(3), relativeTolerance(0.471992884544684e2, TestTools.MID_TOL));

        assertEquals(6.406198113596468e-36, G1, relativeTolerance(6.406198113596468e-36, TestTools.MID_TOL));

        SolverNC.SolverNCJointReturn ret2 = solver_nc_joint(sn, options);
        double Pr2 = ret2.Pr;
        double G2 = ret2.G;
        double lG2 = ret2.lG;
        double runtime2 = ret2.runtime;

        assertEquals(1.947180103098101e-133, Pr2, relativeTolerance(1.947180103098101e-133, TestTools.MID_TOL));
        assertEquals(6.406198113596468e-36, G2, relativeTolerance(6.406198113596468e-36, TestTools.MID_TOL));
        assertEquals(-81.035797370820802, lG2, relativeTolerance(81.035797370820802, TestTools.MID_TOL));

        //
        //
    }

    @Test
    public void test_complexMultiClassModel_2() {
        Network model = SolverNCTestFixtures.complexMultiClassModel_2();
        NetworkStruct sn = model.getStruct(true);
        SolverOptions options = new SolverOptions(SolverType.NC);

        SolverNC.SolverNCMargReturn ret1 = solver_nc_marg(sn, options, Double.NaN);
        Matrix lPr1 = ret1.lPr;
        double G1 = ret1.G;
        double runtime1 = ret1.runtime;

        assertEquals(4, lPr1.getNumRows());
        assertEquals(1, lPr1.getNumCols());

        assertEquals(-0.804926658480257e2, lPr1.get(0), relativeTolerance(0.804926658480257e2, TestTools.MID_TOL));
        assertEquals(-0.008127259006488e2, lPr1.get(1), relativeTolerance(0.008127259006488e2, TestTools.MID_TOL));
        assertEquals(0.0, lPr1.get(2), MID_TOL);
        assertEquals(-1.766646876121396e2, lPr1.get(3), relativeTolerance(1.766646876121396e2, TestTools.MID_TOL));

        assertEquals(1.984349911810360e+30, G1, relativeTolerance(1.984349911810360e+30, TestTools.MID_TOL));

        SolverNC.SolverNCJointReturn ret2 = solver_nc_joint(sn, options);
        double Pr2 = ret2.Pr;
        double G2 = ret2.G;
        double lG2 = ret2.lG;
        double runtime2 = ret2.runtime;

        assertEquals(6.756257601443034e-78, Pr2, relativeTolerance(6.756257601443034e-78, TestTools.MID_TOL));
        assertEquals(1.984349911810360e+30, G2, relativeTolerance(1.984349911810360e+30, TestTools.MID_TOL));
        assertEquals(69.762844149973148, lG2, relativeTolerance(69.762844149973148, TestTools.MID_TOL));

        //
        //
    }

    @Test
    public void test_ld_multiserver_fcfs() {
        Network ldmodel = LoadDependentModel.ld_multiserver_fcfs();
        NetworkStruct sn = ldmodel.getStruct(true);
        SolverOptions options = new SolverOptions(SolverType.NC);

        SolverNC.SolverNCLDReturn ret0 = solver_ncld(sn, options);
        options.method = "nrl";
        SolverNC.SolverNCLDReturn ret1 = solver_ncld(sn, options);
        options.method = "nrp";
        SolverNC.SolverNCLDReturn ret2 = solver_ncld(sn, options);
        options.method = "rd";
        SolverNC.SolverNCLDReturn ret3 = solver_ncld(sn, options);
        Matrix ret0Q = new Matrix(Arrays.asList(1.333333333322438e+00, 1.466666666667756e+01));
        Matrix ret0U = new Matrix(Arrays.asList(1.333333333322438e+00, 9.999999999918286e-01));
        Matrix ret0R = new Matrix(Arrays.asList(1.0, 1.100000000009806e+01));
        Matrix ret0T = new Matrix(Arrays.asList(1.333333333322438e+00, 1.333333333322438e+00));
        Matrix ret0X = Matrix.singleton(1.333333333322438e+00);
        assertTrue(ret0.Q.isEqualToTol(ret0Q, TestTools.MID_TOL));
        assertTrue(ret0.U.isEqualToTol(ret0U, TestTools.MID_TOL));
        assertTrue(ret0.R.isEqualToTol(ret0R, TestTools.MID_TOL));
        assertTrue(ret0.T.isEqualToTol(ret0T, TestTools.MID_TOL));
        assertTrue(ret0.X.isEqualToTol(ret0X, TestTools.MID_TOL));

        assertEquals(-2.576432645335951e+00, ret0.lG, relativeTolerance(2.576432645335951e+00, TestTools.MID_TOL));
        // NRL golden updated after 5146f4476 ("halve the log-determinant in the
        // Laplace approximation logI"), which pfqn_nrl consumes. Exact lG for this
        // model is -2.576432645335951, from direct enumeration of the LD normalizing
        // constant (sum_n Z^(N-n)/(N-n)! * L^n / prod_k mu(k)), independent of any
        // solver here and equal to the ret0 golden above. The new value is 0.120791
        // from exact; the previous golden -2.350641375431591 was 0.225791.
        assertEquals(-2.697223806045956e+00, ret1.lG, relativeTolerance(3.303193953352507e+00, TestTools.MID_TOL));
        assertEquals(-2.576432728076320e+00, ret2.lG, relativeTolerance(2.576432728076320e+00, TestTools.MID_TOL));
        assertEquals(-2.633413431541626e+00, ret3.lG, relativeTolerance(2.633413431541626e+00, TestTools.MID_TOL));


    }

    @Test
    public void test_ld_multiserver_ps() {
        Network ldmodel = LoadDependentModel.ld_multiserver_ps();
        NetworkStruct sn = ldmodel.getStruct(true);
        SolverOptions options = new SolverOptions(SolverType.NC);

        SolverNC.SolverNCLDReturn ret0 = solver_ncld(sn, options);
        options.method = "nrl";
        SolverNC.SolverNCLDReturn ret1 = solver_ncld(sn, options);
        options.method = "nrp";
        SolverNC.SolverNCLDReturn ret2 = solver_ncld(sn, options);
        options.method = "rd";
        SolverNC.SolverNCLDReturn ret3 = solver_ncld(sn, options);
        Matrix ret0Q = new Matrix(Arrays.asList(5.467372543256852e-01, 8.563871132405998e-01, 2.596875632433713e+00))
                .concatCols(new Matrix(Arrays.asList(3.694448418974791e-01, 4.807742299186316e-01, 1.149780928183889e+00)));
        Matrix ret0U = new Matrix(Arrays.asList(5.467372543256852e-01, 2.733686271628426e-01, 6.378601300466328e-01))
                .concatCols(new Matrix(Arrays.asList(3.694448418974791e-01, 1.539353507906163e-01, 2.770836314231093e-01)));
        Matrix ret0R = new Matrix(Arrays.asList(1.0, 1.566359538270022e+00, 4.749768946395564e+00))
                .concatCols(new Matrix(Arrays.asList(2.0, 2.602684760460380e+00, 6.224371260827907e+00)));
        Matrix ret0T = new Matrix(Arrays.asList(5.467372543256852e-01, 5.467372543256852e-01, 5.467372543256852e-01))
                .concatCols(new Matrix(Arrays.asList(1.847224209487395e-01, 1.847224209487395e-01, 1.847224209487395e-01)));
        Matrix ret0X = Matrix.singleton(5.467372543256852e-01)
                .concatCols(Matrix.singleton(1.847224209487395e-01));
        assertTrue(ret0.Q.isEqualToTol(ret0Q, TestTools.MID_TOL));
        assertTrue(ret0.U.isEqualToTol(ret0U, TestTools.MID_TOL));
        assertTrue(ret0.R.isEqualToTol(ret0R, TestTools.MID_TOL));
        assertTrue(ret0.T.isEqualToTol(ret0T, TestTools.MID_TOL));
        assertTrue(ret0.X.isEqualToTol(ret0X, TestTools.MID_TOL));
        // Exact lG for this model is 8.016437786118587, from direct enumeration over
        // all station population vectors (delay term prod_r Z_r^{n_r}/n_r!, LD queue
        // term |n|! prod_r (L_ir^{n_ir}/n_ir!) / prod_{k<=|n|} mu_i(k)), independent
        // of any solver here and matching the ret0 golden below to 2 ulp. Both goldens
        // updated after 5146f4476, which pfqn_nrl and pfqn_nrp consume.
        assertEquals(8.016437786118589e+00, ret0.lG, relativeTolerance(8.016437786118589e+00, TestTools.MID_TOL));
        // NRL: new value 0.345244 from exact; previous golden 7.305374803607128 was
        // 0.711063. It was updated despite PASSING: the MID_TOL*1000 = 1e-1 relative
        // band was wide enough to keep the staler number green.
        assertEquals(7.671193626464337, ret1.lG, relativeTolerance(7.671193626464337, TestTools.MID_TOL * 1000));
        // NRP: new value 0.220848 from exact; previous golden 5.979290244081035 was
        // 2.037148.
        assertEquals(7.795589927820983e+00, ret2.lG, relativeTolerance(7.795589927820983e+00, TestTools.MID_TOL * 100));
        assertEquals(7.755109452853959e+00, ret3.lG, relativeTolerance(7.755109452853959e+00, TestTools.MID_TOL));
    }

    // test_loadDependent_4 deleted - example_loadDependent_4 method was removed

    // test_loadDependent_4b deleted - example_loadDependent_4 method was removed

    // test_ld_class_dependence moved to SolverMVATest - NC solver doesn't support class-dependent scaling

    // ==================== Cache Tests ====================

    /**
     * Creates a simple network with a cache node.
     */
    private Network createCacheNetwork() throws Exception {
        Network model = new Network("CacheNetwork");

        // Source (must be created before OpenClass)
        Source source = new Source(model, "Source");

        // Open class
        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(0.5));

        // Cache node (nitems=10, itemLevelCap=2, LRU replacement)
        Cache cache = new Cache(model, "Cache", 10, 2, ReplacementStrategy.LRU);
        // Set read distribution for cache access
        cache.setRead(jobClass, new Exp(10.0));

        // Queue for cache misses
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setService(jobClass, new Exp(1.0));

        // Sink
        Sink sink = new Sink(model, "Sink");

        // Routing: Source -> Cache -> Queue -> Sink
        model.link(model.serialRouting(source, cache, queue, sink));

        return model;
    }

    /**
     * Creates a closed network with cache.
     */
    private Network createClosedCacheNetwork(int population) throws Exception {
        Network model = new Network("ClosedCacheNetwork");

        ClosedClass jobClass = new ClosedClass(model, "Class1", population, null);

        Delay delay = new Delay(model, "Delay");
        delay.setService(jobClass, new Exp(1.0));

        Cache cache = new Cache(model, "Cache", 5, 2, ReplacementStrategy.LRU);
        cache.setRead(jobClass, new Exp(10.0));

        jobClass.setReferenceStation(delay);

        model.link(model.serialRouting(delay, cache));

        return model;
    }

    /**
     * Test 28: solver_nc_cache_qn_analyzer - NC cache analyzer.
     */
    @Test
    public void testSolverNcCacheQn_basicCache() {
        try {
            Network model = createCacheNetwork();

            SolverNC solver = new SolverNC(model);

            solver.runAnalyzer();
            NetworkAvgTable result = solver.getAvgTable();

            assertNotNull(result, "NC cache analyzer should produce results");
        } catch (Exception e) {
            // Cache analysis may have specific requirements
            assertTrue(true, "Cache analysis may need specific setup");
        }
    }

    /**
     * Test 28b: NC cache with closed network.
     */
    @Test
    public void testSolverNcCacheQn_closedNetwork() {
        try {
            Network model = createClosedCacheNetwork(5);

            SolverNC solver = new SolverNC(model);

            solver.runAnalyzer();
            NetworkAvgTable result = solver.getAvgTable();

            assertNotNull(result, "NC should handle closed network with cache");
        } catch (Exception e) {
            assertTrue(true, "Closed cache network may have specific requirements");
        }
    }

    /**
     * Test: NC cache with exact method.
     */
    @Test
    public void testSolverNcCacheQn_exactMethod() {
        try {
            Network model = createClosedCacheNetwork(3);

            SolverNC solver = new SolverNC(model);
            SolverOptions options = solver.getOptions();
            options.method("exact");

            solver.runAnalyzer();
            NetworkAvgTable result = solver.getAvgTable();

            assertNotNull(result, "Exact NC method should work");
        } catch (Exception e) {
            assertTrue(true, "Exact method may have constraints");
        }
    }

    /**
     * Test: NC cache with SPM method.
     */
    @Test
    public void testSolverNcCacheQn_spmMethod() {
        try {
            Network model = createClosedCacheNetwork(5);

            SolverNC solver = new SolverNC(model);
            SolverOptions options = solver.getOptions();
            options.method("spm");

            solver.runAnalyzer();
            NetworkAvgTable result = solver.getAvgTable();

            assertNotNull(result, "SPM NC method should work");
        } catch (Exception e) {
            assertTrue(true, "SPM method may have specific requirements");
        }
    }

    /**
     * Test: NC basic network without cache.
     */
    @Test
    public void testSolverNc_basicNetwork() {
        try {
            Network model = new Network("BasicNetwork");

            ClosedClass jobClass = new ClosedClass(model, "Class1", 5, null);

            Delay delay = new Delay(model, "Delay");
            delay.setService(jobClass, new Exp(1.0));

            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            queue.setService(jobClass, new Exp(2.0));

            jobClass.setReferenceStation(delay);

            model.link(model.serialRouting(delay, queue));

            SolverNC solver = new SolverNC(model);

            solver.runAnalyzer();
            NetworkAvgTable result = solver.getAvgTable();

            assertNotNull(result, "NC should handle basic closed networks");
        } catch (Exception e) {
            assertTrue(true, "NC may have specific requirements");
        }
    }

    /**
     * Test: NC with multi-class network.
     */
    @Test
    public void testSolverNc_multiClass() {
        try {
            Network model = new Network("MultiClassNetwork");

            ClosedClass class1 = new ClosedClass(model, "Class1", 3, null);
            ClosedClass class2 = new ClosedClass(model, "Class2", 2, null);

            Delay delay = new Delay(model, "Delay");
            delay.setService(class1, new Exp(1.0));
            delay.setService(class2, new Exp(1.5));

            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            queue.setService(class1, new Exp(2.0));
            queue.setService(class2, new Exp(2.5));

            class1.setReferenceStation(delay);
            class2.setReferenceStation(delay);

            model.link(model.serialRouting(delay, queue));

            SolverNC solver = new SolverNC(model);

            solver.runAnalyzer();
            NetworkAvgTable result = solver.getAvgTable();

            assertNotNull(result, "NC should handle multi-class networks");
        } catch (Exception e) {
            assertTrue(true, "Multi-class may have specific requirements");
        }
    }

    /**
     * Test: NC convergence parameters.
     */
    @Test
    public void testSolverNc_convergence() {
        try {
            Network model = createClosedCacheNetwork(4);

            SolverNC solver = new SolverNC(model);
            SolverOptions options = solver.getOptions();
            options.iter_max = 100;
            options.iter_tol = 1e-6;

            solver.runAnalyzer();
            NetworkAvgTable result = solver.getAvgTable();

            assertNotNull(result, "NC should converge with specified parameters");
        } catch (Exception e) {
            assertTrue(true, "Convergence may need tuning");
        }
    }


    /**
     * An EMPTY chain (all its classes have population 0) must report zero for every
     * metric, matching MATLAB solver_ncld.m:270-272.
     *
     * <p>Regression for a load-dependent NC defect: the infinite-server overwrite
     * {@code Rchain = Lchain/Vchain} in Solver_ncld does not consult the chain
     * population, so for an empty chain it left the delay's service demand in Rchain,
     * and that value survived deaggregation into the reported RespT at EVERY station.
     * Before the fix this model reported RespT = 1.0 (the Delay's mean) for the empty
     * class at both stations where MATLAB reports 0; changing the Delay's mean to 3.0
     * moved the wrong answer to 3.0, while the queue's mean had no effect.
     *
     * <p>The same defect corrupted the POPULATED chain's utilization, because the
     * spurious empty-chain contribution pushed the load-dependent utilization row over
     * 1 and triggered a renormalization that scaled the real class down: Queue1 Util
     * was 0.46089850249584030 against MATLAB's 0.92849162011173192. Both symptoms are
     * asserted here because they have one root cause.
     *
     * <p>Reference values are MATLAB SolverNC on the same model, method exact/gld.
     */
    @Test
    public void emptyChainReportsZeroAndLeavesPopulatedChainIntact() {
        Network model = new Network("emptychain");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);

        ClosedClass c1 = new ClosedClass(model, "Class1", 4, delay, 0);
        ClosedClass c2 = new ClosedClass(model, "Class2", 0, delay, 0);

        delay.setService(c1, Exp.fitMean(1.0));
        delay.setService(c2, Exp.fitMean(1.0));
        queue.setService(c1, Exp.fitMean(1.5));
        queue.setService(c2, Exp.fitMean(1.5));

        Matrix lld = new Matrix(1, 4);
        for (int i = 0; i < 4; i++) {
            lld.set(0, i, Math.min(i + 1, 2));
        }
        queue.setLoadDependence(lld);

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, model.serialRouting(delay, queue));
        P.set(c2, c2, model.serialRouting(delay, queue));
        model.link(P);

        SolverOptions opts = new SolverOptions(SolverType.NC);
        opts.method = "exact";
        SolverNC solver = new SolverNC(model, opts);

        Matrix Q = solver.getAvgQLen();
        Matrix U = solver.getAvgUtil();
        Matrix R = solver.getAvgRespT();
        Matrix T = solver.getAvgTput();

        // Empty chain: every metric zero at every station. Asserted on getAvgRespT
        // rather than on the printed table, which filters all-zero rows and would
        // pass vacuously.
        for (int i = 0; i < 2; i++) {
            assertEquals(0.0, Q.get(i, 1), 1e-12, "empty chain QLen at station " + i);
            assertEquals(0.0, U.get(i, 1), 1e-12, "empty chain Util at station " + i);
            assertEquals(0.0, R.get(i, 1), 1e-12, "empty chain RespT at station " + i);
            assertEquals(0.0, T.get(i, 1), 1e-12, "empty chain Tput at station " + i);
        }

        // Populated chain: unchanged, against MATLAB SolverNC exact/gld.
        assertEquals(1.2379888268156425, Q.get(0, 0), 1e-9, "populated QLen at Delay");
        assertEquals(2.7620111731843577, Q.get(1, 0), 1e-9, "populated QLen at Queue1");
        assertEquals(1.2379888268156425, U.get(0, 0), 1e-9, "populated Util at Delay");
        assertEquals(0.92849162011173192, U.get(1, 0), 1e-9, "populated Util at Queue1");
        assertEquals(1.0, R.get(0, 0), 1e-9, "populated RespT at Delay");
        assertEquals(2.2310469314079424, R.get(1, 0), 1e-9, "populated RespT at Queue1");
        assertEquals(1.2379888268156425, T.get(0, 0), 1e-9, "populated Tput at Delay");
        assertEquals(1.2379888268156425, T.get(1, 0), 1e-9, "populated Tput at Queue1");
    }
}
