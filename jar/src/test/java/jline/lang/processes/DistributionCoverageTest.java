package jline.lang.processes;

import org.junit.jupiter.api.Test;
import jline.util.matrix.Matrix;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Coverage tests for distribution classes with 0% coverage.
 * Each test exercises the primary functionality of one distribution class.
 */
public class DistributionCoverageTest {

    private static final double TOL = 1e-6;

    // ========== EmpiricalCDF ==========

    @Test
    void testEmpiricalCDF() {
        // Create CDF data: F(x) values and corresponding x values
        Matrix cdfData = new Matrix(4, 1);
        cdfData.set(0, 0, 0.0);
        cdfData.set(1, 0, 0.25);
        cdfData.set(2, 0, 0.75);
        cdfData.set(3, 0, 1.0);

        Matrix xData = new Matrix(4, 1);
        xData.set(0, 0, 0.0);
        xData.set(1, 0, 1.0);
        xData.set(2, 0, 2.0);
        xData.set(3, 0, 3.0);

        EmpiricalCDF ecdf = new EmpiricalCDF(cdfData, xData);

        assertNotNull(ecdf);
        assertTrue(ecdf.getMean() > 0);
        assertNotNull(ecdf.getData());

        // Test CDF evaluation
        double cdfVal = ecdf.evalCDF(1.5);
        assertTrue(cdfVal >= 0 && cdfVal <= 1);

        // Test sampling
        double[] samples = ecdf.sample(10, new Random(42));
        assertEquals(10, samples.length);
    }

    // ========== MarkovProcess ==========

    @Test
    void testMarkovProcess() {
        // Create a simple 2-state CTMC generator matrix
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -1.0);
        Q.set(0, 1, 1.0);
        Q.set(1, 0, 2.0);
        Q.set(1, 1, -2.0);

        MarkovProcess mp = new MarkovProcess(Q);

        assertNotNull(mp);
        assertEquals(2, mp.getGenerator().getNumRows());
        assertNotNull(mp.getGenerator());
        assertTrue(mp.isFinite());
    }

    // ========== MarkedMMPP ==========

    @Test
    void testMarkedMMPP() {
        // Create MMPP parameters: D0 (transitions without arrivals), D1 (arrivals - must be diagonal)
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -2.0);
        D0.set(0, 1, 1.0);
        D0.set(1, 0, 1.0);
        D0.set(1, 1, -2.0);

        // D1 must be diagonal for MMPP
        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, 1.0);
        D1.set(1, 1, 1.0);

        // Create MatrixCell with D0 and D1
        jline.util.matrix.MatrixCell mmap = new jline.util.matrix.MatrixCell(2);
        mmap.set(0, D0);
        mmap.set(1, D1);

        MarkedMMPP mmpp = new MarkedMMPP(mmap);

        assertNotNull(mmpp);
        assertTrue(mmpp.getMean() > 0);
    }

    // ========== Prior ==========

    @Test
    void testPrior() {
        // Create a prior with two alternatives: Exp(1) and Exp(2)
        java.util.List<Distribution> dists = java.util.Arrays.asList(
            new Exp(1.0),
            new Exp(2.0)
        );
        double[] probs = {0.6, 0.4};

        Prior prior = new Prior(dists, probs);

        assertNotNull(prior);
        assertTrue(prior.getMean() > 0);
        assertEquals(2, prior.getNumAlternatives());
        assertEquals(0.6, prior.getProbability(0), TOL);
    }

    // ========== Cox2 ==========

    @Test
    void testCox2() {
        // Two-phase Coxian with rates mu1=2, mu2=1, phi1=0.5
        Cox2 cox = new Cox2(2.0, 1.0, 0.5);

        assertNotNull(cox);
        assertTrue(cox.getMean() > 0);
        assertTrue(cox.getSCV() > 0);

        // Test fitting
        Cox2 fitted = Cox2.fitMeanAndSCV(1.0, 2.0);
        assertNotNull(fitted);
        assertEquals(1.0, fitted.getMean(), 0.1);
    }

    // ========== Pareto ==========

    @Test
    void testPareto() {
        // Pareto with shape=3, scale=1
        Pareto pareto = new Pareto(3.0, 1.0);

        assertNotNull(pareto);
        assertEquals(1.5, pareto.getMean(), TOL); // mean = alpha*k/(alpha-1) = 3*1/2 = 1.5
        assertTrue(pareto.getSCV() > 0);
        assertTrue(pareto.getVar() > 0);

        // Test CDF
        double cdf = pareto.evalCDF(2.0);
        assertTrue(cdf > 0 && cdf < 1);

        // Test sampling
        double[] samples = pareto.sample(10, new Random(42));
        assertEquals(10, samples.length);
    }

    // ========== Weibull ==========

    @Test
    void testWeibull() {
        // Weibull with shape=2, scale=1
        Weibull weibull = new Weibull(2.0, 1.0);

        assertNotNull(weibull);
        assertTrue(weibull.getMean() > 0);
        assertTrue(weibull.getSCV() > 0);

        // Test CDF
        double cdf = weibull.evalCDF(1.0);
        assertTrue(cdf > 0 && cdf < 1);
    }

    // ========== Lognormal ==========

    @Test
    void testLognormal() {
        // Lognormal with mu=0, sigma=1
        Lognormal lognorm = new Lognormal(0.0, 1.0);

        assertNotNull(lognorm);
        assertTrue(lognorm.getMean() > 0);
        assertTrue(lognorm.getSCV() > 0);

        // Test CDF
        double cdf = lognorm.evalCDF(1.0);
        assertTrue(cdf > 0 && cdf < 1);
    }

    // ========== MarkovChain ==========

    @Test
    void testMarkovChain() {
        // Create a simple 2-state transition matrix
        Matrix P = new Matrix(2, 2);
        P.set(0, 0, 0.7);
        P.set(0, 1, 0.3);
        P.set(1, 0, 0.4);
        P.set(1, 1, 0.6);

        MarkovChain mc = new MarkovChain(P);

        assertNotNull(mc);
        assertEquals(2, mc.getTransMat().getNumRows());
        assertNotNull(mc.getTransMat());
        assertTrue(mc.isFinite());
    }

    // ========== Gamma ==========

    @Test
    void testGamma() {
        // Gamma with shape=2, scale=1
        Gamma gamma = new Gamma(2.0, 1.0);

        assertNotNull(gamma);
        assertEquals(2.0, gamma.getMean(), TOL); // mean = shape * scale
        assertEquals(0.5, gamma.getSCV(), TOL);  // SCV = 1/shape
        assertTrue(gamma.getVar() > 0);

        // Test CDF
        double cdf = gamma.evalCDF(2.0);
        assertTrue(cdf > 0 && cdf < 1);

        // Test sampling
        double[] samples = gamma.sample(10);
        assertEquals(10, samples.length);
    }

    // ========== Normal ==========

    @Test
    void testNormal() {
        // Normal with mu=0, sigma=1
        Normal normal = new Normal(0.0, 1.0);

        assertNotNull(normal);
        assertEquals(0.0, normal.getMean(), TOL);
        assertEquals(1.0, normal.getVar(), TOL);
        assertEquals(0.0, normal.getSkewness(), TOL);

        // Test CDF at mean
        assertEquals(0.5, normal.evalCDF(0.0), TOL);

        // Test sampling
        double[] samples = normal.sample(100, new Random(42));
        assertEquals(100, samples.length);
    }

    // ========== Binomial ==========

    @Test
    void testBinomial() {
        // Binomial with n=10, p=0.5
        Binomial binomial = new Binomial(10, 0.5);

        assertNotNull(binomial);
        assertEquals(5.0, binomial.getMean(), TOL); // mean = n*p
        assertTrue(binomial.getVar() > 0);

        // Test PMF
        double pmf = binomial.evalPMF(5);
        assertTrue(pmf > 0 && pmf < 1);
    }

    // ========== Uniform ==========

    @Test
    void testUniform() {
        // Uniform on [0, 10]
        Uniform uniform = new Uniform(0.0, 10.0);

        assertNotNull(uniform);
        assertEquals(5.0, uniform.getMean(), TOL);
        assertEquals(0.0, uniform.getSkewness(), TOL);

        // Test CDF
        assertEquals(0.5, uniform.evalCDF(5.0), TOL);
        assertEquals(0.0, uniform.evalCDF(-1.0), TOL);
        assertEquals(1.0, uniform.evalCDF(11.0), TOL);

        // Test sampling
        double[] samples = uniform.sample(10, new Random(42));
        assertEquals(10, samples.length);
        for (double s : samples) {
            assertTrue(s >= 0.0 && s <= 10.0);
        }
    }

    // ========== DiscreteUniform ==========

    @Test
    void testDiscreteUniform() {
        // DiscreteUniform on {1, 2, 3, 4, 5}
        DiscreteUniform dunif = new DiscreteUniform(1, 5);

        assertNotNull(dunif);
        assertEquals(3.0, dunif.getMean(), TOL); // mean = (min + max) / 2
        assertTrue(dunif.getVar() > 0);
        assertEquals(0.0, dunif.getSkewness(), TOL);

        // Test CDF
        double cdf = dunif.evalCDF(3);
        assertTrue(cdf > 0 && cdf <= 1);

        // Test sampling
        double[] samples = dunif.sample(10, new Random(42));
        assertEquals(10, samples.length);
    }

    // ========== Bernoulli ==========

    @Test
    void testBernoulli() {
        // Bernoulli with p=0.3
        Bernoulli bernoulli = new Bernoulli(0.3);

        assertNotNull(bernoulli);
        assertEquals(0.3, bernoulli.getMean(), TOL);

        // Test PMF
        assertEquals(0.7, bernoulli.evalPMF(0), TOL);
        assertEquals(0.3, bernoulli.evalPMF(1), TOL);
    }

    // ========== Poisson ==========

    @Test
    void testPoisson() {
        // Poisson with lambda=5
        Poisson poisson = new Poisson(5.0);

        assertNotNull(poisson);
        assertEquals(5.0, poisson.getMean(), TOL);
        assertEquals(5.0, poisson.getVar(), TOL);

        // Test PDF at mean
        double pdf = poisson.evalPDF(5);
        assertTrue(pdf > 0);

        // Test CDF
        double cdf = poisson.evalCDF(5.0);
        assertTrue(cdf > 0 && cdf < 1);
    }

    // ========== MarkedMarkovProcess ==========

    @Test
    void testMarkedMarkovProcess() {
        // Create a simple marked process with generator, event filters, and event list
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -1.0);
        Q.set(0, 1, 1.0);
        Q.set(1, 0, 2.0);
        Q.set(1, 1, -2.0);

        // Create event filter (marks transitions)
        Matrix eventFilter = new Matrix(2, 2);
        eventFilter.set(0, 1, 0.5);  // Some transitions are marked events

        jline.util.matrix.MatrixCell eventFilt = new jline.util.matrix.MatrixCell(1);
        eventFilt.set(0, eventFilter);

        // Empty event list for basic test
        java.util.List<java.util.Map<String, Object>> eventList = new java.util.ArrayList<>();

        MarkedMarkovProcess mmp = new MarkedMarkovProcess(Q, eventFilt, eventList);

        assertNotNull(mmp);
        assertEquals(2, mmp.getGenerator().getNumRows());
        assertNotNull(mmp.getEventFilt());
    }
}
