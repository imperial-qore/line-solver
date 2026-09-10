/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.api.mam.Map_pdf;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.mam.SolverMAM;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static jline.lib.butools.ph.CheckMERepresentation.checkMERepresentation;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for the CME (Concentrated Matrix Exponential) distribution.
 */
public class CMEDistributionTest {

    @BeforeAll
    public static void setUp() {
        Maths.setRandomNumbersMatlab(true);
        RandomManager.setMasterSeed(23000);
    }

    @Test
    public void testIsAnME() {
        CME cme = new CME(1.0, 7);
        assertTrue(cme instanceof ME);
        assertEquals("ME", cme.getName());
        assertEquals(7, cme.getNumberOfPhases());
        assertEquals(7, cme.getOrder());
    }

    @Test
    public void testRepresentationIsValid() {
        for (int order : new int[]{3, 7, 21, 101}) {
            Matrix[] rep = CME.representation(order);
            assertTrue(checkMERepresentation(rep[0], rep[1], 1e-14),
                    "CME of order " + order + " must be a valid ME representation");
            double sum = 0.0;
            for (int i = 0; i < rep[0].getNumElements(); i++) {
                sum += rep[0].get(0, i);
            }
            assertEquals(1.0, sum, 1e-10, "alpha must sum to one at order " + order);
        }
    }

    @Test
    public void testMeanAndSCV() {
        for (int order : new int[]{3, 7, 21, 101}) {
            CME cme = new CME(2.5, order);
            assertEquals(2.5, cme.getMean(), 1e-9);
            // Relative 1e-6: the SCV goes through two inversions of a (2n+1)-square
            // matrix, and at order 101 Matrix.inv loses about eight digits on this
            // spectrum. Native Python, which inverts with LAPACK, reproduces the
            // tabulated cv2 to 1e-13.
            assertEquals(CME.getMinSCV(order), cme.getSCV(), 1e-6 * CME.getMinSCV(order));
        }
    }

    @Test
    public void testSCVBelowErlangBound() {
        // The point of a CME: it beats the phase-type bound 1/order at equal order.
        for (int order : new int[]{7, 21, 101}) {
            assertTrue(CME.getMinSCV(order) < 1.0 / order,
                    "CME of order " + order + " must be more concentrated than Erlang-" + order);
        }
    }

    @Test
    public void testDensityMatchesTabulatedForm() {
        // f(x) = mu1*exp(-mu1*x)*(c + sum_k [a_k cos(k w mu1 x) + b_k sin(k w mu1 x)])
        jline.lib.lti.iltcme.CmeEntry entry = CME.tableEntry(21);
        CME cme = new CME(1.0, 21);
        for (double x : new double[]{0.3, 0.8, 1.0, 1.4}) {
            double closed = 0.0;
            for (int k = 1; k <= entry.n; k++) {
                closed += entry.a.get(k - 1) * Math.cos(k * entry.omega * entry.mu1 * x)
                        + entry.b.get(k - 1) * Math.sin(k * entry.omega * entry.mu1 * x);
            }
            closed = entry.mu1 * Math.exp(-entry.mu1 * x) * (entry.c + closed);
            double pdf = Map_pdf.map_pdf(cme.getProcess(), new double[]{x})[0];
            assertEquals(closed, pdf, 1e-8 * Math.max(1.0, Math.abs(closed)));
        }
    }

    @Test
    public void testDensityIsNonNegative() {
        CME cme = new CME(1.0, 21);
        double[] tset = new double[401];
        for (int i = 0; i < tset.length; i++) {
            tset[i] = i * 0.01;
        }
        double[] pdf = Map_pdf.map_pdf(cme.getProcess(), tset);
        for (int i = 0; i < tset.length; i++) {
            assertTrue(pdf[i] > -1e-10, "CME density must be nonnegative at x = " + tset[i]);
        }
    }

    @Test
    public void testScalingWithMean() {
        // The representation scales as A/mean, so the CDF is a pure time rescaling.
        CME unit = new CME(1.0, 11);
        CME scaled = new CME(3.0, 11);
        for (double x : new double[]{0.2, 0.7, 1.3, 2.0}) {
            assertEquals(unit.evalCDF(x), scaled.evalCDF(3.0 * x), 1e-9 * Math.max(1.0, unit.evalCDF(x)));
        }
        assertEquals(unit.getSCV(), scaled.getSCV(), 1e-12);
    }

    @Test
    public void testFitMeanAndSCV() {
        CME cme = CME.fitMeanAndSCV(4.0, 1e-3);
        assertEquals(4.0, cme.getMean(), 1e-9);
        assertTrue(cme.getSCV() <= 1e-3);
        // Minimality: no lower tabulated order reaches the target.
        for (int order : CME.getSupportedOrders()) {
            if (order < cme.getOrder()) {
                assertTrue(CME.getMinSCV(order) > 1e-3);
            }
        }
    }

    @Test
    public void testSampling() {
        CME cme = new CME(1.0, 11);
        double[] samples = cme.sample(20000, new Random(42));
        double mean = 0.0;
        for (double s : samples) {
            assertTrue(s >= 0.0, "CME samples must be nonnegative");
            mean += s;
        }
        mean /= samples.length;
        double var = 0.0;
        for (double s : samples) {
            var += (s - mean) * (s - mean);
        }
        var /= samples.length;
        assertEquals(1.0, mean, 0.02);
        assertEquals(cme.getSCV(), var / (mean * mean), 1e-3);
    }

    @Test
    public void testInvalidArguments() {
        assertThrows(IllegalArgumentException.class, () -> new CME(1.0, 4));
        assertThrows(IllegalArgumentException.class, () -> new CME(1.0, 1));
        assertThrows(IllegalArgumentException.class, () -> new CME(0.0, 7));
        assertThrows(IllegalArgumentException.class, () -> new CME(-1.0, 7));
        // 2003 phases is beyond the last tabulated entry (n = 1000).
        assertThrows(IllegalArgumentException.class, () -> new CME(1.0, 2003));
        assertThrows(IllegalArgumentException.class, () -> CME.fitMeanAndSCV(1.0, 1e-12));
    }

    @Test
    public void testMG1Queue() {
        // M/CME/1 at rho = 0.5 with a near-deterministic service: the
        // Pollaczek-Khinchine mean queue length is rho^2*(1+scv)/(2*(1-rho)) plus rho.
        // Mirrors the same check in the MATLAB and Python test suites.
        Network model = new Network("M/CME/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        CME service = new CME(1.0, 11);
        source.setArrival(oclass, new Exp(0.5));
        queue.setService(oclass, service);
        model.link(model.serialRouting(source, queue, sink));

        // sn.isph marks the CME station as non-Markovian: its (D0,D1) pair has
        // negative off-diagonal entries, so mu/phi/pie carry no probabilistic
        // reading there, while the exponential source stays Markovian.
        NetworkStruct sn = model.getStruct();
        assertFalse(sn.isph.get(queue).get(oclass));
        assertTrue(sn.isph.get(source).get(oclass));

        NetworkAvgTable avgTable = new SolverMAM(model).getAvgTable();
        double qlen = avgTable.getQLen().get(1);

        double rho = 0.5;
        double pk = rho + rho * rho * (1.0 + service.getSCV()) / (2.0 * (1.0 - rho));
        assertEquals(pk, qlen, 1e-3 * pk);
    }

    private static Network mme1(double rho, int order) {
        Network model = new Network("M/CME/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(rho));
        queue.setService(oclass, new CME(1.0, order));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    @Test
    public void testCtmcMG1IsExact() {
        // The ME embeds in the generator exactly as a phase-type does, keeping the
        // negative off-diagonal entries of A. The stationary vector is then a signed
        // measure, but every aggregate over a phase block is exact, so the mean queue
        // length must reproduce Pollaczek-Khinchine to solver precision.
        int[] orders = {3, 7, 11};
        double[] rhos = {0.3, 0.5};
        for (int order : orders) {
            for (double rho : rhos) {
                SolverOptions options = SolverCTMC.defaultOptions();
                options.cutoff = new Matrix(1, 1, 1);
                options.cutoff.set(0, 0, 30);
                options.verbose = jline.VerboseLevel.SILENT;
                NetworkAvgTable avg = new SolverCTMC(mme1(rho, order), options).getAvgTable();
                double scv = new CME(1.0, order).getSCV();
                double pk = rho + rho * rho * (1.0 + scv) / (2.0 * (1.0 - rho));
                assertEquals(pk, avg.getQLen().get(1), 1e-8 * pk);
                assertEquals(rho, avg.getUtil().get(1), 1e-8);
                assertEquals(rho, avg.getTput().get(1), 1e-8);
            }
        }
    }

    @Test
    public void testCtmcRefusesPerStateProbabilities() {
        // Per-state probabilities and uniformization-based transients do not exist for a
        // signed stationary vector, so they are refused rather than returned.
        SolverOptions options = SolverCTMC.defaultOptions();
        options.cutoff = new Matrix(1, 1, 1);
        options.cutoff.set(0, 0, 5);
        options.verbose = jline.VerboseLevel.SILENT;
        SolverCTMC solver = new SolverCTMC(mme1(0.5, 3), options);
        solver.getAvgTable();
        assertThrows(RuntimeException.class, solver::getProbSysAggr);
        assertThrows(RuntimeException.class, solver::getTranProbSysAggr);
    }

    @Test
    public void testFitMeanAndSCVIsExact() {
        // The CME-plus-exponential convolution matches both moments exactly over the whole
        // range (sY/(1+sY), 1), which is where an Erlang needs ceil(1/scv) phases and a
        // two-phase Coxian cannot go at all.
        double[] scvs = {0.9, 0.5, 0.2, 0.05, 0.01, 1e-3, 1e-4};
        for (double scv : scvs) {
            ME fitted = MEFit.fitMeanAndSCV(2.0, scv);
            assertEquals(2.0, fitted.getMean(), 1e-9);
            assertEquals(scv, fitted.getSCV(), 1e-6 * scv);
            if (scv <= 0.05) {
                // O(1/n^2) instead of the Erlang O(1/n).
                assertTrue(fitted.getNumberOfPhases() < Math.ceil(1.0 / scv));
            }
        }
    }

    @Test
    public void testFitMeanAndSCVRespectsPhaseBudget() {
        // Under a budget the fit returns the closest achievable SCV from below rather than
        // silently truncating an Erlang: 20 phases reach 5.7e-3, Erlang-20 stops at 0.05.
        ME fitted = MEFit.fitMeanAndSCV(1.0, 1e-4, 20);
        assertTrue(fitted.getNumberOfPhases() <= 20);
        assertEquals(1.0, fitted.getMean(), 1e-9);
        assertTrue(fitted.getSCV() < 1.0 / 20);
    }

    @Test
    public void testFitMeanAndSCVRejectsOutOfRange() {
        assertThrows(IllegalArgumentException.class, () -> MEFit.fitMeanAndSCV(1.0, 0.0));
        assertThrows(IllegalArgumentException.class, () -> MEFit.fitMeanAndSCV(1.0, 1.0));
        assertThrows(IllegalArgumentException.class, () -> MEFit.fitMeanAndSCV(1.0, 1.5));
        assertThrows(IllegalArgumentException.class, () -> MEFit.fitMeanAndSCV(0.0, 0.5));
    }

    /**
     * An ME with a negative entry in alpha whose density touches zero in the interior
     * (f(1.2) = 0 against a peak of 2.79), so it admits NO phase-type representation of
     * any order. Spectrum {-1, -2 +- 2i}, mean 0.5329, SCV 2.6843. Same instance as
     * line-test.git/test/testsAPI/test_mam_me_warning.m and the native Python test suite.
     */
    private static ME nonPhaseTypeME() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.61058991931158258);
        alpha.set(0, 1, -0.15547146730086722);
        alpha.set(0, 2, 0.54488154798928464);
        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0);
        A.set(1, 1, -2.0);
        A.set(1, 2, 2.0);
        A.set(2, 1, -2.0);
        A.set(2, 2, -2.0);
        return new ME(alpha, A);
    }

    @Test
    public void testNonPhaseTypeMEHasAnInteriorDensityZero() {
        // The interior zero is what rules out a phase-type representation of ANY order: a
        // PH density is strictly positive throughout the interior of its support.
        // Round-off can make it slightly negative, hence the tolerance. The zero sits at
        // x = 1.2 exactly, so the grid has to land on it.
        ME dist = nonPhaseTypeME();
        double[] fine = new double[2001];
        for (int i = 0; i < fine.length; i++) {
            fine[i] = 1.15 + i * (0.10 / (fine.length - 1));
        }
        double[] pdf = Map_pdf.map_pdf(dist.getProcess(), fine);
        double min = Double.POSITIVE_INFINITY;
        for (double v : pdf) {
            min = Math.min(min, v);
        }
        assertTrue(Math.abs(min) < 1e-10, "the density must reach zero at x = 1.2, got " + min);
        double[] bulk = Map_pdf.map_pdf(dist.getProcess(), new double[]{0.3, 0.6, 2.0, 3.0});
        double max = 0.0;
        for (double v : bulk) {
            max = Math.max(max, v);
        }
        assertTrue(max > 0.1, "the density must be positive away from the zero");
        assertTrue(dist.getAlpha().get(0, 1) < 0, "alpha must have a negative entry");
    }

    @Test
    public void testCtmcSolvesANonPhaseTypeME() {
        // The CTMC branch must not depend on the ME happening to be a phase-type in
        // disguise. Pollaczek-Khinchine still applies, since the service law is a genuine
        // distribution and the arrivals are Poisson.
        double[] rhos = {0.3, 0.5};
        int[] cutoffs = {30, 60};
        for (int i = 0; i < rhos.length; i++) {
            double rho = rhos[i];
            ME service = nonPhaseTypeME();
            double mean = service.getMean();
            double scv = service.getSCV();

            Network model = new Network("M/ME/1 non-PH");
            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");
            OpenClass oclass = new OpenClass(model, "Class1");
            source.setArrival(oclass, new Exp(rho / mean));
            queue.setService(oclass, nonPhaseTypeME());
            model.link(model.serialRouting(source, queue, sink));

            assertFalse(model.getStruct().isph.get(queue).get(oclass));

            SolverOptions options = SolverCTMC.defaultOptions();
            options.cutoff = new Matrix(1, 1, 1);
            options.cutoff.set(0, 0, cutoffs[i]);
            options.verbose = jline.VerboseLevel.SILENT;
            NetworkAvgTable avg = new SolverCTMC(model, options).getAvgTable();
            double pk = rho + rho * rho * (1.0 + scv) / (2.0 * (1.0 - rho));
            assertEquals(pk, avg.getQLen().get(1), 1e-7 * pk);
            assertEquals(rho, avg.getUtil().get(1), 1e-6);
        }
    }

    @Test
    public void testSupportedOrders() {
        int[] orders = CME.getSupportedOrders();
        assertTrue(orders.length > 100);
        assertEquals(3, orders[0]);
        for (int i = 1; i < orders.length; i++) {
            assertTrue(orders[i] > orders[i - 1]);
            assertEquals(1, orders[i] % 2);
        }
    }
}
