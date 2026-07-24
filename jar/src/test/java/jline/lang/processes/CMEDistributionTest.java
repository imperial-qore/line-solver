/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.api.mam.Map_pdf;
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
        assertEquals(cme.getSCV(), var / (mean * mean), 0.01);
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
