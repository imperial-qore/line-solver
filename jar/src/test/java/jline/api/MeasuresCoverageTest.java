package jline.api;

import jline.api.measures.Ms_anderson_darlingKt;
import jline.api.measures.Ms_cramer_von_misesKt;
import jline.api.measures.Ms_kolmogorov_smirnovKt;
import jline.api.measures.Ms_kuiperKt;
import jline.api.measures.Ms_wassersteinKt;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Coverage tests for statistical measures with 0% coverage.
 */
public class MeasuresCoverageTest {

    private static final double TOL = 1e-6;

    // ========== Anderson-Darling ==========

    @Test
    void testAndersonDarling() {
        // Create two identical samples - should have distance 0
        Matrix sample1 = new Matrix(5, 1);
        sample1.set(0, 0, 1.0);
        sample1.set(1, 0, 2.0);
        sample1.set(2, 0, 3.0);
        sample1.set(3, 0, 4.0);
        sample1.set(4, 0, 5.0);

        Matrix sample2 = new Matrix(5, 1);
        sample2.set(0, 0, 1.0);
        sample2.set(1, 0, 2.0);
        sample2.set(2, 0, 3.0);
        sample2.set(3, 0, 4.0);
        sample2.set(4, 0, 5.0);

        double ad = Ms_anderson_darlingKt.ms_anderson_darling(sample1, sample2);
        assertEquals(0.0, ad, TOL);

        // Create two different samples
        Matrix sample3 = new Matrix(5, 1);
        sample3.set(0, 0, 10.0);
        sample3.set(1, 0, 20.0);
        sample3.set(2, 0, 30.0);
        sample3.set(3, 0, 40.0);
        sample3.set(4, 0, 50.0);

        double adDiff = Ms_anderson_darlingKt.ms_anderson_darling(sample1, sample3);
        assertTrue(adDiff > 0.0);
    }

    // ========== Kolmogorov-Smirnov ==========

    @Test
    void testKolmogorovSmirnov() {
        // Create two identical samples - KS distance should be 0
        Matrix sample1 = new Matrix(5, 1);
        sample1.set(0, 0, 1.0);
        sample1.set(1, 0, 2.0);
        sample1.set(2, 0, 3.0);
        sample1.set(3, 0, 4.0);
        sample1.set(4, 0, 5.0);

        Matrix sample2 = new Matrix(5, 1);
        sample2.set(0, 0, 1.0);
        sample2.set(1, 0, 2.0);
        sample2.set(2, 0, 3.0);
        sample2.set(3, 0, 4.0);
        sample2.set(4, 0, 5.0);

        double ks = Ms_kolmogorov_smirnovKt.ms_kolmogorov_smirnov(sample1, sample2);
        assertEquals(0.0, ks, TOL);

        // Create two completely different samples
        Matrix sample3 = new Matrix(5, 1);
        sample3.set(0, 0, 100.0);
        sample3.set(1, 0, 200.0);
        sample3.set(2, 0, 300.0);
        sample3.set(3, 0, 400.0);
        sample3.set(4, 0, 500.0);

        double ksDiff = Ms_kolmogorov_smirnovKt.ms_kolmogorov_smirnov(sample1, sample3);
        // KS distance should be 1.0 for non-overlapping distributions
        assertEquals(1.0, ksDiff, TOL);
    }

    // ========== Cramer-von Mises ==========

    @Test
    void testCramerVonMises() {
        // Create two identical samples
        Matrix sample1 = new Matrix(5, 1);
        sample1.set(0, 0, 1.0);
        sample1.set(1, 0, 2.0);
        sample1.set(2, 0, 3.0);
        sample1.set(3, 0, 4.0);
        sample1.set(4, 0, 5.0);

        Matrix sample2 = new Matrix(5, 1);
        sample2.set(0, 0, 1.0);
        sample2.set(1, 0, 2.0);
        sample2.set(2, 0, 3.0);
        sample2.set(3, 0, 4.0);
        sample2.set(4, 0, 5.0);

        double cvm = Ms_cramer_von_misesKt.ms_cramer_von_mises(sample1, sample2);
        assertEquals(0.0, cvm, TOL);

        // Create two different samples
        Matrix sample3 = new Matrix(5, 1);
        sample3.set(0, 0, 10.0);
        sample3.set(1, 0, 20.0);
        sample3.set(2, 0, 30.0);
        sample3.set(3, 0, 40.0);
        sample3.set(4, 0, 50.0);

        double cvmDiff = Ms_cramer_von_misesKt.ms_cramer_von_mises(sample1, sample3);
        assertTrue(cvmDiff > 0.0);
    }

    // ========== Wasserstein ==========

    @Test
    void testWasserstein() {
        // Create two identical samples
        Matrix sample1 = new Matrix(5, 1);
        sample1.set(0, 0, 1.0);
        sample1.set(1, 0, 2.0);
        sample1.set(2, 0, 3.0);
        sample1.set(3, 0, 4.0);
        sample1.set(4, 0, 5.0);

        Matrix sample2 = new Matrix(5, 1);
        sample2.set(0, 0, 1.0);
        sample2.set(1, 0, 2.0);
        sample2.set(2, 0, 3.0);
        sample2.set(3, 0, 4.0);
        sample2.set(4, 0, 5.0);

        double ws = Ms_wassersteinKt.ms_wasserstein(sample1, sample2);
        assertEquals(0.0, ws, TOL);

        // Create a shifted sample
        Matrix sample3 = new Matrix(5, 1);
        sample3.set(0, 0, 2.0);  // shifted by 1
        sample3.set(1, 0, 3.0);
        sample3.set(2, 0, 4.0);
        sample3.set(3, 0, 5.0);
        sample3.set(4, 0, 6.0);

        double wsDiff = Ms_wassersteinKt.ms_wasserstein(sample1, sample3);
        // Wasserstein distance for a shift of 1 should be close to 1
        assertTrue(wsDiff > 0.0);
    }

    // ========== Kuiper ==========

    @Test
    void testKuiper() {
        // Create two identical samples
        Matrix sample1 = new Matrix(5, 1);
        sample1.set(0, 0, 1.0);
        sample1.set(1, 0, 2.0);
        sample1.set(2, 0, 3.0);
        sample1.set(3, 0, 4.0);
        sample1.set(4, 0, 5.0);

        Matrix sample2 = new Matrix(5, 1);
        sample2.set(0, 0, 1.0);
        sample2.set(1, 0, 2.0);
        sample2.set(2, 0, 3.0);
        sample2.set(3, 0, 4.0);
        sample2.set(4, 0, 5.0);

        double kuiper = Ms_kuiperKt.ms_kuiper(sample1, sample2);
        assertEquals(0.0, kuiper, TOL);

        // Create two different samples
        Matrix sample3 = new Matrix(5, 1);
        sample3.set(0, 0, 10.0);
        sample3.set(1, 0, 20.0);
        sample3.set(2, 0, 30.0);
        sample3.set(3, 0, 40.0);
        sample3.set(4, 0, 50.0);

        double kuiperDiff = Ms_kuiperKt.ms_kuiper(sample1, sample3);
        assertTrue(kuiperDiff > 0.0);
    }

    // ========== Multi-column (2D) tests ==========

    @Test
    void testMultiColumnSamples() {
        // Create 2D matrices (5 rows x 2 columns)
        Matrix sample1 = new Matrix(5, 2);
        sample1.set(0, 0, 1.0); sample1.set(0, 1, 10.0);
        sample1.set(1, 0, 2.0); sample1.set(1, 1, 20.0);
        sample1.set(2, 0, 3.0); sample1.set(2, 1, 30.0);
        sample1.set(3, 0, 4.0); sample1.set(3, 1, 40.0);
        sample1.set(4, 0, 5.0); sample1.set(4, 1, 50.0);

        Matrix sample2 = new Matrix(5, 2);
        sample2.set(0, 0, 1.0); sample2.set(0, 1, 10.0);
        sample2.set(1, 0, 2.0); sample2.set(1, 1, 20.0);
        sample2.set(2, 0, 3.0); sample2.set(2, 1, 30.0);
        sample2.set(3, 0, 4.0); sample2.set(3, 1, 40.0);
        sample2.set(4, 0, 5.0); sample2.set(4, 1, 50.0);

        // All measures should return 0 for identical samples
        assertEquals(0.0, Ms_kolmogorov_smirnovKt.ms_kolmogorov_smirnov(sample1, sample2), TOL);
        assertEquals(0.0, Ms_wassersteinKt.ms_wasserstein(sample1, sample2), TOL);
        assertEquals(0.0, Ms_kuiperKt.ms_kuiper(sample1, sample2), TOL);
        assertEquals(0.0, Ms_cramer_von_misesKt.ms_cramer_von_mises(sample1, sample2), TOL);
        assertEquals(0.0, Ms_anderson_darlingKt.ms_anderson_darling(sample1, sample2), TOL);
    }
}
