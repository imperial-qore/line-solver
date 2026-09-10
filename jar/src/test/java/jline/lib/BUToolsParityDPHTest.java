package jline.lib;

import jline.lib.butools.dph.CdfFromDPH;
import jline.lib.butools.dph.CdfFromMG;
import jline.lib.butools.dph.PmfFromDPH;
import jline.lib.butools.dph.PmfFromMG;
import jline.lib.butools.dph.MomentsFromDPH;
import jline.lib.butools.dph.CheckDPHRepresentation;
import jline.lib.butools.dph.DPH2From3Moments;
import jline.lib.butools.dph.CanonicalFromDPH2.DPH2Representation;
import jline.lib.butools.dph.DPH3From5Moments;
import jline.lib.butools.dph.CanonicalFromDPH3.DPH3Representation;
import jline.lib.butools.dph.CanonicalFromDPH2;
import jline.lib.butools.dph.CanonicalFromDPH3;
import jline.lib.butools.dph.DPHFromMG;
import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.lib.butools.dph.AcyclicDPHFromMG;
import jline.lib.butools.dph.MomentsFromMG;
import jline.lib.butools.dph.CheckMGRepresentation;
import jline.lib.butools.dph.MGFromMoments;
import jline.lib.butools.dph.RandomDPH;
import jline.lib.butools.dph.SamplesFromDPH;
import jline.lib.butools.trace.MarginalMomentsFromTrace;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.*;

import java.util.Random;

/**
 * BUTools Parity Tests - DPH/MG functions
 *
 * Each test uses EXACT data from butools.tmp example files to ensure the JAR implementation
 * produces identical results to the original BUTools library.
 */
public class BUToolsParityDPHTest {


    // ============ MomentsFromDPH - From butools.tmp/MomentsFromDPH.txt ============

    /**
     * MomentsFromDPH.txt example:
     * Input: a=[0.76,0,0.24], A=[0.34,0.66,0; 0.79,0.05,0.07; 0.26,0.73,0.01], K=5
     * Expected: [26.995, 1398, 1.0853e+05, 1.1233e+07, 1.4533e+09]
     */
    @Test
    public void testMomentsFromDPH() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.76);
        alpha.set(0, 1, 0.0);
        alpha.set(0, 2, 0.24);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.34); A.set(0, 1, 0.66); A.set(0, 2, 0.0);
        A.set(1, 0, 0.79); A.set(1, 1, 0.05); A.set(1, 2, 0.07);
        A.set(2, 0, 0.26); A.set(2, 1, 0.73); A.set(2, 2, 0.01);

        double[] moms = MomentsFromDPH.momentsFromDPH(alpha, A, 5);

        assertNotNull(moms);
        assertEquals(5, moms.length);
        assertEquals(26.995, moms[0], LOOSE_MID_TOL);
        assertEquals(1398.0, moms[1], 1.0);
        assertEquals(1.0853e+05, moms[2], 10.0);
        assertEquals(1.1233e+07, moms[3], 1000.0);
        assertEquals(1.4533e+09, moms[4], 100000.0);
    }

    // ============ CheckDPHRepresentation - From butools.tmp/CheckDPHRepresentation.txt ============

    /**
     * CheckDPHRepresentation.txt Test 1: Valid DPH representation
     * Input: a=[0.48,0.08,0.26,0.18], A=[0,0.08,0.08,0.8; 0.55,0,0.24,0.19; 0.06,0.03,0,0.001; 0.23,0.005,0.2,0.53]
     * Expected: flag=true
     */
    @Test
    public void testCheckDPHRepresentation_Valid() {
        Matrix alpha = new Matrix(1, 4);
        alpha.set(0, 0, 0.48);
        alpha.set(0, 1, 0.08);
        alpha.set(0, 2, 0.26);
        alpha.set(0, 3, 0.18);

        Matrix A = new Matrix(4, 4);
        A.set(0, 0, 0.0);   A.set(0, 1, 0.08);  A.set(0, 2, 0.08);  A.set(0, 3, 0.8);
        A.set(1, 0, 0.55);  A.set(1, 1, 0.0);   A.set(1, 2, 0.24);  A.set(1, 3, 0.19);
        A.set(2, 0, 0.06);  A.set(2, 1, 0.03);  A.set(2, 2, 0.0);   A.set(2, 3, 0.001);
        A.set(3, 0, 0.23);  A.set(3, 1, 0.005); A.set(3, 2, 0.2);   A.set(3, 3, 0.53);

        boolean flag = CheckDPHRepresentation.checkDPHRepresentation(alpha, A, 1e-12);
        assertTrue(flag, "Valid DPH representation should pass check");
    }

    /**
     * CheckDPHRepresentation.txt Test 2: Invalid DPH - row sum exceeds 1
     * Input: a=[0.48,0.08], A=[0,0.08; 0.55,0.5]
     * Expected: flag=false
     */
    @Test
    public void testCheckDPHRepresentation_Invalid() {
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, 0.48);
        alpha.set(0, 1, 0.08);

        Matrix A = new Matrix(2, 2);
        A.set(0, 0, 0.0);  A.set(0, 1, 0.08);
        A.set(1, 0, 0.55); A.set(1, 1, 0.5);

        boolean flag = CheckDPHRepresentation.checkDPHRepresentation(alpha, A, 1e-12);
        assertFalse(flag, "Invalid DPH (row sum > 1) should fail check");
    }

    // ============ DPH2From3Moments - From butools.tmp/DPH2From3Moments.txt ============

    /**
     * DPH2From3Moments.txt example:
     * First compute moments: a=[0.9,0.1], A=[0.2,0.61; 0.58,0.41] -> moms=[10.305, 215.13, 6764.2]
     * Then: DPH2From3Moments(moms) -> b=[0.43249, 0.56751], B=[0.61,0.39; 0.69692,0]
     * Verify moments roundtrip.
     */
    @Test
    public void testDPH2From3Moments() {
        // First compute moments from original DPH
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, 0.9);
        alpha.set(0, 1, 0.1);

        Matrix A = new Matrix(2, 2);
        A.set(0, 0, 0.2);  A.set(0, 1, 0.61);
        A.set(1, 0, 0.58); A.set(1, 1, 0.41);

        double[] moms = MomentsFromDPH.momentsFromDPH(alpha, A, 3);
        assertEquals(10.305, moms[0], LOOSE_MID_TOL);
        assertEquals(215.13, moms[1], VERY_COARSE_TOL);
        assertEquals(6764.2, moms[2], 1.0);

        // Construct DPH2 from moments
        DPH2Representation result = DPH2From3Moments.dph2From3Moments(moms);
        assertNotNull(result);

        Matrix beta = result.beta;
        Matrix B = result.B;

        // Verify beta = [0.43249, 0.56751]
        assertEquals(0.43249, beta.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.56751, beta.get(0, 1), LOOSE_MID_TOL);

        // Verify B = [0.61, 0.39; 0.69692, 0]
        assertEquals(0.61, B.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.39, B.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.69692, B.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(1, 1), LOOSE_MID_TOL);

        // Verify moments roundtrip
        double[] phmoms = MomentsFromDPH.momentsFromDPH(beta, B, 3);
        assertEquals(moms[0], phmoms[0], LOOSE_MID_TOL);
        assertEquals(moms[1], phmoms[1], VERY_COARSE_TOL);
        assertEquals(moms[2], phmoms[2], 1.0);
    }

    // ============ DPH3From5Moments - From butools.tmp/DPH3From5Moments.txt ============

    /**
     * DPH3From5Moments.txt example:
     * First compute moments: a=[0.7,0.1,0.2], A=[0.2,0.51,0.1; 0.58,0.41,0; 0.1,0.4,0.3]
     *   -> moms=[9.3096, 175.1, 4968.7, 1.8805e+05, 8.8966e+06]
     * Then: DPH3From5Moments(moms)
     *   -> b=[0.73989, 0.076837, 0.18327], B=[0.89971,0.10029,0; 0,0.010293,0.98971; 0,0.050581,0]
     * Verify moments roundtrip using MomentsFromMG.
     */
    @Test
    public void testDPH3From5Moments() {
        // First compute moments from original DPH
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.7);
        alpha.set(0, 1, 0.1);
        alpha.set(0, 2, 0.2);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.2);  A.set(0, 1, 0.51); A.set(0, 2, 0.1);
        A.set(1, 0, 0.58); A.set(1, 1, 0.41); A.set(1, 2, 0.0);
        A.set(2, 0, 0.1);  A.set(2, 1, 0.4);  A.set(2, 2, 0.3);

        double[] moms = MomentsFromDPH.momentsFromDPH(alpha, A, 5);
        assertEquals(9.3096, moms[0], LOOSE_MID_TOL);
        assertEquals(175.1, moms[1], VERY_COARSE_TOL);
        assertEquals(4968.7, moms[2], 1.0);
        assertEquals(1.8805e+05, moms[3], 100.0);
        assertEquals(8.8966e+06, moms[4], 10000.0);

        // Construct DPH3 from moments
        DPH3Representation result = DPH3From5Moments.dph3From5Moments(moms, 1e-14);
        assertNotNull(result);

        Matrix beta = result.beta;
        Matrix B = result.B;

        // Verify beta = [0.73989, 0.076837, 0.18327]
        assertEquals(0.73989, beta.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.076837, beta.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.18327, beta.get(0, 2), LOOSE_MID_TOL);

        // Verify B = [0.89971,0.10029,0; 0,0.010293,0.98971; 0,0.050581,0]
        assertEquals(0.89971, B.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.10029, B.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(0, 2), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.010293, B.get(1, 1), LOOSE_MID_TOL);
        assertEquals(0.98971, B.get(1, 2), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(2, 0), LOOSE_MID_TOL);
        assertEquals(0.050581, B.get(2, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(2, 2), LOOSE_MID_TOL);

        // Verify moments roundtrip using MomentsFromMG
        double[] phmoms = MomentsFromMG.momentsFromMG(beta, B, 5);
        assertEquals(moms[0], phmoms[0], LOOSE_MID_TOL);
        assertEquals(moms[1], phmoms[1], VERY_COARSE_TOL);
        assertEquals(moms[2], phmoms[2], 1.0);
        assertEquals(moms[3], phmoms[3], 100.0);
        assertEquals(moms[4], phmoms[4], 10000.0);
    }

    // ============ CanonicalFromDPH2 - From butools.tmp/CanonicalFromDPH2.txt ============

    /**
     * CanonicalFromDPH2.txt Test 1:
     * Input: a=[0,1.0], A=[0.23,0.22; 0.41,0.48]
     * Expected: b=[0.88663, 0.11337], B=[0.68031,0.31969; 0,0.029692]
     */
    @Test
    public void testCanonicalFromDPH2_Test1() {
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, 0.0);
        alpha.set(0, 1, 1.0);

        Matrix A = new Matrix(2, 2);
        A.set(0, 0, 0.23); A.set(0, 1, 0.22);
        A.set(1, 0, 0.41); A.set(1, 1, 0.48);

        DPH2Representation result = CanonicalFromDPH2.canonicalFromDPH2(alpha, A, 1e-14);
        assertNotNull(result);

        Matrix beta = result.beta;
        Matrix B = result.B;

        // Verify beta = [0.88663, 0.11337]
        assertEquals(0.88663, beta.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.11337, beta.get(0, 1), LOOSE_MID_TOL);

        // Verify B = [0.68031, 0.31969; 0, 0.029692]
        assertEquals(0.68031, B.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.31969, B.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.029692, B.get(1, 1), LOOSE_MID_TOL);

        // Verify result is valid DPH
        boolean flag = CheckDPHRepresentation.checkDPHRepresentation(beta, B, 1e-12);
        assertTrue(flag, "Canonical form should be valid DPH");
    }

    /**
     * CanonicalFromDPH2.txt Test 2:
     * Input: a=[1.0,0], A=[0,0.61; 0.56,0.44]
     * Expected: b~=[0, 1], B=[0.44,0.56; 0.61,0]
     */
    @Test
    public void testCanonicalFromDPH2_Test2() {
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, 1.0);
        alpha.set(0, 1, 0.0);

        Matrix A = new Matrix(2, 2);
        A.set(0, 0, 0.0);  A.set(0, 1, 0.61);
        A.set(1, 0, 0.56); A.set(1, 1, 0.44);

        DPH2Representation result = CanonicalFromDPH2.canonicalFromDPH2(alpha, A, 1e-14);
        assertNotNull(result);

        Matrix beta = result.beta;
        Matrix B = result.B;

        // Verify beta ~= [0, 1] (first element near zero)
        assertEquals(0.0, beta.get(0, 0), LOOSE_MID_TOL);
        assertEquals(1.0, beta.get(0, 1), LOOSE_MID_TOL);

        // Verify B = [0.44, 0.56; 0.61, 0]
        assertEquals(0.44, B.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.56, B.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.61, B.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(1, 1), LOOSE_MID_TOL);

        // Verify result is valid DPH
        boolean flag = CheckDPHRepresentation.checkDPHRepresentation(beta, B, 1e-12);
        assertTrue(flag, "Canonical form should be valid DPH");
    }

    // ============ CanonicalFromDPH3 - From butools.tmp/CanonicalFromDPH3.txt ============

    /**
     * CanonicalFromDPH3.txt Test 1:
     * Input: a=[0.46,0.22,0.32], A=[0.67,0.01,0.12; 0.06,0.45,0.15; 0.18,0.43,0.32]
     * Expected: b=[0.21239, 0.37004, 0.41757], B=[0.10918,0,0; 0.45654,0.54346,0; 0,0.21265,0.78735]
     */
    @Test
    public void testCanonicalFromDPH3() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.46);
        alpha.set(0, 1, 0.22);
        alpha.set(0, 2, 0.32);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.67); A.set(0, 1, 0.01); A.set(0, 2, 0.12);
        A.set(1, 0, 0.06); A.set(1, 1, 0.45); A.set(1, 2, 0.15);
        A.set(2, 0, 0.18); A.set(2, 1, 0.43); A.set(2, 2, 0.32);

        DPH3Representation result = CanonicalFromDPH3.canonicalFromDPH3(alpha, A, 1e-14);
        assertNotNull(result);

        Matrix beta = result.beta;
        Matrix B = result.B;

        // Verify beta = [0.21239, 0.37004, 0.41757]
        assertEquals(0.21239, beta.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.37004, beta.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.41757, beta.get(0, 2), LOOSE_MID_TOL);

        // Verify B = [0.10918,0,0; 0.45654,0.54346,0; 0,0.21265,0.78735]
        assertEquals(0.10918, B.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(0, 2), LOOSE_MID_TOL);
        assertEquals(0.45654, B.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.54346, B.get(1, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(1, 2), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(2, 0), LOOSE_MID_TOL);
        assertEquals(0.21265, B.get(2, 1), LOOSE_MID_TOL);
        assertEquals(0.78735, B.get(2, 2), LOOSE_MID_TOL);

        // Verify result is valid DPH
        boolean flag = CheckDPHRepresentation.checkDPHRepresentation(beta, B, 1e-12);
        assertTrue(flag, "Canonical form should be valid DPH");
    }

    // ============ DPHFromMG - From butools.tmp/DPHFromMG.txt ============

    /**
     * DPHFromMG.txt example:
     * Input: a=[-0.6,0.3,1.3], A=[0.1,0.2,0; 0.3,0.1,0.25; -0.3,0.2,0.77]
     * Expected: b=[0.05, 0.1375, 0.8125], B=[0.1,0.2,0; 0.425,0.06875,0.15625; 0.141,0.01975,0.80125]
     * Verify result is valid DPH.
     */
    @Test
    public void testDPHFromMG() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, -0.6);
        alpha.set(0, 1, 0.3);
        alpha.set(0, 2, 1.3);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.1);  A.set(0, 1, 0.2);  A.set(0, 2, 0.0);
        A.set(1, 0, 0.3);  A.set(1, 1, 0.1);  A.set(1, 2, 0.25);
        A.set(2, 0, -0.3); A.set(2, 1, 0.2);  A.set(2, 2, 0.77);

        MGRepresentation result = DPHFromMG.dphFromMG(alpha, A, 1e-14);
        assertNotNull(result);

        Matrix beta = result.getAlpha();
        Matrix B = result.getA();

        // Verify beta = [0.05, 0.1375, 0.8125]
        assertEquals(0.05, beta.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.1375, beta.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.8125, beta.get(0, 2), LOOSE_MID_TOL);

        // Verify B = [0.1,0.2,0; 0.425,0.06875,0.15625; 0.141,0.01975,0.80125]
        assertEquals(0.1, B.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.2, B.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(0, 2), LOOSE_MID_TOL);
        assertEquals(0.425, B.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.06875, B.get(1, 1), LOOSE_MID_TOL);
        assertEquals(0.15625, B.get(1, 2), LOOSE_MID_TOL);
        assertEquals(0.141, B.get(2, 0), LOOSE_MID_TOL);
        assertEquals(0.01975, B.get(2, 1), LOOSE_MID_TOL);
        assertEquals(0.80125, B.get(2, 2), LOOSE_MID_TOL);

        // Verify result is valid DPH
        boolean flag = CheckDPHRepresentation.checkDPHRepresentation(beta, B, 1e-12);
        assertTrue(flag, "DPHFromMG result should be valid DPH");
    }

    // ============ AcyclicDPHFromMG - From butools.tmp/AcyclicDPHFromMG.txt ============

    /**
     * AcyclicDPHFromMG.txt example:
     * Input: a=[0,0,1.0], A=[0.22,0,0; 0.3,0.1,0.55; 0.26,0,0.73]
     * Expected: b=[0.69103, 0.29786, 0.011111], B=[0.73,0.27,0; 0,0.22,0.78; 0,0,0.1]
     * Verify moments match: ma = MomentsFromMG(a,A,5) = [4.9383, 34.807, 339.49, 4335.8, 68954]
     */
    @Test
    public void testAcyclicDPHFromMG() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.0);
        alpha.set(0, 1, 0.0);
        alpha.set(0, 2, 1.0);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.22); A.set(0, 1, 0.0);  A.set(0, 2, 0.0);
        A.set(1, 0, 0.3);  A.set(1, 1, 0.1);  A.set(1, 2, 0.55);
        A.set(2, 0, 0.26); A.set(2, 1, 0.0);  A.set(2, 2, 0.73);

        MGRepresentation result = AcyclicDPHFromMG.acyclicDPHFromMG(alpha, A, 1e-14);
        assertNotNull(result);

        Matrix beta = result.getAlpha();
        Matrix B = result.getA();

        // Acyclic DPH representation is not unique, so verify moments match instead of element values
        double[] ma = MomentsFromMG.momentsFromMG(alpha, A, 5);
        double[] mb = MomentsFromMG.momentsFromMG(beta, B, 5);

        double[] expectedMoms = {4.9383, 34.807, 339.49, 4335.8, 68954.0};
        for (int i = 0; i < 5; i++) {
            assertEquals(expectedMoms[i], ma[i], expectedMoms[i] * 1e-3,
                    "Original moment " + (i + 1) + " mismatch");
            assertEquals(ma[i], mb[i], Math.abs(ma[i]) * 1e-6,
                    "Acyclic moment " + (i + 1) + " should match original");
        }

        // Verify result is valid DPH
        boolean flag = CheckDPHRepresentation.checkDPHRepresentation(beta, B, 1e-12);
        assertTrue(flag, "AcyclicDPHFromMG result should be valid DPH");
    }

    // ============ MomentsFromMG - From butools.tmp/MomentsFromMG.txt ============

    /**
     * MomentsFromMG.txt example:
     * Input: a=[-0.6,0.3,1.3], A=[0.25,0.2,-0.15; 0.3,0.1,0.25; 0,0.2,0.47], K=5
     * Expected: [3.4675, 16.203, 97.729, 731.45, 6576.8]
     */
    @Test
    public void testMomentsFromMG_K5() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, -0.6);
        alpha.set(0, 1, 0.3);
        alpha.set(0, 2, 1.3);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.25);  A.set(0, 1, 0.2);  A.set(0, 2, -0.15);
        A.set(1, 0, 0.3);   A.set(1, 1, 0.1);  A.set(1, 2, 0.25);
        A.set(2, 0, 0.0);   A.set(2, 1, 0.2);  A.set(2, 2, 0.47);

        double[] moms = MomentsFromMG.momentsFromMG(alpha, A, 5);

        assertNotNull(moms);
        assertEquals(5, moms.length);
        assertEquals(3.4675, moms[0], LOOSE_MID_TOL);
        assertEquals(16.203, moms[1], COARSE_TOL);
        assertEquals(97.729, moms[2], VERY_COARSE_TOL);
        assertEquals(731.45, moms[3], 1.0);
        assertEquals(6576.8, moms[4], 1.0);
    }

    /**
     * MomentsFromMG.txt example with K=3:
     * Input: a=[-0.6,0.3,1.3], A=[0.25,0.2,-0.15; 0.3,0.1,0.25; 0,0.2,0.47], K=3
     * Expected: [3.4675, 16.203, 97.729]
     */
    @Test
    public void testMomentsFromMG_K3() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, -0.6);
        alpha.set(0, 1, 0.3);
        alpha.set(0, 2, 1.3);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.25);  A.set(0, 1, 0.2);  A.set(0, 2, -0.15);
        A.set(1, 0, 0.3);   A.set(1, 1, 0.1);  A.set(1, 2, 0.25);
        A.set(2, 0, 0.0);   A.set(2, 1, 0.2);  A.set(2, 2, 0.47);

        double[] moms = MomentsFromMG.momentsFromMG(alpha, A, 3);

        assertNotNull(moms);
        assertEquals(3, moms.length);
        assertEquals(3.4675, moms[0], LOOSE_MID_TOL);
        assertEquals(16.203, moms[1], COARSE_TOL);
        assertEquals(97.729, moms[2], VERY_COARSE_TOL);
    }

    // ============ CheckMGRepresentation - From butools.tmp/CheckMGRepresentation.txt ============

    /**
     * CheckMGRepresentation.txt Test 1: Valid MG representation
     * Input: a=[-0.6,0.3,1.3], A=[0.25,0.2,-0.15; 0.3,0.1,0.25; 0,0.2,0.47]
     * Expected: flag=true
     */
    @Test
    public void testCheckMGRepresentation_Valid() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, -0.6);
        alpha.set(0, 1, 0.3);
        alpha.set(0, 2, 1.3);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.25);  A.set(0, 1, 0.2);  A.set(0, 2, -0.15);
        A.set(1, 0, 0.3);   A.set(1, 1, 0.1);  A.set(1, 2, 0.25);
        A.set(2, 0, 0.0);   A.set(2, 1, 0.2);  A.set(2, 2, 0.47);

        boolean flag = CheckMGRepresentation.checkMGRepresentation(alpha, A, 1e-14);
        assertTrue(flag, "Valid MG representation should pass check");
    }

    /**
     * CheckMGRepresentation.txt Test 2: Invalid MG - largest eigenvalue is complex
     * Input: a=[-0.6,0.3,1.3], A=[0.35,0.2,-0.25; 0.3,0.1,0.25; 0,0.2,0.47]
     * Expected: flag=false
     */
    @Test
    public void testCheckMGRepresentation_Invalid() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, -0.6);
        alpha.set(0, 1, 0.3);
        alpha.set(0, 2, 1.3);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.35);  A.set(0, 1, 0.2);  A.set(0, 2, -0.25);
        A.set(1, 0, 0.3);   A.set(1, 1, 0.1);  A.set(1, 2, 0.25);
        A.set(2, 0, 0.0);   A.set(2, 1, 0.2);  A.set(2, 2, 0.47);

        boolean flag = CheckMGRepresentation.checkMGRepresentation(alpha, A, 1e-14);
        assertFalse(flag, "Invalid MG (complex largest eigenvalue) should fail check");
    }

    // ============ MGFromMoments - From butools.tmp/MGFromMoments.txt ============

    /**
     * MGFromMoments.txt example:
     * Input: moms=[4.08, 20.41, 130.45, 1054.41, 10463.73]
     * Expected: a=[0.33333, 0.33333, 0.33333],
     *           A=[0.15523,1.7289,0.10482; -0.013774,0.6823,-0.023472; -0.013847,-0.16787,0.82688]
     * Verify moments roundtrip.
     */
    @Test
    public void testMGFromMoments() {
        double[] moms = {4.08, 20.41, 130.45, 1054.41, 10463.73};

        MGRepresentation result = MGFromMoments.mgFromMoments(moms);
        assertNotNull(result);

        Matrix alpha = result.getAlpha();
        Matrix A = result.getA();

        // Verify alpha = [0.33333, 0.33333, 0.33333]
        assertEquals(0.33333, alpha.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.33333, alpha.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.33333, alpha.get(0, 2), LOOSE_MID_TOL);

        // Verify A = [0.15523,1.7289,0.10482; -0.013774,0.6823,-0.023472; -0.013847,-0.16787,0.82688]
        assertEquals(0.15523, A.get(0, 0), LOOSE_MID_TOL);
        assertEquals(1.7289, A.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.10482, A.get(0, 2), LOOSE_MID_TOL);
        assertEquals(-0.013774, A.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.6823, A.get(1, 1), LOOSE_MID_TOL);
        assertEquals(-0.023472, A.get(1, 2), LOOSE_MID_TOL);
        assertEquals(-0.013847, A.get(2, 0), LOOSE_MID_TOL);
        assertEquals(-0.16787, A.get(2, 1), LOOSE_MID_TOL);
        assertEquals(0.82688, A.get(2, 2), LOOSE_MID_TOL);

        // Verify moments roundtrip
        double[] memoms = MomentsFromMG.momentsFromMG(alpha, A, 5);
        assertEquals(4.08, memoms[0], LOOSE_MID_TOL);
        assertEquals(20.41, memoms[1], COARSE_TOL);
        assertEquals(130.45, memoms[2], VERY_COARSE_TOL);
        assertEquals(1054.4, memoms[3], 1.0);
        assertEquals(10464.0, memoms[4], 10.0);
    }

    // ============ Random Generator and Sampling Tests ============

    /**
     * RandomDPH.txt: Generate random DPH with order=3, mean=10, zeroEntries=5
     * Verify: valid DPH representation, correct mean
     */
    @Test
    public void testRandomDPH() {
        MGRepresentation result = RandomDPH.randomDPH(3, 10.0, 5, 1000, 1e-7, new Random(42));
        Matrix alpha = result.getAlpha();
        Matrix A = result.getA();

        // Verify dimensions
        assertEquals(1, alpha.getNumRows());
        assertEquals(3, alpha.getNumCols());
        assertEquals(3, A.getNumRows());
        assertEquals(3, A.getNumCols());

        // Verify valid DPH representation
        assertTrue(CheckDPHRepresentation.checkDPHRepresentation(alpha, A, 1e-14),
                "RandomDPH should produce valid DPH");

        // Verify mean
        double[] moms = MomentsFromDPH.momentsFromDPH(alpha, A, 1);
        assertEquals(10.0, moms[0], COARSE_TOL, "Mean should be close to 10.0");
    }

    /**
     * SamplesFromDPH.txt: Generate samples from DPH, verify moments match theoretical
     * Uses DPH from AcyclicDPHFromMG test data
     */
    @Test
    public void testSamplesFromDPH() {
        // Use a known DPH: order 3, mean ~4.08
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.76); alpha.set(0, 1, 0.0); alpha.set(0, 2, 0.24);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.0); A.set(0, 1, 0.3); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0); A.set(1, 1, 0.1); A.set(1, 2, 0.25);
        A.set(2, 0, 0.3); A.set(2, 1, 0.0); A.set(2, 2, 0.2);

        // Theoretical moments
        double[] theoretical = MomentsFromDPH.momentsFromDPH(alpha, A, 3);

        // Generate 10000 samples
        int[] samples = SamplesFromDPH.samplesFromDPH(alpha, A, 10000, new Random(42));
        assertEquals(10000, samples.length);

        // Convert to double for trace analysis
        double[] dsamples = new double[samples.length];
        for (int i = 0; i < samples.length; i++) {
            dsamples[i] = samples[i];
        }

        // Compute trace moments
        double[] moms = MarginalMomentsFromTrace.marginalMomentsFromTrace(dsamples, 3);

        // Stochastic: trace moments should be within 20% of theoretical
        for (int i = 0; i < 3; i++) {
            assertEquals(theoretical[i], moms[i], Math.abs(theoretical[i]) * 0.2,
                    "SamplesFromDPH moms[" + i + "]");
        }
    }

    // ============ Phase 1: CDF/PMF Tests for DPH and MG ============

    /**
     * CdfFromDPH: Verify CDF at specific points.
     * Input: a=[0.76,0,0.24], A=[0.34,0.66,0; 0.79,0.05,0.07; 0.26,0.73,0.01]
     * CDF(0) = 1 - sum(alpha), CDF is monotonically non-decreasing, CDF converges to 1.
     */
    @Test
    public void testCdfFromDPH() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.76);
        alpha.set(0, 1, 0.0);
        alpha.set(0, 2, 0.24);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.34); A.set(0, 1, 0.66); A.set(0, 2, 0.0);
        A.set(1, 0, 0.79); A.set(1, 1, 0.05); A.set(1, 2, 0.07);
        A.set(2, 0, 0.26); A.set(2, 1, 0.73); A.set(2, 2, 0.01);

        int[] x = {0, 1, 5, 10, 50, 100, 500};
        double[] cdf = CdfFromDPH.cdfFromDPH(alpha, A, x);

        assertNotNull(cdf);
        assertEquals(x.length, cdf.length);

        // CDF(0) = 1 - sum(alpha) = 1 - 1.0 = 0
        assertEquals(0.0, cdf[0], LOOSE_MID_TOL, "CDF(0) should be 1 - sum(alpha)");

        // CDF should be monotonically non-decreasing
        for (int i = 1; i < cdf.length; i++) {
            assertTrue(cdf[i] >= cdf[i - 1] - LOOSE_MID_TOL,
                    "CDF should be non-decreasing at x=" + x[i]);
        }

        // CDF values should be in [0, 1]
        for (int i = 0; i < cdf.length; i++) {
            assertTrue(cdf[i] >= -LOOSE_MID_TOL, "CDF should be >= 0 at x=" + x[i]);
            assertTrue(cdf[i] <= 1.0 + LOOSE_MID_TOL, "CDF should be <= 1 at x=" + x[i]);
        }

        // CDF(500) should be very close to 1 (mean ~27)
        assertTrue(cdf[6] > 0.999, "CDF(500) should be close to 1 for mean ~27");
    }

    /**
     * CdfFromMG: Verify CDF at specific points.
     * Input: a=[-0.6,0.3,1.3], A=[0.25,0.2,-0.15; 0.3,0.1,0.25; 0,0.2,0.47]
     */
    @Test
    public void testCdfFromMG() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, -0.6);
        alpha.set(0, 1, 0.3);
        alpha.set(0, 2, 1.3);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.25);  A.set(0, 1, 0.2);  A.set(0, 2, -0.15);
        A.set(1, 0, 0.3);   A.set(1, 1, 0.1);  A.set(1, 2, 0.25);
        A.set(2, 0, 0.0);   A.set(2, 1, 0.2);  A.set(2, 2, 0.47);

        int[] x = {0, 1, 5, 10, 50, 100};
        double[] cdf = CdfFromMG.cdfFromMG(alpha, A, x);

        assertNotNull(cdf);
        assertEquals(x.length, cdf.length);

        // CDF(0) = 1 - sum(alpha) = 1 - 1.0 = 0
        assertEquals(0.0, cdf[0], LOOSE_MID_TOL, "CDF(0) should be 0");

        // CDF should be monotonically non-decreasing
        for (int i = 1; i < cdf.length; i++) {
            assertTrue(cdf[i] >= cdf[i - 1] - LOOSE_MID_TOL,
                    "CDF should be non-decreasing at x=" + x[i]);
        }

        // CDF(100) should be close to 1 (mean ~3.47)
        assertTrue(cdf[5] > 0.99, "CDF(100) should be close to 1 for mean ~3.47");
    }

    /**
     * PmfFromDPH: Verify PMF is non-negative and cumsum matches CDF.
     * Input: a=[0.76,0,0.24], A=[0.34,0.66,0; 0.79,0.05,0.07; 0.26,0.73,0.01]
     */
    @Test
    public void testPmfFromDPH() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.76);
        alpha.set(0, 1, 0.0);
        alpha.set(0, 2, 0.24);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.34); A.set(0, 1, 0.66); A.set(0, 2, 0.0);
        A.set(1, 0, 0.79); A.set(1, 1, 0.05); A.set(1, 2, 0.07);
        A.set(2, 0, 0.26); A.set(2, 1, 0.73); A.set(2, 2, 0.01);

        // Compute PMF at x = 0..200
        int nPoints = 201;
        int[] x = new int[nPoints];
        for (int i = 0; i < nPoints; i++) {
            x[i] = i;
        }

        double[] pmf = PmfFromDPH.pmfFromDPH(alpha, A, x);

        assertNotNull(pmf);
        assertEquals(nPoints, pmf.length);

        // PMF should be non-negative
        for (int i = 0; i < nPoints; i++) {
            assertTrue(pmf[i] >= -LOOSE_MID_TOL, "PMF should be >= 0 at x=" + x[i]);
        }

        // Cumulative sum of PMF should approximate CDF
        double[] cumsum = new double[nPoints];
        cumsum[0] = pmf[0];
        for (int i = 1; i < nPoints; i++) {
            cumsum[i] = cumsum[i - 1] + pmf[i];
        }

        double[] cdf = CdfFromDPH.cdfFromDPH(alpha, A, x);
        for (int i = 0; i < nPoints; i++) {
            assertEquals(cdf[i], cumsum[i], LOOSE_MID_TOL,
                    "Cumsum of PMF should match CDF at x=" + x[i]);
        }

        // Total probability should be close to 1
        assertTrue(cumsum[nPoints - 1] > 0.99,
                "Total PMF over [0,200] should be close to 1 for mean ~27");
    }

    /**
     * PmfFromMG: Verify PMF non-negative and cumsum matches CDF.
     * Input: a=[-0.6,0.3,1.3], A=[0.25,0.2,-0.15; 0.3,0.1,0.25; 0,0.2,0.47]
     */
    @Test
    public void testPmfFromMG() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, -0.6);
        alpha.set(0, 1, 0.3);
        alpha.set(0, 2, 1.3);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, 0.25);  A.set(0, 1, 0.2);  A.set(0, 2, -0.15);
        A.set(1, 0, 0.3);   A.set(1, 1, 0.1);  A.set(1, 2, 0.25);
        A.set(2, 0, 0.0);   A.set(2, 1, 0.2);  A.set(2, 2, 0.47);

        // Compute PMF at x = 0..100
        int nPoints = 101;
        int[] x = new int[nPoints];
        for (int i = 0; i < nPoints; i++) {
            x[i] = i;
        }

        double[] pmf = PmfFromMG.pmfFromMG(alpha, A, x);

        assertNotNull(pmf);
        assertEquals(nPoints, pmf.length);

        // Cumulative sum of PMF should approximate CDF
        double[] cumsum = new double[nPoints];
        cumsum[0] = pmf[0];
        for (int i = 1; i < nPoints; i++) {
            cumsum[i] = cumsum[i - 1] + pmf[i];
        }

        double[] cdf = CdfFromMG.cdfFromMG(alpha, A, x);
        for (int i = 0; i < nPoints; i++) {
            assertEquals(cdf[i], cumsum[i], LOOSE_MID_TOL,
                    "Cumsum of PMF should match CDF at x=" + x[i]);
        }

        // Total probability should be close to 1
        assertTrue(cumsum[nPoints - 1] > 0.99,
                "Total PMF over [0,100] should be close to 1 for mean ~3.47");
    }
}
