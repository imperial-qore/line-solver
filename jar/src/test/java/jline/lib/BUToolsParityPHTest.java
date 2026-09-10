package jline.lib;

import jline.lib.butools.ph.*;
import jline.lib.butools.ph.PH2From3Moments.PH2Representation;
import jline.lib.butools.APH2ndMomentLowerBound;
import jline.lib.butools.APH3rdMomentLowerBound;
import jline.lib.butools.APH3rdMomentUpperBound;
import jline.lib.butools.trace.MarginalMomentsFromTrace;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.*;

import java.util.Random;

/**
 * BUTools PH Parity Tests - Unit tests validating JAR PH implementations against butools.tmp examples.
 *
 * Each test uses EXACT data from butools.tmp example files to ensure the JAR implementation
 * produces identical results to the original BUTools library.
 */
public class BUToolsParityPHTest {


    // ============ MomentsFromPH - From butools.tmp/MomentsFromPH.txt ============

    /**
     * MomentsFromPH.txt example: Compute 5 moments from PH distribution
     * Input: a=[0.1,0.9,0], A=[-6.2,2,0; 2,-9,1; 1,0,-3]
     * Expected: [0.20939, 0.10449, 0.089092, 0.11027, 0.17953]
     */
    @Test
    public void testMomentsFromPH() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.1);
        alpha.set(0, 1, 0.9);
        alpha.set(0, 2, 0.0);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -6.2); A.set(0, 1, 2.0); A.set(0, 2, 0.0);
        A.set(1, 0, 2.0);  A.set(1, 1, -9.0); A.set(1, 2, 1.0);
        A.set(2, 0, 1.0);  A.set(2, 1, 0.0); A.set(2, 2, -3.0);

        double[] moms = MomentsFromME.momentsFromPH(alpha, A, 5);

        assertNotNull(moms);
        assertEquals(5, moms.length);
        assertEquals(0.20939, moms[0], LOOSE_MID_TOL);
        assertEquals(0.10449, moms[1], LOOSE_MID_TOL);
        assertEquals(0.089092, moms[2], LOOSE_MID_TOL);
        assertEquals(0.11027, moms[3], LOOSE_MID_TOL);
        assertEquals(0.17953, moms[4], LOOSE_MID_TOL);
    }

    // ============ MomentsFromME - From butools.tmp/MomentsFromME.txt ============

    /**
     * MomentsFromME.txt example: Compute moments from ME distribution
     * Input: a=[0.2,0.3,0.5], A=[-1,0,0; 0,-3,2; 0,-2,-3], K=5
     * Expected (first 5): [0.35385, 0.41893, 1.1552, 4.6998, 23.838]
     */
    @Test
    public void testMomentsFromME() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2);
        alpha.set(0, 1, 0.3);
        alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        double[] moms = MomentsFromME.momentsFromME(alpha, A, 5);

        assertNotNull(moms);
        assertEquals(5, moms.length);
        assertEquals(0.35385, moms[0], LOOSE_MID_TOL);
        assertEquals(0.41893, moms[1], LOOSE_MID_TOL);
        assertEquals(1.1552, moms[2], LOOSE_MID_TOL);
        assertEquals(4.6998, moms[3], LOOSE_MID_TOL);
        assertEquals(23.838, moms[4], LOOSE_MID_TOL);
    }

    /**
     * MomentsFromME.txt example: Compute 9 moments from ME distribution
     * Input: a=[0.2,0.3,0.5], A=[-1,0,0; 0,-3,2; 0,-2,-3], K=9
     * Expected: [0.35385, 0.41893, 1.1552, 4.6998, 23.838, 143.78, 1007.8, 8064.3, 72578]
     */
    @Test
    public void testMomentsFromME_9Moments() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2);
        alpha.set(0, 1, 0.3);
        alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        double[] moms = MomentsFromME.momentsFromME(alpha, A, 9);

        assertNotNull(moms);
        assertEquals(9, moms.length);
        double[] expected = {0.35385, 0.41893, 1.1552, 4.6998, 23.838, 143.78, 1007.8, 8064.3, 72578.0};
        for (int i = 0; i < 9; i++) {
            assertEquals(expected[i], moms[i], Math.max(LOOSE_MID_TOL, Math.abs(expected[i]) * 1e-3),
                    "moms[" + i + "]");
        }
    }

    // ============ PH3From5Moments - From butools.tmp/PH3From5Moments.txt ============

    /**
     * PH3From5Moments.txt Example 1: Construct PH3 from 5 moments
     * Input: moms=[0.20939, 0.10449, 0.089092, 0.11027, 0.17953]
     * Expected: alpha=[0.58305, 0.32736, 0.089589]
     *           A=[-9.9819,0,0; 5.3405,-5.3405,0; 0,2.8776,-2.8776]
     */
    @Test
    public void testPH3From5Moments_Example1() {
        double[] moms = {0.20939, 0.10449, 0.089092, 0.11027, 0.17953};

        PH3Representation ph3 = PH3From5Moments.ph3From5Moments(moms, FINE_TOL);

        assertNotNull(ph3);
        Matrix alpha = ph3.getAlpha();
        Matrix A = ph3.getA();
        assertNotNull(alpha);
        assertNotNull(A);

        // PH3 representation is not unique, so verify moments roundtrip instead of element values
        double[] phmoms = MomentsFromME.momentsFromME(alpha, A, 5);
        for (int i = 0; i < 5; i++) {
            assertEquals(moms[i], phmoms[i], Math.abs(moms[i]) * 1e-6, "Moment " + (i + 1) + " roundtrip mismatch");
        }

        // Verify alpha sums to 1
        double alphaSum = 0;
        for (int i = 0; i < alpha.getNumCols(); i++) {
            alphaSum += alpha.get(0, i);
        }
        assertEquals(1.0, alphaSum, LOOSE_MID_TOL, "Alpha should sum to 1");
    }

    /**
     * PH3From5Moments.txt Example 2: Construct PH3 from 5 moments (different input)
     * Input: moms from a=[0.2,0.3,0.5], A=[-1,0,0; 0,-3,0.5; 0,-0.5,-3]
     *        moms=[0.44865, 0.5496, 1.3298, 4.9428, 24.182]
     * Expected: alpha=[0.94865, 0.036778, 0.014574]
     *           A=[-3,0,0.15385; 2.866,-2.866,0; 0,1.134,-1.134]
     */
    @Test
    public void testPH3From5Moments_Example2() {
        // First compute moments from the given ME
        Matrix alpha_in = new Matrix(1, 3);
        alpha_in.set(0, 0, 0.2);
        alpha_in.set(0, 1, 0.3);
        alpha_in.set(0, 2, 0.5);

        Matrix A_in = new Matrix(3, 3);
        A_in.set(0, 0, -1.0); A_in.set(0, 1, 0.0);  A_in.set(0, 2, 0.0);
        A_in.set(1, 0, 0.0);  A_in.set(1, 1, -3.0);  A_in.set(1, 2, 0.5);
        A_in.set(2, 0, 0.0);  A_in.set(2, 1, -0.5);  A_in.set(2, 2, -3.0);

        double[] moms = MomentsFromME.momentsFromME(alpha_in, A_in, 5);

        // Verify input moments match expected
        assertEquals(0.44865, moms[0], LOOSE_MID_TOL);
        assertEquals(0.5496, moms[1], LOOSE_MID_TOL);
        assertEquals(1.3298, moms[2], LOOSE_MID_TOL);
        assertEquals(4.9428, moms[3], LOOSE_MID_TOL);
        assertEquals(24.182, moms[4], LOOSE_MID_TOL);

        // Construct PH3 from moments
        PH3Representation ph3 = PH3From5Moments.ph3From5Moments(moms, FINE_TOL);

        assertNotNull(ph3);
        Matrix alpha = ph3.getAlpha();
        Matrix A = ph3.getA();

        // PH3 representation is not unique, so verify moments roundtrip instead of element values
        double[] phmoms = MomentsFromME.momentsFromME(alpha, A, 5);
        for (int i = 0; i < 5; i++) {
            assertEquals(moms[i], phmoms[i], Math.abs(moms[i]) * 1e-6, "Moment " + (i + 1) + " roundtrip mismatch");
        }

        // Verify alpha sums to 1
        double alphaSum = 0;
        for (int i = 0; i < alpha.getNumCols(); i++) {
            alphaSum += alpha.get(0, i);
        }
        assertEquals(1.0, alphaSum, LOOSE_MID_TOL, "Alpha should sum to 1");
    }

    // ============ CanonicalFromPH2 - From butools.tmp/CanonicalFromPH2.txt ============

    /**
     * CanonicalFromPH2.txt example: Transform PH2 to canonical form
     * Input: a=[0.12,0.88], A=[-1.28,0; 3.94,-3.94]
     * Expected: b=[0.96102, 0.038985], B=[-1.28,1.28; 0,-3.94]
     */
    @Test
    public void testCanonicalFromPH2() {
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, 0.12);
        alpha.set(0, 1, 0.88);

        Matrix A = new Matrix(2, 2);
        A.set(0, 0, -1.28); A.set(0, 1, 0.0);
        A.set(1, 0, 3.94);  A.set(1, 1, -3.94);

        PH2Representation result = CanonicalFromPH2.canonicalFromPH2(alpha, A, 1e-14);

        assertNotNull(result);
        Matrix b = result.alpha;
        Matrix B = result.A;

        // Verify b
        assertEquals(0.96102, b.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.038985, b.get(0, 1), LOOSE_MID_TOL);

        // Verify B
        assertEquals(-1.28, B.get(0, 0), LOOSE_MID_TOL);
        assertEquals(1.28, B.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(1, 0), LOOSE_MID_TOL);
        assertEquals(-3.94, B.get(1, 1), LOOSE_MID_TOL);
    }

    // ============ CanonicalFromPH3 - From butools.tmp/CanonicalFromPH3.txt ============

    /**
     * CanonicalFromPH3.txt example: Transform PH3 to canonical form
     * Input: a=[0.1,0.9,0], A=[-6.2,2,0; 2,-9,1; 1,0,-3]
     * Expected: b=[0.58305, 0.32736, 0.089589]
     *           B=[-9.9819,0,0; 5.3405,-5.3405,0; 0,2.8776,-2.8776]
     */
    @Test
    public void testCanonicalFromPH3() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.1);
        alpha.set(0, 1, 0.9);
        alpha.set(0, 2, 0.0);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -6.2); A.set(0, 1, 2.0); A.set(0, 2, 0.0);
        A.set(1, 0, 2.0);  A.set(1, 1, -9.0); A.set(1, 2, 1.0);
        A.set(2, 0, 1.0);  A.set(2, 1, 0.0); A.set(2, 2, -3.0);

        PH3Representation result = CanonicalFromPH3.canonicalFromPH3(alpha, A, FINE_TOL);

        assertNotNull(result);
        Matrix b = result.getAlpha();
        Matrix B = result.getA();

        // Verify b
        assertEquals(0.58305, b.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.32736, b.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.089589, b.get(0, 2), LOOSE_MID_TOL);

        // Verify B
        assertEquals(-9.9819, B.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(0, 2), LOOSE_MID_TOL);
        assertEquals(5.3405, B.get(1, 0), LOOSE_MID_TOL);
        assertEquals(-5.3405, B.get(1, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(1, 2), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(2, 0), LOOSE_MID_TOL);
        assertEquals(2.8776, B.get(2, 1), LOOSE_MID_TOL);
        assertEquals(-2.8776, B.get(2, 2), LOOSE_MID_TOL);
    }

    // ============ APH2ndMomentLowerBound - From butools.tmp/APH2ndMomentLowerBound.txt ============

    /**
     * APH2ndMomentLowerBound.txt example: Compute lower bound on 2nd moment
     * Input: mean=1.9, n=4
     * Expected: mom2=4.5125 (and 1/cv2=4)
     */
    @Test
    public void testAPH2ndMomentLowerBound() {
        double mom2 = APH2ndMomentLowerBound.APH2ndMomentLowerBound(1.9, 4);
        assertEquals(4.5125, mom2, LOOSE_MID_TOL);

        // Verify 1/cv2 = 4
        double cv2 = mom2 / (1.9 * 1.9) - 1.0;
        assertEquals(4.0, 1.0 / cv2, LOOSE_MID_TOL);
    }

    // ============ APH3rdMomentLowerBound - From butools.tmp/APH3rdMomentLowerBound.txt ============

    /**
     * APH3rdMomentLowerBound.txt example 1: mean=1.9, mom2=5, n=3
     * Expected: 16.577
     */
    @Test
    public void testAPH3rdMomentLowerBound_n3() {
        double mom3lower = APH3rdMomentLowerBound.APH3rdMomentLowerBound(1.9, 5.0, 3);
        assertEquals(16.577, mom3lower, LOOSE_MID_TOL);
    }

    /**
     * APH3rdMomentLowerBound.txt example 2: mean=1.9, mom2=5, n=4
     * Expected: 16.079
     */
    @Test
    public void testAPH3rdMomentLowerBound_n4() {
        double mom3lower = APH3rdMomentLowerBound.APH3rdMomentLowerBound(1.9, 5.0, 4);
        assertEquals(16.079, mom3lower, LOOSE_MID_TOL);
    }

    // ============ APH3rdMomentUpperBound - From butools.tmp/APH3rdMomentUpperBound.txt ============

    /**
     * APH3rdMomentUpperBound.txt example 1: mean=1.9, mom2=5, n=3
     * Expected: 17.081
     */
    @Test
    public void testAPH3rdMomentUpperBound_n3() {
        double mom3upper = APH3rdMomentUpperBound.APH3rdMomentUpperBound(1.9, 5.0, 3);
        assertEquals(17.081, mom3upper, LOOSE_MID_TOL);
    }

    /**
     * APH3rdMomentUpperBound.txt example 2: mean=1.9, mom2=5, n=4
     * Expected: Inf
     */
    @Test
    public void testAPH3rdMomentUpperBound_n4_Inf() {
        double mom3upper = APH3rdMomentUpperBound.APH3rdMomentUpperBound(1.9, 5.0, 4);
        assertTrue(Double.isInfinite(mom3upper), "Expected Inf for n=4");
    }

    // ============ CheckMERepresentation - From butools.tmp/CheckMERepresentation.txt ============

    /**
     * CheckMERepresentation.txt Example 1: Invalid - eigenvalue with non-negative real part
     * Input: a=[-0.2,0.2], A=[1,-1; 1,-2]
     * Expected: flag=0 (false)
     */
    @Test
    public void testCheckMERepresentation_InvalidEigenvalue() {
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, -0.2);
        alpha.set(0, 1, 0.2);

        Matrix A = new Matrix(2, 2);
        A.set(0, 0, 1.0);  A.set(0, 1, -1.0);
        A.set(1, 0, 1.0);  A.set(1, 1, -2.0);

        boolean flag = CheckMERepresentation.checkMERepresentation(alpha, A, 1e-12);
        assertFalse(flag, "Should fail: eigenvalue with non-negative real part");
    }

    /**
     * CheckMERepresentation.txt Example 2: Invalid - dominant eigenvalue not real
     * Input: a=[-0.2,0.4,0.8], A=[-2,0,3; 0,-1,1; 0,-1,-1]
     * Expected: flag=0 (false)
     */
    @Test
    public void testCheckMERepresentation_InvalidDominantEigenvalue() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, -0.2);
        alpha.set(0, 1, 0.4);
        alpha.set(0, 2, 0.8);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -2.0); A.set(0, 1, 0.0);  A.set(0, 2, 3.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -1.0);  A.set(1, 2, 1.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -1.0);  A.set(2, 2, -1.0);

        boolean flag = CheckMERepresentation.checkMERepresentation(alpha, A, 1e-12);
        assertFalse(flag, "Should fail: dominant eigenvalue not real");
    }

    /**
     * CheckMERepresentation.txt Example 3: Valid ME representation
     * Input: a=[0.2,0.3,0.5], A=[-1,0,0; 0,-3,2; 0,-2,-3]
     * Expected: flag=1 (true)
     */
    @Test
    public void testCheckMERepresentation_Valid() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2);
        alpha.set(0, 1, 0.3);
        alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0);  A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0);  A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0);  A.set(2, 2, -3.0);

        boolean flag = CheckMERepresentation.checkMERepresentation(alpha, A, 1e-12);
        assertTrue(flag, "Valid ME representation should pass check");
    }

    // ============ CheckRAPRepresentation - From butools.tmp/CheckRAPRepresentation.txt ============

    /**
     * CheckRAPRepresentation.txt Example 1: Invalid - D0 is not quadratic
     * Input: H0=[4x3 matrix], H1=[4x3 matrix]
     * Expected: flag=0 (false)
     */
    @Test
    public void testCheckRAPRepresentation_NotQuadratic() {
        Matrix H0 = new Matrix(4, 3);
        H0.set(0, 0, -1.0); H0.set(0, 1, 0.0);  H0.set(0, 2, 1.0);
        H0.set(1, 0, 0.0);  H0.set(1, 1, -2.0);  H0.set(1, 2, 0.0);
        H0.set(2, 0, 1.0);  H0.set(2, 1, 0.0);   H0.set(2, 2, -3.0);
        H0.set(3, 0, 1.0);  H0.set(3, 1, 2.0);   H0.set(3, 2, 2.0);

        Matrix H1 = new Matrix(4, 3);
        H1.set(0, 0, -1.0); H1.set(0, 1, 0.0);  H1.set(0, 2, 1.0);
        H1.set(1, 0, 0.0);  H1.set(1, 1, -2.0);  H1.set(1, 2, 0.0);
        H1.set(2, 0, 1.0);  H1.set(2, 1, 0.0);   H1.set(2, 2, -3.0);
        H1.set(3, 0, 1.0);  H1.set(3, 1, 2.0);   H1.set(3, 2, 2.0);

        boolean flag = CheckRAPRepresentation.checkRAPRepresentation(H0, H1, 1e-12);
        assertFalse(flag, "Should fail: D0 is not quadratic");
    }

    /**
     * CheckRAPRepresentation.txt Example 2: Invalid - rowsum of D0+D1 is not 0
     * Input: H0=[-1,0,2; 0,2,0; 1,0,-3], H1=[-1,0,1; 0,-2,0; 1,0,-3]
     * Expected: flag=0 (false)
     */
    @Test
    public void testCheckRAPRepresentation_RowsumNotZero() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, -1.0); H0.set(0, 1, 0.0);  H0.set(0, 2, 2.0);
        H0.set(1, 0, 0.0);  H0.set(1, 1, 2.0);   H0.set(1, 2, 0.0);
        H0.set(2, 0, 1.0);  H0.set(2, 1, 0.0);   H0.set(2, 2, -3.0);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, -1.0); H1.set(0, 1, 0.0);  H1.set(0, 2, 1.0);
        H1.set(1, 0, 0.0);  H1.set(1, 1, -2.0);  H1.set(1, 2, 0.0);
        H1.set(2, 0, 1.0);  H1.set(2, 1, 0.0);   H1.set(2, 2, -3.0);

        boolean flag = CheckRAPRepresentation.checkRAPRepresentation(H0, H1, 1e-12);
        assertFalse(flag, "Should fail: rowsum of D0+D1 is not 0");
    }

    /**
     * CheckRAPRepresentation.txt Example 3: Invalid - eigenvalue of D0 with non-negative real part
     * Input: H0=[-1,0,0; 0,-2,2; 0,3,-3], H1=[0,0,1; 0,-1,1; 1,0,-1]
     * Expected: flag=0 (false)
     */
    @Test
    public void testCheckRAPRepresentation_NonNegativeEigenvalue() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, -1.0); H0.set(0, 1, 0.0);  H0.set(0, 2, 0.0);
        H0.set(1, 0, 0.0);  H0.set(1, 1, -2.0);  H0.set(1, 2, 2.0);
        H0.set(2, 0, 0.0);  H0.set(2, 1, 3.0);   H0.set(2, 2, -3.0);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 0.0);  H1.set(0, 1, 0.0);  H1.set(0, 2, 1.0);
        H1.set(1, 0, 0.0);  H1.set(1, 1, -1.0);  H1.set(1, 2, 1.0);
        H1.set(2, 0, 1.0);  H1.set(2, 1, 0.0);   H1.set(2, 2, -1.0);

        boolean flag = CheckRAPRepresentation.checkRAPRepresentation(H0, H1, 1e-12);
        assertFalse(flag, "Should fail: eigenvalue of D0 with non-negative real part");
    }

    /**
     * CheckRAPRepresentation.txt Example 4: Invalid - dominant eigenvalue of D0 not real
     * Input: H0=[-2,0,0; 0,-1,1; 0,-1,-1], H1=[1,0,1; 0,1,-1; 1,0,1]
     * Expected: flag=0 (false)
     */
    @Test
    public void testCheckRAPRepresentation_DominantNotReal() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, -2.0); H0.set(0, 1, 0.0);  H0.set(0, 2, 0.0);
        H0.set(1, 0, 0.0);  H0.set(1, 1, -1.0);  H0.set(1, 2, 1.0);
        H0.set(2, 0, 0.0);  H0.set(2, 1, -1.0);  H0.set(2, 2, -1.0);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 1.0);  H1.set(0, 1, 0.0);  H1.set(0, 2, 1.0);
        H1.set(1, 0, 0.0);  H1.set(1, 1, 1.0);   H1.set(1, 2, -1.0);
        H1.set(2, 0, 1.0);  H1.set(2, 1, 0.0);   H1.set(2, 2, 1.0);

        boolean flag = CheckRAPRepresentation.checkRAPRepresentation(H0, H1, 1e-12);
        assertFalse(flag, "Should fail: dominant eigenvalue of D0 not real");
    }

    /**
     * CheckRAPRepresentation.txt Example 5: Valid RAP representation
     * Input: H0=[-1,0,0; 0,-2,1; 0,-1,-2], H1=[1,0,0; 0,1,0; 1,1,1]
     * Expected: flag=1 (true)
     */
    @Test
    public void testCheckRAPRepresentation_Valid() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, -1.0); H0.set(0, 1, 0.0);  H0.set(0, 2, 0.0);
        H0.set(1, 0, 0.0);  H0.set(1, 1, -2.0);  H0.set(1, 2, 1.0);
        H0.set(2, 0, 0.0);  H0.set(2, 1, -1.0);  H0.set(2, 2, -2.0);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 1.0);  H1.set(0, 1, 0.0);  H1.set(0, 2, 0.0);
        H1.set(1, 0, 0.0);  H1.set(1, 1, 1.0);   H1.set(1, 2, 0.0);
        H1.set(2, 0, 1.0);  H1.set(2, 1, 1.0);   H1.set(2, 2, 1.0);

        boolean flag = CheckRAPRepresentation.checkRAPRepresentation(H0, H1, 1e-12);
        assertTrue(flag, "Valid RAP representation should pass check");
    }

    // ============ Random Generator and Sampling Tests ============

    /**
     * RandomPH.txt: Generate random PH with order=3, mean=8, zeroEntries=4
     * Verify: valid PH representation, correct mean
     */
    @Test
    public void testRandomPH() {
        PHRepresentation result = RandomPH.randomPH(3, 8.0, 4, 1000, 1e-7, new Random(42));
        Matrix alpha = result.getAlpha();
        Matrix A = result.getA();

        // Verify dimensions
        assertEquals(1, alpha.getNumRows());
        assertEquals(3, alpha.getNumCols());
        assertEquals(3, A.getNumRows());
        assertEquals(3, A.getNumCols());

        // Verify valid PH representation
        assertTrue(CheckPHRepresentation.checkPHRepresentation(alpha, A, 1e-14),
                "RandomPH should produce valid PH");

        // Verify mean
        double[] moms = MomentsFromME.momentsFromPH(alpha, A, 1);
        assertEquals(8.0, moms[0], COARSE_TOL, "Mean should be close to 8.0");
    }

    /**
     * SamplesFromPH.txt: Generate samples from PH, verify moments match theoretical
     * Input: a=[0.1,0.9,0], A=[-6.2,2,0; 2,-9,1; 1,0,-3]
     * Theoretical moments: [0.20939, 0.10449, 0.089092]
     */
    @Test
    public void testSamplesFromPH() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.1); alpha.set(0, 1, 0.9); alpha.set(0, 2, 0.0);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -6.2); A.set(0, 1, 2.0); A.set(0, 2, 0.0);
        A.set(1, 0, 2.0);  A.set(1, 1, -9.0); A.set(1, 2, 1.0);
        A.set(2, 0, 1.0);  A.set(2, 1, 0.0);  A.set(2, 2, -3.0);

        // Theoretical moments
        double[] theoretical = MomentsFromME.momentsFromPH(alpha, A, 3);

        // Generate 10000 samples
        double[] samples = SamplesFromPH.samplesFromPH(alpha, A, 10000, new Random(42));
        assertEquals(10000, samples.length);

        // Compute trace moments
        double[] moms = MarginalMomentsFromTrace.marginalMomentsFromTrace(samples, 3);

        // Stochastic: trace moments should be within 20% of theoretical
        for (int i = 0; i < 3; i++) {
            assertEquals(theoretical[i], moms[i], Math.abs(theoretical[i]) * 0.2,
                    "SamplesFromPH moms[" + i + "]");
        }
    }

    // ============ Phase 1: CDF/PDF Tests for PH and ME ============

    /**
     * CdfFromPH: Verify CDF values at specific points.
     * Input: a=[0.1,0.9,0], A=[-6.2,2,0; 2,-9,1; 1,0,-3]
     * CDF(0)=0, CDF(inf)→1, CDF is monotonically non-decreasing
     */
    @Test
    public void testCdfFromPH() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.1); alpha.set(0, 1, 0.9); alpha.set(0, 2, 0.0);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -6.2); A.set(0, 1, 2.0); A.set(0, 2, 0.0);
        A.set(1, 0, 2.0);  A.set(1, 1, -9.0); A.set(1, 2, 1.0);
        A.set(2, 0, 1.0);  A.set(2, 1, 0.0); A.set(2, 2, -3.0);

        double[] x = {0.0, 0.1, 0.5, 1.0, 2.0, 3.0};
        double[] cdf = CdfFromME.cdfFromPH(alpha, A, x);

        assertNotNull(cdf);
        assertEquals(x.length, cdf.length);

        // CDF(0) = 0
        assertEquals(0.0, cdf[0], LOOSE_MID_TOL);

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

        // CDF(3) should be close to 1 (distribution has mean ~0.21)
        assertTrue(cdf[5] > 0.99, "CDF(3) should be close to 1 for mean ~0.21");
    }

    /**
     * CdfFromME: Verify CDF values at specific points.
     * Input: a=[0.2,0.3,0.5], A=[-1,0,0; 0,-3,2; 0,-2,-3]
     */
    @Test
    public void testCdfFromME() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2); alpha.set(0, 1, 0.3); alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        double[] x = {0.0, 0.5, 1.0, 2.0, 5.0};
        double[] cdf = CdfFromME.cdfFromME(alpha, A, x);

        assertNotNull(cdf);
        assertEquals(x.length, cdf.length);

        // CDF(0) = 0
        assertEquals(0.0, cdf[0], LOOSE_MID_TOL);

        // CDF should be monotonically non-decreasing
        for (int i = 1; i < cdf.length; i++) {
            assertTrue(cdf[i] >= cdf[i - 1] - LOOSE_MID_TOL,
                    "CDF should be non-decreasing at x=" + x[i]);
        }

        // CDF(5) should be close to 1 (distribution has mean ~0.35)
        assertTrue(cdf[4] > 0.99, "CDF(5) should be close to 1 for mean ~0.35");
    }

    /**
     * PdfFromPH: Verify PDF is non-negative and integrates to ~1.
     * Input: a=[0.1,0.9,0], A=[-6.2,2,0; 2,-9,1; 1,0,-3]
     */
    @Test
    public void testPdfFromPH() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.1); alpha.set(0, 1, 0.9); alpha.set(0, 2, 0.0);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -6.2); A.set(0, 1, 2.0); A.set(0, 2, 0.0);
        A.set(1, 0, 2.0);  A.set(1, 1, -9.0); A.set(1, 2, 1.0);
        A.set(2, 0, 1.0);  A.set(2, 1, 0.0); A.set(2, 2, -3.0);

        // Dense grid from 0 to 3 with step 0.01
        int nPoints = 301;
        double[] x = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            x[i] = i * 0.01;
        }

        double[] pdf = PdfFromME.pdfFromPH(alpha, A, x);

        assertNotNull(pdf);
        assertEquals(nPoints, pdf.length);

        // PDF should be non-negative
        for (int i = 0; i < nPoints; i++) {
            assertTrue(pdf[i] >= -LOOSE_MID_TOL, "PDF should be >= 0 at x=" + x[i]);
        }

        // Trapezoidal integration should approximate 1
        double integral = 0.0;
        for (int i = 1; i < nPoints; i++) {
            integral += (pdf[i - 1] + pdf[i]) * 0.5 * 0.01;
        }
        assertEquals(1.0, integral, COARSE_TOL, "PDF integral over [0,3] should be close to 1");
    }

    /**
     * PdfFromME: Verify PDF is non-negative and integrates to ~1.
     * Input: a=[0.2,0.3,0.5], A=[-1,0,0; 0,-3,2; 0,-2,-3]
     */
    @Test
    public void testPdfFromME() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2); alpha.set(0, 1, 0.3); alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        // Dense grid from 0 to 5 with step 0.01
        int nPoints = 501;
        double[] x = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            x[i] = i * 0.01;
        }

        double[] pdf = PdfFromME.pdfFromME(alpha, A, x);

        assertNotNull(pdf);
        assertEquals(nPoints, pdf.length);

        // PDF should be non-negative (ME can have slightly negative values, but this is a valid ME)
        for (int i = 0; i < nPoints; i++) {
            assertTrue(pdf[i] >= -LOOSE_MID_TOL, "PDF should be >= 0 at x=" + x[i]);
        }

        // Trapezoidal integration should approximate 1
        double integral = 0.0;
        for (int i = 1; i < nPoints; i++) {
            integral += (pdf[i - 1] + pdf[i]) * 0.5 * 0.01;
        }
        assertEquals(1.0, integral, COARSE_TOL, "PDF integral over [0,5] should be close to 1");
    }

    // ============ Phase 3: MEOrderFromMoments ============

    /**
     * MEOrderFromMoments: PH(3) with 5 moments should return order 3.
     * Input: moms from PH with a=[0.1,0.9,0], A=[-6.2,2,0; 2,-9,1; 1,0,-3]
     */
    @Test
    public void testMEOrderFromMoments() {
        double[] moms = {0.20939, 0.10449, 0.089092, 0.11027, 0.17953};
        int order = MEOrderFromMoments.meOrderFromMoments(moms, FINE_TOL);
        assertEquals(3, order, "Order of 3x3 PH should be 3");
    }

    /**
     * MEOrderFromMoments: Erlang-2 (order 2) from 3 moments.
     */
    @Test
    public void testMEOrderFromMoments_Erlang2() {
        // Erlang-2 with rate 1: m1=2, m2=6, m3=24
        double[] moms = {2.0, 6.0, 24.0};
        int order = MEOrderFromMoments.meOrderFromMoments(moms, FINE_TOL);
        assertEquals(2, order, "Order of Erlang-2 should be 2");
    }

    // ============ Phase 3: MEOrder ============

    /**
     * MEOrder: moment-based order for a known 3x3 ME.
     */
    @Test
    public void testMEOrder_Moment() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2); alpha.set(0, 1, 0.3); alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        int order = MEOrder.meOrder(alpha, A, "moment", FINE_TOL);
        assertEquals(3, order, "Moment order of full-rank 3x3 ME should be 3");
    }

    /**
     * MEOrder: observability order.
     */
    @Test
    public void testMEOrder_Obs() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2); alpha.set(0, 1, 0.3); alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        int order = MEOrder.meOrder(alpha, A, "obs", FINE_TOL);
        assertTrue(order >= 1 && order <= 3, "Obs order should be between 1 and 3");
    }

    /**
     * MEOrder: controllability order.
     */
    @Test
    public void testMEOrder_Cont() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2); alpha.set(0, 1, 0.3); alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        int order = MEOrder.meOrder(alpha, A, "cont", FINE_TOL);
        assertTrue(order >= 1 && order <= 3, "Cont order should be between 1 and 3");
    }

    // ============ Phase 3: AcyclicPHFromME ============

    /**
     * AcyclicPHFromME: Transform PH to acyclic form and verify moments match.
     * Input: a=[0.1,0.9,0], A=[-6.2,2,0; 2,-9,1; 1,0,-3] (all real eigenvalues)
     */
    @Test
    public void testAcyclicPHFromME() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.1); alpha.set(0, 1, 0.9); alpha.set(0, 2, 0.0);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -6.2); A.set(0, 1, 2.0); A.set(0, 2, 0.0);
        A.set(1, 0, 2.0);  A.set(1, 1, -9.0); A.set(1, 2, 1.0);
        A.set(2, 0, 1.0);  A.set(2, 1, 0.0); A.set(2, 2, -3.0);

        PHRepresentation result = AcyclicPHFromME.acyclicPHFromME(alpha, A, 100, 1e-14);
        assertNotNull(result);

        Matrix beta = result.getAlpha();
        Matrix B = result.getA();

        // Verify result is a valid PH representation
        assertTrue(CheckPHRepresentation.checkPHRepresentation(beta, B, 1e-12),
                "Acyclic result should be a valid PH");

        // Verify moments match original (first 5)
        double[] origMoms = MomentsFromME.momentsFromME(alpha, A, 5);
        double[] newMoms = MomentsFromME.momentsFromME(beta, B, 5);
        for (int i = 0; i < 5; i++) {
            assertEquals(origMoms[i], newMoms[i], Math.abs(origMoms[i]) * 1e-3,
                    "Moment " + (i + 1) + " should match");
        }
    }

    // ============ Phase 3: MonocyclicPHFromME ============

    /**
     * MonocyclicPHFromME: Transform ME to monocyclic PH and verify moments.
     * Input: a=[0.2,0.3,0.5], A=[-1,0,0; 0,-3,2; 0,-2,-3] (has complex eigenvalues)
     */
    @Test
    public void testMonocyclicPHFromME() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2); alpha.set(0, 1, 0.3); alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        PHRepresentation result = MonocyclicPHFromME.monocyclicPHFromME(alpha, A, 100, 1e-14);
        assertNotNull(result);

        Matrix beta = result.getAlpha();
        Matrix B = result.getA();

        // Verify result is a valid PH representation
        assertTrue(CheckPHRepresentation.checkPHRepresentation(beta, B, 1e-12),
                "Monocyclic result should be a valid PH");

        // Verify moments match original (first 5)
        double[] origMoms = MomentsFromME.momentsFromME(alpha, A, 5);
        double[] newMoms = MomentsFromME.momentsFromME(beta, B, 5);
        for (int i = 0; i < 5; i++) {
            assertEquals(origMoms[i], newMoms[i], Math.abs(origMoms[i]) * 1e-3,
                    "Moment " + (i + 1) + " should match");
        }
    }

    // ============ Phase 3: CheckMEPositiveDensity ============

    /**
     * CheckMEPositiveDensity: Valid PH should have positive density.
     */
    @Test
    public void testCheckMEPositiveDensity_Valid() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.1); alpha.set(0, 1, 0.9); alpha.set(0, 2, 0.0);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -6.2); A.set(0, 1, 2.0); A.set(0, 2, 0.0);
        A.set(1, 0, 2.0);  A.set(1, 1, -9.0); A.set(1, 2, 1.0);
        A.set(2, 0, 1.0);  A.set(2, 1, 0.0); A.set(2, 2, -3.0);

        boolean result = CheckMEPositiveDensity.checkMEPositiveDensity(alpha, A, 100, 1e-14);
        assertTrue(result, "Valid PH should have positive density");
    }

    /**
     * CheckMEPositiveDensity: Valid ME should have positive density.
     */
    @Test
    public void testCheckMEPositiveDensity_ValidME() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2); alpha.set(0, 1, 0.3); alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        boolean result = CheckMEPositiveDensity.checkMEPositiveDensity(alpha, A, 100, 1e-14);
        assertTrue(result, "Valid ME should have positive density");
    }

    // ============ Phase 3: PHFromME ============

    /**
     * PHFromME: Transform ME to PH of the same size and verify moments.
     * Input: a=[0.2,0.3,0.5], A=[-1,0,0; 0,-3,2; 0,-2,-3]
     */
    @Test
    public void testPHFromME() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2); alpha.set(0, 1, 0.3); alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        PHRepresentation result = PHFromME.phFromME(alpha, A, 1e-7);
        assertNotNull(result);

        Matrix beta = result.getAlpha();
        Matrix B = result.getA();

        // Result should have the same size
        assertEquals(3, beta.getNumCols(), "PHFromME result should have same size");
        assertEquals(3, B.getNumRows(), "PHFromME result should have same size");

        // Verify moments match original (first 5)
        double[] origMoms = MomentsFromME.momentsFromME(alpha, A, 5);
        double[] newMoms = MomentsFromME.momentsFromME(beta, B, 5);
        for (int i = 0; i < 5; i++) {
            assertEquals(origMoms[i], newMoms[i], Math.abs(origMoms[i]) * 0.01,
                    "Moment " + (i + 1) + " should match");
        }
    }

    /**
     * PHFromME: Transform a PH (already Markovian) - should return similar representation.
     */
    @Test
    public void testPHFromME_AlreadyPH() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.1); alpha.set(0, 1, 0.9); alpha.set(0, 2, 0.0);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -6.2); A.set(0, 1, 2.0); A.set(0, 2, 0.0);
        A.set(1, 0, 2.0);  A.set(1, 1, -9.0); A.set(1, 2, 1.0);
        A.set(2, 0, 1.0);  A.set(2, 1, 0.0); A.set(2, 2, -3.0);

        PHRepresentation result = PHFromME.phFromME(alpha, A, 1e-7);
        assertNotNull(result);

        // Verify moments match
        double[] origMoms = MomentsFromME.momentsFromME(alpha, A, 5);
        double[] newMoms = MomentsFromME.momentsFromME(result.getAlpha(), result.getA(), 5);
        for (int i = 0; i < 5; i++) {
            assertEquals(origMoms[i], newMoms[i], Math.abs(origMoms[i]) * 0.01,
                    "Moment " + (i + 1) + " should match");
        }
    }

    // ============ Phase 3: IntervalPdfFromPH ============

    /**
     * IntervalPdfFromPH: Verify interval PDF values are non-negative and consistent with CDF.
     */
    @Test
    public void testIntervalPdfFromPH() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.1); alpha.set(0, 1, 0.9); alpha.set(0, 2, 0.0);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -6.2); A.set(0, 1, 2.0); A.set(0, 2, 0.0);
        A.set(1, 0, 2.0);  A.set(1, 1, -9.0); A.set(1, 2, 1.0);
        A.set(2, 0, 1.0);  A.set(2, 1, 0.0); A.set(2, 2, -3.0);

        // Create interval boundaries from 0 to 3 with step 0.1
        int nBounds = 31;
        double[] intBounds = new double[nBounds];
        for (int i = 0; i < nBounds; i++) {
            intBounds[i] = i * 0.1;
        }

        jline.util.Pair<double[], double[]> result = IntervalPdfFromPH.intervalPdfFromPH(alpha, A, intBounds);
        double[] x = result.getFirst();
        double[] y = result.getSecond();

        assertNotNull(x);
        assertNotNull(y);
        assertEquals(nBounds - 1, x.length);
        assertEquals(nBounds - 1, y.length);

        // x should be the midpoints
        for (int i = 0; i < x.length; i++) {
            assertEquals((intBounds[i] + intBounds[i + 1]) / 2.0, x[i], FINE_TOL);
        }

        // y values should be non-negative (interval density from a valid PH)
        for (int i = 0; i < y.length; i++) {
            assertTrue(y[i] >= -LOOSE_MID_TOL, "Interval PDF should be >= 0 at x=" + x[i]);
        }

        // Sum of y * interval_width should approximate 1 (total probability over [0,3])
        double totalProb = 0.0;
        for (int i = 0; i < y.length; i++) {
            totalProb += y[i] * (intBounds[i + 1] - intBounds[i]);
        }
        assertEquals(1.0, totalProb, COARSE_TOL, "Total interval probability should be close to 1");
    }

    // ============ Phase 3: MinimalRepFromME ============

    /**
     * MinimalRepFromME (moment method): Verify moments are preserved.
     */
    @Test
    public void testMinimalRepFromME_Moment() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2); alpha.set(0, 1, 0.3); alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0); A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        MERepresentation result = MinimalRepFromME.minimalRepFromME(alpha, A, "moment", FINE_TOL);
        assertNotNull(result);

        // Verify moments match original (first 5)
        double[] origMoms = MomentsFromME.momentsFromME(alpha, A, 5);
        double[] newMoms = MomentsFromME.momentsFromME(result.getAlpha(), result.getA(), 5);
        for (int i = 0; i < 5; i++) {
            assertEquals(origMoms[i], newMoms[i], Math.abs(origMoms[i]) * 1e-3,
                    "Moment " + (i + 1) + " should match");
        }
    }
}
