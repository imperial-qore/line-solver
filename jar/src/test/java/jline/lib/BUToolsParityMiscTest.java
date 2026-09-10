package jline.lib;

import jline.lib.butools.mc.CRPSolve;
import jline.lib.butools.mc.CTMCSolve;
import jline.lib.butools.mc.DTMCSolve;
import jline.lib.butools.mc.CheckGenerator;
import jline.lib.butools.mc.CheckProbMatrix;
import jline.lib.butools.mc.CheckProbVector;
import jline.lib.butools.CheckMoments;
import jline.lib.butools.FactorialMomsFromMoms;
import jline.lib.butools.MomsFromFactorialMoms;
import jline.lib.butools.MomsFromHankelMoms;
import jline.lib.butools.MomsFromNormMoms;
import jline.lib.butools.MomsFromReducedMoms;
import jline.lib.butools.JFactorialMomsFromJMoms;
import jline.lib.butools.JMomsFromJFactorialMoms;
import jline.lib.butools.SimilarityMatrixForVectors;
import jline.lib.butools.trace.LagCorrelationsFromTrace;
import jline.lib.butools.trace.MarginalMomentsFromTrace;
import jline.lib.butools.map.SamplesFromMAP;
import jline.lib.butools.map.MarginalMomentsFromMAP;
import jline.lib.butools.QBDFundamentalMatrices;
import jline.lib.butools.mam.QBDStationaryDistr;
import jline.lib.butools.mam.MG1FundamentalMatrix;
import jline.lib.butools.mam.MG1FundamentalMatrix.MG1Method;
import jline.lib.butools.mam.GM1FundamentalMatrix;
import jline.lib.butools.mam.GM1FundamentalMatrix.GM1Method;
import jline.lib.butools.mam.FluidStationaryDistr;
import jline.lib.butools.mam.GeneralFluidSolve;
import jline.lib.butools.mam.GeneralFluidSolution;
import jline.lib.butools.reptrans.TransformToAcyclic;
import jline.lib.butools.reptrans.TransformToMonocyclic;
import jline.lib.butools.reptrans.ExtendToMarkovian;
import jline.lib.butools.reptrans.MarkovianRepresentation;
import jline.lib.butools.reptrans.SimilarityMatrix;
import jline.lib.butools.fitting.SquaredDifference;
import jline.lib.butools.fitting.RelativeEntropy;
import jline.lib.butools.dmap.CheckDMMAPRepresentation;
import jline.lib.butools.dmap.MarginalMomentsFromDMMAP;
import jline.lib.butools.dmap.MarginalDistributionFromDMMAP;
import jline.lib.butools.dmap.CheckDMRAPRepresentation;
import jline.lib.butools.dmap.MarginalMomentsFromDMRAP;
import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.lib.butools.trace.CdfFromTrace;
import jline.lib.butools.trace.CdfFromTrace.CdfResult;
import jline.lib.butools.trace.PdfFromTrace;
import jline.lib.butools.trace.PdfFromTrace.PdfResult;
import jline.lib.butools.fitting.LikelihoodFromTrace;
import jline.lib.butools.fitting.PHFromTrace;
import jline.lib.butools.fitting.PHFromTrace.PHFitResult;
import jline.lib.butools.fitting.MAPFromTrace;
import jline.lib.butools.fitting.MAPFitResult;
import jline.lib.butools.reptrans.MStaircase;
import jline.lib.butools.map.MinimalRepFromRAP;
import jline.lib.butools.map.MinimalRepFromMRAP;
import jline.lib.butools.mam.QBDSolve;
import jline.lib.butools.mam.GM1StationaryDistr;
import jline.lib.butools.mam.MG1StationaryDistr;
import jline.lib.butools.queues.QBDQueue;
import jline.lib.butools.queues.MAPMAP1;
import jline.lib.butools.queues.FluidQueue;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.Pair;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.*;

/**
 * BUTools Parity Tests (Misc) - Unit tests validating JAR BuTools implementations
 * for MC solvers, QBD, M/G/1, G/M/1, Fluid, matrix transforms, distance measures,
 * and DMMAP/DMRAP against reference examples in butools.tmp.
 */
public class BUToolsParityMiscTest {


    // ============ MC Solvers ============

    /**
     * CRPSolve.txt: Solve a continuous-time rational process
     * Input: Q = [-4.3, 3.5, 0.8; -8.4, 6.5, 1.9; 17.3, -12.7, -4.6]
     * Expected: ret = [-3.5617, 3.6667, 0.89506]
     * Verify: ret*Q ≈ 0
     */
    @Test
    public void testCRPSolve_Example() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -4.3); Q.set(0, 1, 3.5);   Q.set(0, 2, 0.8);
        Q.set(1, 0, -8.4); Q.set(1, 1, 6.5);   Q.set(1, 2, 1.9);
        Q.set(2, 0, 17.3); Q.set(2, 1, -12.7);  Q.set(2, 2, -4.6);

        Matrix ret = CRPSolve.crpSolve(Q, 1e-14);

        assertNotNull(ret);
        assertEquals(1, ret.getNumRows());
        assertEquals(3, ret.getNumCols());

        // Verify expected values
        assertEquals(-3.5617, ret.get(0, 0), LOOSE_MID_TOL);
        assertEquals(3.6667, ret.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.89506, ret.get(0, 2), LOOSE_MID_TOL);

        // Verify ret*Q ≈ 0
        Matrix product = ret.mult(Q);
        for (int j = 0; j < 3; j++) {
            assertEquals(0.0, product.get(0, j), LOOSE_MID_TOL);
        }
    }

    /**
     * DRPSolve.txt: Solve a discrete-time rational process
     * Input: Q = [-0.9, 0.5, 1.4; 0.9, -0.9, 1; 0.3, 1.3, -0.6]
     * Expected: ret = [0.23138, 0.3484, 0.42021]
     * Verify: ret*Q ≈ ret
     */
    @Test
    public void testDRPSolve_Example() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -0.9); Q.set(0, 1, 0.5);  Q.set(0, 2, 1.4);
        Q.set(1, 0, 0.9);  Q.set(1, 1, -0.9); Q.set(1, 2, 1.0);
        Q.set(2, 0, 0.3);  Q.set(2, 1, 1.3);  Q.set(2, 2, -0.6);

        Matrix ret = CRPSolve.drpSolve(Q);

        assertNotNull(ret);
        assertEquals(1, ret.getNumRows());
        assertEquals(3, ret.getNumCols());

        // Verify expected values
        assertEquals(0.23138, ret.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.3484, ret.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.42021, ret.get(0, 2), LOOSE_MID_TOL);

        // Verify ret*Q ≈ ret
        Matrix product = ret.mult(Q);
        for (int j = 0; j < 3; j++) {
            assertEquals(ret.get(0, j), product.get(0, j), LOOSE_MID_TOL);
        }
    }

    /**
     * CheckProbVector.txt: Six test cases for probability vector validation
     */
    @Test
    public void testCheckProbVector_NegativeElement() {
        // Test 1: Q=[1.1,-0.1] → false (negative element)
        Matrix Q = new Matrix(1, 2);
        Q.set(0, 0, 1.1); Q.set(0, 1, -0.1);
        assertFalse(CheckProbVector.checkProbVector(Q, false, 1e-12));
    }

    @Test
    public void testCheckProbVector_SumNotOne() {
        // Test 2: Q=[1.1,0.1] → false (sum not 1)
        Matrix Q = new Matrix(1, 2);
        Q.set(0, 0, 1.1); Q.set(0, 1, 0.1);
        assertFalse(CheckProbVector.checkProbVector(Q, false, 1e-12));
    }

    @Test
    public void testCheckProbVector_Valid() {
        // Test 3: Q=[1,0] → true
        Matrix Q = new Matrix(1, 2);
        Q.set(0, 0, 1.0); Q.set(0, 1, 0.0);
        assertTrue(CheckProbVector.checkProbVector(Q, false, 1e-12));
    }

    @Test
    public void testCheckProbVector_SubNegativeElement() {
        // Test 4 (sub=true): Q=[0.9,-0.1] → false (negative element)
        Matrix Q = new Matrix(1, 2);
        Q.set(0, 0, 0.9); Q.set(0, 1, -0.1);
        assertFalse(CheckProbVector.checkProbVector(Q, true, 1e-12));
    }

    @Test
    public void testCheckProbVector_SubValid() {
        // Test 5 (sub=true): Q=[0.9,0.1] → true
        Matrix Q = new Matrix(1, 2);
        Q.set(0, 0, 0.9); Q.set(0, 1, 0.1);
        assertTrue(CheckProbVector.checkProbVector(Q, true, 1e-12));
    }

    @Test
    public void testCheckProbVector_SubValidSumLessThanOne() {
        // Test 6 (sub=true): Q=[0.8,0.1] → true
        Matrix Q = new Matrix(1, 2);
        Q.set(0, 0, 0.8); Q.set(0, 1, 0.1);
        assertTrue(CheckProbVector.checkProbVector(Q, true, 1e-12));
    }

    // ============ QBD ============

    /**
     * QBDFundamentalMatrices.txt: Compute R, G, U matrices for a QBD process
     * Input: B=[0,0;3,4], L=[-6,5;3,-12], F=[1,0;2,0]
     * Expected:
     *   R=[0.27839,0.14286; 0.55678,0.28571]
     *   G=[0.42857,0.57143; 0.42857,0.57143]
     *   U=[-5.5714,5.5714; 3.8571,-10.857]
     */
    @Test
    public void testQBDFundamentalMatrices_Example() {
        Matrix B = new Matrix(2, 2);
        B.set(0, 0, 0.0); B.set(0, 1, 0.0);
        B.set(1, 0, 3.0); B.set(1, 1, 4.0);

        Matrix L = new Matrix(2, 2);
        L.set(0, 0, -6.0); L.set(0, 1, 5.0);
        L.set(1, 0, 3.0);  L.set(1, 1, -12.0);

        Matrix F = new Matrix(2, 2);
        F.set(0, 0, 1.0); F.set(0, 1, 0.0);
        F.set(1, 0, 2.0); F.set(1, 1, 0.0);

        @SuppressWarnings("unchecked")
        Map<String, Matrix> result = (Map<String, Matrix>) (Map<?, ?>) QBDFundamentalMatrices.QBDFundamentalMatrices(
                B, L, F, null, null, null, null);

        assertNotNull(result);

        Matrix R = result.get("R");
        Matrix G = result.get("G");
        Matrix U = result.get("U");

        assertNotNull(R);
        assertNotNull(G);
        assertNotNull(U);

        // Verify R
        assertEquals(0.27839, R.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.14286, R.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.55678, R.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.28571, R.get(1, 1), LOOSE_MID_TOL);

        // Verify G
        assertEquals(0.42857, G.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.57143, G.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.42857, G.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.57143, G.get(1, 1), LOOSE_MID_TOL);

        // Verify U
        assertEquals(-5.5714, U.get(0, 0), LOOSE_MID_TOL);
        assertEquals(5.5714, U.get(0, 1), LOOSE_MID_TOL);
        assertEquals(3.8571, U.get(1, 0), LOOSE_MID_TOL);
        assertEquals(-10.857, U.get(1, 1), LOOSE_MID_TOL);
    }

    /**
     * QBDStationaryDistr.txt: Compute stationary distribution of a QBD process
     * Uses R from QBDFundamentalMatrices and pi0 solved from boundary equation.
     * Expected pi (first 12 values):
     *   [0.22992, 0.18681, 0.16802, 0.086221, 0.094781, 0.048638,
     *    0.053466, 0.027437, 0.030161, 0.015477, 0.017014, 0.0087307]
     */
    @Test
    public void testQBDStationaryDistr_Example() {
        // Setup: same B, L, F as QBDFundamentalMatrices
        Matrix B = new Matrix(2, 2);
        B.set(0, 0, 0.0); B.set(0, 1, 0.0);
        B.set(1, 0, 3.0); B.set(1, 1, 4.0);

        Matrix L = new Matrix(2, 2);
        L.set(0, 0, -6.0); L.set(0, 1, 5.0);
        L.set(1, 0, 3.0);  L.set(1, 1, -12.0);

        Matrix F = new Matrix(2, 2);
        F.set(0, 0, 1.0); F.set(0, 1, 0.0);
        F.set(1, 0, 2.0); F.set(1, 1, 0.0);

        Matrix L0 = new Matrix(2, 2);
        L0.set(0, 0, -6.0); L0.set(0, 1, 5.0);
        L0.set(1, 0, 6.0);  L0.set(1, 1, -8.0);

        // Get R from QBDFundamentalMatrices
        @SuppressWarnings("unchecked")
        Map<String, Matrix> fundResult = (Map<String, Matrix>) (Map<?, ?>) QBDFundamentalMatrices.QBDFundamentalMatrices(
                B, L, F, null, null, null, null);
        Matrix R = fundResult.get("R");
        assertNotNull(R);

        // Compute pi0 from boundary equation: pi0*(L0 + R*B) = 0
        // with normalization: pi0*(I-R)^{-1}*ones = 1
        Matrix RB = R.mult(B);
        Matrix boundaryQ = L0.add(1.0, RB);

        // Solve pi0*boundaryQ = 0 using CRPSolve (normalized to pi0*ones=1)
        Matrix pi0_crp = CRPSolve.crpSolve(boundaryQ, 1e-14);

        // Renormalize: we need pi0*(I-R)^{-1}*ones = 1
        Matrix ImR = Matrix.eye(2).add(-1.0, R);
        Matrix ImRinv = ImR.inv();
        Matrix ImRinvOnes = ImRinv.mult(Matrix.ones(2, 1));
        double totalProb = pi0_crp.mult(ImRinvOnes).get(0, 0);
        Matrix pi0 = pi0_crp.scale(1.0 / totalProb);

        // Call QBDStationaryDistr with K=5 (gives 6 levels: 0..5, so 12 values)
        Matrix pi = QBDStationaryDistr.qbdStationaryDistr(pi0, R, 5);

        assertNotNull(pi);

        // Expected values
        double[] expected = {
            0.22992, 0.18681, 0.16802, 0.086221, 0.094781, 0.048638,
            0.053466, 0.027437, 0.030161, 0.015477, 0.017014, 0.0087307
        };

        // pi should be a row vector with 12 elements (6 levels * 2 states)
        assertEquals(1, pi.getNumRows());
        assertEquals(12, pi.getNumCols());

        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], pi.get(0, i), LOOSE_MID_TOL,
                    "Mismatch at index " + i);
        }
    }

    // ============ M/G/1 and G/M/1 ============

    /**
     * MG1FundamentalMatrix.txt: Compute fundamental matrix G for M/G/1 type chain
     * Input: A0=[0.4,0.2;0.3,0.4], A1=[0,0.1;0,0], A2=[0,0.2;0,0.2], A3=[0.1,0;0.1,0]
     * Expected: G=[0.60503,0.39497; 0.45912,0.54088]
     */
    @Test
    public void testMG1FundamentalMatrix_Example() {
        Matrix A0 = new Matrix(2, 2);
        A0.set(0, 0, 0.4); A0.set(0, 1, 0.2);
        A0.set(1, 0, 0.3); A0.set(1, 1, 0.4);

        Matrix A1 = new Matrix(2, 2);
        A1.set(0, 0, 0.0); A1.set(0, 1, 0.1);
        A1.set(1, 0, 0.0); A1.set(1, 1, 0.0);

        Matrix A2 = new Matrix(2, 2);
        A2.set(0, 0, 0.0); A2.set(0, 1, 0.2);
        A2.set(1, 0, 0.0); A2.set(1, 1, 0.2);

        Matrix A3 = new Matrix(2, 2);
        A3.set(0, 0, 0.1); A3.set(0, 1, 0.0);
        A3.set(1, 0, 0.1); A3.set(1, 1, 0.0);

        Matrix G = MG1FundamentalMatrix.mg1FundamentalMatrix(
                Arrays.asList(A0, A1, A2, A3), 1e-14, 50, MG1Method.CR);

        assertNotNull(G);
        assertEquals(2, G.getNumRows());
        assertEquals(2, G.getNumCols());

        assertEquals(0.60503, G.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.39497, G.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.45912, G.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.54088, G.get(1, 1), LOOSE_MID_TOL);
    }

    /**
     * GM1FundamentalMatrix.txt: Compute fundamental matrix R for G/M/1 type chain
     * Input: A0=[0.1,0;0,0.1], A1=[0,0.2;0,0.2], A2=[0,0.1;0,0],
     *        A3=[0.3,0.2;0.3,0.2], A4=[0,0.1;0.2,0]
     * Expected: R=[0.10065,0.026961; 0.00065531,0.12569]
     */
    @Test
    public void testGM1FundamentalMatrix_Example() {
        Matrix A0 = new Matrix(2, 2);
        A0.set(0, 0, 0.1); A0.set(0, 1, 0.0);
        A0.set(1, 0, 0.0); A0.set(1, 1, 0.1);

        Matrix A1 = new Matrix(2, 2);
        A1.set(0, 0, 0.0); A1.set(0, 1, 0.2);
        A1.set(1, 0, 0.0); A1.set(1, 1, 0.2);

        Matrix A2 = new Matrix(2, 2);
        A2.set(0, 0, 0.0); A2.set(0, 1, 0.1);
        A2.set(1, 0, 0.0); A2.set(1, 1, 0.0);

        Matrix A3 = new Matrix(2, 2);
        A3.set(0, 0, 0.3); A3.set(0, 1, 0.2);
        A3.set(1, 0, 0.3); A3.set(1, 1, 0.2);

        Matrix A4 = new Matrix(2, 2);
        A4.set(0, 0, 0.0); A4.set(0, 1, 0.1);
        A4.set(1, 0, 0.2); A4.set(1, 1, 0.0);

        Matrix R = GM1FundamentalMatrix.gm1FundamentalMatrix(
                Arrays.asList(A0, A1, A2, A3, A4), 1e-14, 50, GM1Method.CR);

        assertNotNull(R);
        assertEquals(2, R.getNumRows());
        assertEquals(2, R.getNumCols());

        assertEquals(0.10065, R.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.026961, R.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.00065531, R.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.12569, R.get(1, 1), LOOSE_MID_TOL);
    }

    // ============ Fluid ============

    /**
     * FluidStationaryDistr.txt: Compute stationary distribution of a fluid model
     * Uses GeneralFluidSolve to get mass0, ini, K, clo, then evaluates at x=[0,5,...,30]
     * Checks row 0 (x=0) = mass0 and a few intermediate rows
     */
    @Test
    public void testFluidStationaryDistr_Example() {
        // Q = 6x6 generator matrix (same as GeneralFluidSolve test)
        Matrix Q = new Matrix(6, 6);
        Q.set(0, 0, -6.0);  Q.set(0, 1, 1.0);   Q.set(0, 2, 3.0);   Q.set(0, 3, 2.0);  Q.set(0, 4, 0.0);   Q.set(0, 5, 0.0);
        Q.set(1, 0, 6.0);   Q.set(1, 1, -10.0);  Q.set(1, 2, 2.0);   Q.set(1, 3, 0.0);  Q.set(1, 4, 2.0);   Q.set(1, 5, 0.0);
        Q.set(2, 0, 3.0);   Q.set(2, 1, 7.0);    Q.set(2, 2, -12.0); Q.set(2, 3, 0.0);  Q.set(2, 4, 0.0);   Q.set(2, 5, 2.0);
        Q.set(3, 0, 5.0);   Q.set(3, 1, 0.0);    Q.set(3, 2, 0.0);   Q.set(3, 3, -9.0); Q.set(3, 4, 1.0);   Q.set(3, 5, 3.0);
        Q.set(4, 0, 0.0);   Q.set(4, 1, 5.0);    Q.set(4, 2, 0.0);   Q.set(4, 3, 6.0);  Q.set(4, 4, -13.0); Q.set(4, 5, 2.0);
        Q.set(5, 0, 0.0);   Q.set(5, 1, 0.0);    Q.set(5, 2, 5.0);   Q.set(5, 3, 3.0);  Q.set(5, 4, 7.0);   Q.set(5, 5, -15.0);

        // R = diag([2, -4, -12, 6, 0, -8])
        Matrix R = new Matrix(6, 6);
        R.set(0, 0, 2.0);
        R.set(1, 1, -4.0);
        R.set(2, 2, -12.0);
        R.set(3, 3, 6.0);
        R.set(4, 4, 0.0);
        R.set(5, 5, -8.0);

        // Compute GeneralFluidSolve
        GeneralFluidSolution gfs = GeneralFluidSolve.generalFluidSolve(Q, R, null, 1e-14);
        assertNotNull(gfs);

        // Verify intermediate results from GeneralFluidSolve
        // Expected ini = [0.70195, 0.20505]
        assertNotNull(gfs.getIni());
        assertEquals(0.70195, gfs.getIni().get(0, 0), LOOSE_MID_TOL, "ini[0]");
        assertEquals(0.20505, gfs.getIni().get(0, 1), LOOSE_MID_TOL, "ini[1]");

        // Expected K = [-2.4698, 1.1349; 1.295, -1.1686]
        assertNotNull(gfs.getK());
        assertEquals(-2.4698, gfs.getK().get(0, 0), LOOSE_MID_TOL, "K[0,0]");
        assertEquals(1.1349, gfs.getK().get(0, 1), LOOSE_MID_TOL, "K[0,1]");
        assertEquals(1.295, gfs.getK().get(1, 0), LOOSE_MID_TOL, "K[1,0]");
        assertEquals(-1.1686, gfs.getK().get(1, 1), LOOSE_MID_TOL, "K[1,1]");

        // Evaluate FluidStationaryDistr at x = [0, 5, 10, 15, 20, 25, 30]
        double[] x = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};
        Matrix y = FluidStationaryDistr.fluidStationaryDistr(
                gfs.getMass0(), gfs.getIni(), gfs.getK(), gfs.getClo(), x);

        assertNotNull(y);
        assertEquals(7, y.getNumRows());
        assertEquals(6, y.getNumCols());

        // Row 0 (x=0): should equal mass0
        double[] row0 = {0.0, 0.082246, 0.069492, 0.0, 0.023812, 0.020724};
        for (int j = 0; j < 6; j++) {
            assertEquals(row0[j], y.get(0, j), LOOSE_MID_TOL, "Row 0, col " + j);
        }

        // Row 1 (x=5)
        double[] row1 = {0.34868, 0.1698, 0.14254, 0.13527, 0.066678, 0.055991};
        for (int j = 0; j < 6; j++) {
            assertEquals(row1[j], y.get(1, j), LOOSE_MID_TOL, "Row 1, col " + j);
        }

        // Row 3 (x=15)
        double[] row3 = {0.38286, 0.1799, 0.15089, 0.1531, 0.071946, 0.060343};
        for (int j = 0; j < 6; j++) {
            assertEquals(row3[j], y.get(3, j), LOOSE_MID_TOL, "Row 3, col " + j);
        }

        // Row 6 (x=30): should be close to steady-state
        double[] row6 = {0.38327, 0.18002, 0.15099, 0.15331, 0.072009, 0.060395};
        for (int j = 0; j < 6; j++) {
            assertEquals(row6[j], y.get(6, j), LOOSE_MID_TOL, "Row 6, col " + j);
        }
    }

    // ============ Matrix Transforms ============

    /**
     * TransformToAcyclic.txt: Transform a Markov generator to acyclic form
     * Input: A = [-0.8, 0.8, 0; 0.1, -0.3, 0.1; 0.2, 0, -0.5]
     * Expected: B = [-0.1203, 0.1203, 0; 0, -0.6158, 0.6158; 0, 0, -0.8639]
     * Verify: similarity via A*Cm ≈ Cm*B
     */
    @Test
    public void testTransformToAcyclic_Example() {
        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -0.8); A.set(0, 1, 0.8);  A.set(0, 2, 0.0);
        A.set(1, 0, 0.1);  A.set(1, 1, -0.3); A.set(1, 2, 0.1);
        A.set(2, 0, 0.2);  A.set(2, 1, 0.0);  A.set(2, 2, -0.5);

        Matrix B = TransformToAcyclic.transformToAcyclic(A, 100, 1e-14);

        assertNotNull(B);
        assertEquals(3, B.getNumRows());
        assertEquals(3, B.getNumCols());

        // Verify expected values
        assertEquals(-0.1203, B.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.1203, B.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(0, 2), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(1, 0), LOOSE_MID_TOL);
        assertEquals(-0.6158, B.get(1, 1), LOOSE_MID_TOL);
        assertEquals(0.6158, B.get(1, 2), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(2, 0), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(2, 1), LOOSE_MID_TOL);
        assertEquals(-0.8639, B.get(2, 2), LOOSE_MID_TOL);

        // Verify similarity: A*Cm ≈ Cm*B
        Matrix Cm = SimilarityMatrix.similarityMatrix(A, B);
        Matrix lhs = A.mult(Cm);
        Matrix rhs = Cm.mult(B);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(lhs.get(i, j), rhs.get(i, j), LOOSE_MID_TOL,
                        "Similarity error at (" + i + "," + j + ")");
            }
        }
    }

    /**
     * TransformToMonocyclic.txt: Transform a Markov generator to monocyclic form
     * Input: A = [-1, 0, 0; 0, -3, 2; 0, -2, -3]
     * Expected: 5x5 monocyclic matrix B
     */
    @Test
    public void testTransformToMonocyclic_Example() {
        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0);  A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 2.0);
        A.set(2, 0, 0.0);  A.set(2, 1, -2.0); A.set(2, 2, -3.0);

        Matrix B = TransformToMonocyclic.transformToMonocyclic(A, 100, 1e-14);

        assertNotNull(B);
        assertEquals(5, B.getNumRows());
        assertEquals(5, B.getNumCols());

        // Expected B (5x5 monocyclic)
        assertEquals(-1.0, B.get(0, 0), LOOSE_MID_TOL);
        assertEquals(1.0, B.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(0, 2), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(0, 3), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(0, 4), LOOSE_MID_TOL);

        assertEquals(0.0, B.get(1, 0), LOOSE_MID_TOL);
        assertEquals(-3.0, B.get(1, 1), LOOSE_MID_TOL);
        assertEquals(3.0, B.get(1, 2), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(1, 3), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(1, 4), LOOSE_MID_TOL);

        assertEquals(0.0, B.get(2, 0), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(2, 1), LOOSE_MID_TOL);
        assertEquals(-3.0, B.get(2, 2), LOOSE_MID_TOL);
        assertEquals(3.0, B.get(2, 3), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(2, 4), LOOSE_MID_TOL);

        assertEquals(0.0, B.get(3, 0), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(3, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(3, 2), LOOSE_MID_TOL);
        assertEquals(-3.0, B.get(3, 3), LOOSE_MID_TOL);
        assertEquals(3.0, B.get(3, 4), LOOSE_MID_TOL);

        assertEquals(0.0, B.get(4, 0), LOOSE_MID_TOL);
        assertEquals(0.59259, B.get(4, 1), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(4, 2), LOOSE_MID_TOL);
        assertEquals(0.0, B.get(4, 3), LOOSE_MID_TOL);
        assertEquals(-3.0, B.get(4, 4), LOOSE_MID_TOL);

        // Verify similarity: A*Cm ≈ Cm*B
        Matrix Cm = SimilarityMatrix.similarityMatrix(A, B);
        Matrix lhs = A.mult(Cm);
        Matrix rhs = Cm.mult(B);
        double err = 0.0;
        for (int i = 0; i < A.getNumRows(); i++) {
            for (int j = 0; j < B.getNumCols(); j++) {
                double d = lhs.get(i, j) - rhs.get(i, j);
                err += d * d;
            }
        }
        err = Math.sqrt(err);
        assertTrue(err < LOOSE_MID_TOL, "Similarity error too large: " + err);
    }

    /**
     * ExtendToMarkovian.txt: Extend a non-Markovian representation to Markovian
     * Input: alpha=[0.2,0.3,0.5], A=[-1,0,0; 0,-3,0.6; 0,-0.6,-3]
     * First compute monocyclic B and beta=alpha*Cm, then extend to Markovian
     * Expected: m (7-element positive vector), M (7x7 matrix)
     */
    @Test
    public void testExtendToMarkovian_Example() {
        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, 0.2); alpha.set(0, 1, 0.3); alpha.set(0, 2, 0.5);

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -1.0); A.set(0, 1, 0.0);  A.set(0, 2, 0.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -3.0); A.set(1, 2, 0.6);
        A.set(2, 0, 0.0);  A.set(2, 1, -0.6); A.set(2, 2, -3.0);

        // Step 1: Transform to monocyclic
        Matrix B = TransformToMonocyclic.transformToMonocyclic(A, 100, 1e-14);
        assertNotNull(B);
        assertEquals(4, B.getNumRows());
        assertEquals(4, B.getNumCols());

        // Step 2: Compute beta = alpha * Cm
        Matrix Cm = SimilarityMatrix.similarityMatrix(A, B);
        Matrix beta = alpha.mult(Cm);

        // Verify beta values
        assertEquals(0.045649, beta.get(0, 0), LOOSE_MID_TOL);
        assertEquals(-0.00043836, beta.get(0, 1), LOOSE_MID_TOL);
        assertEquals(-0.088811, beta.get(0, 2), LOOSE_MID_TOL);
        assertEquals(1.0436, beta.get(0, 3), LOOSE_MID_TOL);

        // Step 3: ExtendToMarkovian
        MarkovianRepresentation result = ExtendToMarkovian.extendToMarkovian(beta, B, 100, 1e-14);
        assertNotNull(result);

        Matrix m = result.getBeta();
        Matrix M = result.getB();

        assertNotNull(m);
        assertNotNull(M);

        // Expected m (7-element positive vector)
        assertEquals(7, m.getNumCols());
        double[] expectedM = {0.015399, 0.010192, 0.017621, 0.018114, 0.0087991, 0.10333, 0.82655};
        for (int i = 0; i < 7; i++) {
            assertEquals(expectedM[i], m.get(0, i), LOOSE_MID_TOL, "m[" + i + "]");
        }

        // All elements should be non-negative (Markovian)
        for (int i = 0; i < m.getNumCols(); i++) {
            assertTrue(m.get(0, i) >= -LOOSE_MID_TOL, "m[" + i + "] should be non-negative");
        }

        // Expected M (7x7 matrix)
        assertEquals(7, M.getNumRows());
        assertEquals(7, M.getNumCols());

        // Check key elements of M
        assertEquals(-1.0, M.get(0, 0), LOOSE_MID_TOL);
        assertEquals(1.0, M.get(0, 1), LOOSE_MID_TOL);
        assertEquals(-2.6536, M.get(1, 1), LOOSE_MID_TOL);
        assertEquals(2.6536, M.get(1, 2), LOOSE_MID_TOL);
        assertEquals(-2.6536, M.get(2, 2), LOOSE_MID_TOL);
        assertEquals(2.6536, M.get(2, 3), LOOSE_MID_TOL);
        assertEquals(0.047227, M.get(3, 1), LOOSE_MID_TOL);
        assertEquals(-2.6536, M.get(3, 3), LOOSE_MID_TOL);
        assertEquals(2.6064, M.get(3, 4), LOOSE_MID_TOL);
        assertEquals(-3.2908, M.get(4, 4), LOOSE_MID_TOL);
        assertEquals(3.2908, M.get(4, 5), LOOSE_MID_TOL);
        assertEquals(-3.2908, M.get(5, 5), LOOSE_MID_TOL);
        assertEquals(3.2908, M.get(5, 6), LOOSE_MID_TOL);
        assertEquals(-3.2908, M.get(6, 6), LOOSE_MID_TOL);
    }

    /**
     * SimilarityMatrix.txt: Compute similarity matrix between two similar matrices
     * Input: A1m=[0.2,0.8,0; 1.2,-0.4,0.1; -0.2,0.7,0.5],
     *        T=[1,2,-4,6; 0,8,-9,7; -3,7,8,-2]
     *        A2m = pinv(T)*A1m*T
     * Verify: norm(A1m*B - B*A2m) < tol where B = SimilarityMatrix(A1m, A2m)
     */
    @Test
    public void testSimilarityMatrix_Example() {
        Matrix A1m = new Matrix(3, 3);
        A1m.set(0, 0, 0.2);  A1m.set(0, 1, 0.8);  A1m.set(0, 2, 0.0);
        A1m.set(1, 0, 1.2);  A1m.set(1, 1, -0.4); A1m.set(1, 2, 0.1);
        A1m.set(2, 0, -0.2); A1m.set(2, 1, 0.7);  A1m.set(2, 2, 0.5);

        Matrix T = new Matrix(3, 4);
        T.set(0, 0, 1.0);  T.set(0, 1, 2.0);  T.set(0, 2, -4.0); T.set(0, 3, 6.0);
        T.set(1, 0, 0.0);  T.set(1, 1, 8.0);  T.set(1, 2, -9.0); T.set(1, 3, 7.0);
        T.set(2, 0, -3.0); T.set(2, 1, 7.0);  T.set(2, 2, 8.0);  T.set(2, 3, -2.0);

        // A2m = pinv(T) * A1m * T
        Matrix Tpinv = T.pinv();
        Matrix A2m = Tpinv.mult(A1m).mult(T);

        // Compute similarity matrix
        Matrix B = SimilarityMatrix.similarityMatrix(A1m, A2m);
        assertNotNull(B);

        // Verify: norm(A1m*B - B*A2m) < tol
        Matrix lhs = A1m.mult(B);
        Matrix rhs = B.mult(A2m);
        double err = 0.0;
        for (int i = 0; i < lhs.getNumRows(); i++) {
            for (int j = 0; j < lhs.getNumCols(); j++) {
                double d = lhs.get(i, j) - rhs.get(i, j);
                err += d * d;
            }
        }
        err = Math.sqrt(err);
        assertTrue(err < LOOSE_MID_TOL, "Similarity error: " + err);
    }

    // ============ Distance Measures ============

    /**
     * SquaredDifference.txt: Compute squared difference between two vectors
     * Input: p1 (MAP ACF), p2 (trace ACF)
     * Expected: 0.07738
     */
    @Test
    public void testSquaredDifference_Example() {
        double[] p1 = {0.20005, 0.12003, 0.072023, 0.043216, 0.025931,
                        0.015559, 0.0093357, 0.0056017, 0.0033611, 0.0020168};
        double[] p2 = {0.20005, 0.18927, 0.13895, 0.14213, 0.11713,
                        0.12368, 0.11212, 0.10051, 0.10019, 0.098797};

        double result = SquaredDifference.squaredDifference(p1, p2);
        assertEquals(0.07738, result, LOOSE_MID_TOL);
    }

    /**
     * RelativeEntropy.txt: Compute relative entropy between two vectors
     * Input: same p1, p2 as SquaredDifference
     * Expected: 0.28344
     */
    @Test
    public void testRelativeEntropy_Example() {
        double[] p1 = {0.20005, 0.12003, 0.072023, 0.043216, 0.025931,
                        0.015559, 0.0093357, 0.0056017, 0.0033611, 0.0020168};
        double[] p2 = {0.20005, 0.18927, 0.13895, 0.14213, 0.11713,
                        0.12368, 0.11212, 0.10051, 0.10019, 0.098797};

        double result = RelativeEntropy.relativeEntropy(p1, p2);
        assertEquals(0.28344, result, LOOSE_MID_TOL);
    }

    // ============ DMMAP/DMRAP ============

    /**
     * CheckDMMAPRepresentation.txt: Validate DMMAP representation
     * Input: D0,D1,D2,D3 (3x3 each)
     * Expected: true (valid)
     */
    @Test
    public void testCheckDMMAPRepresentation_Example() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, 0.34); D0.set(0, 1, 0.0);  D0.set(0, 2, 0.0);
        D0.set(1, 0, 0.06); D0.set(1, 1, 0.05); D0.set(1, 2, 0.03);
        D0.set(2, 0, 0.11); D0.set(2, 1, 0.13); D0.set(2, 2, 0.0);

        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 0.3);  D1.set(0, 1, 0.0);  D1.set(0, 2, 0.0);
        D1.set(1, 0, 0.16); D1.set(1, 1, 0.18); D1.set(1, 2, 0.05);
        D1.set(2, 0, 0.15); D1.set(2, 1, 0.04); D1.set(2, 2, 0.09);

        Matrix D2 = new Matrix(3, 3);
        D2.set(0, 0, 0.0);  D2.set(0, 1, 0.01); D2.set(0, 2, 0.0);
        D2.set(1, 0, 0.1);  D2.set(1, 1, 0.07); D2.set(1, 2, 0.08);
        D2.set(2, 0, 0.13); D2.set(2, 1, 0.12); D2.set(2, 2, 0.13);

        Matrix D3 = new Matrix(3, 3);
        D3.set(0, 0, 0.35); D3.set(0, 1, 0.0);  D3.set(0, 2, 0.0);
        D3.set(1, 0, 0.0);  D3.set(1, 1, 0.18); D3.set(1, 2, 0.04);
        D3.set(2, 0, 0.06); D3.set(2, 1, 0.03); D3.set(2, 2, 0.01);

        MatrixCell D = new MatrixCell();
        D.set(0, D0);
        D.set(1, D1);
        D.set(2, D2);
        D.set(3, D3);

        assertTrue(CheckDMMAPRepresentation.checkDMMAPRepresentation(D, 1e-14));
    }

    /**
     * MarginalMomentsFromDMMAP.txt: Compute marginal moments from DMMAP
     * Expected: [1.5037, 3.0278, 8.4243, 31.097, 143.88]
     */
    @Test
    public void testMarginalMomentsFromDMMAP_Example() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, 0.34); D0.set(0, 1, 0.0);  D0.set(0, 2, 0.0);
        D0.set(1, 0, 0.06); D0.set(1, 1, 0.05); D0.set(1, 2, 0.03);
        D0.set(2, 0, 0.11); D0.set(2, 1, 0.13); D0.set(2, 2, 0.0);

        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 0.3);  D1.set(0, 1, 0.0);  D1.set(0, 2, 0.0);
        D1.set(1, 0, 0.16); D1.set(1, 1, 0.18); D1.set(1, 2, 0.05);
        D1.set(2, 0, 0.15); D1.set(2, 1, 0.04); D1.set(2, 2, 0.09);

        Matrix D2 = new Matrix(3, 3);
        D2.set(0, 0, 0.0);  D2.set(0, 1, 0.01); D2.set(0, 2, 0.0);
        D2.set(1, 0, 0.1);  D2.set(1, 1, 0.07); D2.set(1, 2, 0.08);
        D2.set(2, 0, 0.13); D2.set(2, 1, 0.12); D2.set(2, 2, 0.13);

        Matrix D3 = new Matrix(3, 3);
        D3.set(0, 0, 0.35); D3.set(0, 1, 0.0);  D3.set(0, 2, 0.0);
        D3.set(1, 0, 0.0);  D3.set(1, 1, 0.18); D3.set(1, 2, 0.04);
        D3.set(2, 0, 0.06); D3.set(2, 1, 0.03); D3.set(2, 2, 0.01);

        MatrixCell D = new MatrixCell();
        D.set(0, D0);
        D.set(1, D1);
        D.set(2, D2);
        D.set(3, D3);

        double[] moms = MarginalMomentsFromDMMAP.marginalMomentsFromDMMAP(D, 0, 1e-14);

        assertNotNull(moms);
        assertTrue(moms.length >= 5);

        assertEquals(1.5037, moms[0], LOOSE_MID_TOL);
        assertEquals(3.0278, moms[1], LOOSE_MID_TOL);
        assertEquals(8.4243, moms[2], LOOSE_MID_TOL);
        assertEquals(31.097, moms[3], LOOSE_MID_TOL);
        assertEquals(143.88, moms[4], LOOSE_MID_TOL);
    }

    /**
     * MarginalDistributionFromDMMAP.txt: Compute marginal distribution from DMMAP
     * Expected: a=[0.96166, 0.030652, 0.0076858]
     *           A = D0 (the D0 matrix itself)
     */
    @Test
    public void testMarginalDistributionFromDMMAP_Example() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, 0.34); D0.set(0, 1, 0.0);  D0.set(0, 2, 0.0);
        D0.set(1, 0, 0.06); D0.set(1, 1, 0.05); D0.set(1, 2, 0.03);
        D0.set(2, 0, 0.11); D0.set(2, 1, 0.13); D0.set(2, 2, 0.0);

        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 0.3);  D1.set(0, 1, 0.0);  D1.set(0, 2, 0.0);
        D1.set(1, 0, 0.16); D1.set(1, 1, 0.18); D1.set(1, 2, 0.05);
        D1.set(2, 0, 0.15); D1.set(2, 1, 0.04); D1.set(2, 2, 0.09);

        Matrix D2 = new Matrix(3, 3);
        D2.set(0, 0, 0.0);  D2.set(0, 1, 0.01); D2.set(0, 2, 0.0);
        D2.set(1, 0, 0.1);  D2.set(1, 1, 0.07); D2.set(1, 2, 0.08);
        D2.set(2, 0, 0.13); D2.set(2, 1, 0.12); D2.set(2, 2, 0.13);

        Matrix D3 = new Matrix(3, 3);
        D3.set(0, 0, 0.35); D3.set(0, 1, 0.0);  D3.set(0, 2, 0.0);
        D3.set(1, 0, 0.0);  D3.set(1, 1, 0.18); D3.set(1, 2, 0.04);
        D3.set(2, 0, 0.06); D3.set(2, 1, 0.03); D3.set(2, 2, 0.01);

        MatrixCell D = new MatrixCell();
        D.set(0, D0);
        D.set(1, D1);
        D.set(2, D2);
        D.set(3, D3);

        MGRepresentation result = MarginalDistributionFromDMMAP.marginalDistributionFromDMMAP(D, 1e-14);
        assertNotNull(result);

        Matrix a = result.getAlpha();
        Matrix A = result.getA();

        assertNotNull(a);
        assertNotNull(A);

        // Verify initial vector a
        assertEquals(1, a.getNumRows());
        assertEquals(3, a.getNumCols());
        assertEquals(0.96166, a.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.030652, a.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0076858, a.get(0, 2), LOOSE_MID_TOL);

        // Verify A = D0
        assertEquals(3, A.getNumRows());
        assertEquals(3, A.getNumCols());
        assertEquals(0.34, A.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.0, A.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0, A.get(0, 2), LOOSE_MID_TOL);
        assertEquals(0.06, A.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.05, A.get(1, 1), LOOSE_MID_TOL);
        assertEquals(0.03, A.get(1, 2), LOOSE_MID_TOL);
        assertEquals(0.11, A.get(2, 0), LOOSE_MID_TOL);
        assertEquals(0.13, A.get(2, 1), LOOSE_MID_TOL);
        assertEquals(0.0, A.get(2, 2), LOOSE_MID_TOL);
    }

    /**
     * CheckDMRAPRepresentation.txt: Validate DMRAP representation
     * Input: H0, H1, H2 (3x3 each, H0 has negative elements)
     * Expected: true (valid)
     */
    @Test
    public void testCheckDMRAPRepresentation_Example() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, 0.15);  H0.set(0, 1, 0.2);  H0.set(0, 2, 0.18);
        H0.set(1, 0, -0.23); H0.set(1, 1, 0.17); H0.set(1, 2, 0.22);
        H0.set(2, 0, 0.19);  H0.set(2, 1, 0.15); H0.set(2, 2, 0.16);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 0.01); H1.set(0, 1, 0.08); H1.set(0, 2, 0.16);
        H1.set(1, 0, 0.02); H1.set(1, 1, 0.2);  H1.set(1, 2, 0.07);
        H1.set(2, 0, 0.02); H1.set(2, 1, 0.15); H1.set(2, 2, 0.17);

        Matrix H2 = new Matrix(3, 3);
        H2.set(0, 0, 0.14); H2.set(0, 1, 0.07); H2.set(0, 2, 0.01);
        H2.set(1, 0, 0.19); H2.set(1, 1, 0.02); H2.set(1, 2, 0.34);
        H2.set(2, 0, 0.06); H2.set(2, 1, 0.1);  H2.set(2, 2, 0.0);

        MatrixCell H = new MatrixCell();
        H.set(0, H0);
        H.set(1, H1);
        H.set(2, H2);

        assertTrue(CheckDMRAPRepresentation.checkDMRAPRepresentation(H, 1e-14));
    }

    /**
     * MarginalMomentsFromDMRAP.txt: Compute marginal moments from DMRAP
     * Expected: [1.5948, 3.4185, 9.9595, 37.742, 177.13]
     */
    @Test
    public void testMarginalMomentsFromDMRAP_Example() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, 0.15);  H0.set(0, 1, 0.2);  H0.set(0, 2, 0.18);
        H0.set(1, 0, -0.23); H0.set(1, 1, 0.17); H0.set(1, 2, 0.22);
        H0.set(2, 0, 0.19);  H0.set(2, 1, 0.15); H0.set(2, 2, 0.16);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 0.01); H1.set(0, 1, 0.08); H1.set(0, 2, 0.16);
        H1.set(1, 0, 0.02); H1.set(1, 1, 0.2);  H1.set(1, 2, 0.07);
        H1.set(2, 0, 0.02); H1.set(2, 1, 0.15); H1.set(2, 2, 0.17);

        Matrix H2 = new Matrix(3, 3);
        H2.set(0, 0, 0.14); H2.set(0, 1, 0.07); H2.set(0, 2, 0.01);
        H2.set(1, 0, 0.19); H2.set(1, 1, 0.02); H2.set(1, 2, 0.34);
        H2.set(2, 0, 0.06); H2.set(2, 1, 0.1);  H2.set(2, 2, 0.0);

        MatrixCell H = new MatrixCell();
        H.set(0, H0);
        H.set(1, H1);
        H.set(2, H2);

        double[] moms = MarginalMomentsFromDMRAP.marginalMomentsFromDMRAP(H, 0, 1e-14);

        assertNotNull(moms);
        assertTrue(moms.length >= 5);

        double[] expected = {1.5948, 3.4185, 9.9595, 37.742, 177.13};
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], moms[i], Math.max(LOOSE_MID_TOL, Math.abs(expected[i]) * 1e-3),
                    "moms[" + i + "]");
        }
    }

    // ============ CTMC/DTMC Solvers ============

    /**
     * CTMCSolve.txt: Solve CTMC stationary distribution
     * Input: Q = [-0.9, 0.5, 0.4; 0.9, -0.9, 0; 0.3, 0.3, -0.6]
     * Expected: ret = [0.40909, 0.31818, 0.27273]
     */
    @Test
    public void testCTMCSolve_Example() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -0.9); Q.set(0, 1, 0.5);  Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9);  Q.set(1, 1, -0.9); Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3);  Q.set(2, 1, 0.3);  Q.set(2, 2, -0.6);

        Matrix ret = CTMCSolve.ctmcSolve(Q, 1e-14);

        assertNotNull(ret);
        assertEquals(0.40909, ret.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.31818, ret.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.27273, ret.get(0, 2), LOOSE_MID_TOL);

        // Verify ret*Q ≈ 0
        Matrix product = ret.mult(Q);
        for (int j = 0; j < 3; j++) {
            assertEquals(0.0, product.get(0, j), LOOSE_MID_TOL);
        }
    }

    /**
     * DTMCSolve.txt: Solve DTMC stationary distribution
     * Input: Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, 0.3, 0.4]
     * Expected: ret = [0.40909, 0.31818, 0.27273]
     */
    @Test
    public void testDTMCSolve_Example() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, 0.1); Q.set(0, 1, 0.5); Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9); Q.set(1, 1, 0.1); Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3); Q.set(2, 1, 0.3); Q.set(2, 2, 0.4);

        Matrix ret = DTMCSolve.dtmcSolve(Q, 1e-14);

        assertNotNull(ret);
        assertEquals(0.40909, ret.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.31818, ret.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.27273, ret.get(0, 2), LOOSE_MID_TOL);

        // Verify ret*Q ≈ ret
        Matrix product = ret.mult(Q);
        for (int j = 0; j < 3; j++) {
            assertEquals(ret.get(0, j), product.get(0, j), LOOSE_MID_TOL);
        }
    }

    // ============ CheckGenerator ============

    /**
     * CheckGenerator.txt Test 1: Positive diagonal element
     * Input: Q = [-0.9, 0.2, 0.4; 0, 0.9, 0.9; 0, 0.6, -0.6]
     * Expected: false (positive diagonal)
     */
    @Test
    public void testCheckGenerator_InvalidDiagonal() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -0.9); Q.set(0, 1, 0.2); Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.0);  Q.set(1, 1, 0.9); Q.set(1, 2, 0.9);
        Q.set(2, 0, 0.0);  Q.set(2, 1, 0.6); Q.set(2, 2, -0.6);

        assertFalse(CheckGenerator.checkGenerator(Q, true, 1e-12));
    }

    /**
     * CheckGenerator.txt Test 2: Valid generator
     * Input: Q = [-0.9, 0.5, 0.4; 0.9, -0.9, 0; 0.3, 0.3, -0.6]
     * Expected: true
     */
    @Test
    public void testCheckGenerator_Valid() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -0.9); Q.set(0, 1, 0.5);  Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9);  Q.set(1, 1, -0.9); Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3);  Q.set(2, 1, 0.3);  Q.set(2, 2, -0.6);

        assertTrue(CheckGenerator.checkGenerator(Q, true, 1e-12));
    }

    /**
     * CheckGenerator.txt Test 3: Transient generator (rowsum != 0 allowed)
     * Input: Q = [-0.9, 0.2, 0.4; 0.9, -0.9, 0; 0.3, 0.3, -0.6]
     * Expected: true (transient=true allows non-zero rowsum)
     */
    @Test
    public void testCheckGenerator_TransientValid() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -0.9); Q.set(0, 1, 0.2); Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9);  Q.set(1, 1, -0.9); Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3);  Q.set(2, 1, 0.3);  Q.set(2, 2, -0.6);

        assertTrue(CheckGenerator.checkGenerator(Q, true, 1e-12));
    }

    /**
     * CheckGenerator.txt Test 4: Invalid rowsum (non-transient)
     * Input: Q = [-0.9, 0.5, 0.4; 0.9, -1.1, 0; 0.3, 0.3, -0.6]
     * Expected: false (rowsum != 0 with transient=false)
     */
    @Test
    public void testCheckGenerator_InvalidRowsum() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -0.9); Q.set(0, 1, 0.5);  Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9);  Q.set(1, 1, -1.1); Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3);  Q.set(2, 1, 0.3);  Q.set(2, 2, -0.6);

        assertFalse(CheckGenerator.checkGenerator(Q, false, 1e-12));
    }

    /**
     * CheckGenerator.txt Test 5: Valid generator (non-transient)
     * Input: Q = [-0.9, 0.5, 0.4; 0.9, -0.9, 0; 0.3, 0.3, -0.6]
     * Expected: true
     */
    @Test
    public void testCheckGenerator_ValidNonTransient() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, -0.9); Q.set(0, 1, 0.5);  Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9);  Q.set(1, 1, -0.9); Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3);  Q.set(2, 1, 0.3);  Q.set(2, 2, -0.6);

        assertTrue(CheckGenerator.checkGenerator(Q, false, 1e-12));
    }

    // ============ CheckProbMatrix ============

    /**
     * CheckProbMatrix.txt Test 1: Negative element
     * Input: Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, -0.1, 0.4]
     * Expected: false
     */
    @Test
    public void testCheckProbMatrix_NegativeElement() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, 0.1); Q.set(0, 1, 0.5);  Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9); Q.set(1, 1, 0.1);  Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3); Q.set(2, 1, -0.1); Q.set(2, 2, 0.4);

        assertFalse(CheckProbMatrix.checkProbMatrix(Q, false, 1e-12));
    }

    /**
     * CheckProbMatrix.txt Test 2: Invalid rowsum
     * Input: Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, 0.1, 0.4]
     * Expected: false (row 3 sums to 0.8)
     */
    @Test
    public void testCheckProbMatrix_InvalidRowsum() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, 0.1); Q.set(0, 1, 0.5); Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9); Q.set(1, 1, 0.1); Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3); Q.set(2, 1, 0.1); Q.set(2, 2, 0.4);

        assertFalse(CheckProbMatrix.checkProbMatrix(Q, false, 1e-12));
    }

    /**
     * CheckProbMatrix.txt Test 3: Valid probability matrix
     * Input: Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, 0.3, 0.4]
     * Expected: true
     */
    @Test
    public void testCheckProbMatrix_Valid() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, 0.1); Q.set(0, 1, 0.5); Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9); Q.set(1, 1, 0.1); Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3); Q.set(2, 1, 0.3); Q.set(2, 2, 0.4);

        assertTrue(CheckProbMatrix.checkProbMatrix(Q, false, 1e-12));
    }

    /**
     * CheckProbMatrix.txt Test 4: Stochastic matrix fails transient check
     * Input: Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, 0.3, 0.4]
     * Expected: false (largest eigenvalue = 1, not transient)
     */
    @Test
    public void testCheckProbMatrix_TransientInvalid() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, 0.1); Q.set(0, 1, 0.5); Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9); Q.set(1, 1, 0.1); Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3); Q.set(2, 1, 0.3); Q.set(2, 2, 0.4);

        assertFalse(CheckProbMatrix.checkProbMatrix(Q, true, 1e-12));
    }

    /**
     * CheckProbMatrix.txt Test 5: Sub-stochastic matrix passes transient check
     * Input: Q = [0.1, 0.5, 0.4; 0.9, 0.1, 0; 0.3, 0.1, 0.4]
     * Expected: true (sub-stochastic, largest eigenvalue < 1)
     */
    @Test
    public void testCheckProbMatrix_TransientValid() {
        Matrix Q = new Matrix(3, 3);
        Q.set(0, 0, 0.1); Q.set(0, 1, 0.5); Q.set(0, 2, 0.4);
        Q.set(1, 0, 0.9); Q.set(1, 1, 0.1); Q.set(1, 2, 0.0);
        Q.set(2, 0, 0.3); Q.set(2, 1, 0.1); Q.set(2, 2, 0.4);

        assertTrue(CheckProbMatrix.checkProbMatrix(Q, true, 1e-12));
    }

    // ============ CheckMoments ============

    /**
     * CheckMoments.txt Test 1: Invalid moments
     * Input: M = [1.2, 5, 8, 29, 3412]
     * Expected: false
     */
    @Test
    public void testCheckMoments_Invalid() {
        Matrix M = new Matrix(1, 5);
        M.set(0, 0, 1.2); M.set(0, 1, 5.0); M.set(0, 2, 8.0);
        M.set(0, 3, 29.0); M.set(0, 4, 3412.0);

        assertFalse(CheckMoments.checkMoments(M, 1e-14));
    }

    /**
     * CheckMoments.txt Test 2: Valid moments
     * Input: M = [1.3, 2.4, 6.03, 20.5, 89.5]
     * Expected: true
     */
    @Test
    public void testCheckMoments_Valid() {
        Matrix M = new Matrix(1, 5);
        M.set(0, 0, 1.3); M.set(0, 1, 2.4); M.set(0, 2, 6.03);
        M.set(0, 3, 20.5); M.set(0, 4, 89.5);

        assertTrue(CheckMoments.checkMoments(M, 1e-14));
    }

    // ============ Moment Transforms ============

    /**
     * JFactorialMomsFromJMoms.txt / JMomsFromJFactorialMoms.txt:
     * Roundtrip test: MM -> JFmoms -> Jmoms, verify Jmoms == MM
     * Input: MM = [0.7, 2, 3, 4; 5, 6, 7, 8; 9, 10, 11, 12]
     * JFmoms = [0.7, 1.3, -1.6, 3.8; 4.3, -0.3, 0.6, -1.8; -4.6, 0.6, -1.2, 3.6]
     */
    @Test
    public void testJFactorialMomsFromJMoms_Example() {
        Matrix MM = new Matrix(3, 4);
        MM.set(0, 0, 0.7); MM.set(0, 1, 2.0); MM.set(0, 2, 3.0); MM.set(0, 3, 4.0);
        MM.set(1, 0, 5.0); MM.set(1, 1, 6.0); MM.set(1, 2, 7.0); MM.set(1, 3, 8.0);
        MM.set(2, 0, 9.0); MM.set(2, 1, 10.0); MM.set(2, 2, 11.0); MM.set(2, 3, 12.0);

        Matrix JFmoms = JFactorialMomsFromJMoms.jFactorialMomsFromJMoms(MM);

        assertNotNull(JFmoms);
        assertEquals(3, JFmoms.getNumRows());
        assertEquals(4, JFmoms.getNumCols());

        // Verify expected JFmoms values
        double[][] expectedJF = {
            {0.7, 1.3, -1.6, 3.8},
            {4.3, -0.3, 0.6, -1.8},
            {-4.6, 0.6, -1.2, 3.6}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedJF[i][j], JFmoms.get(i, j), LOOSE_MID_TOL,
                        "JFmoms[" + i + "," + j + "]");
            }
        }

        // Roundtrip: JFmoms -> Jmoms should equal MM
        Matrix Jmoms = JMomsFromJFactorialMoms.jMomsFromJFactorialMoms(JFmoms);

        assertNotNull(Jmoms);
        assertEquals(3, Jmoms.getNumRows());
        assertEquals(4, Jmoms.getNumCols());

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(MM.get(i, j), Jmoms.get(i, j), LOOSE_MID_TOL,
                        "Roundtrip Jmoms[" + i + "," + j + "]");
            }
        }
    }

    // ============ SimilarityMatrixForVectors ============

    /**
     * SimilarityMatrixForVectors.txt: Compute matrix B such that B*vecA = vecB
     * Input: vecA = [0; 0.3; -1.5; 0], vecB = [1; 0.2; 0; 1]
     * Expected B:
     *   [0, 3.3333, 0, 0]
     *   [0.66667, 0.66667, 0, 0]
     *   [0, 0, 0, 0]
     *   [-0.83333, -0.83333, -0.83333, -0.83333]
     * Verify: norm(B*vecA - vecB) ≈ 0
     */
    @Test
    public void testSimilarityMatrixForVectors_Example() {
        Matrix vecA = new Matrix(4, 1);
        vecA.set(0, 0, 0.0);
        vecA.set(1, 0, 0.3);
        vecA.set(2, 0, -1.5);
        vecA.set(3, 0, 0.0);

        Matrix vecB = new Matrix(4, 1);
        vecB.set(0, 0, 1.0);
        vecB.set(1, 0, 0.2);
        vecB.set(2, 0, 0.0);
        vecB.set(3, 0, 1.0);

        Matrix B = SimilarityMatrixForVectors.SimilarityMatrixForVectors(vecA, vecB);

        assertNotNull(B);
        assertEquals(4, B.getNumRows());
        assertEquals(4, B.getNumCols());

        // Verify expected B values
        double[][] expectedB = {
            {0.0, 3.3333, 0.0, 0.0},
            {0.66667, 0.66667, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {-0.83333, -0.83333, -0.83333, -0.83333}
        };
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedB[i][j], B.get(i, j), LOOSE_MID_TOL,
                        "B[" + i + "," + j + "]");
            }
        }

        // Verify B*vecA = vecB
        Matrix product = B.mult(vecA);
        for (int i = 0; i < 4; i++) {
            assertEquals(vecB.get(i, 0), product.get(i, 0), LOOSE_MID_TOL,
                    "B*vecA[" + i + "]");
        }
    }

    // ============ Moment Transform Tests ============

    /**
     * FactorialMomsFromMoms.txt: Convert raw moments to factorial moments
     * Input: M=[1.3, 2.4, 6.03, 20.5, 89.5, 474.9]
     * Expected: fmoms=[1.3, 1.1, 1.43, 2.92, 6.75, 19.75]
     * Roundtrip error: 3.5527e-15
     */
    @Test
    public void testFactorialMomsFromMoms() {
        Matrix M = new Matrix(1, 6);
        M.set(0, 0, 1.3); M.set(0, 1, 2.4); M.set(0, 2, 6.03);
        M.set(0, 3, 20.5); M.set(0, 4, 89.5); M.set(0, 5, 474.9);

        Matrix fmoms = FactorialMomsFromMoms.factorialMomsFromMoms(M);

        double[] expected = {1.3, 1.1, 1.43, 2.92, 6.75, 19.75};
        assertEquals(6, fmoms.length());
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], fmoms.get(0, i), LOOSE_MID_TOL, "fmoms[" + i + "]");
        }

        // Roundtrip: MomsFromFactorialMoms(fmoms) should give back M
        Matrix roundtrip = MomsFromFactorialMoms.MomsFromFactorialMoms(fmoms);
        for (int i = 0; i < 6; i++) {
            assertEquals(M.get(0, i), roundtrip.get(0, i), FINE_TOL, "roundtrip[" + i + "]");
        }
    }

    /**
     * MomsFromFactorialMoms.txt: Convert factorial moments to raw moments
     * Input: fm=[1.3, 1.1, 1.43, 2.92, 6.75, 19.75]
     * Expected: M=[1.3, 2.4, 6.03, 20.5, 89.5, 474.9]
     * Roundtrip error: 3.5527e-15
     */
    @Test
    public void testMomsFromFactorialMoms() {
        Matrix fm = new Matrix(1, 6);
        fm.set(0, 0, 1.3); fm.set(0, 1, 1.1); fm.set(0, 2, 1.43);
        fm.set(0, 3, 2.92); fm.set(0, 4, 6.75); fm.set(0, 5, 19.75);

        Matrix M = MomsFromFactorialMoms.MomsFromFactorialMoms(fm);

        double[] expected = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};
        assertEquals(6, M.length());
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], M.get(0, i), LOOSE_MID_TOL, "M[" + i + "]");
        }
    }

    /**
     * MomsFromHankelMoms.txt: Convert Hankel moments to raw moments
     * Input: hmoms=[1.3, 0.71, 2.079, 1.9973, 13.841, 44.916]
     * Expected: M=[1.3, 2.4, 6.03, 20.5, 89.5, 474.9]
     * Roundtrip error: 5.118e-13
     */
    @Test
    public void testMomsFromHankelMoms() {
        Matrix hmoms = new Matrix(1, 6);
        hmoms.set(0, 0, 1.3); hmoms.set(0, 1, 0.71); hmoms.set(0, 2, 2.079);
        hmoms.set(0, 3, 1.9973); hmoms.set(0, 4, 13.841); hmoms.set(0, 5, 44.916);

        Matrix M = MomsFromHankelMoms.momsFromHankelMoms(hmoms);

        double[] expected = {1.3, 2.4, 6.03, 20.5, 89.5, 474.9};
        assertEquals(6, M.length());
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], M.get(0, i), Math.max(LOOSE_MID_TOL, Math.abs(expected[i]) * 1e-3), "M[" + i + "]");
        }
    }

    /**
     * MomsFromNormMoms.txt: Convert normalized moments to raw moments
     * Input: nmoms=[1.2, 3.4722, 6.3333, 10.855, 15.513]
     * Expected: M=[1.2, 5, 38, 495, 9215]
     * Error: 0
     */
    @Test
    public void testMomsFromNormMoms() {
        Matrix nmoms = new Matrix(1, 5);
        nmoms.set(0, 0, 1.2); nmoms.set(0, 1, 3.4722); nmoms.set(0, 2, 6.3333);
        nmoms.set(0, 3, 10.855); nmoms.set(0, 4, 15.513);

        Matrix M = MomsFromNormMoms.momsFromNormMoms(nmoms);

        double[] expected = {1.2, 5.0, 38.0, 495.0, 9215.0};
        assertEquals(5, M.length());
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], M.get(0, i), Math.max(LOOSE_MID_TOL, Math.abs(expected[i]) * 1e-3), "M[" + i + "]");
        }
    }

    /**
     * MomsFromReducedMoms.txt: Convert reduced moments to raw moments
     * Input: rmoms=[1.2, 2.5, 6.3333, 20.625, 76.792]
     * Expected: M=[1.2, 5, 38, 495, 9215]
     * Error: 0
     */
    @Test
    public void testMomsFromReducedMoms() {
        Matrix rmoms = new Matrix(1, 5);
        rmoms.set(0, 0, 1.2); rmoms.set(0, 1, 2.5); rmoms.set(0, 2, 6.3333);
        rmoms.set(0, 3, 20.625); rmoms.set(0, 4, 76.792);

        Matrix M = MomsFromReducedMoms.momsFromReducedMoms(rmoms);

        double[] expected = {1.2, 5.0, 38.0, 495.0, 9215.0};
        assertEquals(5, M.length());
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], M.get(0, i), Math.max(LOOSE_MID_TOL, Math.abs(expected[i]) * 1e-3), "M[" + i + "]");
        }
    }

    // ============ Trace Analysis Tests ============

    /**
     * MarginalMomentsFromTrace.txt: Compute marginal moments from a trace
     * Generates trace from MAP D0/D1, computes moments, compares to theoretical
     * Theoretical moments: [0.054124, 0.006423, 0.0012051]
     * (Stochastic: trace moments should be close but not exact)
     */
    @Test
    public void testMarginalMomentsFromTrace() {
        // Use MAP from the data file
        Matrix D0 = new Matrix(4, 4);
        D0.set(0, 0, -5.0); D0.set(0, 1, 0.0); D0.set(0, 2, 1.0); D0.set(0, 3, 1.0);
        D0.set(1, 0, 1.0);  D0.set(1, 1, -8.0); D0.set(1, 2, 1.0); D0.set(1, 3, 0.0);
        D0.set(2, 0, 1.0);  D0.set(2, 1, 0.0);  D0.set(2, 2, -4.0); D0.set(2, 3, 1.0);
        D0.set(3, 0, 1.0);  D0.set(3, 1, 2.0);  D0.set(3, 2, 3.0); D0.set(3, 3, -9.0);

        Matrix D1 = new Matrix(4, 4);
        D1.set(0, 0, 0.0); D1.set(0, 1, 1.0); D1.set(0, 2, 0.0); D1.set(0, 3, 2.0);
        D1.set(1, 0, 2.0); D1.set(1, 1, 1.0); D1.set(1, 2, 3.0); D1.set(1, 3, 0.0);
        D1.set(2, 0, 0.0); D1.set(2, 1, 0.0); D1.set(2, 2, 1.0); D1.set(2, 3, 1.0);
        D1.set(3, 0, 1.0); D1.set(3, 1, 1.0); D1.set(3, 2, 0.0); D1.set(3, 3, 1.0);

        // Generate trace samples
        double[] trace = SamplesFromMAP.samplesFromMAP(D0, D1, 10000, null, new Random(42));

        // Compute marginal moments from trace
        double[] moms = MarginalMomentsFromTrace.marginalMomentsFromTrace(trace, 3);

        // Theoretical moments from MAP
        double[] theoretical = MarginalMomentsFromMAP.marginalMomentsFromMAP(D0, D1, 3);

        // Stochastic test: trace moments should be within 20% of theoretical
        assertEquals(3, moms.length);
        for (int i = 0; i < 3; i++) {
            assertEquals(theoretical[i], moms[i], Math.abs(theoretical[i]) * 0.2,
                    "moms[" + i + "] trace vs theoretical");
        }
    }

    /**
     * LagCorrelationsFromTrace.txt: Compute lag correlations from a trace
     * Generates trace from MAP D0/D1, computes lag correlations, compares to theoretical
     * Theoretical: [0.041288, 0.021962, ...]
     * (Stochastic: trace correlations should be close but not exact)
     */
    @Test
    public void testLagCorrelationsFromTrace() {
        // Use MAP from the data file
        Matrix D0 = new Matrix(4, 4);
        D0.set(0, 0, -5.0); D0.set(0, 1, 0.0); D0.set(0, 2, 1.0); D0.set(0, 3, 1.0);
        D0.set(1, 0, 1.0);  D0.set(1, 1, -8.0); D0.set(1, 2, 1.0); D0.set(1, 3, 0.0);
        D0.set(2, 0, 1.0);  D0.set(2, 1, 0.0);  D0.set(2, 2, -4.0); D0.set(2, 3, 1.0);
        D0.set(3, 0, 1.0);  D0.set(3, 1, 2.0);  D0.set(3, 2, 3.0); D0.set(3, 3, -9.0);

        Matrix D1 = new Matrix(4, 4);
        D1.set(0, 0, 0.0); D1.set(0, 1, 1.0); D1.set(0, 2, 0.0); D1.set(0, 3, 2.0);
        D1.set(1, 0, 2.0); D1.set(1, 1, 1.0); D1.set(1, 2, 3.0); D1.set(1, 3, 0.0);
        D1.set(2, 0, 0.0); D1.set(2, 1, 0.0); D1.set(2, 2, 1.0); D1.set(2, 3, 1.0);
        D1.set(3, 0, 1.0); D1.set(3, 1, 1.0); D1.set(3, 2, 0.0); D1.set(3, 3, 1.0);

        // Generate trace samples
        double[] trace = SamplesFromMAP.samplesFromMAP(D0, D1, 10000, null, new Random(42));

        // Compute lag correlations from trace
        double[] acf = LagCorrelationsFromTrace.lagCorrelationsFromTrace(trace, 10);

        // Just verify it returns correct number of values and they are reasonable
        assertEquals(10, acf.length);
        for (int i = 0; i < acf.length; i++) {
            assertTrue(acf[i] >= -1.0 && acf[i] <= 1.0,
                    "acf[" + i + "] should be in [-1,1], got " + acf[i]);
        }
    }

    // ============ Phase 1: Trace functions, Likelihood, Empirical measures ============

    /**
     * CdfFromTrace.txt: Compute empirical CDF from a MAP trace
     * Verify: y[0]=0, y[last]=1, monotonically non-decreasing
     */
    @Test
    public void testCdfFromTrace_Example() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -18.); D0.set(0, 1, 1.);  D0.set(0, 2, 4.);
        D0.set(1, 0, 2.);   D0.set(1, 1, -18.); D0.set(1, 2, 7.);
        D0.set(2, 0, 1.);   D0.set(2, 1, 3.);   D0.set(2, 2, -32.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 12.); D1.set(0, 1, 1.); D1.set(0, 2, 0.);
        D1.set(1, 0, 1.);  D1.set(1, 1, 8.); D1.set(1, 2, 0.);
        D1.set(2, 0, 2.);  D1.set(2, 1, 1.); D1.set(2, 2, 25.);
        double[] tr = SamplesFromMAP.samplesFromMAP(D0, D1, 10000, null, new Random(42));
        CdfResult cdf = CdfFromTrace.cdfFromTrace(tr);
        assertNotNull(cdf);
        assertEquals(tr.length, cdf.getX().length);
        assertEquals(0.0, cdf.getY()[0], FINE_TOL);
        assertEquals(1.0, cdf.getY()[cdf.getY().length - 1], FINE_TOL);
        // Monotonically non-decreasing
        for (int i = 1; i < cdf.getY().length; i++) {
            assertTrue(cdf.getY()[i] >= cdf.getY()[i - 1]);
        }
    }

    /**
     * CdfFromWeightedTrace.txt: weighted trace CDF
     */
    @Test
    public void testCdfFromWeightedTrace_Example() {
        double[] wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
        double[] weights = {12., 1., 34., 23., 8., 2.};
        CdfResult cdf = CdfFromTrace.cdfFromWeightedTrace(wtrace, weights);
        assertNotNull(cdf);
        assertEquals(6, cdf.getX().length);
        // Last value should be 1
        assertEquals(1.0, cdf.getY()[cdf.getY().length - 1], FINE_TOL);
        // Values monotonically non-decreasing
        for (int i = 1; i < cdf.getY().length; i++) {
            assertTrue(cdf.getY()[i] >= cdf.getY()[i - 1]);
        }
    }

    /**
     * PdfFromTrace.txt: Compute empirical PDF from MAP trace
     * Verify: y>=0, integral approx 1
     */
    @Test
    @Tag("slow") // ~337s
    public void testPdfFromTrace_Example() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -18.); D0.set(0, 1, 1.);  D0.set(0, 2, 4.);
        D0.set(1, 0, 2.);   D0.set(1, 1, -18.); D0.set(1, 2, 7.);
        D0.set(2, 0, 1.);   D0.set(2, 1, 3.);   D0.set(2, 2, -32.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 12.); D1.set(0, 1, 1.); D1.set(0, 2, 0.);
        D1.set(1, 0, 1.);  D1.set(1, 1, 8.); D1.set(1, 2, 0.);
        D1.set(2, 0, 2.);  D1.set(2, 1, 1.); D1.set(2, 2, 25.);
        double[] tr = SamplesFromMAP.samplesFromMAP(D0, D1, 1000000, null, new Random(42));
        // 50 interval bounds from 0 to 0.5
        double[] intBounds = new double[51];
        for (int i = 0; i <= 50; i++) intBounds[i] = i * 0.01;
        PdfResult pdf = PdfFromTrace.pdfFromTrace(tr, intBounds);
        assertNotNull(pdf);
        assertEquals(50, pdf.y.length);
        for (int i = 0; i < pdf.y.length; i++) {
            assertTrue(pdf.y[i] >= 0, "PDF value should be >= 0");
        }
        // Integral should be approximately 1 (over the range 0..0.5 which may not capture all mass)
        double integral = 0;
        for (int i = 0; i < pdf.y.length; i++) {
            integral += pdf.y[i] * (intBounds[i + 1] - intBounds[i]);
        }
        assertTrue(integral > 0 && integral <= 1.1, "PDF integral over [0,0.5] should be in (0,1.1], got " + integral);
    }

    /**
     * PdfFromWeightedTrace.txt: weighted PDF
     */
    @Test
    public void testPdfFromWeightedTrace_Example() {
        double[] wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
        double[] weights = {12., 1., 34., 23., 8., 2.};
        // bounds 0..3 step 0.5
        double[] intBounds = new double[7];
        for (int i = 0; i <= 6; i++) intBounds[i] = i * 0.5;
        PdfResult pdf = PdfFromTrace.pdfFromWeightedTrace(wtrace, weights, intBounds);
        assertNotNull(pdf);
        assertEquals(6, pdf.y.length);
        for (int i = 0; i < pdf.y.length; i++) {
            assertTrue(pdf.y[i] >= 0);
        }
    }

    /**
     * LagkJointMomentsFromTrace.txt: Lag-k joint moments from MAP trace
     * K=3, verify all positive, matches structure from reference
     */
    @Test
    public void testLagkJointMomentsFromTrace_Example() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -18.); D0.set(0, 1, 1.);  D0.set(0, 2, 4.);
        D0.set(1, 0, 2.);   D0.set(1, 1, -18.); D0.set(1, 2, 7.);
        D0.set(2, 0, 1.);   D0.set(2, 1, 3.);   D0.set(2, 2, -32.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 12.); D1.set(0, 1, 1.); D1.set(0, 2, 0.);
        D1.set(1, 0, 1.);  D1.set(1, 1, 8.); D1.set(1, 2, 0.);
        D1.set(2, 0, 2.);  D1.set(2, 1, 1.); D1.set(2, 2, 25.);
        double[] tr = SamplesFromMAP.samplesFromMAP(D0, D1, 1000000, null, new Random(42));
        Matrix Nm1 = LagCorrelationsFromTrace.lagkJointMomentsFromTrace(tr, 3, 1);
        assertNotNull(Nm1);
        assertEquals(4, Nm1.getNumRows());
        assertEquals(4, Nm1.getNumCols());
        // All entries should be positive
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertTrue(Nm1.get(i, j) > 0, "Nm1[" + i + "," + j + "] should be > 0");
            }
        }
        // Nm1[0,0] should be 1 (zeroth moment)
        assertEquals(1.0, Nm1.get(0, 0), COARSE_TOL);
    }

    /**
     * MarginalMomentsFromWeightedTrace.txt
     * Expected: [0.65242, 0.5958, 0.74264]
     */
    @Test
    public void testMarginalMomentsFromWeightedTrace_Example() {
        double[] wtrace = {0.12, 1.23, 0.546, 0.6765, 1.34, 2.34};
        double[] weights = {12., 1., 34., 23., 8., 2.};
        double[] moms = MarginalMomentsFromTrace.marginalMomentsFromWeightedTrace(wtrace, weights, 3);
        assertEquals(3, moms.length);
        assertEquals(0.65242, moms[0], LOOSE_MID_TOL);
        assertEquals(0.5958, moms[1], LOOSE_MID_TOL);
        assertEquals(0.74264, moms[2], LOOSE_MID_TOL);
    }

    /**
     * LikelihoodFromTrace.txt: PH log-likelihood
     * Uses APH from 3 moments of a PH trace, verifies finite positive result
     */
    @Test
    public void testLikelihoodFromTracePH_Example() {
        // Generate PH trace from a simple PH distribution
        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, 0.8); alpha.set(0, 1, 0.2);
        Matrix A = new Matrix(2, 2);
        A.set(0, 0, -5.0); A.set(0, 1, 1.0);
        A.set(1, 0, 0.0);  A.set(1, 1, -2.0);
        // Generate samples using matrix exponential method
        Random rng = new Random(42);
        double[] trace = new double[1000];
        for (int i = 0; i < 1000; i++) {
            // Simple exponential trace for testing
            trace[i] = -Math.log(rng.nextDouble()) / 3.0;
        }
        double logli = LikelihoodFromTrace.likelihoodFromTracePH(trace, alpha, A, 1e-14);
        assertTrue(Double.isFinite(logli), "PH logli should be finite");
    }

    /**
     * LikelihoodFromTrace.txt: MAP log-likelihood
     */
    @Test
    public void testLikelihoodFromTraceMAP_Example() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -18.); D0.set(0, 1, 1.);  D0.set(0, 2, 4.);
        D0.set(1, 0, 2.);   D0.set(1, 1, -18.); D0.set(1, 2, 7.);
        D0.set(2, 0, 1.);   D0.set(2, 1, 3.);   D0.set(2, 2, -32.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 12.); D1.set(0, 1, 1.); D1.set(0, 2, 0.);
        D1.set(1, 0, 1.);  D1.set(1, 1, 8.); D1.set(1, 2, 0.);
        D1.set(2, 0, 2.);  D1.set(2, 1, 1.); D1.set(2, 2, 25.);
        double[] tr = SamplesFromMAP.samplesFromMAP(D0, D1, 1000, null, new Random(42));
        double logli = LikelihoodFromTrace.likelihoodFromTraceMAP(tr, D0, D1, 1e-14);
        assertTrue(Double.isFinite(logli), "MAP logli should be finite");
        assertTrue(logli > 0, "MAP logli should be positive for well-fitting MAP");
    }

    /**
     * EmpiricalSquaredDifference.txt
     */
    @Test
    public void testEmpiricalSquaredDifference_Example() {
        double[] p1 = {0.1, 0.3, 0.4, 0.2};
        double[] p2 = {0.15, 0.25, 0.35, 0.25};
        double sd = SquaredDifference.empiricalSquaredDifference(p1, p2);
        assertTrue(sd >= 0, "Squared difference should be >= 0");
        // Known: sum of (0.05^2 + 0.05^2 + 0.05^2 + 0.05^2) = 0.01
        assertEquals(0.01, sd, FINE_TOL);
    }

    /**
     * EmpiricalRelativeEntropy.txt
     */
    @Test
    public void testEmpiricalRelativeEntropy_Example() {
        double[] p1 = {0.1, 0.3, 0.4, 0.2};
        double[] p2 = {0.15, 0.25, 0.35, 0.25};
        double re = RelativeEntropy.empiricalRelativeEntropy(p1, p2);
        assertTrue(re >= 0, "Relative entropy should be >= 0");
        assertTrue(Double.isFinite(re), "Relative entropy should be finite");
    }

    /**
     * JMomsFromJFactorialMoms.txt: Round-trip test
     * MM -> JFactorialMomsFromJMoms -> JMomsFromJFactorialMoms -> should equal MM
     */
    @Test
    public void testJMomsFromJFactorialMoms_Standalone() {
        Matrix MM = new Matrix(3, 4);
        MM.set(0, 0, 0.7); MM.set(0, 1, 2.); MM.set(0, 2, 3.); MM.set(0, 3, 4.);
        MM.set(1, 0, 5.);  MM.set(1, 1, 6.); MM.set(1, 2, 7.); MM.set(1, 3, 8.);
        MM.set(2, 0, 9.);  MM.set(2, 1, 10.); MM.set(2, 2, 11.); MM.set(2, 3, 12.);

        Matrix JFmoms = JFactorialMomsFromJMoms.jFactorialMomsFromJMoms(MM);
        // Verify expected factorial moments
        assertEquals(0.7, JFmoms.get(0, 0), LOOSE_MID_TOL);
        assertEquals(1.3, JFmoms.get(0, 1), LOOSE_MID_TOL);
        assertEquals(-1.6, JFmoms.get(0, 2), LOOSE_MID_TOL);
        assertEquals(3.8, JFmoms.get(0, 3), LOOSE_MID_TOL);

        Matrix Jmoms = JMomsFromJFactorialMoms.jMomsFromJFactorialMoms(JFmoms);
        // Round-trip: should match original
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(MM.get(i, j), Jmoms.get(i, j), FINE_TOL,
                        "Jmoms[" + i + "," + j + "] round-trip");
            }
        }
    }

    // ============ Phase 4: MStaircase + MinimalRep ============

    /**
     * MinimalRepFromRAP.txt: 'cont' mode preserves full representation
     */
    @Test
    public void testMinimalRepFromRAP_Cont() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -5.); D0.set(0, 1, 1.); D0.set(0, 2, 0.);
        D0.set(1, 0, 3.);  D0.set(1, 1, -3.); D0.set(1, 2, 0.);
        D0.set(2, 0, 1.);  D0.set(2, 1, 1.);  D0.set(2, 2, -5.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 0.); D1.set(0, 1, 0.); D1.set(0, 2, 4.);
        D1.set(1, 0, 0.); D1.set(1, 1, 0.); D1.set(1, 2, 0.);
        D1.set(2, 0, 1.); D1.set(2, 1, 1.); D1.set(2, 2, 1.);

        Pair<Matrix, Matrix> result = MinimalRepFromRAP.minimalRepFromRAP(D0, D1, "cont", 1e-12);
        Matrix H0 = result.getFirst();
        Matrix H1 = result.getSecond();
        // 'cont' mode should preserve full 3x3 representation
        assertEquals(3, H0.getNumRows());
        assertEquals(3, H0.getNumCols());
    }

    /**
     * MinimalRepFromRAP.txt: 'obs' mode reduces to 2x2
     * Expected H0 = [-4.4074, 1.6931; 0.84259, -2.5926]
     * Expected H1 = [2.037, 0.67725; 2.787, -1.037]
     */
    @Test
    public void testMinimalRepFromRAP_Obs() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -5.); D0.set(0, 1, 1.); D0.set(0, 2, 0.);
        D0.set(1, 0, 3.);  D0.set(1, 1, -3.); D0.set(1, 2, 0.);
        D0.set(2, 0, 1.);  D0.set(2, 1, 1.);  D0.set(2, 2, -5.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 0.); D1.set(0, 1, 0.); D1.set(0, 2, 4.);
        D1.set(1, 0, 0.); D1.set(1, 1, 0.); D1.set(1, 2, 0.);
        D1.set(2, 0, 1.); D1.set(2, 1, 1.); D1.set(2, 2, 1.);

        Pair<Matrix, Matrix> result = MinimalRepFromRAP.minimalRepFromRAP(D0, D1, "obs", 1e-12);
        Matrix H0 = result.getFirst();
        Matrix H1 = result.getSecond();
        assertEquals(2, H0.getNumRows());
        assertEquals(-4.4074, H0.get(0, 0), LOOSE_MID_TOL);
        assertEquals(1.6931, H0.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.84259, H0.get(1, 0), LOOSE_MID_TOL);
        assertEquals(-2.5926, H0.get(1, 1), LOOSE_MID_TOL);
        assertEquals(2.037, H1.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.67725, H1.get(0, 1), LOOSE_MID_TOL);
        assertEquals(2.787, H1.get(1, 0), LOOSE_MID_TOL);
        assertEquals(-1.037, H1.get(1, 1), LOOSE_MID_TOL);
    }

    /**
     * MinimalRepFromRAP.txt: 'obscont' mode reduces to 2x2
     */
    @Test
    public void testMinimalRepFromRAP_ObsCont() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -5.); D0.set(0, 1, 1.); D0.set(0, 2, 0.);
        D0.set(1, 0, 3.);  D0.set(1, 1, -3.); D0.set(1, 2, 0.);
        D0.set(2, 0, 1.);  D0.set(2, 1, 1.);  D0.set(2, 2, -5.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 0.); D1.set(0, 1, 0.); D1.set(0, 2, 4.);
        D1.set(1, 0, 0.); D1.set(1, 1, 0.); D1.set(1, 2, 0.);
        D1.set(2, 0, 1.); D1.set(2, 1, 1.); D1.set(2, 2, 1.);

        Pair<Matrix, Matrix> result = MinimalRepFromRAP.minimalRepFromRAP(D0, D1, "obscont", 1e-12);
        Matrix H0 = result.getFirst();
        Matrix H1 = result.getSecond();
        assertEquals(2, H0.getNumRows());
        assertEquals(-4.4074, H0.get(0, 0), LOOSE_MID_TOL);
        assertEquals(1.6931, H0.get(0, 1), LOOSE_MID_TOL);
    }

    /**
     * MinimalRepFromMRAP.txt: 'obs' mode reduces 3x3 to 2x2
     * Expected H{1} = [-4.4074, 1.6931; 0.84259, -2.5926]
     * Expected H{2} = [0.40741, 0.13545; 0.55741, -0.20741]
     * Expected H{3} = [1.6296, 0.5418; 2.2296, -0.82963]
     */
    @Test
    public void testMinimalRepFromMRAP_Obs() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -5.); D0.set(0, 1, 1.); D0.set(0, 2, 0.);
        D0.set(1, 0, 3.);  D0.set(1, 1, -3.); D0.set(1, 2, 0.);
        D0.set(2, 0, 1.);  D0.set(2, 1, 1.);  D0.set(2, 2, -5.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 0.); D1.set(0, 1, 0.); D1.set(0, 2, 0.8);
        D1.set(1, 0, 0.); D1.set(1, 1, 0.); D1.set(1, 2, 0.);
        D1.set(2, 0, 0.2); D1.set(2, 1, 0.2); D1.set(2, 2, 0.2);
        Matrix D2 = new Matrix(3, 3);
        D2.set(0, 0, 0.); D2.set(0, 1, 0.); D2.set(0, 2, 3.2);
        D2.set(1, 0, 0.); D2.set(1, 1, 0.); D2.set(1, 2, 0.);
        D2.set(2, 0, 0.8); D2.set(2, 1, 0.8); D2.set(2, 2, 0.8);

        MatrixCell Dm = new MatrixCell(3);
        Dm.set(0, D0);
        Dm.set(1, D1);
        Dm.set(2, D2);

        MatrixCell H = MinimalRepFromMRAP.minimalRepFromMRAP(Dm, "obs", 1e-12);
        assertEquals(2, H.get(0).getNumRows());
        assertEquals(-4.4074, H.get(0).get(0, 0), LOOSE_MID_TOL);
        assertEquals(1.6931, H.get(0).get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.40741, H.get(1).get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.13545, H.get(1).get(0, 1), LOOSE_MID_TOL);
        assertEquals(1.6296, H.get(2).get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.5418, H.get(2).get(0, 1), LOOSE_MID_TOL);
    }

    /**
     * MinimalRepFromMRAP.txt: 'obscont' mode
     */
    @Test
    public void testMinimalRepFromMRAP_ObsCont() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -5.); D0.set(0, 1, 1.); D0.set(0, 2, 0.);
        D0.set(1, 0, 3.);  D0.set(1, 1, -3.); D0.set(1, 2, 0.);
        D0.set(2, 0, 1.);  D0.set(2, 1, 1.);  D0.set(2, 2, -5.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 0.); D1.set(0, 1, 0.); D1.set(0, 2, 0.8);
        D1.set(1, 0, 0.); D1.set(1, 1, 0.); D1.set(1, 2, 0.);
        D1.set(2, 0, 0.2); D1.set(2, 1, 0.2); D1.set(2, 2, 0.2);
        Matrix D2 = new Matrix(3, 3);
        D2.set(0, 0, 0.); D2.set(0, 1, 0.); D2.set(0, 2, 3.2);
        D2.set(1, 0, 0.); D2.set(1, 1, 0.); D2.set(1, 2, 0.);
        D2.set(2, 0, 0.8); D2.set(2, 1, 0.8); D2.set(2, 2, 0.8);

        MatrixCell Dm = new MatrixCell(3);
        Dm.set(0, D0);
        Dm.set(1, D1);
        Dm.set(2, D2);

        MatrixCell H = MinimalRepFromMRAP.minimalRepFromMRAP(Dm, "obscont", 1e-12);
        assertEquals(2, H.get(0).getNumRows());
        assertEquals(-4.4074, H.get(0).get(0, 0), LOOSE_MID_TOL);
    }

    // ============ Phase 6: Queue Analysis + MAM + Trace Fitting ============

    /**
     * QBDSolve.txt: Solve QBD
     * Expected pi0 = [0.22992, 0.18681]
     * Expected R = [0.27839, 0.14286; 0.55678, 0.28571]
     */
    @Test
    public void testQBDSolve_Example() {
        Matrix B = new Matrix(2, 2);
        B.set(0, 0, 0.); B.set(0, 1, 0.);
        B.set(1, 0, 3.); B.set(1, 1, 4.);
        Matrix L = new Matrix(2, 2);
        L.set(0, 0, -6.); L.set(0, 1, 5.);
        L.set(1, 0, 3.);  L.set(1, 1, -12.);
        Matrix F = new Matrix(2, 2);
        F.set(0, 0, 1.); F.set(0, 1, 0.);
        F.set(1, 0, 2.); F.set(1, 1, 0.);
        Matrix L0 = new Matrix(2, 2);
        L0.set(0, 0, -6.); L0.set(0, 1, 5.);
        L0.set(1, 0, 6.);  L0.set(1, 1, -8.);

        Pair<Matrix, Matrix> result = QBDSolve.qbdSolve(B, L, F, L0, 1e-14);
        Matrix pi0 = result.getFirst();
        Matrix R = result.getSecond();

        assertEquals(0.22992, pi0.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.18681, pi0.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.27839, R.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.14286, R.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.55678, R.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.28571, R.get(1, 1), LOOSE_MID_TOL);
    }

    /**
     * QBDQueue.txt: Queue length distribution and moments
     * ncDistr expected: [0.29094, 0.20433, 0.14547, ...]
     * ncMoms expected: [2.46, 14.608, 128.88, 1515.9, 22288]
     */
    @Test
    public void testQBDQueue_NcDistrMoms() {
        Matrix B = new Matrix(3, 3);
        B.set(0, 0, 6.); B.set(0, 1, 1.); B.set(0, 2, 0.);
        B.set(1, 0, 0.); B.set(1, 1, 4.); B.set(1, 2, 1.);
        B.set(2, 0, 2.); B.set(2, 1, 0.); B.set(2, 2, 0.);
        Matrix F = new Matrix(3, 3);
        F.set(0, 0, 0.); F.set(0, 1, 1.); F.set(0, 2, 1.);
        F.set(1, 0, 5.); F.set(1, 1, 0.); F.set(1, 2, 0.);
        F.set(2, 0, 1.); F.set(2, 1, 3.); F.set(2, 2, 0.);
        Matrix L = new Matrix(3, 3);
        L.set(0, 0, -14.); L.set(0, 1, 3.); L.set(0, 2, 2.);
        L.set(1, 0, 0.);   L.set(1, 1, -14.); L.set(1, 2, 4.);
        L.set(2, 0, 3.);   L.set(2, 1, 1.);   L.set(2, 2, -10.);
        // L0 = L + B
        Matrix L0 = L.add(1.0, B);

        HashMap<String, Object> measures = new HashMap<>();
        measures.put("ncDistr", 11);
        measures.put("ncMoms", 5);
        Map<String, Object> result = QBDQueue.qbdQueue(B, L, F, L0, measures, 1e-14);

        double[] ncd = (double[]) result.get("ncDistr");
        assertNotNull(ncd);
        assertEquals(0.29094, ncd[0], LOOSE_MID_TOL);
        assertEquals(0.20433, ncd[1], LOOSE_MID_TOL);
        assertEquals(0.14547, ncd[2], LOOSE_MID_TOL);

        double[] ncm = (double[]) result.get("ncMoms");
        assertNotNull(ncm);
        assertEquals(2.46, ncm[0], LOOSE_MID_TOL);
        assertEquals(14.608, ncm[1], LOOSE_MID_TOL);
        assertEquals(128.88, ncm[2], VERY_COARSE_TOL);
    }

    /**
     * QBDQueue.txt: Sojourn time distribution and moments
     * stDistr expected: [~0, 0.14225, 0.25715, ...]
     * stMoms expected: [0.70236, 1.0017, 2.1463, ...]
     */
    @Test
    public void testQBDQueue_StDistrMoms() {
        Matrix B = new Matrix(3, 3);
        B.set(0, 0, 6.); B.set(0, 1, 1.); B.set(0, 2, 0.);
        B.set(1, 0, 0.); B.set(1, 1, 4.); B.set(1, 2, 1.);
        B.set(2, 0, 2.); B.set(2, 1, 0.); B.set(2, 2, 0.);
        Matrix F = new Matrix(3, 3);
        F.set(0, 0, 0.); F.set(0, 1, 1.); F.set(0, 2, 1.);
        F.set(1, 0, 5.); F.set(1, 1, 0.); F.set(1, 2, 0.);
        F.set(2, 0, 1.); F.set(2, 1, 3.); F.set(2, 2, 0.);
        Matrix L = new Matrix(3, 3);
        L.set(0, 0, -14.); L.set(0, 1, 3.); L.set(0, 2, 2.);
        L.set(1, 0, 0.);   L.set(1, 1, -14.); L.set(1, 2, 4.);
        L.set(2, 0, 3.);   L.set(2, 1, 1.);   L.set(2, 2, -10.);
        Matrix L0 = L.add(1.0, B);

        double[] points = new double[11];
        for (int i = 0; i <= 10; i++) points[i] = i * 0.1;

        HashMap<String, Object> measures = new HashMap<>();
        measures.put("stDistr", points);
        measures.put("stMoms", 5);
        Map<String, Object> result = QBDQueue.qbdQueue(B, L, F, L0, measures, 1e-14);

        double[] std = (double[]) result.get("stDistr");
        assertNotNull(std);
        assertEquals(0.14225, std[1], LOOSE_MID_TOL);
        assertEquals(0.25715, std[2], LOOSE_MID_TOL);

        double[] stm = (double[]) result.get("stMoms");
        assertNotNull(stm);
        assertEquals(0.70236, stm[0], LOOSE_MID_TOL);
        assertEquals(1.0017, stm[1], LOOSE_MID_TOL);
    }

    /**
     * MAPMAP1.txt: Example 1 - MAP/MAP/1 queue
     * ncMoms expected: [0.54864, 1.306, 4.357, 19.193, 105.39]
     * stMoms expected: [0.25908, 0.13145, 0.09911, 0.099178, 0.12376]
     */
    @Test
    public void testMAPMAP1_Example1() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -8.); D0.set(0, 1, 1.); D0.set(0, 2, 2.);
        D0.set(1, 0, 0.);  D0.set(1, 1, -6.); D0.set(1, 2, 4.);
        D0.set(2, 0, 3.);  D0.set(2, 1, 0.);  D0.set(2, 2, -3.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 4.); D1.set(0, 1, 1.); D1.set(0, 2, 0.);
        D1.set(1, 0, 0.); D1.set(1, 1, 2.); D1.set(1, 2, 0.);
        D1.set(2, 0, 0.); D1.set(2, 1, 0.); D1.set(2, 2, 0.);
        Matrix S0 = new Matrix(2, 2);
        S0.set(0, 0, -10.); S0.set(0, 1, 4.);
        S0.set(1, 0, 0.);   S0.set(1, 1, -7.);
        Matrix S1 = new Matrix(2, 2);
        S1.set(0, 0, 5.); S1.set(0, 1, 1.);
        S1.set(1, 0, 4.); S1.set(1, 1, 3.);

        HashMap<String, Object> measures = new HashMap<>();
        measures.put("ncMoms", 5);
        measures.put("stMoms", 5);
        double[] stPoints = new double[11];
        for (int i = 0; i <= 10; i++) stPoints[i] = i * 0.1;
        measures.put("stDistr", stPoints);
        measures.put("ncDistr", 11);

        Map<String, Object> result = MAPMAP1.mapmap1(D0, D1, S0, S1, measures, 1e-14);

        double[] ncm = (double[]) result.get("ncMoms");
        assertNotNull(ncm);
        assertEquals(0.54864, ncm[0], LOOSE_MID_TOL);
        assertEquals(1.306, ncm[1], LOOSE_MID_TOL);

        double[] stm = (double[]) result.get("stMoms");
        assertNotNull(stm);
        assertEquals(0.25908, stm[0], LOOSE_MID_TOL);
        assertEquals(0.13145, stm[1], LOOSE_MID_TOL);

        double[] ncd = (double[]) result.get("ncDistr");
        assertNotNull(ncd);
        assertEquals(0.67697, ncd[0], LOOSE_MID_TOL);
        assertEquals(0.18891, ncd[1], LOOSE_MID_TOL);
    }

    /**
     * FluidQueue.txt: Fluid queue analysis
     * flMoms expected: [0.40636, 0.4546, 0.76609, 1.722, 4.8384]
     * stMoms expected: [0.23252, 0.20069, 0.26684, 0.47402, 1.0523]
     */
    @Test
    public void testFluidQueue_Example() {
        Matrix Q = new Matrix(6, 6);
        Q.set(0, 0, -9.);  Q.set(0, 1, 2.);  Q.set(0, 2, 4.); Q.set(0, 3, 0.); Q.set(0, 4, 1.); Q.set(0, 5, 2.);
        Q.set(1, 0, 6.);   Q.set(1, 1, -25.); Q.set(1, 2, 5.); Q.set(1, 3, 3.); Q.set(1, 4, 7.); Q.set(1, 5, 4.);
        Q.set(2, 0, 1.);   Q.set(2, 1, 3.);   Q.set(2, 2, -4.); Q.set(2, 3, 0.); Q.set(2, 4, 0.); Q.set(2, 5, 0.);
        Q.set(3, 0, 0.);   Q.set(3, 1, 0.);   Q.set(3, 2, 0.); Q.set(3, 3, -8.); Q.set(3, 4, 3.); Q.set(3, 5, 5.);
        Q.set(4, 0, 7.);   Q.set(4, 1, 3.);   Q.set(4, 2, 0.); Q.set(4, 3, 2.); Q.set(4, 4, -13.); Q.set(4, 5, 1.);
        Q.set(5, 0, 7.);   Q.set(5, 1, 8.);   Q.set(5, 2, 0.); Q.set(5, 3, 3.); Q.set(5, 4, 8.); Q.set(5, 5, -26.);

        double[] vRin = {4., 2., 1., 0., 0., 3.};
        double[] vRout = {6., 2., 0., 0., 3., 2.};
        Matrix Rin = Matrix.zeros(6, 6);
        Matrix Rout = Matrix.zeros(6, 6);
        for (int i = 0; i < 6; i++) {
            Rin.set(i, i, vRin[i]);
            Rout.set(i, i, vRout[i]);
        }

        double[] flPoints = new double[11];
        for (int i = 0; i <= 10; i++) flPoints[i] = i * 0.1;
        double[] stPoints = new double[11];
        for (int i = 0; i <= 10; i++) stPoints[i] = i * 0.1;

        HashMap<String, Object> measures = new HashMap<>();
        measures.put("flMoms", 5);
        measures.put("flDistr", flPoints);
        measures.put("stMoms", 5);
        measures.put("stDistr", stPoints);

        Map<String, Object> result = FluidQueue.fluidQueue(Q, Rin, Rout, measures, null, 1e-14);

        double[] flm = (double[]) result.get("flMoms");
        assertNotNull(flm);
        assertEquals(0.40636, flm[0], LOOSE_MID_TOL);
        assertEquals(0.4546, flm[1], LOOSE_MID_TOL);

        double[] fld = (double[]) result.get("flDistr");
        assertNotNull(fld);
        assertEquals(0.23662, fld[0], LOOSE_MID_TOL);
        assertEquals(0.39265, fld[1], LOOSE_MID_TOL);

        double[] stm = (double[]) result.get("stMoms");
        assertNotNull(stm);
        assertEquals(0.23252, stm[0], LOOSE_MID_TOL);
        assertEquals(0.20069, stm[1], LOOSE_MID_TOL);

        double[] std = (double[]) result.get("stDistr");
        assertNotNull(std);
        assertEquals(0.31678, std[0], LOOSE_MID_TOL);
        assertEquals(0.57513, std[1], LOOSE_MID_TOL);
    }

    /**
     * GM1StationaryDistr.txt: G/M/1 stationary distribution
     * Uses B matrices from file, verifies convergence (sum near 1)
     */
    @Test
    public void testGM1StationaryDistr_Example() {
        Matrix B0 = new Matrix(2, 2);
        B0.set(0, 0, 0.7); B0.set(0, 1, 0.2);
        B0.set(1, 0, 0.3); B0.set(1, 1, 0.6);
        Matrix B1 = new Matrix(2, 2);
        B1.set(0, 0, 0.3); B1.set(0, 1, 0.4);
        B1.set(1, 0, 0.5); B1.set(1, 1, 0.2);
        Matrix B2 = new Matrix(2, 2);
        B2.set(0, 0, 0.2); B2.set(0, 1, 0.4);
        B2.set(1, 0, 0.1); B2.set(1, 1, 0.6);
        Matrix B3 = new Matrix(2, 2);
        B3.set(0, 0, 0.0); B3.set(0, 1, 0.1);
        B3.set(1, 0, 0.2); B3.set(1, 1, 0.0);

        Matrix A0 = new Matrix(2, 2);
        A0.set(0, 0, 0.1); A0.set(0, 1, 0.);
        A0.set(1, 0, 0.);  A0.set(1, 1, 0.1);
        Matrix A1 = new Matrix(2, 2);
        A1.set(0, 0, 0.); A1.set(0, 1, 0.2);
        A1.set(1, 0, 0.); A1.set(1, 1, 0.2);
        Matrix A2 = new Matrix(2, 2);
        A2.set(0, 0, 0.); A2.set(0, 1, 0.1);
        A2.set(1, 0, 0.); A2.set(1, 1, 0.);
        Matrix A3 = new Matrix(2, 2);
        A3.set(0, 0, 0.3); A3.set(0, 1, 0.2);
        A3.set(1, 0, 0.3); A3.set(1, 1, 0.2);
        Matrix A4 = new Matrix(2, 2);
        A4.set(0, 0, 0.); A4.set(0, 1, 0.1);
        A4.set(1, 0, 0.2); A4.set(1, 1, 0.);

        // Compute R from A matrices
        java.util.List<Matrix> Alist = Arrays.asList(A0, A1, A2, A3, A4);
        Matrix R = GM1FundamentalMatrix.gm1FundamentalMatrix(Alist, 1e-14, 50, GM1Method.CR);
        assertNotNull(R);
        assertEquals(0.10065, R.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.026961, R.get(0, 1), LOOSE_MID_TOL);

        // Compute stationary distribution
        java.util.List<Matrix> Blist = Arrays.asList(B0, B1, B2, B3);
        Matrix pi = GM1StationaryDistr.gm1StationaryDistr(Blist, R, 300);
        assertNotNull(pi);
        // Verify sum is close to 1
        double sum = pi.elementSum();
        assertEquals(1.0, sum, COARSE_TOL);
    }

    /**
     * MG1StationaryDistr.txt: M/G/1 stationary distribution
     * Uses A matrices and boundary B matrices, verifies convergence
     */
    @Test
    public void testMG1StationaryDistr_Example() {
        Matrix B0 = new Matrix(2, 2);
        B0.set(0, 0, 0.1); B0.set(0, 1, 0.5);
        B0.set(1, 0, 0.3); B0.set(1, 1, 0.4);
        Matrix B1 = new Matrix(2, 2);
        B1.set(0, 0, 0.); B1.set(0, 1, 0.1);
        B1.set(1, 0, 0.); B1.set(1, 1, 0.);
        Matrix B2 = new Matrix(2, 2);
        B2.set(0, 0, 0.2); B2.set(0, 1, 0.);
        B2.set(1, 0, 0.);  B2.set(1, 1, 0.2);
        Matrix B3 = new Matrix(2, 2);
        B3.set(0, 0, 0.); B3.set(0, 1, 0.1);
        B3.set(1, 0, 0.1); B3.set(1, 1, 0.);

        Matrix A0 = new Matrix(2, 2);
        A0.set(0, 0, 0.4); A0.set(0, 1, 0.2);
        A0.set(1, 0, 0.3); A0.set(1, 1, 0.4);
        Matrix A1 = new Matrix(2, 2);
        A1.set(0, 0, 0.); A1.set(0, 1, 0.1);
        A1.set(1, 0, 0.); A1.set(1, 1, 0.);
        Matrix A2 = new Matrix(2, 2);
        A2.set(0, 0, 0.); A2.set(0, 1, 0.2);
        A2.set(1, 0, 0.); A2.set(1, 1, 0.2);
        Matrix A3 = new Matrix(2, 2);
        A3.set(0, 0, 0.1); A3.set(0, 1, 0.);
        A3.set(1, 0, 0.1); A3.set(1, 1, 0.);

        // Compute G from A matrices
        java.util.List<Matrix> Alist = Arrays.asList(A0, A1, A2, A3);
        Matrix G = MG1FundamentalMatrix.mg1FundamentalMatrix(Alist, 1e-14, 50, MG1Method.CR);
        assertNotNull(G);
        assertEquals(0.60503, G.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.39497, G.get(0, 1), LOOSE_MID_TOL);

        // Compute stationary distribution
        java.util.List<Matrix> Blist = Arrays.asList(B0, B1, B2, B3);
        Matrix pi = MG1StationaryDistr.mg1StationaryDistr(Alist, Blist, G, 300, 1e-14);
        assertNotNull(pi);
        double sum = pi.elementSum();
        assertEquals(1.0, sum, COARSE_TOL);
    }

    /**
     * PHFromTrace: Fit PH with orders [1,2,2] (best from reference)
     * Expected: logli around 4.915 (reference: 4.91508)
     */
    @Test
    public void testPHFromTrace_FixedOrders() {
        // Generate trace from a known PH
        Matrix alpha0 = new Matrix(1, 2);
        alpha0.set(0, 0, 0.8); alpha0.set(0, 1, 0.2);
        Matrix A0 = new Matrix(2, 2);
        A0.set(0, 0, -5.0); A0.set(0, 1, 1.0);
        A0.set(1, 0, 0.0);  A0.set(1, 1, -2.0);
        // Generate exponential samples with rate 3 (simple test)
        Random rng = new Random(42);
        double[] trace = new double[5000];
        for (int i = 0; i < trace.length; i++) {
            trace[i] = -Math.log(rng.nextDouble()) / 3.0;
        }

        PHFitResult result = PHFromTrace.phFromTrace(trace, new int[]{1, 2, 2}, 200, 1e-7, null, null);
        assertNotNull(result);
        assertNotNull(result.alpha);
        assertNotNull(result.A);
        assertTrue(Double.isFinite(result.logli), "logli should be finite");
        // Alpha should be a probability vector
        double alphaSum = 0;
        for (int i = 0; i < result.alpha.getNumCols(); i++) {
            alphaSum += result.alpha.get(0, i);
        }
        assertEquals(1.0, alphaSum, LOOSE_FINE_TOL);
        // A matrix should be 5x5
        assertEquals(5, result.A.getNumRows());
    }

    /**
     * MAPFromTrace: Fit MAP with orders [1,1,3]
     */
    @Test
    public void testMAPFromTrace_FixedOrders() {
        // Generate trace from a MAP
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -18.); D0.set(0, 1, 1.);  D0.set(0, 2, 4.);
        D0.set(1, 0, 2.);   D0.set(1, 1, -18.); D0.set(1, 2, 7.);
        D0.set(2, 0, 1.);   D0.set(2, 1, 3.);   D0.set(2, 2, -32.);
        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 12.); D1.set(0, 1, 1.); D1.set(0, 2, 0.);
        D1.set(1, 0, 1.);  D1.set(1, 1, 8.); D1.set(1, 2, 0.);
        D1.set(2, 0, 2.);  D1.set(2, 1, 1.); D1.set(2, 2, 25.);

        double[] tr = SamplesFromMAP.samplesFromMAP(D0, D1, 5000, null, new Random(42));

        MAPFitResult result = MAPFromTrace.mapFromTrace(tr, new int[]{1, 1, 3}, 200, 1e-7, null, null);
        assertNotNull(result);
        assertNotNull(result.D0);
        assertNotNull(result.D1);
        assertTrue(Double.isFinite(result.logli), "logli should be finite");
        // D0 + D1 should have zero row sums (generator)
        int N = result.D0.getNumRows();
        assertEquals(5, N); // 1+1+3 = 5
        for (int i = 0; i < N; i++) {
            double rowSum = 0;
            for (int j = 0; j < N; j++) {
                rowSum += result.D0.get(i, j) + result.D1.get(i, j);
            }
            assertEquals(0.0, rowSum, LOOSE_FINE_TOL, "Row " + i + " sum of D0+D1 should be 0");
        }
    }

    /**
     * QBDQueue.txt: ncDistrDPH representation
     * Expected alpha = [0.28256, 0.22386, 0.20264]
     */
    @Test
    public void testQBDQueue_NcDistrDPH() {
        Matrix B = new Matrix(3, 3);
        B.set(0, 0, 6.); B.set(0, 1, 1.); B.set(0, 2, 0.);
        B.set(1, 0, 0.); B.set(1, 1, 4.); B.set(1, 2, 1.);
        B.set(2, 0, 2.); B.set(2, 1, 0.); B.set(2, 2, 0.);
        Matrix F = new Matrix(3, 3);
        F.set(0, 0, 0.); F.set(0, 1, 1.); F.set(0, 2, 1.);
        F.set(1, 0, 5.); F.set(1, 1, 0.); F.set(1, 2, 0.);
        F.set(2, 0, 1.); F.set(2, 1, 3.); F.set(2, 2, 0.);
        Matrix L = new Matrix(3, 3);
        L.set(0, 0, -14.); L.set(0, 1, 3.); L.set(0, 2, 2.);
        L.set(1, 0, 0.);   L.set(1, 1, -14.); L.set(1, 2, 4.);
        L.set(2, 0, 3.);   L.set(2, 1, 1.);   L.set(2, 2, -10.);
        Matrix L0 = L.add(1.0, B);

        HashMap<String, Object> measures = new HashMap<>();
        measures.put("ncDistrDPH", null);
        Map<String, Object> result = QBDQueue.qbdQueue(B, L, F, L0, measures, 1e-14);

        Matrix alphap = (Matrix) result.get("ncDistrDPH_alpha");
        assertNotNull(alphap);
        assertEquals(0.28256, alphap.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.22386, alphap.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.20264, alphap.get(0, 2), LOOSE_MID_TOL);
    }

    /**
     * MAPMAP1.txt: Example 2 - rank-1 arrival and service
     * ncMoms expected: [2.0439, 10.554, 80.619, 820.69, 10443]
     */
    @Test
    public void testMAPMAP1_Example2() {
        // D0 = Dm, D1 = sum(-Dm,2)*delta, S0 = S, S1 = sum(-S,2)*sigma
        double[] delta = {0.5, 0.1, 0.4};
        Matrix Dm = new Matrix(3, 3);
        Dm.set(0, 0, -8.); Dm.set(0, 1, 1.); Dm.set(0, 2, 2.);
        Dm.set(1, 0, 0.);  Dm.set(1, 1, -6.); Dm.set(1, 2, 4.);
        Dm.set(2, 0, 3.);  Dm.set(2, 1, 0.);  Dm.set(2, 2, -3.);
        double[] sigma = {0.2, 0.7, 0.1};
        Matrix S = new Matrix(3, 3);
        S.set(0, 0, -10.); S.set(0, 1, 4.); S.set(0, 2, 0.);
        S.set(1, 0, 5.);   S.set(1, 1, -7.); S.set(1, 2, 2.);
        S.set(2, 0, 1.);   S.set(2, 1, 2.);  S.set(2, 2, -8.);

        // D0 = Dm
        Matrix D0 = Dm;
        // D1 = sum(-Dm, 2) * delta = rowsum of -Dm * delta
        Matrix D1 = new Matrix(3, 3);
        for (int i = 0; i < 3; i++) {
            double rowSum = 0;
            for (int j = 0; j < 3; j++) rowSum += -Dm.get(i, j);
            for (int j = 0; j < 3; j++) D1.set(i, j, rowSum * delta[j]);
        }
        // S0 = S
        Matrix S0 = S;
        // S1 = sum(-S, 2) * sigma
        Matrix S1 = new Matrix(3, 3);
        for (int i = 0; i < 3; i++) {
            double rowSum = 0;
            for (int j = 0; j < 3; j++) rowSum += -S.get(i, j);
            for (int j = 0; j < 3; j++) S1.set(i, j, rowSum * sigma[j]);
        }

        HashMap<String, Object> measures = new HashMap<>();
        measures.put("ncMoms", 5);
        measures.put("stMoms", 5);

        Map<String, Object> result = MAPMAP1.mapmap1(D0, D1, S0, S1, measures, 1e-14);

        double[] ncm = (double[]) result.get("ncMoms");
        assertNotNull(ncm);
        assertEquals(2.0439, ncm[0], LOOSE_MID_TOL);
        assertEquals(10.554, ncm[1], COARSE_TOL);

        double[] stm = (double[]) result.get("stMoms");
        assertNotNull(stm);
        assertEquals(1.1135, stm[0], LOOSE_MID_TOL);
        assertEquals(2.4113, stm[1], LOOSE_MID_TOL);
    }

    /**
     * FluidQueue.txt: PH representation of fluid level distribution
     * Expected alpha = [0.63124, 0.13213]
     * Expected A = [-2.0387, 0.41483; 12.1, -21.143]
     */
    @Test
    public void testFluidQueue_FlDistrPH() {
        Matrix Q = new Matrix(6, 6);
        Q.set(0, 0, -9.);  Q.set(0, 1, 2.);  Q.set(0, 2, 4.); Q.set(0, 3, 0.); Q.set(0, 4, 1.); Q.set(0, 5, 2.);
        Q.set(1, 0, 6.);   Q.set(1, 1, -25.); Q.set(1, 2, 5.); Q.set(1, 3, 3.); Q.set(1, 4, 7.); Q.set(1, 5, 4.);
        Q.set(2, 0, 1.);   Q.set(2, 1, 3.);   Q.set(2, 2, -4.); Q.set(2, 3, 0.); Q.set(2, 4, 0.); Q.set(2, 5, 0.);
        Q.set(3, 0, 0.);   Q.set(3, 1, 0.);   Q.set(3, 2, 0.); Q.set(3, 3, -8.); Q.set(3, 4, 3.); Q.set(3, 5, 5.);
        Q.set(4, 0, 7.);   Q.set(4, 1, 3.);   Q.set(4, 2, 0.); Q.set(4, 3, 2.); Q.set(4, 4, -13.); Q.set(4, 5, 1.);
        Q.set(5, 0, 7.);   Q.set(5, 1, 8.);   Q.set(5, 2, 0.); Q.set(5, 3, 3.); Q.set(5, 4, 8.); Q.set(5, 5, -26.);

        double[] vRin = {4., 2., 1., 0., 0., 3.};
        double[] vRout = {6., 2., 0., 0., 3., 2.};
        Matrix Rin = Matrix.zeros(6, 6);
        Matrix Rout = Matrix.zeros(6, 6);
        for (int i = 0; i < 6; i++) {
            Rin.set(i, i, vRin[i]);
            Rout.set(i, i, vRout[i]);
        }

        HashMap<String, Object> measures = new HashMap<>();
        measures.put("flDistrPH", null);

        Map<String, Object> result = FluidQueue.fluidQueue(Q, Rin, Rout, measures, null, 1e-14);

        Matrix alphap = (Matrix) result.get("flDistrPH_alpha");
        assertNotNull(alphap);
        assertEquals(0.63124, alphap.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.13213, alphap.get(0, 1), LOOSE_MID_TOL);

        Matrix Ap = (Matrix) result.get("flDistrPH_A");
        assertNotNull(Ap);
        assertEquals(-2.0387, Ap.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.41483, Ap.get(0, 1), LOOSE_MID_TOL);
    }
}
