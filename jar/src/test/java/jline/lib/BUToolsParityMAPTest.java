package jline.lib;

import jline.lib.butools.map.MarginalDistributionFromMAP;
import jline.lib.butools.map.MarginalMomentsFromMAP;
import jline.lib.butools.map.LagCorrelationsFromMAP;
import jline.lib.butools.map.MAP2CorrelationBounds;
import jline.lib.butools.map.LagkJointMomentsFromMAP;
import jline.lib.butools.map.LagkJointMomentsFromRAP;
import jline.lib.butools.map.LagkJointMomentsFromMMAP;
import jline.lib.butools.map.LagkJointMomentsFromMRAP;
import jline.lib.butools.map.MarginalDistributionFromMRAP;
import jline.lib.butools.map.MarginalDistributionFromMMAP;
import jline.lib.butools.map.MarginalMomentsFromMRAP;
import jline.lib.butools.map.MarginalMomentsFromMMAP;
import jline.lib.butools.map.CheckMMAPRepresentation;
import jline.lib.butools.map.CheckMRAPRepresentation;
import jline.lib.butools.map.MAPFromRAP;
import jline.lib.butools.map.MMAPFromMRAP;
import jline.lib.butools.map.RandomMMAP;
import jline.lib.butools.map.SamplesFromMMAP;
import jline.lib.butools.map.MRAPFromMoments;
import jline.lib.butools.map.RAPFromMomentsAndCorrelations;
import jline.lib.butools.map.MAPFromFewMomentsAndCorrelations;
import jline.lib.butools.dmap.MarginalMomentsFromDMAP;
import jline.lib.butools.dmap.MarginalDistributionFromDMAP;
import jline.lib.butools.dmap.LagCorrelationsFromDMAP;
import jline.lib.butools.dmap.LagkJointMomentsFromDMAP;
import jline.lib.butools.dmap.CheckDMAPRepresentation;
import jline.lib.butools.dmap.DMAP2FromMoments;
import jline.lib.butools.dmap.CanonicalFromDMAP2;
import jline.lib.butools.dmap.DMAPFromDRAP;
import jline.lib.butools.dmap.MarginalMomentsFromDRAP;
import jline.lib.butools.dmap.MarginalDistributionFromDRAP;
import jline.lib.butools.dmap.LagCorrelationsFromDRAP;
import jline.lib.butools.dmap.LagkJointMomentsFromDRAP;
import jline.lib.butools.dmap.CheckDRAPRepresentation;
import jline.lib.butools.dmap.CheckDMMAPRepresentation;
import jline.lib.butools.dmap.DMMAPFromDMRAP;
import jline.lib.butools.dmap.DMRAPFromMoments;
import jline.lib.butools.dmap.DRAPFromMoments;
import jline.lib.butools.dmap.LagkJointMomentsFromDMMAP;
import jline.lib.butools.dmap.LagkJointMomentsFromDMRAP;
import jline.lib.butools.dmap.MarginalDistributionFromDMRAP;
import jline.lib.butools.dmap.MarginalMomentsFromDMRAP;
import jline.lib.butools.dmap.MarginalMomentsFromDMMAP;
import jline.lib.butools.ph.PHRepresentation;
import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.lib.butools.map.RandomMAP;
import jline.lib.butools.map.MAPRepresentation;
import jline.lib.butools.map.CheckMAPRepresentation;
import jline.lib.butools.map.SamplesFromMAP;
import jline.lib.butools.dmap.RandomDMAP;
import jline.lib.butools.dmap.RandomDMMAP;
import jline.lib.butools.dmap.SamplesFromDMAP;
import jline.lib.butools.dmap.SamplesFromDMMAP;
import jline.lib.butools.trace.MarginalMomentsFromTrace;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.Pair;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.*;

import java.util.Random;

/**
 * BUTools Parity Tests for MAP/RAP/DMAP/DRAP functions.
 * Each test uses EXACT data from butools.tmp example files.
 */
public class BUToolsParityMAPTest {


    // ============ Continuous MAP Tests ============

    /**
     * MarginalDistributionFromMAP.txt: Compute marginal distribution from MAP
     * Input: D0 4x4, D1 4x4
     * Expected: a=[0.14438, 0.23571, 0.33794, 0.28196], A=D0
     */
    @Test
    public void testMarginalDistributionFromMAP() {
        Matrix D0 = new Matrix(4, 4);
        D0.set(0, 0, -0.17); D0.set(0, 1, 0.0);   D0.set(0, 2, 0.0);  D0.set(0, 3, 0.07);
        D0.set(1, 0, 0.01);  D0.set(1, 1, -0.78);  D0.set(1, 2, 0.03); D0.set(1, 3, 0.08);
        D0.set(2, 0, 0.22);  D0.set(2, 1, 0.17);   D0.set(2, 2, -1.1); D0.set(2, 3, 0.02);
        D0.set(3, 0, 0.04);  D0.set(3, 1, 0.12);   D0.set(3, 2, 0.0);  D0.set(3, 3, -0.42);

        Matrix D1 = new Matrix(4, 4);
        D1.set(0, 0, 0.0);  D1.set(0, 1, 0.06); D1.set(0, 2, 0.0);  D1.set(0, 3, 0.04);
        D1.set(1, 0, 0.04); D1.set(1, 1, 0.19); D1.set(1, 2, 0.21); D1.set(1, 3, 0.22);
        D1.set(2, 0, 0.22); D1.set(2, 1, 0.13); D1.set(2, 2, 0.15); D1.set(2, 3, 0.19);
        D1.set(3, 0, 0.05); D1.set(3, 1, 0.0);  D1.set(3, 2, 0.17); D1.set(3, 3, 0.04);

        PHRepresentation result = MarginalDistributionFromMAP.marginalDistributionFromMAP(D0, D1);
        Matrix a = result.getAlpha();
        Matrix A = result.getA();

        // Verify alpha
        double[] expectedAlpha = {0.14438, 0.23571, 0.33794, 0.28196};
        for (int i = 0; i < 4; i++) {
            assertEquals(expectedAlpha[i], a.get(0, i), LOOSE_MID_TOL, "alpha[" + i + "]");
        }

        // Verify A == D0
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(D0.get(i, j), A.get(i, j), LOOSE_MID_TOL, "A[" + i + "," + j + "]");
            }
        }
    }

    /**
     * MarginalMomentsFromMAP.txt: Compute marginal moments from MAP
     * Input: same D0/D1 as MarginalDistributionFromMAP
     * Expected: [3.4433, 34.03, 592.08, 14548, 4.559e+05, 1.727e+07, 7.6526e+08]
     */
    @Test
    public void testMarginalMomentsFromMAP() {
        Matrix D0 = new Matrix(4, 4);
        D0.set(0, 0, -0.17); D0.set(0, 1, 0.0);   D0.set(0, 2, 0.0);  D0.set(0, 3, 0.07);
        D0.set(1, 0, 0.01);  D0.set(1, 1, -0.78);  D0.set(1, 2, 0.03); D0.set(1, 3, 0.08);
        D0.set(2, 0, 0.22);  D0.set(2, 1, 0.17);   D0.set(2, 2, -1.1); D0.set(2, 3, 0.02);
        D0.set(3, 0, 0.04);  D0.set(3, 1, 0.12);   D0.set(3, 2, 0.0);  D0.set(3, 3, -0.42);

        Matrix D1 = new Matrix(4, 4);
        D1.set(0, 0, 0.0);  D1.set(0, 1, 0.06); D1.set(0, 2, 0.0);  D1.set(0, 3, 0.04);
        D1.set(1, 0, 0.04); D1.set(1, 1, 0.19); D1.set(1, 2, 0.21); D1.set(1, 3, 0.22);
        D1.set(2, 0, 0.22); D1.set(2, 1, 0.13); D1.set(2, 2, 0.15); D1.set(2, 3, 0.19);
        D1.set(3, 0, 0.05); D1.set(3, 1, 0.0);  D1.set(3, 2, 0.17); D1.set(3, 3, 0.04);

        double[] moms = MarginalMomentsFromMAP.marginalMomentsFromMAP(D0, D1, 7);

        double[] expected = {3.4433, 34.03, 592.08, 14548.0, 455900.0, 17270000.0, 765260000.0};
        assertEquals(expected.length, moms.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], moms[i], Math.abs(expected[i]) * 1e-3, "moms[" + i + "]");
        }
    }

    /**
     * LagCorrelationsFromMAP.txt: Compute lag correlations from MAP
     * Input: D0 4x4, D1 4x4
     * Expected: [0.00012012, 0.00086176, -0.00022001]
     */
    @Test
    public void testLagCorrelationsFromMAP() {
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

        double[] corr = LagCorrelationsFromMAP.lagCorrelationsFromMAP(D0, D1, 3);

        assertEquals(3, corr.length);
        assertEquals(0.00012012, corr[0], LOOSE_MID_TOL);
        assertEquals(0.00086176, corr[1], LOOSE_MID_TOL);
        assertEquals(-0.00022001, corr[2], LOOSE_MID_TOL);
    }

    /**
     * MAP2CorrelationBounds.txt: Compute MAP2 correlation bounds
     * Input: moms from D0=[-14,1;1,-25], D1=[6,7;3,21]
     * Expected: lb=-0.030588, ub=0.074506
     */
    @Test
    public void testMAP2CorrelationBounds() {
        // First compute moms from the MAP
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -14.0); D0.set(0, 1, 1.0);
        D0.set(1, 0, 1.0);   D0.set(1, 1, -25.0);

        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, 6.0); D1.set(0, 1, 7.0);
        D1.set(1, 0, 3.0); D1.set(1, 1, 21.0);

        double[] moms = MarginalMomentsFromMAP.marginalMomentsFromMAP(D0, D1, 3);
        assertEquals(0.04918, moms[0], LOOSE_MID_TOL);
        assertEquals(0.0052609, moms[1], LOOSE_MID_TOL);
        assertEquals(0.00091819, moms[2], LOOSE_MID_TOL);

        jline.util.Pair<Double, Double> bounds = MAP2CorrelationBounds.map2CorrelationBounds(moms);
        assertEquals(-0.030588, bounds.getFirst(), LOOSE_MID_TOL);
        assertEquals(0.074506, bounds.getSecond(), LOOSE_MID_TOL);
    }

    // ============ Continuous RAP Tests ============

    /**
     * MarginalDistributionFromRAP.txt: Compute marginal distribution from RAP
     * Input: H0 3x3, H1 3x3
     * Expected: a=[0.44444, 0.44444, 0.11111], A=H0
     */
    @Test
    public void testMarginalDistributionFromRAP() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, -2.0); H0.set(0, 1, 0.0);  H0.set(0, 2, 0.0);
        H0.set(1, 0, 0.0);  H0.set(1, 1, -3.0);  H0.set(1, 2, 1.0);
        H0.set(2, 0, 0.0);  H0.set(2, 1, -1.0);  H0.set(2, 2, -2.0);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 1.8); H1.set(0, 1, 0.2); H1.set(0, 2, 0.0);
        H1.set(1, 0, 0.2); H1.set(1, 1, 1.8); H1.set(1, 2, 0.0);
        H1.set(2, 0, 0.2); H1.set(2, 1, 1.8); H1.set(2, 2, 1.0);

        PHRepresentation result = MarginalDistributionFromMAP.marginalDistributionFromRAP(H0, H1);
        Matrix a = result.getAlpha();
        Matrix A = result.getA();

        double[] expectedAlpha = {0.44444, 0.44444, 0.11111};
        for (int i = 0; i < 3; i++) {
            assertEquals(expectedAlpha[i], a.get(0, i), LOOSE_MID_TOL, "alpha[" + i + "]");
        }

        // A should equal H0
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(H0.get(i, j), A.get(i, j), LOOSE_MID_TOL, "A[" + i + "," + j + "]");
            }
        }
    }

    /**
     * MarginalMomentsFromRAP.txt: Compute marginal moments from RAP
     * Input: H0 3x3, H1 3x3
     * Expected: [0.44444, 0.38095, 0.48299, 0.82216, 1.7944]
     */
    @Test
    public void testMarginalMomentsFromRAP() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, -2.0); H0.set(0, 1, 0.0);  H0.set(0, 2, 0.0);
        H0.set(1, 0, 0.0);  H0.set(1, 1, -3.0);  H0.set(1, 2, 1.0);
        H0.set(2, 0, 0.0);  H0.set(2, 1, -1.0);  H0.set(2, 2, -2.0);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 1.8); H1.set(0, 1, 0.2); H1.set(0, 2, 0.0);
        H1.set(1, 0, 0.2); H1.set(1, 1, 1.8); H1.set(1, 2, 0.0);
        H1.set(2, 0, 0.2); H1.set(2, 1, 1.8); H1.set(2, 2, 1.0);

        double[] moms = MarginalMomentsFromMAP.marginalMomentsFromRAP(H0, H1, 5);

        double[] expected = {0.44444, 0.38095, 0.48299, 0.82216, 1.7944};
        assertEquals(expected.length, moms.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], moms[i], LOOSE_MID_TOL, "moms[" + i + "]");
        }
    }

    /**
     * LagCorrelationsFromRAP.txt: Compute lag correlations from RAP
     * Input: H0 3x3, H1 3x3
     * Expected: [-0.0038462, 0.0045604, 0.0058956]
     */
    @Test
    public void testLagCorrelationsFromRAP() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, -2.0); H0.set(0, 1, 0.0);  H0.set(0, 2, 0.0);
        H0.set(1, 0, 0.0);  H0.set(1, 1, -3.0);  H0.set(1, 2, 1.0);
        H0.set(2, 0, 0.0);  H0.set(2, 1, -1.0);  H0.set(2, 2, -2.0);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 1.8); H1.set(0, 1, 0.2); H1.set(0, 2, 0.0);
        H1.set(1, 0, 0.2); H1.set(1, 1, 1.8); H1.set(1, 2, 0.0);
        H1.set(2, 0, 0.2); H1.set(2, 1, 1.8); H1.set(2, 2, 1.0);

        double[] corr = LagCorrelationsFromMAP.lagCorrelationsFromRAP(H0, H1, 3);

        assertEquals(3, corr.length);
        assertEquals(-0.0038462, corr[0], LOOSE_MID_TOL);
        assertEquals(0.0045604, corr[1], LOOSE_MID_TOL);
        assertEquals(0.0058956, corr[2], LOOSE_MID_TOL);
    }

    // ============ Discrete DMAP Tests ============

    /** Helper: create the standard 4x4 DMAP D0 used in multiple DMAP tests */
    private Matrix createDmapD0() {
        Matrix D0 = new Matrix(4, 4);
        D0.set(0, 0, 0.0);  D0.set(0, 1, 0.02); D0.set(0, 2, 0.0);  D0.set(0, 3, 0.0);
        D0.set(1, 0, 0.0);  D0.set(1, 1, 0.17); D0.set(1, 2, 0.2);  D0.set(1, 3, 0.14);
        D0.set(2, 0, 0.16); D0.set(2, 1, 0.17); D0.set(2, 2, 0.02); D0.set(2, 3, 0.18);
        D0.set(3, 0, 0.0);  D0.set(3, 1, 0.0);  D0.set(3, 2, 0.0);  D0.set(3, 3, 0.12);
        return D0;
    }

    /** Helper: create the standard 4x4 DMAP D1 used in multiple DMAP tests */
    private Matrix createDmapD1() {
        Matrix D1 = new Matrix(4, 4);
        D1.set(0, 0, 0.0);  D1.set(0, 1, 0.88); D1.set(0, 2, 0.1);  D1.set(0, 3, 0.0);
        D1.set(1, 0, 0.18); D1.set(1, 1, 0.07); D1.set(1, 2, 0.14); D1.set(1, 3, 0.1);
        D1.set(2, 0, 0.13); D1.set(2, 1, 0.15); D1.set(2, 2, 0.15); D1.set(2, 3, 0.04);
        D1.set(3, 0, 0.31); D1.set(3, 1, 0.18); D1.set(3, 2, 0.12); D1.set(3, 3, 0.27);
        return D1;
    }

    /**
     * MarginalMomentsFromDMAP.txt: Compute marginal moments from DMAP
     * Expected: [1.4955, 2.9542, 7.8852, 27.282, 116.17, 587.04, 3437]
     */
    @Test
    public void testMarginalMomentsFromDMAP() {
        Matrix D0 = createDmapD0();
        Matrix D1 = createDmapD1();

        double[] moms = MarginalMomentsFromDMAP.marginalMomentsFromDMAP(D0, D1, 7);

        double[] expected = {1.4955, 2.9542, 7.8852, 27.282, 116.17, 587.04, 3437.0};
        assertEquals(expected.length, moms.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], moms[i], Math.max(LOOSE_MID_TOL, Math.abs(expected[i]) * 1e-3), "moms[" + i + "]");
        }
    }

    /**
     * MarginalDistributionFromDMAP.txt: Compute marginal distribution from DMAP
     * Expected: a=[0.24388, 0.40412, 0.1941, 0.1579], A=D0
     */
    @Test
    public void testMarginalDistributionFromDMAP() {
        Matrix D0 = createDmapD0();
        Matrix D1 = createDmapD1();

        MGRepresentation result = MarginalDistributionFromDMAP.marginalDistributionFromDMAP(D0, D1, 1e-14);
        Matrix a = result.getAlpha();
        Matrix A = result.getA();

        double[] expectedAlpha = {0.24388, 0.40412, 0.1941, 0.1579};
        for (int i = 0; i < 4; i++) {
            assertEquals(expectedAlpha[i], a.get(0, i), LOOSE_MID_TOL, "alpha[" + i + "]");
        }

        // A should equal D0
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(D0.get(i, j), A.get(i, j), LOOSE_MID_TOL, "A[" + i + "," + j + "]");
            }
        }
    }

    /**
     * LagCorrelationsFromDMAP.txt: Compute lag correlations from DMAP
     * Expected: [-0.045859, 0.010753, -0.0047996]
     */
    @Test
    public void testLagCorrelationsFromDMAP() {
        Matrix D0 = createDmapD0();
        Matrix D1 = createDmapD1();

        double[] corr = LagCorrelationsFromDMAP.lagCorrelationsFromDMAP(D0, D1, 3, 1e-14);

        assertEquals(3, corr.length);
        assertEquals(-0.045859, corr[0], LOOSE_MID_TOL);
        assertEquals(0.010753, corr[1], LOOSE_MID_TOL);
        assertEquals(-0.0047996, corr[2], LOOSE_MID_TOL);
    }

    /**
     * LagkJointMomentsFromDMAP.txt: Compute lag-k joint moments from DMAP
     * K=4, L=1 => 5x5 matrix
     * First row: [1, 1.4955, 2.9542, 7.8852, 27.282]
     */
    @Test
    public void testLagkJointMomentsFromDMAP() {
        Matrix D0 = createDmapD0();
        Matrix D1 = createDmapD1();

        Matrix Nm = LagkJointMomentsFromDMAP.lagkJointMomentsFromDMAP(D0, D1, 4, 1, 1e-14);

        assertEquals(5, Nm.getNumRows());
        assertEquals(5, Nm.getNumCols());

        double[][] expected = {
            {1.0,     1.4955,  2.9542,  7.8852,  27.282},
            {1.4955,  2.2037,  4.2827,  11.293,  38.822},
            {2.9542,  4.2875,  8.1899,  21.315,  72.753},
            {7.8852,  11.326,  21.379,  55.129,  187.21},
            {27.282,  38.993,  73.17,   187.82,  636.23}
        };

        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                assertEquals(expected[i][j], Nm.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expected[i][j]) * 1e-3),
                        "Nm[" + i + "," + j + "]");
            }
        }
    }

    /**
     * CheckDMAPRepresentation.txt Test 1: Size mismatch (D0 3x3, D1 4x4)
     * Expected: false
     */
    @Test
    public void testCheckDMAPRepresentation_SizeMismatch() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, 0.0);  D0.set(0, 1, 0.02); D0.set(0, 2, 0.0);
        D0.set(1, 0, 0.0);  D0.set(1, 1, 0.17); D0.set(1, 2, 0.2);
        D0.set(2, 0, 0.16); D0.set(2, 1, 0.17); D0.set(2, 2, 0.02);

        Matrix D1 = new Matrix(4, 4);
        D1.set(0, 0, 0.0);  D1.set(0, 1, 0.88); D1.set(0, 2, 0.1);  D1.set(0, 3, 0.0);
        D1.set(1, 0, 0.18); D1.set(1, 1, 0.07); D1.set(1, 2, 0.14); D1.set(1, 3, 0.1);
        D1.set(2, 0, 0.13); D1.set(2, 1, 0.15); D1.set(2, 2, 0.15); D1.set(2, 3, 0.04);
        D1.set(3, 0, 0.31); D1.set(3, 1, 0.18); D1.set(3, 2, 0.12); D1.set(3, 3, 0.27);

        boolean flag = CheckDMAPRepresentation.checkDMAPRepresentation(D0, D1, 1e-14);
        assertFalse(flag, "Size mismatch should be detected");
    }

    /**
     * CheckDMAPRepresentation.txt Test 2: Invalid rowsum (D0+D1 rowsum != 1)
     * Expected: false
     */
    @Test
    public void testCheckDMAPRepresentation_InvalidRowsum() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, 0.0);  D0.set(0, 1, 0.02); D0.set(0, 2, 0.0);
        D0.set(1, 0, 0.0);  D0.set(1, 1, 0.17); D0.set(1, 2, 0.2);
        D0.set(2, 0, 0.16); D0.set(2, 1, 0.17); D0.set(2, 2, 0.02);

        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 0.0);  D1.set(0, 1, 0.88); D1.set(0, 2, 0.1);
        D1.set(1, 0, 0.18); D1.set(1, 1, 0.07); D1.set(1, 2, 0.14);
        D1.set(2, 0, 0.13); D1.set(2, 1, 0.15); D1.set(2, 2, 0.15);

        boolean flag = CheckDMAPRepresentation.checkDMAPRepresentation(D0, D1, 1e-14);
        assertFalse(flag, "Invalid rowsum should be detected");
    }

    /**
     * CheckDMAPRepresentation.txt Test 3: Valid DMAP
     * Expected: true
     */
    @Test
    public void testCheckDMAPRepresentation_Valid() {
        Matrix D0 = createDmapD0();
        Matrix D1 = createDmapD1();

        boolean flag = CheckDMAPRepresentation.checkDMAPRepresentation(D0, D1, 1e-14);
        assertTrue(flag, "Valid DMAP should pass check");
    }

    /**
     * DMAP2FromMoments.txt: Construct DMAP2 from moments and correlation
     * Input: moms=[5.1536, 46.587, 626.45], corr=-0.00080286
     * Expected: D0=[0.3,0.65; 0.61538,0], D1=[0.05,0; 0.24462,0.14]
     * Validate moments roundtrip
     */
    @Test
    public void testDMAP2FromMoments() {
        // Compute moms and corr from original DMAP
        Matrix D0orig = new Matrix(2, 2);
        D0orig.set(0, 0, 0.2); D0orig.set(0, 1, 0.7);
        D0orig.set(1, 0, 0.6); D0orig.set(1, 1, 0.1);

        Matrix D1orig = new Matrix(2, 2);
        D1orig.set(0, 0, 0.09); D1orig.set(0, 1, 0.01);
        D1orig.set(1, 0, 0.2);  D1orig.set(1, 1, 0.1);

        double[] moms = MarginalMomentsFromDMAP.marginalMomentsFromDMAP(D0orig, D1orig, 3);
        assertEquals(5.1536, moms[0], LOOSE_MID_TOL);
        assertEquals(46.587, moms[1], LOOSE_MID_TOL);
        assertEquals(626.45, moms[2], LOOSE_MID_TOL);

        double[] corrArr = LagCorrelationsFromDMAP.lagCorrelationsFromDMAP(D0orig, D1orig, 1, 1e-14);
        assertEquals(-0.00080286, corrArr[0], LOOSE_MID_TOL);

        // Reconstruct from moments
        jline.util.Pair<Matrix, Matrix> result = DMAP2FromMoments.dmap2FromMoments(moms, corrArr[0]);
        Matrix D0 = result.getFirst();
        Matrix D1 = result.getSecond();

        // Verify expected values
        assertEquals(0.3, D0.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.65, D0.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.61538, D0.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.0, D0.get(1, 1), LOOSE_MID_TOL);

        assertEquals(0.05, D1.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.0, D1.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.24462, D1.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.14, D1.get(1, 1), LOOSE_MID_TOL);

        // Validate moments roundtrip
        double[] rmoms = MarginalMomentsFromDMAP.marginalMomentsFromDMAP(D0, D1, 3);
        for (int i = 0; i < 3; i++) {
            assertEquals(moms[i], rmoms[i], LOOSE_MID_TOL, "roundtrip moms[" + i + "]");
        }

        double[] rcorr = LagCorrelationsFromDMAP.lagCorrelationsFromDMAP(D0, D1, 1, 1e-14);
        assertEquals(corrArr[0], rcorr[0], LOOSE_MID_TOL, "roundtrip corr");
    }

    /**
     * CanonicalFromDMAP2.txt: Transform DMAP2 to canonical form (first test case)
     * Input: D0=[0.46,0.28; 0.35,0.23], D1=[0.08,0.18; 0.14,0.28]
     * Expected: H0=[0.6785,0.31704; 0,0.011496], H1=[0,0.004455; 0.6285,0.36]
     */
    @Test
    public void testCanonicalFromDMAP2() {
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, 0.46); D0.set(0, 1, 0.28);
        D0.set(1, 0, 0.35); D0.set(1, 1, 0.23);

        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, 0.08); D1.set(0, 1, 0.18);
        D1.set(1, 0, 0.14); D1.set(1, 1, 0.28);

        jline.util.Pair<Matrix, Matrix> result = CanonicalFromDMAP2.canonicalFromDMAP2(D0, D1, 1e-14);
        Matrix H0 = result.getFirst();
        Matrix H1 = result.getSecond();

        assertEquals(0.6785, H0.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.31704, H0.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.0, H0.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.011496, H0.get(1, 1), LOOSE_MID_TOL);

        assertEquals(0.0, H1.get(0, 0), LOOSE_MID_TOL);
        assertEquals(0.004455, H1.get(0, 1), LOOSE_MID_TOL);
        assertEquals(0.6285, H1.get(1, 0), LOOSE_MID_TOL);
        assertEquals(0.36, H1.get(1, 1), LOOSE_MID_TOL);
    }

    /**
     * DMAPFromDRAP.txt: Convert DRAP to DMAP
     * Input: H0 3x3, H1 3x3 (with negative entries)
     * Expected: D0/D1 with all non-negative entries
     */
    @Test
    public void testDMAPFromDRAP() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, 0.0);  H0.set(0, 1, 0.0);  H0.set(0, 2, 0.13);
        H0.set(1, 0, 0.0);  H0.set(1, 1, 0.6);   H0.set(1, 2, 0.18);
        H0.set(2, 0, 0.31); H0.set(2, 1, 0.26);  H0.set(2, 2, 0.02);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 0.0);  H1.set(0, 1, 1.0);   H1.set(0, 2, -0.13);
        H1.set(1, 0, 0.0);  H1.set(1, 1, 0.18);  H1.set(1, 2, 0.04);
        H1.set(2, 0, 0.03); H1.set(2, 1, 0.09);  H1.set(2, 2, 0.29);

        jline.util.Pair<Matrix, Matrix> result = DMAPFromDRAP.dmapFromDRAP(H0, H1, 1e-14);
        Matrix D0 = result.getFirst();
        Matrix D1 = result.getSecond();

        // DMAP representation from DRAP is not unique, so verify:
        // 1. Result is valid DMAP (all entries non-negative)
        boolean validDMAP = CheckDMAPRepresentation.checkDMAPRepresentation(D0, D1, 1e-12);
        assertTrue(validDMAP, "DMAPFromDRAP result should be valid DMAP");

        // 2. Marginal moments match original DRAP
        double[] origMoms = MarginalMomentsFromDRAP.marginalMomentsFromDRAP(H0, H1, 5, 1e-14);
        double[] dmapMoms = MarginalMomentsFromDMAP.marginalMomentsFromDMAP(D0, D1, 5);
        for (int i = 0; i < 5; i++) {
            assertEquals(origMoms[i], dmapMoms[i], Math.abs(origMoms[i]) * 1e-6,
                    "Moment " + (i + 1) + " mismatch between DRAP and DMAP");
        }

        // 3. Lag correlations match
        double[] origCorr = LagCorrelationsFromDRAP.lagCorrelationsFromDRAP(H0, H1, 3, 1e-14);
        double[] dmapCorr = LagCorrelationsFromDMAP.lagCorrelationsFromDMAP(D0, D1, 3, 1e-14);
        for (int i = 0; i < 3; i++) {
            assertEquals(origCorr[i], dmapCorr[i], LOOSE_MID_TOL,
                    "Lag correlation " + (i + 1) + " mismatch between DRAP and DMAP");
        }
    }

    // ============ Discrete DRAP Tests ============

    /** Helper: create the standard 3x3 DRAP H0 */
    private Matrix createDrapH0() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, 0.0);  H0.set(0, 1, 0.0);  H0.set(0, 2, 0.13);
        H0.set(1, 0, 0.0);  H0.set(1, 1, 0.6);   H0.set(1, 2, 0.18);
        H0.set(2, 0, 0.31); H0.set(2, 1, 0.26);  H0.set(2, 2, 0.02);
        return H0;
    }

    /** Helper: create the standard 3x3 DRAP H1 */
    private Matrix createDrapH1() {
        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 0.0);  H1.set(0, 1, 1.0);   H1.set(0, 2, -0.13);
        H1.set(1, 0, 0.0);  H1.set(1, 1, 0.18);  H1.set(1, 2, 0.04);
        H1.set(2, 0, 0.03); H1.set(2, 1, 0.09);  H1.set(2, 2, 0.29);
        return H1;
    }

    /**
     * MarginalMomentsFromDRAP.txt: Compute marginal moments from DRAP
     * Expected: [3.207, 16.898, 130.77, 1347.1, 17343]
     */
    @Test
    public void testMarginalMomentsFromDRAP() {
        Matrix H0 = createDrapH0();
        Matrix H1 = createDrapH1();

        double[] moms = MarginalMomentsFromDRAP.marginalMomentsFromDRAP(H0, H1, 5, 1e-14);

        double[] expected = {3.207, 16.898, 130.77, 1347.1, 17343.0};
        assertEquals(expected.length, moms.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], moms[i], Math.max(LOOSE_MID_TOL, Math.abs(expected[i]) * 1e-3), "moms[" + i + "]");
        }
    }

    /**
     * MarginalDistributionFromDRAP.txt: Compute marginal distribution from DRAP
     * Expected: a=[0.021493, 0.71253, 0.26598], A=H0
     */
    @Test
    public void testMarginalDistributionFromDRAP() {
        Matrix H0 = createDrapH0();
        Matrix H1 = createDrapH1();

        MGRepresentation result = MarginalDistributionFromDRAP.marginalDistributionFromDRAP(H0, H1, 1e-14);
        Matrix a = result.getAlpha();
        Matrix A = result.getA();

        double[] expectedAlpha = {0.021493, 0.71253, 0.26598};
        for (int i = 0; i < 3; i++) {
            assertEquals(expectedAlpha[i], a.get(0, i), LOOSE_MID_TOL, "alpha[" + i + "]");
        }

        // A should equal H0
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(H0.get(i, j), A.get(i, j), LOOSE_MID_TOL, "A[" + i + "," + j + "]");
            }
        }
    }

    /**
     * LagCorrelationsFromDRAP.txt: Compute lag correlations from DRAP
     * Expected: [0.014303, 0.0012424, 7.5868e-06]
     */
    @Test
    public void testLagCorrelationsFromDRAP() {
        Matrix H0 = createDrapH0();
        Matrix H1 = createDrapH1();

        double[] corr = LagCorrelationsFromDRAP.lagCorrelationsFromDRAP(H0, H1, 3, 1e-14);

        assertEquals(3, corr.length);
        assertEquals(0.014303, corr[0], LOOSE_MID_TOL);
        assertEquals(0.0012424, corr[1], LOOSE_MID_TOL);
        assertEquals(7.5868e-06, corr[2], LOOSE_MID_TOL);
    }

    /**
     * LagkJointMomentsFromDRAP.txt: Compute lag-k joint moments from DRAP
     * K=4, L=1 => 5x5 matrix
     * Verify via correlation roundtrip:
     *   corr[i] = (Nx(2,2) - moms(1)^2) / (moms(2) - moms(1)^2)
     */
    @Test
    public void testLagkJointMomentsFromDRAP() {
        Matrix H0 = createDrapH0();
        Matrix H1 = createDrapH1();

        Matrix Nm = LagkJointMomentsFromDRAP.lagkJointMomentsFromDRAP(H0, H1, 4, 1, 1e-14);

        assertEquals(5, Nm.getNumRows());
        assertEquals(5, Nm.getNumCols());

        // Verify correlation roundtrip
        double[] moms = MarginalMomentsFromDRAP.marginalMomentsFromDRAP(H0, H1, 4, 1e-14);
        double[] expectedCorr = {0.014303, 0.0012424, 7.5868e-06};
        double variance = moms[1] - moms[0] * moms[0];

        for (int i = 0; i < 3; i++) {
            Matrix Nx = LagkJointMomentsFromDRAP.lagkJointMomentsFromDRAP(H0, H1, 1, i + 1, 1e-14);
            double cjm = (Nx.get(1, 1) - moms[0] * moms[0]) / variance;
            assertEquals(expectedCorr[i], cjm, LOOSE_MID_TOL, "corr via joint moments[" + i + "]");
        }
    }

    /**
     * CheckDRAPRepresentation.txt Test 1: Size mismatch (D0 4x3, D1 4x3)
     * Expected: false
     */
    @Test
    public void testCheckDRAPRepresentation_SizeMismatch() {
        Matrix H0 = new Matrix(4, 3);
        H0.set(0, 0, 0.0);  H0.set(0, 1, 0.0);  H0.set(0, 2, 0.13);
        H0.set(1, 0, 0.0);  H0.set(1, 1, 0.6);   H0.set(1, 2, 0.18);
        H0.set(2, 0, 0.31); H0.set(2, 1, 0.26);  H0.set(2, 2, 0.02);
        H0.set(3, 0, 0.2);  H0.set(3, 1, 0.0);   H0.set(3, 2, 0.0);

        Matrix H1 = new Matrix(4, 3);
        H1.set(0, 0, 0.0);  H1.set(0, 1, 1.0);   H1.set(0, 2, -0.13);
        H1.set(1, 0, 0.0);  H1.set(1, 1, 0.18);  H1.set(1, 2, 0.04);
        H1.set(2, 0, 0.03); H1.set(2, 1, 0.09);  H1.set(2, 2, 0.29);
        H1.set(3, 0, 0.0);  H1.set(3, 1, 0.8);   H1.set(3, 2, 0.0);

        boolean flag = CheckDRAPRepresentation.checkDRAPRepresentation(H0, H1, 1e-14);
        assertFalse(flag, "Size mismatch should be detected");
    }

    /**
     * CheckDRAPRepresentation.txt Test 2: Invalid rowsum (D0+D1 rowsum != 1)
     * Expected: false
     */
    @Test
    public void testCheckDRAPRepresentation_InvalidRowsum() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, 0.2);  H0.set(0, 1, 0.0);  H0.set(0, 2, 0.13);
        H0.set(1, 0, 0.0);  H0.set(1, 1, 0.6);   H0.set(1, 2, 0.18);
        H0.set(2, 0, 0.31); H0.set(2, 1, 0.26);  H0.set(2, 2, 0.02);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 0.0);  H1.set(0, 1, 1.0);   H1.set(0, 2, -0.13);
        H1.set(1, 0, 0.0);  H1.set(1, 1, 0.18);  H1.set(1, 2, 0.04);
        H1.set(2, 0, 0.03); H1.set(2, 1, 0.09);  H1.set(2, 2, 0.29);

        boolean flag = CheckDRAPRepresentation.checkDRAPRepresentation(H0, H1, 1e-14);
        assertFalse(flag, "Invalid rowsum should be detected");
    }

    /**
     * CheckDRAPRepresentation.txt Test 3: Dominant eigenvalue > 1
     * Expected: false
     */
    @Test
    public void testCheckDRAPRepresentation_DominantEigenvalueGreaterThan1() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, 0.0);  H0.set(0, 1, 0.0);  H0.set(0, 2, 15.0);
        H0.set(1, 0, 0.0);  H0.set(1, 1, 0.6);   H0.set(1, 2, 0.18);
        H0.set(2, 0, 0.31); H0.set(2, 1, 0.26);  H0.set(2, 2, 0.02);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 0.0);  H1.set(0, 1, 1.0);   H1.set(0, 2, -15.0);
        H1.set(1, 0, 0.0);  H1.set(1, 1, 0.18);  H1.set(1, 2, 0.04);
        H1.set(2, 0, 0.03); H1.set(2, 1, 0.09);  H1.set(2, 2, 0.29);

        boolean flag = CheckDRAPRepresentation.checkDRAPRepresentation(H0, H1, 1e-14);
        assertFalse(flag, "Dominant eigenvalue > 1 should be detected");
    }

    /**
     * CheckDRAPRepresentation.txt Test 4: Complex dominant eigenvalue
     * Expected: false
     */
    @Test
    public void testCheckDRAPRepresentation_ComplexEigenvalue() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, 0.0);  H0.set(0, 1, 0.5);  H0.set(0, 2, 0.1);
        H0.set(1, 0, 0.0);  H0.set(1, 1, -1.4);  H0.set(1, 2, 3.1);
        H0.set(2, 0, 0.67); H0.set(2, 1, 0.0);   H0.set(2, 2, 0.4);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 0.0);  H1.set(0, 1, 0.4);   H1.set(0, 2, 0.0);
        H1.set(1, 0, 0.0);  H1.set(1, 1, -0.2);  H1.set(1, 2, -0.5);
        H1.set(2, 0, 0.3);  H1.set(2, 1, -0.7);  H1.set(2, 2, 0.33);

        boolean flag = CheckDRAPRepresentation.checkDRAPRepresentation(H0, H1, 1e-14);
        assertFalse(flag, "Complex dominant eigenvalue should be detected");
    }

    /**
     * CheckDRAPRepresentation.txt Test 5: Valid DRAP
     * Expected: true
     */
    @Test
    public void testCheckDRAPRepresentation_Valid() {
        Matrix H0 = createDrapH0();
        Matrix H1 = createDrapH1();

        boolean flag = CheckDRAPRepresentation.checkDRAPRepresentation(H0, H1, 1e-14);
        assertTrue(flag, "Valid DRAP should pass check");
    }

    // ============ DMMAP/DMRAP Construction and Conversion ============

    /** Helper: create the standard 3x3 DMRAP H0 used in DMMAP/DMRAP tests */
    private Matrix createDmrapH0() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, 0.15);  H0.set(0, 1, 0.2);  H0.set(0, 2, 0.18);
        H0.set(1, 0, -0.23); H0.set(1, 1, 0.17); H0.set(1, 2, 0.22);
        H0.set(2, 0, 0.19);  H0.set(2, 1, 0.15); H0.set(2, 2, 0.16);
        return H0;
    }

    /** Helper: create the standard 3x3 DMRAP H1 used in DMMAP/DMRAP tests */
    private Matrix createDmrapH1() {
        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 0.01); H1.set(0, 1, 0.08); H1.set(0, 2, 0.16);
        H1.set(1, 0, 0.02); H1.set(1, 1, 0.2);  H1.set(1, 2, 0.07);
        H1.set(2, 0, 0.02); H1.set(2, 1, 0.15); H1.set(2, 2, 0.17);
        return H1;
    }

    /** Helper: create the standard 3x3 DMRAP H2 used in DMMAP/DMRAP tests */
    private Matrix createDmrapH2() {
        Matrix H2 = new Matrix(3, 3);
        H2.set(0, 0, 0.14); H2.set(0, 1, 0.07); H2.set(0, 2, 0.01);
        H2.set(1, 0, 0.19); H2.set(1, 1, 0.02); H2.set(1, 2, 0.34);
        H2.set(2, 0, 0.06); H2.set(2, 1, 0.1);  H2.set(2, 2, 0.0);
        return H2;
    }

    /** Helper: create the standard DMMAP D0-D3 cell array */
    private MatrixCell createDmmapCell() {
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
        return D;
    }

    /**
     * DMMAPFromDMRAP.txt: Convert DMRAP to DMMAP
     * Input: H0,H1,H2 (with H0(1,0)=-0.20, H2(1,2)=0.31)
     * Verify: moments and joint moments of DMMAP match those of original DMRAP
     */
    @Test
    public void testDMMAPFromDMRAP_Example() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, 0.15);  H0.set(0, 1, 0.2);  H0.set(0, 2, 0.18);
        H0.set(1, 0, -0.20); H0.set(1, 1, 0.17); H0.set(1, 2, 0.22);
        H0.set(2, 0, 0.19);  H0.set(2, 1, 0.15); H0.set(2, 2, 0.16);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 0.01); H1.set(0, 1, 0.08); H1.set(0, 2, 0.16);
        H1.set(1, 0, 0.02); H1.set(1, 1, 0.2);  H1.set(1, 2, 0.07);
        H1.set(2, 0, 0.02); H1.set(2, 1, 0.15); H1.set(2, 2, 0.17);

        Matrix H2 = new Matrix(3, 3);
        H2.set(0, 0, 0.14); H2.set(0, 1, 0.07); H2.set(0, 2, 0.01);
        H2.set(1, 0, 0.19); H2.set(1, 1, 0.02); H2.set(1, 2, 0.31);
        H2.set(2, 0, 0.06); H2.set(2, 1, 0.1);  H2.set(2, 2, 0.0);

        MatrixCell H = new MatrixCell();
        H.set(0, H0);
        H.set(1, H1);
        H.set(2, H2);

        // Compute original moments for comparison
        double[] origMoms = MarginalMomentsFromDMRAP.marginalMomentsFromDMRAP(H, 5, 1e-14);
        assertEquals(1.6264, origMoms[0], LOOSE_MID_TOL, "origMoms[0]");
        assertEquals(3.6055, origMoms[1], LOOSE_MID_TOL, "origMoms[1]");

        // Convert to DMMAP
        MatrixCell G = DMMAPFromDMRAP.dmmapFromDMRAP(H, 1e-14);
        assertNotNull(G);

        // Verify result is valid DMMAP (non-negative entries)
        assertTrue(CheckDMMAPRepresentation.checkDMMAPRepresentation(G, 1e-10),
                "DMMAPFromDMRAP result should be valid DMMAP");

        // Verify moments match
        double[] dmmapMoms = MarginalMomentsFromDMMAP.marginalMomentsFromDMMAP(G, 5, 1e-14);
        for (int i = 0; i < 5; i++) {
            assertEquals(origMoms[i], dmmapMoms[i], Math.max(LOOSE_MID_TOL, Math.abs(origMoms[i]) * 1e-3),
                    "DMMAP moms[" + i + "]");
        }

        // Verify joint moments match
        MatrixCell origJmom = LagkJointMomentsFromDMRAP.lagkJointMomentsFromDMRAP(H, 3, 1, 1e-14);
        MatrixCell dmmapJmom = LagkJointMomentsFromDMMAP.lagkJointMomentsFromDMMAP(G, 3, 1, 1e-14);
        for (int t = 0; t < 2; t++) {
            Matrix origNm = origJmom.get(t);
            Matrix dmmapNm = dmmapJmom.get(t);
            for (int i = 0; i < origNm.getNumRows(); i++) {
                for (int j = 0; j < origNm.getNumCols(); j++) {
                    assertEquals(origNm.get(i, j), dmmapNm.get(i, j),
                            Math.max(LOOSE_MID_TOL, Math.abs(origNm.get(i, j)) * 1e-3),
                            "JointMom[" + t + "][" + i + "," + j + "]");
                }
            }
        }
    }

    /**
     * DMRAPFromMoments.txt: Reconstruct DMRAP from moments
     * Uses DMMAP D0-D3 to compute moments, then reconstructs DMRAP.
     * Verify output matrices match expected values from data file.
     */
    @Test
    public void testDMRAPFromMoments_Example() {
        MatrixCell D = createDmmapCell();

        // Compute moments from DMRAP (treated as DMRAP with 4 matrices)
        double[] moms = MarginalMomentsFromDMRAP.marginalMomentsFromDMRAP(D, 5, 1e-14);
        assertEquals(1.5037, moms[0], LOOSE_MID_TOL, "moms[0]");
        assertEquals(3.0278, moms[1], LOOSE_MID_TOL, "moms[1]");

        // Compute lag-k joint moments
        MatrixCell Nm = LagkJointMomentsFromDMRAP.lagkJointMomentsFromDMRAP(D, 2, 1, 1e-14);
        assertNotNull(Nm);

        // Reconstruct DMRAP from moments
        MatrixCell H = DMRAPFromMoments.dmrapFromMoments(moms, Nm);
        assertNotNull(H);
        assertEquals(4, H.size(), "Expected 4 matrices in reconstructed DMRAP");

        // Verify H{1} values
        double[][] expectedH1 = {
            {0.067795, 0.67922, -0.42509},
            {-0.018716, 0.0035507, 0.35386},
            {-0.019902, 0.039925, 0.31865}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expectedH1[i][j], H.get(0).get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedH1[i][j]) * 1e-2),
                        "H{1}[" + i + "," + j + "]");
            }
        }

        // Verify H{2} values
        double[][] expectedH2 = {
            {0.19512, 1.6571, -1.5453},
            {-0.0081004, 0.18194, 0.12658},
            {-0.0076905, 0.11526, 0.19294}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expectedH2[i][j], H.get(1).get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedH2[i][j]) * 1e-2),
                        "H{2}[" + i + "," + j + "]");
            }
        }

        // Verify H{3} values
        double[][] expectedH3 = {
            {0.25787, 3.8341, -4.0555},
            {0.12827, -1.3508, 1.2347},
            {0.12666, -1.4076, 1.2929}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expectedH3[i][j], H.get(2).get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedH3[i][j]) * 1e-2),
                        "H{3}[" + i + "," + j + "]");
            }
        }

        // Verify H{4} values
        double[][] expectedH4 = {
            {0.17441, -0.11574, 0.27606},
            {-0.014006, -0.013671, 0.37636},
            {-0.012887, -0.0175, 0.37926}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expectedH4[i][j], H.get(3).get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedH4[i][j]) * 1e-2),
                        "H{4}[" + i + "," + j + "]");
            }
        }
    }

    /**
     * DRAPFromMoments.txt: Reconstruct DRAP from moments
     * Uses DMAP D0/D1 to compute moments, then reconstructs DRAP.
     * Verify output matrices match expected values from data file.
     */
    @Test
    public void testDRAPFromMoments_Example() {
        Matrix D0 = createDmapD0();
        Matrix D1 = createDmapD1();

        // Compute moments
        double[] moms = MarginalMomentsFromDRAP.marginalMomentsFromDRAP(D0, D1, 5, 1e-14);
        assertEquals(1.4955, moms[0], LOOSE_MID_TOL, "moms[0]");
        assertEquals(2.9542, moms[1], LOOSE_MID_TOL, "moms[1]");

        // Compute lag-k joint moments (K=2, L=1 → 3x3 matrix)
        Matrix Nm = LagkJointMomentsFromDRAP.lagkJointMomentsFromDRAP(D0, D1, 2, 1, 1e-14);
        assertNotNull(Nm);

        // Expected Nm values
        double[][] expectedNm = {
            {1.0, 1.4955, 2.9542},
            {1.4955, 2.2037, 4.2827},
            {2.9542, 4.2875, 8.1899}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expectedNm[i][j], Nm.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm[i][j]) * 1e-3),
                        "Nm[" + i + "," + j + "]");
            }
        }

        // Reconstruct DRAP from moments
        jline.util.Pair<Matrix, Matrix> result = DRAPFromMoments.drapFromMoments(moms, Nm);
        Matrix H0 = result.getFirst();
        Matrix H1 = result.getSecond();

        assertNotNull(H0);
        assertNotNull(H1);
        assertEquals(3, H0.getNumRows());
        assertEquals(3, H0.getNumCols());

        // Verify expected H0 values
        double[][] expectedH0 = {
            {0.56447, 0.47188, -0.69474},
            {-0.50857, -0.10551, 0.95921},
            {0.18477, 0.26121, -0.13431}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expectedH0[i][j], H0.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedH0[i][j]) * 1e-2),
                        "H0[" + i + "," + j + "]");
            }
        }

        // Verify expected H1 values
        double[][] expectedH1 = {
            {2.3994, 1.1243, -2.8653},
            {-1.7535, -0.59009, 2.9984},
            {0.95074, 0.51879, -0.7812}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expectedH1[i][j], H1.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedH1[i][j]) * 1e-2),
                        "H1[" + i + "," + j + "]");
            }
        }
    }

    // ============ Lag-k Joint Moments from DMMAP/DMRAP ============

    /**
     * LagkJointMomentsFromDMMAP.txt: Compute lag-k joint moments from DMMAP
     * Input: D0-D3 (same as existing DMMAP tests), K=3, L=1
     * Expected: 3 matrices (Nm{1}, Nm{2}, Nm{3}), each 4x4
     */
    @Test
    public void testLagkJointMomentsFromDMMAP_Example() {
        MatrixCell D = createDmmapCell();

        MatrixCell Nm = LagkJointMomentsFromDMMAP.lagkJointMomentsFromDMMAP(D, 3, 1, 1e-14);
        assertNotNull(Nm);

        // Expected Nm{1} (4x4)
        double[][] expectedNm1 = {
            {0.45395, 0.68525, 1.3856, 3.8671},
            {0.68283, 1.0318, 2.0887, 5.8339},
            {1.3755, 2.0807, 4.2171, 11.789},
            {3.828, 5.7954, 11.756, 32.887}
        };
        Matrix Nm1 = Nm.get(0);
        assertEquals(4, Nm1.getNumRows());
        assertEquals(4, Nm1.getNumCols());
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedNm1[i][j], Nm1.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm1[i][j]) * 1e-3),
                        "Nm{1}[" + i + "," + j + "]");
            }
        }

        // Expected Nm{2} (4x4)
        double[][] expectedNm2 = {
            {0.026281, 0.03323, 0.053055, 0.11925},
            {0.035051, 0.043866, 0.068917, 0.15222},
            {0.060653, 0.07477, 0.11464, 0.24631},
            {0.1482, 0.17996, 0.26901, 0.56074}
        };
        Matrix Nm2 = Nm.get(1);
        assertEquals(4, Nm2.getNumRows());
        assertEquals(4, Nm2.getNumCols());
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedNm2[i][j], Nm2.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm2[i][j]) * 1e-3),
                        "Nm{2}[" + i + "," + j + "]");
            }
        }

        // Expected Nm{3} (4x4)
        double[][] expectedNm3 = {
            {0.51977, 0.78522, 1.5891, 4.438},
            {0.78582, 1.1881, 2.4067, 6.7254},
            {1.5917, 2.4087, 4.8838, 13.657},
            {4.4481, 6.7354, 13.666, 38.235}
        };
        Matrix Nm3 = Nm.get(2);
        assertEquals(4, Nm3.getNumRows());
        assertEquals(4, Nm3.getNumCols());
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedNm3[i][j], Nm3.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm3[i][j]) * 1e-3),
                        "Nm{3}[" + i + "," + j + "]");
            }
        }
    }

    /**
     * LagkJointMomentsFromDMRAP.txt: Compute lag-k joint moments from DMRAP
     * Input: H0,H1,H2 (standard DMRAP), K=3, L=2
     * Expected: 2 matrices (Nm{1}, Nm{2}), each 4x4
     */
    @Test
    public void testLagkJointMomentsFromDMRAP_Example() {
        MatrixCell H = new MatrixCell();
        H.set(0, createDmrapH0());
        H.set(1, createDmrapH1());
        H.set(2, createDmrapH2());

        MatrixCell Nm = LagkJointMomentsFromDMRAP.lagkJointMomentsFromDMRAP(H, 3, 2, 1e-14);
        assertNotNull(Nm);

        // Expected Nm{1} (4x4)
        double[][] expectedNm1 = {
            {0.48798, 0.78047, 1.6785, 4.9029},
            {0.77458, 1.2395, 2.6673, 7.7945},
            {1.6539, 2.6481, 5.7016, 16.669},
            {4.8092, 7.7033, 16.593, 48.526}
        };
        Matrix Nm1 = Nm.get(0);
        assertEquals(4, Nm1.getNumRows());
        assertEquals(4, Nm1.getNumCols());
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedNm1[i][j], Nm1.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm1[i][j]) * 1e-3),
                        "Nm{1}[" + i + "," + j + "]");
            }
        }

        // Expected Nm{2} (4x4)
        double[][] expectedNm2 = {
            {0.51202, 0.81429, 1.7401, 5.0566},
            {0.82019, 1.3036, 2.7837, 8.0853},
            {1.7647, 2.8029, 5.9814, 17.365},
            {5.1503, 8.177, 17.442, 50.619}
        };
        Matrix Nm2 = Nm.get(1);
        assertEquals(4, Nm2.getNumRows());
        assertEquals(4, Nm2.getNumCols());
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedNm2[i][j], Nm2.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm2[i][j]) * 1e-3),
                        "Nm{2}[" + i + "," + j + "]");
            }
        }
    }

    /**
     * MarginalDistributionFromDMRAP.txt: Compute marginal distribution from DMRAP
     * Input: H0,H1,H2 (standard DMRAP)
     * Expected: a=[0.22615, 0.35424, 0.41962], A=H0
     */
    @Test
    public void testMarginalDistributionFromDMRAP_Example() {
        MatrixCell H = new MatrixCell();
        H.set(0, createDmrapH0());
        H.set(1, createDmrapH1());
        H.set(2, createDmrapH2());

        MGRepresentation result = MarginalDistributionFromDMRAP.marginalDistributionFromDMRAP(H, 1e-14);
        assertNotNull(result);

        Matrix a = result.getAlpha();
        Matrix A = result.getA();

        assertNotNull(a);
        assertNotNull(A);

        // Verify initial vector a
        assertEquals(1, a.getNumRows());
        assertEquals(3, a.getNumCols());
        assertEquals(0.22615, a.get(0, 0), LOOSE_MID_TOL, "a[0]");
        assertEquals(0.35424, a.get(0, 1), LOOSE_MID_TOL, "a[1]");
        assertEquals(0.41962, a.get(0, 2), LOOSE_MID_TOL, "a[2]");

        // Verify A = H0
        Matrix H0 = createDmrapH0();
        assertEquals(3, A.getNumRows());
        assertEquals(3, A.getNumCols());
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(H0.get(i, j), A.get(i, j), LOOSE_MID_TOL, "A[" + i + "," + j + "]");
            }
        }
    }

    // ============ Random Generator Tests ============

    /**
     * RandomMAP.txt: Generate random MAP with order=4, mean=1.62, zeroEntries=10
     * Verify: valid MAP representation, correct mean
     */
    @Test
    public void testRandomMAP() {
        MAPRepresentation result = RandomMAP.randomMAP(4, 1.62, 10, 1000, 1e-7, new Random(42));
        Matrix D0 = result.getD0();
        Matrix D1 = result.getD1();

        // Verify dimensions
        assertEquals(4, D0.getNumRows());
        assertEquals(4, D0.getNumCols());
        assertEquals(4, D1.getNumRows());
        assertEquals(4, D1.getNumCols());

        // Verify valid MAP representation
        assertTrue(CheckMAPRepresentation.checkMAPRepresentation(D0, D1, 1e-12),
                "RandomMAP should produce valid MAP");

        // Verify mean
        double[] moms = MarginalMomentsFromMAP.marginalMomentsFromMAP(D0, D1, 1);
        assertEquals(1.62, moms[0], COARSE_TOL, "Mean should be close to 1.62");
    }

    /**
     * RandomDMAP.txt: Generate random DMAP with order=4, mean=5.62, zeroEntries=10
     * Verify: valid DMAP representation, correct mean
     */
    @Test
    public void testRandomDMAP() {
        Pair<Matrix, Matrix> result = RandomDMAP.randomDMAP(4, 5.62, 10, 1000, 1e-7, new Random(42));
        Matrix D0 = result.getFirst();
        Matrix D1 = result.getSecond();

        // Verify dimensions
        assertEquals(4, D0.getNumRows());
        assertEquals(4, D0.getNumCols());
        assertEquals(4, D1.getNumRows());
        assertEquals(4, D1.getNumCols());

        // Verify valid DMAP representation
        assertTrue(CheckDMAPRepresentation.checkDMAPRepresentation(D0, D1, 1e-14),
                "RandomDMAP should produce valid DMAP");

        // Verify mean
        double[] moms = MarginalMomentsFromDMAP.marginalMomentsFromDMAP(D0, D1, 1);
        assertEquals(5.62, moms[0], COARSE_TOL, "Mean should be close to 5.62");
    }

    /**
     * RandomDMMAP.txt: Generate random DMMAP with order=4, types=3, mean=5.62, zeroEntries=10
     * Verify: valid DMMAP representation, correct mean
     */
    @Test
    public void testRandomDMMAP() {
        MatrixCell result = RandomDMMAP.randomDMMAP(4, 3, 5.62, 10, 1000, 1e-7, new Random(42));

        // Verify result has types+1 matrices (D0, D1, ..., Dtypes)
        assertNotNull(result);

        // Verify D0 dimensions
        Matrix D0 = result.get(0);
        assertEquals(4, D0.getNumRows());
        assertEquals(4, D0.getNumCols());

        // Verify valid DMMAP representation
        assertTrue(CheckDMMAPRepresentation.checkDMMAPRepresentation(result, 1e-14),
                "RandomDMMAP should produce valid DMMAP");

        // Verify mean via marginal moments
        double[] moms = MarginalMomentsFromDMMAP.marginalMomentsFromDMMAP(result, 1, 1e-14);
        assertEquals(5.62, moms[0], COARSE_TOL, "Mean should be close to 5.62");
    }

    // ============ Sampling Tests ============

    /**
     * SamplesFromMAP.txt: Generate samples from MAP, verify moments match theoretical
     * Input: D0/D1 (4x4 MAP)
     * Theoretical moments: [3.4433, 34.03, 592.08]
     */
    @Test
    public void testSamplesFromMAP() {
        Matrix D0 = new Matrix(4, 4);
        D0.set(0, 0, -0.17); D0.set(0, 1, 0.0);   D0.set(0, 2, 0.0);  D0.set(0, 3, 0.07);
        D0.set(1, 0, 0.01);  D0.set(1, 1, -0.78);  D0.set(1, 2, 0.03); D0.set(1, 3, 0.08);
        D0.set(2, 0, 0.22);  D0.set(2, 1, 0.17);   D0.set(2, 2, -1.1); D0.set(2, 3, 0.02);
        D0.set(3, 0, 0.04);  D0.set(3, 1, 0.12);   D0.set(3, 2, 0.0);  D0.set(3, 3, -0.42);

        Matrix D1 = new Matrix(4, 4);
        D1.set(0, 0, 0.0);  D1.set(0, 1, 0.06); D1.set(0, 2, 0.0);  D1.set(0, 3, 0.04);
        D1.set(1, 0, 0.04); D1.set(1, 1, 0.19); D1.set(1, 2, 0.21); D1.set(1, 3, 0.22);
        D1.set(2, 0, 0.22); D1.set(2, 1, 0.13); D1.set(2, 2, 0.15); D1.set(2, 3, 0.19);
        D1.set(3, 0, 0.05); D1.set(3, 1, 0.0);  D1.set(3, 2, 0.17); D1.set(3, 3, 0.04);

        // Generate 10000 samples
        double[] samples = SamplesFromMAP.samplesFromMAP(D0, D1, 10000, null, new Random(42));
        assertEquals(10000, samples.length);

        // Compute trace moments
        double[] moms = MarginalMomentsFromTrace.marginalMomentsFromTrace(samples, 3);

        // Theoretical moments
        double[] theoretical = {3.4433, 34.03, 592.08};

        // Stochastic: trace moments should be within 20% of theoretical
        for (int i = 0; i < 3; i++) {
            assertEquals(theoretical[i], moms[i], Math.abs(theoretical[i]) * 0.2,
                    "SamplesFromMAP moms[" + i + "]");
        }
    }

    /**
     * SamplesFromDMAP.txt: Generate samples from DMAP, verify moments match theoretical
     * Input: D0/D1 (4x4 DMAP)
     * Theoretical moments: [1.4955, 2.9542, 7.8852]
     */
    @Test
    public void testSamplesFromDMAP() {
        Matrix D0 = createDmapD0();
        Matrix D1 = createDmapD1();

        // Generate 10000 samples
        int[] samples = SamplesFromDMAP.samplesFromDMAP(D0, D1, 10000, null, 1e-14, new Random(42));
        assertEquals(10000, samples.length);

        // Convert to double array for trace analysis
        double[] dsamples = new double[samples.length];
        for (int i = 0; i < samples.length; i++) {
            dsamples[i] = samples[i];
        }

        // Compute trace moments
        double[] moms = MarginalMomentsFromTrace.marginalMomentsFromTrace(dsamples, 3);

        // Theoretical moments
        double[] theoretical = {1.4955, 2.9542, 7.8852};

        // Stochastic: trace moments should be within 25% of theoretical
        for (int i = 0; i < 3; i++) {
            assertEquals(theoretical[i], moms[i], Math.abs(theoretical[i]) * 0.25,
                    "SamplesFromDMAP moms[" + i + "]");
        }
    }

    /**
     * SamplesFromDMMAP.txt: Generate samples from DMMAP, verify moments match theoretical
     * Input: D0,D1,D2,D3 (3x3 DMMAP)
     * Theoretical moments: [1.5037, 3.0278, 8.4243]
     */
    @Test
    public void testSamplesFromDMMAP() {
        MatrixCell D = createDmmapCell();

        // Generate 10000 samples
        Object result = SamplesFromDMMAP.samplesFromDMMAP(D, 10000, null, 1e-14, new Random(42));

        // Result is int[] for single-type or int[][] for multi-type
        // For DMMAP, result is typically multi-type: Array<IntArray> = int[][]
        // The first column (inter-arrival times) should be used for moment checking
        int[] arrivals;
        if (result instanceof int[]) {
            arrivals = (int[]) result;
        } else {
            // Multi-type: first dimension is samples, second is [interarrival, type]
            int[][] multi = (int[][]) result;
            arrivals = new int[multi.length];
            for (int i = 0; i < multi.length; i++) {
                arrivals[i] = multi[i][0];
            }
        }

        // Convert to double for trace analysis
        double[] dsamples = new double[arrivals.length];
        for (int i = 0; i < arrivals.length; i++) {
            dsamples[i] = arrivals[i];
        }

        // Compute trace moments
        double[] moms = MarginalMomentsFromTrace.marginalMomentsFromTrace(dsamples, 3);

        // Theoretical moments
        double[] theoretical = {1.5037, 3.0278, 8.4243};

        // Stochastic: trace moments should be within 25% of theoretical
        for (int i = 0; i < 3; i++) {
            assertEquals(theoretical[i], moms[i], Math.abs(theoretical[i]) * 0.25,
                    "SamplesFromDMMAP moms[" + i + "]");
        }
    }

    // ============ Continuous MMAP/MRAP Tests ============

    /** Helper: create the standard 4x4 MRAP H0 used in MMAP/MRAP tests (x=0.18) */
    private Matrix createMrapH0() {
        Matrix H0 = new Matrix(4, 4);
        H0.set(0, 0, -5.0);  H0.set(0, 1, 0.28); H0.set(0, 2, 0.9);  H0.set(0, 3, 1.0);
        H0.set(1, 0, 1.0);   H0.set(1, 1, -8.0); H0.set(1, 2, 0.9);  H0.set(1, 3, 0.1);
        H0.set(2, 0, 0.9);   H0.set(2, 1, 0.1);  H0.set(2, 2, -4.0); H0.set(2, 3, 1.0);
        H0.set(3, 0, 1.0);   H0.set(3, 1, 2.0);  H0.set(3, 2, 3.0);  H0.set(3, 3, -9.0);
        return H0;
    }

    /** Helper: create the standard 4x4 MRAP H1 used in MMAP/MRAP tests (x=0.18) */
    private Matrix createMrapH1() {
        // H1 = [0.1-x, 0.7, 0.1, 0.1; 0.1, 1, 1.8, 0.1; 0.1, 0.1, 0.1, 0.7; 0.7, 0.1, 0.1, 0.1]
        Matrix H1 = new Matrix(4, 4);
        H1.set(0, 0, -0.08); H1.set(0, 1, 0.7);  H1.set(0, 2, 0.1);  H1.set(0, 3, 0.1);
        H1.set(1, 0, 0.1);   H1.set(1, 1, 1.0);  H1.set(1, 2, 1.8);  H1.set(1, 3, 0.1);
        H1.set(2, 0, 0.1);   H1.set(2, 1, 0.1);  H1.set(2, 2, 0.1);  H1.set(2, 3, 0.7);
        H1.set(3, 0, 0.7);   H1.set(3, 1, 0.1);  H1.set(3, 2, 0.1);  H1.set(3, 3, 0.1);
        return H1;
    }

    /** Helper: create the standard 4x4 MRAP H2 used in MMAP/MRAP tests */
    private Matrix createMrapH2() {
        Matrix H2 = new Matrix(4, 4);
        H2.set(0, 0, 0.1);  H2.set(0, 1, 0.1);  H2.set(0, 2, 0.1);  H2.set(0, 3, 1.7);
        H2.set(1, 0, 1.8);  H2.set(1, 1, 0.1);  H2.set(1, 2, 1.0);  H2.set(1, 3, 0.1);
        H2.set(2, 0, 0.1);  H2.set(2, 1, 0.1);  H2.set(2, 2, 0.7);  H2.set(2, 3, 0.1);
        H2.set(3, 0, 0.1);  H2.set(3, 1, 1.0);  H2.set(3, 2, 0.1);  H2.set(3, 3, 0.8);
        return H2;
    }

    /** Helper: create the standard MMAP D0-D3 used in continuous MMAP tests */
    private MatrixCell createMmapCell() {
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -1.78); D0.set(0, 1, 0.29);
        D0.set(1, 0, 0.07);  D0.set(1, 1, -0.92);

        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, 0.15); D1.set(0, 1, 0.49);
        D1.set(1, 0, 0.23); D1.set(1, 1, 0.36);

        Matrix D2 = new Matrix(2, 2);
        D2.set(0, 0, 0.11); D2.set(0, 1, 0.2);
        D2.set(1, 0, 0.01); D2.set(1, 1, 0.0);

        Matrix D3 = new Matrix(2, 2);
        D3.set(0, 0, 0.14); D3.set(0, 1, 0.4);
        D3.set(1, 0, 0.11); D3.set(1, 1, 0.14);

        MatrixCell D = new MatrixCell();
        D.set(0, D0);
        D.set(1, D1);
        D.set(2, D2);
        D.set(3, D3);
        return D;
    }

    /**
     * CheckMMAPRepresentation.txt: Validate continuous MMAP
     * Input: D0-D3 (3x3), all valid
     * Expected: true
     */
    @Test
    public void testCheckMMAPRepresentation_Valid() {
        Matrix D0 = new Matrix(3, 3);
        D0.set(0, 0, -1.05); D0.set(0, 1, 0.03); D0.set(0, 2, 0.07);
        D0.set(1, 0, 0.19);  D0.set(1, 1, -1.63); D0.set(1, 2, 0.06);
        D0.set(2, 0, 0.0);   D0.set(2, 1, 0.2);  D0.set(2, 2, -1.03);

        Matrix D1 = new Matrix(3, 3);
        D1.set(0, 0, 0.16); D1.set(0, 1, 0.11); D1.set(0, 2, 0.0);
        D1.set(1, 0, 0.1);  D1.set(1, 1, 0.16); D1.set(1, 2, 0.0);
        D1.set(2, 0, 0.27); D1.set(2, 1, 0.0);  D1.set(2, 2, 0.19);

        Matrix D2 = new Matrix(3, 3);
        D2.set(0, 0, 0.01); D2.set(0, 1, 0.09); D2.set(0, 2, 0.13);
        D2.set(1, 0, 0.26); D2.set(1, 1, 0.21); D2.set(1, 2, 0.05);
        D2.set(2, 0, 0.0);  D2.set(2, 1, 0.16); D2.set(2, 2, 0.07);

        Matrix D3 = new Matrix(3, 3);
        D3.set(0, 0, 0.19); D3.set(0, 1, 0.06); D3.set(0, 2, 0.2);
        D3.set(1, 0, 0.17); D3.set(1, 1, 0.16); D3.set(1, 2, 0.27);
        D3.set(2, 0, 0.0);  D3.set(2, 1, 0.0);  D3.set(2, 2, 0.14);

        MatrixCell D = new MatrixCell();
        D.set(0, D0);
        D.set(1, D1);
        D.set(2, D2);
        D.set(3, D3);

        boolean flag = CheckMMAPRepresentation.checkMMAPRepresentation(D, 1e-14);
        assertTrue(flag, "Valid MMAP should pass check");
    }

    /**
     * CheckMRAPRepresentation.txt: Validate continuous MRAP
     * Input: H0,H1,H2 (4x4), valid (with negative entries in H1)
     * Expected: true
     */
    @Test
    public void testCheckMRAPRepresentation_Valid() {
        MatrixCell H = new MatrixCell();
        H.set(0, createMrapH0());
        H.set(1, createMrapH1());
        H.set(2, createMrapH2());

        boolean flag = CheckMRAPRepresentation.checkMRAPRepresentation(H, 1e-14);
        assertTrue(flag, "Valid MRAP should pass check");
    }

    /**
     * MarginalDistributionFromMRAP.txt: Compute marginal distribution from MRAP
     * Expected: a=[0.17159, 0.21695, 0.27936, 0.3321], A=H0
     */
    @Test
    public void testMarginalDistributionFromMRAP() {
        MatrixCell H = new MatrixCell();
        H.set(0, createMrapH0());
        H.set(1, createMrapH1());
        H.set(2, createMrapH2());

        PHRepresentation result = MarginalDistributionFromMRAP.marginalDistributionFromMRAP(H, 1e-14);
        Matrix a = result.getAlpha();
        Matrix A = result.getA();

        double[] expectedAlpha = {0.17159, 0.21695, 0.27936, 0.3321};
        for (int i = 0; i < 4; i++) {
            assertEquals(expectedAlpha[i], a.get(0, i), LOOSE_MID_TOL, "alpha[" + i + "]");
        }

        // A should equal H0
        Matrix H0 = createMrapH0();
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(H0.get(i, j), A.get(i, j), LOOSE_MID_TOL, "A[" + i + "," + j + "]");
            }
        }
    }

    /**
     * MarginalDistributionFromMMAP.txt: Compute marginal distribution from MMAP
     * Expected: a=[0.36191, 0.63809], A=D0
     */
    @Test
    public void testMarginalDistributionFromMMAP() {
        MatrixCell D = createMmapCell();

        PHRepresentation result = MarginalDistributionFromMMAP.marginalDistributionFromMMAP(D, 1e-14);
        Matrix a = result.getAlpha();
        Matrix A = result.getA();

        double[] expectedAlpha = {0.36191, 0.63809};
        for (int i = 0; i < 2; i++) {
            assertEquals(expectedAlpha[i], a.get(0, i), LOOSE_MID_TOL, "alpha[" + i + "]");
        }

        // A should equal D0
        Matrix D0 = D.get(0);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(D0.get(i, j), A.get(i, j), LOOSE_MID_TOL, "A[" + i + "," + j + "]");
            }
        }
    }

    /**
     * MarginalMomentsFromMRAP.txt: Compute marginal moments from MRAP
     * Expected: [0.33951, 0.24583, 0.27424, 0.41206, 0.77677, 1.7594, 4.6515]
     */
    @Test
    public void testMarginalMomentsFromMRAP() {
        MatrixCell H = new MatrixCell();
        H.set(0, createMrapH0());
        H.set(1, createMrapH1());
        H.set(2, createMrapH2());

        double[] moms = MarginalMomentsFromMRAP.marginalMomentsFromMRAP(H, 7);

        double[] expected = {0.33951, 0.24583, 0.27424, 0.41206, 0.77677, 1.7594, 4.6515};
        assertEquals(expected.length, moms.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], moms[i], Math.max(LOOSE_MID_TOL, Math.abs(expected[i]) * 1e-3), "moms[" + i + "]");
        }
    }

    /**
     * MarginalMomentsFromMMAP.txt: Compute marginal moments from MMAP
     * Expected: [1.0007, 2.1045, 6.8277]
     */
    @Test
    public void testMarginalMomentsFromMMAP() {
        MatrixCell D = createMmapCell();

        double[] moms = MarginalMomentsFromMMAP.marginalMomentsFromMMAP(D, 3);

        double[] expected = {1.0007, 2.1045, 6.8277};
        assertEquals(expected.length, moms.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], moms[i], Math.max(LOOSE_MID_TOL, Math.abs(expected[i]) * 1e-3), "moms[" + i + "]");
        }
    }

    /**
     * LagkJointMomentsFromMAP.txt: Compute lag-k joint moments from MAP
     * Input: D0 4x4, D1 4x4 (same as LagCorrelationsFromMAP test)
     * K=4, L=1 => 5x5 matrix
     */
    @Test
    public void testLagkJointMomentsFromMAP() {
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

        Matrix Nm = LagkJointMomentsFromMAP.lagkJointMomentsFromMAP(D0, D1, 4, 1, 1e-14);

        assertEquals(5, Nm.getNumRows());
        assertEquals(5, Nm.getNumCols());

        double[][] expected = {
            {1.0,     0.34247, 0.25054, 0.28271, 0.42984},
            {0.34247, 0.1173,  0.085789, 0.096807, 0.14721},
            {0.25054, 0.0857,  0.062633, 0.07066, 0.10744},
            {0.28271, 0.096627, 0.070589, 0.079623, 0.12107},
            {0.42984, 0.14686, 0.10727, 0.12099, 0.18396}
        };
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                assertEquals(expected[i][j], Nm.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expected[i][j]) * 1e-3),
                        "Nm[" + i + "," + j + "]");
            }
        }
    }

    /**
     * LagkJointMomentsFromMRAP.txt: Compute lag-k joint moments from MRAP
     * Input: H0,H1,H2 (4x4), K=3, L=2
     * Expected: 2 matrices Nm{1}, Nm{2}, each 4x4
     */
    @Test
    public void testLagkJointMomentsFromMRAP() {
        MatrixCell H = new MatrixCell();
        H.set(0, createMrapH0());
        H.set(1, createMrapH1());
        H.set(2, createMrapH2());

        MatrixCell Nm = LagkJointMomentsFromMRAP.lagkJointMomentsFromMRAP(H, 3, 2);
        assertNotNull(Nm);

        // Expected Nm{1} (4x4)
        double[][] expectedNm1 = {
            {0.41974, 0.14337, 0.1041, 0.11625},
            {0.14138, 0.048248, 0.035017, 0.0391},
            {0.10186, 0.034737, 0.025205, 0.02814},
            {0.11338, 0.038655, 0.028044, 0.031308}
        };
        Matrix Nm1 = Nm.get(0);
        assertEquals(4, Nm1.getNumRows());
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedNm1[i][j], Nm1.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm1[i][j]) * 1e-3),
                        "Nm{1}[" + i + "," + j + "]");
            }
        }

        // Expected Nm{2} (4x4)
        double[][] expectedNm2 = {
            {0.58026, 0.19614, 0.14173, 0.15799},
            {0.19813, 0.066994, 0.048418, 0.053974},
            {0.14397, 0.048697, 0.035199, 0.03924},
            {0.16086, 0.054419, 0.039338, 0.043855}
        };
        Matrix Nm2 = Nm.get(1);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedNm2[i][j], Nm2.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm2[i][j]) * 1e-3),
                        "Nm{2}[" + i + "," + j + "]");
            }
        }
    }

    /**
     * LagkJointMomentsFromMMAP.txt: Compute lag-k joint moments from continuous MMAP
     * Input: D0-D3 (2x2), K=3, L=1
     * Expected: 3 matrices, each 4x4
     */
    @Test
    public void testLagkJointMomentsFromMMAP() {
        MatrixCell D = createMmapCell();

        MatrixCell Nm = LagkJointMomentsFromMMAP.lagkJointMomentsFromMMAP(D, 3, 1);
        assertNotNull(Nm);

        // Expected Nm{1} (4x4)
        double[][] expectedNm1 = {
            {0.60207, 0.60501, 1.2755, 4.1438},
            {0.62913, 0.62913, 1.3226, 4.2901},
            {1.3561, 1.3524, 2.8387, 9.1998},
            {4.4576, 4.4395, 9.3105, 30.16}
        };
        Matrix Nm1 = Nm.get(0);
        assertEquals(4, Nm1.getNumRows());
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedNm1[i][j], Nm1.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm1[i][j]) * 1e-3),
                        "Nm{1}[" + i + "," + j + "]");
            }
        }

        // Expected Nm{2} (4x4)
        double[][] expectedNm2 = {
            {0.080053, 0.078372, 0.16268, 0.52401},
            {0.06033, 0.058276, 0.11997, 0.38467},
            {0.10244, 0.097662, 0.1994, 0.63637},
            {0.28923, 0.27293, 0.5536, 1.7601}
        };
        Matrix Nm2 = Nm.get(1);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedNm2[i][j], Nm2.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm2[i][j]) * 1e-3),
                        "Nm{2}[" + i + "," + j + "]");
            }
        }

        // Expected Nm{3} (4x4)
        double[][] expectedNm3 = {
            {0.31788, 0.31729, 0.66629, 2.1599},
            {0.31121, 0.30821, 0.64424, 2.0831},
            {0.646, 0.63672, 1.3271, 4.2844},
            {2.0808, 2.0455, 4.2565, 13.73}
        };
        Matrix Nm3 = Nm.get(2);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(expectedNm3[i][j], Nm3.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedNm3[i][j]) * 1e-3),
                        "Nm{3}[" + i + "," + j + "]");
            }
        }
    }

    /**
     * MAPFromRAP.txt: Convert RAP to MAP (test case 2)
     * Input: D0=[-2.4,2; 2,-9], D1=[-1.6,2; 3,4]
     * Expected: H0=[-1.8414,0.079468; 0.012509,-9.5586]
     * Verify: output matches expected values and joint moments preserved
     */
    @Test
    public void testMAPFromRAP() {
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -2.4); D0.set(0, 1, 2.0);
        D0.set(1, 0, 2.0);  D0.set(1, 1, -9.0);

        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, -1.6); D1.set(0, 1, 2.0);
        D1.set(1, 0, 3.0);  D1.set(1, 1, 4.0);

        Pair<Matrix, Matrix> result = MAPFromRAP.mapFromRAP(D0, D1, 1e-14);
        Matrix H0 = result.getFirst();
        Matrix H1 = result.getSecond();

        // Verify valid MAP (this test case should produce a valid MAP)
        assertTrue(CheckMAPRepresentation.checkMAPRepresentation(H0, H1, 1e-7),
                "MAPFromRAP result should be valid MAP");

        // Verify joint moments are preserved
        Matrix origNm = LagkJointMomentsFromRAP.lagkJointMomentsFromRAP(D0, D1, 3, 1);
        Matrix mapNm = LagkJointMomentsFromMAP.lagkJointMomentsFromMAP(H0, H1, 3, 1, 1e-14);
        for (int i = 0; i < origNm.getNumRows(); i++) {
            for (int j = 0; j < origNm.getNumCols(); j++) {
                assertEquals(origNm.get(i, j), mapNm.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(origNm.get(i, j)) * 1e-3),
                        "Joint moment [" + i + "," + j + "] mismatch");
            }
        }
    }

    /**
     * MMAPFromMRAP.txt: Convert MRAP to MMAP
     * Verify: moments and joint moments preserved
     */
    @Test
    public void testMMAPFromMRAP() {
        MatrixCell H = new MatrixCell();
        H.set(0, createMrapH0());
        H.set(1, createMrapH1());
        H.set(2, createMrapH2());

        // Compute original moments
        double[] origMoms = MarginalMomentsFromMRAP.marginalMomentsFromMRAP(H, 7);
        MatrixCell origJmom = LagkJointMomentsFromMRAP.lagkJointMomentsFromMRAP(H, 3, 1);

        // Convert
        MatrixCell G = MMAPFromMRAP.mmapFromMRAP(H);
        assertNotNull(G);

        // Verify valid MMAP
        assertTrue(CheckMMAPRepresentation.checkMMAPRepresentation(G, 1e-10),
                "MMAPFromMRAP result should be valid MMAP");

        // Verify moments preserved
        double[] mmapMoms = MarginalMomentsFromMMAP.marginalMomentsFromMMAP(G, 7);
        for (int i = 0; i < origMoms.length; i++) {
            assertEquals(origMoms[i], mmapMoms[i], Math.max(LOOSE_MID_TOL, Math.abs(origMoms[i]) * 1e-3),
                    "MMAP moms[" + i + "]");
        }

        // Verify joint moments preserved
        MatrixCell mmapJmom = LagkJointMomentsFromMMAP.lagkJointMomentsFromMMAP(G, 3, 1);
        for (int t = 0; t < 2; t++) {
            Matrix origNm = origJmom.get(t);
            Matrix mmapNm = mmapJmom.get(t);
            for (int i = 0; i < origNm.getNumRows(); i++) {
                for (int j = 0; j < origNm.getNumCols(); j++) {
                    assertEquals(origNm.get(i, j), mmapNm.get(i, j),
                            Math.max(LOOSE_MID_TOL, Math.abs(origNm.get(i, j)) * 1e-3),
                            "JointMom[" + t + "][" + i + "," + j + "]");
                }
            }
        }
    }

    /**
     * RandomMMAP.txt: Generate random MMAP with order=4, types=3, mean=1.62
     * Verify: valid MMAP, correct mean
     */
    @Test
    public void testRandomMMAP() {
        MatrixCell result = RandomMMAP.randomMMAP(4, 3, 1.62, 10, 1000, 1e-7, new Random(42));
        assertNotNull(result);

        // Verify valid MMAP
        assertTrue(CheckMMAPRepresentation.checkMMAPRepresentation(result, 1e-10),
                "RandomMMAP should produce valid MMAP");

        // Verify mean
        double[] moms = MarginalMomentsFromMMAP.marginalMomentsFromMMAP(result, 1);
        assertEquals(1.62, moms[0], COARSE_TOL, "Mean should be close to 1.62");
    }

    /**
     * SamplesFromMMAP.txt: Generate samples from MMAP
     * Input: D0-D3 (2x2 MMAP)
     * Theoretical moments: [1.0007, 2.1045, 6.8277]
     */
    @Test
    public void testSamplesFromMMAP() {
        MatrixCell D = createMmapCell();

        // Generate 10000 samples
        Object result = SamplesFromMMAP.samplesFromMMAP(D, 10000, null, 1e-14, new Random(42));

        // Result should be an array; extract inter-arrival times
        double[] arrivals;
        if (result instanceof double[]) {
            arrivals = (double[]) result;
        } else {
            // Multi-type: double[][] with [time, type] columns
            double[][] multi = (double[][]) result;
            arrivals = new double[multi.length];
            for (int i = 0; i < multi.length; i++) {
                arrivals[i] = multi[i][0];
            }
        }

        // Compute trace moments
        double[] moms = MarginalMomentsFromTrace.marginalMomentsFromTrace(arrivals, 3);

        // Theoretical moments
        double[] theoretical = {1.0007, 2.1045, 6.8277};

        // Stochastic: trace moments should be within 25% of theoretical
        for (int i = 0; i < 3; i++) {
            assertEquals(theoretical[i], moms[i], Math.abs(theoretical[i]) * 0.25,
                    "SamplesFromMMAP moms[" + i + "]");
        }
    }

    // ============ Phase 5: MAP/RAP Fitting Tests ============

    /**
     * MRAPFromMoments.txt: Reconstruct MRAP from moments
     * Input: G0-G3 (3x3 MMAP-like), compute moms+joint moments, reconstruct
     * Verify H{1} output values match expected from butools.tmp data file
     */
    @Test
    public void testMRAPFromMoments() {
        // Build input from CheckMMAPRepresentation test data
        Matrix G0 = new Matrix(3, 3);
        G0.set(0, 0, -1.05); G0.set(0, 1, 0.03); G0.set(0, 2, 0.07);
        G0.set(1, 0, 0.19);  G0.set(1, 1, -1.63); G0.set(1, 2, 0.06);
        G0.set(2, 0, 0.0);   G0.set(2, 1, 0.2);  G0.set(2, 2, -1.03);

        Matrix G1 = new Matrix(3, 3);
        G1.set(0, 0, 0.16); G1.set(0, 1, 0.11); G1.set(0, 2, 0.0);
        G1.set(1, 0, 0.1);  G1.set(1, 1, 0.16); G1.set(1, 2, 0.0);
        G1.set(2, 0, 0.27); G1.set(2, 1, 0.0);  G1.set(2, 2, 0.19);

        Matrix G2 = new Matrix(3, 3);
        G2.set(0, 0, 0.01); G2.set(0, 1, 0.09); G2.set(0, 2, 0.13);
        G2.set(1, 0, 0.26); G2.set(1, 1, 0.21); G2.set(1, 2, 0.05);
        G2.set(2, 0, 0.0);  G2.set(2, 1, 0.16); G2.set(2, 2, 0.07);

        Matrix G3 = new Matrix(3, 3);
        G3.set(0, 0, 0.19); G3.set(0, 1, 0.06); G3.set(0, 2, 0.2);
        G3.set(1, 0, 0.17); G3.set(1, 1, 0.16); G3.set(1, 2, 0.27);
        G3.set(2, 0, 0.0);  G3.set(2, 1, 0.0);  G3.set(2, 2, 0.14);

        MatrixCell G = new MatrixCell();
        G.set(0, G0);
        G.set(1, G1);
        G.set(2, G2);
        G.set(3, G3);

        // Compute moments and joint moments
        double[] moms = MarginalMomentsFromMRAP.marginalMomentsFromMRAP(G, 5);
        assertEquals(0.99805, moms[0], LOOSE_MID_TOL, "moms[0]");

        MatrixCell Nm = LagkJointMomentsFromMRAP.lagkJointMomentsFromMRAP(G, 2, 1);

        // Reconstruct MRAP
        MatrixCell H = MRAPFromMoments.mrapFromMoments(moms, Nm);
        assertNotNull(H);
        assertEquals(4, H.size(), "Expected 4 matrices");

        // Verify H{1} values match expected from data file
        double[][] expectedH1 = {
            {-1.9473, 3.0344, -2.1704},
            {-0.33434, -0.88118, 0.21355},
            {-0.33363, 0.21321, -0.88152}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expectedH1[i][j], H.get(0).get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedH1[i][j]) * 1e-2),
                        "H{1}[" + i + "," + j + "]");
            }
        }

        // Verify D0+D1+...+DK rowsums are zero (basic MRAP property)
        int N = H.get(0).getNumRows();
        for (int i = 0; i < N; i++) {
            double rowSum = 0.0;
            for (int m = 0; m < H.size(); m++) {
                for (int j = 0; j < N; j++) {
                    rowSum += H.get(m).get(i, j);
                }
            }
            assertEquals(0.0, rowSum, FINE_TOL, "row " + i + " sum should be zero");
        }
    }

    /**
     * RAPFromMoments.txt: Reconstruct RAP from moments (test case 1)
     * Input: G0=[-6.2,2,0; 2,-9,1; 1,0,-3], G1=[2.2,-2,4; 2,2,2; 1,0,1]
     * Verify output matrices match expected from data file
     */
    @Test
    public void testRAPFromMoments() {
        Matrix G0 = new Matrix(3, 3);
        G0.set(0, 0, -6.2); G0.set(0, 1, 2.0); G0.set(0, 2, 0.0);
        G0.set(1, 0, 2.0);  G0.set(1, 1, -9.0); G0.set(1, 2, 1.0);
        G0.set(2, 0, 1.0);  G0.set(2, 1, 0.0);  G0.set(2, 2, -3.0);

        Matrix G1 = new Matrix(3, 3);
        G1.set(0, 0, 2.2);  G1.set(0, 1, -2.0); G1.set(0, 2, 4.0);
        G1.set(1, 0, 2.0);  G1.set(1, 1, 2.0);  G1.set(1, 2, 2.0);
        G1.set(2, 0, 1.0);  G1.set(2, 1, 0.0);  G1.set(2, 2, 1.0);

        double[] moms = MarginalMomentsFromMAP.marginalMomentsFromRAP(G0, G1, 5);
        assertEquals(0.36585, moms[0], LOOSE_MID_TOL, "moms[0]");

        Matrix Nm = LagkJointMomentsFromRAP.lagkJointMomentsFromRAP(G0, G1, 2, 1);

        // Reconstruct
        Pair<Matrix, Matrix> result = MRAPFromMoments.rapFromMoments(moms, Nm);
        Matrix H0 = result.getFirst();
        Matrix H1 = result.getSecond();

        // Verify H0 values match expected from data file
        double[][] expectedH0 = {
            {-12.949, 36.78, -24.817},
            {-1.1102, -2.5113, 0.91705},
            {-0.71205, 0.68912, -2.7393}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expectedH0[i][j], H0.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedH0[i][j]) * 1e-2),
                        "H0[" + i + "," + j + "]");
            }
        }

        // Verify H1 values match expected from data file
        double[][] expectedH1 = {
            {9.2672, -99.958, 91.678},
            {1.1693, -2.1771, 3.7123},
            {0.65292, 3.9994, -1.8901}
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expectedH1[i][j], H1.get(i, j),
                        Math.max(LOOSE_MID_TOL, Math.abs(expectedH1[i][j]) * 1e-2),
                        "H1[" + i + "," + j + "]");
            }
        }

        // Verify H0+H1 rowsums are zero
        for (int i = 0; i < 3; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < 3; j++) {
                rowSum += H0.get(i, j) + H1.get(i, j);
            }
            assertEquals(0.0, rowSum, FINE_TOL, "row " + i + " sum should be zero");
        }
    }

    /**
     * RAPFromMomentsAndCorrelations.txt: Reconstruct RAP from moments and correlations
     * Input: H0=[-6.2,2,0; 2,-9,1; 1,0,-3], H1=[2.2,0,2; 0,4,2; 0,1,1]
     * Verify moments and correlations roundtrip
     */
    @Test
    public void testRAPFromMomentsAndCorrelations() {
        Matrix H0 = new Matrix(3, 3);
        H0.set(0, 0, -6.2); H0.set(0, 1, 2.0);  H0.set(0, 2, 0.0);
        H0.set(1, 0, 2.0);  H0.set(1, 1, -9.0);  H0.set(1, 2, 1.0);
        H0.set(2, 0, 1.0);  H0.set(2, 1, 0.0);   H0.set(2, 2, -3.0);

        Matrix H1 = new Matrix(3, 3);
        H1.set(0, 0, 2.2);  H1.set(0, 1, 0.0);   H1.set(0, 2, 2.0);
        H1.set(1, 0, 0.0);  H1.set(1, 1, 4.0);   H1.set(1, 2, 2.0);
        H1.set(2, 0, 0.0);  H1.set(2, 1, 1.0);   H1.set(2, 2, 1.0);

        double[] mom = MarginalMomentsFromMAP.marginalMomentsFromRAP(H0, H1, 5);
        assertEquals(0.29774, mom[0], LOOSE_MID_TOL, "mom[0]");

        double[] corr = LagCorrelationsFromMAP.lagCorrelationsFromRAP(H0, H1, 3);
        assertEquals(0.012394, corr[0], LOOSE_MID_TOL, "corr[0]");

        // Reconstruct
        Pair<Matrix, Matrix> result = RAPFromMomentsAndCorrelations.rapFromMomentsAndCorrelations(mom, corr);
        Matrix G0 = result.getFirst();
        Matrix G1 = result.getSecond();

        // Verify moments roundtrip
        double[] rmom = MarginalMomentsFromMAP.marginalMomentsFromRAP(G0, G1, 5);
        for (int i = 0; i < 5; i++) {
            assertEquals(mom[i], rmom[i], Math.max(LOOSE_MID_TOL, Math.abs(mom[i]) * 1e-3),
                    "roundtrip moms[" + i + "]");
        }

        // Verify correlations roundtrip
        double[] rcorr = LagCorrelationsFromMAP.lagCorrelationsFromRAP(G0, G1, 3);
        for (int i = 0; i < 3; i++) {
            assertEquals(corr[i], rcorr[i], LOOSE_MID_TOL, "roundtrip corr[" + i + "]");
        }
    }

    /**
     * MAPFromFewMomentsAndCorrelations.txt Test 1: 2 moments, negative correlation
     * Input: moms=[1.1, 6.05], corr1=-0.17
     * Expected: 4x4 MAP, moments and correlation roundtrip
     */
    @Test
    public void testMAPFromFewMomentsAndCorrelations_2Mom_NegCorr() {
        double[] moms = {1.1, 6.05};
        double corr1 = -0.17;

        Matrix[] result = MAPFromFewMomentsAndCorrelations.mapFromFewMomentsAndCorrelations(moms, corr1, null);
        Matrix D0 = result[0];
        Matrix D1 = result[1];

        // Verify valid MAP
        assertTrue(CheckMAPRepresentation.checkMAPRepresentation(D0, D1, 1e-10),
                "Result should be valid MAP");

        // Verify moments roundtrip
        double[] rmoms = MarginalMomentsFromMAP.marginalMomentsFromMAP(D0, D1, 2);
        assertEquals(1.1, rmoms[0], LOOSE_MID_TOL, "roundtrip moms[0]");
        assertEquals(6.05, rmoms[1], LOOSE_MID_TOL, "roundtrip moms[1]");

        // Verify correlation roundtrip
        double[] rcorr = LagCorrelationsFromMAP.lagCorrelationsFromMAP(D0, D1, 1);
        assertEquals(-0.17, rcorr[0], LOOSE_MID_TOL, "roundtrip corr1");
    }

    /**
     * MAPFromFewMomentsAndCorrelations.txt Test 2: 3 moments, negative correlation
     * Input: moms=[1.2, 4.32, 20.0], corr1=-0.4
     * Expected: 6x6 MAP, moments and correlation roundtrip
     */
    @Test
    public void testMAPFromFewMomentsAndCorrelations_3Mom_NegCorr() {
        double[] moms = {1.2, 4.32, 20.0};
        double corr1 = -0.4;

        Matrix[] result = MAPFromFewMomentsAndCorrelations.mapFromFewMomentsAndCorrelations(moms, corr1, null);
        Matrix D0 = result[0];
        Matrix D1 = result[1];

        // Verify valid MAP
        assertTrue(CheckMAPRepresentation.checkMAPRepresentation(D0, D1, 1e-10),
                "Result should be valid MAP");

        // Verify moments roundtrip
        double[] rmoms = MarginalMomentsFromMAP.marginalMomentsFromMAP(D0, D1, 3);
        assertEquals(1.2, rmoms[0], LOOSE_MID_TOL, "roundtrip moms[0]");
        assertEquals(4.32, rmoms[1], LOOSE_MID_TOL, "roundtrip moms[1]");
        assertEquals(20.0, rmoms[2], LOOSE_MID_TOL, "roundtrip moms[2]");

        // Verify correlation roundtrip
        double[] rcorr = LagCorrelationsFromMAP.lagCorrelationsFromMAP(D0, D1, 1);
        assertEquals(-0.4, rcorr[0], LOOSE_MID_TOL, "roundtrip corr1");
    }

    /**
     * MAPFromFewMomentsAndCorrelations.txt Test 3: 3 moments, positive correlation
     * Input: moms=[1.2, 4.32, 20.0], corr1=0.4
     * Expected: 9x9 MAP, moments and correlation roundtrip
     */
    @Test
    public void testMAPFromFewMomentsAndCorrelations_3Mom_PosCorr() {
        double[] moms = {1.2, 4.32, 20.0};
        double corr1 = 0.4;

        Matrix[] result = MAPFromFewMomentsAndCorrelations.mapFromFewMomentsAndCorrelations(moms, corr1, null);
        Matrix D0 = result[0];
        Matrix D1 = result[1];

        // Verify valid MAP
        assertTrue(CheckMAPRepresentation.checkMAPRepresentation(D0, D1, 1e-10),
                "Result should be valid MAP");

        // Verify moments roundtrip
        double[] rmoms = MarginalMomentsFromMAP.marginalMomentsFromMAP(D0, D1, 3);
        assertEquals(1.2, rmoms[0], LOOSE_MID_TOL, "roundtrip moms[0]");
        assertEquals(4.32, rmoms[1], LOOSE_MID_TOL, "roundtrip moms[1]");
        assertEquals(20.0, rmoms[2], LOOSE_MID_TOL, "roundtrip moms[2]");

        // Verify correlation roundtrip
        double[] rcorr = LagCorrelationsFromMAP.lagCorrelationsFromMAP(D0, D1, 1);
        assertEquals(0.4, rcorr[0], LOOSE_MID_TOL, "roundtrip corr1");
    }
}
