/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.aoi;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Analytic validation of the Age-of-Information APIs (jline.api.aoi):
 * - the M/M/1 FCFS mean AoI closed form (Inoue et al.) is recomputed inline;
 * - LST-based general routines must reduce to the specialized closed forms
 *   when instantiated with the matching distribution;
 * - LCFS-PR dominance over FCFS for M/M/1;
 * - the aoi_lst_* factory transforms are valid LSTs.
 */
public class AoiApiTest {

    private static final double TOL = 1e-9;
    private static final double NUM_TOL = 1e-6;

    private static final double LAMBDA = 1.0;
    private static final double MU = 2.0;
    private static final double RHO = LAMBDA / MU;

    @Test
    public void fcfsMm1MatchesClosedForm() {
        AoiResult r = Aoi_fcfs_mm1.aoi_fcfs_mm1(LAMBDA, MU);
        double meanExact = (1.0 / MU) * (1.0 + 1.0 / RHO + RHO * RHO / (1.0 - RHO));
        double peakExact = (1.0 / MU) * (1.0 + 1.0 / RHO + RHO / (1.0 - RHO));
        assertEquals(meanExact, r.getMeanAoI(), TOL, "M/M/1 FCFS mean AoI");
        assertEquals(peakExact, r.getPeakAoI(), TOL, "M/M/1 FCFS peak AoI");
    }

    @Test
    public void fcfsMgi1ReducesToMm1ForExponentialService() {
        AoiResult exact = Aoi_fcfs_mm1.aoi_fcfs_mm1(LAMBDA, MU);
        LstFunction expLst = Aoi_lst.aoi_lst_exp(MU);
        AoiLstResult general = Aoi_fcfs_mgi1.aoi_fcfs_mgi1(
                LAMBDA, expLst, 1.0 / MU, 2.0 / (MU * MU));
        assertEquals(exact.getMeanAoI(), general.getMeanAoI(), NUM_TOL,
                "M/GI/1 with exponential H must equal M/M/1 mean AoI");
        assertEquals(exact.getPeakAoI(), general.getPeakAoI(), NUM_TOL,
                "M/GI/1 with exponential H must equal M/M/1 peak AoI");
    }

    @Test
    public void fcfsGim1ReducesToMm1ForExponentialArrivals() {
        AoiResult exact = Aoi_fcfs_mm1.aoi_fcfs_mm1(LAMBDA, MU);
        LstFunction expLst = Aoi_lst.aoi_lst_exp(LAMBDA);
        AoiLstResult general = Aoi_fcfs_gim1.aoi_fcfs_gim1(
                expLst, MU, 1.0 / LAMBDA, 2.0 / (LAMBDA * LAMBDA));
        assertEquals(exact.getMeanAoI(), general.getMeanAoI(), NUM_TOL,
                "GI/M/1 with exponential Y must equal M/M/1 mean AoI");
    }

    @Test
    public void fcfsMd1ReducesToMgi1WithDeterministicLst() {
        double d = 0.4;
        AoiResult md1 = Aoi_fcfs_md1.aoi_fcfs_md1(LAMBDA, d);
        LstFunction detLst = Aoi_lst.aoi_lst_det(d);
        AoiLstResult general = Aoi_fcfs_mgi1.aoi_fcfs_mgi1(LAMBDA, detLst, d, d * d);
        assertEquals(md1.getMeanAoI(), general.getMeanAoI(), NUM_TOL,
                "M/D/1 closed form vs M/GI/1 with deterministic LST");
    }

    @Test
    public void lcfsprMm1ReducesFromGeneralRoutine() {
        AoiResult exact = Aoi_lcfspr_mm1.aoi_lcfspr_mm1(LAMBDA, MU);
        LstFunction expLst = Aoi_lst.aoi_lst_exp(MU);
        AoiLstResult general = Aoi_lcfspr_mgi1.aoi_lcfspr_mgi1(
                LAMBDA, expLst, 1.0 / MU, 2.0 / (MU * MU));
        assertEquals(exact.getMeanAoI(), general.getMeanAoI(), NUM_TOL,
                "LCFS-PR M/GI/1 with exponential H must equal LCFS-PR M/M/1");
    }

    @Test
    public void lcfsprDominatesFcfsInMm1() {
        double fcfs = Aoi_fcfs_mm1.aoi_fcfs_mm1(LAMBDA, MU).getMeanAoI();
        double lcfspr = Aoi_lcfspr_mm1.aoi_lcfspr_mm1(LAMBDA, MU).getMeanAoI();
        assertTrue(lcfspr <= fcfs + TOL,
                "LCFS-PR mean AoI must not exceed FCFS mean AoI in M/M/1");
    }

    @Test
    public void fcfsDm1ReducesToGim1WithDeterministicLst() {
        double tau = 1.0; // deterministic interarrival
        AoiResult dm1 = Aoi_fcfs_dm1.aoi_fcfs_dm1(tau, MU);
        LstFunction detLst = Aoi_lst.aoi_lst_det(tau);
        AoiLstResult general = Aoi_fcfs_gim1.aoi_fcfs_gim1(detLst, MU, tau, tau * tau);
        assertEquals(dm1.getMeanAoI(), general.getMeanAoI(), NUM_TOL,
                "D/M/1 closed form vs GI/M/1 with deterministic LST");
    }

    @Test
    public void lstFactoriesProduceValidTransforms() {
        LstFunction[] lsts = {
                Aoi_lst.aoi_lst_exp(2.0),
                Aoi_lst.aoi_lst_det(0.5),
                Aoi_lst.aoi_lst_erlang(2, 4.0)
        };
        for (LstFunction lst : lsts) {
            assertEquals(1.0, lst.evaluate(0.0), 1e-9, "LST(0) must equal 1");
            double prev = 1.0;
            for (double s : new double[]{0.1, 0.5, 1.0, 2.0}) {
                double v = lst.evaluate(s);
                assertTrue(v > 0 && v <= 1.0 + 1e-12, "LST value outside (0,1]");
                assertTrue(v <= prev + 1e-12, "LST must be nonincreasing");
                prev = v;
            }
        }
    }

    @Test
    public void lcfsVariantsAreOrderedByPreemptionDegree() {
        // For M/M/1, discarding (lcfsd) and serving (lcfss) variants are valid
        // sampling disciplines: their mean AoI is finite and positive, and the
        // preemptive variant is the best of the three LCFS policies.
        LstFunction expLst = Aoi_lst.aoi_lst_exp(MU);
        double eH = 1.0 / MU, eH2 = 2.0 / (MU * MU);
        double lcfspr = Aoi_lcfspr_mgi1.aoi_lcfspr_mgi1(LAMBDA, expLst, eH, eH2).getMeanAoI();
        double lcfsd = Aoi_lcfsd_mgi1.aoi_lcfsd_mgi1(LAMBDA, expLst, eH, eH2).getMeanAoI();
        double lcfss = Aoi_lcfss_mgi1.aoi_lcfss_mgi1(LAMBDA, expLst, eH, eH2).getMeanAoI();
        assertTrue(Double.isFinite(lcfsd) && lcfsd > 0, "LCFS-D mean AoI invalid");
        assertTrue(Double.isFinite(lcfss) && lcfss > 0, "LCFS-S mean AoI invalid");
        assertTrue(lcfspr <= lcfsd + TOL, "LCFS-PR must dominate LCFS-D in M/M/1");
        assertTrue(lcfspr <= lcfss + TOL, "LCFS-PR must dominate LCFS-S in M/M/1");
    }

    private static jline.util.matrix.Matrix scalar(double v) {
        jline.util.matrix.Matrix m = new jline.util.matrix.Matrix(1, 1);
        m.set(0, 0, v);
        return m;
    }

    @Test
    public void singleBufferMfqMatchesMatlabReference() {
        // Matrix-fluid-queue AoI solver, single-buffer, exponential service
        // (sigma=[1], S=[-2]); reference values from MATLAB solveSingleBuffer.
        AoiMfqResult r = Aoi_solve_singlebuffer.aoi_solve_singlebuffer(
                1.0, scalar(1.0), scalar(-2.0), 0.5);
        assertEquals(1.609524, r.getAoiMean(), 1e-5, "single-buffer mean AoI");
        assertEquals(1.251020, r.getAoiVar(), 1e-5, "single-buffer AoI variance");
    }

    @Test
    public void bufferlessMfqMatchesMatlabReference() {
        // Matrix-fluid-queue AoI solver, bufferless; reference from MATLAB
        // solveBufferless with exponential arrival and service.
        AoiMfqResult r = Aoi_solve_bufferless.aoi_solve_bufferless(
                scalar(1.0), scalar(-1.0), scalar(1.0), scalar(-2.0), 0.5);
        assertEquals(1.566667, r.getAoiMean(), 1e-5, "bufferless mean AoI");
        assertEquals(1.298889, r.getAoiVar(), 1e-5, "bufferless AoI variance");
    }
}
