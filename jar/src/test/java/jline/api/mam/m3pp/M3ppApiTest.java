/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.mam.m3pp;

import jline.api.mam.Map_isfeasible;
import jline.api.mam.Mmap_count_mean;
import jline.api.mam.Mmap_lambda;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the marked-MMPP (M3PP) generators and fitters
 * (jline.api.mam.m3pp) via structural and fit-consistency invariants:
 * - random M3PPs must be feasible marked MAPs whose class matrices sum to
 *   the aggregate D1;
 * - count fitting must preserve the per-class rates it was given;
 * - the counting process mean must grow as lambda*t.
 */
public class M3ppApiTest {

    private static final double TOL = 1e-8;

    private static void assertValidMmap(MatrixCell mmap, String context) {
        assertNotNull(mmap, context + ": null MMAP");
        assertTrue(mmap.size() >= 3, context + ": marked MAP needs D0, D1 and class matrices");
        Matrix d1 = mmap.get(1);
        int n = d1.getNumRows();
        // Class matrices are nonnegative and sum to the aggregate D1
        Matrix sum = new Matrix(n, n);
        for (int k = 2; k < mmap.size(); k++) {
            Matrix dk = mmap.get(k);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    assertTrue(dk.get(i, j) >= -TOL,
                            context + ": class matrix D1" + (k - 1) + " has negative entry");
                    sum.set(i, j, sum.get(i, j) + dk.get(i, j));
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                assertEquals(d1.get(i, j), sum.get(i, j), 1e-8,
                        context + ": class matrices must sum to D1 at (" + i + "," + j + ")");
            }
        }
        // The aggregate (D0, D1) must be a feasible MAP
        assertTrue(Map_isfeasible.map_isfeasible(new MatrixCell(mmap.get(0), mmap.get(1))),
                context + ": aggregate MAP infeasible");
    }

    @Test
    public void randomM3ppIsAFeasibleMarkedMap() {
        MatrixCell mmap = M3pp_rand.m3pp_rand(2, 2, 23000L);
        assertValidMmap(mmap, "m3pp_rand(2,2)");
        Matrix lambda = Mmap_lambda.mmap_lambda(mmap);
        for (int k = 0; k < lambda.getNumElements(); k++) {
            assertTrue(lambda.get(k) > 0, "class rate " + k + " must be positive");
        }
    }

    // Reference M3PP(2,2): D0 = [-3.5 0.5; 0.2 -1.2], D11 = diag(2.0, 0.6),
    // D12 = diag(1.0, 0.4). The statistics below are its exact counting
    // characteristics at t1 = 1, t2 = t3 = 10, taken from MATLAB; the fit must
    // reproduce the process, which is what the reference implementation does.
    private static final double REF_A = 1.5714285714285716;
    private static final double REF_BT1 = 1.4168254519828447;
    private static final double REF_BT2 = 2.2723905395103197;
    private static final double REF_BINF = 2.4840180227935313;
    private static final double REF_M3T2 = 136.99785175249326;
    private static final double REF_T1 = 1.0;
    private static final double REF_T2 = 10.0;
    private static final double REF_T3 = 10.0;
    private static final double[] REF_AI = {1.0, 0.57142857142857151};
    private static final double[] REF_DVT3 = {3.8858201161539156, -20.681375237690119};

    private static MatrixCell fitReference() {
        Matrix[] fitted = M3pp2m_fitc.m3pp2m_fitc(REF_A, REF_BT1, REF_BT2, REF_BINF, REF_M3T2,
                REF_T1, REF_T2, REF_AI, REF_DVT3, REF_T3);
        MatrixCell cell = new MatrixCell(fitted.length);
        for (int i = 0; i < fitted.length; i++) {
            cell.set(i, fitted[i]);
        }
        return cell;
    }

    @Test
    public void fitcPreservesPerClassRates() {
        MatrixCell fitted = fitReference();
        assertValidMmap(fitted, "m3pp2m_fitc");
        Matrix lambda = Mmap_lambda.mmap_lambda(fitted);
        double aggregate = 0;
        for (int k = 0; k < lambda.getNumElements(); k++) {
            aggregate += lambda.get(k);
        }
        assertEquals(REF_A, aggregate, 1e-6, "fitted aggregate rate must match the target rate");
        assertEquals(REF_AI[0], lambda.get(0), 1e-6, "fitted class 1 rate");
        assertEquals(REF_AI[1], lambda.get(1), 1e-6, "fitted class 2 rate");
    }

    @Test
    public void fitcMatchesTheReferenceImplementation() {
        // Pinned against MATLAB m3pp2m_fitc.m; native Python agrees to 8e-12.
        // Guards the closed form against silent drift in any of the codebases.
        MatrixCell fitted = fitReference();
        double[][] expectedD0 = {{-3.5015668436003242, 0.50033168008885631},
                                 {0.19979451688103933, -1.2002668033751962}};
        double[][] expectedD11 = {{1.7006036136346727, 0.0}, {0.0, 0.72023206588399402}};
        double[][] expectedD12 = {{1.3006315498767951, 0.0}, {0.0, 0.28024022061016285}};
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(expectedD0[i][j], fitted.get(0).get(i, j), 1e-9, "D0(" + i + "," + j + ")");
                assertEquals(expectedD11[i][j], fitted.get(2).get(i, j), 1e-9, "D11(" + i + "," + j + ")");
                assertEquals(expectedD12[i][j], fitted.get(3).get(i, j), 1e-9, "D12(" + i + "," + j + ")");
            }
        }
    }

    @Test
    public void countMeanGrowsLinearlyWithRate() {
        MatrixCell fitted = fitReference();
        double t = 10.0;
        Matrix mean = Mmap_count_mean.mmap_count_mean(fitted, t);
        assertEquals(REF_AI[0] * t, mean.get(0), 1e-5, "class 1 count mean must be lambda1*t");
        assertEquals(REF_AI[1] * t, mean.get(1), 1e-5, "class 2 count mean must be lambda2*t");
    }

    @org.junit.jupiter.api.Test
    public void poissonCountProcessHasVarianceEqualToMean() {
        // For a Poisson process the count over [0,t] has Var[N(t)] = E[N(t)]
        // = lambda*t and index of dispersion 1.
        Matrix d0 = new Matrix(1, 1);
        d0.set(0, 0, -2.0);
        Matrix d1 = new Matrix(1, 1);
        d1.set(0, 0, 2.0);
        MatrixCell poisson = new MatrixCell(new Matrix[]{d0, d1, d1});
        double t = 5.0;
        double mean = jline.api.mam.Mmap_count_mean.mmap_count_mean(poisson, t).get(0);
        double var = jline.api.mam.Mmap_count_var.mmap_count_var(poisson, t).get(0);
        assertEquals(10.0, mean, 1e-6, "Poisson(2) count mean over t=5");
        assertEquals(10.0, var, 1e-6, "Poisson(2) count variance equals the mean");
        assertEquals(1.0, jline.api.mam.Mmap_count_idc.mmap_count_idc(poisson, t).get(0),
                1e-6, "Poisson index of dispersion is 1");
    }
}
