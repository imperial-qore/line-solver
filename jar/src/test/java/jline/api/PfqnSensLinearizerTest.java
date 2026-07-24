/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import java.util.Random;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.mva.Pfqn_linearizer;
import jline.api.pfqn.sens.Pfqn_sens_linearizer;
import jline.api.pfqn.sens.Pfqn_sens_mom;
import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates {@link Pfqn_sens_linearizer#pfqn_sens_linearizer}, the LINEARIZER-2 /
 * LINEARIZER-3 moment approximation of Strelen (1990), Section 5, mirroring the
 * MATLAB harness {@code pfqn_sens_linearizer_validate.m}.
 *
 * <p>This routine is an APPROXIMATION, so it must not be held to machine
 * precision. The checks are therefore of three kinds:</p>
 *
 * <ul>
 *   <li>A. accuracy bands against the exact {@link Pfqn_sens_mom}, on models small
 *       enough for the exact lattice. The reference reports relative errors below
 *       2.1% on E[Q], 4.1% on E[Q^2] and 6.2% on E[Q^3] over its own 51 networks;
 *       the bands asserted here are of that order. This is the only meaningful
 *       statement of correctness for an approximation: it must track the exact
 *       answer, not equal it;</li>
 *   <li>B. exactness where the approximation degenerates. At a population of one
 *       job the CORE estimate of the queue lengths one job down is identically zero
 *       whatever the delta terms are, so the Linearizer equations coincide with the
 *       exact MVA and every moment must match {@link Pfqn_sens_mom} to roundoff.
 *       This pins the derivative algebra independently of the heuristic;</li>
 *   <li>C. structural invariants that hold for any population: the mean queue
 *       lengths conserve the population, and they agree with LINE's own
 *       {@link Pfqn_linearizer}, which runs the same heuristic without
 *       derivatives.</li>
 * </ul>
 *
 * <p>Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks and
 * its Linearizer", Performance Evaluation 11:127-142, 1990.</p>
 */
public class PfqnSensLinearizerTest {

    /** B: the one-job case must be exact */
    private static final double TOL_EXACT = 1e-9;
    /** C: population conservation */
    private static final double TOL_POP = 1e-8;
    /** C: agreement with LINE's pfqn_linearizer on the means */
    private static final double TOL_LIN = 5e-2;

    // A: the bands are the accuracy the reference itself claims in Section 5 over
    // its 51 networks, so this asserts that our port reproduces the paper's own
    // accuracy statement rather than some slacker figure. The seed is fixed, so the
    // model set is deterministic and these are hard regression guards.
    private static final double BAND_M = 0.021;    // r(Q)   < 2.1% in the reference
    private static final double BAND_M2 = 0.041;   // r(Q^2) < 4.1% in the reference
    private static final double BAND_M3 = 0.062;   // r(Q^3) < 6.2% in the reference
    // The reference does not report an error on the variance. It is naturally
    // larger than the one on E[Q^2] because Var = E[Q^2] - E[Q]^2 is a difference of
    // larger numbers, so relative error is amplified; banded here only to catch
    // regressions.
    private static final double BAND_VAR = 0.08;

    /** max_k |a_k - b_k| / max(1, |a_k|, |b_k|), matching the MATLAB relerr. */
    private static double relerr(double[] a, double[] b) {
        double e = 0.0;
        for (int k = 0; k < a.length; k++) {
            double scale = Math.max(1.0, Math.max(Math.abs(a[k]), Math.abs(b[k])));
            e = Math.max(e, Math.abs(a[k] - b[k]) / scale);
        }
        return e;
    }

    private static double relerr(Matrix a, Matrix b) {
        return relerr(flat(a), flat(b));
    }

    private static double[] flat(Matrix m) {
        double[] v = new double[m.getNumRows() * m.getNumCols()];
        int k = 0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                v[k++] = m.get(i, j);
            }
        }
        return v;
    }

    private static double[] rowSums(Matrix m) {
        double[] v = new double[m.getNumRows()];
        for (int i = 0; i < m.getNumRows(); i++) {
            double s = 0.0;
            for (int j = 0; j < m.getNumCols(); j++) {
                s += m.get(i, j);
            }
            v[i] = s;
        }
        return v;
    }

    /** A / C. random closed models against the exact moment analysis. */
    @Test
    public void sensLinearizerTracksExactMomentsWithinThePapersBands() {
        Random rng = new Random(7);
        double errPop = 0.0;
        double errLin = 0.0;
        double worstM = 0.0;
        double worstM2 = 0.0;
        double worstM3 = 0.0;
        double worstVar = 0.0;
        int nA = 0;

        for (int trial = 1; trial <= 60; trial++) {
            int M = 2 + rng.nextInt(3);
            int R = 1 + rng.nextInt(2);
            Matrix L = new Matrix(M, R);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    L.set(i, r, 0.2 + rng.nextDouble());
                }
            }
            Matrix N = new Matrix(1, R);
            for (int r = 0; r < R; r++) {
                N.set(0, r, 1 + rng.nextInt(4));
            }
            Matrix Z = new Matrix(1, R);
            if (trial % 2 == 0) {
                for (int r = 0; r < R; r++) {
                    Z.set(0, r, 0.5 + rng.nextDouble());
                }
            }

            Ret.pfqnSensLinearizer app = Pfqn_sens_linearizer.pfqn_sens_linearizer(L, N, Z);
            Ret.pfqnSensMom ex = Pfqn_sens_mom.pfqn_sens_mom(L, N, Z);

            worstM = Math.max(worstM, relerr(app.m, ex.m));
            worstM2 = Math.max(worstM2, relerr(app.M2, ex.M2));
            worstM3 = Math.max(worstM3, relerr(app.M3, ex.M3));
            worstVar = Math.max(worstVar, relerr(app.Var, ex.Var));
            nA++;

            // ---- C. population conservation --------------------------------
            for (int r = 0; r < R; r++) {
                double inNet = 0.0;
                for (int i = 0; i < M; i++) {
                    inNet += app.Q.get(i, r);
                }
                double inDelay = app.X.get(0, r) * Z.get(0, r);
                errPop = Math.max(errPop, Math.abs(inNet + inDelay - N.get(0, r))
                        / Math.max(1.0, N.get(0, r)));
            }

            // ---- C. means against LINE's own Linearizer --------------------
            SchedStrategy[] type = new SchedStrategy[M];
            for (int i = 0; i < M; i++) {
                type[i] = SchedStrategy.PS;
            }
            Ret.pfqnAMVA lin = Pfqn_linearizer.pfqn_linearizer(L, N, Z, type, 1e-10, 500);
            errLin = Math.max(errLin, relerr(rowSums(app.Q), rowSums(lin.Q)));
        }

        assertTrue(nA > 0, "no model was checked against the exact moment analysis");
        assertTrue(worstM <= BAND_M,
                "E[Q] accuracy vs exact pfqn_sens_mom (" + nA + " models): " + (100 * worstM)
                        + "% (band " + (100 * BAND_M) + "%)");
        assertTrue(worstVar <= BAND_VAR,
                "Var[Q] accuracy vs exact pfqn_sens_mom: " + (100 * worstVar)
                        + "% (band " + (100 * BAND_VAR) + "%)");
        assertTrue(worstM2 <= BAND_M2,
                "E[Q^2] accuracy vs exact pfqn_sens_mom: " + (100 * worstM2)
                        + "% (band " + (100 * BAND_M2) + "%)");
        assertTrue(worstM3 <= BAND_M3,
                "E[Q^3] accuracy vs exact pfqn_sens_mom: " + (100 * worstM3)
                        + "% (band " + (100 * BAND_M3) + "%)");
        assertTrue(errPop <= TOL_POP, "population conservation: " + errPop);
        assertTrue(errLin <= TOL_LIN, "means vs pfqn_linearizer: " + errLin);
    }

    /** B. one job: the Linearizer equations degenerate to the exact MVA. */
    @Test
    public void sensLinearizerIsExactWithASingleJob() {
        Random rng = new Random(70);
        double errExact = 0.0;
        for (int trial = 1; trial <= 12; trial++) {
            int M = 2 + rng.nextInt(3);
            Matrix L = new Matrix(M, 1);
            for (int i = 0; i < M; i++) {
                L.set(i, 0, 0.2 + rng.nextDouble());
            }
            Matrix Z = new Matrix(1, 1);
            if (trial % 2 == 0) {
                Z.set(0, 0, 0.4 + rng.nextDouble());
            }
            Matrix N = new Matrix(1, 1);
            N.set(0, 0, 1);

            Ret.pfqnSensLinearizer app = Pfqn_sens_linearizer.pfqn_sens_linearizer(L, N, Z);
            Ret.pfqnSensMom ex = Pfqn_sens_mom.pfqn_sens_mom(L, N, Z);
            errExact = Math.max(errExact, relerr(app.m, ex.m));
            errExact = Math.max(errExact, relerr(app.Var, ex.Var));
            errExact = Math.max(errExact, relerr(app.M2, ex.M2));
            errExact = Math.max(errExact, relerr(app.M3, ex.M3));
            errExact = Math.max(errExact, relerr(app.Cov, ex.Cov));
        }
        assertTrue(errExact <= TOL_EXACT, "one job, exact vs pfqn_sens_mom: " + errExact);
    }
}
