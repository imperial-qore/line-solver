/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sn;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Tests whether a process representation admits a phase-type reading.
 *
 * A representation is Markovian when D0 has nonnegative off-diagonal entries, every D_k
 * with k &gt;= 1 is nonnegative, and the entry vector pie is nonnegative. Exactly under
 * those conditions do sn.mu, sn.phi and sn.pie carry their probabilistic reading
 * (mu_i = -D0(i,i) is a rate, phi_i is a completion probability, pie is a distribution
 * over phases), which is what the CTMC state space, SSA and the fluid ODEs consume.
 *
 * A matrix-exponential (ME) or rational (RAP) process fails the test: its moments,
 * transforms and aggregated stationary measures remain exact, but the per-phase
 * quantities are signed. See _kb/04-networkstruct.md.
 *
 * Mirrors matlab/src/api/sn/sn_is_phasetype.m and the native Python
 * line_solver.api.sn.sn_is_phasetype.
 */
public final class SnIsPhaseType {

    private SnIsPhaseType() {}

    /**
     * Tests a (D0, D1, ...) representation without an entry vector.
     *
     * @param proc the process representation
     * @return true when the representation is Markovian
     */
    public static boolean snIsPhaseType(MatrixCell proc) {
        return snIsPhaseType(proc, null);
    }

    /**
     * Tests a (D0, D1, ...) representation and its entry vector.
     *
     * @param proc the process representation
     * @param pie  the entry vector, or null to skip that test
     * @return true when the representation is Markovian
     */
    public static boolean snIsPhaseType(MatrixCell proc, Matrix pie) {
        double tol = GlobalConstants.Zero;

        // An empty, scalar-parameter or NaN-carrying entry describes a disabled or
        // not-yet-Markovian process; there is no phase decomposition to invalidate, so
        // the caller is not blocked by this test.
        if (proc == null || proc.size() < 2) {
            return true;
        }
        Matrix D0 = proc.get(0);
        if (D0 == null || D0.isEmpty() || D0.hasNaN()) {
            return true;
        }
        int n = D0.getNumRows();
        if (D0.getNumCols() != n) {
            return true;
        }

        // Off-diagonal entries of D0 are transition rates between phases.
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j && D0.get(i, j) < -tol) {
                    return false;
                }
            }
        }

        // D1 and any further D_k are jump matrices and must be nonnegative.
        for (int k = 1; k < proc.size(); k++) {
            Matrix Dk = proc.get(k);
            if (Dk == null || Dk.hasNaN()) {
                continue;
            }
            for (int i = 0; i < Dk.getNumRows(); i++) {
                for (int j = 0; j < Dk.getNumCols(); j++) {
                    if (Dk.get(i, j) < -tol) {
                        return false;
                    }
                }
            }
        }

        if (pie != null && !pie.isEmpty() && !pie.hasNaN()) {
            for (int i = 0; i < pie.getNumElements(); i++) {
                if (pie.get(i) < -tol) {
                    return false;
                }
            }
        }

        return true;
    }
}
