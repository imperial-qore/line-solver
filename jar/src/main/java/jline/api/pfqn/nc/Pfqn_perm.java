/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.nc;

import jline.lib.perm.Permanent;
import jline.util.matrix.Matrix;

/**
 * Permanent of a demand matrix, with optional column multiplicities.
 *
 * <p>The pfqn_ entry point of the permanent library. It exists because the
 * product-form joint queue-length probability of the per-station TOTAL
 * populations is a permanent of the demand matrix replicated once per job,
 * which is a normalizing-constant quantity rather than a general-purpose linear
 * algebra one; see {@link jline.api.pfqn.Pfqn_jointmarg}.</p>
 *
 * <p>Orientation is chosen before repeated lines are grouped, inside
 * {@link Permanent}. That is a correctness concern, not an optimisation:
 * perm(A) is transpose-invariant but the Ryser SUM is not, and exploiting
 * repeated rows silently expands the transpose.</p>
 *
 * <p>Reference: H. J. Ryser, "Combinatorial Mathematics", Carus Mathematical
 * Monographs 14, Mathematical Association of America, 1963.</p>
 */
public final class Pfqn_perm {
    private Pfqn_perm() {}

    /**
     * Permanent of the square matrix A.
     *
     * @param A square matrix
     * @return the permanent
     */
    public static double pfqn_perm(Matrix A) {
        if (A == null || A.getNumRows() == 0) {
            // The permanent of the empty matrix is 1, which is what makes an
            // all-empty occupancy vector free of any special case upstream.
            return 1.0;
        }
        if (A.getNumRows() != A.getNumCols()) {
            throw new IllegalArgumentException("pfqn_perm: the matrix must be square, got "
                    + A.getNumRows() + "x" + A.getNumCols());
        }
        return new Permanent(A, true).value;
    }

    /**
     * Permanent of the matrix whose column j is column j of A repeated m[j]
     * times, so that sum(m) equals the number of rows of A.
     *
     * @param A matrix of distinct columns
     * @param m multiplicity of each column
     * @return the permanent of the expanded matrix
     */
    public static double pfqn_perm(Matrix A, int[] m) {
        if (A == null || m == null) {
            return 1.0;
        }
        if (m.length != A.getNumCols()) {
            throw new IllegalArgumentException("pfqn_perm: the multiplicity vector has " + m.length
                    + " entries but the matrix has " + A.getNumCols() + " columns");
        }
        // Materialize the expansion so that both entry points share the
        // orientation choice made inside Permanent; without this the
        // two-argument form is pinned to whichever orientation the caller
        // happened to build.
        Matrix expanded = null;
        for (int j = 0; j < m.length; j++) {
            if (m[j] <= 0) {
                continue;
            }
            Matrix col = Matrix.extractColumn(A, j, null);
            Matrix rep = col.repmat(1, m[j]);
            expanded = (expanded == null) ? rep : Matrix.concatColumns(expanded, rep, null);
        }
        return pfqn_perm(expanded);
    }
}
