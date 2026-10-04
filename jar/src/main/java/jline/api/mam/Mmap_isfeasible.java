/**
 * @file Marked Markovian Arrival Process feasibility validation
 *
 * Validates mathematical feasibility of MMAP representations including stochastic
 * properties and marking consistency. Essential for ensuring valid multiclass models.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_isfeasible {
    private Mmap_isfeasible() {}

    /**
     * Checks the feasibility of a Markovian Arrival Process with marked arrivals (MMAP).
     *
     * This method verifies whether the given MMAP is feasible, which includes checking:
     * 1. The underlying MAP (formed by matrices D0 and D1) is feasible.
     * 2. The elements of each additional event matrix (D1c) are non-negative.
     * 3. The sum of the matrices representing different event types (D1c) does not exceed the visible transition matrix (D1) in absolute value.
     *
     * @param MMAP the MatrixCell containing the transition matrices of the MMAP, with D0, D1, ..., Dc representing different types of events
     * @return true if the MMAP is feasible, false otherwise
     */
    public static boolean mmap_isfeasible(MatrixCell MMAP) {
        MatrixCell MAP = new MatrixCell();
        MAP.set(0, MMAP.get(0));
        MAP.set(1, MMAP.get(1));
        if (!Map_isfeasible.map_isfeasible(MAP)) {
            return false;
        }
        int C = MMAP.size() - 2;
        // elements of D1c are >= 0
        double smallest;
        for (int c = 0; c < C; c++) {
            smallest = MMAP.get(2 + c).elementMin();
            if (smallest < -GlobalConstants.FineTol) return false;
        }
        // D1 = D11 + D12 + ... + D1C. Work on a copy: subEq would otherwise
        // overwrite the caller's D1 with the residual.
        Matrix S = MMAP.get(1).copy();
        for (int c = 0; c < C; c++) {
            S.subEq(MMAP.get(2 + c));
        }
        return !(S.elementMaxAbs() > GlobalConstants.FineTol);
    }
}
