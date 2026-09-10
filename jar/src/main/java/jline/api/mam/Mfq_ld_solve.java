/**
 * @file First/second-order level-dependent (multi-regime) fluid queue solver
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.List;

import jline.lib.butools.mam.LevelDependentFluidSolution;
import jline.lib.butools.mam.SecondOrderLevelDependentFluidSolve;
import jline.util.matrix.Matrix;

public final class Mfq_ld_solve {
    private Mfq_ld_solve() {}

    /** Matrix-exponential solution of a multi-regime first/second-order fluid queue. */
    public static LevelDependentFluidSolution mfq_ld_solve(List<Matrix> Q, List<Matrix> R, List<Matrix> S,
                                                           double[] T, double[] boundaryL, double[] boundaryU,
                                                           List<Matrix> Qt, double prec) {
        return SecondOrderLevelDependentFluidSolve.solve(Q, R, S, T, boundaryL, boundaryU, Qt, prec);
    }
}
