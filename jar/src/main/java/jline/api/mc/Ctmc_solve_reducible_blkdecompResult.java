/**
 * @file CTMC reducible solve full result type
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.List;

import jline.util.matrix.Matrix;

/**
 * Result of {@link Ctmc_solve_reducible_blkdecomp#ctmc_solve_reducible_blkdecomp_full}.
 */
public final class Ctmc_solve_reducible_blkdecompResult {
    private final Matrix pi;
    private final List<Matrix> pis;
    private final Matrix pi0;
    private final List<List<Integer>> scc;
    private final List<Boolean> isrec;

    public Ctmc_solve_reducible_blkdecompResult(Matrix pi, List<Matrix> pis, Matrix pi0,
                                                List<List<Integer>> scc, List<Boolean> isrec) {
        this.pi = pi;
        this.pis = pis;
        this.pi0 = pi0;
        this.scc = scc;
        this.isrec = isrec;
    }

    public Matrix getPi() { return pi; }
    public List<Matrix> getPis() { return pis; }
    public Matrix getPi0() { return pi0; }
    public List<List<Integer>> getScc() { return scc; }
    public List<Boolean> getIsrec() { return isrec; }
}
