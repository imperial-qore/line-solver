package jline.api.mdd;

/**
 * Knobs of the level iteration in {@link Mdd_mcd}.
 *
 * <p>The defaults are deliberately much tighter than a solver-level fixed-point
 * tolerance: the level iteration is an INNER numerical solve and
 * {@link Mdd_mcd} verifies the population invariant at 1e-6, so a loose
 * tolerance converges short of the fixed point and trips that guard. Do not
 * wire SolverOptions.iter_tol (sized for AMVA outer loops) into these.</p>
 */
public class MddMcdOptions {

    /** Convergence tolerance on the level marginals. */
    public double tol = 1e-12;
    /** Maximum coupled sweeps before the iteration is declared non-convergent. */
    public int maxiter = 500;
    /** Print a summary of the level sizes and the iteration count. */
    public boolean verbose = false;
    /** Optional warm-start level vectors, one per paper level; null to start uniform. */
    public double[][] initpik = null;
}
