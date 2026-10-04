/**
 * @file Chandy-Neuse population-scaled termination cutoff for approximate MVA
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

/**
 * Termination cutoff 1/(4000 + 16*sum(N)) of the Linearizer.
 *
 * Published in K. M. Chandy, D. Neuse, "Linearizer: A Heuristic Algorithm for
 * Queuing Network Models of Computing Systems", Commun. ACM 25(2):126-134,
 * 1982, p.129 and appendix. The iteration continues while
 *
 *   max_{i,r} |Q^I(i,r) - Q^{I-1}(i,r)| / N_r &gt; 1/(4000 + 16*|N|),
 *
 * |N| = sum(N). The paper motivates the scaling with |N|: at large populations
 * removing one job changes the queue lengths very little, so a fixed cutoff
 * would terminate the iteration prematurely. It also notes that the expression
 * stays below 0.00025 even at very small populations.
 *
 * The same expression is what LQNS uses as its termination test, set in the
 * SchweitzerCommon constructor of libmva/src/mva.cc; that code carries no
 * citation, and the paper above is its source.
 *
 * Passing NaN as the tol argument of pfqn_bs / pfqn_egflinearizer selects BOTH
 * this cutoff and the normalized-maximum metric of the paper, which is the
 * published test; passing pfqn_cntol(N) as a plain number selects only the
 * cutoff, with those methods' own convergence metric. NaN is the sentinel
 * because it cannot collide with any legitimate tolerance and it is the one
 * form the MATLAB, Python and C++ twins share (MATLAB and Python additionally
 * accept the string 'cn').
 *
 * @since LINE 3.0
 */
public final class Pfqn_cntol {
    private Pfqn_cntol() {}

    /** Termination cutoff at the given population vector. */
    public static double pfqn_cntol(Matrix N) {
        return 1.0 / (4000.0 + 16.0 * N.elementSum());
    }

    /** Termination cutoff at the given total population. */
    public static double pfqn_cntol(double totalPopulation) {
        return 1.0 / (4000.0 + 16.0 * totalPopulation);
    }

    /** True when tol is the sentinel requesting the Chandy-Neuse test. */
    public static boolean isCntol(double tol) {
        return Double.isNaN(tol);
    }
}
