/**
 * @file Neuse-Chandy SCAT (Self-Correcting Approximation Technique) approximate MVA.
 *
 * SCAT shares the Linearizer fixed point: it carries the mean queue lengths at
 * the target population N and at the R reduced populations N-e_s, and corrects
 * the Bard-Schweitzer proportionality assumption with the fraction difference
 *
 *   Delta(i,r,s) = Q(i,r|N-e_s)/(N-e_s)_r - Q(i,r|N)/N_r,
 *
 * held fixed while an inner MVA fixed point is iterated. It differs from
 * Linearizer in that this correction is refreshed ONCE: SCAT stops after the
 * first pass, where Linearizer performs the fixed three passes of Chandy and
 * Neuse (1982), Sec. 4. Cost is therefore about one third of Linearizer's, and
 * accuracy sits between Bard-Schweitzer (the Delta=0 special case) and
 * Linearizer.
 *
 * SCAT's second departure from Linearizer, fitting a probability mass function
 * centred on the mean queue length at queue-dependent centres instead of
 * propagating the MVA distribution recursion (Krzesinski and Greyling 1984,
 * Sec. 4), does not arise here: this entry point covers single-server and delay
 * stations only, exactly as Pfqn_linearizer does. That mass function is
 * available separately as the "scat" marginal rule of Pfqn_ab_amva.
 *
 * Reference: D. Neuse, K. M. Chandy, "SCAT: A Heuristic Algorithm for Queueing
 * Network Models of Computing Systems", ACM SIGMETRICS Perform. Eval. Rev.
 * 10(3), 1981.
 */
package jline.api.pfqn.mva;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

public final class Pfqn_scat {
    private Pfqn_scat() {}

    public static Ret.pfqnAMVA pfqn_scat(Matrix L,
                                          Matrix N,
                                          Matrix Z,
                                          SchedStrategy[] type,
                                          double tol,
                                          int maxiter) {
        return pfqn_scat(L, N, Z, type, tol, maxiter, null);
    }

    public static Ret.pfqnAMVA pfqn_scat(Matrix L,
                                          Matrix N,
                                          Matrix Z,
                                          SchedStrategy[] type,
                                          double tol,
                                          int maxiter,
                                          Matrix QN0) {
        Matrix alpha = new Matrix(N.getNumRows(), N.getNumCols());
        for (int i = 0; i < alpha.getNumRows(); i++) {
            for (int j = 0; j < alpha.getNumCols(); j++) {
                alpha.set(i, j, 1.0);
            }
        }
        // npasses = 1 is what separates SCAT from Linearizer: one Delta refresh, not three
        return Pfqn_egflinearizer.pfqn_egflinearizer(L, N, Z, type, tol, maxiter, alpha, QN0, 1);
    }

    public static Ret.pfqnAMVA pfqn_scat(Matrix L,
                                          Matrix N,
                                          Matrix Z,
                                          SchedStrategy[] type) {
        return pfqn_scat(L, N, Z, type, 1.0e-8, 1000);
    }
}
