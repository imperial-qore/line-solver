/**
 * @file Iterative BJB(k) balanced job bounds (Zahorjan et al.) for single-class closed networks
 *
 * BJB(k) at iteration count k; BJB(1) recovers the noniterative Balanced Job Bound. Realized
 * through the validated pfqn_pbh recursion. Ported at parity from MATLAB pfqn_bjbk.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_bjbk {
    private Pfqn_bjbk() {}

    /**
     * Iterative balanced job bounds BJB(k).
     *
     * @param L service demand vector (M x 1)
     * @param N population
     * @param Z think time
     * @param k iteration count (>= 1)
     * @return {Xlo, Xhi} lower/upper throughput bounds
     */
    public static double[] pfqn_bjbk(Matrix L, double N, double Z, int k) {
        return Pfqn_pbh.pfqn_pbh(L, N, Z, k);
    }
}
