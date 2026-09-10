/**
 * @file Iterative PB(k) proportional bounds (Eager-Sevcik) for single-class closed networks
 *
 * PB(k) is the level-k member of the Performance Bound Hierarchy; realized through the
 * validated pfqn_pbh recursion. Ported at parity from MATLAB pfqn_pbk.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_pbk {
    private Pfqn_pbk() {}

    /**
     * Iterative proportional bounds PB(k).
     *
     * @param L service demand vector (M x 1)
     * @param N population
     * @param Z think time
     * @param k iteration count (>= 0)
     * @return {Xlo, Xhi} lower/upper throughput bounds
     */
    public static double[] pfqn_pbk(Matrix L, double N, double Z, int k) {
        return Pfqn_pbh.pfqn_pbh(L, N, Z, k);
    }
}
