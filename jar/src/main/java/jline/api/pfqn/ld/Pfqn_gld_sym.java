/**
 * @file Symbolic load-dependent normalizing constant of a closed network.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import jline.util.matrix.SymMatrix;
import jline.util.symbolic.SymContext;
import jline.util.symbolic.SymExpr;

/**
 * The gld recursion over EXACT SYMBOLIC demands and rates.
 *
 * <p>Port of the {@code isSym} arms of {@code matlab/src/api/pfqn/pfqn_gld.m},
 * and the twin of the native python object-array arm.
 *
 * <p>Three arms, as the reference has:
 * <ul>
 *   <li>M == 1, the closed form. The multinomial is formed EXACTLY as a ratio of
 *       factorials rather than through {@code exp(factln(...))}: factln returns
 *       a double, and carrying that into exact arithmetic would leave a rational
 *       approximation of a float where a small integer belongs. A class with no
 *       jobs is dropped by the test {@code N > 0}, which is concrete; a class
 *       WITH jobs is assumed to have a non-zero demand, there being no way to
 *       prove otherwise on a symbol.</li>
 *   <li>R == 1, delegated to {@link Pfqn_gldsingle_sym}. NOT to the lld kernel,
 *       whose threshold is found by comparing rates.</li>
 *   <li>otherwise the station/class recursion, which compares nothing.</li>
 * </ul>
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Pfqn_gld_sym {
    private Pfqn_gld_sym() {}

    /**
     * Normalizing constant of a load-dependent closed network, symbolically.
     *
     * @param L  demands (M x R)
     * @param N  populations (1 x R)
     * @param mu load-dependent scaling factors (M x sum(N))
     * @return the normalizing constant
     */
    public static Ret.pfqnNcSym pfqn_gld_sym(SymMatrix L, Matrix N, SymMatrix mu) {
        if (mu.context() != L.context()) {
            throw new RuntimeException("pfqn_gld_sym: L and mu come from different symbolic contexts.");
        }
        return new Ret.pfqnNcSym(gld(L, N, mu));
    }

    private static SymExpr gld(SymMatrix L, Matrix N, SymMatrix mu) {
        SymContext ctx = L.context();
        int M = L.getNumRows();
        int R = L.getNumCols();

        if (M == 0) {
            return ctx.zero();
        }
        int Nt = 0;
        for (int r = 0; r < R; r++) {
            Nt += (int) N.get(r);
        }
        if (Nt == 0) {
            return ctx.one();
        }

        if (M == 1) {
            // multinomial(N) * prod(L(1,r)^N_r over active r) / prod(mu(1,1..Nt))
            SymExpr acc = factorial(ctx, Nt);
            for (int r = 0; r < R; r++) {
                int nr = (int) N.get(r);
                if (nr > 0) {
                    acc = acc.divide(factorial(ctx, nr));
                    acc = acc.multiply(L.get(0, r).pow(nr));
                }
            }
            for (int k = 0; k < Nt; k++) {
                acc = acc.divide(mu.get(0, k));
            }
            return acc;
        }

        if (R == 1) {
            return Pfqn_gldsingle_sym.pfqn_gldsingle_sym(L, N, mu).G;
        }

        // G = gld(L without the last station, N, mu without it)
        //   + sum_r (L(M,r)/mu(M,0)) * gld(L, N - e_r, mushift(mu, M))
        SymExpr G = gld(L.extractRows(0, M - 1), N, mu.extractRows(0, M - 1));
        for (int r = 0; r < R; r++) {
            if (N.get(r) > 0) {
                Matrix N1 = N.copy();
                N1.set(0, r, N.get(r) - 1);
                SymExpr coeff = L.get(M - 1, r).divide(mu.get(M - 1, 0));
                G = G.add(coeff.multiply(gld(L, N1, mushift(mu, M))));
            }
        }
        return G;
    }

    /**
     * {@code pfqn_mushift}: drop the first rate of row {@code m}, shifting the
     * rest of that row left by one and leaving every other row alone.
     */
    private static SymMatrix mushift(SymMatrix mu, int m) {
        SymContext ctx = mu.context();
        int rows = mu.getNumRows();
        int cols = mu.getNumCols();
        SymMatrix out = new SymMatrix(ctx, rows, Math.max(cols - 1, 0));
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j + 1 <= cols - 1; j++) {
                out.set(i, j, i == m - 1 ? mu.get(i, j + 1) : mu.get(i, j));
            }
        }
        return out;
    }

    /**
     * Exact k! as a symbolic constant.
     *
     * <p>Built as a PRODUCT OF SMALL FACTORS rather than from a precomputed
     * long, so that a population past 20! -- where a long silently wraps --
     * still lands exactly; the ring folds the product into one coefficient.
     */
    private static SymExpr factorial(SymContext ctx, int k) {
        SymExpr acc = ctx.one();
        for (int i = 2; i <= k; i++) {
            acc = acc.multiply(ctx.constant((long) i));
        }
        return acc;
    }
}
