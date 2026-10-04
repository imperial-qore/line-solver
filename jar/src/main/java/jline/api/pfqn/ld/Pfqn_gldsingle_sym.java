/**
 * @file Symbolic single-class load-dependent auxiliary normalizing constant.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.HashMap;
import java.util.Map;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import jline.util.matrix.SymMatrix;
import jline.util.symbolic.SymExpr;

/**
 * The single-class gld recursion over EXACT SYMBOLIC demands and rates.
 *
 * <p>Port of the {@code isSym} arm of
 * {@code matlab/src/api/pfqn/pfqn_gldsingle.m}, and the twin of the native
 * python object-array arm in {@code api/pfqn/ncld.py}.
 *
 * <p>THE RECURSION IS UNCHANGED FROM THE NUMERIC ONE, which is the whole point:
 * {@code g(m,n,tm) = g(m-1,n,1) + L(m)*g(m,n-1,tm+1)/mu(m,tm)} is {@code +},
 * {@code *} and {@code /} throughout, so it is field-agnostic and runs over a
 * rational function exactly as it runs over a double.
 *
 * <p>WHAT IS DROPPED, AND WHY, mirroring the reference: the LOG-DOMAIN path,
 * because there is no log here and none is needed (the log domain exists to stop
 * a double underflowing, and an exact rational cannot); and the ZERO-DEMAND
 * GUARD and the LOAD-INDEPENDENCE SCAN, because {@code L > 0} and
 * {@code min(row) == 1} are comparisons with no truth value on a symbol. The
 * reference resolves the first by selecting on {@code N}, which is concrete, and
 * skips the second in favour of the recursion, which compares nothing. So does
 * this.
 *
 * <p>{@code pfqn_lldsingle} is NOT the kernel here, for the same reason it is not
 * in the reference: its saving rests on locating the index past which a rate row
 * is constant, and that is a comparison.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Pfqn_gldsingle_sym {
    private Pfqn_gldsingle_sym() {}

    /**
     * Normalizing constant of a single-class load-dependent model, symbolically.
     *
     * @param L  demands at all stations (M x 1)
     * @param N  number of jobs, one class
     * @param mu load-dependent scaling factors (M x sum(N))
     * @return the normalizing constant
     */
    public static Ret.pfqnNcSym pfqn_gldsingle_sym(SymMatrix L, Matrix N, SymMatrix mu) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (R > 1) {
            throw new RuntimeException(
                    "pfqn_gldsingle_sym: multiclass model detected. pfqn_gldsingle_sym is for single class models.");
        }
        if (mu.context() != L.context()) {
            throw new RuntimeException(
                    "pfqn_gldsingle_sym: L and mu come from different symbolic contexts.");
        }
        int Nt = (int) N.get(0);
        SymExpr zero = L.context().zero();
        SymExpr one = L.context().one();

        // Sparse, keyed exactly as the complex twin is: the stencil is
        // triangular in (n,tm) and a dense cube would allocate what it never reads.
        Map<Ret.pfqnGldIndex, SymExpr> g = new HashMap<Ret.pfqnGldIndex, SymExpr>();
        g.put(new Ret.pfqnGldIndex(1, 1, 1), zero);
        for (int n = 1; n <= Nt; n++) {
            g.put(new Ret.pfqnGldIndex(1, n + 1, 2), zero);
        }
        for (int m = 1; m <= M; m++) {
            for (int tm = 1; tm <= Nt + 1; tm++) {
                g.put(new Ret.pfqnGldIndex(m + 1, 1, tm + 1), one);
            }
            for (int nn = 1; nn <= Nt; nn++) {
                for (int tmm = 1; tmm <= Nt - nn + 1; tmm++) {
                    SymExpr prev = g.get(new Ret.pfqnGldIndex(m, nn + 1, 2));
                    SymExpr feed = g.get(new Ret.pfqnGldIndex(m + 1, nn, tmm + 2));
                    g.put(new Ret.pfqnGldIndex(m + 1, nn + 1, tmm + 1),
                            prev.add(L.get(m - 1, 0).multiply(feed).divide(mu.get(m - 1, tmm - 1))));
                }
            }
        }
        return new Ret.pfqnNcSym(g.get(new Ret.pfqnGldIndex(M + 1, Nt + 1, 2)));
    }
}
