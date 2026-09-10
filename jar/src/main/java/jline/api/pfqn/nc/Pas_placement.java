/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.nc;

/**
 * Placement-order logic of a pass-and-swap (P&S) / order-independent network
 * with swap graph H. Isolates the check that decides which class orderings are
 * feasible (adhere to the placement partial order) for {@link Pfqn_pas_is}.
 *
 * <p>An ordering c = (c_1, ..., c_ell) is feasible iff it is non-decreasing with
 * respect to H, i.e. class a never appears before class b whenever H(b,a) != 0
 * (Comte and Dorsman, 2021, arXiv:2009.12299). Equivalently H(b,a) != 0 means b
 * must be placed before a. The transitive closure of "must precede" yields the
 * precedence matrix P(i,j) = 1 iff class i must be placed before class j, so an
 * ordering is feasible iff every class is placed only after all its
 * P-predecessors. H may be a mere Hasse diagram; the closure makes it explicit.
 *
 * <p>Port of matlab/src/api/pfqn/pas_placement.m.
 */
public final class Pas_placement {
    private Pas_placement() {}

    /**
     * Transitive closure of the "must precede" relation. Input H is the direct
     * relation (H[i][j] != 0 iff i precedes j); the result P[i][j] = 1 iff i
     * precedes j through any chain. Returns null for an empty/null H.
     */
    public static int[][] closure(int[][] H) {
        if (H == null || H.length == 0) {
            return null;
        }
        int R = H.length;
        int[][] P = new int[R][R];
        for (int i = 0; i < R; i++) {
            for (int j = 0; j < R; j++) {
                P[i][j] = (H[i][j] != 0) ? 1 : 0;
            }
        }
        // Iterate R times so paths of any length are captured (early exit on
        // convergence). Pnext = P OR (P * H).
        for (int it = 0; it < R; it++) {
            boolean changed = false;
            int[][] Pnext = new int[R][R];
            for (int i = 0; i < R; i++) {
                for (int j = 0; j < R; j++) {
                    int v = P[i][j];
                    if (v == 0) {
                        for (int k = 0; k < R; k++) {
                            if (P[i][k] != 0 && H[k][j] != 0) {
                                v = 1;
                                break;
                            }
                        }
                    }
                    Pnext[i][j] = v;
                    if (v != P[i][j]) {
                        changed = true;
                    }
                }
            }
            P = Pnext;
            if (!changed) {
                break;
            }
        }
        return P;
    }

    /**
     * Class j is placeable next iff it is present (x[j] &gt; 0) and no still-present
     * class i must precede it: sum_i x[i] * P[i][j] == 0. Returns the placeable
     * class indices. A null P means no constraints (pure OI): every present class.
     */
    public static int[] placeable(int[] x, int[][] P) {
        int R = x.length;
        int count = 0;
        int[] tmp = new int[R];
        for (int j = 0; j < R; j++) {
            if (x[j] <= 0) {
                continue;
            }
            boolean blocked = false;
            if (P != null) {
                for (int i = 0; i < R; i++) {
                    if (x[i] > 0 && P[i][j] != 0) {
                        blocked = true;
                        break;
                    }
                }
            }
            if (!blocked) {
                tmp[count++] = j;
            }
        }
        int[] out = new int[count];
        System.arraycopy(tmp, 0, out, 0, count);
        return out;
    }
}
