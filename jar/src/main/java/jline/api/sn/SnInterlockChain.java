package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

/**
 * Aggregation of a class-indexed interlock matrix to the chain basis the MVA solvers work in.
 *
 * ILclass(r,s) is the share of the class-s queue that a class-r arrival must not see, the
 * interlocked flow of Franks (1999), Eq. (4.7). Two classes of the same chain belong to the
 * same client, so the diagonal blocks carry no information and the chain diagonal stays zero:
 * an arrival always sees its own chain in full.
 *
 * Reference: G. Franks, "Performance Analysis of Distributed Server Systems", PhD thesis,
 * Carleton University, 1999, Ch. 4.
 */
public final class SnInterlockChain {
    private SnInterlockChain() {}

    public static Matrix snInterlockChain(NetworkStruct sn, Matrix ILclass) {
        if (ILclass == null || ILclass.isEmpty()) {
            return null;
        }
        int R = sn.nclasses;
        if (ILclass.getNumRows() != R || ILclass.getNumCols() != R) {
            throw new RuntimeException(String.format(
                    "snInterlockChain: the interlock matrix is %dx%d but the model has %d classes.",
                    ILclass.getNumRows(), ILclass.getNumCols(), R));
        }
        int K = sn.nchains;
        Matrix ILchain = new Matrix(K, K);
        boolean any = false;
        for (int cr = 0; cr < K; cr++) {
            for (int cs = 0; cs < K; cs++) {
                if (cr == cs) {
                    continue;
                }
                double best = 0.0;
                for (int r = 0; r < R; r++) {
                    if (sn.chains.get(cr, r) == 0) {
                        continue;
                    }
                    for (int s = 0; s < R; s++) {
                        if (sn.chains.get(cs, s) == 0) {
                            continue;
                        }
                        best = Math.max(best, ILclass.get(r, s));
                    }
                }
                if (best > 0) {
                    ILchain.set(cr, cs, best);
                    any = true;
                }
            }
        }
        return any ? ILchain : null;
    }
}
