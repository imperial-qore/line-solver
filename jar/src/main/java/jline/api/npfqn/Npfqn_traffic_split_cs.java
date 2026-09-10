/**
 * @file NPFQN Traffic Splitting with Class Switching
 *
 * Implements traffic splitting algorithms for non-product-form queueing networks
 * with class switching capabilities. Handles the decomposition of traffic streams
 * while accounting for class transitions in NPFQN analysis.
 *
 * @since LINE 3.0
 */
package jline.api.npfqn;

import jline.api.mam.Mmap_normalize;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.util.HashMap;
import java.util.Map;

public final class Npfqn_traffic_split_cs {
    private Npfqn_traffic_split_cs() {}

    /**
     * Splits MMAP traffic flows with class switching.
     *
     * @param MMAP MMAP to be split
     * @param P    Class switching matrix
     * @return Map of split MMAPs
     */
    public static Map<Integer, MatrixCell> npfqn_traffic_split_cs(MatrixCell MMAP, Matrix P) {
        MMAP.size();
        int R = P.getNumRows();
        int J = P.getNumCols();
        long M = FastMath.round(J / (float) R);
        Map<Integer, MatrixCell> SMMAP = new HashMap<Integer, MatrixCell>();
        for (int jst = 0; jst < M; jst++) {
            SMMAP.put(jst, new MatrixCell());
            SMMAP.get(jst).set(0, MMAP.get(0).add(1.0, MMAP.get(1)));
            SMMAP.get(jst).set(1, new Matrix(MMAP.get(1).getNumRows(), MMAP.get(1).getNumCols()));
            for (int s = 0; s < R; s++) {
                SMMAP.get(jst).set(2 + s, new Matrix(SMMAP.get(jst).get(0).getNumRows(),
                        SMMAP.get(jst).get(0).getNumCols()));
                for (int r = 0; r < R; r++) {
                    Matrix a = MMAP.get(2 + r).copy();
                    // see _kb/03-api-layer.md for rationale
                    a.scaleEq(P.get(r, jst * R + s));
                    SMMAP.get(jst).set(2 + s, SMMAP.get(jst).get(2 + s).add(1.0, a));
                    SMMAP.get(jst).set(1, SMMAP.get(jst).get(1).add(1.0, a));
                    // the split flow moves OUT of the hidden part D0 (which was
                    // initialized to D0+D1 with all arrivals folded in)
                    SMMAP.get(jst).set(0, SMMAP.get(jst).get(0).add(-1.0, a));
                }
            }
            SMMAP.put(jst, Mmap_normalize.mmap_normalize(SMMAP.get(jst)));
        }
        return SMMAP;
    }
}
