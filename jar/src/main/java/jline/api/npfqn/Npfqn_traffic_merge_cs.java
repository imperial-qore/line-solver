/**
 * @file NPFQN Traffic Merging with Class Switching
 *
 * Implements traffic merging algorithms for non-product-form queueing networks
 * with class switching capabilities. Handles the aggregation of multiple traffic
 * streams while accounting for class transitions in NPFQN analysis.
 *
 * @since LINE 3.0
 */
package jline.api.npfqn;

import jline.api.mam.Mmap_mark;
import jline.api.mam.Mmap_super;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.HashMap;
import java.util.Map;

public final class Npfqn_traffic_merge_cs {
    private Npfqn_traffic_merge_cs() {}

    /**
     * Merges MMAP traffic flows with class switching.
     *
     * @param MMAPs Map of MMAP traffic flows
     * @param prob  Probability matrix for class switching
     * @param config Configuration string ("default" or "super")
     * @return Merged and normalized MMAP traffic flow
     */
    public static MatrixCell npfqn_traffic_merge_cs(Map<Integer, MatrixCell> MMAPs, Matrix prob, String config) {
        int n = MMAPs.size();
        int R = prob.getNumCols();
        Map<Integer, MatrixCell> MMAP_copy = new HashMap<Integer, MatrixCell>();
        for (Map.Entry<Integer, MatrixCell> entry : MMAPs.entrySet()) {
            MMAP_copy.put(entry.getKey(), new MatrixCell(entry.getValue()));
        }
        for (int i = 0; i < n; i++) {
            Matrix P = new Matrix(R, R, R * R);
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    P.set(r, s, prob.get((i - 1) * R + r, s));
                }
            }
            MMAP_copy.put(i, Mmap_mark.mmap_mark(MMAP_copy.get(i), P));
        }
        MatrixCell SMMAP = new MatrixCell();
        if (n == 1) {
            SMMAP = MMAPs.get(0);
        } else {
            if ("default".equals(config) || "super".equals(config)) {
                SMMAP = MMAPs.get(0);
                for (int j = 1; j < n; j++) {
                    SMMAP = Mmap_super.mmap_super(SMMAP, MMAP_copy.get(j), "match");
                }
            }
        }
        return SMMAP;
    }

    public static MatrixCell npfqn_traffic_merge_cs(Map<Integer, MatrixCell> MMAPs, Matrix prob) {
        return npfqn_traffic_merge_cs(MMAPs, prob, "default");
    }
}
