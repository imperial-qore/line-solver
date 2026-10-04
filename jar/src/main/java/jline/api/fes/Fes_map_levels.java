/**
 * @file Per-level processes of a load-dependent flow-equivalent server
 *
 * @since LINE 3.0
 */
package jline.api.fes;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Expands a MAP into the per-level processes of a load-dependent server.
 *
 * A flow-equivalent server is described by one MAP (F0^k,F1^k) per population level
 * k=1..n. A single MAP is replicated over the levels and scaled by min(k,mi), which
 * reproduces a queue with mi servers and, for mi infinite, a delay station serving at
 * rate k*mu. The scaling is exact for exponential service and is the load-dependent
 * rate approximation otherwise.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Fes_map_levels {
    private Fes_map_levels() {}

    /**
     * Replicates a load independent MAP over n levels.
     *
     * @param MAP process in the form (D0,D1)
     * @param n   number of levels required
     * @return list of per-level processes
     */
    public static List<MatrixCell> fes_map_levels(MatrixCell MAP, int n) {
        return fes_map_levels(MAP, n, 1);
    }

    /**
     * Replicates a load independent MAP over n levels, scaling level k by min(k,mi).
     *
     * @param MAP process in the form (D0,D1)
     * @param n   number of levels required
     * @param mi  number of servers, Double.POSITIVE_INFINITY for a delay
     * @return list of per-level processes
     */
    public static List<MatrixCell> fes_map_levels(MatrixCell MAP, int n, double mi) {
        List<MatrixCell> levels = new ArrayList<MatrixCell>(n);
        for (int k = 1; k <= n; k++) {
            double s = Math.min((double) k, mi);
            levels.add(new MatrixCell(MAP.get(0).scale(s), MAP.get(1).scale(s)));
        }
        return levels;
    }

    /**
     * Validates a load dependent descriptor and returns its first n levels.
     *
     * @param levels per-level processes, at least n of them
     * @param n      number of levels required
     * @return list of per-level processes
     */
    public static List<MatrixCell> fes_map_levels(List<MatrixCell> levels, int n) {
        if (levels.size() < n) {
            throw new IllegalArgumentException("The flow-equivalent server is defined for "
                    + levels.size() + " levels but " + n + " are required.");
        }
        int mf = levels.get(0).get(0).getNumRows();
        List<MatrixCell> out = new ArrayList<MatrixCell>(n);
        for (int k = 0; k < n; k++) {
            if (levels.get(k).get(0).getNumRows() != mf) {
                throw new IllegalArgumentException(
                        "All levels of a flow-equivalent server must have the same number of phases.");
            }
            out.add(levels.get(k));
        }
        return out;
    }
}
