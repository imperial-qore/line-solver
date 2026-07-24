/**
 * @file Marked Markovian Arrival Process class decomposition
 *
 * Extracts individual MAP processes for each marked class from MMAP representations.
 * Used for analyzing per-class behavior and comparing multiclass vs single-class models.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.HashMap;
import java.util.Map;

import jline.util.matrix.MatrixCell;

public final class Mmap_maps {
    private Mmap_maps() {}

    /**
     * Extracts K Markovian Arrival Processes (MAPs) from a given MMAP, one for each class.
     *
     * <p>This method creates a MAP for each class in the MMAP by separating the transitions related to that class.
     * Each resulting MAP has its own set of transition matrices, derived from the original MMAP.
     *
     * @param MMAP the original MMAP
     * @return a map containing K MAPs, each stored in a MatrixCell
     */
    public static Map<Integer, MatrixCell> mmap_maps(MatrixCell MMAP) {
        int K = MMAP.size() - 2;
        Map<Integer, MatrixCell> Maps = new HashMap<Integer, MatrixCell>();
        for (int k = 0; k < K; k++) {
            MatrixCell mc = new MatrixCell();
            mc.set(0, MMAP.get(0).add(1.0, MMAP.get(1)).add(-1.0, MMAP.get(2 + k)));
            mc.set(1, MMAP.get(2 + k));
            Maps.put(k, mc);
        }
        return Maps;
    }
}
