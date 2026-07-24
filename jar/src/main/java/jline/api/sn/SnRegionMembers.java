package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

public final class SnRegionMembers {
    private SnRegionMembers() {}

    /**
     * Station membership mask of finite capacity region f.
     *
     * Membership is read from sn.regionmembers.get(f), which refreshRegions records
     * directly from the region's node list. It cannot be derived from sn.region.get(f):
     * -1 there means "unbounded", which is indistinguishable from "not a member", so a
     * region constrained only by regionlincon (or only by a memory budget) reads as
     * empty and is silently ignored.
     *
     * Rmat and memMat provide the legacy derivation, used only for an sn built before
     * regionmembers existed (for instance one deserialised from an older model file).
     * That derivation carries the ambiguity above and is not equivalent.
     *
     * @param sn     network structure
     * @param f      region index
     * @param Rmat   sn.region.get(f), Matrix(M, K+1), used by the legacy fallback
     * @param memMat sn.regionmaxmem.get(f), Matrix(M, 1) or null, used by the legacy fallback
     * @return boolean[M], true where station i belongs to region f
     */
    public static boolean[] snRegionMembers(NetworkStruct sn, int f, Matrix Rmat, Matrix memMat) {
        int M = Rmat.getNumRows();
        boolean[] mask = new boolean[M];
        if (sn.regionmembers != null && sn.regionmembers.size() > f && sn.regionmembers.get(f) != null
                && sn.regionmembers.get(f).getNumRows() >= M) {
            Matrix mem = sn.regionmembers.get(f);
            for (int i = 0; i < M; i++) {
                mask[i] = mem.get(i, 0) != 0;
            }
            return mask;
        }
        int K = Rmat.getNumCols() - 1;
        for (int i = 0; i < M; i++) {
            boolean m = false;
            for (int c = 0; c <= K; c++) {
                if (Rmat.get(i, c) != -1) { m = true; break; }
            }
            if (!m && memMat != null && memMat.getNumRows() > i && memMat.get(i, 0) != -1) { m = true; }
            mask[i] = m;
        }
        return mask;
    }
}
