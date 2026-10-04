package jline.lang.state;

import jline.lang.constant.HeteroSchedPolicy;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Server-pool bookkeeping of a heterogeneous station (Queue.addServerType) as rewritten
 * for the CTMC by {@code jline.solvers.ctmc.handlers.Solver_ctmc_pools}; port of MATLAB
 * sn.nodeparam{ind}.ctmcpool. The class-r server block holds one sub-block of phases per
 * compatible pool, in ascending pool order, so a phase identifies both the pool and the
 * service phase. Indices are 0-based.
 *
 * <p>ALIS and FAIRNESS keep a global pool order. When some class has two or more
 * compatible pools that order is part of the state, as a 1-based index into
 * {@link #perms} held in the local-variable column nvars(ind, 2R).</p>
 */
public class CtmcPool implements Serializable {
    /** Number of pools. */
    public int ntypes;
    /** Servers per pool. */
    public int[] count;
    /** compat[t][r]: pool t serves class r. */
    public boolean[][] compat;
    /** Pool-selection policy of the station. */
    public HeteroSchedPolicy policy;
    /** pools[r]: compatible pools of class r, ascending; empty when r is not served. */
    public int[][] pools;
    /** off[r][k], len[r][k]: offset and length of block k in the class-r server block. */
    public int[][] off;
    public int[][] len;
    /** alpha[r][k], exit[r][k]: entry vector and exit rates of block k. */
    public double[][][] alpha;
    public double[][][] exit;
    /** D0 of block k of class r. */
    public Matrix[][] D0;
    /** fsfrate[t][r]: 1/mean of the law of r at pool t (FSF). */
    public double[][] fsfrate;
    /** Pools by ascending number of compatible classes, ties by index (ALFS). */
    public int[] alfsorder;
    /** True when the ALIS/FAIRNESS pool order is part of the state. */
    public boolean rotate;
    /** Pool orders, lexicographic (identity first). */
    public int[][] perms;
    /** 0-based position of the pool-order variable inside the local-variable block. */
    public int varpos;

    private Map<String, Integer> permIndex;

    /** 1-based index of a pool order in {@link #perms}. */
    public int permIndexOf(int[] order) {
        if (permIndex == null) {
            permIndex = new HashMap<String, Integer>();
            for (int i = 0; i < perms.length; i++) {
                permIndex.put(Arrays.toString(perms[i]), i + 1);
            }
        }
        return permIndex.get(Arrays.toString(order));
    }

    /** Index of pool t among the compatible pools of class r, or -1. */
    public int blockOf(int r, int t) {
        for (int k = 0; k < pools[r].length; k++) {
            if (pools[r][k] == t) {
                return k;
            }
        }
        return -1;
    }
}
