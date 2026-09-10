/**
 * @file Screen for promotion paths blocked by storage cost caps
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

public final class Cache_cost_pathcheck {
    private Cache_cost_pathcheck() {}

    /**
     * One (item, list) pair that is size-feasible for the list but unreachable
     * because an intermediate list on the promotion path rejects the item.
     */
    public static final class BlockedPair {
        public final int item;
        public final int list;
        public final int blockingList;

        public BlockedPair(int item, int list, int blockingList) {
            this.item = item;
            this.list = list;
            this.blockingList = blockingList;
        }
    }

    /**
     * Detects promotion paths blocked by storage cost caps. The constrained
     * normalizing constant E(m,k) of Casale-Gast, IEEE/ACM Trans. Networking
     * 29(2), 2021, Sec. IX, sums over every size-feasible state, while under
     * RR-C(m) an item only reaches list j one list at a time along the path
     * from the miss list. A cap on an intermediate list therefore makes
     * size-feasible states unreachable and E(m,k) normalizes over states the
     * cache never visits. An empty report is a necessary, not sufficient,
     * condition for the two sets to agree.
     *
     * @param gamma Cache access factors (n x h).
     * @param sigma Item storage costs (sizes), one per item.
     * @param k Per-list storage cost caps, one per list.
     * @param parent Parent list of each list, 0-based with -1 for lists rooted in the miss list.
     * @return the blocked pairs, empty when the screen finds none.
     */
    public static List<BlockedPair> cache_cost_pathcheck(Matrix gamma, Matrix sigma, Matrix k, int[] parent) {
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        List<BlockedPair> viol = new ArrayList<BlockedPair>();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < h; j++) {
                if (gamma.get(i, j) == 0.0 || sigma.get(i) > k.get(j)) {
                    continue; // item i never resides in list j anyway
                }
                int l = parent[j];
                while (l >= 0) {
                    if (sigma.get(i) > k.get(l)) {
                        viol.add(new BlockedPair(i, j, l));
                        break;
                    }
                    l = parent[l];
                }
            }
        }
        return viol;
    }
}
