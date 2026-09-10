package jline.api.mdd;

import java.util.ArrayList;
import java.util.List;

/**
 * Reachability set generation into a decision diagram.
 *
 * <p>After A.S. Miner, G. Ciardo, "Efficient Reachability Set Generation and
 * Storage Using Decision Diagrams", ICATPN 1999, LNCS 1639, pp.6-25.</p>
 */
public class Mdd_reachset {

    private Mdd_reachset() {}

    /**
     * Generate and store the reachability set into a quasi-reduced ordered MDD.
     *
     * <p>This is the basic (explicit-frontier) realisation: a breadth-first
     * search enumerates successors while the MDD provides the O(K) membership
     * test that replaces the usual explicit visited hash. The stored set lives
     * entirely in the MDD (O(#nodes) memory); only the transient BFS frontier is
     * held explicitly. Symbolic image computation / saturation, which removes
     * the explicit frontier too, is the natural next step but is out of scope
     * here.</p>
     *
     * @param domain per-level local-state counts, values 0..domain[k]-1
     * @param init the initial global state, 0-based local values
     * @param nextfun the next-state function
     * @return an MDD holding every state reachable from init
     */
    public static MDD mdd_reachset(int[] domain, int[] init, MddNextState nextfun) {
        MDD mdd = new MDD(domain);
        mdd.insert(init);

        List<int[]> frontier = new ArrayList<int[]>();
        frontier.add(init);
        int head = 0;
        while (head < frontier.size()) {
            int[] s = frontier.get(head);
            head++;
            int[][] successors = nextfun.next(s);
            if (successors != null) {
                for (int r = 0; r < successors.length; r++) {
                    int[] t = successors[r];
                    if (!mdd.member(t)) {
                        mdd.insert(t);
                        frontier.add(t);
                    }
                }
            }
            // drop already-expanded rows periodically to bound frontier memory
            if (head > 1024 && 2 * head > frontier.size()) {
                frontier = new ArrayList<int[]>(frontier.subList(head, frontier.size()));
                head = 0;
            }
        }
        mdd.compact();   // reclaim the dead nodes left by the append-only inserts
        return mdd;
    }
}
