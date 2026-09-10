package jline.api.mdd;

import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Quasi-reduced ordered Multi-valued Decision Diagram.
 *
 * <p>Compact symbolic store for a set of discrete states, after A.S. Miner,
 * G. Ciardo, "Efficient Reachability Set Generation and Storage Using Decision
 * Diagrams", ICATPN 1999, LNCS 1639, pp.6-25.</p>
 *
 * <p>A global state is a K-tuple of <em>local</em> state values (one per
 * level/submodel), state[k] in {0,...,domain[k]-1}. The set is stored as a
 * directed acyclic graph with K variable levels plus a terminal level: level 0
 * is the top (root), a node at level k has domain[k] outgoing arcs to level
 * k+1 nodes, and a state belongs to the set iff its path of arcs reaches the
 * TRUE terminal. Canonicity is enforced by a per-level unique table (no
 * duplicate nodes) and by collapsing the all-FALSE node to the FALSE terminal.
 * Storage is O(#nodes), typically O(K * #local-states), instead of O(|S|) as in
 * an explicit state list.</p>
 *
 * <p>Two constant terminals encode the boolean value of a completed path:
 * TERM_FALSE = 0 (empty subgraph) and TERM_TRUE = -1 (state in set). Node ids
 * are 1-based positive integers so that 0 can serve as TERM_FALSE, matching the
 * MATLAB and python implementations.</p>
 *
 * @see Mdd_reachset
 * @see Mdd_mcd
 */
public class MDD {

    /** Terminal node "1": a completed path is accepted. */
    public static final int TERM_TRUE = -1;
    /** Terminal node "0": empty subgraph. */
    public static final int TERM_FALSE = 0;

    /** Number of variable levels. */
    public final int K;
    /** domain[k] is the number of local states at level k. */
    public final int[] domain;
    /** node.get(k).get(id-1) holds the arc row of level-k node id. */
    private final List<List<int[]>> node;
    /** Per-level unique table, arc row -> node id. */
    private final List<Map<ArcKey, Integer>> uniq;
    /** Id of the top (level-0) node; TERM_FALSE for the empty set. */
    private int root;
    /** Memoized per-node state counts; -1 marks "not yet computed". */
    private long[][] cnt;
    private boolean dirty;

    /** Create an empty set over the given per-level domains. */
    public MDD(int[] domain) {
        this.domain = Arrays.copyOf(domain, domain.length);
        this.K = this.domain.length;
        this.node = new ArrayList<List<int[]>>(K);
        this.uniq = new ArrayList<Map<ArcKey, Integer>>(K);
        for (int k = 0; k < K; k++) {
            this.node.add(new ArrayList<int[]>());
            this.uniq.add(new HashMap<ArcKey, Integer>());
        }
        this.root = TERM_FALSE;
        this.cnt = new long[K][];
        this.dirty = true;
    }

    /** Build an MDD from a set of 0-based state tuples. */
    public static MDD fromStates(int[] domain, int[][] states) {
        MDD obj = new MDD(domain);
        for (int i = 0; i < states.length; i++) {
            obj.insert(states[i]);
        }
        return obj;
    }

    /** Add a K-tuple of 0-based local values to the set. */
    public void insert(int[] state) {
        this.root = addState(0, this.root, state);
        this.dirty = true;
    }

    /** True iff state is in the set; O(K). */
    public boolean member(int[] state) {
        int id = this.root;
        for (int k = 0; k < K; k++) {
            if (id == TERM_FALSE) {
                return false;
            }
            id = node.get(k).get(id - 1)[state[k]];
        }
        return id == TERM_TRUE;
    }

    /** |S|, the number of stored states. */
    public long cardinality() {
        ensureCounts();
        return childCount(0, this.root);
    }

    /** Id of the root node, or TERM_FALSE for the empty set. */
    public int getRoot() {
        return this.root;
    }

    /** Number of live nodes at level k. */
    public int nodeCount(int k) {
        return node.get(k).size();
    }

    /** Arc row of level-k node id (1-based id); the array is the live one, do not mutate. */
    public int[] arcs(int k, int id) {
        return node.get(k).get(id - 1);
    }

    /**
     * 0-based lexicographic rank of state among the stored set, level 0 most
     * significant, or -1 when the state is not stored.
     *
     * <p>This is a bijection S to {0,...,|S|-1}, so a generator matrix can be
     * assembled without an explicit state list.</p>
     */
    public long index(int[] state) {
        ensureCounts();
        long idx = 0;
        int id = this.root;
        for (int k = 0; k < K; k++) {
            if (id == TERM_FALSE) {
                return -1;
            }
            int[] a = node.get(k).get(id - 1);
            int v = state[k];
            for (int vv = 0; vv < v; vv++) {
                idx += childCount(k + 1, a[vv]);
            }
            id = a[v];
        }
        if (id != TERM_TRUE) {
            return -1;
        }
        return idx;
    }

    /** All stored states as rows, in index() order. */
    public int[][] enumerate() {
        if (this.root == TERM_FALSE) {
            return new int[0][K];
        }
        List<int[]> out = new ArrayList<int[]>();
        enumBelow(0, this.root, new int[K], out);
        return out.toArray(new int[out.size()][]);
    }

    /** Export the diagram as plain arrays for downstream algorithms. */
    public MddStruct toStruct() {
        int[] nnodes = new int[K];
        int[][][] tbl = new int[K][][];
        for (int k = 0; k < K; k++) {
            List<int[]> lvl = node.get(k);
            nnodes[k] = lvl.size();
            tbl[k] = new int[lvl.size()][];
            for (int p = 0; p < lvl.size(); p++) {
                tbl[k][p] = Arrays.copyOf(lvl.get(p), domain[k]);
            }
        }
        return new MddStruct(K, Arrays.copyOf(domain, K), this.root, nnodes, tbl);
    }

    /**
     * Reclaim dead nodes left by the append-only build.
     *
     * <p>Membership, index and enumerate are unchanged. A production MDD would
     * reference-count instead and never accumulate dead nodes; this is the
     * basic sweep.</p>
     */
    public void compact() {
        boolean[][] vis = reachableIds();
        List<List<int[]>> newnode = new ArrayList<List<int[]>>(K);
        int[][] remap = new int[K][];
        for (int k = 0; k < K; k++) {
            List<int[]> lvl = node.get(k);
            remap[k] = new int[lvl.size() + 1];
            List<int[]> kept = new ArrayList<int[]>();
            for (int p = 0; p < lvl.size(); p++) {
                if (vis[k][p]) {
                    kept.add(lvl.get(p));
                    remap[k][p + 1] = kept.size();
                }
            }
            newnode.add(kept);
        }
        for (int k = 0; k < K - 1; k++) {
            for (int[] row : newnode.get(k)) {
                for (int v = 0; v < domain[k]; v++) {
                    if (row[v] > 0) {
                        row[v] = remap[k + 1][row[v]];
                    }
                }
            }
        }
        for (int k = 0; k < K; k++) {
            node.set(k, newnode.get(k));
        }
        if (this.root != TERM_FALSE) {
            this.root = remap[0][this.root];
        }
        for (int k = 0; k < K; k++) {
            Map<ArcKey, Integer> table = uniq.get(k);
            table.clear();
            List<int[]> lvl = node.get(k);
            for (int p = 0; p < lvl.size(); p++) {
                table.put(new ArcKey(lvl.get(p)), Integer.valueOf(p + 1));
            }
        }
        this.dirty = true;
    }

    /** Storage description of the current set; only reachable nodes are counted. */
    public MddStats stats() {
        boolean[][] vis = reachableIds();
        int[] perLevel = new int[K];
        int numNodes = 0;
        long mddInts = 0;
        int tableNodes = 0;
        for (int k = 0; k < K; k++) {
            int c = 0;
            for (int p = 0; p < vis[k].length; p++) {
                if (vis[k][p]) {
                    c++;
                }
            }
            perLevel[k] = c;
            numNodes += c;
            mddInts += (long) c * domain[k];
            tableNodes += node.get(k).size();
        }
        long numStates = cardinality();
        return new MddStats(K, perLevel, numNodes, tableNodes, numStates, mddInts,
                numStates * K);
    }

    @Override
    public String toString() {
        MddStats s = stats();
        StringBuilder sb = new StringBuilder();
        sb.append("  MDD  ").append(K).append(" levels, domains [");
        for (int k = 0; k < K; k++) {
            sb.append(domain[k]);
            if (k < K - 1) {
                sb.append(' ');
            }
        }
        sb.append("]\n       ").append(s.numStates).append(" states stored in ")
          .append(s.numNodes).append(" nodes\n       footprint ").append(s.mddInts)
          .append(" ints vs ").append(s.explicitInts).append(" explicit (")
          .append(String.format("%.1f", s.compression())).append("x compression)");
        if (s.tableNodes > s.numNodes) {
            sb.append("\n       (").append(s.tableNodes - s.numNodes)
              .append(" dead nodes in tables; call compact() to reclaim)");
        }
        return sb.toString();
    }

    // -- internals -----------------------------------------------------------

    /** Canonical node creation through the per-level unique table. */
    private int makeNode(int k, int[] arcRow) {
        boolean allFalse = true;
        for (int v = 0; v < arcRow.length; v++) {
            if (arcRow[v] != TERM_FALSE) {
                allFalse = false;
                break;
            }
        }
        if (allFalse) {
            return TERM_FALSE;                       // collapse the empty node
        }
        ArcKey key = new ArcKey(arcRow);
        Integer found = uniq.get(k).get(key);
        if (found != null) {
            return found.intValue();
        }
        node.get(k).add(arcRow);
        int id = node.get(k).size();
        uniq.get(k).put(key, Integer.valueOf(id));
        return id;
    }

    /**
     * Recursively add one state below node id at level k.
     *
     * <p>Nodes are immutable and shared, so this rebuilds the path bottom-up
     * rather than mutating in place.</p>
     */
    private int addState(int k, int id, int[] state) {
        if (k >= K) {
            return TERM_TRUE;
        }
        int[] arcRow;
        if (id == TERM_FALSE) {
            arcRow = new int[domain[k]];             // all TERM_FALSE
        } else {
            arcRow = Arrays.copyOf(node.get(k).get(id - 1), domain[k]);
        }
        int v = state[k];
        arcRow[v] = addState(k + 1, arcRow[v], state);
        return makeNode(k, arcRow);
    }

    private long childCount(int k, int childId) {
        if (k >= K) {
            return childId == TERM_TRUE ? 1L : 0L;
        }
        if (childId == TERM_FALSE) {
            return 0L;
        }
        return countNode(k, childId);
    }

    private long countNode(int k, int id) {
        long c = cnt[k][id - 1];
        if (c >= 0) {
            return c;
        }
        int[] a = node.get(k).get(id - 1);
        c = 0;
        for (int v = 0; v < domain[k]; v++) {
            c += childCount(k + 1, a[v]);
        }
        cnt[k][id - 1] = c;
        return c;
    }

    private void ensureCounts() {
        boolean sized = !dirty;
        if (sized) {
            for (int k = 0; k < K; k++) {
                if (cnt[k] == null || cnt[k].length != node.get(k).size()) {
                    sized = false;
                    break;
                }
            }
        }
        if (sized) {
            return;
        }
        for (int k = 0; k < K; k++) {
            cnt[k] = new long[node.get(k).size()];
            Arrays.fill(cnt[k], -1L);
        }
        dirty = false;
    }

    /** Per-level masks of nodes reachable from the root. */
    private boolean[][] reachableIds() {
        boolean[][] vis = new boolean[K][];
        for (int k = 0; k < K; k++) {
            vis[k] = new boolean[node.get(k).size()];
        }
        if (this.root == TERM_FALSE) {
            return vis;
        }
        vis[0][this.root - 1] = true;
        Deque<int[]> stack = new ArrayDeque<int[]>();
        stack.push(new int[]{0, this.root});
        while (!stack.isEmpty()) {
            int[] top = stack.pop();
            int k = top[0];
            int id = top[1];
            if (k == K - 1) {
                continue;                             // children are terminals
            }
            int[] a = node.get(k).get(id - 1);
            for (int v = 0; v < domain[k]; v++) {
                int ch = a[v];
                if (ch > 0 && !vis[k + 1][ch - 1]) {
                    vis[k + 1][ch - 1] = true;
                    stack.push(new int[]{k + 1, ch});
                }
            }
        }
        return vis;
    }

    private void enumBelow(int k, int id, int[] prefix, List<int[]> out) {
        int[] a = node.get(k).get(id - 1);
        if (k == K - 1) {
            for (int v = 0; v < domain[k]; v++) {
                if (a[v] == TERM_TRUE) {
                    prefix[k] = v;
                    out.add(Arrays.copyOf(prefix, K));
                }
            }
            return;
        }
        for (int v = 0; v < domain[k]; v++) {
            int child = a[v];
            if (child != TERM_FALSE) {
                prefix[k] = v;
                enumBelow(k + 1, child, prefix, out);
            }
        }
    }

    /** Hashable wrapper of an arc row, the key of the per-level unique table. */
    private static final class ArcKey {
        private final int[] arcs;
        private final int hash;

        ArcKey(int[] arcs) {
            this.arcs = arcs;
            this.hash = Arrays.hashCode(arcs);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) {
                return true;
            }
            if (!(o instanceof ArcKey)) {
                return false;
            }
            return Arrays.equals(this.arcs, ((ArcKey) o).arcs);
        }

        @Override
        public int hashCode() {
            return hash;
        }
    }
}
