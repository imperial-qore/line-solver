/**
 * @file Synchronous call DAG carrying reference-task customers into a layer.
 *
 * Port of {@code matlab/src/api/lqn/lqn_ref_routes.m}.
 *
 * @since LINE 3.0
 */
package jline.api.lqn;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import jline.lang.constant.CallType;
import jline.lang.layered.LayeredNetworkStruct;

/**
 * Resolves, for the caller set of one layer, the synchronous call graph along which reference (REF)
 * task customers descend to those callers, and the mean number of times each entry and each call is
 * invoked per REF cycle.
 *
 * <p>SolverLN gives every caller of a layer its own client chain, each carrying its own thread pool
 * as a population, so a layer holds one customer per caller even when those callers are the same REF
 * customers arriving by different routes. The caller set of a group, {@link Group#members}, is what
 * the refpath interlock method turns into ONE chain.
 *
 * <p>Nodes are ENTRIES, because the work a hop charges depends on which entry was called. Visits are
 * computed topologically, {@code v(u) = sum over parents of v(p)*w(a)*callmean}, and routes are
 * COUNTED in the same pass rather than enumerated. Only SYNC calls are followed: an ASYNC call is
 * send-no-reply, so the customer below it is not the REF customer.
 *
 * <p>INDEX BASE. {@link LayeredNetworkStruct} is 0-based, so every index here (tasks, entries,
 * activities, calls, and the positions into {@link Group#entries}) is 0-based as well.
 */
public final class LqnRefRoutes {
    private LqnRefRoutes() {}

    /** Default route cap, config.interlock_maxpaths. */
    public static final int DEFAULT_MAXPATHS = 32;

    /** One reference task and the synchronous DAG below it. */
    public static final class Group {
        /** Task index of the REF task at the root. */
        public int reftask;
        /** True when the REF task is itself a caller of this layer. */
        public boolean headIsCaller;
        /** Callers of this layer on the DAG, in descent order. */
        public int[] members;
        /** Every DAG entry, topologically ordered, root entries first. */
        public int[] entries;
        /** lqn.parent of each entry. */
        public int[] etask;
        /** True where that entry's task is a caller of this layer. */
        public boolean[] ismember;
        /** Mean invocations of each entry per REF cycle. */
        public double[] vEntry;
        /** Per entry, the activity indices reachable from it. */
        public List<int[]> actweightActs;
        /** Per entry, the mean executions of each of those activities per invocation. */
        public List<double[]> actweightW;
        /** Rows {cidx, fromPos, toPos, aidx, vCall}, positions into entries, vCall per REF cycle. */
        public double[][] calls;
        /** Positions in entries forming the prefix, topological order. */
        public int[] prefixPos;
        /** True where that prefix position is a first caller. */
        public boolean[] prefixTerm;
        /** Distinct REF-to-caller routes, counted. */
        public double npaths;
        /** Min of lqn.maxmult over the DAG tasks. A diagnostic: the chain population is never capped by it. */
        public double poolmin;
    }

    /** The groups found, plus the reason for a fallback. */
    public static final class Result {
        /** One element per reference task reaching the callers; empty when why is non-empty. */
        public final List<Group> groups = new ArrayList<Group>();
        /** Non-empty when the layer must fall back to another interlock method. */
        public String why = "";
    }

    public static Result lqn_ref_routes(LayeredNetworkStruct lqn, Collection<Integer> callers) {
        return lqn_ref_routes(lqn, callers, DEFAULT_MAXPATHS, null);
    }

    /**
     * @param lqn       the layered network structure
     * @param callers   task indices that call the layer's server
     * @param maxpaths  refuse the layer above this many reference-task routes into it
     * @param serverSet server elements of the layer, or null; a prefix node whose task is one of them is
     *                  recursion and refuses the layer
     */
    public static Result lqn_ref_routes(LayeredNetworkStruct lqn, Collection<Integer> callers,
                                        double maxpaths, Collection<Integer> serverSet) {
        Result R = new Result();
        if (lqn.ncalls == 0 || callers == null || callers.isEmpty()) {
            return R;
        }
        boolean[] isCaller = new boolean[lqn.nidx];
        for (Integer c : new TreeSet<Integer>(callers)) {
            isCaller[c] = true;
        }

        List<List<double[]>> succ = syncSuccessors(lqn);

        List<Integer> reftasks = new ArrayList<Integer>();
        for (int t = 0; t < lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if (lqn.isref.get(0, tidx) != 0) {
                reftasks.add(tidx);
            }
        }
        if (reftasks.isEmpty()) {
            return R;
        }

        // A caller reachable from two REF tasks is two INDEPENDENT customer pools, and merging them into one
        // chain would invent a correlation the model does not contain, so the whole layer falls back.
        int[] nrefOf = new int[lqn.nidx];
        for (int r : reftasks) {
            boolean[] seen = reachableEntries(lqn, succ, r);
            TreeSet<Integer> mem = new TreeSet<Integer>();
            for (int i = 0; i < lqn.nidx; i++) {
                if (seen[i]) {
                    int p = (int) lqn.parent.get(0, i);
                    if (p >= 0 && isCaller[p]) {
                        mem.add(p);
                    }
                }
            }
            if (isCaller[r]) {
                mem.add(r);
            }
            for (int m : mem) {
                nrefOf[m]++;
            }
        }
        for (int i = 0; i < lqn.nidx; i++) {
            if (nrefOf[i] > 1) {
                R.why = String.format("task '%s' is reachable from %d reference tasks, whose customer pools are independent",
                        nameOf(lqn, i), nrefOf[i]);
                return R;
            }
        }

        for (int r : reftasks) {
            String[] gwhy = new String[]{""};
            Group grp = buildGroup(lqn, succ, r, isCaller, maxpaths, serverSet, gwhy);
            if (!gwhy[0].isEmpty()) {
                R.why = gwhy[0];
                R.groups.clear();
                return R;
            }
            if (grp != null) {
                R.groups.add(grp);
            }
        }
        return R;
    }

    /**
     * Per calling ENTRY, the rows {cidx, called entry, calling activity, callmean} of its synchronous calls. The
     * calling entry is found by inverting actsof over the entry range, since the parent of an activity is its task.
     */
    private static List<List<double[]>> syncSuccessors(LayeredNetworkStruct lqn) {
        List<List<double[]>> succ = new ArrayList<List<double[]>>(lqn.nidx);
        for (int i = 0; i < lqn.nidx; i++) {
            succ.add(new ArrayList<double[]>());
        }
        int[] entryOfAct = new int[lqn.nidx];
        Arrays.fill(entryOfAct, -1);
        for (int e = 0; e < lqn.nentries; e++) {
            int eidx = lqn.eshift + e;
            List<Integer> acts = lqn.actsof == null ? null : lqn.actsof.get(eidx);
            if (acts != null) {
                for (int a : acts) {
                    entryOfAct[a] = eidx;
                }
            }
        }
        for (int cidx = 0; cidx < lqn.ncalls; cidx++) {
            if (lqn.calltype.get(cidx) != CallType.SYNC) {
                continue;
            }
            int aidx = (int) lqn.callpair.get(cidx, 0);
            int eidxTo = (int) lqn.callpair.get(cidx, 1);
            if (aidx < 0 || eidxTo < 0) {
                continue;
            }
            int eidxFrom = entryOfAct[aidx];
            if (eidxFrom < 0) {
                continue;
            }
            succ.get(eidxFrom).add(new double[]{cidx, eidxTo, aidx, callMeanOf(lqn, cidx)});
        }
        return succ;
    }

    /** Mean number of invocations carried by call cidx. */
    private static double callMeanOf(LayeredNetworkStruct lqn, int cidx) {
        double m = Double.NaN;
        if (lqn.callproc_mean != null && lqn.callproc_mean.get(cidx) != null) {
            m = lqn.callproc_mean.get(cidx);
        }
        if (!isFinite(m) && lqn.callproc != null && lqn.callproc.get(cidx) != null) {
            m = lqn.callproc.get(cidx).getMean();
        }
        if (!isFinite(m)) {
            m = 1;
        }
        return m;
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }

    /** Printable name of an LQN element. */
    private static String nameOf(LayeredNetworkStruct lqn, int idx) {
        if (lqn.hashnames != null && lqn.hashnames.get(idx) != null && !lqn.hashnames.get(idx).isEmpty()) {
            return lqn.hashnames.get(idx);
        }
        return "#" + idx;
    }

    private static List<Integer> entriesOf(LayeredNetworkStruct lqn, int tidx) {
        List<Integer> l = lqn.entriesof == null ? null : lqn.entriesof.get(tidx);
        return l == null ? Collections.<Integer>emptyList() : l;
    }

    /** Every entry reachable from task tidx over SYNC calls, without pruning. */
    private static boolean[] reachableEntries(LayeredNetworkStruct lqn, List<List<double[]>> succ, int tidx) {
        boolean[] seen = new boolean[lqn.nidx];
        Deque<Integer> stack = new ArrayDeque<Integer>();
        for (int e : entriesOf(lqn, tidx)) {
            stack.push(e);
        }
        while (!stack.isEmpty()) {
            int eidx = stack.pop();
            if (seen[eidx]) {
                continue;
            }
            seen[eidx] = true;
            for (double[] row : succ.get(eidx)) {
                stack.push((int) row[1]);
            }
        }
        return seen;
    }

    /** One group: the synchronous DAG below REF task r, its per-cycle visit counts, and the prefix above the first callers. */
    private static Group buildGroup(LayeredNetworkStruct lqn, List<List<double[]>> succ, int r, boolean[] isCaller,
                                    double maxpaths, Collection<Integer> serverSet, String[] why) {
        List<Integer> roots = entriesOf(lqn, r);
        if (roots.isEmpty()) {
            return null;
        }

        // Depth-first sweep with an on-stack marker. A back edge is REFUSED: a cycle makes v(u) a geometric series
        // in a quantity the chain would have to express as a self-loop through the layer's own server.
        final int WHITE = 0, GREY = 1, BLACK = 2;
        int[] color = new int[lqn.nidx];
        List<Integer> post = new ArrayList<Integer>();
        for (int e0 : roots) {
            if (color[e0] != WHITE) {
                continue;
            }
            List<int[]> stack = new ArrayList<int[]>();
            stack.add(new int[]{e0, 0});
            while (!stack.isEmpty()) {
                int[] top = stack.get(stack.size() - 1);
                int u = top[0];
                int ci = top[1];
                if (ci == 0) {
                    color[u] = GREY;
                }
                List<double[]> kids = succ.get(u);
                if (ci < kids.size()) {
                    top[1] = ci + 1;
                    int v = (int) kids.get(ci)[1];
                    if (color[v] == GREY) {
                        why[0] = String.format("the synchronous call graph below '%s' is cyclic", nameOf(lqn, v));
                        return null;
                    } else if (color[v] == WHITE) {
                        stack.add(new int[]{v, 0});
                    }
                } else {
                    color[u] = BLACK;
                    post.add(u);
                    stack.remove(stack.size() - 1);
                }
            }
        }

        int n = post.size();
        if (n == 0) {
            return null;
        }
        int[] entries = new int[n];
        for (int i = 0; i < n; i++) {
            entries[i] = post.get(n - 1 - i); // topological order
        }
        int[] pos = new int[lqn.nidx];
        Arrays.fill(pos, -1);
        for (int i = 0; i < n; i++) {
            pos[entries[i]] = i;
        }
        int[] etask = new int[n];
        boolean[] ismem = new boolean[n];
        boolean anyMem = false;
        for (int i = 0; i < n; i++) {
            etask[i] = (int) lqn.parent.get(0, entries[i]);
            ismem[i] = etask[i] >= 0 && isCaller[etask[i]];
            anyMem |= ismem[i];
        }
        if (!anyMem) {
            return null; // this REF reaches none of the callers
        }

        // Per-entry activity weights, then the call list
        List<int[]> awActs = new ArrayList<int[]>();
        List<double[]> awW = new ArrayList<double[]>();
        List<double[]> calls = new ArrayList<double[]>();
        for (int i = 0; i < n; i++) {
            int u = entries[i];
            int[][] actsOut = new int[1][];
            double[] w = entryActWeights(lqn, u, actsOut, why);
            if (!why[0].isEmpty()) {
                return null;
            }
            awActs.add(actsOut[0]);
            awW.add(w);
            for (double[] row : succ.get(u)) {
                int v = (int) row[1];
                if (pos[v] < 0) {
                    continue;
                }
                int aidx = (int) row[2];
                double wa = 0;
                for (int k = 0; k < actsOut[0].length; k++) {
                    if (actsOut[0][k] == aidx) {
                        wa = w[k];
                        break;
                    }
                }
                calls.add(new double[]{row[0], i, pos[v], aidx, wa * row[3]});
            }
        }

        // Topological visits. The REF task selects among its own entries with equal probability, so a root
        // entry is visited 1/nentries times per REF cycle.
        double[] vEntry = new double[n];
        List<Integer> rootpos = new ArrayList<Integer>();
        for (int e : roots) {
            if (pos[e] >= 0) {
                rootpos.add(pos[e]);
            }
        }
        for (int p : rootpos) {
            vEntry[p] = 1.0 / roots.size();
        }
        for (int i = 0; i < n; i++) {
            for (double[] c : calls) {
                if ((int) c[1] == i) {
                    vEntry[(int) c[2]] += vEntry[i] * c[4];
                }
            }
        }
        for (double[] c : calls) {
            c[4] = vEntry[(int) c[1]] * c[4];
        }

        // The prefix: entries reachable from a root without passing THROUGH a caller. A caller entry reached this
        // way is the prefix's terminal, where the layer's own activity subgraph takes over.
        boolean[] inPrefix = new boolean[n];
        boolean[] prefTerm = new boolean[n];
        double[] npath = new double[n];
        for (int p : rootpos) {
            inPrefix[p] = true;
            npath[p] = 1;
        }
        for (int i = 0; i < n; i++) {
            if (!inPrefix[i]) {
                continue;
            }
            if (ismem[i]) {
                prefTerm[i] = true;
                continue; // do not descend past a caller
            }
            // A HOP whose task is a server of this layer would charge at the client Delay work that the layer's
            // own station exists to serve, and the customer would be in two places at once.
            if (serverSet != null && serverSet.contains(etask[i])) {
                why[0] = String.format("task '%s' is both an intermediate on the reference path and a server of this layer",
                        nameOf(lqn, etask[i]));
                return null;
            }
            for (double[] c : calls) {
                if ((int) c[1] == i) {
                    int j = (int) c[2];
                    inPrefix[j] = true;
                    npath[j] += npath[i];
                }
            }
        }

        double np = 0;
        for (int i = 0; i < n; i++) {
            if (prefTerm[i]) {
                np += npath[i];
            }
        }
        if (np > maxpaths) {
            why[0] = String.format("the reference path into this layer carries %s distinct routes, above config.interlock_maxpaths = %s",
                    fmtG(np), fmtG(maxpaths));
            return null;
        }

        Group g = new Group();
        g.reftask = r;
        g.headIsCaller = isCaller[r];
        Set<Integer> members = new LinkedHashSet<Integer>();
        for (int i = 0; i < n; i++) {
            if (ismem[i]) {
                members.add(etask[i]);
            }
        }
        g.members = toArray(members);
        g.entries = entries;
        g.etask = etask;
        g.ismember = ismem;
        g.vEntry = vEntry;
        g.actweightActs = awActs;
        g.actweightW = awW;
        g.calls = calls.toArray(new double[calls.size()][]);
        List<Integer> pp = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            if (inPrefix[i]) {
                pp.add(i);
            }
        }
        g.prefixPos = toArray(pp);
        g.prefixTerm = new boolean[pp.size()];
        for (int k = 0; k < pp.size(); k++) {
            g.prefixTerm[k] = prefTerm[pp.get(k)];
        }
        g.npaths = np;
        g.poolmin = poolMin(lqn, etask);
        return g;
    }

    /** MATLAB %g of an integral or finite count. */
    private static String fmtG(double x) {
        if (x == Math.rint(x) && Math.abs(x) < 1e15) {
            return String.valueOf((long) x);
        }
        return String.valueOf(x);
    }

    private static int[] toArray(Collection<Integer> c) {
        int[] a = new int[c.size()];
        int k = 0;
        for (int x : c) {
            a[k++] = x;
        }
        return a;
    }

    /** Smallest thread pool on the path. DIAGNOSTIC ONLY: capping the chain at it would understate throughput. */
    private static double poolMin(LayeredNetworkStruct lqn, int[] tasks) {
        double p = Double.POSITIVE_INFINITY;
        if (lqn.maxmult == null) {
            return p;
        }
        for (int t : new TreeSet<Integer>(boxed(tasks))) {
            if (t >= 0 && lqn.maxmult.getNumCols() > t) {
                p = Math.min(p, lqn.maxmult.get(0, t));
            }
        }
        return p;
    }

    private static List<Integer> boxed(int[] a) {
        List<Integer> l = new ArrayList<Integer>(a.length);
        for (int x : a) {
            l.add(x);
        }
        return l;
    }

    /**
     * The mean executions of each activity of entry eidx per invocation of that entry, propagated over lqn.graph so
     * an OR-branch splits its successors by the declared probabilities. Returns the weights and sets actsOut[0] to
     * the activities they belong to. An activity-graph loop sets why and returns null.
     */
    private static double[] entryActWeights(LayeredNetworkStruct lqn, int eidx, int[][] actsOut, String[] why) {
        List<Integer> actsL = lqn.actsof == null ? null : lqn.actsof.get(eidx);
        if (actsL == null || actsL.isEmpty()) {
            actsOut[0] = new int[0];
            return new double[0];
        }
        int na = actsL.size();
        int[] nodeset = new int[na + 1];
        nodeset[0] = eidx;
        for (int k = 0; k < na; k++) {
            nodeset[k + 1] = actsL.get(k);
        }
        int m = nodeset.length;
        int[] pos = new int[lqn.nidx];
        Arrays.fill(pos, -1);
        for (int i = 0; i < m; i++) {
            pos[nodeset[i]] = i;
        }
        double[][] A = new double[m][m];
        for (int i = 0; i < m; i++) {
            int u = nodeset[i];
            for (int v = 0; v < lqn.nidx; v++) {
                double g = lqn.graph.get(u, v);
                if (g != 0 && pos[v] >= 0) {
                    A[i][pos[v]] = g;
                }
            }
        }

        // Topological propagation with an explicit cycle test: a loop makes the executions a geometric series the
        // chain cannot express, so the layer falls back rather than silently truncating the count.
        int[] remaining = new int[m];
        for (int j = 0; j < m; j++) {
            for (int i = 0; i < m; i++) {
                if (A[i][j] > 0) {
                    remaining[j]++;
                }
            }
        }
        remaining[0] = 0; // the entry is the source
        double[] w = new double[m];
        w[0] = 1;
        Deque<Integer> queue = new ArrayDeque<Integer>();
        for (int j = 0; j < m; j++) {
            if (remaining[j] == 0) {
                queue.addLast(j);
            }
        }
        boolean[] done = new boolean[m];
        int ndone = 0;
        while (!queue.isEmpty()) {
            int i = queue.pollFirst();
            if (done[i]) {
                continue;
            }
            done[i] = true;
            ndone++;
            for (int j = 0; j < m; j++) {
                if (A[i][j] > 0) {
                    w[j] += w[i] * A[i][j];
                    remaining[j]--;
                    if (remaining[j] <= 0 && !done[j]) {
                        queue.addLast(j);
                    }
                }
            }
        }
        if (ndone < m) {
            why[0] = String.format("the activity graph of entry '%s' contains a loop", nameOf(lqn, eidx));
            return null;
        }
        int[] acts = new int[na];
        double[] wa = new double[na];
        for (int k = 0; k < na; k++) {
            acts[k] = nodeset[k + 1];
            wa[k] = w[k + 1];
        }
        actsOut[0] = acts;
        return wa;
    }
}
