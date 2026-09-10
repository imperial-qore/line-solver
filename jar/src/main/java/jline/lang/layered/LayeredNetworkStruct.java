/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.Copyable;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.CallType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.DiscreteDistribution;
import jline.lang.processes.Distribution;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Class summarizing the characteristics of a LayeredNetwork
 *
 * <p><b>Index base.</b> Element indices are 0-BASED. A host, task, entry,
 * activity or call has a local index in {@code 0..n-1}, and its absolute index
 * is {@code shift + local}, so absolute indices run {@code 0..nidx-1} and calls
 * run {@code 0..ncalls-1}. The shift values themselves are element COUNTS
 * ({@code hshift=0}, {@code tshift=nhosts}, {@code eshift=nhosts+ntasks}, ...)
 * and are the same numbers they were under the former 1-based numbering; only
 * the local range moved. Matrices carry no padding row or column, and an unset
 * or not-found index is -1, never 0 (0 is the first host). Prefer the
 * {@code hostIdx/taskIdx/entryIdx/actIdx} converters and the {@code isHostIdx}
 * family below over open-coded {@code shift + local} arithmetic.
 *
 * <p>MATLAB keeps the 1-based numbering; {@code JLINE.from_jline_struct_layered}
 * translates. See {@code _kb/04-networkstruct.md}.
 */
public class LayeredNetworkStruct implements Copyable {
    public int nidx;
    public int nhosts;
    public int ntasks;
    public int nacts;
    public int nentries;
    public int ncalls;

    public int tshift;
    public int eshift;
    public int ashift;
    public int hshift;
    public int cshift;

    public Map<Integer, List<Integer>> tasksof;
    public Map<Integer, List<Integer>> entriesof;
    public Map<Integer, List<Integer>> actsof;
    public Map<Integer, List<Integer>> callsof;

    // Host demand distribution and primitive representations
    @Deprecated // Use hostdem_type, hostdem_params, hostdem_mean, hostdem_scv, hostdem_proc instead
    public Map<Integer, Distribution> hostdem;
    public Map<Integer, ProcessType> hostdem_type;
    public Map<Integer, Matrix> hostdem_params;
    public Map<Integer, Double> hostdem_mean;
    public Map<Integer, Double> hostdem_scv;
    public Map<Integer, MatrixCell> hostdem_proc;

    // Think time distribution and primitive representations
    @Deprecated // Use think_type, think_params, think_mean, think_scv, think_proc instead
    public Map<Integer, Distribution> think;
    public Map<Integer, ProcessType> think_type;
    public Map<Integer, Matrix> think_params;
    public Map<Integer, Double> think_mean;
    public Map<Integer, Double> think_scv;
    public Map<Integer, MatrixCell> think_proc;

    // Activity think time distribution and primitive representations
    @Deprecated // Use actthink_type, actthink_params, actthink_mean, actthink_scv, actthink_proc instead
    public Map<Integer, Distribution> actthink;
    public Map<Integer, ProcessType> actthink_type;
    public Map<Integer, Matrix> actthink_params;
    public Map<Integer, Double> actthink_mean;
    public Map<Integer, Double> actthink_scv;
    public Map<Integer, MatrixCell> actthink_proc;

    public Map<Integer, int[]> itemcap;  // Changed to int[] to support multi-level caches

    // Item process distribution and primitive representations
    @Deprecated // Use itemproc_type, itemproc_params, itemproc_mean, itemproc_scv, itemproc_proc instead
    public Map<Integer, DiscreteDistribution> itemproc;
    public Map<Integer, ProcessType> itemproc_type;
    public Map<Integer, Matrix> itemproc_params;
    public Map<Integer, Double> itemproc_mean;
    public Map<Integer, Double> itemproc_scv;
    public Map<Integer, MatrixCell> itemproc_proc;

    public Map<Integer, CallType> calltype;

    // Call process distribution and primitive representations
    @Deprecated // Use callproc_type, callproc_params, callproc_mean, callproc_scv, callproc_proc instead
    public Map<Integer, Distribution> callproc;
    public Map<Integer, ProcessType> callproc_type;
    public Map<Integer, Matrix> callproc_params;
    public Map<Integer, Double> callproc_mean;
    public Map<Integer, Double> callproc_scv;
    public Map<Integer, MatrixCell> callproc_proc;

    public Map<Integer, String> callnames;
    public Map<Integer, String> callhashnames;

    public Map<Integer, SchedStrategy> sched;
    /**
     * Scheduling priority of a task, 0 on every other index. Lower is served
     * first, the convention sn.classprio uses in a Network. Read by the priority
     * disciplines (HOL) of a host or a task.
     */
    public Matrix prio;
    public Map<Integer, String> names;
    public Map<Integer, String> hashnames;
    public Matrix mult;
    public Matrix maxmult;
    public Matrix repl;
    public Matrix type;
    public Matrix graph;
    public Matrix dag;

    public Matrix replygraph;
    public Matrix actphase;  // Phase number (1 or 2) for each activity
    public Matrix nitems;

    public Matrix replacestrat;
    public Matrix replacement; // Alias for replacestrat for Python compatibility
    public Matrix iscache;
    public Matrix hasretrieval; // 1 on a CacheTask index whose miss path is a delayed-hit retrieval
    public Matrix parent;

    // Matrix representation of scheduling strategies (for Python compatibility)
    public Matrix schedid;

    public Matrix iscaller;
    public Matrix issynccaller;
    public Matrix isasynccaller;
    public Matrix callpair;
    /**
     * Calls dispatched as a group by a routing strategy. Empty for every model
     * that does not use Activity.synchCallRoundRobin. Requires the squashed
     * layering; see SolverLN.assertCallGroups.
     */
    public List<CallGroupStruct> callgroups = new ArrayList<>();

    /**
     * A caller activity, the strategy, and the entries it dispatches over.
     * Serializable because LayeredNetworkStruct deep-copies through
     * ObjectOutputStream.
     */
    public static class CallGroupStruct implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
        public final int caller;
        public final RoutingStrategy strategy;
        public final List<Integer> targets;

        public CallGroupStruct(int caller, RoutingStrategy strategy, List<Integer> targets) {
            this.caller = caller;
            this.strategy = strategy;
            this.targets = targets;
        }
    }

    public Matrix taskgraph;

    public Matrix actpretype;
    public Matrix actposttype;

    // Quorum matrix: actquorum(0, join_target_idx) = k, the number of PRE_AND predecessors
    // that must complete before the join fires. Equals the predecessor count when the join
    // has no quorum (wait for all). Indexed by the global index of the join target activity.
    public Matrix actquorum;

    public Matrix isref;

    public Matrix hassetup;

    // Fan-out matrix: fanout(source_task_idx, dest_task_idx) = fan-out value
    // Used for replication load distribution (fork-join semantics)
    public Matrix fanout;

    // Setup time distribution and primitive representations
    @Deprecated // Use setuptime_type, setuptime_params, setuptime_mean, setuptime_scv, setuptime_proc instead
    public Map<Integer, Distribution> setuptime;
    public Map<Integer, ProcessType> setuptime_type;
    public Map<Integer, Matrix> setuptime_params;
    public Map<Integer, Double> setuptime_mean;
    public Map<Integer, Double> setuptime_scv;
    public Map<Integer, MatrixCell> setuptime_proc;

    // Delay-off time distribution and primitive representations
    @Deprecated // Use delayofftime_type, delayofftime_params, delayofftime_mean, delayofftime_scv, delayofftime_proc instead
    public Map<Integer, Distribution> delayofftime;
    public Map<Integer, ProcessType> delayofftime_type;
    public Map<Integer, Matrix> delayofftime_params;
    public Map<Integer, Double> delayofftime_mean;
    public Map<Integer, Double> delayofftime_scv;
    public Map<Integer, MatrixCell> delayofftime_proc;

    // Arrival distribution and primitive representations (for entries with open arrivals)
    @Deprecated // Use arrival_type, arrival_params, arrival_mean, arrival_scv, arrival_proc instead
    public Map<Integer, Distribution> arrival;
    public Map<Integer, ProcessType> arrival_type;
    public Map<Integer, Matrix> arrival_params;
    public Map<Integer, Double> arrival_mean;
    public Map<Integer, Double> arrival_scv;
    public Map<Integer, MatrixCell> arrival_proc;

    // Admission constraint A*n <= b on the layer station of host or task i, keyed by global index i.
    // get(i)[0] is Matrix(C_i,K_i) and get(i)[1] is Matrix(C_i,1). Columns are that host's tasks, or
    // that task's entries, in tasksof/entriesof order. Absent where unconstrained -- see _kb/04-networkstruct.md
    public Map<Integer, Matrix[]> lincon;

    // Service-rate dependences on the layer station of host or task i, keyed by global index i.
    // lldscaling.get(i) is the vector alpha(n) applied at total population n; cdscaling and jdscaling
    // are handles over that server's operands, in the tasksof/entriesof order used by lincon columns,
    // with cdscalingpeak/jdscalingpeak their per-operand peak rate scaling -- see _kb/04-networkstruct.md
    public Map<Integer, Matrix> lldscaling;
    public Map<Integer, SerializableFunction<Matrix, Matrix>> cdscaling;
    public Map<Integer, Matrix> cdscalingpeak;
    public Map<Integer, SerializableFunction<Matrix, Matrix>> jdscaling;
    public Map<Integer, Matrix> jdscalingpeak;

    /**
     * Compatibility pools declared on a layer server, by element index, absent
     * where none. SolverLN lowers each to the activated-server rate of
     * {@link jline.api.sn.SnCompatRate}; see _kb/04-networkstruct.md.
     */
    public Map<Integer, ServerPools> pools;

    /** The resolved pool block of one server: compat(t, j) nonzero = pool t may serve operand j. */
    public static class ServerPools implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
        public List<String> names;
        public Matrix counts;
        public Matrix rates;
        public Matrix compat;

        public ServerPools(List<String> names, Matrix counts, Matrix rates, Matrix compat) {
            this.names = names;
            this.counts = counts;
            this.rates = rates;
            this.compat = compat;
        }
    }

    public Matrix conntasks; // 1 * n matrix
    
    public List<Integer> hitmissaidx;
    public Integer hitaidx;
    public Integer missaidx;

    public LayeredNetworkStruct() {

    }

    /** Absolute index of host h, h in 0..nhosts-1. */
    public int hostIdx(int h) {
        return hshift + h;
    }

    /** Absolute index of task t, t in 0..ntasks-1. */
    public int taskIdx(int t) {
        return tshift + t;
    }

    /** Absolute index of entry e, e in 0..nentries-1. */
    public int entryIdx(int e) {
        return eshift + e;
    }

    /** Absolute index of activity a, a in 0..nacts-1. */
    public int actIdx(int a) {
        return ashift + a;
    }

    /** Host-local index of an absolute index, or -1 if it is not a host. */
    public int hostOf(int idx) {
        return isHostIdx(idx) ? idx - hshift : -1;
    }

    /** Task-local index of an absolute index, or -1 if it is not a task. */
    public int taskOf(int idx) {
        return isTaskIdx(idx) ? idx - tshift : -1;
    }

    /** Entry-local index of an absolute index, or -1 if it is not an entry. */
    public int entryOf(int idx) {
        return isEntryIdx(idx) ? idx - eshift : -1;
    }

    /** Activity-local index of an absolute index, or -1 if it is not an activity. */
    public int actOf(int idx) {
        return isActIdx(idx) ? idx - ashift : -1;
    }

    public boolean isHostIdx(int idx) {
        return idx >= hshift && idx < hshift + nhosts;
    }

    public boolean isTaskIdx(int idx) {
        return idx >= tshift && idx < tshift + ntasks;
    }

    public boolean isEntryIdx(int idx) {
        return idx >= eshift && idx < eshift + nentries;
    }

    public boolean isActIdx(int idx) {
        return idx >= ashift && idx < ashift + nacts;
    }

    /**
     * Checks the index-base invariants of this struct.
     *
     * The element index space is 0-based and unpadded, so the shifts must
     * partition 0..nidx-1 and every index-valued field must stay inside it.
     * Throws IllegalStateException on the first violation.
     */
    public void validateIndexSpace() {
        if (hshift != 0 || tshift != nhosts || eshift != nhosts + ntasks
                || ashift != nhosts + ntasks + nentries) {
            throw new IllegalStateException("LayeredNetworkStruct shifts do not partition 0.." + (nidx - 1)
                    + ": hshift=" + hshift + " tshift=" + tshift + " eshift=" + eshift + " ashift=" + ashift);
        }
        if (nidx != nhosts + ntasks + nentries + nacts) {
            throw new IllegalStateException("LayeredNetworkStruct nidx=" + nidx + " does not match its element counts");
        }
        if (type != null && type.length() != nidx) {
            throw new IllegalStateException("LayeredNetworkStruct type must be 1x" + nidx
                    + " (unpadded), got length " + type.length());
        }
        if (graph != null && (graph.getNumRows() != nidx || graph.getNumCols() != nidx)) {
            throw new IllegalStateException("LayeredNetworkStruct graph must be " + nidx + "x" + nidx
                    + " (unpadded), got " + graph.getNumRows() + "x" + graph.getNumCols());
        }
        if (parent != null) {
            for (int i = 0; i < Math.min(nidx, parent.length()); i++) {
                double p = parent.get(i);
                if (p < -1 || p >= nidx) {
                    throw new IllegalStateException("LayeredNetworkStruct parent[" + i + "]=" + p
                            + " is outside -1.." + (nidx - 1));
                }
            }
        }
    }



    public void print() {
        System.out.println("nidx: " + nidx);
        System.out.println("nhosts: " + nhosts);
        System.out.println("ntasks: " + ntasks);
        System.out.println("nentries: " + nentries);
        System.out.println("nacts: " + nacts);
        System.out.println("ncalls: " + ncalls);
        System.out.println("hshift: " + hshift);
        System.out.println("tshift: " + tshift);
        System.out.println("eshift: " + eshift);
        System.out.println("ashift: " + ashift);
        System.out.println("cshift: " + cshift);
        System.out.println("tasksof: " + tasksof);
        System.out.println("entriesof: " + entriesof);
        System.out.println("actof: " + actsof);
        System.out.println("callof: " + callsof);
        System.out.println("hostdem_type: " + hostdem_type);
        System.out.println("hostdem_mean: " + hostdem_mean);
        System.out.println("think_type: " + think_type);
        System.out.println("think_mean: " + think_mean);
        System.out.println("sched: " + sched);
        System.out.println("names: " + names);
        System.out.println("hashnames: " + hashnames);
        System.out.println("mult: ");
        mult.print();
        System.out.println("repl: ");
        repl.print();
        System.out.println("type: ");
        type.print();
        System.out.println("nitems: ");
        nitems.print();
        System.out.println("itemcap: " + itemcap);
        System.out.println("replacement: ");
        replacestrat.print();
        System.out.println("itemproc_type: " + itemproc_type);
        System.out.println("itemproc_mean: " + itemproc_mean);
        System.out.println("calltype: " + calltype);

        System.out.println("callpair: ");
        callpair.print();
        System.out.println("callproc_type: " + callproc_type);
        System.out.println("callproc_mean: " + callproc_mean);

        System.out.println("callnames: " + callnames);
        System.out.println("callhashname: " + callhashnames);
        System.out.println("actpretype: ");
        actpretype.print();
        System.out.println("actposttype: ");
        actposttype.print();
        System.out.println("graph: ");
        graph.printNonZero();
        System.out.println("parent: ");
        parent.print();
        System.out.println("replygraph: ");
        replygraph.printNonZero();
        System.out.println("iscache: ");
        iscache.print();
        System.out.println("iscaller: ");
        iscaller.printNonZero();
        System.out.println("issynccaller: ");
        issynccaller.printNonZero();
        System.out.println("isasynccaller: ");
        isasynccaller.printNonZero();
        System.out.println("isref: ");
        isref.print();
    }

}

