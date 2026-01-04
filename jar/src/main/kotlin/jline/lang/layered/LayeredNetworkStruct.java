/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.Copyable;
import jline.lang.constant.CallType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.DiscreteDistribution;
import jline.lang.processes.Distribution;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Class summarizing the characteristics of a LayeredNetwork
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
    public Map<Integer, DiscreteDistribution> callproc;
    public Map<Integer, ProcessType> callproc_type;
    public Map<Integer, Matrix> callproc_params;
    public Map<Integer, Double> callproc_mean;
    public Map<Integer, Double> callproc_scv;
    public Map<Integer, MatrixCell> callproc_proc;

    public Map<Integer, String> callnames;
    public Map<Integer, String> callhashnames;

    public Map<Integer, SchedStrategy> sched;
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
    public Matrix parent;

    // Matrix representation of scheduling strategies (for Python compatibility)
    public Matrix schedid;

    public Matrix iscaller;
    public Matrix issynccaller;
    public Matrix isasynccaller;
    public Matrix callpair;
    public Matrix taskgraph;

    public Matrix actpretype;
    public Matrix actposttype;

    public Matrix isref;

    public Matrix isfunction;

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

    public Matrix conntasks; // 1 * n matrix
    
    public List<Integer> hitmissaidx;
    public Integer hitaidx;
    public Integer missaidx;

    public LayeredNetworkStruct() {

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

