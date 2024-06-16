package jline.lang.layered;

import jline.lang.constant.CallType;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.DiscreteDistribution;
import jline.lang.distributions.Distribution;
import jline.lang.distributions.Geometric;
import jline.util.Matrix;

import java.util.List;
import java.util.Map;

/**
 * Class summarizing the characteristics of a LayeredNetwork
 */
public class LayeredNetworkStruct {
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
    public Map<Integer, Distribution> hostdem;
    public Map<Integer, Distribution> think;
    public Map<Integer, Integer> itemcap;
    public Map<Integer, DiscreteDistribution> itemproc;
    public Map<Integer, CallType> calltype;
    public Map<Integer, Geometric> callproc;
    public Map<Integer, String> callnames;
    public Map<Integer, String> callhashnames;

    public Map<Integer, SchedStrategy> sched;
    public Matrix schedid;
    public Map<Integer, String> names;
    public Map<Integer, String> hashnames;
    public Matrix mult;
    public Matrix repl;
    public Matrix type;
    public Matrix graph;

    public Matrix replygraph;
    public Matrix nitems;

    public Matrix replacement;
    public Matrix iscache;
    public Matrix parent;

    public Matrix iscaller;
    public Matrix issynccaller;
    public Matrix isasynccaller;
    public Matrix callpair;
    public Matrix taskgraph;

    public Matrix actpretype;
    public Matrix actposttype;

    public Matrix isref;

    public Matrix conntasks; // 1 * n matrix

    public LayeredNetworkStruct () {

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
        System.out.println("hostdem: " + hostdem);
        System.out.println("think: " + think);
        System.out.println("sched: " + sched);
        System.out.println("schedid: ");
        schedid.print();
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
        replacement.print();
        System.out.println("itemproc: " + itemproc);
        System.out.println("calltype: " + calltype);

        System.out.println("callpair: ");
        callpair.print();
        System.out.println("callproc: " + callproc);

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
