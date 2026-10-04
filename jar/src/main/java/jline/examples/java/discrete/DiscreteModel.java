/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.discrete;

import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Geometric;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Discrete-time (slotted) model builders.
 *
 * <p>Time advances one slot at a time; in each slot a job in service completes
 * with probability p and an arrival occurs with probability b, both recorded at
 * the end of the slot with the departure resolved before the arrival (Daduna's
 * LA rule and D/A rule). Every rate is therefore a per-slot probability and
 * every time is a number of slots.</p>
 *
 * <p>Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS
 * 2046, Springer, 2001.</p>
 */
public class DiscreteModel {

    /** Geo/Geo/1 with an unbounded buffer (theorem 2.3, corollary 2.7). */
    public static Network dt_geogeo1(double a, double s) {
        Network model = new Network("GeoGeo1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Geometric(a));
        queue.setService(jobClass, new Geometric(s));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /** Geo/Geo/1/L loss system (corollary 2.8). */
    public static Network dt_geogeo1_loss(double a, double s, int L) {
        Network model = new Network("GeoGeo1L");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Geometric(a));
        queue.setService(jobClass, new Geometric(s));
        queue.setCap(L);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /**
     * Load-dependent Bernoulli server, p(n) = s min(n,c), the discrete-time
     * multiserver approximation of example 2.10.
     */
    public static Network dt_bernoulli_loaddep(double a, double s, int c, int L) {
        Network model = new Network("LoadDepBernoulli");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Geometric(a));
        queue.setService(jobClass, new Geometric(s));
        queue.setCap(L);
        Matrix alpha = new Matrix(1, L);
        for (int n = 1; n <= L; n++) {
            alpha.set(0, n - 1, Math.min(n, c));
        }
        queue.setLoadDependence(alpha);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /** Closed cycle of Bernoulli servers (corollary 3.4). */
    public static Network dt_cycle(double[] p, int N) {
        return cycle(p, N, -1, null);
    }

    /**
     * Closed cycle whose station {@code ldStation} runs at
     * p(n) = p min(n,ldServers) (theorem 3.2).
     */
    public static Network dt_cycle_loaddep(double[] p, int N, int ldStation, int ldServers) {
        return cycle(p, N, ldStation, buildAlpha(N, ldServers));
    }

    private static Matrix buildAlpha(int N, int servers) {
        Matrix alpha = new Matrix(1, N);
        for (int n = 1; n <= N; n++) {
            alpha.set(0, n - 1, Math.min(n, servers));
        }
        return alpha;
    }

    private static Network cycle(double[] p, int N, int ldStation, Matrix alpha) {
        Network model = new Network(ldStation < 0 ? "BernoulliCycle" : "BernoulliCycleLD");
        List<Queue> station = new ArrayList<Queue>();
        for (int j = 0; j < p.length; j++) {
            station.add(new Queue(model, "Queue" + (j + 1), SchedStrategy.FCFS));
        }
        ClosedClass jobClass = new ClosedClass(model, "Jobs", N, station.get(0));
        for (int j = 0; j < p.length; j++) {
            station.get(j).setService(jobClass, new Geometric(p[j]));
        }
        if (ldStation >= 0) {
            station.get(ldStation).setLoadDependence(alpha);
        }
        model.link(Network.serialRouting(new ArrayList<jline.lang.nodes.Node>(station).toArray(new jline.lang.nodes.Node[0])));
        return model;
    }

    /**
     * Multichain closed cycle (section 3.2). Both chains follow the same cycle
     * and never switch class, so the joint queue length law is the unichain one
     * at the aggregate population.
     */
    public static Network dt_cycle_multiclass(double[] p, int[] pops) {
        Network model = new Network("BernoulliCycleMC");
        int J = p.length;
        List<Queue> station = new ArrayList<Queue>();
        for (int j = 0; j < J; j++) {
            station.add(new Queue(model, "Queue" + (j + 1), SchedStrategy.FCFS));
        }
        List<JobClass> classes = new ArrayList<JobClass>();
        for (int g = 0; g < pops.length; g++) {
            ClosedClass cls = new ClosedClass(model, "Chain" + (g + 1), pops[g], station.get(0));
            classes.add(cls);
            for (int j = 0; j < J; j++) {
                station.get(j).setService(cls, new Geometric(p[j]));
            }
        }
        RoutingMatrix P = model.initRoutingMatrix();
        for (int g = 0; g < classes.size(); g++) {
            for (int j = 0; j < J; j++) {
                P.set(classes.get(g), classes.get(g), station.get(j), station.get((j + 1) % J), 1.0);
            }
        }
        model.link(P);
        return model;
    }

    /** Service probabilities used by the cycle examples. */
    public static double[] defaultCycleRates() {
        return Arrays.copyOf(new double[]{0.5, 0.25, 0.7}, 3);
    }
}
