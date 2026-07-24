/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.Network;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

/**
 * Server-cluster models built with the {@link Network#cluster} static factories.
 *
 * <p>The topology is fixed: Source -> Dispatcher (Router) -> Server[1..M] -> Sink for the
 * open variants, and Think (Delay) -> Dispatcher -> Server[1..M] -> Think for the closed
 * variant. The factories collapse ~20 lines of node wiring into one call.
 */
public class ClusterModel {

    /**
     * Single-class open cluster with four PS servers and random dispatching.
     */
    public static Network cl_basic() {
        Matrix lambda = new Matrix(1, 1);
        lambda.set(0, 0, 0.4);
        Matrix D = new Matrix(4, 1);
        for (int i = 0; i < 4; i++) D.set(i, 0, 1.0);  // mean service time = 1
        return Network.clusterPs(lambda, D, RoutingStrategy.RAND);
    }

    /**
     * Heterogeneous open cluster: each of three FCFS servers has a different service time.
     */
    public static Network cl_heterogeneous() {
        Matrix lambda = new Matrix(1, 1);
        lambda.set(0, 0, 0.5);
        Matrix D = new Matrix(3, 1);
        D.set(0, 0, 0.5);  // fast server
        D.set(1, 0, 1.0);
        D.set(2, 0, 2.0);  // slow server
        Matrix S = new Matrix(3, 1, 3);
        S.fill(1.0);
        return Network.clusterFcfs(lambda, D, S, RoutingStrategy.JSQ);
    }

    /**
     * Closed cluster with three PS servers, single class, random dispatching.
     */
    public static Network cl_closed() {
        Matrix N = new Matrix(1, 1);
        N.set(0, 0, 8);
        Matrix Z = new Matrix(1, 1);
        Z.set(0, 0, 1.0);  // think time
        Matrix D = new Matrix(3, 1);
        for (int i = 0; i < 3; i++) D.set(i, 0, 1.0);
        Matrix S = new Matrix(3, 1, 3);
        S.fill(1.0);
        SchedStrategy[] sched = new SchedStrategy[3];
        for (int i = 0; i < 3; i++) sched[i] = SchedStrategy.PS;
        return Network.clusterClosed(N, Z, D, sched, S, RoutingStrategy.RAND);
    }

    /**
     * Two-class open cluster (e.g., interactive vs batch traffic).
     */
    public static Network cl_multiclass() {
        Matrix lambda = new Matrix(1, 2);
        lambda.set(0, 0, 0.3);
        lambda.set(0, 1, 0.2);
        Matrix D = new Matrix(2, 2);
        D.set(0, 0, 1.0); D.set(0, 1, 0.5);
        D.set(1, 0, 1.0); D.set(1, 1, 0.5);
        Matrix S = new Matrix(2, 1, 2);
        S.fill(1.0);
        return Network.clusterPs(lambda, D, S, RoutingStrategy.RAND);
    }
}
