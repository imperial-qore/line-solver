/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVA;
import jline.util.Maths;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

/**
 * Examples of models with load-dependent stations
 */
public class LoadDependentModel {
    /**
     * Basic load-dependent queue with FCFS scheduling.
     * <p>
     * Features:
     * - Closed network with 16 jobs and delay node
     * - FCFS queue with load-dependent service capacity
     * - Service capacity increases linearly up to 2 servers
     * - Alpha function: min(jobs+1, 2) servers available
     * - Demonstrates load-dependent server allocation
     *
     * @return configured load-dependent network model
     */
    public static Network ld_multiserver_fcfs() {
        int N = 16; // number of jobs
        int c = 2; // number of servers
        Network model = new Network("model");
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", N, node1, 0);
        node1.setService(jobclass1, Exp.fitMean(1.00));
        node2.setService(jobclass1, Exp.fitMean(1.500));
        Matrix alpha = new Matrix(1, N);
        for (int i = 0; i < N; i++) {
            alpha.set(0, i, Maths.min(i + 1, c));
        }
        node2.setLoadDependence(alpha);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        model.link(routingMatrix);
        return model;
    }

    /**
     * Multi-class load-dependent network with PS scheduling.
     * <p>
     * Features:
     * - Two closed classes: Class1 (4 jobs), Class2 (2 jobs)
     * - PS queue with total load-dependent capacity
     * - Service capacity based on total population across both classes
     * - Different service rates for each class
     * - Demonstrates multi-class load dependence
     *
     * @return configured multi-class load-dependent model
     */
    public static Network ld_multiserver_ps_twoclasses() {
        int N = 4; // number of jobs
        int c = 2; // number of servers
        Network model = new Network("model");
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", N, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", N / 2.0, node1, 0);
        node1.setService(jobclass1, Exp.fitMean(1.00));
        node1.setService(jobclass2, Exp.fitMean(2.00));
        node2.setService(jobclass1, Exp.fitMean(1.500));
        node2.setService(jobclass2, Exp.fitMean(2.500));
        Matrix alpha = new Matrix(1, (N + N / 2));
        for (int i = 0; i < (N + N / 2); i++) {
            alpha.set(0, i, Maths.min(i + 1, c));
        }
        node2.setLoadDependence(alpha);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.0);
        model.link(routingMatrix);
        return model;
    }

    /**
     * Three-station network with multiple load-dependent queues.
     * <p>
     * Features:
     * - Two closed classes with different populations
     * - Two PS queues (Queue1, Queue2) both with load dependence
     * - Serial routing: Delay → Queue1 → Queue2 → Delay
     * - Each queue has capacity up to 3 servers
     * - Different service rates at each station per class
     *
     * @return configured multi-station load-dependent model
     */
    public static Network ld_multiserver_ps() {
        int N = 4; // number of jobs
        int c = 3; // number of servers
        Network model = new Network("model");
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", N, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", N / 2.0, node1, 0);
        node1.setService(jobclass1, Exp.fitMean(1.00));
        node1.setService(jobclass2, Exp.fitMean(2.00));

        node2.setService(jobclass1, Exp.fitMean(1.500));
        node2.setService(jobclass2, Exp.fitMean(2.500));
        Matrix alpha = new Matrix(1, (N + N / 2));
        for (int i = 0; i < (N + N / 2); i++) {
            alpha.set(0, i, Maths.min(i + 1, c));
        }
        node2.setLoadDependence(alpha);

        node3.setService(jobclass1, Exp.fitMean(3.500));
        node3.setService(jobclass2, Exp.fitMean(4.500));
        alpha = new Matrix(1, (N + N / 2));
        for (int i = 0; i < (N + N / 2); i++) {
            alpha.set(0, i, Maths.min(i + 1, c));
        }
        node3.setLoadDependence(alpha);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.0);

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.0);
        model.link(routingMatrix);
        return model;
    }


    /**
     * Class-dependent service capacity model.
     * <p>
     * Features:
     * - Two closed classes with different populations
     * - PS queue with class-dependent capacity function
     * - Service capacity depends only on Class1 population
     * - Demonstrates selective load dependence by class
     * - Custom SerializableFunction for capacity calculation
     *
     * @return configured class-dependent network model
     */
    public static Network ld_class_dependence() {
        int N = 16;
        int c = 2;

        Network cdmodel = new Network("model");

        Delay node1 = new Delay(cdmodel, "Delay");
        Queue node2 = new Queue(cdmodel, "Queue1", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(cdmodel, "Class1", N, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(cdmodel, "Class2", N / 2.0, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.0));
        node1.setService(jobclass2, Exp.fitMean(2.0));

        node2.setService(jobclass1, Exp.fitMean(1.5));
        node2.setService(jobclass2, Exp.fitMean(2.5));

        SerializableFunction<Matrix, Double> beta =
                ni -> {
                    double class1Jobs = ni.get(0, 0); // Assuming ni is a 1x2 Matrix
                    return Math.min(class1Jobs, c);
                };
        node2.setClassDependence(beta);

        RoutingMatrix routingMatrix = cdmodel.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.0);

        cdmodel.link(routingMatrix);
        return cdmodel;
    }

    /**
     * Main method for testing and demonstrating load-dependent examples.
     *
     * <p>Currently configured to run ld_multiserver_ps_twoclasses() and solve it
     * using the MVA solver with default method settings.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        Network model = ld_multiserver_ps_twoclasses();
        SolverOptions options = new SolverOptions(SolverType.MVA);
        options.method = "default";
        MVA solver = new MVA(model, options);
        NetworkAvgTable t = solver.getAvgTable();
        t.print(options);
    }
}
