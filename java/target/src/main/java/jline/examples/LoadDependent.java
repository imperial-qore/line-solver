package jline.examples;

import jline.util.Maths;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.Exp;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.Matrix;
import jline.solvers.NetworkAvgTable;

import java.util.Arrays;
import java.util.Collections;

/**
 * Examples of models with load-dependent stations
 */
public class LoadDependent {
    public static Network ex1(){
        int N = 16; // number of jobs
        int c = 2; // number of servers
        Network model = new Network("model");
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", N, node1, 0);
        node1.setService(jobclass1, Exp.fitMean(1.000000));
        node2.setService(jobclass1, Exp.fitMean(1.500000));
        Matrix alpha = new Matrix(1,N);
        for(int i = 0; i < N; i++){
            alpha.set(0, i, Maths.min(i+1, c));
        }
        node2.setLoadDependence(alpha);
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Collections.singletonList(jobclass1),
                Arrays.asList(node1, node2));
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        model.link(routingMatrix);
        return model;
    }

    public static Network ex2(){
        int N = 4; // number of jobs
        int c = 2; // number of servers
        Network model = new Network("model");
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", N, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", N/2, node1, 0);
        node1.setService(jobclass1, Exp.fitMean(1.000000));
        node1.setService(jobclass2, Exp.fitMean(2.000000));
        node2.setService(jobclass1, Exp.fitMean(1.500000));
        node2.setService(jobclass2, Exp.fitMean(2.500000));
        Matrix alpha = new Matrix(1,(N+N/2));
        for(int i = 0; i < (N+N/2); i++){
            alpha.set(0, i, Maths.min(i+1, c));
        }
        node2.setLoadDependence(alpha);
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2));
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.0);
        model.link(routingMatrix);
        return model;
    }

    public static Network ex3(){
        int N = 4; // number of jobs
        int c = 3; // number of servers
        Network model = new Network("model");
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", N, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", N/2, node1, 0);
        node1.setService(jobclass1, Exp.fitMean(1.000000));
        node1.setService(jobclass2, Exp.fitMean(2.000000));

        node2.setService(jobclass1, Exp.fitMean(1.500000));
        node2.setService(jobclass2, Exp.fitMean(2.500000));
        Matrix alpha = new Matrix(1,(N+N/2));
        for(int i = 0; i < (N+N/2); i++){
            alpha.set(0, i, Maths.min(i+1, c));
        }
        node2.setLoadDependence(alpha);

        node3.setService(jobclass1, Exp.fitMean(3.500000));
        node3.setService(jobclass2, Exp.fitMean(4.500000));
        alpha = new Matrix(1,(N+N/2));
        for(int i = 0; i < (N+N/2); i++){
            alpha.set(0, i, Maths.min(i+1, c));
        }
        node3.setLoadDependence(alpha);
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3));
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.0);

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.0);
        model.link(routingMatrix);
        return model;
    }

    public static void main(String[] args) {
        Network model = ex2();
        SolverOptions options = new SolverOptions(SolverType.MVA);
        options.method = "default";
        NetworkSolver solver = new SolverMVA(model, options);
        NetworkAvgTable t = solver.getAvgTable();
        t.print(options);
    }
}
