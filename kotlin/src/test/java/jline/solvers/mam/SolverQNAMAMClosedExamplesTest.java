package jline.solvers.mam;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.APH;
import jline.lang.distributions.Exp;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.NetworkAvgTable;

import java.util.Arrays;

public class SolverQNAMAMClosedExamplesTest {
    public static void main(String[] args) {
        Network model1 = example_closedModel_2();


        SolverOptions options = new SolverOptions(SolverType.MAM);
        options.method = "qnamam";

        NetworkSolver solver1 = new SolverMAM(model1, options);
        NetworkAvgTable t1 = solver1.getAvgTable();
        t1.print(options);
    }

    public static Network example_closedModel_1() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model,"Queue2", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 4, node1, 0);

        node1.setService(jobclass1, APH.fitMeanAndSCV(1,1.0/4)); // (Delay,Class1)
        node2.setService(jobclass1, APH.fitMeanAndSCV(2,1.0/5)); // (Queue1,Class1)
        node3.setService(jobclass1, APH.fitMeanAndSCV(4,1.0/6));
        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1),
                Arrays.asList(node1, node2,node3));


        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1);
        model.link(routingMatrix);

        return model;
    }
    public static Network example_closedModel_2() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model,"Queue2", SchedStrategy.FCFS);
        Queue node4 = new Queue(model,"Queue3", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 6, node1, 0);

        node1.setService(jobclass1, APH.fitMeanAndSCV(1,1.0/2)); // (Delay,Class1)
        node2.setService(jobclass1, APH.fitMeanAndSCV(2,1.0/4)); // (Queue1,Class1)
        node3.setService(jobclass1, APH.fitMeanAndSCV(3,1.0/6));
        node4.setService(jobclass1, APH.fitMeanAndSCV(4,1.0/8));
        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1),
                Arrays.asList(node1, node2,node3,node4));


        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 0.5); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node4, 0.5);
        routingMatrix.set(jobclass1, jobclass1, node3, node2, 0.5);
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 0.5);
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1);
        model.link(routingMatrix);

        return model;
    }
}
