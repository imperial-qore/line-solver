package jline.examples;

import jline.lang.*;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.solvers.NetworkSolver;
import jline.solvers.mva.SolverMVA;

import java.util.Arrays;

/**
 * Miscellaneous examples
 */
public class Misc {

    public static Network ex1() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 1, node1, 0);

        node1.setService(jobclass1, Erlang.fitMeanAndOrder(3,2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5,3.0,10.0)); // (Delay,Class2)

        node2.setService(jobclass1, new HyperExp(0.1,1.0,10.0)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue1,Class2)

        node3.setService(jobclass1, new HyperExp(0.1,1.0,10.0)); // (Queue2,Class1)
        node3.setService(jobclass2, new Erlang(1,2)); // (Queue2,Class2)

        model.addLink(node1,node1);
        model.addLink(node1,node2);
        model.addLink(node1,node3);
        model.addLink(node2,node1);
        model.addLink(node3,node1);

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3));

        node1.setProbRouting(jobclass1, node2, 0.300000); // (Delay,Class1) -> (Queue1,Class1)
        node1.setProbRouting(jobclass1, node3, 0.700000); // (Delay,Class1) -> (Queue2,Class1)
        node2.setProbRouting(jobclass1, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        node3.setProbRouting(jobclass1, node1, 1.000000); // (Queue2,Class1) -> (Delay,Class1)

        node1.setRouting(jobclass2, RoutingStrategy.RAND); // (Delay,Class2) -> (Delay,Class2)
        node2.setRouting(jobclass2, RoutingStrategy.RAND); // (Delay,Class2) -> (Delay,Class2)
        node3.setRouting(jobclass2, RoutingStrategy.RAND); // (Delay,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        model.getStruct(false).rtnodes.print();

        return model;
    }
    public static Network ex2() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue node4 = new Queue(model, "Queue3", SchedStrategy.PS);
        Queue node5 = new Queue(model, "Queue4", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(0.333333)); // (Delay,Class1)
        node1.setService(jobclass2, Exp.fitMean(0.100000)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(0.100000)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.100000)); // (Queue1,Class2)
        node3.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(0.100000)); // (Queue2,Class2)
        node4.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue3,Class1)
        node4.setService(jobclass2, Exp.fitMean(0.100000)); // (Queue3,Class2)
        node5.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue4,Class1)
        node5.setService(jobclass2, Exp.fitMean(0.100000)); // (Queue4,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4, node5));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue1,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Queue2,Class1) -> (Queue3,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node5, 1.000000); // (Queue3,Class1) -> (Queue4,Class1)
        routingMatrix.set(jobclass1, jobclass1, node5, node1, 1.000000); // (Queue4,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue1,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node4, 1.000000); // (Queue2,Class2) -> (Queue3,Class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node5, 1.000000); // (Queue3,Class2) -> (Queue4,Class2)
        routingMatrix.set(jobclass2, jobclass2, node5, node1, 1.000000); // (Queue4,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }
    public static Network ex3() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.DPS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 1, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(0.333333)); // (Delay,Class1)
        node1.setService(jobclass2, Exp.fitMean(2.000000)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(10.000000)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue1,Class2)
        node3.setService(jobclass1, Exp.fitMean(10.000000)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue2,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.300000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node3, 0.700000); // (Delay,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000); // (Queue2,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 0.700000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node1, node3, 0.300000); // (Delay,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue1,Class2) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.000000); // (Queue2,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }
    public static Network ex4() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 1, node1, 0);

        node1.setService(jobclass1, APH.fitMeanAndSCV(0.666667,0.500000)); // (Delay,Class1)
        node1.setService(jobclass2, APH.fitMeanAndSCV(0.216667,1.579882)); // (Delay,Class2)
        node2.setService(jobclass1, APH.fitMeanAndSCV(0.190000,5.038781)); // (Queue1,Class1)
        node2.setService(jobclass2, APH.fitMeanAndSCV(0.700000,1.040816)); // (Queue1,Class2)
        node3.setService(jobclass1, APH.fitMeanAndSCV(0.190000,5.038781)); // (Queue2,Class1)
        node3.setService(jobclass2, APH.fitMeanAndSCV(2.000000,0.500000)); // (Queue2,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.300000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node3, 0.700000); // (Delay,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000); // (Queue2,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node1, 0.333333); // (Delay,Class2) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 0.333333); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node1, node3, 0.333333); // (Delay,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue1,Class2) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.000000); // (Queue2,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }


    public static void main(String[] args) throws Exception {
        Network model = ex1();

        NetworkStruct sn = model.getStruct(false);
        sn.rt.print();
        NetworkSolver solver = new SolverMVA(model);
        solver.getAvgTable();
    }
}
