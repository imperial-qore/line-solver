package jline.examples;

import jline.lang.*;
import jline.lang.constant.RoutingStrategy;
import jline.lang.distributions.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.util.Maths;
import jline.util.Matrix;

import java.util.Arrays;

/**
 * Examples of closed queueing networks
 */
public class ClosedModel {

    public static Network example_closedModel_1() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node2.setService(jobclass1, Exp.fitMean(1.500000)); // (Queue1,Class1)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.700000); // (Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.3); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_2() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,Class2)
        node2.setService(jobclass1, new HyperExp(0.1, 1.0, 10.0)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(1.0)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.3);
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.1);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 0.2);

        routingMatrix.set(jobclass1, jobclass2, node1, node1, 0.6);
        routingMatrix.set(jobclass1, jobclass2, node2, node1, 0.8);

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);

        routingMatrix.set(jobclass2, jobclass1, node2, node1, 1.0);

        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_3() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        node2.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 1, node1, 0);

        node1.setService(jobclass1, new Erlang(3,2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,Class2)
        node1.setService(jobclass3, Exp.fitMean(1.000000)); // (Delay,Class3)

        node2.setService(jobclass1, new HyperExp(0.1, 1.0,10.0)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(2.0)); // (Queue1,Class2)
        node2.setService(jobclass3, new Exp(3.0)); // (Queue1,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.3);
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.1);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 0.2);

        routingMatrix.set(jobclass1, jobclass2, node1, node1, 0.6);
        routingMatrix.set(jobclass1, jobclass2, node2, node1, 0.8);

        routingMatrix.set(jobclass2, jobclass1, node2, node1, 1);

        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1);

        Matrix c = Maths.circul(model.getNumberOfStations());

        routingMatrix.set(jobclass3, jobclass3, node1, node1, c.get(0,0));
        routingMatrix.set(jobclass3, jobclass3, node1, node2, c.get(0,1));
        routingMatrix.set(jobclass3, jobclass3, node2, node1, c.get(1,0));
        routingMatrix.set(jobclass3, jobclass3, node2, node2, c.get(1,1));

        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_4() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        node2.setNumberOfServers(3);
        node3.setNumberOfServers(3);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 2, node1, 0);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 1, node1, 0);

        node1.setService(jobclass1, new Exp(1)); // (Delay,Class1)
        node1.setService(jobclass2, new Exp(1)); // (Delay,Class2)
        node1.setService(jobclass3, new Exp(10)); // (Delay,Class3)
        node1.setService(jobclass4, new Exp(1)); // (Delay,Class4)

        node2.setService(jobclass1, new Exp(1)); // (Queue1,Class1)
        node2.setService(jobclass2, new Erlang(1, 2)); // (Queue1,Class2)
        node2.setService(jobclass3, new Exp(10)); // (Queue1,Class3)
        node2.setService(jobclass4, new Exp(1)); // (Queue1,Class4)

        node3.setService(jobclass1, Disabled.getInstance()); // (Queue2,Class1)
        node3.setService(jobclass2, Disabled.getInstance()); // (Queue2,Class2)
        node3.setService(jobclass3, new Erlang(1, 2)); // (Queue2,Class3)
        node3.setService(jobclass4, new Exp(1)); // (Queue2,Class4)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.500000);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000);

        routingMatrix.set(jobclass1, jobclass2, node1, node2, 0.500000);

        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000);
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.000000);

        routingMatrix.set(jobclass2, jobclass1, node1, node2, 1.000000);

        routingMatrix.set(jobclass3, jobclass3, node1, node2, 0.250000);
        routingMatrix.set(jobclass3, jobclass3, node1, node3, 0.250000);
        routingMatrix.set(jobclass3, jobclass3, node2, node1, 1.000000);
        routingMatrix.set(jobclass3, jobclass3, node3, node1, 1.000000);

        routingMatrix.set(jobclass3, jobclass4, node1, node2, 0.500000);

        routingMatrix.set(jobclass4, jobclass4, node2, node1, 1.000000);
        routingMatrix.set(jobclass4, jobclass4, node3, node1, 1.000000);

        routingMatrix.set(jobclass4, jobclass3, node1, node2, 1.000000);

        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_6() {
        Network model = new Network("Model");

        Matrix csMatrix = new Matrix("[0,1;1,0]");

        ClassSwitch node1 = new ClassSwitch(model, "CS", csMatrix);
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue node4 = new Queue(model, "Delay", SchedStrategy.INF);
        node4.setNumberOfServers(Integer.MAX_VALUE);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 15, node4, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 5, node4, 0);

        node2.setService(jobclass1, Exp.fitMean(1.500000)); // (Queue1,Class1)
        node2.setService(jobclass2, Erlang.fitMeanAndOrder(1.5, 2)); // (Queue1,Class2)

        node3.setService(jobclass1, Erlang.fitMeanAndOrder(1.5, 2)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(1.500000)); // (Queue2,Class2)
        node4.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node4.setService(jobclass2, Exp.fitMean(1.000000)); // (Delay,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node4, node2, 0.500000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node3, 0.500000); // (Delay,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass2, node1, node4, 1.000000); // (CS,Class1) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass1, node1, node4, 1.000000); // (CS,Class2) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node4, node2, 0.500000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node3, 0.500000); // (Delay,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (CS,Class1)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue1,Class2) -> (CS,Class2)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000); // (Queue2,Class1) -> (CS,Class1)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.000000); // (Queue2,Class2) -> (CS,Class2)

        node1.setRouting(jobclass1, RoutingStrategy.RAND);
        node2.setRouting(jobclass1, RoutingStrategy.RAND);
        node3.setRouting(jobclass1, RoutingStrategy.RAND);
        node4.setRouting(jobclass1, RoutingStrategy.RROBIN);

        node1.setRouting(jobclass2, RoutingStrategy.RAND);
        node2.setRouting(jobclass2, RoutingStrategy.RAND);
        node3.setRouting(jobclass2, RoutingStrategy.RAND);
        node4.setRouting(jobclass2, RoutingStrategy.RROBIN, node2, 1);
        node4.setRouting(jobclass2, RoutingStrategy.RROBIN, node3, 2);

        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_7ps() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,Class2)
        node2.setService(jobclass1, new Exp(1)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(1)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue1,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_7fcfs() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,Class1)
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue1,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_7lcfspr() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.LCFSPR);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, APH.fitMeanAndSCV(0.666667,0.500000)); // (Delay,Class1)
        node1.setService(jobclass2, APH.fitMeanAndSCV(0.216667,1.579882)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue1,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_8() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        node2.setNumberOfServers(3);
        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 4, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, new Exp(1)); // (Delay,Class1)
        node1.setService(jobclass2, new Exp(1));

        node2.setService(jobclass1, new Exp(1)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(10));


        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue1,Class2) -> (Delay,Class2)
        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_9() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 10, node1, 0);


        node1.setService(jobclass1, Exp.fitMean(1.000000));
        node2.setService(jobclass1, Exp.fitMean(1.500000));
        node3.setService(jobclass1, Exp.fitMean(3.000000));
        node1.setService(jobclass2, Exp.fitMean(1.000000));
        node2.setService(jobclass2, Exp.fitMean(1.500000));
        node3.setService(jobclass2, Exp.fitMean(3.000000));

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue1,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000); // (Queue2,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue1,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.000000); // (Queue2,Class2) -> (Delay,Class2)
        model.link(routingMatrix);

        return model;
    }

    public static void main(String[] args) throws Exception {
        Network model = example_closedModel_4();
        new SolverMVA(model).getAvgTable().print();
        new SolverNC(model).getAvgTable().print();
    }

}
