package jline.examples;

import jline.lang.*;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Exp;
import jline.lang.nodes.*;
import jline.solvers.NetworkSolver;
import jline.solvers.mva.SolverMVA;
import jline.util.Matrix;

import java.util.Map;

/**
 * Examples of state probability computations
 */
public class StateProbabilities {

    public static Network example_stateProbabilities_1() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        node3.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.250000)); // (Queue1,Class2)
        node3.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(0.333333)); // (Queue2,Class2)

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
    public static Network example_stateProbabilities_2() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        node3.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 4, node1, 0);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 0, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node1.setService(jobclass2, Exp.fitMean(0.500000)); // (Delay,Class2)
        node1.setService(jobclass3, Exp.fitMean(1.000000)); // (Delay,Class3)
        node1.setService(jobclass4, Exp.fitMean(1.000000)); // (Delay,Class4)
        node2.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.250000)); // (Queue1,Class2)
        node2.setService(jobclass3, Exp.fitMean(0.200000)); // (Queue1,Class3)
        node2.setService(jobclass4, Exp.fitMean(1.000000)); // (Queue1,Class4)
        node3.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(0.333333)); // (Queue2,Class2)
        node3.setService(jobclass3, Exp.fitMean(0.200000)); // (Queue2,Class3)
        node3.setService(jobclass4, Exp.fitMean(0.500000)); // (Queue2,Class4)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        routingMatrix.set(jobclass1, jobclass2, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        routingMatrix.set(jobclass1, jobclass3, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass1, jobclass4, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));

        routingMatrix.set(jobclass2, jobclass1, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        routingMatrix.set(jobclass2, jobclass2, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        routingMatrix.set(jobclass2, jobclass3, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass2, jobclass4, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));

        routingMatrix.set(jobclass3, jobclass1, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass3, jobclass2, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass3, jobclass3, new Matrix("[0,1,0; 0,0,1; 0,0,0]"));
        routingMatrix.set(jobclass3, jobclass4, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));

        routingMatrix.set(jobclass4, jobclass1, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass4, jobclass2, new Matrix("[0,0,0; 0,0,0; 0,0,0]"));
        routingMatrix.set(jobclass4, jobclass3, new Matrix("[0,0,0; 0,0,0; 1,0,0]"));
        routingMatrix.set(jobclass4, jobclass4, new Matrix("[0,0,1; 0,0,0; 0,0,0]"));

        model.link(routingMatrix);

        return model;
    }

    public static Network example_stateProbabilities_3() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        node3.setNumberOfServers(2);
        Router node4 = new Router(model, "CS_Queue2_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 3, node1, 0);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 0, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node1.setService(jobclass2, Exp.fitMean(0.500000)); // (Delay,Class2)
        node1.setService(jobclass3, Exp.fitMean(1.000000)); // (Delay,Class3)
        node1.setService(jobclass4, Exp.fitMean(1.000000)); // (Delay,Class4)
        node2.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.250000)); // (Queue1,Class2)
        node2.setService(jobclass3, Exp.fitMean(0.200000)); // (Queue1,Class3)
        node2.setService(jobclass4, Exp.fitMean(1.000000)); // (Queue1,Class4)
        node3.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(0.333333)); // (Queue2,Class2)
        node3.setService(jobclass3, Exp.fitMean(0.200000)); // (Queue2,Class3)
        node3.setService(jobclass4, Exp.fitMean(0.500000)); // (Queue2,Class4)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue1,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Queue2,Class1) -> (CS_Queue2_to_Delay,Class1)
        routingMatrix.set(jobclass1, jobclass2, node4, node1, 1.000000); // (CS_Queue2_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass1, node4, node1, 1.000000); // (CS_Queue2_to_Delay,Class2) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue1,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node4, 1.000000); // (Queue2,Class2) -> (CS_Queue2_to_Delay,Class2)
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.000000); // (Delay,Class3) -> (Queue1,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node3, 1.000000); // (Queue1,Class3) -> (Queue2,Class3)
        routingMatrix.set(jobclass3, jobclass3, node3, node4, 1.000000); // (Queue2,Class3) -> (CS_Queue2_to_Delay,Class3)
        routingMatrix.set(jobclass3, jobclass4, node4, node1, 1.000000); // (CS_Queue2_to_Delay,Class3) -> (Delay,Class4)
        routingMatrix.set(jobclass4, jobclass3, node4, node1, 1.000000); // (CS_Queue2_to_Delay,Class4) -> (Delay,Class3)
        routingMatrix.set(jobclass4, jobclass4, node1, node3, 1.000000); // (Delay,Class4) -> (Queue2,Class4)
        routingMatrix.set(jobclass4, jobclass4, node2, node3, 1.000000); // (Queue1,Class4) -> (Queue2,Class4)
        routingMatrix.set(jobclass4, jobclass4, node3, node4, 1.000000); // (Queue2,Class4) -> (CS_Queue2_to_Delay,Class4)

        model.link(routingMatrix);

        return model;
    }
    public static Network example_stateProbabilities_4() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.PS);
        node3.setNumberOfServers(3);
        Router node4 = new Router(model, "CS_Queue3_to_Queue1"); // Dummy node, class switching is embedded in the routing matrix P

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 1, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 1, node1, 0);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 1, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue1,Class1)
        node1.setService(jobclass2, Exp.fitMean(0.500000)); // (Queue1,Class2)
        node1.setService(jobclass3, Exp.fitMean(1.000000)); // (Queue1,Class3)
        node1.setService(jobclass4, Exp.fitMean(1.000000)); // (Queue1,Class4)
        node2.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue2,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.250000)); // (Queue2,Class2)
        node2.setService(jobclass3, Exp.fitMean(0.200000)); // (Queue2,Class3)
        node2.setService(jobclass4, Exp.fitMean(1.000000)); // (Queue2,Class4)
        node3.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue3,Class1)
        node3.setService(jobclass2, Exp.fitMean(0.333333)); // (Queue3,Class2)
        node3.setService(jobclass3, Exp.fitMean(0.200000)); // (Queue3,Class3)
        node3.setService(jobclass4, Exp.fitMean(0.500000)); // (Queue3,Class4)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Queue1,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue2,Class1) -> (Queue3,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Queue3,Class1) -> (CS_Queue3_to_Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass2, node4, node1, 1.000000); // (CS_Queue3_to_Queue1,Class1) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass1, node4, node1, 1.000000); // (CS_Queue3_to_Queue1,Class2) -> (Queue1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Queue1,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue2,Class2) -> (Queue3,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node4, 1.000000); // (Queue3,Class2) -> (CS_Queue3_to_Queue1,Class2)
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.000000); // (Queue1,Class3) -> (Queue2,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node3, 1.000000); // (Queue2,Class3) -> (Queue3,Class3)
        routingMatrix.set(jobclass3, jobclass3, node3, node4, 1.000000); // (Queue3,Class3) -> (CS_Queue3_to_Queue1,Class3)
        routingMatrix.set(jobclass3, jobclass4, node4, node1, 1.000000); // (CS_Queue3_to_Queue1,Class3) -> (Queue1,Class4)
        routingMatrix.set(jobclass4, jobclass3, node4, node1, 1.000000); // (CS_Queue3_to_Queue1,Class4) -> (Queue1,Class3)
        routingMatrix.set(jobclass4, jobclass4, node1, node3, 1.000000); // (Queue1,Class4) -> (Queue3,Class4)
        routingMatrix.set(jobclass4, jobclass4, node2, node3, 1.000000); // (Queue2,Class4) -> (Queue3,Class4)
        routingMatrix.set(jobclass4, jobclass4, node3, node4, 1.000000); // (Queue3,Class4) -> (CS_Queue3_to_Queue1,Class4)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_stateProbabilities_5() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        node3.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.250000)); // (Queue1,Class2)
        node3.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(0.333333)); // (Queue2,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000);
        routingMatrix.set(jobclass1, jobclass2, node2, node3, 1.000000);
        routingMatrix.set(jobclass2, jobclass1, node1, node2, 1.000000);
        routingMatrix.set(jobclass2, jobclass1, node3, node1, 1.000000);
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000);

        model.link(routingMatrix);

        return model;
    }

    public static Network example_stateProbabilities_6() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        node3.setNumberOfServers(2);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.250000)); // (Queue1,Class2)
        node3.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(0.333333)); // (Queue2,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000);
        routingMatrix.set(jobclass1, jobclass2, node2, node3, 1.000000);
        routingMatrix.set(jobclass2, jobclass1, node1, node2, 1.000000);
        routingMatrix.set(jobclass2, jobclass1, node3, node1, 1.000000);
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000);

        model.link(routingMatrix);

        return model;
    }

    public static void main(String[] args) throws Exception {
        Network model = example_stateProbabilities_2();
    }
}
