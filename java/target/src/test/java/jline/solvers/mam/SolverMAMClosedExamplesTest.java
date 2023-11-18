package jline.solvers.mam;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.Erlang;
import jline.lang.distributions.Exp;
import jline.lang.distributions.HyperExp;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lib.M3A;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.Matrix;
import jline.solvers.NetworkAvgTable;

import java.util.Arrays;
import java.util.Collections;
import java.util.Map;

public class SolverMAMClosedExamplesTest {

    public static void main(String[] args) {

        Network model1 = example_closedModel_1();
        SolverMAM solver1 = new SolverMAM(model1);
        solver1.getAvgTable().print(new SolverOptions());
        System.out.println(solver1.result.runtime);

        Network model2 = example_closedModel_2();
        SolverMAM solver2 = new SolverMAM(model2);
        solver2.getAvgTable().print(new SolverOptions());
        System.out.println(solver2.result.runtime);

        Network model3 = example_closedModel_3();
        SolverMAM solver3 = new SolverMAM(model3);
        solver3.getAvgTable().print(new SolverOptions());
        System.out.println(solver3.result.runtime);

        Network model4 = example_closedModel_5();
        SolverMAM solver4 = new SolverMAM(model4);
        solver4.getAvgTable().print(new SolverOptions());
        System.out.println(solver4.result.runtime);

        Network model5 = example_closedModel_8();
        SolverMAM solver5 = new SolverMAM(model5);
        solver5.getAvgTable().print(new SolverOptions());
        System.out.println(solver5.result.runtime);

        Network model6 = example_closedModel_9();
        SolverMAM solver6 = new SolverMAM(model6);
        solver6.getAvgTable().print(new SolverOptions());
        System.out.println(solver6.result.runtime);

    }
    public static Network example_closedModel_1() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node2.setService(jobclass1, Exp.fitMean(1.500000)); // (Queue1,Class1)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Collections.singletonList(jobclass1),
                Arrays.asList(node1, node2));

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.700000); // (Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.300000); // (Delay,Class1) -> (Queue1,Class1)
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
        node1.setService(jobclass2, new HyperExp(0.5, 30.0, 100.0)); // (Delay,Class2)
        node2.setService(jobclass1, new HyperExp(0.1, 10.0, 100.0)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(1.0)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2));

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
        node1.setService(jobclass2, new HyperExp(0.5,3.0,10.0)); // (Delay,Class2)
        node1.setService(jobclass3, new Exp(1)); // (Delay,Class3)
        node2.setService(jobclass1, new HyperExp(0.1,1.0,10.0)); // (Queue1,Class1)
        node2.setService(jobclass2, new Exp(2)); // (Queue1,Class2)
        node2.setService(jobclass3, new Exp(3)); // (Queue1,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2, jobclass3),
                Arrays.asList(node1, node2));

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.3); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.1); // (Delay,Class1) -> (CS_Delay_to_Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 0.2); // (Queue1,Class1) -> (CS_Queue1_to_Delay,Class1)
        routingMatrix.set(jobclass1, jobclass2, node1, node1, 0.6); // (CS_Delay_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node1, 0.8); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass1, node2, node1, 1); // (CS_Delay_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.set(jobclass3, jobclass3, node2, node1, 1.000000); // (CS_Queue1_to_Delay,Class2) -> (Delay,Class1)
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_9() {
        Network model = new Network("Model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 10, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(1.500000)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.500000)); // (Queue1,Class2)
        node3.setService(jobclass1, Exp.fitMean(2.000000)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(2.000000)); // (Queue2,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3));

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.700000); // (Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.300000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000); // (Queue2,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node1, 0.700000); // (Delay,Class2) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass2, node1, node3, 0.300000); // (Delay,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue1,Class2) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.000000); // (Queue2,Class2) -> (Delay,Class2)

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

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Delay,Class2)
        node2.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.100000)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue1,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_closedModel_5() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay1");
        Delay node2 = new Delay(model, "Delay2");
        Queue node3 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node4 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(91.000000)); // (Delay1,Class1)
        node1.setService(jobclass2, Exp.fitMean(92.000000)); // (Delay1,Class2)
        node2.setService(jobclass1, Exp.fitMean(93.000000)); // (Delay2,Class1)
        node2.setService(jobclass2, Exp.fitMean(94.000000)); // (Delay2,Class2)
        node3.setService(jobclass1, Exp.fitMean(10.000000)); // (Queue1,Class1)
        node3.setService(jobclass2, Exp.fitMean(5.000000)); // (Queue1,Class2)
        node4.setService(jobclass1, Exp.fitMean(5.000000)); // (Queue2,Class1)
        node4.setService(jobclass2, Exp.fitMean(9.000000)); // (Queue2,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay1,Class1) -> (Delay2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Delay2,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Queue1,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.000000); // (Queue2,Class1) -> (Delay1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay1,Class2) -> (Delay2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Delay2,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node4, 1.000000); // (Queue1,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node1, 1.000000); // (Queue2,Class2) -> (Delay1,Class2)

        model.link(routingMatrix);

        return model;
    }
}



