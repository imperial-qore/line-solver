package jline.solvers.mva;

import jline.examples.java.basic.ClosedModel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.*;
import jline.lang.processes.APH;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class SolverMVATestFixtures {

    public static Network test_qna_model1() {
        Network model = new Network("model");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");

        OpenClass openClass = new OpenClass(model, "Class1", 0);
        delay.setService(openClass, new HyperExp(0.5, 3.0, 10.0));
        queue.setService(openClass, new Exp(1));
        source.setArrival(openClass, new Exp(0.1));

        RoutingMatrix p = new RoutingMatrix(model, model.getJobClasses(), model.getNodes());
        p.addConnection(source, delay, openClass, 1);
        p.addConnection(delay, queue, openClass, 1);
        p.addConnection(queue, sink, openClass, 1);

        model.link(p);

        return model;
    }

    public static Network test_qna_model5() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");
        Router node4 = new Router(model, "VSink1");
        Router node5 = new Router(model, "VSink2");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class1", 0);

        node1.setArrival(jobclass1, Exp.fitMean(1.00)); // (Source,Class1)
        node2.setService(jobclass1, Exp.fitMean(0.010000)); // (Queue1,Class1)

        node1.setArrival(jobclass2, new Exp(1.0));
        node2.setService(jobclass2, Exp.fitMean(0.010000));

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4, node5));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node4, 0.600); // (Queue1,Class1) -> (VSink1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node5, 0.400); // (Queue1,Class1) -> (VSink2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node3, 1.00); // (VSink1,Class1) -> (Sink,Class1)
        routingMatrix.set(jobclass1, jobclass1, node5, node3, 1.00); // (VSink2,Class1) -> (Sink,Class1)


        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00); // (Source,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node2, node4, 0.100); // (Queue1,Class1) -> (VSink1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node2, node5, 0.900); // (Queue1,Class1) -> (VSink2,Class1)
        routingMatrix.set(jobclass2, jobclass2, node4, node3, 1.00); // (VSink1,Class1) -> (Sink,Class1)
        routingMatrix.set(jobclass2, jobclass2, node5, node3, 1.00); // (VSink2,Class1) -> (Sink,Class1)

        model.link(routingMatrix);

        return model;
    }

    public static Network test_qna_model6() {
        Network model = new Network("model2");
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue node4 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue node5 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Sink node6 = new Sink(model, "Sink");

        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);
        OpenClass jobclass3 = new OpenClass(model, "Class3", 0);

        node1.setArrival(jobclass1, Exp.fitMean(5));
        node1.setArrival(jobclass2, Exp.fitMean(8));
        node1.setArrival(jobclass3, Exp.fitMean(7));
        node2.setService(jobclass1, Exp.fitMean(0.3));
        node2.setService(jobclass2, Exp.fitMean(0.5));
        node2.setService(jobclass3, Exp.fitMean(0.6));
        node3.setService(jobclass1, Exp.fitMean(1.1));
        node3.setService(jobclass2, Exp.fitMean(1.3));
        node3.setService(jobclass3, Exp.fitMean(1.5));
        node4.setService(jobclass1, Exp.fitMean(2.0));
        node4.setService(jobclass2, Exp.fitMean(2.1));
        node4.setService(jobclass3, Exp.fitMean(1.9));
        node5.setService(jobclass1, Exp.fitMean(1.5));
        node5.setService(jobclass2, Exp.fitMean(0.9));
        node5.setService(jobclass3, Exp.fitMean(2.3));

        RoutingMatrix p = new RoutingMatrix(model, model.getJobClasses(), model.getNodes());
        p.addConnection(jobclass1, jobclass1, node1, node2, 1);
        p.addConnection(jobclass1, jobclass1, node2, node3, 0.25);
        p.addConnection(jobclass1, jobclass1, node2, node4, 0.25);
        p.addConnection(jobclass1, jobclass1, node2, node5, 0.25);
        p.addConnection(jobclass1, jobclass1, node2, node6, 0.25);
        p.addConnection(jobclass1, jobclass1, node3, node2, 1);
        p.addConnection(jobclass1, jobclass1, node4, node2, 1);
        p.addConnection(jobclass1, jobclass1, node5, node2, 1);
        p.addConnection(jobclass2, jobclass2, node1, node2, 1);
        p.addConnection(jobclass2, jobclass2, node2, node3, 0.25);
        p.addConnection(jobclass2, jobclass2, node2, node4, 0.25);
        p.addConnection(jobclass2, jobclass2, node2, node5, 0.25);
        p.addConnection(jobclass2, jobclass2, node2, node6, 0.25);
        p.addConnection(jobclass2, jobclass2, node3, node2, 1);
        p.addConnection(jobclass2, jobclass2, node4, node2, 1);
        p.addConnection(jobclass2, jobclass2, node5, node2, 1);
        p.addConnection(jobclass3, jobclass3, node1, node2, 1);
        p.addConnection(jobclass3, jobclass3, node2, node3, 0.25);
        p.addConnection(jobclass3, jobclass3, node2, node4, 0.25);
        p.addConnection(jobclass3, jobclass3, node2, node5, 0.25);
        p.addConnection(jobclass3, jobclass3, node2, node6, 0.25);
        p.addConnection(jobclass3, jobclass3, node3, node2, 1);
        p.addConnection(jobclass3, jobclass3, node4, node2, 1);
        p.addConnection(jobclass3, jobclass3, node5, node2, 1);

        model.link(p);

        return model;
    }

    public static Network test_cqn_twoclass_hyperl() throws IllegalAccessException {
        Network model = new Network("cqn_twoclass_hyperl");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Router node3 = new Router(model, "CS_Delay_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P
        Router node4 = new Router(model, "CS_Queue1_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, APH.fitMeanAndSCV(0.666667, 0.500)); // (Delay,Class1)
        node1.setService(jobclass2, APH.fitMeanAndSCV(0.216667, 1.579882)); // (Delay,Class2)
        node2.setService(jobclass1, APH.fitMeanAndSCV(0.190000, 5.038781)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.00)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4));

        routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 0.100); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node1, node3, 0.900); // (Delay,Class1) -> (CS_Delay_to_Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node2, node4, 1.00); // (Queue1,Class1) -> (CS_Queue1_to_Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node3, node1, 0.333333); // (CS_Delay_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node4, node1, 0.200); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass2, node3, node1, 0.666667); // (CS_Delay_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.addConnection(jobclass1, jobclass2, node4, node1, 0.800); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.addConnection(jobclass2, jobclass1, node4, node1, 1.00); // (CS_Queue1_to_Delay,Class2) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass2, jobclass2, node1, node2, 1.00); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.addConnection(jobclass2, jobclass2, node2, node4, 1.00); // (Queue1,Class2) -> (CS_Queue1_to_Delay,Class2)
        routingMatrix.addConnection(jobclass2, jobclass2, node3, node1, 1.00); // (CS_Delay_to_Delay,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network test_cqn_threeclass_hyperl() throws IllegalAccessException {
        Network model = new Network("cqn_threeclass_hyperl");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        node2.setNumberOfServers(2);
        Router node3 = new Router(model, "CS_Delay_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P
        Router node4 = new Router(model, "CS_Queue1_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 1, node1, 0);

        node1.setService(jobclass1, APH.fitMeanAndSCV(0.666667, 0.500)); // (Delay,Class1)
        node1.setService(jobclass2, APH.fitMeanAndSCV(0.216667, 1.579882)); // (Delay,Class2)
        node1.setService(jobclass3, Exp.fitMean(1.00)); // (Delay,Class3)
        node2.setService(jobclass1, APH.fitMeanAndSCV(0.190000, 5.038781)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.500)); // (Queue1,Class2)
        node2.setService(jobclass3, Exp.fitMean(0.333333)); // (Queue1,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2, jobclass3),
                Arrays.asList(node1, node2, node3, node4));

        routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 0.100); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node1, node3, 0.900); // (Delay,Class1) -> (CS_Delay_to_Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node2, node4, 1.00); // (Queue1,Class1) -> (CS_Queue1_to_Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node3, node1, 0.333333); // (CS_Delay_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node4, node1, 0.200); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass2, node3, node1, 0.666667); // (CS_Delay_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.addConnection(jobclass1, jobclass2, node4, node1, 0.800); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.addConnection(jobclass2, jobclass1, node4, node1, 1.00); // (CS_Queue1_to_Delay,Class2) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass2, jobclass2, node1, node2, 1.00); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.addConnection(jobclass2, jobclass2, node2, node4, 1.00); // (Queue1,Class2) -> (CS_Queue1_to_Delay,Class2)
        routingMatrix.addConnection(jobclass2, jobclass2, node3, node1, 1.00); // (CS_Delay_to_Delay,Class2) -> (Delay,Class2)
        routingMatrix.addConnection(jobclass3, jobclass3, node1, node2, 1.00); // (Delay,Class3) -> (Queue1,Class3)
        routingMatrix.addConnection(jobclass3, jobclass3, node2, node4, 1.00); // (Queue1,Class3) -> (CS_Queue1_to_Delay,Class3)
        routingMatrix.addConnection(jobclass3, jobclass3, node3, node1, 1.00); // (CS_Delay_to_Delay,Class3) -> (Delay,Class3)
        routingMatrix.addConnection(jobclass3, jobclass3, node4, node1, 1.00); // (CS_Queue1_to_Delay,Class3) -> (Delay,Class3)

        model.link(routingMatrix);

        return model;
    }

    public static Network test_example_closedModel_2() throws IllegalAccessException {
        Network model = new Network("example_closedModel_2");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Router node3 = new Router(model, "CS_Delay_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P
        Router node4 = new Router(model, "CS_Queue1_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, APH.fitMeanAndSCV(0.666667, 0.500)); // (Delay,Class1)
        node1.setService(jobclass2, APH.fitMeanAndSCV(0.216667, 1.579882)); // (Delay,Class2)
        node2.setService(jobclass1, APH.fitMeanAndSCV(0.190000, 5.038781)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.00)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4));

        routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 0.100); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node1, node3, 0.900); // (Delay,Class1) -> (CS_Delay_to_Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node2, node4, 1.00); // (Queue1,Class1) -> (CS_Queue1_to_Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node3, node1, 0.333333); // (CS_Delay_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node4, node1, 0.200); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass2, node3, node1, 0.666667); // (CS_Delay_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.addConnection(jobclass1, jobclass2, node4, node1, 0.800); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.addConnection(jobclass2, jobclass1, node4, node1, 1.00); // (CS_Queue1_to_Delay,Class2) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass2, jobclass2, node1, node2, 1.00); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.addConnection(jobclass2, jobclass2, node2, node4, 1.00); // (Queue1,Class2) -> (CS_Queue1_to_Delay,Class2)
        routingMatrix.addConnection(jobclass2, jobclass2, node3, node1, 1.00); // (CS_Delay_to_Delay,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network test_example_closedModel_3() throws IllegalAccessException {
        Network model = new Network("example_closedModel_3");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        node2.setNumberOfServers(2);
        Router node3 = new Router(model, "CS_Delay_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P
        Router node4 = new Router(model, "CS_Queue1_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 0, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 1, node1, 0);

        node1.setService(jobclass1, APH.fitMeanAndSCV(0.666667, 0.500)); // (Delay,Class1)
        node1.setService(jobclass2, APH.fitMeanAndSCV(0.216667, 1.579882)); // (Delay,Class2)
        node1.setService(jobclass3, Exp.fitMean(1.00)); // (Delay,Class3)
        node2.setService(jobclass1, APH.fitMeanAndSCV(0.190000, 5.038781)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.500)); // (Queue1,Class2)
        node2.setService(jobclass3, Exp.fitMean(0.333333)); // (Queue1,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2, jobclass3),
                Arrays.asList(node1, node2, node3, node4));

        routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 0.100); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node1, node3, 0.900); // (Delay,Class1) -> (CS_Delay_to_Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node2, node4, 1.00); // (Queue1,Class1) -> (CS_Queue1_to_Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node3, node1, 0.333333); // (CS_Delay_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass1, node4, node1, 0.200); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass1, jobclass2, node3, node1, 0.666667); // (CS_Delay_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.addConnection(jobclass1, jobclass2, node4, node1, 0.800); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.addConnection(jobclass2, jobclass1, node4, node1, 1.00); // (CS_Queue1_to_Delay,Class2) -> (Delay,Class1)
        routingMatrix.addConnection(jobclass2, jobclass2, node1, node2, 1.00); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.addConnection(jobclass2, jobclass2, node2, node4, 1.00); // (Queue1,Class2) -> (CS_Queue1_to_Delay,Class2)
        routingMatrix.addConnection(jobclass2, jobclass2, node3, node1, 1.00); // (CS_Delay_to_Delay,Class2) -> (Delay,Class2)
        routingMatrix.addConnection(jobclass3, jobclass3, node1, node2, 1.00); // (Delay,Class3) -> (Queue1,Class3)
        routingMatrix.addConnection(jobclass3, jobclass3, node2, node4, 1.00); // (Queue1,Class3) -> (CS_Queue1_to_Delay,Class3)
        routingMatrix.addConnection(jobclass3, jobclass3, node3, node1, 1.00); // (CS_Delay_to_Delay,Class3) -> (Delay,Class3)
        routingMatrix.addConnection(jobclass3, jobclass3, node4, node1, 1.00); // (CS_Queue1_to_Delay,Class3) -> (Delay,Class3)

        model.link(routingMatrix);

        return model;
    }

}
