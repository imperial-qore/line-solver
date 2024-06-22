package jline.examples;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.lang.processes.Replayer;
import jline.util.Matrix;

/**
 * Examples of open queueing networks
 */
public class OpenModel {
    public static Network example_openModel_1(){
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Source node3 = new Source(model, "Source");
        Sink node4 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);

        node1.setService(jobclass1, new HyperExp(0.5, 3.0, 10.0));
        node2.setService(jobclass1, new Exp(1));
        //node2.setService(jobclass1, new Replayer("/home/gcasale/Dropbox/code/line-solver.git/python/gettingstarted/example_trace.txt"));
        node3.setArrival(jobclass1, new Exp(0.1));

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1);
        routingMatrix.set(jobclass1, jobclass1, node2, node4, 1);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1);

        model.link(routingMatrix);

        return model;
    }


    public static Network example_openModel_2() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Delay node2 = new Delay(model, "Station1");
        Queue node3 = new Queue(model, "Station2", SchedStrategy.PS);
        Queue node4 = new Queue(model, "Station3", SchedStrategy.PS);
        Sink node5 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        // TODO: model building quite different from original though equivalent model
        node1.setArrival(jobclass1, Exp.fitMean(50.000000)); // (Source,Class1)
        node1.setArrival(jobclass2, Exp.fitMean(25.000000)); // (Source,Class2)
        node2.setService(jobclass1, Exp.fitMean(91.000000)); // (Station1,Class1)
        node2.setService(jobclass2, Exp.fitMean(92.000000)); // (Station1,Class2)
        node3.setService(jobclass1, Exp.fitMean(10.000000)); // (Station2,Class1)
        node3.setService(jobclass2, Exp.fitMean(5.000000)); // (Station2,Class2)
        node4.setService(jobclass1, Exp.fitMean(5.000000)); // (Station3,Class1)
        node4.setService(jobclass2, Exp.fitMean(9.000000)); // (Station3,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Source,Class1) -> (Station1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Station1,Class1) -> (Station2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Station2,Class1) -> (Station3,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node5, 1.000000); // (Station3,Class1) -> (Sink,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Source,Class2) -> (Station1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Station1,Class2) -> (Station2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node4, 1.000000); // (Station2,Class2) -> (Station3,Class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node5, 1.000000); // (Station3,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_openModel_3(){
        Network model = new Network("model");

        // Block 1: nodes
        Source node1 = new Source(model, "Source 1");
        Queue node2 = new Queue(model, "Queue 1", SchedStrategy.PS);
        ClassSwitch node3 = new ClassSwitch(model, "ClassSwitch 1");
        Sink node4 = new Sink(model, "Sink 1");
        Queue node5 = new Queue(model, "Queue 2", SchedStrategy.PS);

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class A", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class B", 0);
        OpenClass jobclass3 = new OpenClass(model, "Class C", 0);

        node1.setArrival(jobclass1, Exp.fitMean(0.500000)); // (Source 1,Class A)
        node1.setArrival(jobclass2, Exp.fitMean(1.000000)); // (Source 1,Class B)
        node1.setArrival(jobclass3, Disabled.getInstance()); // (Source 1,Class C)
        node2.setService(jobclass1, Exp.fitMean(0.200000)); // (Queue 1,Class A)
        node2.setService(jobclass2, Exp.fitMean(0.300000)); // (Queue 1,Class B)
        node2.setService(jobclass3, Exp.fitMean(0.333333)); // (Queue 1,Class C)
        node5.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue 2,Class A)
        node5.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue 2,Class B)
        node5.setService(jobclass3, Exp.fitMean(0.150000)); // (Queue 2,Class C)

        ClassSwitchMatrix C = node3.initClassSwitchMatrix();
        C = new ClassSwitchMatrix(Matrix.eye(C.getNumRows()));
        node3.setClassSwitchingMatrix(C);

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Source 1,Class A) -> (Queue 1,Class A)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue 1,Class A) -> (ClassSwitch 1,Class A)
        routingMatrix.set(jobclass1, jobclass1, node5, node4, 1.000000); // (Queue 2,Class A) -> (Sink 1,Class A)
        routingMatrix.set(jobclass1, jobclass3, node3, node5, 1.000000); // (CS_ClassSwitch 1_to_Queue 2,Class A) -> (Queue 2,Class C)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Source 1,Class B) -> (Queue 1,Class B)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue 1,Class B) -> (ClassSwitch 1,Class B)
        routingMatrix.set(jobclass2, jobclass2, node5, node4, 1.000000); // (Queue 2,Class B) -> (Sink 1,Class B)
        routingMatrix.set(jobclass2, jobclass3, node3, node5, 1.000000); // (CS_ClassSwitch 1_to_Queue 2,Class B) -> (Queue 2,Class C)
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.000000); // (Source 1,Class C) -> (Queue 1,Class C)
        routingMatrix.set(jobclass3, jobclass3, node2, node3, 1.000000); // (Queue 1,Class C) -> (ClassSwitch 1,Class C)
        routingMatrix.set(jobclass3, jobclass3, node3, node5, 1.000000); // (ClassSwitch 1,Class C) -> (CS_ClassSwitch 1_to_Queue 2,Class C)
        routingMatrix.set(jobclass3, jobclass3, node5, node4, 1.000000); // (Queue 2,Class C) -> (Sink 1,Class C)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_openModel_4() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "OpenClass", 0);

        node1.setArrival(jobclass1, Exp.fitMean(1.000000)); // (Source,OpenClass)
        node2.setService(jobclass1, new Replayer("example_trace.txt")); // (Queue,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Source,OpenClass) -> (Queue,OpenClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue,OpenClass) -> (Sink,OpenClass)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_openModel_5() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");
        Router node4 = new Router(model, "VSink1");
        Router node5 = new Router(model, "VSink2");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);

        node1.setArrival(jobclass1, Exp.fitMean(1.000000)); // (Source,Class1)
        node2.setService(jobclass1, Exp.fitMean(0.010000)); // (Queue1,Class1)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Source,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node4, 0.600000); // (Queue1,Class1) -> (VSink1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node5, 0.400000); // (Queue1,Class1) -> (VSink2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node3, 1.000000); // (VSink1,Class1) -> (Sink,Class1)
        routingMatrix.set(jobclass1, jobclass1, node5, node3, 1.000000); // (VSink2,Class1) -> (Sink,Class1)

        model.link(routingMatrix);

        return model;
    }

}
