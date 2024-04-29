package jline.examples;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.*;
import jline.lang.nodes.*;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.solvers.NetworkAvgTable;

import java.util.Arrays;

/**
 * Examples of mixed queueing networks
 */
public class MixedModel {

    public static Network ex1() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Source node5 = new Source(model, "Source");
        Sink node6 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 100, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue1,ClosedClass)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue1,OpenClass)
        node2.setService(jobclass1, Exp.fitMean(0.500000)); // (Queue2,ClosedClass)
        node2.setService(jobclass2, Exp.fitMean(0.707107)); // (Queue2,OpenClass)
        node3.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue3,ClosedClass)
        node3.setService(jobclass2, Exp.fitMean(0.577350)); // (Queue3,OpenClass)
        node4.setService(jobclass1, Exp.fitMean(0.250000)); // (Queue4,ClosedClass)
        node4.setService(jobclass2, Exp.fitMean(0.500000)); // (Queue4,OpenClass)
        node5.setArrival(jobclass1, Disabled.getInstance()); // (Source,ClosedClass)
        node5.setArrival(jobclass2, APH.fitMeanAndSCV(3.000000,64.000000)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4, node5, node6));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Queue1,ClosedClass) -> (Queue2,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue2,ClosedClass) -> (Queue3,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Queue3,ClosedClass) -> (Queue4,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.000000); // (Queue4,ClosedClass) -> (Queue1,ClosedClass)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Queue1,OpenClass) -> (Queue2,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue2,OpenClass) -> (Queue3,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node6, 1.000000); // (Queue3,OpenClass) -> (Sink,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node4, node1, 1.000000); // (Queue4,OpenClass) -> (Queue1,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node5, node1, 1.000000); // (Source,OpenClass) -> (Queue1,OpenClass)

        model.link(routingMatrix);

        return model;
    }

    public static Network ex1_line() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Source node3 = new Source(model, "Source");
        Sink node4 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 2, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, new Erlang(3, 2)); // (Delay,ClosedClass)
        node1.setService(jobclass2, new HyperExp(0.5,3.0,10.0)); // (Delay,OpenClass)
        node2.setService(jobclass1, new HyperExp(0.1,1.0,10.0)); // (Queue1,ClosedClass)
        node2.setService(jobclass2, new Exp(1)); // (Queue1,OpenClass)
        node3.setArrival(jobclass2, new Exp(0.1)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay,ClosedClass) -> (Queue1,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,ClosedClass) -> (Delay,ClosedClass)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,OpenClass) -> (Queue1,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node4, 1.000000); // (Queue1,OpenClass) -> (Sink,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.000000); // (Source,OpenClass) -> (Delay,OpenClass)

        model.link(routingMatrix);

        return model;
    }

    public static Network ex2() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.PS);
        node2.setNumberOfServers(2);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.PS);
        node3.setNumberOfServers(3);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.PS);
        node4.setNumberOfServers(4);
        Queue node5 = new Queue(model, "Queue5", SchedStrategy.PS);
        node5.setNumberOfServers(5);
        Source node6 = new Source(model, "Source");
        Sink node7 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 3, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue1,ClosedClass)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue1,OpenClass)
        node2.setService(jobclass1, Exp.fitMean(0.500000)); // (Queue2,ClosedClass)
        node2.setService(jobclass2, Exp.fitMean(0.707107)); // (Queue2,OpenClass)
        node3.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue3,ClosedClass)
        node3.setService(jobclass2, Exp.fitMean(0.577350)); // (Queue3,OpenClass)
        node4.setService(jobclass1, Exp.fitMean(0.250000)); // (Queue4,ClosedClass)
        node4.setService(jobclass2, Exp.fitMean(0.500000)); // (Queue4,OpenClass)
        node5.setService(jobclass1, Exp.fitMean(0.200000)); // (Queue5,ClosedClass)
        node5.setService(jobclass2, Exp.fitMean(0.447214)); // (Queue5,OpenClass)
        node6.setArrival(jobclass1, Disabled.getInstance()); // (Source,ClosedClass)
        node6.setArrival(jobclass2, Exp.fitMean(3.333333)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4, node5, node6, node7));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Queue1,ClosedClass) -> (Queue2,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue2,ClosedClass) -> (Queue3,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Queue3,ClosedClass) -> (Queue4,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.000000); // (Queue4,ClosedClass) -> (Queue1,ClosedClass)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Queue1,OpenClass) -> (Queue2,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue2,OpenClass) -> (Queue3,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node5, 1.000000); // (Queue3,OpenClass) -> (Queue5,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node4, node1, 1.000000); // (Queue4,OpenClass) -> (Queue1,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node5, node7, 1.000000); // (Queue5,OpenClass) -> (Sink,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node6, node1, 1.000000); // (Source,OpenClass) -> (Queue1,OpenClass)

        model.link(routingMatrix);

        return model;
    }

    public static Network ex3() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        node2.setNumberOfServers(2);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        node3.setNumberOfServers(3);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        node4.setNumberOfServers(4);
        Queue node5 = new Queue(model, "Queue5", SchedStrategy.FCFS);
        node5.setNumberOfServers(5);
        Source node6 = new Source(model, "Source");
        Sink node7 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 3, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue1,ClosedClass)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue1,OpenClass)
        node2.setService(jobclass1, Exp.fitMean(0.500000)); // (Queue2,ClosedClass)
        node2.setService(jobclass2, Exp.fitMean(0.707107)); // (Queue2,OpenClass)
        node3.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue3,ClosedClass)
        node3.setService(jobclass2, Exp.fitMean(0.577350)); // (Queue3,OpenClass)
        node4.setService(jobclass1, Exp.fitMean(0.250000)); // (Queue4,ClosedClass)
        node4.setService(jobclass2, Exp.fitMean(0.500000)); // (Queue4,OpenClass)
        node5.setService(jobclass1, Exp.fitMean(0.200000)); // (Queue5,ClosedClass)
        node5.setService(jobclass2, Exp.fitMean(0.447214)); // (Queue5,OpenClass)
        node6.setArrival(jobclass1, Disabled.getInstance()); // (Source,ClosedClass)
        node6.setArrival(jobclass2, Exp.fitMean(3.333333)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4, node5, node6, node7));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Queue1,ClosedClass) -> (Queue2,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue2,ClosedClass) -> (Queue3,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Queue3,ClosedClass) -> (Queue4,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.000000); // (Queue4,ClosedClass) -> (Queue1,ClosedClass)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Queue1,OpenClass) -> (Queue2,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue2,OpenClass) -> (Queue3,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node5, 1.000000); // (Queue3,OpenClass) -> (Queue5,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node4, node1, 1.000000); // (Queue4,OpenClass) -> (Queue1,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node5, node7, 1.000000); // (Queue5,OpenClass) -> (Sink,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node6, node1, 1.000000); // (Source,OpenClass) -> (Queue1,OpenClass)

        model.link(routingMatrix);

        return model;
    }

    public static Network ex4() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Source node5 = new Source(model, "Source");
        Sink node6 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 100, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue1,ClosedClass)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue1,OpenClass)
        node2.setService(jobclass1, Exp.fitMean(0.500000)); // (Queue2,ClosedClass)
        node2.setService(jobclass2, Exp.fitMean(0.707107)); // (Queue2,OpenClass)
        node3.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue3,ClosedClass)
        node3.setService(jobclass2, Exp.fitMean(0.577350)); // (Queue3,OpenClass)
        node4.setService(jobclass1, Exp.fitMean(0.250000)); // (Queue4,ClosedClass)
        node4.setService(jobclass2, Exp.fitMean(0.500000)); // (Queue4,OpenClass)
        node5.setArrival(jobclass1, Disabled.getInstance()); // (Source,ClosedClass)
        node5.setArrival(jobclass2, APH.fitMeanAndSCV(3.000000, 64.000000)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4, node5, node6));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Queue1,ClosedClass) -> (Queue2,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue2,ClosedClass) -> (Queue3,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Queue3,ClosedClass) -> (Queue4,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.000000); // (Queue4,ClosedClass) -> (Queue1,ClosedClass)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Queue1,OpenClass) -> (Queue2,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue2,OpenClass) -> (Queue3,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node6, 1.000000); // (Queue3,OpenClass) -> (Sink,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node4, node1, 1.000000); // (Queue4,OpenClass) -> (Queue1,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node5, node1, 1.000000); // (Source,OpenClass) -> (Queue1,OpenClass)


        model.link(routingMatrix);

        return model;
    }


    public static Network ex5() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.PS);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.PS);
        Source node5 = new Source(model, "Source");
        Sink node6 = new Sink(model, "Sink");

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 100, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "OpenClass", 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue1,ClosedClass)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue1,OpenClass)
        node2.setService(jobclass1, Exp.fitMean(0.500000)); // (Queue2,ClosedClass)
        node2.setService(jobclass2, Exp.fitMean(0.707107)); // (Queue2,OpenClass)
        node3.setService(jobclass1, Exp.fitMean(0.333333)); // (Queue3,ClosedClass)
        node3.setService(jobclass2, Exp.fitMean(0.577350)); // (Queue3,OpenClass)
        node4.setService(jobclass1, Exp.fitMean(0.250000)); // (Queue4,ClosedClass)
        node4.setService(jobclass2, Exp.fitMean(0.500000)); // (Queue4,OpenClass)


        node5.setArrival(jobclass1, Disabled.getInstance()); // (Source,ClosedClass)
        node5.setArrival(jobclass2, new Exp(0.3)); // (Source,OpenClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4, node5, node6));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Queue1,ClosedClass) -> (Queue2,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue2,ClosedClass) -> (Queue3,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Queue3,ClosedClass) -> (Queue4,ClosedClass)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.000000); // (Queue4,ClosedClass) -> (Queue1,ClosedClass)


        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Queue1,OpenClass) -> (Queue2,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue2,OpenClass) -> (Queue3,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node3, node6, 1.000000); // (Queue3,OpenClass) -> (Sink,OpenClass)

//        routingMatrix.set(jobclass2, jobclass2, node4, node1, 1.000000); // (Queue4,OpenClass) -> (Queue1,OpenClass)
        routingMatrix.set(jobclass2, jobclass2, node5, node1, 1.000000); // (Source,OpenClass) -> (Queue1,OpenClass)


        model.link(routingMatrix);

        return model;

    }



    public static void main(String[] args) throws Exception {
        Network model = ex2();

        NetworkStruct sn = model.getStruct(false);
        sn.rt.print();
        SolverOptions options = new SolverOptions(SolverType.MVA);
        options.method = "exact";
        NetworkSolver solver = new SolverMVA(model, options);
        NetworkAvgTable t = solver.getAvgTable();
        t.print(options);
    }
}
