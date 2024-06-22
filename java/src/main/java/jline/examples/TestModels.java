package jline.examples;

import jline.lang.*;
import jline.lang.constant.GlobalConstants;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.constant.VerboseLevel;
import jline.lang.distributions.*;
import jline.lang.layered.*;
import jline.lang.nodes.*;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.env.SolverEnv;
import jline.solvers.fluid.SolverFluid;
import jline.util.Matrix;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Examples of caching models
 */
public class TestModels {


    public static void main(String[] args) {

    }

    public static Network test_closedModel_2() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Router node3 = new Router(model, "CS_Delay_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P
        Router node4 = new Router(model, "CS_Queue1_to_Delay"); // Dummy node, class switching is embedded in the routing matrix P

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);

        node1.setService(jobclass1, APH.fitMeanAndSCV(0.666667,0.500000)); // (Delay,Class1)
        node1.setService(jobclass2, APH.fitMeanAndSCV(0.216667,1.579882)); // (Delay,Class2)
        node2.setService(jobclass1, APH.fitMeanAndSCV(0.190000,5.038781)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue1,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.100000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node3, 0.900000); // (Delay,Class1) -> (CS_Delay_to_Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node4, 1.000000); // (Queue1,Class1) -> (CS_Queue1_to_Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 0.333333); // (CS_Delay_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 0.200000); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass2, node3, node1, 0.666667); // (CS_Delay_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.set(jobclass1, jobclass2, node4, node1, 0.800000); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass1, node4, node1, 1.000000); // (CS_Queue1_to_Delay,Class2) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node4, 1.000000); // (Queue1,Class2) -> (CS_Queue1_to_Delay,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.000000); // (CS_Delay_to_Delay,Class2) -> (Delay,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network test_closedModel_3() {
        Network model = new Network("myModel");

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

        node1.setService(jobclass1, APH.fitMeanAndSCV(0.666667,0.500000)); // (Delay,Class1)
        node1.setService(jobclass2, APH.fitMeanAndSCV(0.216667,1.579882)); // (Delay,Class2)
        node1.setService(jobclass3, Exp.fitMean(1.000000)); // (Delay,Class3)
        node2.setService(jobclass1, APH.fitMeanAndSCV(0.190000,5.038781)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.500000)); // (Queue1,Class2)
        node2.setService(jobclass3, Exp.fitMean(0.333333)); // (Queue1,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.100000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node3, 0.900000); // (Delay,Class1) -> (CS_Delay_to_Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node4, 1.000000); // (Queue1,Class1) -> (CS_Queue1_to_Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 0.333333); // (CS_Delay_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 0.200000); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass2, node3, node1, 0.666667); // (CS_Delay_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.set(jobclass1, jobclass2, node4, node1, 0.800000); // (CS_Queue1_to_Delay,Class1) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass1, node4, node1, 1.000000); // (CS_Queue1_to_Delay,Class2) -> (Delay,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node4, 1.000000); // (Queue1,Class2) -> (CS_Queue1_to_Delay,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.000000); // (CS_Delay_to_Delay,Class2) -> (Delay,Class2)
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.000000); // (Delay,Class3) -> (Queue1,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node4, 1.000000); // (Queue1,Class3) -> (CS_Queue1_to_Delay,Class3)
        routingMatrix.set(jobclass3, jobclass3, node3, node1, 1.000000); // (CS_Delay_to_Delay,Class3) -> (Delay,Class3)
        routingMatrix.set(jobclass3, jobclass3, node4, node1, 1.000000); // (CS_Queue1_to_Delay,Class3) -> (Delay,Class3)

        model.link(routingMatrix);

        return model;
    }

    public static Network test_closedModel_4() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        node2.setNumberOfServers(3);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        node3.setNumberOfServers(3);
        Router node4 = new Router(model, "CS_Delay_to_Queue1"); // Dummy node, class switching is embedded in the routing matrix P

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);
        ClosedClass jobclass3 = new ClosedClass(model, "Class3", 2, node1, 0);
        ClosedClass jobclass4 = new ClosedClass(model, "Class4", 1, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node1.setService(jobclass2, Exp.fitMean(1.000000)); // (Delay,Class2)
        node1.setService(jobclass3, Exp.fitMean(0.100000)); // (Delay,Class3)
        node1.setService(jobclass4, Exp.fitMean(1.000000)); // (Delay,Class4)
        node2.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue1,Class1)
        node2.setService(jobclass2, APH.fitMeanAndSCV(2.000000,0.500000)); // (Queue1,Class2)
        node2.setService(jobclass3, Exp.fitMean(0.100000)); // (Queue1,Class3)
        node2.setService(jobclass4, Exp.fitMean(1.000000)); // (Queue1,Class4)
        node3.setService(jobclass1, Disabled.getInstance()); // (Queue2,Class1)
        node3.setService(jobclass2, Disabled.getInstance()); // (Queue2,Class2)
        node3.setService(jobclass3, APH.fitMeanAndSCV(2.000000,0.500000)); // (Queue2,Class3)
        node3.setService(jobclass4, Exp.fitMean(1.000000)); // (Queue2,Class4)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node4, 1.000000); // (Delay,Class1) -> (CS_Delay_to_Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000); // (Queue2,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node2, 0.500000); // (CS_Delay_to_Queue1,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass2, node4, node2, 0.500000); // (CS_Delay_to_Queue1,Class1) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass1, node4, node2, 1.000000); // (CS_Delay_to_Queue1,Class2) -> (Queue1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node4, 1.000000); // (Delay,Class2) -> (CS_Delay_to_Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.000000); // (Queue1,Class2) -> (Delay,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.000000); // (Queue2,Class2) -> (Delay,Class2)
        routingMatrix.set(jobclass3, jobclass3, node1, node3, 0.250000); // (Delay,Class3) -> (Queue2,Class3)
        routingMatrix.set(jobclass3, jobclass3, node1, node4, 0.750000); // (Delay,Class3) -> (CS_Delay_to_Queue1,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node1, 1.000000); // (Queue1,Class3) -> (Delay,Class3)
        routingMatrix.set(jobclass3, jobclass3, node3, node1, 1.000000); // (Queue2,Class3) -> (Delay,Class3)
        routingMatrix.set(jobclass3, jobclass3, node4, node2, 0.333333); // (CS_Delay_to_Queue1,Class3) -> (Queue1,Class3)
        routingMatrix.set(jobclass3, jobclass4, node4, node2, 0.666667); // (CS_Delay_to_Queue1,Class3) -> (Queue1,Class4)
        routingMatrix.set(jobclass4, jobclass3, node4, node2, 1.000000); // (CS_Delay_to_Queue1,Class4) -> (Queue1,Class3)
        routingMatrix.set(jobclass4, jobclass4, node1, node4, 1.000000); // (Delay,Class4) -> (CS_Delay_to_Queue1,Class4)
        routingMatrix.set(jobclass4, jobclass4, node2, node1, 1.000000); // (Queue1,Class4) -> (Delay,Class4)
        routingMatrix.set(jobclass4, jobclass4, node3, node1, 1.000000); // (Queue2,Class4) -> (Delay,Class4)

        model.link(routingMatrix);

        return model;
    }

    public static Network test_closedModel_5() {
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
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

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

    public static Network test_closedModel_8() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node2.setService(jobclass1, Exp.fitMean(1.500000)); // (Queue1,Class1)
        node3.setService(jobclass1, Exp.fitMean(2.000000)); // (Queue2, Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.500000); // (Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.300000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node3, 0.200000); // (Delay,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000); // (Queue2,Class1) -> (Delay,Class1)
        model.link(routingMatrix);

        // TODO
        return model;
    }

    public static Network test_closedModel_9() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 3, node1, 0);

        node1.setService(jobclass1, new Erlang(1, 5)); // (Delay,Class1)
        node2.setService(jobclass1, new Erlang(2, 5)); // (Queue1,Class1)
        node3.setService(jobclass1, new Erlang(0.5, 2)); // (Queue2, Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.500000); // (Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.300000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node3, 0.200000); // (Delay,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000); // (Queue2,Class1) -> (Delay,Class1)
        model.link(routingMatrix);

        return model;
    }

    public static Network test_CQN_JMT_2() {
        Network model = new Network("jmtregression_11");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay 1");
        Queue node2 = new Queue(model, "Queue 1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue 2", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 1, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(0.010000)); // (Delay 1,Class1)
        node2.setService(jobclass1, Exp.fitMean(75.000000)); // (Queue 1,Class1)
        node3.setService(jobclass1, Exp.fitMean(5.000000)); // (Queue 2,Class1)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay 1,Class1) -> (Queue 1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue 1,Class1) -> (Queue 2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 1.000000); // (Queue 2,Class1) -> (Delay 1,Class1)

        model.link(routingMatrix);

        return model;
    }

    public static Network ex1() {
        /*  M/M/1 queue
         */
        Network model = new Network("MM1LowU");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass openClass = new OpenClass(model, "Open Class");
        source.setArrival(openClass, new Exp(2));
        queue.setService(openClass, new Exp(10));
        model.link(model.serialRouting(source, queue, sink));

        return model;
    }

    public static Network ex2() {
        /*  M/M/2 queue
         */
        Network model = new Network("MM2HighU");
        OpenClass openClass = new OpenClass(model, "Open Class");
        Source source = new Source(model, "Source");
        source.setArrival(openClass, new Exp(8));
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setService(openClass, new Exp(10));
        queue.setNumberOfServers(1);
        Sink sink = new Sink(model, "Sink");

        model.link(model.serialRouting(source, queue, sink));

        return model;
    }

    public static Network ex3() {
        /*  3 markovian queues in series
         */
        Network model = new Network("3 Series");
        OpenClass openClass = new OpenClass(model, "Open Class");
        Source source = new Source(model, "Source");
        source.setArrival(openClass, new Exp(8));
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        queue1.setService(openClass, new Exp(12));
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        queue2.setService(openClass, new Exp(11));
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        queue3.setService(openClass, new Exp(10));
        Sink sink = new Sink(model, "Sink");

        model.link(model.serialRouting(source, queue1, queue2, queue3, sink));

        return model;
    }

    public static Network ex4() {
        /*  3 queues in parallel with Erlang-distributed service times
         */
        Network model = new Network("Parallel Erlang");
        OpenClass openClass = new OpenClass(model, "Open Class");
        Source source = new Source(model, "Source");
        source.setArrival(openClass, new Exp(30));
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        queue1.setService(openClass, new Erlang(12, 3));
        queue1.setNumberOfServers(3);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.LCFS);
        queue2.setService(openClass, new Erlang(12, 3));
        queue2.setNumberOfServers(3);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.SIRO);
        queue3.setService(openClass, new Erlang(12, 3));
        queue3.setNumberOfServers(3);
        Sink sink = new Sink(model, "Sink");


        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(source, queue1);
        routingMatrix.set(queue1, sink);
        routingMatrix.set(source, queue2);
        routingMatrix.set(queue2, sink);
        routingMatrix.set(source, queue3);
        routingMatrix.set(queue3, sink);

        model.link(routingMatrix);


        return model;
    }

    public static Network ex5() {
        /* A closed network of 3 queues
         */
        Network model = new Network("3 Closed");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);

        ClosedClass closedClass = new ClosedClass(model, "Closed Class", 3, queue1);

        queue1.setService(closedClass, new Exp(1));
        queue2.setService(closedClass, new Exp(2));
        queue3.setService(closedClass, new Exp(3));

        model.link(model.serialRouting(queue1, queue2, queue3));
        return model;
    }

    public static Network ex6() {
        /*An M/M/1/10 queue
         */
        Network model = new Network("MM1 10");
        OpenClass openClass = new OpenClass(model, "Open Class");
        Source source = new Source(model, "Source");
        source.setArrival(openClass, new Exp(8));
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setService(openClass, new Exp(10));
        queue.setCap(10);
        Sink sink = new Sink(model, "Sink");

        model.link(model.serialRouting(source, queue, sink));

        return model;
    }

    public static Network erlangExample1() {
        Network model = new Network("MRP");
        Delay delay = new Delay(model, "Working State");
        Queue queue = new Queue(model, "RepairQueue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);

        ClosedClass closedClass = new ClosedClass(model, "Machines", 3, delay);
        delay.setService(closedClass, new Exp(0.5));
        queue.setService(closedClass, new Erlang(1, 2));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    public static Network erlangExample2() {
        Network model = new Network("MRP");
        Delay delay = new Delay(model, "Working State");
        Queue queue = new Queue(model, "RepairQueue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);

        ClosedClass closedClass = new ClosedClass(model, "Machines", 3, delay);
        delay.setService(closedClass, new Exp(0.5));
        queue.setService(closedClass, new Erlang(5, 2));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    public static Network erlangExample3() {
        Network model = new Network("MRP");
        Delay delay = new Delay(model, "Working State");
        Queue queue = new Queue(model, "RepairQueue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);

        ClosedClass closedClass = new ClosedClass(model, "Machines", 4, delay);
        delay.setService(closedClass, new Erlang(1, 5));
        queue.setService(closedClass, new Erlang(2, 3));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    public static Network ctmcExample1() {
        Network model = new Network("M/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");

        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, new Exp(1));
        queue.setNumberOfServers(1);
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    public static Network ex7() {
        /* A queue with two different open classes
         */
        Network model = new Network("2CDSDC");
        OpenClass openClass1 = new OpenClass(model, "Open 1");
        OpenClass openClass2 = new OpenClass(model, "Open 2");
        Source source = new Source(model, "Source");
        source.setArrival(openClass1, new Exp(8));
        source.setArrival(openClass2, new Exp(5));
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setService(openClass1, new Exp(12));
        queue.setService(openClass2, new Exp(16));
        queue.setClassCap(openClass1, 5);
        queue.setClassCap(openClass2, 3);
        Sink sink = new Sink(model, "Sink");

        model.link(model.serialRouting(source, queue, sink));

        return model;
    }

    public static Network ex8() {
        Network model = new Network("2CDSDC");
        Delay Node1 = new Delay(model, "Delay");
        Queue Node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        ClosedClass closedClass1 = new ClosedClass(model, "Closed 1", 10, Node1, 0);

        Node1.setService(closedClass1, new Exp(1));
        Node2.setService(closedClass1, new Exp(0.6666667));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(Node1, Node1, 0.7);
        routingMatrix.set(closedClass1, Node1, Node2, 0.3);
        routingMatrix.set(closedClass1, Node2, Node1, 1.0);
        routingMatrix.set(closedClass1, Node2, Node2, 0.0);

        model.link(routingMatrix);
        return model;
    }

    public static Network ex9() {
        Network model = new Network("ForkJoin");

        Source source = new Source(model, "Src");
        Fork fork = new Fork(model);
        Queue upper = new Queue(model, "upper Q");
        Queue lower = new Queue(model, "Lower Q");
        Join join = new Join(model);
        Queue post = new Queue(model, "Post");
        Sink sink = new Sink(model, "Sink");

        OpenClass openClass = new OpenClass(model, "oclass");

        source.setArrival(openClass, new Exp(10));
        upper.setService(openClass, new Exp(6));
        lower.setService(openClass, new Exp(8));
        post.setService(openClass, new Exp(15));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(openClass, source, fork, 1);
        routingMatrix.set(openClass, fork, upper, 0.5);
        routingMatrix.set(openClass, fork, lower, 0.5);
        routingMatrix.set(openClass, upper, join, 1);
        routingMatrix.set(openClass, lower, join, 1);
        routingMatrix.set(openClass, join, post, 1);
        routingMatrix.set(openClass, post, sink, 1);

        model.link(routingMatrix);
        return model;
    }

    public static Network ex10() {
        Network model = new Network("closed net");
        Queue Node1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue Node2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        ClosedClass closedClass1 = new ClosedClass(model, "Closed 1", 6, Node1, 0);

        Node1.setService(closedClass1, new Exp(3));
        Node2.setService(closedClass1, new Exp(5));

        model.link(model.serialRouting(Node1, Node2));
        return model;
    }

    public static Network ex11() {
        /*  3 queues in parallel with Erlang-distributed service times
         */
        Network model = new Network("Parallel Erlang");
        OpenClass openClass = new OpenClass(model, "Open Class 1");
        OpenClass oClass2 = new OpenClass(model, "Open Class 2");
        Source source = new Source(model, "Source");
        source.setArrival(openClass, new Exp(30));
        source.setArrival(oClass2, new Exp(10));
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        queue1.setService(openClass, new Erlang(12, 3));
        queue1.setService(oClass2, new Erlang(15, 2));
        queue1.setNumberOfServers(3);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        queue2.setService(openClass, new Erlang(12, 3));
        queue2.setService(oClass2, new Erlang(15, 2));
        queue2.setNumberOfServers(3);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.SIRO);
        queue3.setService(openClass, new Erlang(12, 3));
        queue3.setService(oClass2, new Erlang(15, 2));
        queue3.setNumberOfServers(3);
        Sink sink = new Sink(model, "Sink");


        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(source, queue1);
        routingMatrix.set(queue1, sink);
        routingMatrix.set(source, queue2);
        routingMatrix.set(queue2, sink);
        routingMatrix.set(source, queue3);
        routingMatrix.set(queue3, sink);

        model.link(routingMatrix);

        return model;
    }

    public static Network aphExample1() {
        Network model = new Network("MRP");
        Delay delay = new Delay(model, "Working State");
        Queue queue = new Queue(model, "RepairQueue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);

        ClosedClass closedClass = new ClosedClass(model, "Machines", 3, delay);
        delay.setService(closedClass, new Exp(0.5));
        queue.setService(closedClass, APH.fitMeanAndSCV(4.0, 0.5));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    public static Network aphExample2() {
        Network model = new Network("M/APH/1");
        Delay delay = new Delay(model, "Working State");
        Queue queue = new Queue(model, "RepairQueue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);
        ClosedClass closedClass1 = new ClosedClass(model, "Machines", 1, delay);
        ClosedClass closedClass2 = new ClosedClass(model, "Machines2", 1, delay);

        delay.setService(closedClass1, APH.fitMeanAndSCV(16.0000, 16.0000));
        delay.setService(closedClass2, APH.fitMeanAndSCV(8.0000, 16.0000));
        queue.setService(closedClass1, new Exp(0.1));
        queue.setService(closedClass2, new Exp(0.1));

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(closedClass1, closedClass1, delay, queue, 1.000000); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(closedClass1, closedClass1, queue, delay, 1.000000); // (Queue1,Class1) -> (Delay,Class1)
        routingMatrix.set(closedClass2, closedClass2, delay, queue, 1.000000); // (Delay,Class2) -> (Queue1,Class2)
        routingMatrix.set(closedClass2, closedClass2, queue, delay, 1.000000); // (Queue1,Class2) -> (Delay,Class2)

        model.link(routingMatrix);
//        model.link(model.serialRouting(delay, queue));

        return model;
    }

    public static Network aphExample3() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Queue node1 = new Queue(model, "queue1", SchedStrategy.FCFS);
        node1.setNumberOfServers(12);
        Queue node2 = new Queue(model, "queue2", SchedStrategy.FCFS);
        node2.setNumberOfServers(40);
        Delay node3 = new Delay(model, "delay1");
        Delay node4 = new Delay(model, "delay2");
        Delay node5 = new Delay(model, "delay3");
        Delay node6 = new Delay(model, "delay4");
        Delay node7 = new Delay(model, "delay5");
        Delay node8 = new Delay(model, "delay6");
        Delay node9 = new Delay(model, "delay7");
        Delay node10 = new Delay(model, "delay8");
        Delay node11 = new Delay(model, "delay9");
        Delay node12 = new Delay(model, "delay10");
        List<ClassSwitch> switchNodes = Arrays.asList(
                new ClassSwitch(model, "cs_queue2_delay9"),
                new ClassSwitch(model, "cs_delay1_queue1"),
                new ClassSwitch(model, "cs_delay2_delay5"),
                new ClassSwitch(model, "cs_delay3_delay2"),
                new ClassSwitch(model, "cs_delay5_delay8"),
                new ClassSwitch(model, "cs_delay6_delay5"),
                new ClassSwitch(model, "cs_delay7_delay8"),
                new ClassSwitch(model, "cs_delay9_queue1"),
                new ClassSwitch(model, "cs_delay9_delay4"),
                new ClassSwitch(model, "cs_delay10_delay3")
        );
        //Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "CClass1", 20, node8, 0);

        node1.setService(jobclass1, Exp.fitMean(64.000000)); // (queue1,CClass1)
        node2.setService(jobclass1, Exp.fitMean(0.015625)); // (queue2,CClass1)
        node3.setService(jobclass1, APH.fitMeanAndSCV(1.000000, 16.000000)); // (delay1,CClass1)
        node4.setService(jobclass1, APH.fitMeanAndSCV(2.000000, 32.000000)); // (delay2,CClass1)
        node5.setService(jobclass1, APH.fitMeanAndSCV(0.015625, 2.000000)); // (delay3,CClass1)
        node6.setService(jobclass1, APH.fitMeanAndSCV(1.000000, 0.500000)); // (delay4,CClass1)
        node7.setService(jobclass1, new Erlang(128.000000, 64)); // (delay5,CClass1)
        node8.setService(jobclass1, APH.fitMeanAndSCV(0.125000, 4.000000)); // (delay6,CClass1)
        node9.setService(jobclass1, Exp.fitMean(0.031250)); // (delay7,CClass1)
        node10.setService(jobclass1, Exp.fitMean(0.250000)); // (delay8,CClass1)
        node11.setService(jobclass1, APH.fitMeanAndSCV(0.031250, 4.000000)); // (delay9,CClass1)
        node12.setService(jobclass1, APH.fitMeanAndSCV(16.000000, 16.000000)); // (delay10,CClass1)
        // Initialise class switch matrix
        for (ClassSwitch switchNode : switchNodes) {
            ClassSwitchMatrix csMatrix = switchNode.initClassSwitchMatrix();
            csMatrix.set(jobclass1.getIndex() - 1, jobclass1.getIndex() - 1, 1);
            switchNode.setClassSwitchingMatrix(csMatrix);
        }


        // Block 3: topology
        RoutingMatrix P = model.initRoutingMatrix();

        P.set(jobclass1, jobclass1, node1, node7, 1);
        P.set(jobclass1, jobclass1, node2, node9, 5.000000e-01);
        P.set(jobclass1, jobclass1, node2, switchNodes.get(0), 5.000000e-01);
        P.set(jobclass1, jobclass1, node3, switchNodes.get(1), 1);
        P.set(jobclass1, jobclass1, node4, switchNodes.get(2), 1);
        P.set(jobclass1, jobclass1, node5, node3, 5.000000e-01);
        P.set(jobclass1, jobclass1, node5, switchNodes.get(3), 5.000000e-01);
        P.set(jobclass1, jobclass1, node6, node1, 1);
        P.set(jobclass1, jobclass1, node7, node8, 5.000000e-01);
        P.set(jobclass1, jobclass1, node7, switchNodes.get(4), 5.000000e-01);
        P.set(jobclass1, jobclass1, node8, switchNodes.get(5), 1);
        P.set(jobclass1, jobclass1, node9, switchNodes.get(6), 1);
        P.set(jobclass1, jobclass1, node10, node2, 8.480000e-01);
        P.set(jobclass1, jobclass1, node10, node12, 1.520000e-01);
        P.set(jobclass1, jobclass1, node11, switchNodes.get(7), 5.000000e-01);
        P.set(jobclass1, jobclass1, node11, switchNodes.get(8), 5.000000e-01);
        P.set(jobclass1, jobclass1, node12, switchNodes.get(9), 1);
        P.set(jobclass1, jobclass1, switchNodes.get(0), node11, 1);
        P.set(jobclass1, jobclass1, switchNodes.get(1), node1, 1);
        P.set(jobclass1, jobclass1, switchNodes.get(2), node7, 1);
        P.set(jobclass1, jobclass1, switchNodes.get(3), node4, 1);
        P.set(jobclass1, jobclass1, switchNodes.get(4), node10, 1);
        P.set(jobclass1, jobclass1, switchNodes.get(5), node7, 1);
        P.set(jobclass1, jobclass1, switchNodes.get(6), node10, 1);
        P.set(jobclass1, jobclass1, switchNodes.get(7), node1, 1);
        P.set(jobclass1, jobclass1, switchNodes.get(8), node6, 1);
        P.set(jobclass1, jobclass1, switchNodes.get(9), node5, 1);
        model.link(P);
        return model;
    }

    public static Network closed2() {
        Network model = new Network("MRP");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);
        ClosedClass closedClass = new ClosedClass(model, "Machines", 1, delay, 0);
        ClosedClass closedClass2 = new ClosedClass(model, "Machines2", 1, delay, 0);
        delay.setService(closedClass, new Exp(1));
        delay.setService(closedClass2, new Exp(1));

        queue.setService(closedClass, new Exp(1));
        queue.setService(closedClass2, new Exp(1));
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(closedClass, closedClass, delay, delay, 0.3);
        routingMatrix.set(closedClass, closedClass, delay, queue, 0.1);
        routingMatrix.set(closedClass, closedClass, queue, delay, 0.2);
        routingMatrix.set(closedClass, closedClass, queue, queue, 0);
        routingMatrix.set(closedClass, closedClass2, delay, delay, 0.6);
        routingMatrix.set(closedClass, closedClass2, delay, queue, 0.0);
        routingMatrix.set(closedClass, closedClass2, queue, delay, 0.8);
        routingMatrix.set(closedClass, closedClass2, queue, queue, 0.0);
        routingMatrix.set(closedClass2, closedClass2, delay, delay, 0.0);
        routingMatrix.set(closedClass2, closedClass2, delay, queue, 1.0);
        routingMatrix.set(closedClass2, closedClass2, queue, delay, 0.0);
        routingMatrix.set(closedClass2, closedClass2, queue, queue, 0.0);
        routingMatrix.set(closedClass2, closedClass, delay, delay, 0.0);
        routingMatrix.set(closedClass2, closedClass, delay, queue, 0.0);
        routingMatrix.set(closedClass2, closedClass, queue, delay, 1.0);
        routingMatrix.set(closedClass2, closedClass, queue, queue, 0.0);
        model.link(routingMatrix);
        return model;
    }

    public static Network mixed1() {
        Network model = new Network("MRP");
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Source node3 = new Source(model, "Source");
        Sink node4 = new Sink(model, "Sink");
        ClosedClass jobclass1 = new ClosedClass(model, "ClosedClass", 2, node1, 0);
        OpenClass jobclass2 = new OpenClass(model, "myClass", 0);
        node1.setService(jobclass1, new Exp(1));
        node1.setService(jobclass2, new Exp(1));
        node2.setService(jobclass1, new Erlang(3, 2));
        node2.setService(jobclass2, new Exp(1));

        node3.setArrival(jobclass2, new Exp(0.1));
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node1, node3, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node1, node4, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node2, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node2, node4, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node3, node1, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node3, node2, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node3, node3, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node4, node2, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node4, node3, 0.0);
        routingMatrix.set(jobclass1, jobclass1, node4, node4, 0.0);

        routingMatrix.set(jobclass1, jobclass2, node1, node1, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node1, node2, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node1, node3, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node1, node4, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node2, node1, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node2, node3, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node2, node4, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node3, node1, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node3, node2, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node3, node3, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node3, node4, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node4, node1, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node4, node2, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node4, node3, 0.0);
        routingMatrix.set(jobclass1, jobclass2, node4, node4, 0.0);

        routingMatrix.set(jobclass2, jobclass1, node1, node1, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node1, node2, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node1, node3, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node1, node4, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node2, node1, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node2, node2, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node2, node3, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node2, node4, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node3, node1, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node3, node2, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node3, node3, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node3, node4, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node4, node1, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node4, node2, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node4, node3, 0.0);
        routingMatrix.set(jobclass2, jobclass1, node4, node4, 0.0);

        routingMatrix.set(jobclass2, jobclass2, node1, node1, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node3, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node1, node4, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node2, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node2, node4, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node3, node1, 1.0);
        routingMatrix.set(jobclass2, jobclass2, node3, node2, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node3, node3, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node3, node4, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node4, node1, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node4, node2, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node4, node3, 0.0);
        routingMatrix.set(jobclass2, jobclass2, node4, node4, 0.0);
        model.link(routingMatrix);
        return model;
    }

    public static LayeredNetwork testSimple() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_simple");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF);
        T1.setThinkTime(new Exp(0.01));
        T1.on(P1);

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp(1));
        A1.on(T1);
        A1.boundTo(E1);
//        T1.addPrecedence(ActivityPrecedence.Sequence("A1","A2"));
        return model;
    }

    public static LayeredNetwork test0() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_1");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF);
        Task T2 = new Task(model, "T2", 50, SchedStrategy.FCFS);
        T1.setThinkTime(new Exp(1));
        T1.on(P1);
        T2.setThinkTime(new Exp(1));
        T2.on(P2);

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);


        Activity A1 = new Activity(model, "A1", new Exp(1));
        A1.on(T1);
        A1.boundTo(E1);
        A1.synchCall(E2, 1);

        Activity A2 = new Activity(model, "A2", new Exp(1));
        A2.on(T2);
        A2.boundTo(E2);
        A2.repliesTo(E2);
//        T1.addPrecedence(ActivityPrecedence.Sequence("A1","A2"));
        return model;
    }

    public static LayeredNetwork test1() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_1");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF);
        T1.setThinkTime(new Exp(1.0 / 100));
        T1.on(P1);

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp(10));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Exp(1.0 / 1.5));
        A2.on(T1);

        T1.addPrecedence(ActivityPrecedence.Sequence("A1", "A2"));

        return model;
    }

    public static LayeredNetwork test2() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_2");

        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Exp(1.0 / 100));

        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS);
        T2.on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Entry E2 = new Entry(model, "E2");
        E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1 / 1.6));
        A1.on(T1);
        A1.boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Immediate());
        A2.on(T1);

        Activity A3 = new Activity(model, "A3", new Immediate());
        A3.on(T1);
        A3.synchCall(E2, 1);

        Activity A4 = new Activity(model, "A4", new Exp(1.0 / 5));
        A4.on(T2);
        A4.boundTo(E2);

        Activity A5 = new Activity(model, "A5", new Exp(1));
        A5.on(T2);
        A5.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Sequence("A1", "A2"));
        T1.addPrecedence(ActivityPrecedence.Sequence("A2", "A3"));
        T2.addPrecedence(ActivityPrecedence.Sequence("A4", "A5"));

        return model;
    }

    public static LayeredNetwork test35() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF); T1.on(P1);
        T1.setThinkTime(Erlang.fitMeanAndSCV(0.0001,0.5));
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF); T2.on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1"); E1.on(T1);
        Entry E2 = new Entry(model, "E2"); E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1)); A1.on(T1); A1.boundTo(E1); A1.synchCall(E2,3);
        Activity A2 = new Activity(model, "A2", APH.fitMeanAndSCV(1,10)); A2.on(T2); A2.boundTo(E2); A2.repliesTo(E2);

        return model;
    }

    public static LayeredNetwork test3() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF); T1.on(P1);
        T1.setThinkTime(Erlang.fitMeanAndSCV(0.0001,0.5));
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF); T2.on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1"); E1.on(T1);
        Entry E2 = new Entry(model, "E2"); E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp(1)); A1.on(T1); A1.boundTo(E1); A1.synchCall(E2,3);
        Activity A2 = new Activity(model, "A2", APH.fitMeanAndSCV(1.0,10.0)); A2.on(T2); A2.boundTo(E2); A2.repliesTo(E2);

        return model;
    }

    public static LayeredNetwork test4() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_LQN_4");
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Exp(100));

        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2);
        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);

        Entry E2 = new Entry(model, "E2").on(T2);

        Activity A1 = new Activity(model, "A1", new Exp( 1.6)).on(T1).boundTo(E1);

        Activity A2 = new Activity(model, "A2", new Immediate()).on(T1).synchCall(E2);

        Activity A3 = new Activity(model, "A3", new Exp(5)).on(T2).boundTo(E2);

        Activity A4 = new Activity(model, "A4", new Exp(1)).on(T2).repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Sequence("A1", "A2"));
        T2.addPrecedence(ActivityPrecedence.Sequence("A3", "A4"));

        return model;
    }

    public static LayeredNetwork testAndForkJoin() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_and_fork_join");
        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp((double) 1/2)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", new Exp((double)1/3)).on(T1);
        Activity A3 = new Activity(model, "A3", new Exp((double)1/4)).on(T1);
        Activity A4 = new Activity(model, "A4", new Exp((double)1/5)).on(T1);
        Activity A5 = new Activity(model, "A5", new Exp((double)1/6)).on(T1);

        T1.addPrecedence(ActivityPrecedence.AndFork("A1", Arrays.asList("A2", "A3", "A4"), new Matrix(0, 0)));
        T1.addPrecedence(ActivityPrecedence.AndJoin(Arrays.asList("A2", "A3", "A4"), "A5", new Matrix(1)));

        return model;
    }

    public static LayeredNetwork testOrForkJoin() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_or_fork_join");
        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.FCFS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1");
        E1.on(T1);

        Activity A1 = new Activity(model, "A1", new Exp((double) 1/2)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", new Exp((double)1/3)).on(T1);
        Activity A3 = new Activity(model, "A3", new Exp((double)1/4)).on(T1);
        Activity A4 = new Activity(model, "A4", new Exp((double)1/5)).on(T1);
        Activity A5 = new Activity(model, "A5", new Exp((double)1/6)).on(T1);

        Matrix params = new Matrix(1, 3);
        params.set(0,0,0.3);
        params.set(0,1,0.3);
        params.set(0,2,0.4);
        T1.addPrecedence(ActivityPrecedence.OrFork("A1", Arrays.asList("A2", "A3", "A4"),  params));
        T1.addPrecedence(ActivityPrecedence.OrJoin(Arrays.asList("A2", "A3", "A4"), "A5"));

        return model;
    }

    public static LayeredNetwork testLoop() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_loop");
        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.INF);
//        Processor P2 = new Processor(model, "P2", 10, SchedStrategy.INF);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF);
        T1.on(P1);
        T1.setThinkTime(new Immediate());
//
//        Task T2 = new Task(model, "T2", 1, SchedStrategy.INF);
//        T2.on(P2);
//        T2.setThinkTime(new Immediate());

        Entry E1 = new Entry(model, "E1").on(T1);

//        Entry E2 = new Entry(model, "E2");
//        E2.on(T2);

        Activity A1 = new Activity(model, "A1", new Exp( (double)1/2)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", new Exp((double) 1/3)).on(T1);
        Activity A3 = new Activity(model, "A3", new Exp((double) 1/4)).on(T1);

//        A3.synchCall(E2);

//        Activity B1 = new Activity(model, "B1", new Exp( 1/2));
//        B1.on(T2);
//        B1.boundTo(E2);

//        Activity B2 = new Activity(model, "B2", new Exp(1/3));
//        B2.on(T2);
//        B2.repliesTo(E2);

        T1.addPrecedence(ActivityPrecedence.Loop("A1", Arrays.asList("A2", "A3"), new Matrix(2)));
//        T1.addPrecedence(ActivityPrecedence.Sequence("B1", "B2"));

        return model;
    }

    public static LayeredNetwork testAllPrecedences() throws Exception {
        LayeredNetwork model = new LayeredNetwork("test_loop_network");
        Processor P1 = new Processor(model, "P1", 10, SchedStrategy.INF);
        Processor P2 = new Processor(model, "P2", 10, SchedStrategy.INF);
        Processor P3 = new Processor(model, "P3", 5, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1).setThinkTime(new Immediate());
        Task T2 = new Task(model, "T2", 1, SchedStrategy.INF).on(P2).setThinkTime(new Immediate());
        Task T3 = new Task(model, "T3", 20, SchedStrategy.INF).on(P3).setThinkTime(new Exp((double) 1 / 10));

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T3);

        Activity A1 = new Activity(model, "A1", new Exp( 1)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", new Exp(1/2)).on(T1);
        Activity A3 = new Activity(model, "A3", new Exp(1/3)).on(T1).synchCall(E2);

        Activity B1 = new Activity(model, "B1", new Exp( 10)).on(T2).boundTo(E2);
        Activity B2 = new Activity(model, "B2", new Exp(5)).on(T2);
        Activity B3 = new Activity(model, "B3", new Exp(10/3)).on(T2);
        Activity B4 = new Activity(model, "B4", new Exp(10/4)).on(T2);
        Activity B5 = new Activity(model, "B5", new Exp(2)).on(T2);
        Activity B6 = new Activity(model, "B6", new Exp(10/6)).on(T2).synchCall(E3).repliesTo(E2);

        Activity C1 = new Activity(model, "C1", new Exp( 10)).on(T3).boundTo(E3);
        Activity C2 = new Activity(model, "C2", new Exp(5)).on(T3);
        Activity C3 = new Activity(model, "C3", new Exp(10/3)).on(T3);
        Activity C4 = new Activity(model, "C4", new Exp(10/4)).on(T3);
        Activity C5 = new Activity(model, "C5", new Exp(2)).on(T3).repliesTo(E3);

        T1.addPrecedence(ActivityPrecedence.Loop("A1", Arrays.asList("A2", "A3"), new Matrix(3)));
        T2.addPrecedence(ActivityPrecedence.Sequence("B4", "B5"));
        T2.addPrecedence(ActivityPrecedence.AndFork("B1", Arrays.asList("B2", "B3", "B4"), new Matrix(0, 0)));
        T2.addPrecedence(ActivityPrecedence.AndJoin(Arrays.asList("B2", "B3", "B5"), "B6", new Matrix(0, 0)));
        T3.addPrecedence(ActivityPrecedence.OrFork("C1", Arrays.asList("C2", "C3", "C4"), new Matrix(Arrays.asList(0.3, 0.3, 0.4))));
        T3.addPrecedence(ActivityPrecedence.OrJoin(Arrays.asList("C2", "C3", "C4"), "C5"));

        return model;
    }

    public static LayeredNetwork  testOfBizFCFS() throws Exception {

        LayeredNetwork model = new LayeredNetwork("myLayeredModel");

        Processor P1 = new Processor(model, "FrontEnd_CPU_Processor", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "USAGE_DELAY_Processor", 1, SchedStrategy.PS);
        Processor P3 = new Processor(model, "UsageScenario_userType1_1_Processor", 1, SchedStrategy.PS);
        Processor P4 = new Processor(model, "RequestHandler_HandlerIF_main_345_Processor", 1, SchedStrategy.PS);
        Processor P5 = new Processor(model, "RequestHandler_HandlerIF_login_345_Processor", 1, SchedStrategy.PS);
        Processor P6 = new Processor(model, "RequestHandler_HandlerIF_checkLogin_345_Processor", 1, SchedStrategy.PS);
        Processor P7 = new Processor(model, "RequestHandler_HandlerIF_logout_345_Processor", 1, SchedStrategy.PS);
        Processor P8 = new Processor(model, "UsageScenario_userType2_7_Processor", 1, SchedStrategy.PS);
        Processor P9 = new Processor(model, "RequestHandler_HandlerIF_quickadd_345_Processor", 1, SchedStrategy.PS);

        Task T1 = new Task(model, "FrontEnd_CPU_Task", 1, SchedStrategy.FCFS).on(P1).setThinkTime(new Exp(1/5.0));
        Task T2 = new Task(model, "USAGE_DELAY_Task", 1, SchedStrategy.FCFS).on(P2).setThinkTime(new Exp(1/5.0));
        Task T3 = new Task(model, "UsageScenario_userType1_1_Task", 10, SchedStrategy.REF).on(P3).setThinkTime(new Exp(0.1));
        Task T4 = new Task(model, "RequestHandler_HandlerIF_main_345_Task", 1, SchedStrategy.FCFS).on(P4).setThinkTime(new Exp(1.0/2));
        Task T5 = new Task(model, "RequestHandler_HandlerIF_login_345_Task", 1, SchedStrategy.FCFS).on(P5).setThinkTime(new Exp(1));
        Task T6 = new Task(model, "RequestHandler_HandlerIF_checkLogin_345_Task", 1, SchedStrategy.FCFS).on(P6).setThinkTime(new Exp(1.0/2));
        Task T7 = new Task(model, "RequestHandler_HandlerIF_logout_345_Task", 1, SchedStrategy.FCFS).on(P7).setThinkTime(new Exp(1.0/5));
        Task T8 = new Task(model, "UsageScenario_userType2_7_Task", 10, SchedStrategy.REF).on(P8).setThinkTime(new Exp(0.0714286));
        Task T9 = new Task(model, "RequestHandler_HandlerIF_quickadd_345_Task", 1, SchedStrategy.FCFS).on(P9).setThinkTime(new Exp(1.0/10));

        Entry E1 = new Entry(model, "FrontEnd_CPU_Entry").on(T1);
        Entry E2 = new Entry(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50_Entry").on(T1);
        Entry E3 = new Entry(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50_Entry").on(T1);
        Entry E4 = new Entry(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50_Entry").on(T1);
        Entry E5 = new Entry(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50_Entry").on(T1);
        Entry E6 = new Entry(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50_Entry").on(T1);
        Entry E7 = new Entry(model, "USAGE_DELAY0_Entry").on(T2);
        Entry E8 = new Entry(model, "UsageScenario_userType1_1_Entry").on(T3);
        Entry E9 = new Entry(model, "RequestHandler_HandlerIF_main_345_Entry").on(T4);
        Entry E10 = new Entry(model, "RequestHandler_HandlerIF_login_345_Entry").on(T5);
        Entry E11 = new Entry(model, "RequestHandler_HandlerIF_checkLogin_345_Entry").on(T6);
        Entry E12 = new Entry(model, "RequestHandler_HandlerIF_logout_345_Entry").on(T7);
        Entry E13 = new Entry(model, "UsageScenario_userType2_7_Entry").on(T8);
        Entry E14 = new Entry(model, "RequestHandler_HandlerIF_quickadd_345_Entry").on(T9);

        Activity A1 = new Activity(model, "FrontEnd_CPU_Activity", new Exp(1.0)).on(T1).boundTo(E1).repliesTo(E1);
        Activity A2 = new Activity(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100)).on(T1).boundTo(E2).repliesTo(E2);
        Activity A3 = new Activity(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50_Activity", new Exp(100)).on(T1).boundTo(E3).repliesTo(E3);
        Activity A4 = new Activity(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100)).on(T1).boundTo(E4).repliesTo(E4);
        Activity A5 = new Activity(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50_Activity", new Exp(100)).on(T1).boundTo(E5).repliesTo(E5);
        Activity A6 = new Activity(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50_Activity", new Exp(100)).on(T1).boundTo(E6).repliesTo(E6);
        Activity A7 = new Activity(model, "USAGE_DELAY0_Activity", new Exp(1.0/5)).on(T2).boundTo(E7).repliesTo(E7);
        Activity A8 = new Activity(model, "Start2", new Exp(1.0/3)).on(T3).boundTo(E8);
        Activity A9 = new Activity(model, "EntryLevelSystemCallcall1main1", new Exp(1.0/2)).on(T3).synchCall(E9,1);
        Activity A10 = new Activity(model, "EntryLevelSystemCallcall1login", new Exp(1.0/5)).on(T3).synchCall(E10,1);
        Activity A11 = new Activity(model, "EntryLevelSystemCallcall1checkLogin1", new Exp(1.0/5)); A11.on(T3); A11.synchCall(E11,1);
        Activity A12 = new Activity(model, "EntryLevelSystemCallcall1checkLogin2", new Exp(1.0/5)); A12.on(T3); A12.synchCall(E11,1);
        Activity A13 = new Activity(model, "EntryLevelSystemCallcall1logout", new Exp(1.0/5)); A13.on(T3); A13.synchCall(E12,1);
        Activity A14 = new Activity(model, "EntryLevelSystemCallcall1main2", new Exp(1.0/5)); A14.on(T3); A14.synchCall(E9,1);
        Activity A15 = new Activity(model, "Stop6", new Exp(1.0/5)); A15.on(T3);
        Activity A16 = new Activity(model, "StartAction_start__EkLVIMhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)).on(T4).boundTo(E9);
        Activity A17 = new Activity(model, "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50", new Exp(1.0/5)); A17.on(T4).synchCall(E2,1);
        Activity A18 = new Activity(model, "StopAction_stop__EkL8MMhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A18.on(T4); A18.repliesTo(E9);
        Activity A19 = new Activity(model, "StartAction_aName__XJKk8g26EeSPwb7XgvxhWQ_34_5", new Exp(1.0/5)); A19.on(T5); A19.boundTo(E10);
        Activity A20 = new Activity(model, "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50", new Exp(1.0/5)); A20.on(T5); A20.synchCall(E3,1);
        Activity A21 = new Activity(model, "StopAction_aName__ae0SIA26EeSPwb7XgvxhWQ_34_5", new Exp(1.0/5)); A21.on(T5); A21.repliesTo(E10);
        Activity A22 = new Activity(model, "StartAction_start__EkI44MhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A22.on(T6); A22.boundTo(E11);
        Activity A23 = new Activity(model, "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50", new Exp(1.0/5)); A23.on(T6); A23.synchCall(E4,1);
        Activity A24 = new Activity(model, "StopAction_stop__EkJf8MhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A24.on(T6); A24.repliesTo(E11);
        Activity A25 = new Activity(model, "StartAction_start__EkVtMMhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A25.on(T7); A25.boundTo(E12);
        Activity A26 = new Activity(model, "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50", new Exp(1.0/5)); A26.on(T7); A26.synchCall(E5,1);
        Activity A27 = new Activity(model, "StopAction_stop__EkVtMchoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A27.on(T7); A27.repliesTo(E12);
        Activity A28 = new Activity(model, "Start8", new Exp(1.0/5)); A28.on(T8); A28.boundTo(E13);
        Activity A29 = new Activity(model, "EntryLevelSystemCallcall2main1", new Exp(1.0/5)); A29.on(T8); A29.synchCall(E9,1);
        Activity A30 = new Activity(model, "EntryLevelSystemCallcall2login", new Exp(1.0/5)); A30.on(T8); A30.synchCall(E10,1);
        Activity A31 = new Activity(model, "EntryLevelSystemCallcall2checkLogin1", new Exp(1.0/5)); A31.on(T8); A31.synchCall(E11,1);
        Activity A32 = new Activity(model, "EntryLevelSystemCallcall2checkLogin2", new Exp(1.0/5)); A32.on(T8); A32.synchCall(E11,1);
        Activity A33 = new Activity(model, "EntryLevelSystemCallcall2main2", new Exp(1.0/5)); A33.on(T8); A33.synchCall(E9,1);
        Activity A34 = new Activity(model, "EntryLevelSystemCallcall2quickadd", new Exp(1.0/5)); A34.on(T8); A34.synchCall(E14,1);
        Activity A35 = new Activity(model, "EntryLevelSystemCallcall2logout", new Exp(1.0/5)); A35.on(T8); A35.synchCall(E12,1);
        Activity A36 = new Activity(model, "EntryLevelSystemCallcall2main3", new Exp(1.0/5)); A36.on(T8); A36.synchCall(E9,1);
        Activity A37 = new Activity(model, "Stop9", new Exp(1.0/5)); A37.on(T8);
        Activity A38 = new Activity(model, "StartAction_start__EkBkIMhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A38.on(T9); A38.boundTo(E14);
        Activity A39 = new Activity(model, "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50", new Exp(1.0/5)); A39.on(T9); A39.synchCall(E6,1);
        Activity A40 = new Activity(model, "StopAction_stop__EkCLMMhoEeKON4DtRoKCMw_34_5", new Exp(1.0/5)); A40.on(T9); A40.repliesTo(E14);

        T3.addPrecedence(ActivityPrecedence.Sequence("Start2", "EntryLevelSystemCallcall1main1"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1main1", "EntryLevelSystemCallcall1login"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1login", "EntryLevelSystemCallcall1checkLogin1"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1checkLogin1", "EntryLevelSystemCallcall1checkLogin2"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1checkLogin2", "EntryLevelSystemCallcall1logout"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1logout", "EntryLevelSystemCallcall1main2"));
        T3.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall1main2", "Stop6"));
        T4.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkLVIMhoEeKON4DtRoKCMw_34_5", "InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50"));
        T4.addPrecedence(ActivityPrecedence.Sequence("InternalAction_main__Iu1-wMhoEeKON4DtRoKCMw_34_50", "StopAction_stop__EkL8MMhoEeKON4DtRoKCMw_34_5"));
        T5.addPrecedence(ActivityPrecedence.Sequence("StartAction_aName__XJKk8g26EeSPwb7XgvxhWQ_34_5", "InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50"));
        T5.addPrecedence(ActivityPrecedence.Sequence("InternalAction_login__YgxRGw26EeSPwb7XgvxhWQ_34_50", "StopAction_aName__ae0SIA26EeSPwb7XgvxhWQ_34_5"));
        T6.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkI44MhoEeKON4DtRoKCMw_34_5", "InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50"));
        T6.addPrecedence(ActivityPrecedence.Sequence("InternalAction_CheckLogin__vHMZsMhoEeKON4DtRoKCMw_34_50", "StopAction_stop__EkJf8MhoEeKON4DtRoKCMw_34_5"));
        T7.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkVtMMhoEeKON4DtRoKCMw_34_5", "InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50"));
        T7.addPrecedence(ActivityPrecedence.Sequence("InternalAction_logout__j-2E4MhrEeKON4DtRoKCMw_34_50", "StopAction_stop__EkVtMchoEeKON4DtRoKCMw_34_5"));
        T8.addPrecedence(ActivityPrecedence.Sequence("Start8", "EntryLevelSystemCallcall2main1"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2main1", "EntryLevelSystemCallcall2login"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2login", "EntryLevelSystemCallcall2checkLogin1"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2checkLogin1", "EntryLevelSystemCallcall2checkLogin2"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2checkLogin2", "EntryLevelSystemCallcall2main2"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2main2", "EntryLevelSystemCallcall2quickadd"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2quickadd", "EntryLevelSystemCallcall2logout"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2logout", "EntryLevelSystemCallcall2main3"));
        T8.addPrecedence(ActivityPrecedence.Sequence("EntryLevelSystemCallcall2main3", "Stop9"));
        T9.addPrecedence(ActivityPrecedence.Sequence("StartAction_start__EkBkIMhoEeKON4DtRoKCMw_34_5", "InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50"));
        T9.addPrecedence(ActivityPrecedence.Sequence("InternalAction_addcart__5JEHQMhoEeKON4DtRoKCMw_34_50", "StopAction_stop__EkCLMMhoEeKON4DtRoKCMw_34_5"));

//        SolverLN solver = new SolverLN(model);
//        solver.getEnsemble().get(7).getConnectionMatrix().print();
//        for(int i=0;i<solver.getEnsemble().size();i++) {
//            System.out.println("Solver result "+i);
//            Network layer = solver.getEnsemble().get(i);
//            SolverMVA layersolver = new SolverMVA(layer);
//            layersolver.getAvgTable();
//        }
//        solver.getEnsembleAvg();
        return model;
    }

    public static Network test_mixedModel_1() {
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
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

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

    // For demonstration of State-Dependent Random Environments
    public static SolverEnv test_randomEnvironment_4() {

      int E = 2;
      Env envModel = new Env("MyEnv", E);
      String[] envName = {"Stage1", "Stage2"};
      String[] envType = {"FAST", "SLOW"};

      // Model in Stage 1
      Network modelStage1 = new Network("model");
      Queue queue1Stage1 = new Queue(modelStage1, "Queue1", SchedStrategy.PS);
      Queue queue2Stage1 = new Queue(modelStage1, "Queue2", SchedStrategy.PS);
      ClosedClass class1Stage1 = new ClosedClass(modelStage1, "Class1", 8, queue1Stage1, 0);
      queue1Stage1.setService(class1Stage1, new Exp(100));
      queue2Stage1.setService(class1Stage1, new Exp(10));
      modelStage1.link(modelStage1.serialRouting(queue1Stage1, queue2Stage1));
      envModel.addStage(0, envName[0], envType[0], modelStage1);

      // Model in Stage 2
      Network modelStage2 = new Network("model");
      Queue queue1Stage2 = new Queue(modelStage2, "Queue1", SchedStrategy.FCFS);
      Queue queue2Stage2 = new Queue(modelStage2, "Queue2", SchedStrategy.FCFS);
      ClosedClass class1Stage2 = new ClosedClass(modelStage2, "Class1", 8, queue1Stage2, 0);
      queue1Stage2.setService(class1Stage2, new Exp(1));
      queue2Stage2.setService(class1Stage2, new Exp(10));
      modelStage2.link(modelStage2.serialRouting(queue1Stage2, queue2Stage2));
      envModel.addStage(1, envName[1], envType[1], modelStage2);

      Matrix envRates = new Matrix(2, 2);
      envRates.set(0, 1, 100);
      envRates.set(1, 0, 0.01);

      Env.ResetEnvRatesFunction resetEnvRatesFunction =
          (originalDist, QExit, UExit, TExit) -> {
            double lambda = originalDist.getRate();
            lambda *= UExit.sumRows(0); // Time-averaged utilisation at Queue1
            return new Exp(Math.max(lambda, GlobalConstants.Zero));
          };

      for (int e = 0; e < E; e++) {
        for (int h = 0; h < E; h++) {
          if (envRates.get(e, h) > 0) {
            envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
          }
        }
      }
      envModel.resetEnvRatesFun[0][1] = resetEnvRatesFunction;

      SolverOptions options = new SolverOptions(SolverType.ENV);
      options.iter_tol = 0.01;
      options.timespan[0] = 0;
      options.method = "statedep";

      SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
      fluidOptions.method = "matrix";
      fluidOptions.stiff = false;
      fluidOptions.setODEMaxStep(0.1);
      fluidOptions.verbose = VerboseLevel.SILENT;

      NetworkSolver[] solvers = new NetworkSolver[E];
      for (int e = 0; e < E; e++) {
        solvers[e] = new SolverFluid(envModel.getModel(e));
        solvers[e].options = fluidOptions;
      }

      return new SolverEnv(envModel, solvers, options);
    }

    // For demonstration of p-Norm Smoothing combined with SolverEnv
    public static SolverEnv test_randomEnvironment_5() {

      int E = 2;
      Env envModel = new Env("MyEnv", E);
      String[] envName = {"Stage1", "Stage2"};
      String[] envType = {"SLOW", "FAST"};

      // Model in Stage 1
      Network modelStage1 = new Network("model");
      Delay delayStage1 = new Delay(modelStage1, "Delay");
      Queue queueStage1 = new Queue(modelStage1, "Queue", SchedStrategy.PS);
      ClosedClass class1Stage1 = new ClosedClass(modelStage1, "Class1", 8, delayStage1, 0);
      delayStage1.setService(class1Stage1, new Exp(1));
      queueStage1.setService(class1Stage1, new Exp(8));
      modelStage1.link(modelStage1.serialRouting(delayStage1, queueStage1));
      envModel.addStage(0, envName[0], envType[0], modelStage1);

      // Model in Stage 2
      Network modelStage2 = new Network("model");
      Delay delayStage2 = new Delay(modelStage2, "Delay");
      Queue queueStage2 = new Queue(modelStage2, "Queue", SchedStrategy.PS);
      ClosedClass class1Stage2 = new ClosedClass(modelStage2, "Class1", 8, delayStage2, 0);
      delayStage2.setService(class1Stage2, new Exp(4));
      queueStage2.setService(class1Stage2, new Exp(8));
      modelStage2.link(modelStage2.serialRouting(delayStage2, queueStage2));
      envModel.addStage(1, envName[1], envType[1], modelStage2);

      Matrix envRates = new Matrix(2, 2);
      envRates.set(0, 1, 0.001); // Replace with 1000 to test Figure 6.11 results
      envRates.set(1, 0, 0.001); // Replace with 1000 to test Figure 6.11 results

      for (int e = 0; e < E; e++) {
        for (int h = 0; h < E; h++) {
          if (envRates.get(e, h) > 0) {
            envModel.addTransition(e, h, new Exp(envRates.get(e, h)));
          }
        }
      }

      SolverOptions options = new SolverOptions(SolverType.ENV);
      options.iter_tol = 0.01;
      options.timespan[0] = 0;

      SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
      fluidOptions.method = "matrix";
      fluidOptions.stiff = false;
      fluidOptions.setODEMaxStep(0.1);
      fluidOptions.verbose = VerboseLevel.SILENT;

      NetworkSolver[] solvers = new NetworkSolver[E];
      for (int e = 0; e < E; e++) {
        solvers[e] = new SolverFluid(envModel.getModel(e));
        solvers[e].options = fluidOptions;
        // Activating p-Norm Smoothing by providing each solver with an initial set of pStar values.
        // Going forwards, the pStar values are refreshed within SolverEnv's 'analyze' method
  //      PStarSearcher searcher = new PStarSearcher();
  //      Matrix targetQueueLengths = searcher.generateTargetQueueLengths(solvers[e].model);
  //      PointValuePair pStarValues = searcher.findPStarValues(solvers[e].model, targetQueueLengths);
  //      for (int i = 0; i < solvers[e].model.getNumberOfNodes(); i++) {
  //        solvers[e].options.config.pstar.add(i, pStarValues.getPoint()[i]);
  //      }
      }

      return new SolverEnv(envModel, solvers, options);
    }

    // node 3 arrival time matches matlab, ex3_line kept as used in other tests
    public static Network test_openModel_1(){
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

    public static Network test_openModel_1b() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source 1");
        Queue node2 = new Queue(model, "Queue 1", SchedStrategy.PS);
        ClassSwitch node3 = new ClassSwitch(model, "ClassSwitch 1"); // Dummy node, class switching is embedded in the routing matrix P
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

        // Block 3: topology
        ClassSwitchMatrix C = node3.initClassSwitchMatrix();
        C.setTo(Matrix.eye(C.length()));
        node3.setClassSwitchingMatrix(C);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Source 1,Class A) -> (Queue 1,Class A)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue 1,Class A) -> (ClassSwitch 1,Class A)
        routingMatrix.set(jobclass1, jobclass1, node5, node4, 1.000000); // (Queue 2,Class A) -> (Sink 1,Class A)
        routingMatrix.set(jobclass1, jobclass3, node3, node5, 1.000000); // (ClassSwitch 1,Class A) -> (Queue 2,Class C)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Source 1,Class B) -> (Queue 1,Class B)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue 1,Class B) -> (ClassSwitch 1,Class B)
        routingMatrix.set(jobclass2, jobclass2, node5, node4, 1.000000); // (Queue 2,Class B) -> (Sink 1,Class B)
        routingMatrix.set(jobclass2, jobclass3, node3, node5, 1.000000); // (ClassSwitch 1,Class B) -> (Queue 2,Class C)
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.000000); // (Source 1,Class C) -> (Queue 1,Class C)
        routingMatrix.set(jobclass3, jobclass3, node2, node3, 1.000000); // (Queue 1,Class C) -> (ClassSwitch 1,Class C)
        routingMatrix.set(jobclass3, jobclass3, node3, node5, 1.000000); // (ClassSwitch 1,Class C) -> (Queue 2,Class C)
        routingMatrix.set(jobclass3, jobclass3, node5, node4, 1.000000); // (Queue 2,Class C) -> (Sink 1,Class C)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_openModel_6() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue node4 = new Queue(model, "Queue3", SchedStrategy.PS);
        Queue node5 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Sink node6 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);
        OpenClass jobclass3 = new OpenClass(model, "Class3", 0);

        node1.setArrival(jobclass1, Exp.fitMean(5.000000)); // (Source,Class1)
        node1.setArrival(jobclass2, Exp.fitMean(8.000000)); // (Source,Class2)
        node1.setArrival(jobclass3, Exp.fitMean(7.000000)); // (Source,Class3)
        node2.setService(jobclass1, Exp.fitMean(0.300000)); // (Queue1,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.500000)); // (Queue1,Class2)
        node2.setService(jobclass3, Exp.fitMean(0.600000)); // (Queue1,Class3)
        node3.setService(jobclass1, Exp.fitMean(1.100000)); // (Queue2,Class1)
        node3.setService(jobclass2, Exp.fitMean(1.300000)); // (Queue2,Class2)
        node3.setService(jobclass3, Exp.fitMean(1.500000)); // (Queue2,Class3)
        node4.setService(jobclass1, Exp.fitMean(2.000000)); // (Queue3,Class1)
        node4.setService(jobclass2, Exp.fitMean(2.100000)); // (Queue3,Class2)
        node4.setService(jobclass3, Exp.fitMean(1.900000)); // (Queue3,Class3)
        node5.setService(jobclass1, Exp.fitMean(1.500000)); // (Queue4,Class1)
        node5.setService(jobclass2, Exp.fitMean(0.900000)); // (Queue4,Class2)
        node5.setService(jobclass3, Exp.fitMean(2.300000)); // (Queue4,Class3)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Source,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 0.250000); // (Queue1,Class1) -> (Queue2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node4, 0.250000); // (Queue1,Class1) -> (Queue3,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node5, 0.250000); // (Queue1,Class1) -> (Queue4,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node6, 0.250000); // (Queue1,Class1) -> (Sink,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node2, 1.000000); // (Queue2,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node2, 1.000000); // (Queue3,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node5, node2, 1.000000); // (Queue4,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Source,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 0.250000); // (Queue1,Class2) -> (Queue2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node4, 0.250000); // (Queue1,Class2) -> (Queue3,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node5, 0.250000); // (Queue1,Class2) -> (Queue4,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node6, 0.250000); // (Queue1,Class2) -> (Sink,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node2, 1.000000); // (Queue2,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node2, 1.000000); // (Queue3,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node5, node2, 1.000000); // (Queue4,Class2) -> (Queue1,Class2)
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.000000); // (Source,Class3) -> (Queue1,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node3, 0.250000); // (Queue1,Class3) -> (Queue2,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node4, 0.250000); // (Queue1,Class3) -> (Queue3,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node5, 0.250000); // (Queue1,Class3) -> (Queue4,Class3)
        routingMatrix.set(jobclass3, jobclass3, node2, node6, 0.250000); // (Queue1,Class3) -> (Sink,Class3)
        routingMatrix.set(jobclass3, jobclass3, node3, node2, 1.000000); // (Queue2,Class3) -> (Queue1,Class3)
        routingMatrix.set(jobclass3, jobclass3, node4, node2, 1.000000); // (Queue3,Class3) -> (Queue1,Class3)
        routingMatrix.set(jobclass3, jobclass3, node5, node2, 1.000000); // (Queue4,Class3) -> (Queue1,Class3)

        model.link(routingMatrix);

        return model;
    }
}
