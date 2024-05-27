package jline.examples;

import jline.lang.nodes.*;
import jline.solvers.ctmc.SolverCTMC;
import jline.lang.nodes.Delay;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.*;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mva.SolverMVA;
import jline.util.Matrix;

import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;

import javax.xml.parsers.ParserConfigurationException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Getting started examples
 */
public class GettingStarted {

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

    public static Network ex1_line() {
        Network model = new Network("M/M/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass openClass = new OpenClass(model, "myClass");
        source.setArrival(openClass, new Exp(1));
        queue.setService(openClass, new Exp(2));

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


        RoutingMatrix routingMatrix = new RoutingMatrix(model, Collections.singletonList(openClass),
                Arrays.asList(source, queue1, queue2, queue3, sink));
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

    public static Network ex5_line() {
        /* A closed network of 3 queues
         */
        Network model = new Network("3 Closed");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

        ClosedClass jobClass1 = new ClosedClass(model, "Class1", 1, queue);
        ClosedClass jobClass2 = new ClosedClass(model, "Class2", 0, queue);
        ClosedClass jobClass3 = new ClosedClass(model, "Class3", 0, queue);

        queue.setService(jobClass1, Erlang.fitMeanAndOrder(1, 2));
        queue.setService(jobClass2, Erlang.fitMeanAndOrder(2, 2));
        queue.setService(jobClass3, Erlang.fitMeanAndOrder(3, 2));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass1, jobClass2, queue, queue, 1.0);
        P.set(jobClass2, jobClass3, queue, queue, 1.0);
        P.set(jobClass3, jobClass1, queue, queue, 1.0);


        //jobClass1.setCompletes(false);
        //jobClass2.setCompletes(false);
        model.link(P);
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

    public static Network matlabExample3() {
        Network model = new Network("MRP");
        Delay delay = new Delay(model, "Working State");
        Queue queue = new Queue(model, "RepairQueue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);

        ClosedClass closedClass = new ClosedClass(model, "Machines", 3, delay);
        delay.setService(closedClass, new Exp(0.5));
        queue.setService(closedClass, new Exp(4.0));
        model.link(model.serialRouting(delay, queue));
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

        RoutingMatrix routingMatrix = new RoutingMatrix(model, Collections.singletonList(closedClass1),
                Arrays.asList(Node1, Node2));
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

        RoutingMatrix routingMatrix = new RoutingMatrix(model, Collections.singletonList(openClass),
                Arrays.asList(source, fork, upper, lower, join, post, sink));
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


        RoutingMatrix routingMatrix = new RoutingMatrix(model, Arrays.asList(openClass, oClass2),
                Arrays.asList(source, queue1, queue2, queue3, sink));
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
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(closedClass1, closedClass2),
                Arrays.asList(delay, queue));

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
        RoutingMatrix P = new RoutingMatrix(model,
                Collections.singletonList(jobclass1),
                Arrays.asList(node1, node2, node3, node4, node5, node6, node7, node8, node9, node10, node11, node12,
                        switchNodes.get(0), switchNodes.get(1), switchNodes.get(2), switchNodes.get(3), switchNodes.get(4), switchNodes.get(5), switchNodes.get(6), switchNodes.get(7), switchNodes.get(8), switchNodes.get(9)));
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
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(closedClass, closedClass2),
                Arrays.asList(delay, queue));
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
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4));
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


    public static void main(String[] args) throws IllegalAccessException, ParserConfigurationException {
        Network model = ex1_line();
        //SolverCTMC solver = new SolverCTMC(model);
        //SolverJMT solver = new SolverJMT(model);
        SolverMVA solver = new SolverMVA(model);
        //solver.applyCutoff(3);
        solver.getAvgTable().print();
        solver.getAvgSysTable().print();
        //solver.jsimgView();
//        solver.getGenerator();
//        solver.getStateSpace();
//        solver.getProbabilityVector();
    }
}
