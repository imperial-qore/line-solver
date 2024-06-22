package jline.examples;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Exp;
import jline.lang.distributions.Immediate;
import jline.lang.nodes.*;

import java.util.Arrays;
import java.util.Collections;

/**
 * Examples of models with fork-join subsystems
 */
public class ForkJoinModel {

    public static Network example_forkJoin_1(){
        Network model = new Network("model");

        Source source = new Source(model,"Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass1 = new OpenClass(model, "class1");
        source.setArrival(jobclass1, new Exp(0.05));
        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(2.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, source, fork, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1);
        routingMatrix.set(jobclass1, jobclass1, queue2, join, 1);
        routingMatrix.set(jobclass1, jobclass1, join, sink, 1);
        model.link(routingMatrix);

        return model;
    }

    public static Network example_forkJoin_2(){
        Network model = new Network("model");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue2", SchedStrategy.PS);
        Fork node4 = new Fork(model, "Fork");
        node4.setTasksPerLink(2);
        Join node5 = new Join(model, "Join", node4);
        Sink node6 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "class2", 0);

        node1.setArrival(jobclass1, new Exp(0.25)); // (Source,class1)
        node1.setArrival(jobclass2, new Exp(0.25)); // (Source,class2)

        node2.setService(jobclass1, new Exp(1.0)); // (Queue1,class1)
        node2.setService(jobclass2, new Immediate()); // (Queue1,class2)
        node3.setService(jobclass1, new Exp(0.75)); // (Queue2,class1)
        node3.setService(jobclass2, new Exp(2.0)); // (Queue2,class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node4, 1.000000); // (Source,class1) -> (Fork,class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node5, 1.000000); // (Queue1,class1) -> (Join,class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node5, 1.000000); // (Queue2,class1) -> (Join,class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node2, 1.000000); // (Fork,class1) -> (Queue1,class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node3, 1.000000); // (Fork,class1) -> (Queue2,class1)
        routingMatrix.set(jobclass1, jobclass1, node5, node6, 1.000000); // (Join,class1) -> (Sink,class1)

        routingMatrix.set(jobclass2, jobclass2, node1, node4, 1.000000); // (Source,class2) -> (Fork,class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node5, 1.000000); // (Queue1,class2) -> (Join,class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node5, 1.000000); // (Queue2,class2) -> (Join,class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node2, 1.000000); // (Fork,class2) -> (Queue1,class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node3, 1.000000); // (Fork,class2) -> (Queue2,class2)
        routingMatrix.set(jobclass2, jobclass2, node5, node6, 1.000000); // (Join,class2) -> (Sink,class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_forkJoin_3() {
        Network model = new Network("model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Fork node2 = new Fork(model, "Fork1");
        node2.setTasksPerLink(1);
        Fork node3 = new Fork(model, "Fork1_1");
        node3.setTasksPerLink(2);
        Join node4 = new Join(model, "Join1", node2);
        Join node5 = new Join(model, "Join1_1", node3);
        Queue node6 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node7 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "class1", 5, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "class2", 2, node1, 0);

        node1.setService(jobclass1, new Exp(0.25)); // (Delay,class1)
        node1.setService(jobclass2, new Exp(0.25)); // (Delay,class2)
        node6.setService(jobclass1, new Exp(1.0)); // (Queue1,class1)
        node6.setService(jobclass2, new Exp(2.0)); // (Queue1,class2)
        node7.setService(jobclass1, new Exp(0.75)); // (Queue2,class1)
        node7.setService(jobclass2, new Exp(2.0)); // (Queue2,class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4, node5, node6, node7));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Delay,class1) -> (Fork1,class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node6, 1.000000); // (Fork1,class1) -> (Queue1,class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node7, 1.000000); // (Fork1,class1) -> (Queue2,class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node1, 1.000000); // (Join1,class1) -> (Delay,class1)
        routingMatrix.set(jobclass1, jobclass1, node6, node4, 1.000000); // (Queue1,class1) -> (Join1,class1)
        routingMatrix.set(jobclass1, jobclass1, node7, node4, 1.000000); // (Queue2,class1) -> (Join1,class1)

        routingMatrix.set(jobclass2, jobclass2, node1, node3, 1.000000); // (Delay,class2) -> (Fork1_1,class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node6, 1.000000); // (Fork1,class2) -> (Queue1,class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node7, 1.000000); // (Fork1,class2) -> (Queue2,class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node2, 1.000000); // (Fork1_1,class2) -> (Fork1,class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node5, 1.000000); // (Join1,class2) -> (Join1_1,class2)
        routingMatrix.set(jobclass2, jobclass2, node5, node1, 1.000000); // (Join1_1,class2) -> (Delay,class2)
        routingMatrix.set(jobclass2, jobclass2, node6, node4, 1.000000); // (Queue1,class2) -> (Join1,class2)
        routingMatrix.set(jobclass2, jobclass2, node7, node4, 1.000000); // (Queue2,class2) -> (Join1,class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network example_forkJoin_4() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network example_forkJoin_5(){
        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 5 ,delay);

        delay.setService(jobclass1, new Exp(1.0));
        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(1.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, delay, fork, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1);
        routingMatrix.set(jobclass1, jobclass1, queue2, join, 1);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1);
        model.link(routingMatrix);

        return model;
    }

    public static Network example_forkJoin_6() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network example_forkJoin_7() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network example_forkJoin_8() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network example_forkJoin_9(){
        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 10,delay, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "class2", 10,delay, 0);

        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(2.0));
        delay.setService(jobclass1, new Exp(0.5));

        queue1.setService(jobclass2, new Exp(1.0));
        queue2.setService(jobclass2, new Exp(2.0));
        delay.setService(jobclass2, new Exp(0.2));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, delay, fork, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1);
        routingMatrix.set(jobclass1, jobclass1, queue2, join, 1);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1);

        routingMatrix.set(jobclass2, jobclass2, delay, fork, 1);
        routingMatrix.set(jobclass2, jobclass2, fork, queue1, 1);
        routingMatrix.set(jobclass2, jobclass2, fork, queue2, 1);
        routingMatrix.set(jobclass2, jobclass2, queue1, join, 1);
        routingMatrix.set(jobclass2, jobclass2, queue2, join, 1);
        routingMatrix.set(jobclass2, jobclass2, join, delay, 1);
        model.link(routingMatrix);

        return model;
    }

    public static Network example_forkJoin_10(){
        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 10,delay, 0);

        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(2.0));
        queue3.setService(jobclass1, new Exp(1.0));
        delay.setService(jobclass1, new Exp(0.5));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, delay, fork, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1);
        routingMatrix.set(jobclass1, jobclass1, queue1, join, 1);
        routingMatrix.set(jobclass1, jobclass1, queue2, queue3, 1);
        routingMatrix.set(jobclass1, jobclass1, queue3, join, 1);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1);

        model.link(routingMatrix);

        return model;
    }

    public static Network example_forkJoin_11() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network example_forkJoin_12() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network example_forkJoin_13() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network example_forkJoin_cs() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network example_forkJoin_nested(){
        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay1");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Fork fork2 = new Fork(model, "Fork2");
        Join join2 = new Join(model, "Join2", fork2);

        ClosedClass jobclass1 = new ClosedClass(model, "class1", 1,delay, 0);

        queue1.setService(jobclass1, new Exp(1.0));
        queue2.setService(jobclass1, new Exp(1.0));
        delay.setService(jobclass1, new Exp(0.5));
        queue3.setService(jobclass1, new Exp(2.0));
        queue4.setService(jobclass1, new Exp(2.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, delay, fork, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue1, 1);
        routingMatrix.set(jobclass1, jobclass1, fork, queue2, 1);
        routingMatrix.set(jobclass1, jobclass1, queue1, fork2, 1);
        routingMatrix.set(jobclass1, jobclass1, fork2, queue3, 1);
        routingMatrix.set(jobclass1, jobclass1, fork2, queue4, 1);
        routingMatrix.set(jobclass1, jobclass1, queue3, join2, 1);
        routingMatrix.set(jobclass1, jobclass1, queue4, join2, 1);
        routingMatrix.set(jobclass1, jobclass1, join2, join, 1);
        routingMatrix.set(jobclass1, jobclass1, queue2, join, 1);
        routingMatrix.set(jobclass1, jobclass1, join, delay, 1);

        model.link(routingMatrix);

        return model;
    }

    public static Network example_forkJoin_series() {
        Network model = new Network("model");
        // TODO
        return model;
    }
}
