// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.io;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Coxian;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.fluid.SolverFluid;

import java.util.Arrays;
import java.util.Collections;

public class M2MTestFixtures {

    public static void main(String[] args) {

        Network model = closed_ex1();

        SolverFluid solver = new SolverFluid(model);
        solver.options.stiff = false;
        solver.getAvgTable();
    }

    public static Network closed_ex1() {

        // Closed Example 1 - 2 Stations
        Network model = new Network("closed_ex1");

        // Block 1: nodes
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 8, delay, 0);
        delay.setService(class1, new Exp(4));
        queue.setService(class1, new Exp(8));

        // Block 3: topology
        model.link(Network.serialRouting(delay, queue));

        return model;
    }

    public static Network closed_ex2() {

        // Closed Example 2 - 4 Stations
        Network model = new Network("closed_ex2");

        // Block 1: nodes
        Delay delay1 = new Delay(model, "Delay1");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Delay delay2 = new Delay(model, "Delay2");
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 16, delay1, 0);
        delay1.setService(class1, new Exp(4));
        queue1.setService(class1, new Exp(8));
        delay2.setService(class1, new Exp(4));
        queue2.setService(class1, new Exp(8));

        // Block 3: topology
        model.link(Network.serialRouting(delay1, queue1, delay2, queue2));

        return model;
    }

    public static Network closed_ex3() {

        // Closed Example 3 - 8 Stations
        Network model = new Network("closed_ex3");

        // Block 1: nodes
        Delay delay1 = new Delay(model, "Delay1");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Delay delay2 = new Delay(model, "Delay2");
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Delay delay3 = new Delay(model, "Delay3");
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Delay delay4 = new Delay(model, "Delay4");
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 32, delay1, 0);
        delay1.setService(class1, new Exp(4));
        queue1.setService(class1, new Exp(8));
        delay2.setService(class1, new Exp(4));
        queue2.setService(class1, new Exp(8));
        delay3.setService(class1, new Exp(4));
        queue3.setService(class1, new Exp(8));
        delay4.setService(class1, new Exp(4));
        queue4.setService(class1, new Exp(8));

        // Block 3: topology
        model.link(Network.serialRouting(delay1, queue1, delay2, queue2, delay3, queue3, delay4, queue4));

        return model;
    }

    public static Network closed_ex4() {

        // Closed Example 4 - 2 JobClasses
        Network model = new Network("closed_ex4");

        // Block 1: nodes
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 8, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 8, delay, 0);
        delay.setService(class1, new Exp(8));
        delay.setService(class2, new Exp(7));
        queue.setService(class1, new Exp(4));
        queue.setService(class2, new Exp(3.5));

        // Block 3: topology
        model.link(Network.serialRouting(delay, queue));

        return model;
    }

    public static Network closed_ex5() {

        // Closed Example 5 - 4 JobClasses
        Network model = new Network("closed_ex5");

        // Block 1: nodes
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 4, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 4, delay, 0);
        ClosedClass class3 = new ClosedClass(model, "Class3", 4, delay, 0);
        ClosedClass class4 = new ClosedClass(model, "Class4", 4, delay, 0);
        delay.setService(class1, new Exp(8));
        delay.setService(class2, new Exp(7));
        delay.setService(class3, new Exp(6));
        delay.setService(class4, new Exp(5));
        queue.setService(class1, new Exp(4));
        queue.setService(class2, new Exp(3.5));
        queue.setService(class3, new Exp(3));
        queue.setService(class4, new Exp(2.5));

        // Block 3: topology
        model.link(Network.serialRouting(delay, queue));

        return model;
    }

    public static Network closed_ex6() {

        // Closed Example 6 - 8 JobClasses
        Network model = new Network("closed_ex6");

        // Block 1: nodes
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 2, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 2, delay, 0);
        ClosedClass class3 = new ClosedClass(model, "Class3", 2, delay, 0);
        ClosedClass class4 = new ClosedClass(model, "Class4", 2, delay, 0);
        ClosedClass class5 = new ClosedClass(model, "Class5", 2, delay, 0);
        ClosedClass class6 = new ClosedClass(model, "Class6", 2, delay, 0);
        ClosedClass class7 = new ClosedClass(model, "Class7", 2, delay, 0);
        ClosedClass class8 = new ClosedClass(model, "Class8", 2, delay, 0);
        delay.setService(class1, new Exp(8));
        delay.setService(class2, new Exp(7));
        delay.setService(class3, new Exp(6));
        delay.setService(class4, new Exp(5));
        delay.setService(class5, new Exp(4));
        delay.setService(class6, new Exp(3));
        delay.setService(class7, new Exp(2));
        delay.setService(class8, new Exp(1));
        queue.setService(class1, new Exp(4));
        queue.setService(class2, new Exp(3.5));
        queue.setService(class3, new Exp(3));
        queue.setService(class4, new Exp(2.5));
        queue.setService(class5, new Exp(2));
        queue.setService(class6, new Exp(1.5));
        queue.setService(class7, new Exp(1));
        queue.setService(class8, new Exp(0.5));

        // Block 3: topology
        model.link(Network.serialRouting(delay, queue));

        return model;
    }

    public static Network closed_ex7() {

        // Closed Example 7 - Balanced Utilisation
        Network model = new Network("closed_ex7");

        // Block 1: nodes
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 16, queue1, 0);
        queue1.setService(class1, new Exp(16));
        queue2.setService(class1, new Exp(16));
        queue3.setService(class1, new Exp(16));
        queue4.setService(class1, new Exp(16));

        // Block 3: topology
        model.link(Network.serialRouting(queue1, queue2, queue3, queue4));

        return model;
    }

    public static Network closed_ex8() {

        // Closed Example 8 - Imbalanced Utilisation
        Network model = new Network("closed_ex8");

        // Block 1: nodes
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.FCFS);

        // Block 2: classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 16, queue1, 0);
        queue1.setService(class1, new Exp(16));
        queue2.setService(class1, new Exp(1));
        queue3.setService(class1, new Exp(32));
        queue4.setService(class1, new Exp(16));

        // Block 3: topology
        model.link(Network.serialRouting(queue1, queue2, queue3, queue4));

        return model;
    }

    public static Network closed_ex9() {

        // Closed Example 9 - Complex
        Network model = new Network("closed_ex9");

        // Block 1: nodes
        Delay delay1 = new Delay(model, "Delay1");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Delay delay2 = new Delay(model, "Delay2");
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Delay delay3 = new Delay(model, "Delay3");
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Delay delay4 = new Delay(model, "Delay4");
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 2, delay1, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 2, delay1, 0);
        ClosedClass class3 = new ClosedClass(model, "Class3", 2, delay2, 0);
        ClosedClass class4 = new ClosedClass(model, "Class4", 2, delay2, 0);
        ClosedClass class5 = new ClosedClass(model, "Class5", 2, delay3, 0);
        ClosedClass class6 = new ClosedClass(model, "Class6", 2, delay3, 0);
        ClosedClass class7 = new ClosedClass(model, "Class7", 2, delay4, 0);
        ClosedClass class8 = new ClosedClass(model, "Class8", 2, delay4, 0);

        delay1.setService(class1, new Exp(8));
        delay1.setService(class2, new Exp(7));
        delay1.setService(class3, new Exp(6));
        delay1.setService(class4, new Exp(5));
        delay1.setService(class5, new Exp(4));
        delay1.setService(class6, new Exp(3));
        delay1.setService(class7, new Exp(2));
        delay1.setService(class8, new Exp(1));
        delay2.setService(class1, new Exp(8));
        delay2.setService(class2, new Exp(7));
        delay2.setService(class3, new Exp(6));
        delay2.setService(class4, new Exp(5));
        delay2.setService(class5, new Exp(4));
        delay2.setService(class6, new Exp(3));
        delay2.setService(class7, new Exp(2));
        delay2.setService(class8, new Exp(1));
        delay3.setService(class1, new Exp(8));
        delay3.setService(class2, new Exp(7));
        delay3.setService(class3, new Exp(6));
        delay3.setService(class4, new Exp(5));
        delay3.setService(class5, new Exp(4));
        delay3.setService(class6, new Exp(3));
        delay3.setService(class7, new Exp(2));
        delay3.setService(class8, new Exp(1));
        delay4.setService(class1, new Exp(8));
        delay4.setService(class2, new Exp(7));
        delay4.setService(class3, new Exp(6));
        delay4.setService(class4, new Exp(5));
        delay4.setService(class5, new Exp(4));
        delay4.setService(class6, new Exp(3));
        delay4.setService(class7, new Exp(2));
        delay4.setService(class8, new Exp(1));

        queue1.setService(class1, new Exp(4));
        queue1.setService(class2, new Exp(3.5));
        queue1.setService(class3, new Exp(3));
        queue1.setService(class4, new Exp(2.5));
        queue1.setService(class5, new Exp(2));
        queue1.setService(class6, new Exp(1.5));
        queue1.setService(class7, new Exp(1));
        queue1.setService(class8, new Exp(0.5));
        queue2.setService(class1, new Exp(4));
        queue2.setService(class2, new Exp(3.5));
        queue2.setService(class3, new Exp(3));
        queue2.setService(class4, new Exp(2.5));
        queue2.setService(class5, new Exp(2));
        queue2.setService(class6, new Exp(1.5));
        queue2.setService(class7, new Exp(1));
        queue2.setService(class8, new Exp(0.5));
        queue3.setService(class1, new Exp(4));
        queue3.setService(class2, new Exp(3.5));
        queue3.setService(class3, new Exp(3));
        queue3.setService(class4, new Exp(2.5));
        queue3.setService(class5, new Exp(2));
        queue3.setService(class6, new Exp(1.5));
        queue3.setService(class7, new Exp(1));
        queue3.setService(class8, new Exp(0.5));
        queue4.setService(class1, new Exp(4));
        queue4.setService(class2, new Exp(3.5));
        queue4.setService(class3, new Exp(3));
        queue4.setService(class4, new Exp(2.5));
        queue4.setService(class5, new Exp(2));
        queue4.setService(class6, new Exp(1.5));
        queue4.setService(class7, new Exp(1));
        queue4.setService(class8, new Exp(0.5));

        // Block 3: topology
        model.link(Network.serialRouting(delay1, queue1, delay2, queue2, delay3, queue3, delay4, queue4));

        return model;
    }

    public static Network open_ex1() {

        // Open Example 1 - 2 Stations
        Network model = new Network("open_ex1");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        source.setArrival(class1, new Exp(2));
        queue1.setService(class1, new Exp(8));
        queue2.setService(class1, new Exp(16));

        // Block 3: topology
        model.link(Network.serialRouting(source, queue1, queue2, sink));

        return model;
    }

    public static Network open_ex2() {

        // Open Example 2 - 4 Stations
        Network model = new Network("open_ex2");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        source.setArrival(class1, new Exp(2));
        queue1.setService(class1, new Exp(8));
        queue2.setService(class1, new Exp(16));
        queue3.setService(class1, new Exp(8));
        queue4.setService(class1, new Exp(16));

        // Block 3: topology
        model.link(Network.serialRouting(source, queue1, queue2, queue3, queue4, sink));

        return model;
    }

    public static Network open_ex3() {

        // Open Example 3 - 8 Stations
        Network model = new Network("open_ex3");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.PS);
        Queue queue5 = new Queue(model, "Queue5", SchedStrategy.FCFS);
        Queue queue6 = new Queue(model, "Queue6", SchedStrategy.PS);
        Queue queue7 = new Queue(model, "Queue7", SchedStrategy.FCFS);
        Queue queue8 = new Queue(model, "Queue8", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        source.setArrival(class1, new Exp(2));
        queue1.setService(class1, new Exp(8));
        queue2.setService(class1, new Exp(16));
        queue3.setService(class1, new Exp(8));
        queue4.setService(class1, new Exp(16));
        queue5.setService(class1, new Exp(8));
        queue6.setService(class1, new Exp(16));
        queue7.setService(class1, new Exp(8));
        queue8.setService(class1, new Exp(16));

        // Block 3: topology
        model.link(
                Network.serialRouting(
                        source, queue1, queue2, queue3, queue4, queue5, queue6, queue7, queue8, sink));

        return model;
    }

    public static Network open_ex4() {

        // Open Example 4 - 2 JobClasses
        Network model = new Network("open_ex4");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        OpenClass class2 = new OpenClass(model, "Class2");
        source.setArrival(class1, new Exp(2));
        source.setArrival(class2, new Exp(2));
        queue.setService(class1, new Exp(24));
        queue.setService(class2, new Exp(20));

        // Block 3: topology
        model.link(Network.serialRouting(source, queue, sink));

        return model;
    }

    public static Network open_ex5() {

        // Open Example 5 - 4 JobClasses
        Network model = new Network("open_ex5");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        OpenClass class2 = new OpenClass(model, "Class2");
        OpenClass class3 = new OpenClass(model, "Class3");
        OpenClass class4 = new OpenClass(model, "Class4");
        source.setArrival(class1, new Exp(2));
        source.setArrival(class2, new Exp(2));
        source.setArrival(class3, new Exp(2));
        source.setArrival(class4, new Exp(2));
        queue.setService(class1, new Exp(24));
        queue.setService(class2, new Exp(20));
        queue.setService(class3, new Exp(24));
        queue.setService(class4, new Exp(20));

        // Block 3: topology
        model.link(Network.serialRouting(source, queue, sink));

        return model;
    }

    public static Network open_ex6() {

        // Open Example 6 - 8 JobClasses
        Network model = new Network("open_ex6");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        OpenClass class2 = new OpenClass(model, "Class2");
        OpenClass class3 = new OpenClass(model, "Class3");
        OpenClass class4 = new OpenClass(model, "Class4");
        OpenClass class5 = new OpenClass(model, "Class5");
        OpenClass class6 = new OpenClass(model, "Class6");
        OpenClass class7 = new OpenClass(model, "Class7");
        OpenClass class8 = new OpenClass(model, "Class8");
        source.setArrival(class1, new Exp(2));
        source.setArrival(class2, new Exp(2));
        source.setArrival(class3, new Exp(2));
        source.setArrival(class4, new Exp(2));
        source.setArrival(class5, new Exp(2));
        source.setArrival(class6, new Exp(2));
        source.setArrival(class7, new Exp(2));
        source.setArrival(class8, new Exp(2));
        queue.setService(class1, new Exp(24));
        queue.setService(class2, new Exp(20));
        queue.setService(class3, new Exp(24));
        queue.setService(class4, new Exp(20));
        queue.setService(class5, new Exp(24));
        queue.setService(class6, new Exp(20));
        queue.setService(class7, new Exp(24));
        queue.setService(class8, new Exp(20));

        // Block 3: topology
        model.link(Network.serialRouting(source, queue, sink));

        return model;
    }

    public static Network open_ex7() {

        // Open Example 7 - Low Utilisation
        Network model = new Network("open_ex7");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        source.setArrival(class1, new Exp(2));
        queue1.setService(class1, new Exp(64));

        // Block 3: topology
        model.link(Network.serialRouting(source, queue1, sink));

        return model;
    }

    public static Network open_ex8() {

        // Open Example 8 - High Utilisation
        Network model = new Network("open_ex8");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        source.setArrival(class1, new Exp(2));
        queue1.setService(class1, new Exp(2.05));

        // Block 3: topology
        model.link(Network.serialRouting(source, queue1, sink));

        return model;
    }

    public static Network open_ex9() {

        // Open Example 9 - Mid Utilisation
        Network model = new Network("open_ex9");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        source.setArrival(class1, new Exp(2));
        queue1.setService(class1, new Exp(4));

        // Block 3: topology
        model.link(Network.serialRouting(source, queue1, sink));

        return model;
    }

    public static Network open_ex10() {

        // Open Example 10 - Complex
        Network model = new Network("open_ex10");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        OpenClass class2 = new OpenClass(model, "Class2");
        OpenClass class3 = new OpenClass(model, "Class3");
        OpenClass class4 = new OpenClass(model, "Class4");
        source.setArrival(class1, new Exp(2));
        source.setArrival(class2, new Exp(2));
        source.setArrival(class3, new Exp(2));
        source.setArrival(class4, new Exp(2));
        queue1.setService(class1, new Exp(32));
        queue1.setService(class2, new Exp(32));
        queue1.setService(class3, new Exp(32));
        queue1.setService(class4, new Exp(32));
        queue2.setService(class1, new Exp(32));
        queue2.setService(class2, new Exp(32));
        queue2.setService(class3, new Exp(32));
        queue2.setService(class4, new Exp(32));
        queue3.setService(class1, new Exp(32));
        queue3.setService(class2, new Exp(32));
        queue3.setService(class3, new Exp(32));
        queue3.setService(class4, new Exp(32));
        queue4.setService(class1, new Exp(32));
        queue4.setService(class2, new Exp(32));
        queue4.setService(class3, new Exp(32));
        queue4.setService(class4, new Exp(32));

        // Block 3: topology
        model.link(Network.serialRouting(source, queue1, queue2, queue3, queue4, sink));

        return model;
    }

    public static Network other_ex1() {

        // Example 1 - 4 Queues in Parallel, different scheduling algorithms, phase-type service
        // distributions
        Network model = new Network("ex1");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Queue queue4 = new Queue(model, "Queue4", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        source.setArrival(class1, new Exp(2));
        queue1.setService(class1, new Exp(4));
        queue2.setService(class1, new Erlang(4, 2));
        queue3.setService(class1, new Exp(8));
        queue4.setService(class1, new Erlang(6, 2));

        // Block 3: topology
        RoutingMatrix routingMatrix =
                new RoutingMatrix(
                        model,
                        Collections.singletonList(class1),
                        Arrays.asList(source, queue1, queue2, queue3, queue4, sink));
        routingMatrix.set(class1, source, queue1, 0.3);
        routingMatrix.set(class1, source, queue2, 0.2);
        routingMatrix.set(class1, source, queue3, 0.3);
        routingMatrix.set(class1, source, queue4, 0.2);
        routingMatrix.set(class1, queue1, sink, 1);
        routingMatrix.set(class1, queue2, sink, 1);
        routingMatrix.set(class1, queue3, sink, 1);
        routingMatrix.set(class1, queue4, sink, 1);

        model.link(routingMatrix);

        return model;
    }

    public static Network other_ex2() {

        // Example 2
        Network model = new Network("ex2");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass class1 = new OpenClass(model, "Class1");
        OpenClass class2 = new OpenClass(model, "Class2");
        OpenClass class3 = new OpenClass(model, "Class3");
        source.setArrival(class1, new Exp(2));
        source.setArrival(class2, new Exp(2));
        source.setArrival(class3, new Exp(1));
        queue1.setService(class1, new Exp(32));
        queue1.setService(class2, new Exp(32));
        queue1.setService(class3, new Exp(32));
        queue2.setService(class1, new Exp(32));
        queue2.setService(class2, new Exp(32));
        queue2.setService(class3, new Exp(16));
        queue3.setService(class1, new Exp(16));
        queue3.setService(class2, new Exp(16));
        queue3.setService(class3, new Exp(8));

        // Block 3: topology
        RoutingMatrix routingMatrix =
                new RoutingMatrix(
                        model,
                        Arrays.asList(class1, class2, class3),
                        Arrays.asList(source, queue1, queue2, queue3, sink));
        routingMatrix.set(class1, class1, source, queue1, 1);
        routingMatrix.set(class1, class1, queue1, queue2, 1);
        routingMatrix.set(class1, class3, queue2, queue1, 0.2);
        routingMatrix.set(class1, class1, queue2, queue3, 0.8);
        routingMatrix.set(class1, class1, queue3, sink, 1);

        routingMatrix.set(class2, class2, source, queue1, 1);
        routingMatrix.set(class2, class2, queue1, queue2, 1);
        routingMatrix.set(class2, class3, queue2, queue1, 0.2);
        routingMatrix.set(class2, class2, queue2, queue3, 0.8);
        routingMatrix.set(class2, class2, queue3, sink, 1);

        routingMatrix.set(class3, class3, source, queue1, 1);
        routingMatrix.set(class3, class3, queue1, queue2, 1);
        routingMatrix.set(class3, class3, queue2, queue1, 0.2);
        routingMatrix.set(class3, class3, queue2, queue3, 0.8);
        routingMatrix.set(class3, class3, queue3, sink, 1);

        model.link(routingMatrix);

        return model;
    }

    public static Network ruuskanenExample2() {

        Network model = new Network("RuuskanenExample2");

        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        queue1.setNumberOfServers(4);
        queue2.setNumberOfServers(8);

        ClosedClass class1 = new ClosedClass(model, "Job1", 50, delay, 0);
        delay.setService(class1, new Exp(0.6));
        queue1.setService(class1, Coxian.fitMeanAndSCV(0.5, 0.5));
        queue2.setService(class1, Coxian.fitMeanAndSCV(1, 10));

        model.link(Network.serialRouting(delay, queue1, queue2));

        return model;
    }
}
