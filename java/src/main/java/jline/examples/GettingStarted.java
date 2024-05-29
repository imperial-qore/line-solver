package jline.examples;

import jline.lang.constant.RoutingStrategy;
import jline.lang.nodes.*;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.*;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mam.SolverMAM;
import jline.lang.processes.Replayer;

import javax.xml.parsers.ParserConfigurationException;
import java.net.URISyntaxException;
import java.nio.file.Paths;

/**
 * Getting started examples
 */
public class GettingStarted {

    public static Network ex1() {
        Network model = new Network("M/M/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, new Exp(2));

        model.link(model.serialRouting(source, queue, sink));

        return model;
    }

    public static Network ex2() {
        Network model = new Network("M/G/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass1 = new OpenClass(model, "Class1");
        OpenClass jobclass2 = new OpenClass(model, "Class2");
        source.setArrival(jobclass1, new Exp(0.5));
        source.setArrival(jobclass2, new Exp(0.5));
        queue.setService(jobclass1, Erlang.fitMeanAndSCV(1, (double) 1 / 3));
        try {
            String fileName = Paths.get(GettingStarted.class.getResource("/example_trace.txt").toURI()).toString();
            queue.setService(jobclass2, new Replayer(fileName));
        } catch (URISyntaxException e) {
            throw new RuntimeException(e);
        }

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass1, jobclass1, Network.serialRouting(source,queue,sink));
        P.set(jobclass2, jobclass2, Network.serialRouting(source,queue,sink));
        model.link(P);

        return model;
    }

    public static Network ex3() {
        Network model = new Network("MRP");

        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Working State");
        Queue queue = new Queue(model, "RepairQueue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);

        ClosedClass closedClass = new ClosedClass(model, "Machines", 3, delay);
        delay.setService(closedClass, new Exp(0.5));
        queue.setService(closedClass, new Exp(4.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    public static Network ex4() {
        Network model = new Network("RRLB");
        Source source = new Source(model, "Source");
        Router lb = new Router(model, "LB");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass = new OpenClass(model, "Class1");
        source.setArrival(jobclass, new Exp(1));
        queue1.setService(jobclass, new Exp(2));
        queue2.setService(jobclass, new Exp(2));

        model.addLinks(new Node[][]{
                {source, lb},
                {lb, queue1},
                {lb, queue2},
                {queue1, sink},
                {queue2, sink}
        });

        return model;
    }

    public static Network ex5() {
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

    public static void main(String[] args) throws IllegalAccessException, ParserConfigurationException {
        Network model = new Network("RRLB");
        Source source = new Source(model, "Source");
        Router lb = new Router(model, "LB");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobclass = new OpenClass(model, "Class1");
        source.setArrival(jobclass, new Exp(1));
        queue1.setService(jobclass, new Exp(2));
        queue2.setService(jobclass, new Exp(2));

        model.addLinks(new Node[][]{
                {source, lb},
                {lb, queue1},
                {lb, queue2},
                {queue1, sink},
                {queue2, sink}
        });
        new SolverJMT(model).getAvgTable().print();
        model.reset();
        lb.setRouting(jobclass, RoutingStrategy.RROBIN);
        new SolverJMT(model).getAvgTable().print();
    }

}
