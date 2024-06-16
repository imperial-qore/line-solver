package jline.examples;

import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SolverType;
import jline.lang.constant.VerboseLevel;
import jline.lang.nodes.*;
import jline.lang.processes.Replayer;
import jline.solvers.*;
import jline.solvers.ctmc.SolverCTMC;
import jline.lang.nodes.Delay;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.*;
import jline.solvers.jmt.JMTOptions;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mva.SolverMVA;

import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;

import javax.xml.parsers.ParserConfigurationException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Paths;
import de.xypron.jcobyla.Calcfc;
import de.xypron.jcobyla.Cobyla;
import jline.solvers.nc.SolverNC;


import java.util.function.Function;

import static jline.solvers.jmt.JMTOptions.*;

/**
 * Getting started examples
 */
public class GettingStarted {

    public static void getting_started_1() {
        Network model = new Network("M/M/1");
        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        // Block 2: classes
        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, new Exp(1.0)); // (source, jobclass)
        queue.setService(jobclass, new Exp(2.0)); // (queue, jobclass)
        // Block 3: topology
        model.link(model.serialRouting(source, queue, sink));
        // Block 4: solution
        NetworkAvgTable avgTable = new SolverJMT(model, "seed", 23000).getAvgTable();
        avgTable.print();
        //model.jsimgView();
        avgTable.tget(queue, jobclass).print();
        avgTable.tget("Queue", "Class1").print();
    }

    public static void getting_started_2() {
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
            URI fileURI = GettingStarted.class.getResource("/example_trace.txt").toURI();
            String fileName = Paths.get(fileURI).toString();
            queue.setService(jobclass2, new Replayer(fileName));
            // queue.setService(jobclass2, new Replayer(fileName).fitAPH());
        } catch (URISyntaxException e) {
            throw new RuntimeException(e);
        }
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass1, jobclass1, Network.serialRouting(source,queue,sink));
        P.set(jobclass2, jobclass2, Network.serialRouting(source,queue,sink));
        model.link(P);
        NetworkAvgTable avgTable = new SolverJMT(model, "seed", 23000).getAvgTable();
        avgTable.print();
    }

    public static void getting_started_3() {
        Network model = new Network("MRP");
        Delay delay = new Delay(model, "WorkingState");
        Queue queue = new Queue(model, "RepairQueue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);
        ClosedClass closedClass = new ClosedClass(model, "Machines", 3, delay);
        delay.setService(closedClass, new Exp(0.5));
        queue.setService(closedClass, new Exp(4.0));
        model.link(model.serialRouting(delay, queue));
        SolverCTMC solver = new SolverCTMC(model);
        NetworkAvgTable ctmcAvgTable = solver.getAvgTable();
        ctmcAvgTable.print();
//        ArrayList<SSAStateMatrix> stateSpace = solver.getStateSpace();
//        for (SSAStateMatrix matrix : stateSpace) {
//            matrix.print();
//        }
//        Matrix infGen = solver.getGenerator();
//        infGen.print();
    }

    public static void getting_started_4() {
        Network model = new Network("RRLB");
        Source source = new Source(model, "Source");
        Router lb = new Router(model, "LB");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "Class1");
        source.setArrival(jobclass, new Exp(1.0));
        queue1.setService(jobclass, new Exp(2.0));
        queue2.setService(jobclass, new Exp(2.0));
        model.addLinks(new Node[][]{
                {source, lb},
                {lb, queue1},
                {lb, queue2},
                {queue1, sink},
                {queue2, sink}
        });
        new SolverJMT(model, "seed", 23000).getAvgTable().print();
        model.reset();
        lb.setRouting(jobclass, RoutingStrategy.RROBIN);
        new SolverJMT(model, "seed", 23000).getAvgTable().print();
    }

    public static void getting_started_5() {
        Network model = new Network("RL");
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
        model.link(P);
        new SolverNC(model).getAvgTable().print();
        new SolverNC(model).getAvgSysTable().print();
        jobClass1.setCompletes(false);
        jobClass2.setCompletes(false);
        new SolverNC(model).getAvgSysTable().print();
    }


    public static void getting_started_8() {
        Network model = new Network("LoadBalCQN");

        // Block 1: nodes
        Delay delay = new Delay(model, "Think");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass cclass = new ClosedClass(model, "Job1", 16, delay);
        delay.setService(cclass, new Exp(1));
        queue1.setService(cclass, new Exp(0.75));
        queue2.setService(cclass, new Exp(0.50));

        // Block 3: topology
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(cclass, queue1, delay, 1.0);
        P.set(cclass, queue2, delay, 1.0);
        model.link(P);

        Function<Double, Double> routingModel = (p) -> {
            // Block 4: solution
            P.set(cclass, delay, queue1, p);
            P.set(cclass, delay, queue2, 1 - p);
            model.reset();
            model.link(P);
            SolverMVA solver = new SolverMVA(model, "exact");
            return (Double) solver.getAvgSysRespT().get(0);
        };

        Calcfc objFun = new Calcfc() {
            @Override
            public double compute(int n, int m, double[] x, double[] con) {
                return routingModel.apply(x[0]) ;
            }
        };

        double[] p = {0.5};
        Cobyla.findMinimum(objFun, 1, 0, p, 0.5, 1.0e-8, 0, 10000);

        System.out.println("Optimal p: " + p[0]);
    }

    public static void main(String[] args) throws IllegalAccessException, ParserConfigurationException {
        getting_started_1();
    }
}