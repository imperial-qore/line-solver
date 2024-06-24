package jline.examples;

import jline.lang.constant.*;
import jline.lang.nodes.*;
import jline.lang.processes.Replayer;
import jline.solvers.*;
import jline.solvers.ctmc.SolverCTMC;
import jline.lang.nodes.Delay;
import jline.lang.*;
import jline.lang.distributions.*;
import jline.solvers.fluid.SolverFluid;
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
import jline.solvers.ssa.SolverSSA;


import java.util.function.Function;

import static jline.solvers.jmt.JMTOptions.*;

/**
 * Getting started examples
 */
public class GettingStarted {

    public static Network getting_started_1() {
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

        return model;
    }

    public static Network getting_started_2() {
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
        P.set(jobclass1, Network.serialRouting(source,queue,sink));
        P.set(jobclass2, Network.serialRouting(source,queue,sink));
        model.link(P);
        NetworkAvgTable avgTable = new SolverJMT(model, "seed", 23000).getAvgTable();
        avgTable.print();

        return model;
    }

    public static Network getting_started_3() {
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
        return model;
    }

    public static Network getting_started_4() {
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
        return model;
    }

    public static Network getting_started_5() {
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
        return model;
    }

    public static Network getting_started_6() {
            Network model = new Network("Model");

            // Block 1: nodes
            Delay clientDelay = new Delay(model, "Client");
            Cache cacheNode = new Cache(model, "Cache", 1000, 50, ReplacementStrategy.LRU);
            Delay cacheDelay = new Delay(model, "CacheDelay");

            // Block 2: classes
            ClosedClass clientClass = new ClosedClass(model, "ClientClass", 1, clientDelay, 0);
            ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, clientDelay, 0);
            ClosedClass missClass = new ClosedClass(model, "MissClass", 0, clientDelay, 0);

            clientDelay.setService(clientClass, new Immediate()); // (Client,ClientClass)
            cacheDelay.setService(hitClass, Exp.fitMean(0.200000)); // (CacheDelay,HitClass)
            cacheDelay.setService(missClass, Exp.fitMean(1.000000)); // (CacheDelay,MissClass)

            cacheNode.setRead(clientClass, new Zipf(1.4,1000));
            cacheNode.setHitClass(clientClass, hitClass);
            cacheNode.setMissClass(clientClass, missClass);

            // Block 3: topology
            RoutingMatrix routingMatrix = model.initRoutingMatrix();

            routingMatrix.set(clientClass, clientClass, clientDelay, cacheNode, 1.000000); // (Client,ClientClass) -> (Cache,ClientClass)

            routingMatrix.set(hitClass, hitClass, cacheNode, cacheDelay, 1.000000); // (Client,HitClass) -> (Cache,HitClass)
            routingMatrix.set(missClass, missClass, cacheNode, cacheDelay, 1.000000); // (Cache,MissClass) -> (CacheDelay,MissClass)

            routingMatrix.set(hitClass, clientClass, cacheDelay, clientDelay, 1.000000); // (Cache,HitClass) -> (CacheDelay,HitClass)
            routingMatrix.set(missClass, clientClass, cacheDelay, clientDelay, 1.000000); // (Client,MissClass) -> (Cache,MissClass)

            model.link(routingMatrix);
            //TODO
            //new SolverSSA(model,"samples",2e4,"seed",1,"verbose",true).getAvgTable().print();
            return model;
        }

    public static Network getting_started_7() {
        Network model = new Network("Model");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 5, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
        node2.setService(jobclass1, Exp.fitMean(2.000000)); // (Queue1,Class1)

        // Block 3: topology
        model.link(Network.serialRouting(node1,node2));

        //TODO
        //RDfluid = new SolverFluid(model).getCdfRespT();
        //RDsim = new SolverJMT(model,"seed",23000,"samples",1e4).getCdfRespT();

        return model;
    }

    public static Network getting_started_8() {
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
        return model;
    }

    public static void getting_started_9() {
        // TODO
    }

    public static void getting_started_10() {
        // TODO
    }

    public static void main(String[] args) throws IllegalAccessException, ParserConfigurationException {
        getting_started_6();
    }
}