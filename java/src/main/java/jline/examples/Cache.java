package jline.examples;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.Disabled;
import jline.lang.distributions.DiscreteSampler;
import jline.lang.distributions.Exp;
import jline.lang.nodes.*;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.Matrix;
import jline.solvers.NetworkAvgTable;

import java.util.Arrays;

/**
 * Examples of caching models
 */
public class Cache {
    // TODO: missing NC methods for caches

    public static Network ex1(){
        Network model = new Network("model");

        int n = 5; // Number of items
        Matrix m = new Matrix(1,1); // cache capacity
        m.set(0, 2);

        Source source = new Source(model, "Source");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.FIFO, null);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(2));

        Matrix p = new Matrix(1, n);
        p.fill(1.0/n);
        DiscreteSampler pAccess = new DiscreteSampler(p);

        cacheNode.setRead(jobClass, pAccess);

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobClass, hitClass, missClass),
                Arrays.asList(source, cacheNode, sink));

        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1);

        model.link(routingMatrix);

        return model;
    }

    public static Network ex2(){
        Network model = new Network("model");

        int n = 5; // Number of items
        Matrix m = new Matrix(1,1); // cache capacity
        m.set(0, 2);

        Delay delay = new Delay(model, "Delay");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.FIFO, null);


        ClosedClass jobClass = new ClosedClass(model, "JobClass", 1, delay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, delay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, delay, 0);

        delay.setService(jobClass, new Exp(1));

        Matrix p = new Matrix(1, n);
        p.fill(1.0/n);
        DiscreteSampler pAccess = new DiscreteSampler(p);

        cacheNode.setRead(jobClass, pAccess);

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobClass, hitClass, missClass),
                Arrays.asList(delay, cacheNode));

        routingMatrix.set(jobClass, jobClass, delay, cacheNode, 1);
        routingMatrix.set(hitClass, jobClass, cacheNode, delay, 1);
        routingMatrix.set(missClass, jobClass, cacheNode, delay, 1);

        model.link(routingMatrix);

        return model;
    }

    public static Network ex5(){
        Network model = new Network("model");

        int n = 5; // Number of items
        Matrix m = new Matrix(1,1); // cache capacity
        m.set(0, 2);

        Source source = new Source(model, "Source");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.FIFO, null);
        Router routerNode = new Router(model, "Router");
        Delay delay1 = new Delay(model, "Delay1");
        Delay delay2 = new Delay(model, "Delay2");
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(2));
        source.setArrival(hitClass, new Disabled());
        source.setArrival(missClass, new Disabled());

        delay1.setService(hitClass, new Exp(10));
        delay1.setService(missClass, new Exp(1));
        delay2.setService(hitClass, new Exp(20));
        delay2.setService(missClass, new Exp(2));


        Matrix p = new Matrix(1, n);
        p.fill(1.0/n);
        DiscreteSampler pAccess = new DiscreteSampler(p);

        cacheNode.setRead(jobClass, pAccess);

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        model.addLink(source, cacheNode);
        model.addLink(cacheNode, routerNode);
        model.addLink(routerNode, delay1);
        model.addLink(routerNode, delay2);
        model.addLink(delay1, sink);
        model.addLink(delay2, sink);

        source.setProbRouting(jobClass, cacheNode, 1);

        cacheNode.setRouting(jobClass, RoutingStrategy.RAND);
        cacheNode.setProbRouting(hitClass, routerNode, 1);
        cacheNode.setProbRouting(missClass, routerNode, 1);

        routerNode.setRouting(hitClass, RoutingStrategy.RAND);
        routerNode.setRouting(missClass, RoutingStrategy.RAND);

        delay1.setProbRouting(hitClass, sink, 1);
        delay1.setProbRouting(missClass, sink, 1);

        delay2.setProbRouting(hitClass, sink, 1);
        delay2.setProbRouting(missClass, sink, 1);

        return model;
    }

    public static void main(String[] args) {
        Network model = ex5();
        SolverOptions options = new SolverOptions(SolverType.MVA);
        options.seed = 1;
        options.iter_max = 100;
        NetworkSolver solver = new SolverMVA(model, options);
        try{
            NetworkAvgTable t = solver.getAvgTable();
            t.print(options);
        } catch (IllegalArgumentException e){
            // Only throws an exception because getAvgTable is called instead of getNodeAvgTable.
            // getAvgTable is used just to test the call sequence and that everything is working properly.
        }
        System.out.println("Hit ratio: " + ((jline.lang.nodes.Cache) model.getNodes().get(1)).getHitRatio().get(0));
        System.out.println("Miss ratio: " + ((jline.lang.nodes.Cache) model.getNodes().get(1)).getMissRatio().get(0));
    }
}
