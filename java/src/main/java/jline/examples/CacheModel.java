package jline.examples;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.distributions.Disabled;
import jline.lang.distributions.DiscreteSampler;
import jline.lang.distributions.Exp;
import jline.lang.nodes.*;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.mva.MVAOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.Matrix;
import jline.solvers.NetworkAvgTable;

import java.util.Arrays;

/**
 * Examples of caching models
 */
public class CacheModel {

    public static Network example_cacheModel_1(){
        Network model = new Network("model");

        int n = 5; // Number of items
        int m = 2;

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

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);

        model.link(routingMatrix);

        return model;
    }

    public static Network example_cacheModel_2(){
        Network model = new Network("model");

        int n = 5; // Number of items
        int m = 2; // Cache capacity

        Delay delay = new Delay(model, "Delay");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.FIFO);

        ClosedClass jobClass = new ClosedClass(model, "JobClass", 1, delay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, delay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, delay, 0);

        delay.setService(jobClass, new Exp(1.0));

        Matrix p = new Matrix(1, n); p.fill(1.0/n);
        DiscreteSampler pAccess = new DiscreteSampler(p);
        cacheNode.setRead(jobClass, pAccess);

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, delay, cacheNode, 1.0);
        routingMatrix.set(hitClass, jobClass, cacheNode, delay, 1.0);
        routingMatrix.set(missClass, jobClass, cacheNode, delay, 1.0);

        model.link(routingMatrix);

        return model;
    }

    public static Network example_cacheModel_3() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network example_cacheModel_4() {
        Network model = new Network("model");
        // TODO
        return model;
    }

    public static Network example_cacheModel_5(){
        Network model = new Network("model");

        int n = 5; // Number of items
        int m = 2;

        Source source = new Source(model, "Source");
        Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.FIFO);
        Router routerNode = new Router(model, "Router");
        Delay delay1 = new Delay(model, "Delay1");
        Delay delay2 = new Delay(model, "Delay2");
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(2));
        delay1.setService(hitClass, new Exp(10));
        delay1.setService(missClass, new Exp(1.0));
        delay2.setService(hitClass, new Exp(20));
        delay2.setService(missClass, new Exp(2));

        Matrix p = new Matrix(1, n);
        p.fill(1.0/n);
        DiscreteSampler pAccess = new DiscreteSampler(p);

        cacheNode.setRead(jobClass, pAccess);

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        model.addLink(source, cacheNode);
        model.addLink(delay1, sink);
        model.addLink(delay2, sink);

        source.setProbRouting(jobClass, cacheNode, 1.0);

        cacheNode.setRouting(jobClass, RoutingStrategy.RAND);
        cacheNode.setProbRouting(hitClass, routerNode, 1.0);
        cacheNode.setProbRouting(missClass, routerNode, 1.0);

        routerNode.setRouting(hitClass, RoutingStrategy.RAND);
        routerNode.setRouting(missClass, RoutingStrategy.RAND);

        delay1.setProbRouting(hitClass, sink, 1.0);
        delay1.setProbRouting(missClass, sink, 1.0);

        delay2.setProbRouting(hitClass, sink, 1.0);
        delay2.setProbRouting(missClass, sink, 1.0);

        // TODO
        // new SolverCTMC(model,"keep",false,"cutoff",1,"seed",1).getAvgNodeTable().print();
        // new SolverSSA(model,"samples",1e4,"verbose",true,"method","serial","seed",1).getAvgNodeTable().print();
        // new SolverMVA(model,"seed",1).getAvgNodeTable().print();
        // new SolverNC(model,"seed",1).getAvgNodeTable().print();

        System.out.printf("hitRatio = %d\n", cacheNode.getHitRatio().get(0,0));
        System.out.printf("hitRatio = %d\n", cacheNode.getMissRatio().get(0,0));

        return model;
    }

    public static void main(String[] args) {
        example_cacheModel_5();
    }
}
