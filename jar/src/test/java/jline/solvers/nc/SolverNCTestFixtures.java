package jline.solvers.nc;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.DiscreteSampler;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

public class SolverNCTestFixtures {

    public static Network complexMultiClassModel_1() {
        Network model = new Network("example_complexMultiClassModel");

        Queue node1 = new Queue(model, "Queue1", SchedStrategy.INF);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.INF);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.PS);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.PS);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 100, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 33, node4, 1);

        node1.setService(jobclass1, new Exp(1.0));
        node1.setService(jobclass2, new Exp(Math.sqrt(1.0)));
        node2.setService(jobclass1, new Exp(2.0));
        node2.setService(jobclass2, new Exp(Math.sqrt(2.0)));
        node3.setService(jobclass1, new Exp(3.0));
        node3.setService(jobclass2, new Exp(Math.sqrt(3.0)));
        node4.setService(jobclass1, new Exp(4.0));
        node4.setService(jobclass2, new Exp(Math.sqrt(4.0)));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 1.00);
        routingMatrix.addConnection(jobclass1, jobclass1, node2, node3, 1.00);
        routingMatrix.addConnection(jobclass1, jobclass1, node3, node4, 1.00);
        routingMatrix.addConnection(jobclass1, jobclass1, node4, node1, 1.00);

        routingMatrix.addConnection(jobclass2, jobclass2, node1, node2, 1.00);
        routingMatrix.addConnection(jobclass2, jobclass2, node2, node4, 1.00);
        routingMatrix.addConnection(jobclass2, jobclass2, node4, node1, 1.00);

        model.link(routingMatrix);

        return model;
    }
    
    public static Network complexMultiClassModel_2() {
        Network model = new Network("example_complexMultiClassModel");

        Queue node1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue node2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue node3 = new Queue(model, "Queue3", SchedStrategy.INF);
        Queue node4 = new Queue(model, "Queue4", SchedStrategy.INF);

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 88, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 33, node4, 1);

        node1.setService(jobclass1, new Exp(1.0));
        node1.setService(jobclass2, new Exp(Math.sqrt(1.0)));
        node2.setService(jobclass1, new Exp(2.0));
        node2.setService(jobclass2, new Exp(Math.sqrt(2.0)));
        node3.setService(jobclass1, new Exp(3.0));
        node3.setService(jobclass2, new Exp(Math.sqrt(3.0)));
        node4.setService(jobclass1, new Exp(4.0));
        node4.setService(jobclass2, new Exp(Math.sqrt(4.0)));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 1.00);
        routingMatrix.addConnection(jobclass1, jobclass1, node2, node3, 1.00);
        routingMatrix.addConnection(jobclass1, jobclass1, node3, node4, 1.00);
        routingMatrix.addConnection(jobclass1, jobclass1, node4, node1, 1.00);

        routingMatrix.addConnection(jobclass2, jobclass2, node1, node2, 1.00);
        routingMatrix.addConnection(jobclass2, jobclass2, node2, node4, 1.00);
        routingMatrix.addConnection(jobclass2, jobclass2, node4, node1, 1.00);

        model.link(routingMatrix);

        return model;
    }

    public static Network cacheModel_1(){
        Network model = new Network("model");

        int n = 10; // Number of items
        Matrix m = new Matrix(1,2); // cache capacity
        m.set(0, 2);
        m.set(1, 2);

        Source source = new Source(model, "Source");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.LRU, null);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(2));

        Matrix p = new Matrix(1, n);
        p.fill(1.0/n);
        DiscreteSampler pAccess = new DiscreteSampler(p);

        cacheNode.setRead(jobClass, pAccess);
        //cacheNode.setRead(jobClass, new Zipf(0.8, n));

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);

        model.link(routingMatrix);
        return model;
    }

    public static Network cacheModel_2(){
        Network model = new Network("model");

        int n = 100; // Number of items
        Matrix m = new Matrix(1,2); // cache capacity
        m.set(0, 20);
        m.set(1, 20);

        Source source = new Source(model, "Source");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.LRU, null);

        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(2));

        Matrix p = new Matrix(1, n);
        p.fill(1.0/n);
        DiscreteSampler pAccess = new DiscreteSampler(p);

        cacheNode.setRead(jobClass, pAccess);
        //cacheNode.setRead(jobClass, new Zipf(0.8, n));

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);

        model.link(routingMatrix);
        return model;
    }

//    public static Network cacheModel_3(){
//        Network model = new Network("model");
//
//        int n = 100; // Number of items
//        Matrix m = new Matrix(1,2); // cache capacity
//        m.set(0, 20);
//        m.set(1, 20);
//
//        Source source = new Source(model, "Source");
//        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.HLRU, null);
//        Sink sink = new Sink(model, "Sink");
//
//        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
//        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
//        OpenClass missClass = new OpenClass(model, "MissClass", 0);
//
//        source.setArrival(jobClass, new Exp(2));
//
//        Matrix p = new Matrix(1, n);
//        p.fill(1.0/n);
//        DiscreteSampler pAccess = new DiscreteSampler(p);
//
//        cacheNode.setRead(jobClass, pAccess);
//        //cacheNode.setRead(jobClass, new Zipf(0.8, n));
//
//        cacheNode.setHitClass(jobClass, hitClass);
//        cacheNode.setMissClass(jobClass, missClass);
//
//        RoutingMatrix routingMatrix = model.initRoutingMatrix();
//
//        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
//        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
//        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);
//
//        model.link(routingMatrix);
//        return model;
//    }

    public static Network cacheModel_4(){
        Network model = new Network("model");

        int n = 100; // Number of items
        Matrix m = new Matrix(1,2); // cache capacity
        m.set(0, 20);
        m.set(1, 20);

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
        //cacheNode.setRead(jobClass, new Zipf(0.8, n));

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);

        model.link(routingMatrix);
        return model;
    }

    public static Network cacheModel_5(){
        Network model = new Network("model");

        int n = 100; // Number of items
        Matrix m = new Matrix(1,2); // cache capacity
        m.set(0, 20);
        m.set(1, 20);

        Source source = new Source(model, "Source");
        jline.lang.nodes.Cache cacheNode = new jline.lang.nodes.Cache(model, "Cache", n, m, ReplacementStrategy.RR, null);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);

        source.setArrival(jobClass, new Exp(2));

        Matrix p = new Matrix(1, n);
        p.fill(1.0/n);
        DiscreteSampler pAccess = new DiscreteSampler(p);

        cacheNode.setRead(jobClass, pAccess);
        //cacheNode.setRead(jobClass, new Zipf(0.8, n));

        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobClass, jobClass, source, cacheNode, 1.0);
        routingMatrix.set(hitClass, hitClass, cacheNode, sink, 1.0);
        routingMatrix.set(missClass, missClass, cacheNode, sink, 1.0);

        model.link(routingMatrix);
        return model;
    }


}
