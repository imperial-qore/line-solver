package jline.solvers.mva;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.Exp;
import jline.lang.distributions.HyperExp;
import jline.lang.nodes.*;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.solvers.NetworkAvgTable;

import java.util.Arrays;

public class SolverQNAOpenExamplesTest {
    public static void main(String[] args) {

        Network model1 = model1();
        Network model2 = model5();
        Network model3 = model6();

        SolverOptions options = new SolverOptions(SolverType.MVA);
        options.method = "qna";

        NetworkSolver solver1 = new SolverMVA(model1, options);
        NetworkAvgTable t1 = solver1.getAvgTable();
        t1.print(options);

        NetworkSolver solver2 = new SolverMVA(model2, options);
        NetworkAvgTable t2 = solver2.getAvgTable();
        t2.print(options);

        NetworkSolver solver3 = new SolverMVA(model3, options);
        NetworkAvgTable t3 = solver3.getAvgTable();
        t3.print(options);

    }
    public static Network model6(){
        Network model = new Network("model2");
        Source node1 = new Source(model,"Source");
        Queue node2 = new Queue(model,"Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model,"Queue2",SchedStrategy.FCFS);
        Queue node4 =  new Queue(model,"Queue3", SchedStrategy.FCFS);
        Queue node5 = new Queue(model,"Queue4",SchedStrategy.FCFS);
        Sink node6 = new Sink(model,"Sink");

        OpenClass jobclass1 = new OpenClass(model,"Class1",0);
        OpenClass jobclass2 = new OpenClass(model,"Class2",0);
        OpenClass jobclass3 = new OpenClass(model,"Class3",0);

        node1.setArrival(jobclass1, Exp.fitMean(5));
        node1.setArrival(jobclass2,Exp.fitMean(8));
        node1.setArrival(jobclass3,Exp.fitMean(7));
        node2.setService(jobclass1,Exp.fitMean(0.3));
        node2.setService(jobclass2,Exp.fitMean(0.5));
        node2.setService(jobclass3,Exp.fitMean(0.6));
        node3.setService(jobclass1,Exp.fitMean(1.1));
        node3.setService(jobclass2,Exp.fitMean(1.3));
        node3.setService(jobclass3,Exp.fitMean(1.5));
        node4.setService(jobclass1,Exp.fitMean(2.0));
        node4.setService(jobclass2,Exp.fitMean(2.1));
        node4.setService(jobclass3,Exp.fitMean(1.9));
        node5.setService(jobclass1,Exp.fitMean(1.5));
        node5.setService(jobclass2,Exp.fitMean(0.9));
        node5.setService(jobclass3,Exp.fitMean(2.3));

        RoutingMatrix p = new RoutingMatrix(model, model.getJobClass(), model.getNodes());
        p.addConnection(jobclass1,jobclass1,node1,node2,1);
        p.addConnection(jobclass1,jobclass1,node2,node3,0.25);
        p.addConnection(jobclass1,jobclass1,node2,node4,0.25);
        p.addConnection(jobclass1,jobclass1,node2,node5,0.25);
        p.addConnection(jobclass1,jobclass1,node2,node6,0.25);
        p.addConnection(jobclass1,jobclass1,node3,node2,1);
        p.addConnection(jobclass1,jobclass1,node4,node2,1);
        p.addConnection(jobclass1,jobclass1,node5,node2,1);
        p.addConnection(jobclass2,jobclass2,node1,node2,1);
        p.addConnection(jobclass2,jobclass2,node2,node3,0.25);
        p.addConnection(jobclass2,jobclass2,node2,node4,0.25);
        p.addConnection(jobclass2,jobclass2,node2,node5,0.25);
        p.addConnection(jobclass2,jobclass2,node2,node6,0.25);
        p.addConnection(jobclass2,jobclass2,node3,node2,1);
        p.addConnection(jobclass2,jobclass2,node4,node2,1);
        p.addConnection(jobclass2,jobclass2,node5,node2,1);
        p.addConnection(jobclass3,jobclass3,node1,node2,1);
        p.addConnection(jobclass3,jobclass3,node2,node3,0.25);
        p.addConnection(jobclass3,jobclass3,node2,node4,0.25);
        p.addConnection(jobclass3,jobclass3,node2,node5,0.25);
        p.addConnection(jobclass3,jobclass3,node2,node6,0.25);
        p.addConnection(jobclass3,jobclass3,node3,node2,1);
        p.addConnection(jobclass3,jobclass3,node4,node2,1);
        p.addConnection(jobclass3,jobclass3,node5,node2,1);

        model.link(p);

        return model;
    }

    public static Network model1(){
        Network model = new Network("model");

        Delay delay = new Delay(model,"Delay");
        Queue queue = new Queue(model,"Queue1", SchedStrategy.FCFS);
        Source source = new Source(model,"Source");
        Sink sink = new Sink(model,"Sink");

        OpenClass openClass = new OpenClass(model,"Class1",0);
        delay.setService(openClass, new HyperExp(0.5,3.0,10.0));
        queue.setService(openClass, new Exp(1));
        source.setArrival(openClass, new Exp(0.1));

        RoutingMatrix p = new RoutingMatrix(model, model.getJobClass(), model.getNodes());
        p.addConnection(source,delay,openClass,1);
        p.addConnection(delay,queue,openClass,1);
        p.addConnection(queue,sink,openClass,1);

        model.link(p);

        return model;
    }

    public static Network model5() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");
        Router node4 = new Router(model, "VSink1");
        Router node5 = new Router(model, "VSink2");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class1", 0);

        node1.setArrival(jobclass1, Exp.fitMean(1.000000)); // (Source,Class1)
        node2.setService(jobclass1, Exp.fitMean(0.010000)); // (Queue1,Class1)

        node1.setArrival(jobclass2,new Exp(1.0));
        node2.setService(jobclass2, Exp.fitMean(0.010000));

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1,jobclass2),
                Arrays.asList(node1, node2, node3, node4, node5));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Source,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node4, 0.600000); // (Queue1,Class1) -> (VSink1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node5, 0.400000); // (Queue1,Class1) -> (VSink2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node3, 1.000000); // (VSink1,Class1) -> (Sink,Class1)
        routingMatrix.set(jobclass1, jobclass1, node5, node3, 1.000000); // (VSink2,Class1) -> (Sink,Class1)


        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Source,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node2, node4, 0.100000); // (Queue1,Class1) -> (VSink1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node2, node5, 0.900000); // (Queue1,Class1) -> (VSink2,Class1)
        routingMatrix.set(jobclass2, jobclass2, node4, node3, 1.000000); // (VSink1,Class1) -> (Sink,Class1)
        routingMatrix.set(jobclass2, jobclass2, node5, node3, 1.000000); // (VSink2,Class1) -> (Sink,Class1)

        model.link(routingMatrix);

        return model;
    }

}
