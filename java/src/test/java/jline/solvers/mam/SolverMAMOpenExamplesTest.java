package jline.solvers.mam;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Disabled;
import jline.lang.distributions.Exp;
import jline.lang.distributions.HyperExp;
import jline.lang.nodes.*;
import jline.solvers.SolverOptions;
import jline.util.Matrix;

import java.util.Arrays;

public class SolverMAMOpenExamplesTest {
    public static void main(String[] args) {

        Network model1 = model1();
        Network model2 = model2();
        Network model3 = model3();
        Network model4 = model5();
        Network model5 = model6();


        SolverMAM solver1 = new SolverMAM(model1);
        solver1.getAvgTable().print(new SolverOptions());
        System.out.println(solver1.result.runtime);

        SolverMAM solver2 = new SolverMAM(model2);
        solver2.getAvgTable().print(new SolverOptions());
        System.out.println(solver2.result.runtime);

        SolverMAM solver3 = new SolverMAM(model3);
        solver3.getAvgTable().print(new SolverOptions());
        System.out.println(solver3.result.runtime);

        SolverMAM solver4 = new SolverMAM(model4);
        solver4.getAvgTable().print(new SolverOptions());
        System.out.println(solver4.result.runtime);

        SolverMAM solver5 = new SolverMAM(model5);
        solver5.getAvgTable().print(new SolverOptions());
        System.out.println(solver5.result.runtime);
    }

    public static Network model6(){
        Network model = new Network("model2");
        Source node1 = new Source(model,"Source");
        Queue node2 = new Queue(model,"Queue1",SchedStrategy.FCFS);
        Queue node3 = new Queue(model,"Queue2",SchedStrategy.FCFS);
        Queue node4 =  new Queue(model,"Queue3", SchedStrategy.FCFS);
        Queue node5 = new Queue(model,"Queue4",SchedStrategy.FCFS);
        Sink node6 = new Sink(model,"Sink");

        OpenClass jobclass1 = new OpenClass(model,"Class1",0);
        OpenClass jobclass2 = new OpenClass(model,"Class2",0);
        OpenClass jobclass3 = new OpenClass(model,"Class3",0);

        node1.setArrival(jobclass1,Exp.fitMean(5));
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

    public static Network model2(){
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Delay node2 = new Delay(model, "Station1");
        Queue node3 = new Queue(model, "Station2", SchedStrategy.PS);
        Queue node4 = new Queue(model, "Station3", SchedStrategy.PS);
        Sink node5 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, Exp.fitMean(50.000000)); // (Source,Class1)
        node1.setArrival(jobclass2, Exp.fitMean(25.000000)); // (Source,Class2)
        node2.setService(jobclass1, Exp.fitMean(91.000000)); // (Station1,Class1)
        node2.setService(jobclass2, Exp.fitMean(92.000000)); // (Station1,Class2)
        node3.setService(jobclass1, Exp.fitMean(10.000000)); // (Station2,Class1)
        node3.setService(jobclass2, Exp.fitMean(5.000000)); // (Station2,Class2)
        node4.setService(jobclass1, Exp.fitMean(5.000000)); // (Station3,Class1)
        node4.setService(jobclass2, Exp.fitMean(9.000000)); // (Station3,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4, node5));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Source,Class1) -> (Station1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Station1,Class1) -> (Station2,Class1)
        routingMatrix.set(jobclass1, jobclass1, node3, node4, 1.000000); // (Station2,Class1) -> (Station3,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node5, 1.000000); // (Station3,Class1) -> (Sink,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Source,Class2) -> (Station1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Station1,Class2) -> (Station2,Class2)
        routingMatrix.set(jobclass2, jobclass2, node3, node4, 1.000000); // (Station2,Class2) -> (Station3,Class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node5, 1.000000); // (Station3,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;


    }

    public static Network model_compare(){
        Network model = new Network("model");
        Source source = new Source(model,"Source");
        Queue queue1 = new Queue(model,"Queue1",SchedStrategy.FCFS);
        Queue queue2 = new Queue(model,"Queue2",SchedStrategy.FCFS);
        Sink sink = new Sink(model,"Sink");

        OpenClass openClass1 = new OpenClass(model,"Class1",0);
        OpenClass openClass2 = new OpenClass(model,"Class2",0);

        source.setArrival(openClass1,Exp.fitMean(10));
        source.setArrival(openClass2,Exp.fitMean(4));

        queue1.setService(openClass1,Exp.fitMean(0.1));
        queue1.setService(openClass2,Exp.fitMean(0.2));

        queue2.setService(openClass1,Exp.fitMean(0.3));
        queue2.setService(openClass2,Exp.fitMean(0.2));

        RoutingMatrix p = new RoutingMatrix(model, model.getJobClass(), model.getNodes());
        p.addConnection(openClass1,openClass1,source,queue1,1);
        p.addConnection(openClass1,openClass1,queue1,queue2,1);
        p.addConnection(openClass1,openClass1,queue2,queue1,0.1);
        p.addConnection(openClass1,openClass1,queue2,sink,0.9);
        p.addConnection(openClass2,openClass2,source,queue1,1);
        p.addConnection(openClass2,openClass2,queue1,queue2,1);
        p.addConnection(openClass2,openClass2,queue2,queue1,0.1);
        p.addConnection(openClass2,openClass2,queue2,sink,0.9);

        model.link(p);

        return model;
    }

    public static Network model3() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source 1");
        Queue node2 = new Queue(model, "Queue 1", SchedStrategy.PS);
        ClassSwitch node3 = new ClassSwitch(model, "ClassSwitch 1"); // Dummy node, class switching is embedded in the routing matrix P
        Sink node4 = new Sink(model, "Sink 1");
        Queue node5 = new Queue(model, "Queue 2", SchedStrategy.PS);

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class A", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class B", 0);
        OpenClass jobclass3 = new OpenClass(model, "Class C", 0);

        node1.setArrival(jobclass1, Exp.fitMean(0.500000)); // (Source 1,Class A)
        node1.setArrival(jobclass2, Exp.fitMean(1.000000)); // (Source 1,Class B)
        node1.setArrival(jobclass3, Disabled.getInstance()); // (Source 1,Class C)
        node2.setService(jobclass1, Exp.fitMean(0.200000)); // (Queue 1,Class A)
        node2.setService(jobclass2, Exp.fitMean(0.300000)); // (Queue 1,Class B)
        node2.setService(jobclass3, Exp.fitMean(0.333333)); // (Queue 1,Class C)
        node5.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue 2,Class A)
        node5.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue 2,Class B)
        node5.setService(jobclass3, Exp.fitMean(0.150000)); // (Queue 2,Class C)

        // Block 3: topology
        Matrix C = node3.initClassSwitchMatrix();
        C.setTo(Matrix.eye(C.length()));
        node3.setClassSwitchingMatrix(C);

        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2, jobclass3),
                Arrays.asList(node1, node2, node3, node4, node5));

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.000000); // (Source 1,Class A) -> (Queue 1,Class A)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue 1,Class A) -> (ClassSwitch 1,Class A)
        routingMatrix.set(jobclass1, jobclass1, node5, node4, 1.000000); // (Queue 2,Class A) -> (Sink 1,Class A)
        routingMatrix.set(jobclass1, jobclass3, node3, node5, 1.000000); // (ClassSwitch 1,Class A) -> (Queue 2,Class C)
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.000000); // (Source 1,Class B) -> (Queue 1,Class B)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue 1,Class B) -> (ClassSwitch 1,Class B)
        routingMatrix.set(jobclass2, jobclass2, node5, node4, 1.000000); // (Queue 2,Class B) -> (Sink 1,Class B)
        routingMatrix.set(jobclass2, jobclass3, node3, node5, 1.000000); // (ClassSwitch 1,Class B) -> (Queue 2,Class C)
        routingMatrix.set(jobclass3, jobclass3, node1, node2, 1.000000); // (Source 1,Class C) -> (Queue 1,Class C)
        routingMatrix.set(jobclass3, jobclass3, node2, node3, 1.000000); // (Queue 1,Class C) -> (ClassSwitch 1,Class C)
        routingMatrix.set(jobclass3, jobclass3, node3, node5, 1.000000); // (ClassSwitch 1,Class C) -> (Queue 2,Class C)
        routingMatrix.set(jobclass3, jobclass3, node5, node4, 1.000000); // (Queue 2,Class C) -> (Sink 1,Class C)

        model.link(routingMatrix);

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
