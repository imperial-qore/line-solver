package jline.solvers.mam;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.APH;
import jline.lang.distributions.Exp;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.NetworkAvgTable;

public class SolverQNAMAMOpenExamplesTest {
    public static void main(String[] args) {

        Network model1 = model3();

        SolverOptions options = new SolverOptions(SolverType.MAM);
        options.method = "qnamam";

        NetworkSolver solver1 = new SolverMAM(model1, options);
        NetworkAvgTable t1 = solver1.getAvgTable();
        t1.print(options);


    }
    public static Network model1(){
        Network model = new Network("model1");
        Source node1 = new Source(model,"Source");
        Queue node2 = new Queue(model,"Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model,"Queue2",SchedStrategy.FCFS);
        Queue node4 =  new Queue(model,"Queue3", SchedStrategy.FCFS);
        Queue node5 = new Queue(model,"Queue4",SchedStrategy.FCFS);
        Sink node6 = new Sink(model,"Sink");

        OpenClass jobclass1 = new OpenClass(model,"Class1",0);
        OpenClass jobclass2 = new OpenClass(model,"Class2",0);


        node1.setArrival(jobclass1, Exp.fitMean(5));
        node1.setArrival(jobclass2,Exp.fitMean(8));

        node2.setService(jobclass1,Exp.fitMean(0.3));
        node2.setService(jobclass2,Exp.fitMean(0.5));

        node3.setService(jobclass1,Exp.fitMean(1.1));
        node3.setService(jobclass2,Exp.fitMean(1.3));

        node4.setService(jobclass1,Exp.fitMean(2.0));
        node4.setService(jobclass2,Exp.fitMean(2.1));

        node5.setService(jobclass1,Exp.fitMean(1.5));
        node5.setService(jobclass2,Exp.fitMean(0.9));


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


        model.link(p);

        return model;
    }
    public static Network model2(){
        Network model = new Network("model2");
        Source node1 = new Source(model,"Source");
        Queue node2 = new Queue(model,"Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model,"Queue2",SchedStrategy.FCFS);

        Sink node4 = new Sink(model,"Sink");

        OpenClass jobclass1 = new OpenClass(model,"Class1",0);


        node1.setArrival(jobclass1, APH.fitMeanAndSCV(2,1.0/6));
        node2.setService(jobclass1,APH.fitMeanAndSCV(0.3,1.0/5));
        node3.setService(jobclass1,APH.fitMeanAndSCV(0.2,1.0/10));



        RoutingMatrix p = new RoutingMatrix(model, model.getJobClass(), model.getNodes());
        p.addConnection(jobclass1,jobclass1,node1,node2,1);
        p.addConnection(jobclass1,jobclass1,node2,node3,1);
        p.addConnection(jobclass1,jobclass1,node3,node4,0.5);
        p.addConnection(jobclass1,jobclass1,node3,node2,0.5);


        model.link(p);

        return model;
    }
    public static Network model3(){
        Network model = new Network("model2");
        Source node1 = new Source(model,"Source");
        Queue node2 = new Queue(model,"Queue1", SchedStrategy.FCFS);
        Queue node3 = new Queue(model,"Queue2",SchedStrategy.FCFS);
        Queue node4 = new Queue(model,"Queue3",SchedStrategy.FCFS);

        Sink node5 = new Sink(model,"Sink");

        OpenClass jobclass1 = new OpenClass(model,"Class1",0);


        node1.setArrival(jobclass1, APH.fitMeanAndSCV(20,1.0/4));
        node2.setService(jobclass1,APH.fitMeanAndSCV(1,4));
        node3.setService(jobclass1,APH.fitMeanAndSCV(2,1.0/6));
        node4.setService(jobclass1,APH.fitMeanAndSCV(4,6));


        RoutingMatrix p = new RoutingMatrix(model, model.getJobClass(), model.getNodes());
        p.addConnection(jobclass1,jobclass1,node1,node2,1);
        p.addConnection(jobclass1,jobclass1,node2,node3,1);
        p.addConnection(jobclass1,jobclass1,node3,node4,0.2);
        p.addConnection(jobclass1,jobclass1,node3,node2,0.5);
        p.addConnection(jobclass1,jobclass1,node3,node5,0.3);
        p.addConnection(jobclass1,jobclass1,node4,node5,0.9);
        p.addConnection(jobclass1,jobclass1,node4,node2,0.1);
        model.link(p);

        return model;
    }


}
