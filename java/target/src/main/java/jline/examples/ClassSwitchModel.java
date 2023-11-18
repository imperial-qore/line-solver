package jline.examples;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Exp;
import jline.lang.nodes.*;
import jline.solvers.NetworkSolver;
import jline.solvers.mva.SolverMVA;
import jline.util.Matrix;

import java.util.Arrays;

/**
 * Examples of models with class switching
 */
public class ClassSwitchModel {

    public static Network ex1() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Source node1 = new Source(model, "Source 1");
        Queue node2 = new Queue(model, "Queue 1", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink 1");
        jline.lang.nodes.ClassSwitch node4 = new jline.lang.nodes.ClassSwitch(model, "ClassSwitch 1"); // Dummy node, class switching is embedded in the routing matrix P

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, Exp.fitMean(10.000000)); // (Source 1,Class1)
        node1.setArrival(jobclass2, Exp.fitMean(2.000000)); // (Source 1,Class2)
        node2.setService(jobclass1, Exp.fitMean(1.000000)); // (Queue 1,Class1)
        node2.setService(jobclass2, Exp.fitMean(1.000000)); // (Queue 1,Class2)

        Matrix csMatrix = node4.initClassSwitchMatrix();
        csMatrix.set(jobclass1.getIndex()-1, jobclass1.getIndex()-1, 0.3);
        csMatrix.set(jobclass1.getIndex()-1, jobclass2.getIndex()-1, 0.7);
        csMatrix.set(jobclass2.getIndex()-1, jobclass1.getIndex()-1, 1.0);
        csMatrix.set(jobclass2.getIndex()-1, jobclass2.getIndex()-1, 0.0);
        node4.setClassSwitchingMatrix(csMatrix);

        // Block 3: topology
        RoutingMatrix routingMatrix = new RoutingMatrix(model,
                Arrays.asList(jobclass1, jobclass2),
                Arrays.asList(node1, node2, node3, node4));

        routingMatrix.set(jobclass1, jobclass1, node1, node4, 1.000000); // (Source 1,Class1) -> (ClassSwitch 1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.000000); // (Queue 1,Class1) -> (Sink 1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node4, node2, 1.000000); // (ClassSwitch 1,Class1) -> (Queue 1,Class1)
        routingMatrix.set(jobclass2, jobclass2, node1, node4, 1.000000); // (Source 1,Class2) -> (ClassSwitch 1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.000000); // (Queue 1,Class2) -> (Sink 1,Class2)
        routingMatrix.set(jobclass2, jobclass2, node4, node2, 1.000000); // (Queue 1,Class2) -> (Sink 1,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static void main(String[] args) throws Exception {
        Network model = ex1();

        NetworkStruct sn = model.getStruct(false);
        sn.rt.print();
        NetworkSolver solver = new SolverMVA(model);
        solver.getAvgTable();
    }

}
