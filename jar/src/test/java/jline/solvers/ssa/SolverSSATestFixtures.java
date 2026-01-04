package jline.solvers.ssa;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;

public class SolverSSATestFixtures {

    public static Network cqn_repairmenps() {
        Network model = new Network("myModel");

        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.PS);

        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);

        node1.setService(jobclass1, Exp.fitMean(1.00)); // (Delay,Class1)
        node2.setService(jobclass1, Exp.fitMean(1.50)); // (Queue1,Class1)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node1, 0.70); // (Delay,Class1) -> (Delay,Class1)
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 0.3); // (Delay,Class1) -> (Queue1,Class1)
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.00); // (Queue1,Class1) -> (Delay,Class1)

        model.link(routingMatrix);

        return model;
    }

}
