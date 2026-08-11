/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import java.util.List;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;

/**
 * Extended test fixtures for the Schmidt MVA algorithm.
 */
public class SolverMVATestSchmidtExtensionFixtures {

    public static Network network_2Class_1Node() {
        Network network = new Network("network");

        Queue q1 = new Queue(network, "q1", SchedStrategy.FCFS);
        q1.setNumberOfServers(2);

        ClosedClass c1 = new ClosedClass(network, "c1", 2, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 2, q1);

        q1.setService(c1, Exp.fitMean(1.0));
        q1.setService(c2, Exp.fitMean(2.0));

        RoutingMatrix routingMatrix = new RoutingMatrix(network, List.of(c1, c2), List.of(q1));

        routingMatrix.addConnection(c1, c1, q1, q1, 1.0);
        routingMatrix.addConnection(c2, c2, q1, q1, 1.0);

        network.link(routingMatrix);

        return network;
    }

    public static Network network_multiclass_1Node() {
        Network network = new Network("network");

        Queue q1 = new Queue(network, "q1", SchedStrategy.FCFS);
        q1.setNumberOfServers(2);

        ClosedClass c1 = new ClosedClass(network, "c1", 1, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 1, q1);
        ClosedClass c3 = new ClosedClass(network, "c3", 1, q1);
        ClosedClass c4 = new ClosedClass(network, "c4", 1, q1);

        q1.setService(c1, Exp.fitMean(1.0));
        q1.setService(c2, Exp.fitMean(8.0));
        q1.setService(c3, Exp.fitMean(16.0));
        q1.setService(c4, Exp.fitMean(128.0));

        RoutingMatrix routingMatrix = new RoutingMatrix(network, List.of(c1, c2, c3, c4), List.of(q1));

        routingMatrix.addConnection(c1, c1, q1, q1, 1.0);
        routingMatrix.addConnection(c2, c2, q1, q1, 1.0);
        routingMatrix.addConnection(c3, c3, q1, q1, 1.0);
        routingMatrix.addConnection(c4, c4, q1, q1, 1.0);

        network.link(routingMatrix);

        return network;
    }

    public static Network fcfs_varies_3_classes_large_difference() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.PS);
        q1.setNumberOfServers(3);
        Queue q2 = new Queue(network, "q2", SchedStrategy.FCFS);
        q2.setNumberOfServers(4);

        Queue q3 = new Queue(network, "q3", SchedStrategy.FCFS);
        q3.setNumberOfServers(3);

        ClosedClass c1 = new ClosedClass(network, "c1", 3, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 3, q1);
        ClosedClass c3 = new ClosedClass(network, "c3", 3, q1);

        q1.setService(c1, Exp.fitMean(64.5));
        q1.setService(c2, Exp.fitMean(64.5));
        q1.setService(c3, Exp.fitMean(64.5));

        q2.setService(c1, Exp.fitMean(18.84));
        q2.setService(c2, Exp.fitMean(4.90));
        q2.setService(c3, Exp.fitMean(117.03));

        q3.setService(c1, Exp.fitMean(3.74));
        q3.setService(c2, Exp.fitMean(127.02));
        q3.setService(c3, Exp.fitMean(10.10));

        RoutingMatrix routingMatrix = network.initRoutingMatrix();

        routingMatrix.addConnection(c1, c1, q1, q2, 0.5);
        routingMatrix.addConnection(c2, c2, q1, q2, 0.5);
        routingMatrix.addConnection(c3, c3, q1, q2, 0.5);

        routingMatrix.addConnection(c1, c1, q1, q3, 0.5);
        routingMatrix.addConnection(c2, c2, q1, q3, 0.5);
        routingMatrix.addConnection(c3, c3, q1, q3, 0.5);

        routingMatrix.addConnection(c1, c1, q2, q1, 0.5);
        routingMatrix.addConnection(c2, c2, q2, q1, 0.5);
        routingMatrix.addConnection(c3, c3, q2, q1, 0.5);

        routingMatrix.addConnection(c1, c1, q2, q3, 0.5);
        routingMatrix.addConnection(c2, c2, q2, q3, 0.5);
        routingMatrix.addConnection(c3, c3, q2, q3, 0.5);

        routingMatrix.addConnection(c1, c1, q3, q1, 0.5);
        routingMatrix.addConnection(c2, c2, q3, q1, 0.5);
        routingMatrix.addConnection(c3, c3, q3, q1, 0.5);

        routingMatrix.addConnection(c1, c1, q3, q2, 0.5);
        routingMatrix.addConnection(c2, c2, q3, q2, 0.5);
        routingMatrix.addConnection(c3, c3, q3, q2, 0.5);

        network.link(routingMatrix);
        return network;
    }

    public static Network performanceValidation() {
        Network network = new Network("performanceValidation");

        Queue q1 = new Queue(network, "q1", SchedStrategy.PS);
        q1.setNumberOfServers(4);
        Queue q2 = new Queue(network, "q2", SchedStrategy.FCFS);
        q2.setNumberOfServers(3);
        Queue q3 = new Queue(network, "q3", SchedStrategy.FCFS);
        q3.setNumberOfServers(2);

        ClosedClass c1 = new ClosedClass(network, "c1", 4, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 5, q1);
        ClosedClass c3 = new ClosedClass(network, "c3", 5, q1);
        ClosedClass c4 = new ClosedClass(network, "c4", 3, q1);

        q1.setService(c1, Exp.fitMean(84.01040080312384));
        q1.setService(c2, Exp.fitMean(94.77038866304343));
        q1.setService(c3, Exp.fitMean(41.35757935978896));
        q1.setService(c4, Exp.fitMean(9.019203818986776));

        q2.setService(c1, Exp.fitMean(115.89436272155214));
        q2.setService(c2, Exp.fitMean(30.20993010827056));
        q2.setService(c3, Exp.fitMean(85.18681829894625));
        q2.setService(c4, Exp.fitMean(87.816460168896));

        q3.setService(c1, Exp.fitMean(72.55502570493324));
        q3.setService(c2, Exp.fitMean(67.09795222853062));
        q3.setService(c3, Exp.fitMean(67.2631380345578));
        q3.setService(c4, Exp.fitMean(3.396292720283628));

        RoutingMatrix routingMatrix = network.initRoutingMatrix();

        routingMatrix.addConnection(c1, c1, q1, q2, 0.5);
        routingMatrix.addConnection(c2, c2, q1, q2, 0.5);
        routingMatrix.addConnection(c3, c3, q1, q2, 0.5);
        routingMatrix.addConnection(c4, c4, q1, q2, 0.5);

        routingMatrix.addConnection(c1, c1, q1, q3, 0.5);
        routingMatrix.addConnection(c2, c2, q1, q3, 0.5);
        routingMatrix.addConnection(c3, c3, q1, q3, 0.5);
        routingMatrix.addConnection(c4, c4, q1, q3, 0.5);

        routingMatrix.addConnection(c1, c1, q2, q1, 0.5);
        routingMatrix.addConnection(c2, c2, q2, q1, 0.5);
        routingMatrix.addConnection(c3, c3, q2, q1, 0.5);
        routingMatrix.addConnection(c4, c4, q2, q1, 0.5);

        routingMatrix.addConnection(c1, c1, q2, q3, 0.5);
        routingMatrix.addConnection(c2, c2, q2, q3, 0.5);
        routingMatrix.addConnection(c3, c3, q2, q3, 0.5);
        routingMatrix.addConnection(c4, c4, q2, q3, 0.5);

        routingMatrix.addConnection(c1, c1, q3, q1, 0.5);
        routingMatrix.addConnection(c2, c2, q3, q1, 0.5);
        routingMatrix.addConnection(c3, c3, q3, q1, 0.5);
        routingMatrix.addConnection(c4, c4, q3, q1, 0.5);

        routingMatrix.addConnection(c1, c1, q3, q2, 0.5);
        routingMatrix.addConnection(c2, c2, q3, q2, 0.5);
        routingMatrix.addConnection(c3, c3, q3, q2, 0.5);
        routingMatrix.addConnection(c4, c4, q3, q2, 0.5);

        network.link(routingMatrix);
        return network;
    }

    public static Network fcfs_varies_3_classes_small_difference() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.PS);
        q1.setNumberOfServers(2);
        Queue q2 = new Queue(network, "q2", SchedStrategy.FCFS);
        q2.setNumberOfServers(2);

        Queue q3 = new Queue(network, "q3", SchedStrategy.FCFS);
        q3.setNumberOfServers(2);

        ClosedClass c1 = new ClosedClass(network, "c1", 4, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 4, q1);
        ClosedClass c3 = new ClosedClass(network, "c3", 4, q1);

        q1.setService(c1, Exp.fitMean(1));
        q1.setService(c2, Exp.fitMean(1));
        q1.setService(c3, Exp.fitMean(1));

        q2.setService(c1, Exp.fitMean(0.9));
        q2.setService(c2, Exp.fitMean(1.0));
        q2.setService(c3, Exp.fitMean(1.5));

        q3.setService(c1, Exp.fitMean(0.9));
        q3.setService(c2, Exp.fitMean(1.0));
        q3.setService(c3, Exp.fitMean(1.1));

        RoutingMatrix routingMatrix = network.initRoutingMatrix();

        routingMatrix.addConnection(c1, c1, q1, q2, 0.5);
        routingMatrix.addConnection(c2, c2, q1, q2, 0.5);
        routingMatrix.addConnection(c3, c3, q1, q2, 0.5);

        routingMatrix.addConnection(c1, c1, q1, q3, 0.5);
        routingMatrix.addConnection(c2, c2, q1, q3, 0.5);
        routingMatrix.addConnection(c3, c3, q1, q3, 0.5);

        routingMatrix.addConnection(c1, c1, q2, q1, 0.5);
        routingMatrix.addConnection(c2, c2, q2, q1, 0.5);
        routingMatrix.addConnection(c3, c3, q2, q1, 0.5);

        routingMatrix.addConnection(c1, c1, q2, q3, 0.5);
        routingMatrix.addConnection(c2, c2, q2, q3, 0.5);
        routingMatrix.addConnection(c3, c3, q2, q3, 0.5);

        routingMatrix.addConnection(c1, c1, q3, q1, 0.5);
        routingMatrix.addConnection(c2, c2, q3, q1, 0.5);
        routingMatrix.addConnection(c3, c3, q3, q1, 0.5);

        routingMatrix.addConnection(c1, c1, q3, q2, 0.5);
        routingMatrix.addConnection(c2, c2, q3, q2, 0.5);
        routingMatrix.addConnection(c3, c3, q3, q2, 0.5);

        network.link(routingMatrix);
        return network;
    }

    public static Network fcfs_varies_3_paper_example() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.PS);
        q1.setNumberOfServers(4);
        Queue q2 = new Queue(network, "q2", SchedStrategy.FCFS);
        q2.setNumberOfServers(2);

        ClosedClass c1 = new ClosedClass(network, "c1", 5, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 2, q1);
        ClosedClass c3 = new ClosedClass(network, "c3", 5, q1);

        q1.setService(c1, Exp.fitMean(119));
        q1.setService(c2, Exp.fitMean(86));
        q1.setService(c3, Exp.fitMean(120));

        q2.setService(c1, Exp.fitMean(122));
        q2.setService(c2, Exp.fitMean(22));
        q2.setService(c3, Exp.fitMean(26));

        RoutingMatrix routingMatrix = network.initRoutingMatrix();

        routingMatrix.addConnection(c1, c1, q1, q2, 1.0);
        routingMatrix.addConnection(c2, c2, q1, q2, 1.0);
        routingMatrix.addConnection(c3, c3, q1, q2, 1.0);

        routingMatrix.addConnection(c1, c1, q2, q1, 1.0);
        routingMatrix.addConnection(c2, c2, q2, q1, 1.0);
        routingMatrix.addConnection(c3, c3, q2, q1, 1.0);

        network.link(routingMatrix);
        return network;
    }

    public static Network single_server() {
        Network network = new Network("network");

        Queue q1 = new Queue(network, "q1", SchedStrategy.PS);
        q1.setNumberOfServers(1);
        Queue q2 = new Queue(network, "q2", SchedStrategy.FCFS);
        q2.setNumberOfServers(1);
        Queue q3 = new Queue(network, "q3", SchedStrategy.INF);
        Queue q4 = new Queue(network, "q4", SchedStrategy.PS);
        q4.setNumberOfServers(1);
        Queue q5 = new Queue(network, "q5", SchedStrategy.PS);
        q5.setNumberOfServers(1);

        ClosedClass c1 = new ClosedClass(network, "c1", 1, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 1, q1);
        ClosedClass c3 = new ClosedClass(network, "c3", 1, q1);

        q1.setService(c1, Exp.fitMean(1.0277551665592826));
        q1.setService(c2, Exp.fitMean(1.1070999601244844));
        q1.setService(c3, Exp.fitMean(1.335475451851549));

        q2.setService(c1, Exp.fitMean(1.102991390769498));
        q2.setService(c2, Exp.fitMean(1.0856009820123402));
        q2.setService(c3, Exp.fitMean(1.2732187446445709));

        q3.setService(c1, Exp.fitMean(1.3505594895470652));
        q3.setService(c2, Exp.fitMean(1.0624555324122023));
        q3.setService(c3, Exp.fitMean(1.3202447131722856));

        q4.setService(c1, Exp.fitMean(1.2877117673316287));
        q4.setService(c2, Exp.fitMean(1.0948244659274469));
        q4.setService(c3, Exp.fitMean(1.3809838125228326));

        q5.setService(c1, Exp.fitMean(1.3104381608843825));
        q5.setService(c2, Exp.fitMean(1.0791504871680548));
        q5.setService(c3, Exp.fitMean(1.3480546987209776));

        RoutingMatrix routingMatrix = network.initRoutingMatrix();

        routingMatrix.addConnection(c1, c1, q1, q2, 0.25);
        routingMatrix.addConnection(c1, c1, q1, q3, 0.25);
        routingMatrix.addConnection(c1, c1, q1, q4, 0.25);
        routingMatrix.addConnection(c1, c1, q1, q5, 0.25);

        routingMatrix.addConnection(c2, c2, q1, q2, 0.25);
        routingMatrix.addConnection(c2, c2, q1, q3, 0.25);
        routingMatrix.addConnection(c2, c2, q1, q4, 0.25);
        routingMatrix.addConnection(c2, c2, q1, q5, 0.25);

        routingMatrix.addConnection(c3, c3, q1, q2, 0.25);
        routingMatrix.addConnection(c3, c3, q1, q3, 0.25);
        routingMatrix.addConnection(c3, c3, q1, q4, 0.25);
        routingMatrix.addConnection(c3, c3, q1, q5, 0.25);

        routingMatrix.addConnection(c1, c1, q2, q1, 0.25);
        routingMatrix.addConnection(c1, c1, q2, q3, 0.25);
        routingMatrix.addConnection(c1, c1, q2, q4, 0.25);
        routingMatrix.addConnection(c1, c1, q2, q5, 0.25);

        routingMatrix.addConnection(c2, c2, q2, q1, 0.25);
        routingMatrix.addConnection(c2, c2, q2, q3, 0.25);
        routingMatrix.addConnection(c2, c2, q2, q4, 0.25);
        routingMatrix.addConnection(c2, c2, q2, q5, 0.25);

        routingMatrix.addConnection(c3, c3, q2, q1, 0.25);
        routingMatrix.addConnection(c3, c3, q2, q3, 0.25);
        routingMatrix.addConnection(c3, c3, q2, q4, 0.25);
        routingMatrix.addConnection(c3, c3, q2, q5, 0.25);

        routingMatrix.addConnection(c1, c1, q3, q1, 0.25);
        routingMatrix.addConnection(c1, c1, q3, q2, 0.25);
        routingMatrix.addConnection(c1, c1, q3, q4, 0.25);
        routingMatrix.addConnection(c1, c1, q3, q5, 0.25);

        routingMatrix.addConnection(c2, c2, q3, q1, 0.25);
        routingMatrix.addConnection(c2, c2, q3, q2, 0.25);
        routingMatrix.addConnection(c2, c2, q3, q4, 0.25);
        routingMatrix.addConnection(c2, c2, q3, q5, 0.25);

        routingMatrix.addConnection(c3, c3, q3, q1, 0.25);
        routingMatrix.addConnection(c3, c3, q3, q2, 0.25);
        routingMatrix.addConnection(c3, c3, q3, q4, 0.25);
        routingMatrix.addConnection(c3, c3, q3, q5, 0.25);

        routingMatrix.addConnection(c1, c1, q4, q1, 0.25);
        routingMatrix.addConnection(c1, c1, q4, q2, 0.25);
        routingMatrix.addConnection(c1, c1, q4, q3, 0.25);
        routingMatrix.addConnection(c1, c1, q4, q5, 0.25);

        routingMatrix.addConnection(c2, c2, q4, q1, 0.25);
        routingMatrix.addConnection(c2, c2, q4, q2, 0.25);
        routingMatrix.addConnection(c2, c2, q4, q3, 0.25);
        routingMatrix.addConnection(c2, c2, q4, q5, 0.25);

        routingMatrix.addConnection(c3, c3, q4, q1, 0.25);
        routingMatrix.addConnection(c3, c3, q4, q2, 0.25);
        routingMatrix.addConnection(c3, c3, q4, q3, 0.25);
        routingMatrix.addConnection(c3, c3, q4, q5, 0.25);

        routingMatrix.addConnection(c1, c1, q5, q1, 0.25);
        routingMatrix.addConnection(c1, c1, q5, q2, 0.25);
        routingMatrix.addConnection(c1, c1, q5, q3, 0.25);
        routingMatrix.addConnection(c1, c1, q5, q4, 0.25);

        routingMatrix.addConnection(c2, c2, q5, q1, 0.25);
        routingMatrix.addConnection(c2, c2, q5, q2, 0.25);
        routingMatrix.addConnection(c2, c2, q5, q3, 0.25);
        routingMatrix.addConnection(c2, c2, q5, q4, 0.25);

        routingMatrix.addConnection(c3, c3, q5, q1, 0.25);
        routingMatrix.addConnection(c3, c3, q5, q2, 0.25);
        routingMatrix.addConnection(c3, c3, q5, q3, 0.25);
        routingMatrix.addConnection(c3, c3, q5, q4, 0.25);

        network.link(routingMatrix);
        return network;
    }

    public static Network small_example() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.PS);
        q1.setNumberOfServers(1);
        Queue q2 = new Queue(network, "q2", SchedStrategy.FCFS);
        q2.setNumberOfServers(1);

        ClosedClass c1 = new ClosedClass(network, "c1", 2, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 2, q1);

        q1.setService(c1, Exp.fitMean(1.3));
        q1.setService(c2, Exp.fitMean(1.0));

        q2.setService(c1, Exp.fitMean(1.2));
        q2.setService(c2, Exp.fitMean(1.0));

        RoutingMatrix routingMatrix = network.initRoutingMatrix();

        routingMatrix.addConnection(c1, c1, q1, q2, 1.0);
        routingMatrix.addConnection(c2, c2, q1, q2, 1.0);

        routingMatrix.addConnection(c1, c1, q2, q1, 1.0);
        routingMatrix.addConnection(c2, c2, q2, q1, 1.0);

        network.link(routingMatrix);
        return network;
    }
}
