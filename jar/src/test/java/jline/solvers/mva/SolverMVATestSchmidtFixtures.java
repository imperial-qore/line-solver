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
 * Test fixtures for the Schmidt MVA algorithm.
 */
public class SolverMVATestSchmidtFixtures {

    /**
     * Network from FES model
     *
     * @return Network
     */
    public static Network test_textbook_model1() {
        Network network = new Network("test_textbook_model1");

        Queue q1 = new Queue(network, "q1");
        Queue q2 = new Queue(network, "q2");
        Queue q3 = new Queue(network, "q3");
        Queue q4 = new Queue(network, "q4");

        ClosedClass c1 = new ClosedClass(network, "c1", 2, q1);

        q1.setService(c1, Exp.fitRate(1.0));
        q2.setService(c1, Exp.fitRate(2.0));
        q3.setService(c1, Exp.fitRate(3.0));
        q4.setService(c1, Exp.fitRate(4.0));

        RoutingMatrix routingMatrix = new RoutingMatrix(network, List.of(c1), List.of(q1, q2, q3, q4));

        routingMatrix.addConnection(c1, c1, q1, q2, 0.5);
        routingMatrix.addConnection(c1, c1, q1, q3, 0.5);

        routingMatrix.addConnection(c1, c1, q2, q1, 0.5);
        routingMatrix.addConnection(c1, c1, q2, q4, 0.5);

        routingMatrix.addConnection(c1, c1, q3, q1, 0.5);
        routingMatrix.addConnection(c1, c1, q3, q4, 0.5);

        routingMatrix.addConnection(c1, c1, q4, q2, 0.4);
        routingMatrix.addConnection(c1, c1, q4, q3, 0.4);

        routingMatrix.addConnection(c1, c1, q4, q4, 0.2);

        network.link(routingMatrix);

        return network;

    }

    public static Network fcfs_varies() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.FCFS);
        q1.setNumberOfServers(2);
        Queue q2 = new Queue(network, "q2", SchedStrategy.FCFS);
        q2.setNumberOfServers(2);

        ClosedClass c1 = new ClosedClass(network, "c1", 5, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 5, q1);

        q1.setService(c1, Exp.fitMean(1.1));
        q1.setService(c2, Exp.fitMean(1.0));

        q2.setService(c1, Exp.fitMean(1.0));
        q2.setService(c2, Exp.fitMean(1.0));


        RoutingMatrix routingMatrix = new RoutingMatrix(network, List.of(c1, c2), List.of(q1, q2));

        routingMatrix.addConnection(c1, c1, q1, q2, 1.0);
        routingMatrix.addConnection(c2, c2, q1, q2, 1.0);

        routingMatrix.addConnection(c1, c1, q2, q1, 1.0);
        routingMatrix.addConnection(c2, c2, q2, q1, 1.0);

        network.link(routingMatrix);
        return network;
    }

    public static Network fcfs_varies_3_classes() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.FCFS);
        q1.setNumberOfServers(2);
        Queue q2 = new Queue(network, "q2", SchedStrategy.PS);
        q2.setNumberOfServers(2);

        ClosedClass c1 = new ClosedClass(network, "c1", 4, q2);
        ClosedClass c2 = new ClosedClass(network, "c2", 4, q2);
        ClosedClass c3 = new ClosedClass(network, "c3", 4, q2);

        q1.setService(c1, Exp.fitMean(16.1));
        q1.setService(c2, Exp.fitMean(8.0));
        q1.setService(c3, Exp.fitMean(0.9));

        q2.setService(c1, Exp.fitMean(1.0));
        q2.setService(c2, Exp.fitMean(1.0));
        q2.setService(c3, Exp.fitMean(1.0));


        RoutingMatrix routingMatrix = new RoutingMatrix(network, List.of(c1, c2, c3), List.of(q1, q2));

        routingMatrix.addConnection(c1, c1, q1, q2, 1.0);
        routingMatrix.addConnection(c2, c2, q1, q2, 1.0);
        routingMatrix.addConnection(c3, c3, q1, q2, 1.0);

        routingMatrix.addConnection(c1, c1, q2, q1, 1.0);
        routingMatrix.addConnection(c2, c2, q2, q1, 1.0);
        routingMatrix.addConnection(c3, c3, q2, q1, 1.0);

        network.link(routingMatrix);
        return network;
    }

    public static Network fcfs_varies_1_node() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.FCFS);
        q1.setNumberOfServers(2);

        ClosedClass c1 = new ClosedClass(network, "c1", 3, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 3, q1);
        ClosedClass c3 = new ClosedClass(network, "c3", 3, q1);

        q1.setService(c1, Exp.fitMean(1.1));
        q1.setService(c2, Exp.fitMean(1.0));
        q1.setService(c3, Exp.fitMean(0.9));

        RoutingMatrix routingMatrix = new RoutingMatrix(network, List.of(c1, c2, c3), List.of(q1));

        routingMatrix.addConnection(c1, c1, q1, q1, 1.0);
        routingMatrix.addConnection(c2, c2, q1, q1, 1.0);
        routingMatrix.addConnection(c3, c3, q1, q1, 1.0);

        network.link(routingMatrix);
        return network;
    }


    public static Network fcfs_1_class_1_node() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.FCFS);
        q1.setNumberOfServers(2);

        ClosedClass c1 = new ClosedClass(network, "c1", 20, q1);

        q1.setService(c1, Exp.fitMean(1.0));

        RoutingMatrix routingMatrix = new RoutingMatrix(network, List.of(c1), List.of(q1));

        routingMatrix.addConnection(c1, c1, q1, q1, 1.0);

        network.link(routingMatrix);
        return network;
    }

    public static Network fcfs_2_class_2_nodes() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.PS);
        q1.setNumberOfServers(2);
        Queue q2 = new Queue(network, "q2", SchedStrategy.FCFS);

        ClosedClass c1 = new ClosedClass(network, "c1", 4, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 4, q1);


        q1.setService(c1, Exp.fitMean(1.0));
        q1.setService(c2, Exp.fitMean(1.0));

        q2.setService(c1, Exp.fitMean(1.5));
        q2.setService(c2, Exp.fitMean(0.9));

        RoutingMatrix routingMatrix = new RoutingMatrix(network, List.of(c1, c2), List.of(q1, q2));

        routingMatrix.addConnection(c1, c1, q1, q1, 0.0);
        routingMatrix.addConnection(c2, c2, q1, q2, 1.0);

        routingMatrix.addConnection(c1, c1, q2, q1, 1.0);
        routingMatrix.addConnection(c1, c1, q2, q2, 0.0);

        network.link(routingMatrix);
        return network;
    }

    public static Network fcfs_2_class_1_node() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.FCFS);
        q1.setNumberOfServers(2);

        ClosedClass c1 = new ClosedClass(network, "c1", 6, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 6, q1);

        q1.setService(c1, Exp.fitMean(1.1));
        q1.setService(c2, Exp.fitMean(1.0));

        RoutingMatrix routingMatrix = new RoutingMatrix(network, List.of(c1, c2), List.of(q1));

        routingMatrix.addConnection(c1, c1, q1, q1, 1.0);
        routingMatrix.addConnection(c2, c2, q1, q1, 1.0);

        network.link(routingMatrix);
        return network;

    }

    public static Network fcfs_3_nodes_3_classes_varied_visit_ratio() {
        Network network = new Network("network");
        Queue q1 = new Queue(network, "q1", SchedStrategy.PS);
        q1.setNumberOfServers(2);
        Queue q2 = new Queue(network, "q2", SchedStrategy.FCFS);
        q2.setNumberOfServers(2);
        Queue q3 = new Queue(network, "q3", SchedStrategy.FCFS);
        q3.setNumberOfServers(3);

        ClosedClass c1 = new ClosedClass(network, "c1", 2, q1);
        ClosedClass c2 = new ClosedClass(network, "c2", 3, q1);
        ClosedClass c3 = new ClosedClass(network, "c3", 2, q1);

        q1.setService(c1, Exp.fitMean(30));
        q1.setService(c2, Exp.fitMean(30));
        q1.setService(c3, Exp.fitMean(30));

        q2.setService(c1, Exp.fitMean(114));
        q2.setService(c2, Exp.fitMean(117));
        q2.setService(c3, Exp.fitMean(56));

        q3.setService(c1, Exp.fitMean(50));
        q3.setService(c2, Exp.fitMean(37));
        q3.setService(c3, Exp.fitMean(119));

        RoutingMatrix routingMatrix = network.initRoutingMatrix();

        routingMatrix.addConnection(c1, c1, q1, q2, 0.3);
        routingMatrix.addConnection(c1, c1, q1, q3, 0.7);

        routingMatrix.addConnection(c2, c2, q1, q2, 0.5);
        routingMatrix.addConnection(c2, c2, q1, q3, 0.5);

        routingMatrix.addConnection(c3, c3, q1, q2, 0.9);
        routingMatrix.addConnection(c3, c3, q1, q3, 0.1);


        routingMatrix.addConnection(c1, c1, q2, q1, 0.4);
        routingMatrix.addConnection(c1, c1, q2, q3, 0.6);

        routingMatrix.addConnection(c2, c2, q2, q1, 0.7);
        routingMatrix.addConnection(c2, c2, q2, q3, 0.3);

        routingMatrix.addConnection(c3, c3, q2, q1, 0.5);
        routingMatrix.addConnection(c3, c3, q2, q3, 0.5);


        routingMatrix.addConnection(c1, c1, q3, q1, 1);
        routingMatrix.addConnection(c1, c1, q3, q2, 0);

        routingMatrix.addConnection(c2, c2, q3, q1, 0.5);
        routingMatrix.addConnection(c2, c2, q3, q2, 0.5);

        routingMatrix.addConnection(c3, c3, q3, q1, 0.5);
        routingMatrix.addConnection(c3, c3, q3, q2, 0.5);

        network.link(routingMatrix);
        return network;
    }

}
