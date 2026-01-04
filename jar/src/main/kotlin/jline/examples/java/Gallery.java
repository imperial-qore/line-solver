/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.models;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.VerboseLevel;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.solvers.jmt.JMT;
import org.apache.commons.math3.util.FastMath;

import javax.xml.parsers.ParserConfigurationException;
import java.net.URI;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static jline.util.Maths.linSpace;

/**
 * Gallery of simple and classical queueing models
 */
public class Gallery {
    /**
     * APH/M/1 queue with Acyclic Phase-Type arrivals.
     * <p>
     * Features:
     * - APH arrival process with high variability (SCV ≈ 1.999)
     * - FCFS queue with exponential service (rate 2.0)
     * - Demonstrates acyclic phase-type modeling
     *
     * @return configured APH/M/1 network model
     */
    public static Network gallery_aphm1() {
        Network model = new Network("APH/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, APH.fitCentral(1, 0.99, 1.999));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /**
     * Cox/M/1 queue with Coxian arrivals.
     * <p>
     * Features:
     * - Coxian arrival process with high variability (SCV ≈ 1.999)
     * - FCFS queue with exponential service (rate 2.0)
     * - Demonstrates Coxian distribution modeling
     *
     * @return configured Cox/M/1 network model
     */
    public static Network gallery_coxm1() {
        Network model = new Network("Cox/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Coxian.fitCentral(1, 0.99, 1.999));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /**
     * Default closed queueing network with 2 stations.
     *
     * @return closed queueing network with M=2 stations, no delay
     */
    public static Network gallery_cqn() {
        return gallery_cqn(2, false, 2300);
    }

    /**
     * Closed queueing network with specified number of stations.
     *
     * @param M number of stations
     * @return closed queueing network with M stations, no delay
     */
    public static Network gallery_cqn(int M) {
        return gallery_cqn(M, false, 2300);
    }

    /**
     * Closed queueing network with optional delay station.
     *
     * @param M        number of PS queue stations
     * @param useDelay whether to include an additional delay station
     * @return closed queueing network with M+1 stations if delay is used
     */
    public static Network gallery_cqn(int M, boolean useDelay) {
        return gallery_cqn(M, useDelay, 2300);
    }

    /**
     * Parameterized closed queueing network generator.
     * <p>
     * Features:
     * - M PS queue stations with random service rates
     * - Optional additional delay station with fixed rate
     * - Single closed class with random population (3 to 10*M+3 jobs)
     * - Serial routing through all stations
     * - Configurable random seed for reproducibility
     *
     * @param M        number of PS queue stations
     * @param useDelay whether to include an additional delay station
     * @param seed     random seed for reproducible parameter generation
     * @return configured closed queueing network
     */
    public static Network gallery_cqn(int M, boolean useDelay, long seed) {
        Network model = new Network("model");

        Random random = new Random();
        random.setSeed(seed);

        ServiceStation[] station = useDelay ? new ServiceStation[M + 1] : new ServiceStation[M];

        for (int i = 0; i < M; i++) {
            station[i] = new Queue(model, "Queue" + i, SchedStrategy.PS);
        }

        if (useDelay) {
            station[M] = new Delay(model, "Delay1");
        }
        JobClass jobclass = new ClosedClass(model, "Class1", FastMath.round(random.nextDouble() * 10 * M + 3), station[0], 0);

        for (int i = 0; i < M; i++) {
            station[i].setService(jobclass, Exp.fitMean(random.nextDouble() + i));
        }

        if (useDelay) {
            station[M].setService(jobclass, Exp.fitMean(2.0));
        }

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, Network.serialRouting(station));
        model.link(P);

        return model;
    }

    public static Network gallery_cqn_multiclass() {
        return gallery_cqn_multiclass(2, 2, false, 2300);
    }

    public static Network gallery_cqn_multiclass(int M, int R) {
        return gallery_cqn_multiclass(M, R, false, 2300);
    }

    public static Network gallery_cqn_multiclass(int M, int R, boolean useDelay) {
        return gallery_cqn_multiclass(M, R, useDelay, 2300);
    }

    public static Network gallery_cqn_multiclass(int M, int R, boolean useDelay, long seed) {
        Network model = new Network("model");

        Random random = new Random();
        random.setSeed(seed);

        ServiceStation[] station = useDelay ? new ServiceStation[M + 1] : new ServiceStation[M];

        for (int i = 0; i < M; i++) {
            station[i] = new Queue(model, "Queue" + i, SchedStrategy.PS);
        }

        if (useDelay) {
            station[M] = new Delay(model, "Delay1");
        }

        JobClass[] jobclass = new JobClass[R];

        for (int r = 0; r < R; r++) {
            jobclass[r] = new ClosedClass(model, "Class" + r, 5, station[0], 0);
        }

        for (int r = 0; r < R; r++) {
            for (int i = 0; i < M; i++) {
                station[i].setService(jobclass[r], Exp.fitMean(Math.round(random.nextDouble() * 50)));
            }

            if (useDelay) {
                station[M].setService(jobclass[r], Exp.fitMean(Math.round(random.nextDouble() * 100)));
            }
        }

        RoutingMatrix P = model.initRoutingMatrix();
        for (int r = 0; r < R; r++) {
            P.set(jobclass[r], Network.serialRouting(station));
        }
        model.link(P);

        return model;
    }

    public static Network gallery_dm1() {
        Network model = new Network("Det/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Det(1));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_erldk() {
        return gallery_erldk(2);
    }

    public static Network gallery_erldk(int k) {
        Network model = new Network("Erl/Det/k");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Erlang.fitMeanAndOrder(1.0, 5));
        queue.setService(oclass, new Det(2.0 / k));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_erlerl1_reentrant() {
        return gallery_erlerl1_reentrant(5);
    }

    public static Network gallery_erlerl1_reentrant(int n) {
        Network model = new Network("Erl/Erl/1-Reentrant");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, Erlang.fitMeanAndOrder(1, n)); // (Source,Class1)
        node1.setArrival(jobclass2, Disabled.getInstance()); // (Source,Class2)
        node2.setService(jobclass1, Erlang.fitMeanAndOrder(0.5, n)); // (Queue,Class1)
        node2.setService(jobclass2, new Exp(3)); // (Queue,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 1.00); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_erlm1() {
        Network model = new Network("Erl/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Erlang.fitMeanAndOrder(1.0, 5));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_erlm1_ps() {
        Network model = new Network("Erl/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.PS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Erlang.fitMeanAndOrder(1.0, 5));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_erlm1_reentrant() {
        Network model = new Network("Erl/M/1-Reentrant");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, Erlang.fitMeanAndOrder(1, 5)); // (Source,Class1)
        node1.setArrival(jobclass2, Disabled.getInstance()); // (Source,Class2)
        node2.setService(jobclass1, new Exp(2)); // (Queue,Class1)
        node2.setService(jobclass2, new Exp(3)); // (Queue,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 1.00); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_gamm1() {
        Network model = new Network("Gam/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Gamma.fitMeanAndSCV(1.0, 1.0 / 5));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_hyperl1_feedback() {
        Network model = new Network("Hyper/Erl/1-Feedback");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass1 = new OpenClass(model, "Class1");
        source.setArrival(oclass1, HyperExp.fitMeanAndSCV(1, 64));
        queue.setService(oclass1, Erlang.fitMeanAndOrder(0.5, 5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass1, oclass1, source, queue, 1.0);
        P.set(oclass1, oclass1, queue, queue, 0.9);
        P.set(oclass1, oclass1, queue, sink, 0.1);
        model.link(P);
        return model;
    }

    public static Network gallery_hyperl1_reentrant() {
        Network model = new Network("Hyper/Erl/1-Reentrant");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, HyperExp.fitMeanAndSCV(1, 64)); // (Source,Class1)
        node1.setArrival(jobclass2, Disabled.getInstance()); // (Source,Class2)
        node2.setService(jobclass1, Erlang.fitMeanAndOrder(0.5, 5)); // (Queue,Class1)
        node2.setService(jobclass2, new Exp(3)); // (Queue,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 1.00); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_hyperlk() {
        return gallery_hyperlk(2);
    }

    public static Network gallery_hyperlk(int k) {
        Network model = new Network("Hyper/Erl/k");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, HyperExp.fitMeanAndSCVBalanced(1 / 1.8, 4));
        queue.setService(oclass, Erlang.fitMeanAndSCV(1, 0.25));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_hyphyp1_linear() {
        return gallery_hyphyp1_linear(2, 0.9);
    }

    public static Network gallery_hyphyp1_linear(int n) {
        return gallery_hyphyp1_linear(n, 0.9);
    }

    public static Network gallery_hyphyp1_linear(int n, Double Umax) {
        Network model = new Network("Hyp/Hyp/1-Linear");

        // Block 1: nodes
        List<Node> nodes = new ArrayList<>();
        nodes.add(new Source(model, "mySource"));
        for (int i = 1; i <= n; i++) {
            nodes.add(new Queue(model, "Queue" + i, SchedStrategy.FCFS));
        }
        nodes.add(new Sink(model, "mySink"));

        // Block 2: classes
        OpenClass oclass = new OpenClass(model, "myClass");
        ((Source) nodes.get(0)).setArrival(oclass, HyperExp.fitMeanAndSCV(1, 2));

        double[] firstHalf = linSpace(0.1, Umax, n / 2);
        double[] means;

        if (n % 2 == 0) {
            // Even case
            means = new double[n];
            for (int i = 0; i < firstHalf.length; i++) {
                means[i] = firstHalf[i];
                means[n - 1 - i] = firstHalf[i];
            }
        } else {
            // Odd case
            means = new double[n];
            for (int i = 0; i < firstHalf.length; i++) {
                means[i] = firstHalf[i];
                means[n - 1 - i] = firstHalf[i];
            }
            means[firstHalf.length] = Umax;
        }

        for (int i = 1; i <= n; i++) {
            ((Queue) nodes.get(i)).setService(oclass, HyperExp.fitMeanAndSCV(means[i - 1], 1 + i));
        }

        // Block 3: topology
        model.link(Network.serialRouting(nodes));

        return model;
    }

    public static Network gallery_hyphyp1_reentrant() {
        Network model = new Network("Hyper/Hyper/1-Reentrant");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, HyperExp.fitMeanAndSCV(1, 64)); // (Source,Class1)
        node1.setArrival(jobclass2, Disabled.getInstance()); // (Source,Class2)
        node2.setService(jobclass1, HyperExp.fitMeanAndSCV(0.5, 4)); // (Queue,Class1)
        node2.setService(jobclass2, new Exp(3)); // (Queue,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 1.00); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_hyphyp1_tandem() {
        return gallery_hyphyp1_linear(2);
    }

    public static Network gallery_hypm1() {
        Network model = new Network("H2/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, HyperExp.fitMeanAndSCV(1, 64));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_hypm1_reentrant() {
        Network model = new Network("Hyper/Erl/1-Reentrant");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, HyperExp.fitMeanAndSCV(1, 4)); // (Source,Class1)
        node1.setArrival(jobclass2, Disabled.getInstance()); // (Source,Class2)
        node2.setService(jobclass1, new Exp(2)); // (Queue,Class1)
        node2.setService(jobclass2, new Exp(3)); // (Queue,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 1.00); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_mapm1() {
        return gallery_mapm1(MAP.rand(2));
    }

    public static Network gallery_mapm1(MAP map) {
        Network model = new Network("MAP/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, map);
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_mapmk() {
        return gallery_mapmk(MAP.rand(2), 2);
    }

    public static Network gallery_mapmk(MAP map) {
        return gallery_mapmk(map, 2);
    }

    public static Network gallery_mapmk(MAP map, int k) {
        Network model = new Network("MAP/M/k");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, map);
        queue.setService(oclass, new Exp(2));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_mdk() {
        return gallery_mdk(2);
    }

    public static Network gallery_mdk(int k) {
        Network model = new Network("M/D/k");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Exp.fitMean(1));
        queue.setService(oclass, new Det(2.0 / k));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_merl1() {
        Network model = new Network("M/E/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, Erlang.fitMeanAndOrder(0.5, 2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_merl1_linear() {
        return gallery_merl1_linear(2, 0.9);
    }

    public static Network gallery_merl1_linear(int n) {
        return gallery_merl1_linear(n, 0.9);
    }

    public static Network gallery_merl1_linear(int n, Double Umax) {
        Network model = new Network("M/Erl/1-Linear");

        // Block 1: nodes
        List<Node> nodes = new ArrayList<>();
        nodes.add(new Source(model, "mySource"));
        for (int i = 1; i <= n; i++) {
            nodes.add(new Queue(model, "Queue" + i, SchedStrategy.FCFS));
        }
        nodes.add(new Sink(model, "mySink"));

        // Block 2: classes
        OpenClass oclass = new OpenClass(model, "myClass");
        ((Source) nodes.get(0)).setArrival(oclass, new Exp(1));

        double[] firstHalf = linSpace(0.1, Umax, n / 2);
        double[] means;

        if (n % 2 == 0) {
            // Even case
            means = new double[n];
            for (int i = 0; i < firstHalf.length; i++) {
                means[i] = firstHalf[i];
                means[n - 1 - i] = firstHalf[i];
            }
        } else {
            // Odd case
            means = new double[n];
            for (int i = 0; i < firstHalf.length; i++) {
                means[i] = firstHalf[i];
                means[n - 1 - i] = firstHalf[i];
            }
            means[firstHalf.length] = Umax;
        }

        for (int i = 1; i <= n; i++) {
            ((Queue) nodes.get(i)).setService(oclass, Erlang.fitMeanAndOrder(means[i - 1], i));
        }

        // Block 3: topology
        model.link(Network.serialRouting(nodes));

        return model;
    }

    public static Network gallery_merl1_reentrant() {
        Network model = new Network("M/Erl/1-Reentrant");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, Exp.fitMean(1.00)); // (Source,Class1)
        node1.setArrival(jobclass2, Disabled.getInstance()); // (Source,Class2)
        node2.setService(jobclass1, Erlang.fitMeanAndOrder(0.5, 5)); // (Queue,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.333333)); // (Queue,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 1.00); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_merl1_tandem() {
        return gallery_merl1_linear(2);
    }

    public static Network gallery_merlk() {
        return gallery_merlk(2);
    }

    public static Network gallery_merlk(int k) {
        Network model = new Network("M/E/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, Erlang.fitMeanAndOrder(0.5, 2));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_mhyp1() {
        Network model = new Network("M/H/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, HyperExp.fitMeanAndSCV(0.5, 4));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_mhyp1_linear() {
        return gallery_mhyp1_linear(2, 0.9);
    }

    public static Network gallery_mhyp1_linear(int n) {
        return gallery_mhyp1_linear(n, 0.9);
    }

    public static Network gallery_mhyp1_linear(int n, Double Umax) {
        Network model = new Network("M/Hyp/1-Linear");

        // Block 1: nodes
        List<Node> nodes = new ArrayList<>();
        nodes.add(new Source(model, "mySource"));
        for (int i = 1; i <= n; i++) {
            nodes.add(new Queue(model, "Queue" + i, SchedStrategy.FCFS));
        }
        nodes.add(new Sink(model, "mySink"));

        // Block 2: classes
        OpenClass oclass = new OpenClass(model, "myClass");
        ((Source) nodes.get(0)).setArrival(oclass, new Exp(1));

        double[] firstHalf = linSpace(0.1, Umax, n / 2);
        double[] means;

        if (n % 2 == 0) {
            // Even case
            means = new double[n];
            for (int i = 0; i < firstHalf.length; i++) {
                means[i] = firstHalf[i];
                means[n - 1 - i] = firstHalf[i];
            }
        } else {
            // Odd case
            means = new double[n];
            for (int i = 0; i < firstHalf.length; i++) {
                means[i] = firstHalf[i];
                means[n - 1 - i] = firstHalf[i];
            }
            means[firstHalf.length] = Umax;
        }

        for (int i = 1; i <= n; i++) {
            ((Queue) nodes.get(i)).setService(oclass, HyperExp.fitMeanAndSCV(means[i - 1], n));
        }

        // Block 3: topology
        model.link(Network.serialRouting(nodes));

        return model;
    }

    public static Network gallery_mhyp1_reentrant() {
        Network model = new Network("M/Hyper/1-Reentrant");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, Exp.fitMean(1.00)); // (Source,Class1)
        node1.setArrival(jobclass2, Disabled.getInstance()); // (Source,Class2)
        node2.setService(jobclass1, Coxian.fitMeanAndSCV(0.5, 4)); // (Queue,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.333333)); // (Queue,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 1.00); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_mhyp1_tandem() {
        return gallery_mhyp1_linear(2);
    }

    public static Network gallery_mhypk() {
        return gallery_mhypk(2);
    }

    public static Network gallery_mhypk(int k) {
        Network model = new Network("M/H/k");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, HyperExp.fitMeanAndSCV(0.5, 4));
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /**
     * Classic M/M/1 queue.
     * <p>
     * Features:
     * - Exponential arrivals (rate 1.0) and service (rate 2.0)
     * - FCFS scheduling
     * - Single server queue
     * - Fundamental queueing model
     *
     * @return configured M/M/1 network model
     */
    public static Network gallery_mm1() {
        Network model = new Network("M/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_mm1_feedback() {
        return gallery_mm1_feedback(1.0 / 3);
    }

    /**
     * M/M/1 queue with feedback routing.
     * <p>
     * Features:
     * - Jobs can return to the queue with probability p
     * - Exit to sink with probability (1-p)
     * - Demonstrates feedback queueing systems
     * - Higher effective service rate due to feedback
     *
     * @param p feedback probability (jobs returning to queue)
     * @return configured M/M/1 feedback network model
     */
    public static Network gallery_mm1_feedback(double p) {
        Network model = new Network("M/M/1-Feedback");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass1 = new OpenClass(model, "Class1");
        source.setArrival(oclass1, Exp.fitMean(1));
        queue.setService(oclass1, Exp.fitMean(0.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass1, oclass1, source, queue, 1.0);
        P.set(oclass1, oclass1, queue, queue, p);
        P.set(oclass1, oclass1, queue, sink, 1 - p);
        model.link(P);
        return model;
    }

    public static Network gallery_mm1_linear() {
        return gallery_mm1_linear(2, 0.9);
    }

    public static Network gallery_mm1_linear(Integer n) {
        return gallery_mm1_linear(n, 0.9);
    }

    public static Network gallery_mm1_linear(Integer n, Double Umax) {
        Network model = new Network("M/M/1-Linear");

        // Block 1: nodes
        List<Node> nodes = new ArrayList<>();
        nodes.add(new Source(model, "mySource"));
        for (int i = 1; i <= n; i++) {
            nodes.add(new Queue(model, "Queue" + i, SchedStrategy.FCFS));
        }
        nodes.add(new Sink(model, "mySink"));

        // Block 2: classes
        OpenClass oclass = new OpenClass(model, "myClass");
        ((Source) nodes.get(0)).setArrival(oclass, new Exp(1));

        double[] firstHalf = linSpace(0.1, Umax, n / 2);
        double[] means;

        if (n % 2 == 0) {
            // Even case
            means = new double[n];
            for (int i = 0; i < firstHalf.length; i++) {
                means[i] = firstHalf[i];
                means[n - 1 - i] = firstHalf[i];
            }
        } else {
            // Odd case
            means = new double[n];
            for (int i = 0; i < firstHalf.length; i++) {
                means[i] = firstHalf[i];
                means[n - 1 - i] = firstHalf[i];
            }
            means[firstHalf.length] = Umax;
        }


        for (int i = 1; i <= n; i++) {
            ((Queue) nodes.get(i)).setService(oclass, Exp.fitMean(means[i - 1]));
        }

        // Block 3: topology
        model.link(Network.serialRouting(nodes));

        return model;
    }

    public static Network gallery_mm1_multiclass() {
        Network model = new Network("M[2]/M[2]/1");

        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");

        OpenClass oclass1 = new OpenClass(model, "myClass1");
        source.setArrival(oclass1, new Exp(1));
        queue.setService(oclass1, new Exp(4));

        OpenClass oclass2 = new OpenClass(model, "myClass2");
        source.setArrival(oclass2, new Exp(0.5));
        queue.setService(oclass2, new Exp(4));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass1, Network.serialRouting(source, queue, sink));
        P.set(oclass2, Network.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    public static Network gallery_mm1_prio() {
        Network model = new Network("M[2]/M[2]/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.HOL);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass1 = new OpenClass(model, "myClass1", 1);
        source.setArrival(oclass1, new Exp(1));
        queue.setService(oclass1, new Exp(4));
        OpenClass oclass2 = new OpenClass(model, "myClass2", 0);
        source.setArrival(oclass2, new Exp(0.5));
        queue.setService(oclass2, new Exp(4));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass1, Network.serialRouting(source, queue, sink));
        P.set(oclass2, Network.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    public static Network gallery_mm1_ps() {
        Network model = new Network("M/M/1-PS");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.PS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_mm1_ps_feedback() {
        return gallery_mm1_ps_feedback(1.0 / 3.0);
    }

    public static Network gallery_mm1_ps_feedback(double p) {
        Network model = new Network("M/M/1-PS-Feedback");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.PS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Exp.fitMean(1));
        queue.setService(oclass, Exp.fitMean(0.5));
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(oclass, oclass, source, queue, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(oclass, oclass, queue, queue, p); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(oclass, oclass, queue, sink, 1 - p); // (Queue,Class2) -> (Sink,Class2)

        model.link(routingMatrix);
        return model;
    }

    public static Network gallery_mm1_ps_multiclass() {
        Network model = new Network("M[2]/M[2]/1");

        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.PS);
        Sink sink = new Sink(model, "mySink");

        OpenClass oclass1 = new OpenClass(model, "myClass1");
        source.setArrival(oclass1, new Exp(1));
        queue.setService(oclass1, new Exp(4));

        OpenClass oclass2 = new OpenClass(model, "myClass2");
        source.setArrival(oclass2, new Exp(0.5));
        queue.setService(oclass2, new Exp(4));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass1, Network.serialRouting(source, queue, sink));
        P.set(oclass2, Network.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    public static Network gallery_mm1_ps_reentrant() {
        Network model = new Network("M/M/1-PS-Reentrant");

        // Block 1: nodes
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.PS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, Exp.fitMean(1.00)); // (Source,Class1)
        node1.setArrival(jobclass2, Disabled.getInstance()); // (Source,Class2)
        node2.setService(jobclass1, Exp.fitMean(0.50)); // (Queue,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.333333)); // (Queue,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 1.00); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_mm1_reentrant() {
        Network model = new Network("M/M/1-Reentrant");

        // Block 1: nodes			
        Source node1 = new Source(model, "Source");
        Queue node2 = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink node3 = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
        OpenClass jobclass2 = new OpenClass(model, "Class2", 0);

        node1.setArrival(jobclass1, Exp.fitMean(1.00)); // (Source,Class1)
        node1.setArrival(jobclass2, Disabled.getInstance()); // (Source,Class2)
        node2.setService(jobclass1, Exp.fitMean(0.50)); // (Queue,Class1)
        node2.setService(jobclass2, Exp.fitMean(0.333333)); // (Queue,Class2)

        // Block 3: topology	
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 1.00); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue,Class2) -> (Sink,Class2)

        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_mm1_tandem() {
        return gallery_mm1_linear(2, 0.9);
    }

    public static Network gallery_mm1_tandem(Double Umax) {
        return gallery_mm1_linear(2, Umax);
    }

    public static Network gallery_mm1_tandem_multiclass() {
        Network model = new Network("M[2]/M[2]/1 -> -/M[2]/1");

        Source source = new Source(model, "mySource");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");

        OpenClass oclass1 = new OpenClass(model, "myClass1");
        source.setArrival(oclass1, new Exp(1));
        queue1.setService(oclass1, new Exp(4));
        queue1.setService(oclass1, new Exp(6));

        OpenClass oclass2 = new OpenClass(model, "myClass2");
        source.setArrival(oclass2, new Exp(0.5));
        queue2.setService(oclass2, new Exp(2));
        queue2.setService(oclass2, new Exp(6));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass1, Network.serialRouting(source, queue1, queue2, sink));
        P.set(oclass2, Network.serialRouting(source, queue1, queue2, sink));
        model.link(P);

        return model;
    }

    public static Network gallery_mmap1() {
        MAP map = MAP.rand();
        map.setMean(0.5);
        return gallery_mmap1(map);
    }

    public static Network gallery_mmap1(MAP map) {
        Network model = new Network("M/MAP/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.PS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, map);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_mmap1_multiclass() {
        MAP map1 = MAP.rand();
        map1.setMean(0.5);
        MAP map2 = MAP.rand();
        map2.setMean(0.5);
        return gallery_mmap1_multiclass(map1, map2);
    }

    public static Network gallery_mmap1_multiclass(MAP map1, MAP map2) {
        Network model = new Network("M/MAP/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.PS);
        Sink sink = new Sink(model, "mySink");

        OpenClass oclass1 = new OpenClass(model, "myClass1");
        source.setArrival(oclass1, new Exp(0.35 / map1.getMean()));
        queue.setService(oclass1, map1);
        OpenClass oclass2 = new OpenClass(model, "myClass2");
        source.setArrival(oclass2, new Exp(0.15 / map2.getMean()));
        queue.setService(oclass2, map2);

        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_mmapk() {
        return gallery_mmapk(MAP.rand(), 2);
    }

    public static Network gallery_mmapk(MAP map, int k) {
        Network model = new Network("M/MAP/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.PS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, map);
        queue.setNumberOfServers(k);
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_mmk() {
        return gallery_mmk(2);

    }

    /**
     * M/M/k queue with multiple servers.
     * <p>
     * Features:
     * - Exponential arrivals (rate 1.0) and service (rate 0.5 per server)
     * - k identical servers with FCFS scheduling
     * - Demonstrates multi-server queueing
     * - Higher capacity than single server systems
     *
     * @param k number of servers
     * @return configured M/M/k network model
     */
    public static Network gallery_mmk(int k) {
        Network model = new Network("M/M/k");

        // Block 1: nodes
        Source node1 = new Source(model, "mySource");
        Queue node2 = new Queue(model, "myQueue", SchedStrategy.FCFS);
        node2.setNumberOfServers(k);
        Sink node3 = new Sink(model, "mySink");

        // Block 2: classes
        OpenClass jobclass1 = new OpenClass(model, "myClass", 0);

        node1.setArrival(jobclass1, Exp.fitMean(1.00)); // (mySource,myClass)
        node2.setService(jobclass1, Exp.fitMean(0.500000)); // (myQueue,myClass)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (mySource,myClass) -> (myQueue,myClass)
        routingMatrix.set(jobclass1, jobclass1, node2, node3, 1.00); // (myQueue,myClass) -> (mySink,myClass)

        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_mpar1() {
        Network model = new Network("Par/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, Pareto.fitMeanAndSCV(0.5, 64));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_parm1() {
        Network model = new Network("Par/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Pareto.fitMeanAndSCV(1, 64));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_repairmen() {
        return gallery_repairmen(2, 23000);
    }

    public static Network gallery_repairmen(int seed) {
        return gallery_repairmen(1, seed);
    }


    /**
     * Machine repair model (repairmen problem).
     * <p>
     * Features:
     * - Closed network representing machines and repair facility
     * - nServers repair servers with PS scheduling
     * - Random number of machines (3 to 13)
     * - Working state modeled as delay station
     * - Random service rates for realistic modeling
     *
     * @param nServers number of repair servers
     * @param seed     random seed for reproducible parameters
     * @return configured machine repair network model
     */
    public static Network gallery_repairmen(int nServers, long seed) {
        Network model = new Network("model");

        Random random = new Random();
        random.setSeed(seed);

        ServiceStation[] station = new ServiceStation[2];
        station[0] = new Queue(model, "Queue0", SchedStrategy.PS);
        station[0].setNumberOfServers(nServers);
        station[1] = new Delay(model, "Delay1");

        JobClass jobclass = new ClosedClass(model, "Class1", FastMath.round(random.nextDouble() * 10 + 3), station[0], 0);

        station[0].setService(jobclass, Exp.fitMean(random.nextDouble()));
        station[1].setService(jobclass, Exp.fitMean(2.0));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, Network.serialRouting(station));
        model.link(P);

        return model;
    }

    public static Network gallery_replayerm1() {
        URI fileURI = null;
        try {
            fileURI = Gallery.class.getResource("/example_trace.txt").toURI();
        } catch (Exception e) {
            e.printStackTrace();
        }
        assert fileURI != null;
        String fileName = Paths.get(fileURI).toString();
        return gallery_replayerm1(fileName);
    }

    public static Network gallery_replayerm1(String fileName) {
        Network model = new Network("Trace/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        Replayer replayer = new Replayer(fileName);
        source.setArrival(oclass, replayer);
        queue.setService(oclass, new Exp(3 / replayer.getMean()));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_um1() {
        Network model = new Network("U/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Uniform(1, 2));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /**
     * Main method for testing and demonstrating gallery examples.
     *
     * <p>Currently configured to:
     * - Run gallery_replayerm1() with trace-driven arrivals
     * - Solve using JMT solver with debug verbosity
     * - Display model visualization
     *
     * @param args command line arguments (not used)
     * @throws IllegalAccessException       if model access is restricted
     * @throws ParserConfigurationException if XML parsing fails
     */
    public static void main(String[] args) throws IllegalAccessException, ParserConfigurationException {
        Network model = gallery_replayerm1();
        new JMT(model, "verbose", VerboseLevel.DEBUG).getAvgTable().print();
        model.view();
    }
}
