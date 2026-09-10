/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.models;

import jline.lang.*;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.layered.*;
import jline.gen.NetworkGenerator;
import jline.gen.LayeredNetworkGenerator;
import jline.util.matrix.Matrix;
import jline.VerboseLevel;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.solvers.wrappers.jmt.JMT;
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
        Network model = new Network("D/M/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
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
        node2.setService(jobclass1, Erlang.fitMeanAndOrder(0.1, n)); // (Queue,Class1)
        node2.setService(jobclass2, new Exp(10)); // (Queue,Class2)

        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00); // (Source,Class1) -> (Queue,Class1)
        routingMatrix.set(jobclass1, jobclass2, node2, node2, 0.50); // (CS_Queue_to_Queue,Class1) -> (Queue,Class2)
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 0.50); // (Queue,Class2) -> (Sink,Class2)

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
        queue.setService(oclass1, Erlang.fitMeanAndOrder(0.05, 5));
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
        node2.setNumberOfServers(2);
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
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
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

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");

        OpenClass oclass1 = new OpenClass(model, "myClass1");
        source.setArrival(oclass1, new Exp(1));
        queue1.setService(oclass1, new Exp(4));
        queue2.setService(oclass1, new Exp(6));

        OpenClass oclass2 = new OpenClass(model, "myClass2");
        source.setArrival(oclass2, new Exp(0.5));
        queue1.setService(oclass2, new Exp(2));
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
        Network model = new Network("M/Par/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(1));
        queue.setService(oclass, Pareto.fitMeanAndSCV(0.5, 64));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_parm1() {
        Network model = new Network("Par/M/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
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
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
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

    public static Network gallery_detm1() {
        Network model = new Network("D/M/1");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, new Det(1));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_erlerl1() {
        return gallery_erlerl1(5);
    }

    public static Network gallery_erlerl1(int n) {
        Network model = new Network("Erl/Erl/1");

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
        routingMatrix.set(jobclass2, jobclass2, node2, node3, 1.00); // (Queue,Class2) -> (Sink,Class2)
        model.link(routingMatrix);

        return model;
    }

    public static Network gallery_erlm1ps() {
        Network model = new Network("Er/M/1-PS");
        Source source = new Source(model, "mySource");
        Queue queue = new Queue(model, "myQueue", SchedStrategy.PS);
        Sink sink = new Sink(model, "mySink");
        OpenClass oclass = new OpenClass(model, "myClass");
        source.setArrival(oclass, Erlang.fitMeanAndOrder(1.0, 5));
        queue.setService(oclass, new Exp(2));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    public static Network gallery_lukumar_reentrant() {
        return gallery_lukumar_reentrant("FCFS");
    }

    public static Network gallery_lukumar_reentrant(String schedStrategy) {
        SchedStrategy sched1, sched2;
        switch (schedStrategy.toUpperCase()) {
            case "HOL":
                sched1 = SchedStrategy.HOL;
                sched2 = SchedStrategy.HOL;
                break;
            case "PS":
                sched1 = SchedStrategy.PS;
                sched2 = SchedStrategy.PS;
                break;
            default:
                sched1 = SchedStrategy.FCFS;
                sched2 = SchedStrategy.FCFS;
                break;
        }

        Network model = new Network("Lu-Kumar-Reentrant");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue station1 = new Queue(model, "Station1", sched1);
        Queue station2 = new Queue(model, "Station2", sched2);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes (priority is the third constructor argument)
        OpenClass class1 = new OpenClass(model, "Class1", 1);
        OpenClass class2 = new OpenClass(model, "Class2", 0);
        OpenClass class3 = new OpenClass(model, "Class3", 1);
        OpenClass class4 = new OpenClass(model, "Class4", 0);

        // Block 3: arrivals
        double arrivalRate = 0.08;
        source.setArrival(class1, new Exp(arrivalRate));
        source.setArrival(class2, Disabled.getInstance());
        source.setArrival(class3, new Exp(arrivalRate));
        source.setArrival(class4, Disabled.getInstance());

        // Block 4: service times (asymmetric Kumar-Seidman configuration)
        double m1 = 10.0, m2 = 1.0, m3 = 10.0, m4 = 1.0;
        station1.setService(class1, new Exp(1 / m1));
        station1.setService(class2, Disabled.getInstance());
        station1.setService(class3, Disabled.getInstance());
        station1.setService(class4, new Exp(1 / m4));
        station2.setService(class1, Disabled.getInstance());
        station2.setService(class2, new Exp(1 / m2));
        station2.setService(class3, new Exp(1 / m3));
        station2.setService(class4, Disabled.getInstance());

        // Block 5: routing (two-chain topology)
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, source, station1, 1); // Source -> Class1@Station1
        P.set(class1, class2, station1, station2, 1); // Class1@Station1 -> Class2@Station2
        P.set(class2, class2, station2, sink, 1); // Class2@Station2 -> Sink
        P.set(class3, class3, source, station2, 1); // Source -> Class3@Station2
        P.set(class3, class4, station2, station1, 1); // Class3@Station2 -> Class4@Station1
        P.set(class4, class4, station1, sink, 1); // Class4@Station1 -> Sink
        model.link(P);

        return model;
    }

    public static LayeredNetwork gallery_multitier() {
        LayeredNetwork model = new LayeredNetwork("testLQN3");

        // Layer 1: client
        Processor P0 = new Processor(model, "P0", 1, SchedStrategy.PS);
        Task T0 = new Task(model, "T0", 1, SchedStrategy.REF).on(P0);
        Entry E0 = new Entry(model, "E0").on(T0);

        // Layer 2: application server
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", 1, SchedStrategy.FCFS).on(P1);
        Entry E10 = new Entry(model, "E10").on(T1);
        Entry E11 = new Entry(model, "E11").on(T1);
        Entry E12 = new Entry(model, "E12").on(T1);
        Entry E13 = new Entry(model, "E13").on(T1);

        // Layer 3: database
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2);
        Entry E20 = new Entry(model, "E20").on(T2);
        Entry E21 = new Entry(model, "E21").on(T2);
        Entry E22 = new Entry(model, "E22").on(T2);
        Entry E23 = new Entry(model, "E23").on(T2);

        // Client activities
        Activity A0 = new Activity(model, "A0", new Exp(1.0)).on(T0).boundTo(E0).synchCall(E12, 1.0);
        Activity A1 = new Activity(model, "A1", new Exp(1.0)).on(T0).synchCall(E10, 1.0);
        Activity A2 = new Activity(model, "A2", new Exp(1.0)).on(T0).synchCall(E11, 1.0);
        Activity A3 = new Activity(model, "A3", new Exp(1.0)).on(T0).synchCall(E13, 1.0);

        // Application activities
        Activity B0 = new Activity(model, "B0", new Exp(1.0)).on(T1).boundTo(E10);
        Activity B1 = new Activity(model, "B1", new Exp(1.0)).on(T1).repliesTo(E10);
        Activity B2 = new Activity(model, "B2", new Exp(1.0)).on(T1).boundTo(E11);
        Activity B3 = new Activity(model, "B3", new Exp(1.0)).on(T1).synchCall(E21, 1.0).repliesTo(E11);
        Activity B4 = new Activity(model, "B4", new Exp(1.0)).on(T1).boundTo(E12).synchCall(E20, 1.0).repliesTo(E12);
        Activity B5 = new Activity(model, "B5", new Exp(1.0)).on(T1).boundTo(E13);
        Activity B6 = new Activity(model, "B6", new Exp(1.0)).on(T1);
        Activity B7 = new Activity(model, "B7", new Exp(1.0)).on(T1).synchCall(E22, 1.0);
        Activity B7a = new Activity(model, "B7a", new Exp(1.0)).on(T1);
        Activity B7b = new Activity(model, "B7b", new Exp(1.0)).on(T1).synchCall(E23, 1.0);
        Activity B8 = new Activity(model, "B8", new Exp(1.0)).on(T1).repliesTo(E13);

        // Database activities
        Activity C0 = new Activity(model, "C0", new Exp(1.0)).on(T2).boundTo(E20);
        Activity C1 = new Activity(model, "C1", new Exp(1.0)).on(T2).repliesTo(E20);
        Activity C2 = new Activity(model, "C2", new Exp(1.0)).on(T2).boundTo(E21).repliesTo(E21);
        Activity C3 = new Activity(model, "C3", new Exp(1.0)).on(T2).boundTo(E22);
        Activity C4 = new Activity(model, "C4", new Exp(1.0)).on(T2);
        Activity C5 = new Activity(model, "C5", new Exp(1.0)).on(T2).repliesTo(E22);
        Activity C6 = new Activity(model, "C6", new Exp(1.0)).on(T2).boundTo(E23).repliesTo(E23);

        // Precedences
        T0.addPrecedence(ActivityPrecedence.Serial(A0, A1, A2, A3));
        T1.addPrecedence(ActivityPrecedence.Serial(B0, B1));
        T1.addPrecedence(ActivityPrecedence.Serial(B2, B3));
        T1.addPrecedence(ActivityPrecedence.Serial(B5, B6, B7));
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.7);
        probs.set(0, 1, 0.3);
        List<Activity> b7targets = new ArrayList<Activity>();
        b7targets.add(B7a);
        b7targets.add(B7b);
        T1.addPrecedence(ActivityPrecedence.OrFork(B7, b7targets, probs));
        List<Activity> b7sources = new ArrayList<Activity>();
        b7sources.add(B7a);
        b7sources.add(B7b);
        T1.addPrecedence(ActivityPrecedence.OrJoin(b7sources, B8));
        T2.addPrecedence(ActivityPrecedence.Serial(C0, C1));
        T2.addPrecedence(ActivityPrecedence.Serial(C3, C4, C5));

        return model;
    }

    public static LayeredNetwork gallery_multitier_storage() {
        LayeredNetwork model = new LayeredNetwork("testLQN3_Cache");

        // Layer 1: client
        Processor P0 = new Processor(model, "P0", 1, SchedStrategy.PS);
        Task T0 = new Task(model, "T0", 1, SchedStrategy.REF).on(P0);
        Entry E0 = new Entry(model, "E0").on(T0);

        // Layer 2: application server
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", 1, SchedStrategy.FCFS).on(P1);
        Entry E10 = new Entry(model, "E10").on(T1);
        Entry E11 = new Entry(model, "E11").on(T1);
        Entry E12 = new Entry(model, "E12").on(T1);
        Entry E13 = new Entry(model, "E13").on(T1);

        // Layer 3: database
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2);
        Entry E20 = new Entry(model, "E20").on(T2);
        Entry E21 = new Entry(model, "E21").on(T2);
        Entry E22 = new Entry(model, "E22").on(T2);
        Entry E23 = new Entry(model, "E23").on(T2);

        // Layer 4: cache (10 items, capacity 2, LRU, uniform access)
        int totalitems = 10;
        int cachecapacity = 2;
        Matrix accessProbs = new Matrix(1, totalitems);
        for (int i = 0; i < totalitems; i++) {
            accessProbs.set(0, i, 1.0 / totalitems);
        }
        DiscreteSampler pAccess = new DiscreteSampler(accessProbs);
        Processor P3 = new Processor(model, "P3", 1, SchedStrategy.PS);
        CacheTask T3 = new CacheTask(model, "T3", totalitems, cachecapacity, ReplacementStrategy.LRU, 1);
        T3.on(P3);
        ItemEntry E3 = new ItemEntry(model, "E3", totalitems, pAccess).on(T3);

        // Client activities
        Activity A0 = new Activity(model, "A0", new Exp(1.0)).on(T0).boundTo(E0).synchCall(E12, 1.0);
        Activity A1 = new Activity(model, "A1", new Exp(1.0)).on(T0).synchCall(E10, 1.0);
        Activity A2 = new Activity(model, "A2", new Exp(1.0)).on(T0).synchCall(E11, 1.0);
        Activity A3 = new Activity(model, "A3", new Exp(1.0)).on(T0).synchCall(E13, 1.0);

        // Application activities
        Activity B0 = new Activity(model, "B0", new Exp(1.0)).on(T1).boundTo(E10);
        Activity B1 = new Activity(model, "B1", new Exp(1.0)).on(T1).repliesTo(E10);
        Activity B2 = new Activity(model, "B2", new Exp(1.0)).on(T1).boundTo(E11);
        Activity B3 = new Activity(model, "B3", new Exp(1.0)).on(T1).synchCall(E21, 1.0).repliesTo(E11);
        Activity B4 = new Activity(model, "B4", new Exp(1.0)).on(T1).boundTo(E12).synchCall(E20, 1.0).repliesTo(E12);
        Activity B5 = new Activity(model, "B5", new Exp(1.0)).on(T1).boundTo(E13);
        Activity B6 = new Activity(model, "B6", new Exp(1.0)).on(T1);
        Activity B7 = new Activity(model, "B7", new Exp(1.0)).on(T1).synchCall(E22, 1.0);
        Activity B7a = new Activity(model, "B7a", new Exp(1.0)).on(T1).repliesTo(E13);
        Activity B7b = new Activity(model, "B7b", new Exp(1.0)).on(T1).synchCall(E23, 1.0).repliesTo(E13);

        // Database activities
        Activity C0 = new Activity(model, "C0", new Exp(1.0)).on(T2).boundTo(E20);
        Activity C1 = new Activity(model, "C1", new Exp(1.0)).on(T2).synchCall(E3, 1.0).repliesTo(E20);
        Activity C2 = new Activity(model, "C2", new Exp(1.0)).on(T2).boundTo(E21).synchCall(E3, 1.0).repliesTo(E21);
        Activity C3 = new Activity(model, "C3", new Exp(1.0)).on(T2).boundTo(E22);
        Activity C4 = new Activity(model, "C4", new Exp(1.0)).on(T2);
        Activity C5 = new Activity(model, "C5", new Exp(1.0)).on(T2).repliesTo(E22);
        Activity C6 = new Activity(model, "C6", new Exp(1.0)).on(T2).boundTo(E23).repliesTo(E23);

        // Cache activities
        Activity D0 = new Activity(model, "D0", new Immediate()).on(T3).boundTo(E3);
        Activity D1a = new Activity(model, "D1a", new Exp(1.0)).on(T3).repliesTo(E3);
        Activity D1b = new Activity(model, "D1b", new Exp(0.5)).on(T3).repliesTo(E3);

        // Precedences
        T0.addPrecedence(ActivityPrecedence.Serial(A0, A1, A2, A3));
        T1.addPrecedence(ActivityPrecedence.Serial(B0, B1));
        T1.addPrecedence(ActivityPrecedence.Serial(B2, B3));
        T1.addPrecedence(ActivityPrecedence.Serial(B5, B6, B7));
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.7);
        probs.set(0, 1, 0.3);
        List<Activity> b7targets = new ArrayList<Activity>();
        b7targets.add(B7a);
        b7targets.add(B7b);
        T1.addPrecedence(ActivityPrecedence.OrFork(B7, b7targets, probs));
        T2.addPrecedence(ActivityPrecedence.Serial(C0, C1));
        T2.addPrecedence(ActivityPrecedence.Serial(C3, C4, C5));
        List<Activity> cacheTargets = new ArrayList<Activity>();
        cacheTargets.add(D1a);
        cacheTargets.add(D1b);
        T3.addPrecedence(ActivityPrecedence.CacheAccess(D0, cacheTargets));

        return model;
    }

    public static Network gallery_fj_open() {
        Network model = new Network("Fork-Join-Open");
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "class1");
        source.setArrival(oclass, new Exp(0.05));
        queue1.setService(oclass, new Exp(1.0));
        queue2.setService(oclass, new Exp(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, source, fork, 1.0);
        P.set(oclass, oclass, fork, queue1, 1.0);
        P.set(oclass, oclass, fork, queue2, 1.0);
        P.set(oclass, oclass, queue1, join, 1.0);
        P.set(oclass, oclass, queue2, join, 1.0);
        P.set(oclass, oclass, join, sink, 1.0);
        model.link(P);
        return model;
    }

    public static Network gallery_fj_closed() {
        Network model = new Network("Fork-Join-Closed");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        ClosedClass oclass = new ClosedClass(model, "class1", 5, delay);
        delay.setService(oclass, new Exp(1.0));
        queue1.setService(oclass, new Exp(1.0));
        queue2.setService(oclass, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, delay, fork, 1.0);
        P.set(oclass, oclass, fork, queue1, 1.0);
        P.set(oclass, oclass, fork, queue2, 1.0);
        P.set(oclass, oclass, queue1, join, 1.0);
        P.set(oclass, oclass, queue2, join, 1.0);
        P.set(oclass, oclass, join, delay, 1.0);
        model.link(P);
        return model;
    }

    /**
     * A closed fork-join whose join fires on a 2-of-3 QUORUM: the third sibling is discarded when
     * it arrives. SolverLDES and SolverJMT reproduce it exactly; SolverMVA and SolverNC charge the
     * second order statistic of the branch completion times (FJ_ordstat_exp).
     */
    public static Network gallery_fj_quorum() {
        Network model = new Network("Fork-Join-Quorum");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.PS);
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        ClosedClass oclass = new ClosedClass(model, "class1", 5, delay);
        delay.setService(oclass, new Exp(1.0));
        queue1.setService(oclass, new Exp(2.0));
        queue2.setService(oclass, new Exp(2.0));
        queue3.setService(oclass, new Exp(2.0));
        join.setStrategy(oclass, JoinStrategy.PARTIAL);
        join.setRequired(oclass, 2);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, delay, fork, 1.0);
        P.set(oclass, oclass, fork, queue1, 1.0);
        P.set(oclass, oclass, fork, queue2, 1.0);
        P.set(oclass, oclass, fork, queue3, 1.0);
        P.set(oclass, oclass, queue1, join, 1.0);
        P.set(oclass, oclass, queue2, join, 1.0);
        P.set(oclass, oclass, queue3, join, 1.0);
        P.set(oclass, oclass, join, delay, 1.0);
        model.link(P);
        return model;
    }

    public static Network gallery_cache_lru() {
        Network model = new Network("Cache-LRU");
        int n = 5; // number of items
        int m = 2; // cache capacity
        Delay delay = new Delay(model, "Delay");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.LRU);
        ClosedClass jobClass = new ClosedClass(model, "JobClass", 1, delay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, delay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, delay, 0);
        delay.setService(jobClass, new Exp(1));
        Matrix pmatrix = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            pmatrix.set(0, i, 1.0 / n);
        }
        DiscreteSampler pAccess = new DiscreteSampler(pmatrix);
        cacheNode.setRead(jobClass, pAccess);
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, delay, cacheNode, 1.0);
        P.set(hitClass, jobClass, cacheNode, delay, 1.0);
        P.set(missClass, jobClass, cacheNode, delay, 1.0);
        model.link(P);
        return model;
    }

    public static Network gallery_cache_routing() {
        Network model = new Network("Cache-Routing");
        int n = 4; // number of items
        int m = 2; // cache capacity
        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, m, ReplacementStrategy.LRU);
        Queue hitQueue = new Queue(model, "HitQueue", SchedStrategy.FCFS);
        Queue missQueue = new Queue(model, "MissQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "InitClass", 0);
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);
        source.setArrival(jobClass, new Exp(1));
        hitQueue.setService(hitClass, new Exp(2.0));
        missQueue.setService(missClass, new Exp(1.0));
        Matrix pmatrix = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            pmatrix.set(0, i, 1.0 / n);
        }
        DiscreteSampler pAccess = new DiscreteSampler(pmatrix);
        cacheNode.setRead(jobClass, pAccess);
        cacheNode.setHitClass(jobClass, hitClass);
        cacheNode.setMissClass(jobClass, missClass);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, cacheNode, 1.0);
        P.set(hitClass, hitClass, cacheNode, hitQueue, 1.0);
        P.set(hitClass, hitClass, hitQueue, sink, 1.0);
        P.set(missClass, missClass, cacheNode, missQueue, 1.0);
        P.set(missClass, missClass, missQueue, sink, 1.0);
        model.link(P);
        return model;
    }

    public static LayeredNetwork gallery_lqn_basic() {
        LayeredNetwork model = new LayeredNetwork("LQN-Basic");
        Processor P1 = new Processor(model, "P1", 2, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 3, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF).on(P1).setThinkTime(new Exp(1.0 / 2));
        Task T2 = new Task(model, "T2", 50, SchedStrategy.FCFS).on(P1).setThinkTime(new Exp(1.0 / 3));
        Task T3 = new Task(model, "T3", 25, SchedStrategy.FCFS).on(P2).setThinkTime(new Exp(1.0 / 4));
        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T3);
        Activity A1 = new Activity(model, "AS1", new Exp(10)).on(T1).boundTo(E1).synchCall(E2, 1);
        Activity A2 = new Activity(model, "AS2", new Exp(20)).on(T2).boundTo(E2).synchCall(E3, 5).repliesTo(E2);
        Activity A3 = new Activity(model, "AS3", new Exp(50)).on(T3).boundTo(E3).repliesTo(E3);
        return model;
    }

    public static LayeredNetwork gallery_lqn_workflows() {
        LayeredNetwork model = new LayeredNetwork("LQN-Workflows");
        Processor P1 = new Processor(model, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        Task T1 = new Task(model, "T1", 1, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Immediate());
        Entry E1 = new Entry(model, "Entry").on(T1);

        Processor P2 = new Processor(model, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        Task T2 = new Task(model, "T2", Integer.MAX_VALUE, SchedStrategy.INF).on(P2).setThinkTime(new Immediate());
        Entry E2 = new Entry(model, "E2").on(T2);

        Processor P3 = new Processor(model, "P3", 5, SchedStrategy.PS);
        Task T3 = new Task(model, "T3", Integer.MAX_VALUE, SchedStrategy.INF).on(P3);
        T3.setThinkTime(Exp.fitMean(10));
        Entry E3 = new Entry(model, "E3").on(T3);

        Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1);
        Activity A2 = new Activity(model, "A2", Exp.fitMean(2)).on(T1);
        Activity A3 = new Activity(model, "A3", Exp.fitMean(3)).on(T1).synchCall(E2, 1.0);

        Activity B1 = new Activity(model, "B1", Exp.fitMean(0.1)).on(T2).boundTo(E2);
        Activity B2 = new Activity(model, "B2", Exp.fitMean(0.2)).on(T2);
        Activity B3 = new Activity(model, "B3", Exp.fitMean(0.3)).on(T2);
        Activity B4 = new Activity(model, "B4", Exp.fitMean(0.4)).on(T2);
        Activity B5 = new Activity(model, "B5", Exp.fitMean(0.5)).on(T2);
        Activity B6 = new Activity(model, "B6", Exp.fitMean(0.6)).on(T2).synchCall(E3, 1.0).repliesTo(E2);

        Activity C1 = new Activity(model, "C1", Exp.fitMean(0.1)).on(T3).boundTo(E3);
        Activity C2 = new Activity(model, "C2", Exp.fitMean(0.2)).on(T3);
        Activity C3 = new Activity(model, "C3", Exp.fitMean(0.3)).on(T3);
        Activity C4 = new Activity(model, "C4", Exp.fitMean(0.4)).on(T3);
        Activity C5 = new Activity(model, "C5", Exp.fitMean(0.5)).on(T3).repliesTo(E3);

        List<Activity> loopActs = new ArrayList<Activity>();
        loopActs.add(A2);
        loopActs.add(A3);
        T1.addPrecedence(ActivityPrecedence.Loop(A1, loopActs, 3));
        T2.addPrecedence(ActivityPrecedence.Serial(B4, B5));
        List<Activity> andForkTargets = new ArrayList<Activity>();
        andForkTargets.add(B2);
        andForkTargets.add(B3);
        andForkTargets.add(B4);
        T2.addPrecedence(ActivityPrecedence.AndFork(B1, andForkTargets));
        List<Activity> andJoinSources = new ArrayList<Activity>();
        andJoinSources.add(B2);
        andJoinSources.add(B3);
        andJoinSources.add(B5);
        T2.addPrecedence(ActivityPrecedence.AndJoin(andJoinSources, B6));
        List<Activity> orForkTargets = new ArrayList<Activity>();
        orForkTargets.add(C2);
        orForkTargets.add(C3);
        orForkTargets.add(C4);
        Matrix orProbs = new Matrix(1, 3);
        orProbs.set(0, 0, 0.3);
        orProbs.set(0, 1, 0.3);
        orProbs.set(0, 2, 0.4);
        T3.addPrecedence(ActivityPrecedence.OrFork(C1, orForkTargets, orProbs));
        List<Activity> orJoinSources = new ArrayList<Activity>();
        orJoinSources.add(C2);
        orJoinSources.add(C3);
        orJoinSources.add(C4);
        T3.addPrecedence(ActivityPrecedence.OrJoin(orJoinSources, C5));
        return model;
    }

    public static Network gallery_mm1k() {
        return gallery_mm1k(3);
    }

    public static Network gallery_mm1k(int K) {
        Network model = new Network("M/M/1/K");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(1);
        queue.setCapacity(K);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1", 0);
        source.setArrival(oclass, new Exp(0.8));
        queue.setService(oclass, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, source, queue, 1.0);
        P.set(oclass, oclass, queue, sink, 1.0);
        model.link(P);
        return model;
    }

    public static Network gallery_fcr() {
        return gallery_fcr(3);
    }

    public static Network gallery_fcr(int K) {
        Network model = new Network("FCR-Dropping");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1", 0);
        source.setArrival(oclass, new Exp(0.8));
        queue.setService(oclass, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oclass, oclass, source, queue, 1.0);
        P.set(oclass, oclass, queue, sink, 1.0);
        model.link(P);
        List<Node> regionNodes = new ArrayList<Node>();
        regionNodes.add(queue);
        Region fcr = model.addRegion(regionNodes);
        fcr.setGlobalMaxJobs(K);
        fcr.setDropRule(oclass, true);
        return model;
    }

    public static Environment gallery_renv_breakdown() {
        Network model = new Network("ServerWithFailures");
        Source source = new Source(model, "Arrivals");
        Queue queue = new Queue(model, "Server", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Departures");
        OpenClass jobclass = new OpenClass(model, "Jobs");
        source.setArrival(jobclass, new Exp(0.8));
        queue.setService(jobclass, new Exp(2.0));
        queue.setNumberOfServers(1);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, source, queue, 1.0);
        P.set(jobclass, jobclass, queue, sink, 1.0);
        model.link(P);
        Environment env = new Environment("ServerEnv");
        env.addNodeFailureRepair(model, queue, new Exp(0.1), new Exp(1.0), new Exp(0.5));
        env.init();
        return env;
    }

    public static Network gallery_qn_random() {
        return gallery_qn_random(23000L);
    }

    public static Network gallery_qn_random(long seed) {
        // Closed queueing network: 3 queues, 1 delay, 2 closed classes (always stable).
        // FCFS/Exp with a deterministic (cyclic) topology so output is reproducible
        // under the given seed.
        NetworkGenerator generator = new NetworkGenerator(
            "fcfs", "Probabilities", "Exp", "medium",
            false, false, false, false, NetworkGenerator::cyclicGraph);
        generator.setSeed(seed);
        return generator.generate(3, 1, 0, 2);
    }

    public static LayeredNetwork gallery_lqn_random() {
        return gallery_lqn_random(23000L);
    }

    public static LayeredNetwork gallery_lqn_random(long seed) {
        // Random layered network: 1 client, 2 levels, 4 tasks, 2 processors.
        LayeredNetworkGenerator generator = new LayeredNetworkGenerator();
        generator.setSeed(seed);
        return generator.generate(1, 2, 4, 2);
    }
}
