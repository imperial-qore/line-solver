/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.*;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Place;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Transition;
import jline.lang.processes.Coxian;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.lang.processes.Immediate;
import jline.lang.processes.Pareto;
import jline.solvers.wrappers.jmt.JMT;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Examples of stochastic Petri net models
 */
public class StochPetriNetModel {

    /**
     * Basic open stochastic Petri net with single transition.
     * <p>
     * Features:
     * - Open network: Source → Place → Transition → Sink
     * - Single transition T1 with exponential firing time (rate 4.0)
     * - Transition requires 1 token from P1 to fire
     * - Infinite server capacity for transition
     * - Demonstrates basic Petri net structure in LINE
     *
     * @return configured basic stochastic Petri net model
     */
    public static Network spn_basic_open() {
        Network model = new Network("model");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");
        Place P = new Place(model, "P1");
        Transition T = new Transition(model, "T1");

        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, Exp.fitMean(1.0));

        Mode mode1 = T.addMode("Mode1");
        T.setNumberOfServers(mode1, Integer.MAX_VALUE);
        T.setDistribution(mode1, new Exp(4));
        T.setEnablingConditions(mode1, jobclass, P, 1);
        T.setFiringOutcome(mode1, jobclass, sink, 1);

        model.link(Network.serialRouting(source, P, T, sink));

        NetworkStruct sn = model.getStruct();
        TransitionNodeParam np = (TransitionNodeParam) sn.nodeparam.get(T);

        return model;
    }

    /**
     * Complex open stochastic Petri net with immediate transitions.
     * <p>
     * Features:
     * - 7 places and 8 transitions with mixed timing strategies
     * - Immediate transitions (T2, T3, T4, T5) with priorities and weights
     * - Timed transitions with Exp and Erlang distributions
     * - Inhibiting conditions (T5 inhibited by P6)
     * - Multiple enabling conditions and firing outcomes per transition
     * - Initial state configuration with tokens in P1 and P5
     *
     * @return configured complex stochastic Petri net model
     */
    public static Network spn_open_sevenplaces() {
        Network model = new Network("model");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");
        List<Place> P = new ArrayList<>();
        for (int i = 0; i < 7; i++) {
            P.add(new Place(model, "P" + (i + 1)));
        }
        List<Transition> T = new ArrayList<>();
        for (int i = 0; i < 8; i++) {
            T.add(new Transition(model, "T" + i + 1));
        }
        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, Exp.fitMean(1.0));

        // T1
        Mode mode1 = T.get(0).addMode("Mode1");
        T.get(0).setNumberOfServers(mode1, Integer.MAX_VALUE);
        T.get(0).setDistribution(mode1, new Exp(4));
        T.get(0).setEnablingConditions(mode1, jobclass, P.get(0), 1);
        T.get(0).setFiringOutcome(mode1, jobclass, P.get(1), 1);

        // T2
        Mode mode2 = T.get(1).addMode("Mode1");
        T.get(1).setNumberOfServers(mode2, Integer.MAX_VALUE);
        T.get(1).setEnablingConditions(mode2, jobclass, P.get(1), 1);
        T.get(1).setFiringOutcome(mode2, jobclass, P.get(2), 1);
        T.get(1).setTimingStrategy(mode2, TimingStrategy.IMMEDIATE);
        T.get(1).setFiringPriorities(mode2, 1);
        T.get(1).setFiringWeights(mode2, 1);

        // T3
        Mode mode3 = T.get(2).addMode("Mode1");
        T.get(2).setNumberOfServers(mode3, Integer.MAX_VALUE);
        T.get(2).setEnablingConditions(mode3, jobclass, P.get(1), 1);
        T.get(2).setFiringOutcome(mode3, jobclass, P.get(3), 1);
        T.get(2).setTimingStrategy(mode3, TimingStrategy.IMMEDIATE);
        T.get(2).setFiringPriorities(mode3, 1);
        // T.get(2).setFiringWeights(1,0.6);
        T.get(2).setFiringPriorities(mode3, 1);

        // T4
        Mode mode4 = T.get(3).addMode("Mode1");
        T.get(3).setNumberOfServers(mode4, Integer.MAX_VALUE);
        T.get(3).setEnablingConditions(mode4, jobclass, P.get(2), 1);
        T.get(3).setEnablingConditions(mode4, jobclass, P.get(4), 1);
        T.get(3).setFiringOutcome(mode4, jobclass, P.get(4), 1);
        T.get(3).setFiringOutcome(mode4, jobclass, P.get(5), 1);
        T.get(3).setTimingStrategy(mode4, TimingStrategy.IMMEDIATE);
        T.get(3).setFiringPriorities(mode4, 1);

        // T5
        Mode mode5 = T.get(4).addMode("Mode1");
        T.get(4).setNumberOfServers(mode5, Integer.MAX_VALUE);
        T.get(4).setEnablingConditions(mode5, jobclass, P.get(3), 1);
        T.get(4).setEnablingConditions(mode5, jobclass, P.get(4), 1);
        T.get(4).setFiringOutcome(mode5, jobclass, P.get(6), 1);
        T.get(4).setInhibitingConditions(mode5, jobclass, P.get(5), 1);
        T.get(4).setTimingStrategy(mode5, TimingStrategy.IMMEDIATE);
        T.get(4).setFiringPriorities(mode5, 1);

        // T6
        Mode mode6 = T.get(5).addMode("Mode1");
        T.get(5).setNumberOfServers(mode6, Integer.MAX_VALUE);
        T.get(5).setDistribution(mode6, new Erlang(2, 2));
        T.get(5).setEnablingConditions(mode6, jobclass, P.get(5), 1);
        T.get(5).setFiringOutcome(mode6, jobclass, P.get(0), 1);

        // T7
        Mode mode7 = T.get(6).addMode("Mode1");
        T.get(6).setNumberOfServers(mode7, Integer.MAX_VALUE);
        T.get(6).setDistribution(mode7, new Exp(2));
        T.get(6).setEnablingConditions(mode7, jobclass, P.get(6), 1);
        T.get(6).setFiringOutcome(mode7, jobclass, P.get(0), 1);
        T.get(6).setFiringOutcome(mode7, jobclass, P.get(4), 1);

        // T8
        Mode mode8 = T.get(7).addMode("Mode1");
        T.get(7).setNumberOfServers(mode8, Integer.MAX_VALUE);
        T.get(7).setDistribution(mode8, new Exp(2));
        T.get(7).setEnablingConditions(mode8, jobclass, P.get(3), 1);
        T.get(7).setFiringOutcome(mode8, jobclass, sink, 1);

        RoutingMatrix routingMatrix = model.initRoutingMatrix(); // initialize routing matrix
        routingMatrix.set(jobclass, jobclass, source, P.get(0), 1.0); // (Source,Class1) -> (P1,Class1)

        routingMatrix.set(jobclass, jobclass, P.get(0), T.get(0), 1.0); // (P1,Class1) -> (T1,Class1)
        routingMatrix.set(jobclass, jobclass, P.get(1), T.get(1), 1.0); // (P2,Class1) -> (T2,Class1)
        routingMatrix.set(jobclass, jobclass, P.get(1), T.get(2), 1.0); // (P2,Class1) -> (T3,Class1)
        routingMatrix.set(jobclass, jobclass, P.get(2), T.get(3), 1.0); // (P3,Class1) -> (T4,Class1)
        routingMatrix.set(jobclass, jobclass, P.get(3), T.get(4), 1.0); // (P4,Class1) -> (T5,Class1)
        routingMatrix.set(jobclass, jobclass, P.get(4), T.get(3), 1.0); // (P5,Class1) -> (T4,Class1)
        routingMatrix.set(jobclass, jobclass, P.get(4), T.get(4), 1.0); // (P5,Class1) -> (T5,Class1)
        routingMatrix.set(jobclass, jobclass, P.get(5), T.get(4), 1.0); // (P6,Class1) -> (T5,Class1)
        routingMatrix.set(jobclass, jobclass, P.get(5), T.get(5), 1.0); // (P6,Class1) -> (T6,Class1)
        routingMatrix.set(jobclass, jobclass, P.get(6), T.get(6), 1.0); // (P7,Class1) -> (T7,Class1)

        routingMatrix.set(jobclass, jobclass, T.get(0), P.get(1), 1.0); // (T1,Class1) -> (P2,Class1)
        routingMatrix.set(jobclass, jobclass, T.get(1), P.get(2), 1.0); // (T2,Class1) -> (P3,Class1)
        routingMatrix.set(jobclass, jobclass, T.get(2), P.get(3), 1.0); // (T3,Class1) -> (P4,Class1)
        routingMatrix.set(jobclass, jobclass, T.get(3), P.get(4), 1.0); // (T4,Class1) -> (P5,Class1)
        routingMatrix.set(jobclass, jobclass, T.get(3), P.get(5), 1.0); // (T4,Class1) -> (P6,Class1)
        routingMatrix.set(jobclass, jobclass, T.get(4), P.get(6), 1.0); // (T5,Class1) -> (P7,Class1)
        routingMatrix.set(jobclass, jobclass, T.get(5), P.get(0), 1.0); // (T6,Class1) -> (P1,Class1)
        routingMatrix.set(jobclass, jobclass, T.get(6), sink, 1.0); // (T7,Class1) -> (Sink,Class1)
        routingMatrix.set(jobclass, jobclass, T.get(6), P.get(0), 1.0); // (T7,Class1) -> (P1,Class1)
        routingMatrix.set(jobclass, jobclass, T.get(6), P.get(4), 1.0); // (T7,Class1) -> (P5,Class1)

        routingMatrix.set(jobclass, jobclass, P.get(3), T.get(7), 1.0); // (P4,Class1) -> (T5,Class1)
        routingMatrix.set(jobclass, jobclass, T.get(7), sink, 1.0); // (T8,Class1) -> (Sink,Class1)

        model.link(routingMatrix);
        // Set Initial State
        source.setState(Matrix.singleton(0));
        P.get(0).setState(Matrix.singleton(2));
        P.get(1).setState(Matrix.singleton(0));
        P.get(2).setState(Matrix.singleton(0));
        P.get(3).setState(Matrix.singleton(0));
        P.get(4).setState(Matrix.singleton(1));
        P.get(5).setState(Matrix.singleton(0));
        P.get(6).setState(Matrix.singleton(0));
        return model;
    }

    /**
     * Closed stochastic Petri net with batch processing.
     * <p>
     * Features:
     * - Closed system with 10 tokens circulating between 2 places
     * - T1 requires 4 tokens to fire, produces 4 tokens
     * - T2 requires 2 tokens to fire, produces 2 tokens
     * - Demonstrates batch token processing in Petri nets
     * - All tokens initially placed in P1
     *
     * @return configured batch processing Petri net model
     */
    public static Network spn_twomodes() {
        Network model = new Network("model");

        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");

        ClosedClass jobclass = new ClosedClass(model, "Class1", 10, P1, 0);

        // T1
        Mode mode1 = T1.addMode("Mode1");
        T1.setDistribution(mode1, new Exp(2));
        T1.setEnablingConditions(mode1, jobclass, P1, 4);
        T1.setFiringOutcome(mode1, jobclass, P2, 4);

        // T2
        Mode mode2 = T2.addMode("Mode2");
        T2.setDistribution(mode2, new Exp(3));
        T2.setEnablingConditions(mode2, jobclass, P2, 2);
        T2.setFiringOutcome(mode2, jobclass, P1, 2);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass, jobclass, P1, T1, 1.0);
        routingMatrix.set(jobclass, jobclass, P2, T2, 1.0);
        routingMatrix.set(jobclass, jobclass, T1, P2, 1.0);
        routingMatrix.set(jobclass, jobclass, T2, P1, 1.0);

        model.link(routingMatrix);

        P1.setState(Matrix.singleton(jobclass.getPopulation()));

        return model;
    }

    /**
     * Closed stochastic Petri net with competing transitions.
     * <p>
     * Features:
     * - 8 tokens in closed system with 3 places
     * - T1 and T2 compete for tokens from P1 (require 2 and 3 tokens respectively)
     * - T3 and T4 return tokens to P1 from P2 and P3
     * - Different firing rates create resource competition
     * - Demonstrates resource contention in Petri nets
     *
     * @return configured competing transitions Petri net model
     */
    public static Network spn_fourmodes() {
        Network model = new Network("model");

        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Place P3 = new Place(model, "P3");

        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        Transition T3 = new Transition(model, "T3");
        Transition T4 = new Transition(model, "T4");

        ClosedClass jobclass = new ClosedClass(model, "Class1", 8, P1, 0);

        Mode mode1 = T1.addMode("Mode1");
        T1.setDistribution(mode1, new Exp(2));
        T1.setEnablingConditions(mode1, jobclass, P1, 2);
        T1.setFiringOutcome(mode1, jobclass, P2, 2);

        Mode mode2 = T2.addMode("Mode2");
        T2.setDistribution(mode2, new Exp(1));
        T2.setEnablingConditions(mode2, jobclass, P1, 3);
        T2.setFiringOutcome(mode2, jobclass, P3, 3);

        Mode mode3 = T3.addMode("Mode3");
        T3.setDistribution(mode3, new Exp(4));
        T3.setEnablingConditions(mode3, jobclass, P2, 1);
        T3.setFiringOutcome(mode3, jobclass, P1, 1);

        Mode mode4 = T4.addMode("Mode4");
        T4.setDistribution(mode4, new Exp(2));
        T4.setEnablingConditions(mode4, jobclass, P3, 2);
        T4.setFiringOutcome(mode4, jobclass, P1, 2);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass, jobclass, P1, T1, 1.0);
        routingMatrix.set(jobclass, jobclass, P1, T2, 1.0);
        routingMatrix.set(jobclass, jobclass, P2, T3, 1.0);
        routingMatrix.set(jobclass, jobclass, P3, T4, 1.0);

        routingMatrix.set(jobclass, jobclass, T1, P2, 1.0);
        routingMatrix.set(jobclass, jobclass, T2, P3, 1.0);
        routingMatrix.set(jobclass, jobclass, T3, P1, 1.0);
        routingMatrix.set(jobclass, jobclass, T4, P1, 1.0);

        model.link(routingMatrix);

        P1.setState(Matrix.singleton(jobclass.getPopulation()));

        return model;
    }

    /**
     * Closed stochastic Petri net with multiple firing modes and inhibition.
     * <p>
     * Features:
     * - 4 tokens in closed system with 3 places
     * - T1 has two firing modes: Mode1 (2 tokens → P2), Mode2 (1 token → P3)
     * - T3 has inhibiting condition: fires only when P2 has no tokens
     * - Demonstrates mode-based firing and inhibiting arcs
     * - Complex token flow patterns with conditional transitions
     *
     * @return configured multi-mode Petri net with inhibition
     */
    public static Network spn_inhibiting() {
        // Closed model with multiple firing modes and inhibiting conditions
        Network model = new Network("model");

        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Place P3 = new Place(model, "P3");

        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        Transition T3 = new Transition(model, "T3");

        ClosedClass jobclass = new ClosedClass(model, "Class1", 4, P1, 0);

        Mode mode1 = T1.addMode("Mode1");
        T1.setDistribution(mode1, new Exp(2));
        T1.setEnablingConditions(mode1, jobclass, P1, 2);
        T1.setFiringOutcome(mode1, jobclass, P2, 2);

        Mode mode2 = T1.addMode("Mode2");
        T1.setDistribution(mode2, new Exp(1));
        T1.setEnablingConditions(mode2, jobclass, P1, 1);
        T1.setFiringOutcome(mode2, jobclass, P3, 1);

        Mode mode3 = T2.addMode("Mode3");
        T2.setDistribution(mode3, new Exp(4));
        T2.setEnablingConditions(mode3, jobclass, P2, 1);
        T2.setFiringOutcome(mode3, jobclass, P1, 1);

        Mode mode4 = T3.addMode("Mode4");
        T3.setDistribution(mode4, new Exp(1));
        T3.setEnablingConditions(mode4, jobclass, P3, 3);
        T3.setInhibitingConditions(mode4, jobclass, P2, 1);
        T3.setFiringOutcome(mode4, jobclass, P1, 3);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass, jobclass, P1, T1, 1.0);
        routingMatrix.set(jobclass, jobclass, P2, T2, 1.0);
        routingMatrix.set(jobclass, jobclass, P2, T3, 1.0);
        routingMatrix.set(jobclass, jobclass, P3, T3, 1.0);

        routingMatrix.set(jobclass, jobclass, T1, P2, 1.0);
        routingMatrix.set(jobclass, jobclass, T1, P3, 1.0);
        routingMatrix.set(jobclass, jobclass, T2, P1, 1.0);
        routingMatrix.set(jobclass, jobclass, T3, P1, 1.0);

        model.link(routingMatrix);

        P1.setState(Matrix.singleton(jobclass.getPopulation()));

        return model;
    }

    /**
     * Closed stochastic Petri net with diverse service distributions.
     * <p>
     * Features:
     * - 2 tokens circulating through 4 places in series
     * - Different firing distributions: Exp, Erlang, HyperExp, Coxian
     * - T4 uses custom Coxian distribution with specified phases
     * - Demonstrates various probability distributions in Petri nets
     * - All transitions require and produce 2 tokens (synchronous firing)
     *
     * @return configured Petri net with diverse distributions
     */
    public static Network spn_closed_fourplaces() {
        // Closed model with Erlang, HyperExp, Coxian distributions
        Network model = new Network("model");

        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Place P3 = new Place(model, "P3");
        Place P4 = new Place(model, "P4");

        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        Transition T3 = new Transition(model, "T3");
        Transition T4 = new Transition(model, "T4");

        ClosedClass jobclass = new ClosedClass(model, "Class1", 2, P1, 0);

        Mode mode1 = T1.addMode("Mode1");
        T1.setDistribution(mode1, new Exp(2));
        T1.setEnablingConditions(mode1, jobclass, P1, 2);
        T1.setFiringOutcome(mode1, jobclass, P2, 2);

        Mode mode2 = T2.addMode("Mode2");
        T2.setDistribution(mode2, new Erlang(3, 4));
        T2.setEnablingConditions(mode2, jobclass, P2, 2);
        T2.setFiringOutcome(mode2, jobclass, P3, 2);

        Mode mode3 = T3.addMode("Mode3");
        T3.setDistribution(mode3, new HyperExp(0.7, 3, 1.5));
        T3.setEnablingConditions(mode3, jobclass, P3, 2);
        T3.setFiringOutcome(mode3, jobclass, P4, 2);

        Matrix mu0 = new Matrix(2, 1);
        mu0.set(0, 0, 1.0);
        mu0.set(1, 0, 2.0);

        Matrix phi0 = new Matrix(2, 1);
        phi0.set(0, 0, 0.6);
        phi0.set(1, 0, 1.0);

        Mode mode4 = T4.addMode("Mode4");
        T4.setDistribution(mode4, new Coxian(mu0, phi0));
        T4.setEnablingConditions(mode4, jobclass, P4, 2);
        T4.setFiringOutcome(mode4, jobclass, P1, 2);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass, jobclass, P1, T1, 1.0);
        routingMatrix.set(jobclass, jobclass, P2, T2, 1.0);
        routingMatrix.set(jobclass, jobclass, P3, T3, 1.0);
        routingMatrix.set(jobclass, jobclass, P4, T4, 1.0);

        routingMatrix.set(jobclass, jobclass, T1, P2, 1.0);
        routingMatrix.set(jobclass, jobclass, T2, P3, 1.0);
        routingMatrix.set(jobclass, jobclass, T3, P4, 1.0);
        routingMatrix.set(jobclass, jobclass, T4, P1, 1.0);

        model.link(routingMatrix);

        P1.setState(Matrix.singleton(jobclass.getPopulation()));

        return model;
    }

    /**
     * Multi-class closed stochastic Petri net.
     * <p>
     * Features:
     * - Two job classes: Class1 (10 tokens), Class2 (7 tokens)
     * - T1 has different modes for each class with different requirements
     * - Class1: 2 tokens required, Class2: 1 token required
     * - T2 and T3 handle different classes with different batch sizes
     * - Demonstrates multi-class token management in Petri nets
     *
     * @return configured multi-class stochastic Petri net model
     */
    public static Network spn_closed_twoplaces() {
        // Closed model with multiple job classes
        Network model = new Network("model");

        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        Transition T3 = new Transition(model, "T3");

        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, P1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 7, P1, 0);

        // T1
        Mode mode1 = T1.addMode("Mode1");
        T1.setDistribution(mode1, new Exp(2));
        T1.setEnablingConditions(mode1, jobclass1, P1, 2);
        T1.setFiringOutcome(mode1, jobclass1, P2, 2);

        // T1
        Mode mode2 = T1.addMode("Mode2");
        T1.setDistribution(mode2, new Exp(3));
        T1.setEnablingConditions(mode2, jobclass2, P1, 1);
        T1.setFiringOutcome(mode2, jobclass2, P2, 1);

        // T2
        Mode mode3 = T2.addMode("Mode3");
        T2.setDistribution(mode3, new Erlang(1.5, 2));
        T2.setEnablingConditions(mode3, jobclass1, P2, 1);
        T2.setFiringOutcome(mode3, jobclass1, P1, 1);

        // T3
        Mode mode4 = T3.addMode("Mode4");
        T3.setDistribution(mode4, new Exp(0.5));
        T3.setEnablingConditions(mode4, jobclass2, P2, 4);
        T3.setFiringOutcome(mode4, jobclass2, P1, 4);

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, P1, T1, 1.0);
        routingMatrix.set(jobclass2, jobclass2, P1, T1, 1.0);
        routingMatrix.set(jobclass1, jobclass1, P2, T2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, P2, T3, 1.0);
        routingMatrix.set(jobclass1, jobclass1, T1, P2, 1.0);
        routingMatrix.set(jobclass2, jobclass2, T1, P2, 1.0);
        routingMatrix.set(jobclass1, jobclass1, T2, P1, 1.0);
        routingMatrix.set(jobclass2, jobclass2, T3, P1, 1.0);

        model.link(routingMatrix);

        P1.setState(new Matrix("[10,7]"));

        return model;
    }

    /**
     * Single-place Petri net with multiple firing modes.
     * <p>
     * Features:
     * - Single place P1 with 1 token
     * - Single transition T1 with 3 different firing modes:
     * - Mode1: Exponential distribution (mean 1.0)
     * - Mode2: Erlang distribution (mean 1.0, order 2)
     * - Mode3: HyperExponential distribution (mean 1.0, SCV 4.0)
     * - Self-loop: transition fires and returns token to same place
     * - Demonstrates multiple stochastic modes in single transition
     *
     * @return configured multi-mode single-place Petri net
     */
    public static Network spn_basic_closed() {
        Network model = new Network("model");

        // Places
        Place P1 = new Place(model, "P1");

        // Transition
        Transition T1 = new Transition(model, "T1");

        // Job Class
        ClosedClass jobclass = new ClosedClass(model, "Class1", 1, P1, 0);

        // Mode 1: Exponential with mean 1
        Mode mode1 = T1.addMode("Mode1");
        T1.setDistribution(mode1, Exp.fitMean(1.0)); // mean = 1
        T1.setEnablingConditions(mode1, jobclass, P1, 1);
        T1.setFiringOutcome(mode1, jobclass, P1, 1);

        // Mode 2: Erlang with mean 1 and order 2
        Mode mode2 = T1.addMode("Mode2");
        T1.setDistribution(mode2, Erlang.fitMeanAndOrder(1, 2)); // mean = 1 -> rate = 2, k = 2
        T1.setEnablingConditions(mode2, jobclass, P1, 1);
        T1.setFiringOutcome(mode2, jobclass, P1, 1);

        // Mode 3: HyperExponential with mean 1 and SCV = 4
        Mode mode3 = T1.addMode("Mode3");
        T1.setDistribution(mode3, HyperExp.fitMeanAndSCV(1.0, 4.0)); // mean = 1, SCV = 4
        T1.setEnablingConditions(mode3, jobclass, P1, 1);
        T1.setFiringOutcome(mode3, jobclass, P1, 1);

        // Routing Matrix
        RoutingMatrix routingMatrix = model.initRoutingMatrix();

        routingMatrix.set(jobclass, jobclass, P1, T1, 1.0);
        routingMatrix.set(jobclass, jobclass, T1, P1, 1.0);

        model.link(routingMatrix);

        // Set initial state
        P1.setState(Matrix.singleton(jobclass.getPopulation()));

        return model;
    }


/**
     * Open stochastic Petri net with Pareto service time.
     * <p>
     * Features:
     * - Open network: Source → Place → Transition → Sink
     * - Single transition T1 with Pareto firing time (shape=3, scale=1)
     * - Pareto distribution has mean = shape*scale/(shape-1) = 3*1/(3-1) = 1.5
     * - Demonstrates non-Markovian (heavy-tailed) firing times in Petri nets
     * - Single server capacity for transition
     *
     * @return configured stochastic Petri net model with Pareto service
     */
    public static Network spn_pareto_service() {
        Network model = new Network("model");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");
        Place P = new Place(model, "P1");
        Transition T = new Transition(model, "T1");

        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, new Exp(0.5)); // arrival rate 0.5

        // T1 with Pareto service time
        // Pareto(shape, scale) - shape must be >= 2
        // With shape=3 and scale=1, mean service time = 3*1/(3-1) = 1.5
        Mode mode1 = T.addMode("Mode1");
        T.setNumberOfServers(mode1, 1);
        T.setDistribution(mode1, new Pareto(3, 1)); // Pareto with shape=3, scale=1
        T.setEnablingConditions(mode1, jobclass, P, 1);
        T.setFiringOutcome(mode1, jobclass, sink, 1);

        model.link(Network.serialRouting(source, P, T, sink));

        return model;
    }

    /**
     * A product-form net solved ANALYTICALLY, and one that no queueing network
     * expresses.
     *
     * <p>The cycle P0 -&gt; T0 -&gt; P1 -&gt; T1 -&gt; P2 -&gt; T2 -&gt; P0 has a
     * product form, so {@code new NC(model)} solves it exactly through
     * {@link jline.api.spn.Spn_pf} (complex balance),
     * {@link jline.api.mdd.Mdd_rec} (the normalising constant by one walk of the
     * decision diagram holding the reachable set) and
     * {@link jline.api.spn.Spn_metrics}.</p>
     *
     * @return the 3-place cyclic net at N = 4 with rates {1, 1.5, 2}
     */
    public static Network spn_productform_cyclic() {
        Network model = new Network("spn");
        Place[] pl = new Place[3];
        Transition[] tr = new Transition[3];
        double[] rates = {1.0, 1.5, 2.0};
        for (int i = 0; i < 3; i++) {
            pl[i] = new Place(model, "P" + i);
            tr[i] = new Transition(model, "T" + i);
        }
        ClosedClass jc = new ClosedClass(model, "Class1", 4, pl[0], 0);
        for (int i = 0; i < 3; i++) {
            Mode m = tr[i].addMode("fire");
            tr[i].setDistribution(m, new Exp(rates[i]));
            tr[i].setNumberOfServers(m, Integer.valueOf(1));
            tr[i].setEnablingConditions(m, jc, pl[i], 1);
            tr[i].setFiringOutcome(m, jc, pl[(i + 1) % 3], 1);
        }
        RoutingMatrix P = model.initRoutingMatrix();
        for (int i = 0; i < 3; i++) {
            P.set(jc, jc, pl[i], tr[i], 1.0);
            P.set(jc, jc, tr[i], pl[(i + 1) % 3], 1.0);
        }
        model.link(P);
        for (int i = 0; i < 3; i++) {
            pl[i].setState(Matrix.singleton(i == 0 ? 4 : 0));
        }
        return model;
    }

    /**
     * P0 -(Tf)-&gt; P1 + P2 -(Tj)-&gt; P3 -(Tb)-&gt; P0.
     *
     * <p>Tf consumes ONE token and produces TWO, Tj the reverse, so the marking
     * is not a conserved job population and there is no queueing-network
     * counterpart -- the limitation the MDD-rec paper opens with. Its place
     * invariant is 2*m0 + m1 + m2 + 2*m3, not the token count. SolverNC's
     * {@code rec} method solves it exactly.</p>
     *
     * @return the fork-join net with 3 tokens at P0
     */
    public static Network spn_productform_forkjoin() {
        Network model = new Network("fj");
        Place[] pl = new Place[4];
        for (int i = 0; i < 4; i++) {
            pl[i] = new Place(model, "P" + i);
        }
        Transition tf = new Transition(model, "Tf");
        Transition tj = new Transition(model, "Tj");
        Transition tb = new Transition(model, "Tb");
        ClosedClass jc = new ClosedClass(model, "C", 3, pl[0], 0);

        Mode mf = tf.addMode("f");
        tf.setDistribution(mf, new Exp(1.3));
        tf.setNumberOfServers(mf, Integer.valueOf(1));
        tf.setEnablingConditions(mf, jc, pl[0], 1);
        tf.setFiringOutcome(mf, jc, pl[1], 1);
        tf.setFiringOutcome(mf, jc, pl[2], 1);

        Mode mj = tj.addMode("j");
        tj.setDistribution(mj, new Exp(0.7));
        tj.setNumberOfServers(mj, Integer.valueOf(1));
        tj.setEnablingConditions(mj, jc, pl[1], 1);
        tj.setEnablingConditions(mj, jc, pl[2], 1);
        tj.setFiringOutcome(mj, jc, pl[3], 1);

        Mode mb = tb.addMode("b");
        tb.setDistribution(mb, new Exp(1.9));
        tb.setNumberOfServers(mb, Integer.valueOf(1));
        tb.setEnablingConditions(mb, jc, pl[3], 1);
        tb.setFiringOutcome(mb, jc, pl[0], 1);

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jc, jc, pl[0], tf, 1.0);
        P.set(jc, jc, tf, pl[1], 1.0);
        P.set(jc, jc, tf, pl[2], 1.0);
        P.set(jc, jc, pl[1], tj, 1.0);
        P.set(jc, jc, pl[2], tj, 1.0);
        P.set(jc, jc, tj, pl[3], 1.0);
        P.set(jc, jc, pl[3], tb, 1.0);
        P.set(jc, jc, tb, pl[0], 1.0);
        model.link(P);
        for (int i = 0; i < 4; i++) {
            pl[i].setState(Matrix.singleton(i == 0 ? 3 : 0));
        }
        return model;
    }


    /**
     * Liu (1998) Fig. 2b: four servers in a line, blocking before service.
     *
     * <p>{@code t1 -> p5 -> t2 -> p4 -> t3 -> p3 -> t4}, with p2, p1 and p0
     * holding the free slots of the three finite buffers (3, 2 and 4). Each
     * buffer is a conserved pair of places, so the net is a strongly connected
     * marked graph and all four transitions carry the same throughput.</p>
     *
     * @param mu the four firing rates
     * @return the production line of the paper's Table 2
     */
    public static Network spn_lpbounds_prodline(double[] mu) {
        Network model = new Network("liu98");
        Place p5 = new Place(model, "p5");
        Place p4 = new Place(model, "p4");
        Place p3 = new Place(model, "p3");
        Place p2 = new Place(model, "p2");
        Place p1 = new Place(model, "p1");
        Place p0 = new Place(model, "p0");
        Transition t1 = new Transition(model, "t1");
        Transition t2 = new Transition(model, "t2");
        Transition t3 = new Transition(model, "t3");
        Transition t4 = new Transition(model, "t4");
        ClosedClass jc = new ClosedClass(model, "Class1", 9, p2, 0);

        Mode m1 = t1.addMode("m1");
        t1.setDistribution(m1, new Exp(mu[0]));
        t1.setEnablingConditions(m1, jc, p2, 1);
        t1.setFiringOutcome(m1, jc, p5, 1);

        Mode m2 = t2.addMode("m2");
        t2.setDistribution(m2, new Exp(mu[1]));
        t2.setEnablingConditions(m2, jc, p5, 1);
        t2.setEnablingConditions(m2, jc, p1, 1);
        t2.setFiringOutcome(m2, jc, p4, 1);
        t2.setFiringOutcome(m2, jc, p2, 1);

        Mode m3 = t3.addMode("m3");
        t3.setDistribution(m3, new Exp(mu[2]));
        t3.setEnablingConditions(m3, jc, p4, 1);
        t3.setEnablingConditions(m3, jc, p0, 1);
        t3.setFiringOutcome(m3, jc, p3, 1);
        t3.setFiringOutcome(m3, jc, p1, 1);

        Mode m4 = t4.addMode("m4");
        t4.setDistribution(m4, new Exp(mu[3]));
        t4.setEnablingConditions(m4, jc, p3, 1);
        t4.setFiringOutcome(m4, jc, p0, 1);

        RoutingMatrix R = model.initRoutingMatrix();
        R.set(jc, jc, p2, t1, 1.0);
        R.set(jc, jc, t1, p5, 1.0);
        R.set(jc, jc, p5, t2, 1.0);
        R.set(jc, jc, p1, t2, 1.0);
        R.set(jc, jc, t2, p4, 1.0);
        R.set(jc, jc, t2, p2, 1.0);
        R.set(jc, jc, p4, t3, 1.0);
        R.set(jc, jc, p0, t3, 1.0);
        R.set(jc, jc, t3, p3, 1.0);
        R.set(jc, jc, t3, p1, 1.0);
        R.set(jc, jc, p3, t4, 1.0);
        R.set(jc, jc, t4, p0, 1.0);
        model.link(R);

        p5.setState(Matrix.singleton(0));
        p4.setState(Matrix.singleton(0));
        p3.setState(Matrix.singleton(0));
        p2.setState(Matrix.singleton(3));
        p1.setState(Matrix.singleton(2));
        p0.setState(Matrix.singleton(4));
        return model;
    }

    /**
     * Closed queueing Petri net with two queueing places.
     *
     * <p>A CPU (single-server FCFS queueing place) and a think stage (infinite
     * server) exchange a fixed population through two immediate transitions.
     * This is the queueing-Petri-net rendering of a machine-repairman model, so
     * {@link #spn_queueing_place_ref} gives the exact cross-check.</p>
     *
     * @param N the token population
     * @return the queueing Petri net
     */
    public static Network spn_queueing_place(int N) {
        Network model = new Network("QueueingPetriNet");
        // setService turns each Place into a queueing place whose embedded queue
        // is served under the scheduling strategy of its constructor.
        Place cpu = new Place(model, "CPU", SchedStrategy.FCFS);
        Place think = new Place(model, "Think", SchedStrategy.INF);
        ClosedClass jobs = new ClosedClass(model, "Jobs", N, think, 0);
        cpu.setService(jobs, new Exp(1.5));
        think.setService(jobs, new Exp(0.5));

        Transition toCPU = new Transition(model, "toCPU");
        Mode m1 = toCPU.addMode("m1");
        toCPU.setTimingStrategy(m1, TimingStrategy.IMMEDIATE);
        toCPU.setEnablingConditions(m1, jobs, think, 1);
        toCPU.setFiringOutcome(m1, jobs, cpu, 1);

        Transition toThink = new Transition(model, "toThink");
        Mode m2 = toThink.addMode("m2");
        toThink.setTimingStrategy(m2, TimingStrategy.IMMEDIATE);
        toThink.setEnablingConditions(m2, jobs, cpu, 1);
        toThink.setFiringOutcome(m2, jobs, think, 1);

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobs, jobs, think, toCPU, 1.0);
        P.set(jobs, jobs, toCPU, cpu, 1.0);
        P.set(jobs, jobs, cpu, toThink, 1.0);
        P.set(jobs, jobs, toThink, think, 1.0);
        model.link(P);

        think.setState(Matrix.singleton(N));
        cpu.setState(Matrix.singleton(0));
        return model;
    }

    /**
     * The finite-population Delay + M/M/1 network the queueing Petri net of
     * {@link #spn_queueing_place} is equivalent to.
     *
     * @param N the closed population
     * @return the reference queueing network
     */
    public static Network spn_queueing_place_ref(int N) {
        Network ref = new Network("ref");
        Delay delay = new Delay(ref, "Think");
        Queue queue = new Queue(ref, "CPU", SchedStrategy.FCFS);
        ClosedClass cc = new ClosedClass(ref, "Jobs", N, delay, 0);
        delay.setService(cc, new Exp(0.5));
        queue.setService(cc, new Exp(1.5));
        ref.link(Network.serialRouting(delay, queue));
        return ref;
    }

    /**
     * A closed net whose fluid answer is EXACT.
     *
     * <p>Every mode is infinite-server with one input arc, so
     * {@code min(m/w, Inf) = m} and the drift is LINEAR: the fluid mean is then
     * the exact mean and the covariance the exact covariance (a binomial
     * marking, {@code Binomial(4, 3/5)}, variance 0.96).</p>
     *
     * @return the two-place linear-drift net
     */
    public static Network spn_fluid_exact() {
        Network exact = new Network("spn_fluid_exact");
        Place P1 = new Place(exact, "P1");
        Place P2 = new Place(exact, "P2");
        Transition T1 = new Transition(exact, "T1");
        Transition T2 = new Transition(exact, "T2");
        ClosedClass jc = new ClosedClass(exact, "Class1", 4, P1, 0);
        Mode m1 = T1.addMode("Mode1");
        T1.setNumberOfServers(m1, Integer.valueOf(Integer.MAX_VALUE));
        T1.setDistribution(m1, new Exp(2));
        T1.setEnablingConditions(m1, jc, P1, 1);
        T1.setFiringOutcome(m1, jc, P2, 1);
        Mode m2 = T2.addMode("Mode2");
        T2.setNumberOfServers(m2, Integer.valueOf(Integer.MAX_VALUE));
        T2.setDistribution(m2, new Exp(3));
        T2.setEnablingConditions(m2, jc, P2, 1);
        T2.setFiringOutcome(m2, jc, P1, 1);
        RoutingMatrix R = exact.initRoutingMatrix();
        R.set(jc, jc, P1, T1, 1.0);
        R.set(jc, jc, T1, P2, 1.0);
        R.set(jc, jc, P2, T2, 1.0);
        R.set(jc, jc, T2, P1, 1.0);
        exact.link(R);
        P1.setState(Matrix.singleton(jc.getPopulation()));
        P2.setState(Matrix.singleton(0));
        return exact;
    }

    /**
     * {@code P1 -T1-> P2 -(immediate)-> P3 -T3-> P1}.
     *
     * <p>The vanishing place P2 holds exactly zero mass, so the net answers as
     * {@link #spn_fluid_reduced} does. That is what the algebraic flow buys: an
     * approximation of the immediate transition by a large finite rate would
     * only approach it.</p>
     *
     * @return the three-place net with one immediate transition
     */
    public static Network spn_fluid_immediate() {
        Network imm = new Network("spn_fluid_immediate");
        Place Q1 = new Place(imm, "P1");
        Place Q2 = new Place(imm, "P2");
        Place Q3 = new Place(imm, "P3");
        Transition U1 = new Transition(imm, "T1");
        Transition Ui = new Transition(imm, "Ti");
        Transition U3 = new Transition(imm, "T3");
        ClosedClass jq = new ClosedClass(imm, "Class1", 4, Q1, 0);
        Mode a1 = U1.addMode("M1");
        U1.setDistribution(a1, new Exp(2));
        U1.setEnablingConditions(a1, jq, Q1, 1);
        U1.setFiringOutcome(a1, jq, Q2, 1);
        Mode ai = Ui.addMode("Mi");
        Ui.setDistribution(ai, new Immediate());
        Ui.setTimingStrategy(ai, TimingStrategy.IMMEDIATE);
        Ui.setEnablingConditions(ai, jq, Q2, 1);
        Ui.setFiringOutcome(ai, jq, Q3, 1);
        Mode a3 = U3.addMode("M3");
        U3.setDistribution(a3, new Exp(3));
        U3.setEnablingConditions(a3, jq, Q3, 1);
        U3.setFiringOutcome(a3, jq, Q1, 1);
        RoutingMatrix Ri = imm.initRoutingMatrix();
        Ri.set(jq, jq, Q1, U1, 1.0);
        Ri.set(jq, jq, U1, Q2, 1.0);
        Ri.set(jq, jq, Q2, Ui, 1.0);
        Ri.set(jq, jq, Ui, Q3, 1.0);
        Ri.set(jq, jq, Q3, U3, 1.0);
        Ri.set(jq, jq, U3, Q1, 1.0);
        imm.link(Ri);
        Q1.setState(Matrix.singleton(jq.getPopulation()));
        Q2.setState(Matrix.singleton(0));
        Q3.setState(Matrix.singleton(0));
        return imm;
    }

    /**
     * {@link #spn_fluid_immediate} with the immediate transition eliminated by
     * hand, which is the answer the algebraic flow must reproduce.
     *
     * @return the reduced two-place net
     */
    public static Network spn_fluid_reduced() {
        Network red = new Network("spn_fluid_reduced");
        Place W1 = new Place(red, "P1");
        Place W3 = new Place(red, "P3");
        Transition V1 = new Transition(red, "T1");
        Transition V3 = new Transition(red, "T3");
        ClosedClass jr = new ClosedClass(red, "Class1", 4, W1, 0);
        Mode b1 = V1.addMode("M1");
        V1.setDistribution(b1, new Exp(2));
        V1.setEnablingConditions(b1, jr, W1, 1);
        V1.setFiringOutcome(b1, jr, W3, 1);
        Mode b3 = V3.addMode("M3");
        V3.setDistribution(b3, new Exp(3));
        V3.setEnablingConditions(b3, jr, W3, 1);
        V3.setFiringOutcome(b3, jr, W1, 1);
        RoutingMatrix Rr = red.initRoutingMatrix();
        Rr.set(jr, jr, W1, V1, 1.0);
        Rr.set(jr, jr, V1, W3, 1.0);
        Rr.set(jr, jr, W3, V3, 1.0);
        Rr.set(jr, jr, V3, W1, 1.0);
        red.link(Rr);
        W1.setState(Matrix.singleton(jr.getPopulation()));
        W3.setState(Matrix.singleton(0));
        return red;
    }

    /**
     * An M/M/1 queue written as an OPEN Petri net: a Source Exp(lambda) feeds
     * place P1, whose tokens drain through a single-server Transition Exp(mu)
     * to a Sink.
     *
     * @param lambda the arrival rate
     * @param mu     the firing rate
     * @return the open net
     */
    public static Network spn_nrm_mm1(double lambda, double mu) {
        Network model = new Network("mm1spn");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");
        Place P1 = new Place(model, "P1");
        Transition T1 = new Transition(model, "T1");
        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, new Exp(lambda));
        Mode mode = T1.addMode("Mode1");
        T1.setDistribution(mode, new Exp(mu));
        T1.setEnablingConditions(mode, jobclass, P1, 1);
        T1.setFiringOutcome(mode, jobclass, sink, 1);
        RoutingMatrix R = model.initRoutingMatrix();
        R.set(jobclass, jobclass, source, P1, 1.0);
        R.set(jobclass, jobclass, P1, T1, 1.0);
        R.set(jobclass, jobclass, T1, sink, 1.0);
        model.link(R);
        return model;
    }

    /**
     * The open tandem of {@link #spn_nrm_mm1}: two places in series, each an
     * M/M/1 queue at its own rate.
     *
     * @param lambda the arrival rate
     * @param mu1    the firing rate of T1
     * @param mu2    the firing rate of T2
     * @return the open tandem net
     */
    public static Network spn_nrm_tandem(double lambda, double mu1, double mu2) {
        Network model = new Network("tandemspn");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");
        Place P1 = new Place(model, "P1");
        Place P2 = new Place(model, "P2");
        Transition T1 = new Transition(model, "T1");
        Transition T2 = new Transition(model, "T2");
        OpenClass jobclass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobclass, new Exp(lambda));
        Mode m1 = T1.addMode("Mode1");
        T1.setDistribution(m1, new Exp(mu1));
        T1.setEnablingConditions(m1, jobclass, P1, 1);
        T1.setFiringOutcome(m1, jobclass, P2, 1);
        Mode m2 = T2.addMode("Mode1");
        T2.setDistribution(m2, new Exp(mu2));
        T2.setEnablingConditions(m2, jobclass, P2, 1);
        T2.setFiringOutcome(m2, jobclass, sink, 1);
        RoutingMatrix R = model.initRoutingMatrix();
        R.set(jobclass, jobclass, source, P1, 1.0);
        R.set(jobclass, jobclass, P1, T1, 1.0);
        R.set(jobclass, jobclass, T1, P2, 1.0);
        R.set(jobclass, jobclass, P2, T2, 1.0);
        R.set(jobclass, jobclass, T2, sink, 1.0);
        model.link(R);
        return model;
    }

    /**
     * Colored generalized stochastic Petri net (CGSPN).
     * <p>
     * A token color is a job class, so a place holds one marking per color and
     * a transition declares one mode per color it serves. Two colors, Gold and
     * Silver, circulate between a Buffer place and a Server place:
     * <ul>
     * <li>the immediate transition admit has one mode per color, each held off
     * by inhibiting arcs on BOTH colors, so the server is a mutual-exclusion
     * resource and the firing weights arbitrate between the colors;</li>
     * <li>the timed transition serve returns the token at a color-dependent
     * rate.</li>
     * </ul>
     * Admission is immediate, so the server is never idle: the weights split
     * the completions 2:1 in favour of Gold, hence X_Gold = 2*X_Silver,
     * U_c = X_c/mu_c and U_Gold + U_Silver = 1.
     *
     * @return configured colored GSPN
     */
    public static Network spn_colored_gspn() {
        Network model = new Network("ColoredGSPN");

        // Declare every node before parameterizing any of it: the enabling,
        // inhibiting and firing matrices are indexed by (node, class) over the
        // nodes that exist when they are first set.
        Place buf = new Place(model, "Buffer");
        Place srv = new Place(model, "Server");
        Transition admit = new Transition(model, "admit");
        Transition serve = new Transition(model, "serve");

        // Two token colors, both starting in the buffer
        ClosedClass gold = new ClosedClass(model, "Gold", 2, buf, 0);
        ClosedClass silver = new ClosedClass(model, "Silver", 2, buf, 0);

        // Immediate transition: one mode per color, admitting a token only when
        // the server holds no token of either color. The weights set the mix.
        Mode mg = admit.addMode("gold");
        admit.setTimingStrategy(mg, TimingStrategy.IMMEDIATE);
        admit.setEnablingConditions(mg, gold, buf, 1);
        admit.setInhibitingConditions(mg, gold, srv, 1);
        admit.setInhibitingConditions(mg, silver, srv, 1);
        admit.setFiringOutcome(mg, gold, srv, 1);
        admit.setFiringWeights(mg, 2.0);

        Mode ms = admit.addMode("silver");
        admit.setTimingStrategy(ms, TimingStrategy.IMMEDIATE);
        admit.setEnablingConditions(ms, silver, buf, 1);
        admit.setInhibitingConditions(ms, gold, srv, 1);
        admit.setInhibitingConditions(ms, silver, srv, 1);
        admit.setFiringOutcome(ms, silver, srv, 1);
        admit.setFiringWeights(ms, 1.0);

        // Timed transition: one mode per color, color-dependent service rates
        Mode sg = serve.addMode("gold");
        serve.setDistribution(sg, new Exp(3.0));
        serve.setEnablingConditions(sg, gold, srv, 1);
        serve.setFiringOutcome(sg, gold, buf, 1);

        Mode ss = serve.addMode("silver");
        serve.setDistribution(ss, new Exp(1.5));
        serve.setEnablingConditions(ss, silver, srv, 1);
        serve.setFiringOutcome(ss, silver, buf, 1);

        // Topology: both colors follow Buffer -> admit -> Server -> serve -> Buffer
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        for (JobClass c : new JobClass[]{gold, silver}) {
            routingMatrix.set(c, c, buf, admit, 1.0);
            routingMatrix.set(c, c, admit, srv, 1.0);
            routingMatrix.set(c, c, srv, serve, 1.0);
            routingMatrix.set(c, c, serve, buf, 1.0);
        }
        model.link(routingMatrix);

        // Initial marking: two tokens of each color in the buffer, server empty
        Matrix bufMarking = new Matrix(1, 2);
        bufMarking.set(0, 0, 2);
        bufMarking.set(0, 1, 2);
        buf.setState(bufMarking);
        srv.setState(new Matrix(1, 2));

        return model;
    }

    /**
     * Main method for testing and demonstrating stochastic Petri net examples.
     *
     * <p>Currently configured to:
     * - Run spn_basic_closed() with multiple firing modes
     * - Solve using JMT solver with specified seed (23000)
     * - Print average performance metrics
     * - Launch JMT simulation GUI viewer
     * - SSA analysis is commented out
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        Network model = spn_basic_closed();
        
        System.out.println("SSA Analysis");
//        new SSA(model).getAvgTable().print();
//        System.out.println("JMT Analysis");
        new JMT(model, "seed", 23000).getAvgTable().print();
//        new JMT(model, "seed", 23000).jsimgView();
    }
}
