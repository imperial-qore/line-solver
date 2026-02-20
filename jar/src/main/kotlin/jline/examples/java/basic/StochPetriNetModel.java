/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.*;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.nodes.Place;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Transition;
import jline.lang.processes.Coxian;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.lang.processes.Pareto;
import jline.solvers.jmt.JMT;
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
