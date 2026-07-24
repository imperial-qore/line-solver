package jline.lang.state;

import jline.lang.*;
import jline.lang.constant.EventType;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Test models for State event testing, replicating example_stochPetriNet_8
 */
public class TestSPNModels {
    
    /**
     * Creates a test SPN model with 3 modes at T1:
     * - Mode 1: Exponential distribution (mean=1, rate=1)
     * - Mode 2: Erlang distribution (mean=1, order=2, gives 2 phases with rate=2)
     * - Mode 3: HyperExp distribution (mean=1, SCV=4)
     * 
     * This model structure matches example_stochPetriNet_8.m
     */
    public static Network createTestModelWithThreeModes() {
        Network model = new Network("model");
        
        // Create single place (self-loop model)
        Place P1 = new Place(model, "P1");
        
        // Create transition
        Transition T1 = new Transition(model, "T1");
        
        // Create job class with population=1
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, P1, 0);
        
        // Mode 1: Exponential with mean=1 (rate=1)
        Mode mode1 = T1.addMode("Mode1");
        T1.setNumberOfServers(mode1, 1);
        T1.setDistribution(mode1, Exp.fitMean(1.0));
        T1.setEnablingConditions(mode1, jobClass, P1, 1);
        T1.setFiringOutcome(mode1, jobClass, P1, 1); // Self-loop back to P1
        
        // Mode 2: Erlang with mean=1 and order=2 (2 phases, each with rate=2)
        Mode mode2 = T1.addMode("Mode2");
        T1.setNumberOfServers(mode2, 1);
        T1.setDistribution(mode2, Erlang.fitMeanAndOrder(1.0, 2));
        T1.setEnablingConditions(mode2, jobClass, P1, 1);
        T1.setFiringOutcome(mode2, jobClass, P1, 1); // Self-loop back to P1
        
        // Mode 3: HyperExp with mean=1 and SCV=4
        Mode mode3 = T1.addMode("Mode3");
        T1.setNumberOfServers(mode3, 1);
        T1.setDistribution(mode3, HyperExp.fitMeanAndSCV(1.0, 4.0));
        T1.setEnablingConditions(mode3, jobClass, P1, 1);
        T1.setFiringOutcome(mode3, jobClass, P1, 1); // Self-loop back to P1
        
        // Set up routing (self-loop: P1 -> T1 -> P1)
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, P1, T1, 1.0);
        routingMatrix.set(jobClass, jobClass, T1, P1, 1.0);
        model.link(routingMatrix);
        
        // Set initial state
        P1.setState(Matrix.singleton(jobClass.getPopulation()));
        
        return model;
    }
    
    /**
     * Creates a test SPN model with 3 identical Exp(1) modes.
     * Expected throughput: 3.0 (race of 3 Exp(1) gives E[min]=1/3)
     */
    public static Network createTestModelWithThreeExpModes() {
        Network model = new Network("model");

        // Create single place (self-loop model)
        Place P1 = new Place(model, "P1");

        // Create transition
        Transition T1 = new Transition(model, "T1");

        // Create job class with population=1
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, P1, 0);

        // Mode 1: Exponential with mean=1
        Mode mode1 = T1.addMode("Mode1");
        T1.setNumberOfServers(mode1, 1);
        T1.setDistribution(mode1, Exp.fitMean(1.0));
        T1.setEnablingConditions(mode1, jobClass, P1, 1);
        T1.setFiringOutcome(mode1, jobClass, P1, 1);

        // Mode 2: Exponential with mean=1
        Mode mode2 = T1.addMode("Mode2");
        T1.setNumberOfServers(mode2, 1);
        T1.setDistribution(mode2, Exp.fitMean(1.0));
        T1.setEnablingConditions(mode2, jobClass, P1, 1);
        T1.setFiringOutcome(mode2, jobClass, P1, 1);

        // Mode 3: Exponential with mean=1
        Mode mode3 = T1.addMode("Mode3");
        T1.setNumberOfServers(mode3, 1);
        T1.setDistribution(mode3, Exp.fitMean(1.0));
        T1.setEnablingConditions(mode3, jobClass, P1, 1);
        T1.setFiringOutcome(mode3, jobClass, P1, 1);

        // Set up routing (self-loop: P1 -> T1 -> P1)
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, P1, T1, 1.0);
        routingMatrix.set(jobClass, jobClass, T1, P1, 1.0);
        model.link(routingMatrix);

        // Set initial state
        P1.setState(Matrix.singleton(jobClass.getPopulation()));

        return model;
    }

    /**
     * Creates a test SPN model with 3 identical Erlang(1,2) modes.
     * Expected throughput: 3 * 1/E[min(Erlang_1, Erlang_2, Erlang_3)]
     */
    public static Network createTestModelWithThreeErlangModes() {
        Network model = new Network("model");

        Place P1 = new Place(model, "P1");
        Transition T1 = new Transition(model, "T1");
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, P1, 0);

        // 3 Erlang(mean=1, order=2) modes
        for (int i = 1; i <= 3; i++) {
            Mode mode = T1.addMode("Mode" + i);
            T1.setNumberOfServers(mode, 1);
            T1.setDistribution(mode, Erlang.fitMeanAndOrder(1.0, 2));
            T1.setEnablingConditions(mode, jobClass, P1, 1);
            T1.setFiringOutcome(mode, jobClass, P1, 1);
        }

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, P1, T1, 1.0);
        routingMatrix.set(jobClass, jobClass, T1, P1, 1.0);
        model.link(routingMatrix);

        P1.setState(Matrix.singleton(jobClass.getPopulation()));
        return model;
    }

    /**
     * Creates a test SPN model with 3 identical HyperExp(1,4) modes.
     */
    public static Network createTestModelWithThreeHyperExpModes() {
        Network model = new Network("model");

        Place P1 = new Place(model, "P1");
        Transition T1 = new Transition(model, "T1");
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, P1, 0);

        // 3 HyperExp(mean=1, SCV=4) modes
        for (int i = 1; i <= 3; i++) {
            Mode mode = T1.addMode("Mode" + i);
            T1.setNumberOfServers(mode, 1);
            T1.setDistribution(mode, HyperExp.fitMeanAndSCV(1.0, 4.0));
            T1.setEnablingConditions(mode, jobClass, P1, 1);
            T1.setFiringOutcome(mode, jobClass, P1, 1);
        }

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, P1, T1, 1.0);
        routingMatrix.set(jobClass, jobClass, T1, P1, 1.0);
        model.link(routingMatrix);

        P1.setState(Matrix.singleton(jobClass.getPopulation()));
        return model;
    }

    /**
     * Creates a test SPN model with multiple servers for testing global events
     * Similar structure but with 2 servers for mode1 and mode3
     */
    public static Network createTestModelWithMultipleServers() {
        Network model = new Network("model");

        // Create single place (self-loop model)
        Place P1 = new Place(model, "P1");
        
        // Create transition
        Transition T1 = new Transition(model, "T1");
        
        // Create job class with more population to test multiple servers
        ClosedClass jobClass = new ClosedClass(model, "Class1", 10, P1, 0);
        
        // Mode 1: Exponential with 2 servers
        Mode mode1 = T1.addMode("Mode1");
        T1.setNumberOfServers(mode1, 2);
        T1.setDistribution(mode1, Exp.fitMean(1.0));
        T1.setEnablingConditions(mode1, jobClass, P1, 1);
        T1.setFiringOutcome(mode1, jobClass, P1, 1); // Self-loop
        
        // Mode 2: Erlang (2 phases)
        Mode mode2 = T1.addMode("Mode2");
        T1.setNumberOfServers(mode2, 1);
        T1.setDistribution(mode2, Erlang.fitMeanAndOrder(1.0, 2));
        T1.setEnablingConditions(mode2, jobClass, P1, 1);
        T1.setFiringOutcome(mode2, jobClass, P1, 1); // Self-loop
        
        // Mode 3: HyperExp with 2 servers
        Mode mode3 = T1.addMode("Mode3");
        T1.setNumberOfServers(mode3, 2);
        T1.setDistribution(mode3, HyperExp.fitMeanAndSCV(1.0, 4.0));
        T1.setEnablingConditions(mode3, jobClass, P1, 1);
        T1.setFiringOutcome(mode3, jobClass, P1, 1); // Self-loop
        
        // Set up routing (self-loop: P1 -> T1 -> P1)
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobClass, jobClass, P1, T1, 1.0);
        routingMatrix.set(jobClass, jobClass, T1, P1, 1.0);
        model.link(routingMatrix);
        
        // Set initial state
        P1.setState(Matrix.singleton(jobClass.getPopulation()));
        
        return model;
    }
    
    /**
     * Helper to create sync events matching the MATLAB test structure
     */
    public static Sync createPhaseSync(int node, EventType eventType, int jobClass) {
        Sync sync = new Sync();
        sync.active.put(0, new Event(eventType, node, jobClass));
        sync.passive.put(0, new Event(EventType.PHASE, node, jobClass));
        return sync;
    }
    
    /**
     * Helper to create global sync events for enabling/firing.
     * For ENABLE events, passive events point to the enabling place (P1 at index 0).
     * For FIRE events, passive events include both PRE (consume) and POST (produce) events.
     */
    public static GlobalSync createGlobalSync(int transitionNode, EventType eventType, int mode) {
        GlobalSync gsync = new GlobalSync();
        List<ModeEvent> activeEvents = new ArrayList<>();
        activeEvents.add(new ModeEvent(eventType, transitionNode, mode, 1.0));
        gsync.setActive(activeEvents);

        List<ModeEvent> passiveEvents = new ArrayList<>();
        // For test models, P1 (place) is always at node index 0
        int placeNode = 0;

        switch (eventType) {
            case ENABLE:
                // Passive events point to enabling places with LOCAL event type
                passiveEvents.add(new ModeEvent(EventType.LOCAL, placeNode, mode, 1.0));
                break;
            case FIRE:
                // PRE event: consume tokens from input place
                passiveEvents.add(new ModeEvent(EventType.PRE, placeNode, mode, 1.0));
                // POST event: produce tokens to output place (self-loop in test model)
                passiveEvents.add(new ModeEvent(EventType.POST, placeNode, mode, 1.0));
                break;
            default:
                // For other events, use PHASE as before
                passiveEvents.add(new ModeEvent(EventType.PHASE, transitionNode, mode));
                break;
        }
        gsync.setPassive(passiveEvents);

        return gsync;
    }
}