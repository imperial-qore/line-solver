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
public class PetriStateTestFixtures {
    
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
     * Helper to create global sync events for enabling/firing
     */
    public static GlobalSync createGlobalSync(int node, EventType eventType, int mode) {
        GlobalSync gsync = new GlobalSync();
        ModeEvent activeEvent = new ModeEvent(eventType, node, mode);
        
        List<ModeEvent> active = new ArrayList<>();
        active.add(activeEvent);
        gsync.setActive(active);
        
        List<ModeEvent> passive = new ArrayList<>();
        if (eventType == EventType.FIRE) {
            // For FIRE events, we need PRE and POST events
            // PRE event consumes from input place (P1 at index 0)
            ModeEvent preEvent = new ModeEvent(EventType.PRE, 0, mode);
            passive.add(preEvent);
            // POST event produces to output place (P1 at index 0 for self-loop)
            ModeEvent postEvent = new ModeEvent(EventType.POST, 0, mode);
            passive.add(postEvent);
        } else {
            // For other events (ENABLE), use PHASE event
            ModeEvent passiveEvent = new ModeEvent(EventType.PHASE, 0, mode);
            passive.add(passiveEvent);
        }
        gsync.setPassive(passive);
        
        return gsync;
    }
}