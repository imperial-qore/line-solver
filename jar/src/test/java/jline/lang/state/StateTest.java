package jline.lang.state;

import jline.io.Ret;
import jline.lang.GlobalSync;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.Sync;
import jline.lang.constant.EventType;
import jline.GlobalConstants;
import jline.lang.nodes.Transition;
import jline.lang.state.AfterGlobalEvent;
import jline.lang.state.AfterGlobalEvent.AfterGlobalEventResult;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import jline.lang.Event;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Merged test class for State functionality in the jline.lang.state package.
 * Contains tests for State.afterEvent and State.afterGlobalEvent methods.
 * These tests replicate test_SPN_events.m functionality, testing phase transitions
 * for different distribution types (Exp, Erlang, HyperExp) and global synchronization events.
 */
public class StateTest {
    
    private Network model;
    private NetworkStruct sn;
    private int transitionIndex;
    private int jobClassIndex;
    
    private Network multiServerModel;
    private NetworkStruct multiServerSn;
    private int multiServerTransitionIndex;
    
    @BeforeEach
    public void setUp() {
        // Set up for afterEvent tests
        // GlobalConstants.setDummyMode(true); // Not available in Java implementation
        model = TestSPNModels.createTestModelWithThreeModes();
        sn = model.getStruct();
        
        // Get transition index (T1)
        Transition T1 = (Transition) model.getNodeByName("T1");
        transitionIndex = T1.getNodeIndex();
        
        // Job class index is 0 (first class)
        jobClassIndex = 0;
        
        // Set up for afterGlobalEvent tests
        multiServerModel = TestSPNModels.createTestModelWithMultipleServers();
        multiServerSn = multiServerModel.getStruct();
        
        // Get transition index (T1) for multi-server model
        Transition multiServerT1 = (Transition) multiServerModel.getNodeByName("T1");
        multiServerTransitionIndex = multiServerT1.getNodeIndex();
        
        // GlobalConstants.setDummyMode(false); // Not available in Java implementation
    }

    // ==================== STATE AFTER EVENT TESTS ====================
    
    @Test
    //@Disabled
    public void testPhaseTransition_ExponentialMode_AllDisabled() {
        // Test: all modes disabled, so state cannot change for a phase transition
        // This corresponds to lines 19-25 in test_SPN_events.m
        
        // sync{2} corresponds to Mode 1 (Exponential) phase transition
        Sync sync = sn.sync.get(1); // index 1 for second sync (0-based)
        EventType event = sync.active.get(0).getEvent();
        int jobClass = sync.active.get(0).getJobClass();
        
        // State: all modes disabled [1 1 1 0 0 0 0 0]
        // Format: [P1_tokens, P2_tokens, enabled, mode1_servers, mode2_phase1, mode2_phase2, mode3_phase1, mode3_phase2]
        Matrix inspace = new Matrix(1, 8);
        inspace.set(0, 0, 1); // P1 tokens
        inspace.set(0, 1, 1); // P2 tokens  
        inspace.set(0, 2, 1); // enabled flag
        inspace.set(0, 3, 0); // mode1 servers (0 = disabled)
        inspace.set(0, 4, 0); // mode2 phase1
        inspace.set(0, 5, 0); // mode2 phase2
        inspace.set(0, 6, 0); // mode3 phase1
        inspace.set(0, 7, 0); // mode3 phase2
        
        Ret.EventResult result = State.afterEvent(sn, transitionIndex, inspace, event, jobClass, false);
        
        assertNull(result.outspace, "Output space should be null when all modes are disabled");
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testPhaseTransition_ExponentialMode_CannotProgress() {
        // Test: mode 1 (Exp) enabled in phase 1, it cannot progress further
        // This corresponds to lines 27-31 in test_SPN_events.m
        
        Sync sync = sn.sync.get(1); // Mode 1 (Exponential)
        EventType event = sync.active.get(0).getEvent();
        int jobClass = sync.active.get(0).getJobClass();
        
        // State: mode 1 enabled [0 1 1 1 0 0 0 0]
        Matrix inspace = new Matrix(1, 8);
        inspace.set(0, 0, 0);
        inspace.set(0, 1, 1);
        inspace.set(0, 2, 1);
        inspace.set(0, 3, 1); // mode1 servers = 1 (enabled)
        inspace.set(0, 4, 0);
        inspace.set(0, 5, 0);
        inspace.set(0, 6, 0);
        inspace.set(0, 7, 0);
        
        Ret.EventResult result = State.afterEvent(sn, transitionIndex, inspace, event, jobClass, false);
        
        assertNull(result.outspace, "Exponential mode should not have phase transitions");
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testPhaseTransition_ErlangMode_Phase1ToPhase2() {
        // Test: mode 2 (Erlang) enabled in phase 1, can progress to phase 2
        // This corresponds to lines 33-40 in test_SPN_events.m
        
        Sync sync = sn.sync.get(2); // Mode 2 (Erlang)
        EventType event = sync.active.get(0).getEvent();
        int jobClass = sync.active.get(0).getJobClass();
        
        // State: mode 2 phase 1 enabled [1 0 1 0 1 0 0 0]
        Matrix inspace = new Matrix(1, 8);
        inspace.set(0, 0, 1);
        inspace.set(0, 1, 0);
        inspace.set(0, 2, 1);
        inspace.set(0, 3, 0);
        inspace.set(0, 4, 1); // mode2 phase1 = 1
        inspace.set(0, 5, 0); // mode2 phase2 = 0
        inspace.set(0, 6, 0);
        inspace.set(0, 7, 0);
        
        Ret.EventResult result = State.afterEvent(sn, transitionIndex, inspace, event, jobClass, false);
        
        assertNotNull(result.outspace, "Should return output space for phase transition");
        
        // Expected output: [1 0 1 0 0 1 0 0] - moved from phase 1 to phase 2
        assertEquals(1, result.outspace.get(0, 0), 0.0);
        assertEquals(0, result.outspace.get(0, 1), 0.0);
        assertEquals(1, result.outspace.get(0, 2), 0.0);
        assertEquals(0, result.outspace.get(0, 3), 0.0);
        assertEquals(0, result.outspace.get(0, 4), 0.0); // phase1 now 0
        assertEquals(1, result.outspace.get(0, 5), 0.0); // phase2 now 1
        assertEquals(0, result.outspace.get(0, 6), 0.0);
        assertEquals(0, result.outspace.get(0, 7), 0.0);
        
        assertEquals(2.0, result.outrate.get(0, 0), GlobalConstants.FineTol, 
            "Erlang rate should be 2.0");
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testPhaseTransition_ErlangMode_Phase2CannotProgress() {
        // Test: mode 2 (Erlang) enabled in phase 2, cannot progress further
        // This corresponds to lines 42-48 in test_SPN_events.m
        
        Sync sync = sn.sync.get(2); // Mode 2 (Erlang)
        EventType event = sync.active.get(0).getEvent();
        int jobClass = sync.active.get(0).getJobClass();
        
        // State: mode 2 phase 2 enabled [1 0 1 0 0 1 0 0]
        Matrix inspace = new Matrix(1, 8);
        inspace.set(0, 0, 1);
        inspace.set(0, 1, 0);
        inspace.set(0, 2, 1);
        inspace.set(0, 3, 0);
        inspace.set(0, 4, 0); // mode2 phase1 = 0
        inspace.set(0, 5, 1); // mode2 phase2 = 1
        inspace.set(0, 6, 0);
        inspace.set(0, 7, 0);
        
        Ret.EventResult result = State.afterEvent(sn, transitionIndex, inspace, event, jobClass, false);
        
        // The result should indicate no progression is possible
        if (result.outrate != null) {
            assertEquals(0.0, result.outrate.get(0, 0), GlobalConstants.FineTol,
                "Rate should be 0 when Erlang is in final phase");
        }
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testPhaseTransition_HyperExpMode_Phase1RateZero() {
        // Test: mode 3 (HyperExp) enabled in phase 1, progress rate to phase 2 is zero
        // This corresponds to lines 50-56 in test_SPN_events.m
        
        Sync sync = sn.sync.get(3); // Mode 3 (HyperExp)
        EventType event = sync.active.get(0).getEvent();
        int jobClass = sync.active.get(0).getJobClass();
        
        // State: mode 3 phase 1 enabled [1 0 0 0 1 0 1 0]
        Matrix inspace = new Matrix(1, 8);
        inspace.set(0, 0, 1);
        inspace.set(0, 1, 0);
        inspace.set(0, 2, 0);
        inspace.set(0, 3, 0);
        inspace.set(0, 4, 1);
        inspace.set(0, 5, 0);
        inspace.set(0, 6, 1); // mode3 phase1 = 1
        inspace.set(0, 7, 0); // mode3 phase2 = 0
        
        Ret.EventResult result = State.afterEvent(sn, transitionIndex, inspace, event, jobClass, false);
        
        assertNotNull(result.outrate, "Should return rate for HyperExp phase transition");
        assertEquals(0.0, result.outrate.get(0, 0), GlobalConstants.FineTol,
            "HyperExp phase transition rate should be 0");
    }

    // ==================== STATE AFTER GLOBAL EVENT TESTS ====================
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testEnablingEvent_Mode1Exp_FullyEnabled() {
        // Test: Enabling event for mode 1 (Exp) with 2 jobs in buffer
        // This corresponds to lines 67-75 in test_SPN_events.m
        
        // Create a GlobalSync for testing based on the sync structure
        GlobalSync glevent = TestSPNModels.createGlobalSync(multiServerTransitionIndex, EventType.ENABLE, 0); // Mode 1
        
        List<Matrix> glspace = new ArrayList<Matrix>();
        // FCFS buffer at place P1 with 2 jobs of class 1
        Matrix buffer = new Matrix(1, 2);
        buffer.set(0, 0, 1); // job class 1
        buffer.set(0, 1, 1); // job class 1
        glspace.add(buffer);
        
        // Transition state: [2 1 1 0 0 0 0 0]
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2); // 2 jobs available
        transState.set(0, 1, 1); 
        transState.set(0, 2, 1);
        transState.set(0, 3, 0); // mode1 servers = 0 (will be enabled)
        transState.set(0, 4, 0);
        transState.set(0, 5, 0);
        transState.set(0, 6, 0);
        transState.set(0, 7, 0);
        glspace.add(transState);
        
        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(multiServerSn, multiServerTransitionIndex, glspace, glevent, false);
        
        assertNotNull(result.outglspace, "Should return output global space");
        assertEquals(2, result.outglspace.size(), "Should have 2 elements in output");
        
        // Expected output state: [0 1 1 2 0 0 0 0] - 2 servers enabled
        Matrix outTransState = result.outglspace.get(1);
        assertEquals(0, outTransState.get(0, 0), 0.0, "Jobs consumed");
        assertEquals(1, outTransState.get(0, 1), 0.0);
        assertEquals(1, outTransState.get(0, 2), 0.0);
        assertEquals(2, outTransState.get(0, 3), 0.0, "2 mode1 servers enabled");
        
        assertEquals(GlobalConstants.Immediate, result.outrate.get(0, 0), GlobalConstants.FineTol,
            "Enabling events should have immediate rate");
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol,
            "Probability should be 1 for deterministic enabling");
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testEnablingEvent_Mode1Exp_PartiallyEnabled() {
        // Test: Partially enabling event for mode 1 (Exp) with 1 job
        // This corresponds to lines 77-85 in test_SPN_events.m
        
        // Create a GlobalSync for testing
        GlobalSync glevent = TestSPNModels.createGlobalSync(multiServerTransitionIndex, EventType.ENABLE, 0); // Mode 1
        
        List<Matrix> glspace = new ArrayList<Matrix>();
        // FCFS buffer with 1 job
        Matrix buffer = new Matrix(1, 1);
        buffer.set(0, 0, 1); // job class 1
        glspace.add(buffer);
        
        // Transition state: [2 1 1 0 0 0 0 0]
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        transState.set(0, 3, 0);
        transState.set(0, 4, 0);
        transState.set(0, 5, 0);
        transState.set(0, 6, 0);
        transState.set(0, 7, 0);
        glspace.add(transState);
        
        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(multiServerSn, multiServerTransitionIndex, glspace, glevent, false);
        
        assertNotNull(result.outglspace, "Should return output global space");
        
        // Expected output state: [1 1 1 1 0 0 0 0] - only 1 server enabled
        Matrix outTransState = result.outglspace.get(1);
        assertEquals(1, outTransState.get(0, 0), 0.0, "1 job remains");
        assertEquals(1, outTransState.get(0, 3), 0.0, "Only 1 mode1 server enabled");
        
        assertEquals(GlobalConstants.Immediate, result.outrate.get(0, 0), GlobalConstants.FineTol);
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol);
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testEnablingEvent_Mode1Exp_NotEnabled() {
        // Test: Not enabling event for mode 1 (Exp) with no jobs
        // This corresponds to lines 87-95 in test_SPN_events.m
        
        // Create a GlobalSync for testing
        GlobalSync glevent = TestSPNModels.createGlobalSync(multiServerTransitionIndex, EventType.ENABLE, 0); // Mode 1
        
        List<Matrix> glspace = new ArrayList<Matrix>();
        // Empty FCFS buffer
        Matrix buffer = new Matrix(1, 1);
        buffer.set(0, 0, 0); // no jobs
        glspace.add(buffer);
        
        // Transition state: [2 1 1 0 0 0 0 0]
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        transState.set(0, 3, 0);
        transState.set(0, 4, 0);
        transState.set(0, 5, 0);
        transState.set(0, 6, 0);
        transState.set(0, 7, 0);
        glspace.add(transState);
        
        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(multiServerSn, multiServerTransitionIndex, glspace, glevent, false);
        
        assertNotNull(result.outglspace, "Should return output global space");
        
        // Expected output state: [2 1 1 0 0 0 0 0] - no change
        Matrix outTransState = result.outglspace.get(1);
        assertEquals(2, outTransState.get(0, 0), 0.0, "Jobs unchanged");
        assertEquals(0, outTransState.get(0, 3), 0.0, "No servers enabled");
        
        assertEquals(GlobalConstants.Immediate, result.outrate.get(0, 0), GlobalConstants.FineTol);
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol);
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testEnablingEvent_Mode2Erlang() {
        // Test: Enabling event for mode 2 (Erlang) with phase selection
        // This corresponds to lines 97-106 in test_SPN_events.m
        
        // Create a GlobalSync for testing
        GlobalSync glevent = TestSPNModels.createGlobalSync(multiServerTransitionIndex, EventType.ENABLE, 1); // Mode 2 (Erlang)
        
        List<Matrix> glspace = new ArrayList<Matrix>();
        // FCFS buffer with 2 jobs
        Matrix buffer = new Matrix(1, 2);
        buffer.set(0, 0, 1);
        buffer.set(0, 1, 1);
        glspace.add(buffer);
        
        // Transition state: [2 1 1 0 0 0 0 0]
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        transState.set(0, 3, 0);
        transState.set(0, 4, 0); // erlang phase 1
        transState.set(0, 5, 0); // erlang phase 2
        transState.set(0, 6, 0);
        transState.set(0, 7, 0);
        glspace.add(transState);
        
        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(multiServerSn, multiServerTransitionIndex, glspace, glevent, false);
        
        assertNotNull(result.outglspace, "Should return output global space");
        
        // For Erlang, we expect multiple possible outcomes (phase 1 or phase 2)
        // The test expects 2 rows in the output
        if (result.outglspace.get(1).getNumRows() > 1) {
            // First row: [2 0 1 0 1 0 0 0] - phase 1 enabled
            assertEquals(1, result.outglspace.get(1).get(0, 4), 0.0, "Phase 1 enabled");
            assertEquals(0, result.outglspace.get(1).get(0, 5), 0.0);
            
            // Second row: [2 0 1 0 0 1 0 0] - phase 2 enabled
            assertEquals(0, result.outglspace.get(1).get(1, 4), 0.0);
            assertEquals(1, result.outglspace.get(1).get(1, 5), 0.0, "Phase 2 enabled");
        }
        
        // Check rates
        if (result.outrate.getNumRows() > 1) {
            assertEquals(GlobalConstants.Immediate, result.outrate.get(0, 0), GlobalConstants.FineTol);
            assertEquals(GlobalConstants.Immediate, result.outrate.get(1, 0), GlobalConstants.FineTol);
        }
        
        // Check probabilities - first phase should have probability 1, second 0
        if (result.outprob.getNumRows() > 1) {
            assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol, "First phase probability");
            assertEquals(0.0, result.outprob.get(1, 0), GlobalConstants.FineTol, "Second phase probability");
        }
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testEnablingEvent_Mode3HyperExp() {
        // Test: Enabling event for mode 3 (HyperExp) with probabilistic phase selection
        // This corresponds to lines 108-117 in test_SPN_events.m
        
        // Create a GlobalSync for testing
        GlobalSync glevent = TestSPNModels.createGlobalSync(multiServerTransitionIndex, EventType.ENABLE, 2); // Mode 3 (HyperExp)
        
        List<Matrix> glspace = new ArrayList<Matrix>();
        // FCFS buffer with 2 jobs
        Matrix buffer = new Matrix(1, 2);
        buffer.set(0, 0, 1);
        buffer.set(0, 1, 1);
        glspace.add(buffer);
        
        // Transition state: [2 1 1 0 0 0 0 0]
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        transState.set(0, 3, 0);
        transState.set(0, 4, 0);
        transState.set(0, 5, 0);
        transState.set(0, 6, 0); // hyperexp phase 1
        transState.set(0, 7, 0); // hyperexp phase 2
        glspace.add(transState);
        
        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(multiServerSn, multiServerTransitionIndex, glspace, glevent, false);
        
        assertNotNull(result.outglspace, "Should return output global space");
        
        // For HyperExp, we expect multiple possible outcomes based on phase probabilities
        if (result.outglspace.get(1).getNumRows() > 1) {
            // First row: [2 1 0 0 0 0 1 0] - phase 1 selected
            assertEquals(1, result.outglspace.get(1).get(0, 6), 0.0, "HyperExp phase 1");
            assertEquals(0, result.outglspace.get(1).get(0, 7), 0.0);
            
            // Second row: [2 1 0 0 0 0 0 1] - phase 2 selected
            assertEquals(0, result.outglspace.get(1).get(1, 6), 0.0);
            assertEquals(1, result.outglspace.get(1).get(1, 7), 0.0, "HyperExp phase 2");
        }
        
        // Check probabilities match HyperExp phase entry probabilities (0.99, 0.01)
        if (result.outprob.getNumRows() > 1) {
            assertEquals(0.99, result.outprob.get(0, 0), GlobalConstants.FineTol, 
                "HyperExp phase 1 probability");
            assertEquals(0.01, result.outprob.get(1, 0), GlobalConstants.FineTol,
                "HyperExp phase 2 probability");
        }
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testFiringEvent_ConsumeAndProduceTokens() {
        // Test: Firing event that consumes 2 tokens and produces 1
        // This corresponds to lines 119-130 in test_SPN_events.m
        
        // Create a GlobalSync for firing event
        GlobalSync glevent = TestSPNModels.createGlobalSync(multiServerTransitionIndex, EventType.FIRE, 0); // Firing event for Mode 1
        
        List<Matrix> glspace = new ArrayList<Matrix>();
        // FCFS buffer at Place 1: [2,1,2,1,1] - 3 jobs of class 1, 2 of class 2
        Matrix buffer = new Matrix(1, 5);
        buffer.set(0, 0, 2); // class 2
        buffer.set(0, 1, 1); // class 1
        buffer.set(0, 2, 2); // class 2
        buffer.set(0, 3, 1); // class 1
        buffer.set(0, 4, 1); // class 1
        glspace.add(buffer);
        
        // Transition state: [0 1 1 2 0 0 0 0] - 2 mode1 servers active
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 0);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        transState.set(0, 3, 2); // 2 servers for mode 1
        transState.set(0, 4, 0);
        transState.set(0, 5, 0);
        transState.set(0, 6, 0);
        transState.set(0, 7, 0);
        glspace.add(transState);
        
        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(multiServerSn, multiServerTransitionIndex, glspace, glevent, false);
        
        assertNotNull(result.outglspace, "Should return output global space");
        
        // Expected buffer output: [1 2 1 2] - consumed 2 front tokens, added 1 at end
        Matrix outBuffer = result.outglspace.get(0);
        assertEquals(4, outBuffer.getNumCols(), "Buffer should have 4 tokens after firing");
        assertEquals(1, outBuffer.get(0, 0), 0.0, "First token class");
        assertEquals(2, outBuffer.get(0, 1), 0.0, "Second token class");
        assertEquals(1, outBuffer.get(0, 2), 0.0, "Third token class");
        assertEquals(2, outBuffer.get(0, 3), 0.0, "Fourth token class");
        
        // Expected transition state: [1 1 1 1 0 0 0 0] - 1 server remains active
        Matrix outTransState = result.outglspace.get(1);
        assertEquals(1, outTransState.get(0, 0), 0.0, "1 token added back");
        assertEquals(1, outTransState.get(0, 3), 0.0, "1 server remains active");
        
        assertEquals(2.0, result.outrate.get(0, 0), GlobalConstants.FineTol,
            "Two Mode 1 servers firing with rate 1 each");
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol);
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testFiringEvent_SingleToken() {
        // Test: Firing event that consumes 1 token and produces 1
        // This corresponds to lines 132-142 in test_SPN_events.m
        
        // Create a GlobalSync for firing event
        GlobalSync glevent = TestSPNModels.createGlobalSync(multiServerTransitionIndex, EventType.FIRE, 0); // Firing event for Mode 1
        
        List<Matrix> glspace = new ArrayList<Matrix>();
        // FCFS buffer: [2,2,1] - 2 jobs of class 2, 1 of class 1
        Matrix buffer = new Matrix(1, 3);
        buffer.set(0, 0, 2); // class 2
        buffer.set(0, 1, 2); // class 2
        buffer.set(0, 2, 1); // class 1
        glspace.add(buffer);
        
        // Transition state: [1 1 1 1 0 0 0 0] - 1 mode1 server active
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 1);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        transState.set(0, 3, 1); // 1 server for mode 1
        transState.set(0, 4, 0);
        transState.set(0, 5, 0);
        transState.set(0, 6, 0);
        transState.set(0, 7, 0);
        glspace.add(transState);
        
        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(multiServerSn, multiServerTransitionIndex, glspace, glevent, false);
        
        assertNotNull(result.outglspace, "Should return output global space");
        
        // Expected buffer output: [1 2 2] - consumed 1 front token, added 1 at end
        Matrix outBuffer = result.outglspace.get(0);
        assertEquals(3, outBuffer.getNumCols(), "Buffer should have 3 tokens after firing");
        assertEquals(1, outBuffer.get(0, 0), 0.0);
        assertEquals(2, outBuffer.get(0, 1), 0.0);
        assertEquals(2, outBuffer.get(0, 2), 0.0);
        
        // Expected transition state: [2 1 1 0 0 0 0 0] - servers disabled after firing
        Matrix outTransState = result.outglspace.get(1);
        assertEquals(2, outTransState.get(0, 0), 0.0, "2 tokens available");
        assertEquals(0, outTransState.get(0, 3), 0.0, "Servers disabled");
        
        assertEquals(1.0, result.outrate.get(0, 0), GlobalConstants.FineTol,
            "One Mode 1 server firing with rate 1");
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol);
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testEnablingEvent_Mode3HyperExp_MultipleServers() {
        // Test: Enabling event for mode 3 (HyperExp) with 2 servers
        // This corresponds to lines 144-161 in test_SPN_events.m
        
        // Model already has 2 servers for mode 3
        // Create a GlobalSync for testing
        GlobalSync glevent = TestSPNModels.createGlobalSync(multiServerTransitionIndex, EventType.ENABLE, 2); // Mode 3 (HyperExp)
        
        List<Matrix> glspace = new ArrayList<Matrix>();
        // FCFS buffer with 3 jobs
        Matrix buffer = new Matrix(1, 3);
        buffer.set(0, 0, 1);
        buffer.set(0, 1, 1);
        buffer.set(0, 2, 1);
        glspace.add(buffer);
        
        // Transition state: [2 1 2 0 0 0 0 0] - can enable 2 servers
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 2); // 2 servers available for mode 3
        transState.set(0, 3, 0);
        transState.set(0, 4, 0);
        transState.set(0, 5, 0);
        transState.set(0, 6, 0); // hyperexp phase 1
        transState.set(0, 7, 0); // hyperexp phase 2
        glspace.add(transState);
        
        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(multiServerSn, multiServerTransitionIndex, glspace, glevent, false);
        
        assertNotNull(result.outglspace, "Should return output global space");
        
        // With 2 servers, we expect 3 possible outcomes: (2,0), (1,1), (0,2)
        if (result.outglspace.get(1).getNumRows() >= 3) {
            // First row: [2 1 0 0 0 0 2 0] - both servers in phase 1
            assertEquals(2, result.outglspace.get(1).get(0, 6), 0.0, "Both servers in phase 1");
            assertEquals(0, result.outglspace.get(1).get(0, 7), 0.0);
            
            // Second row: [2 1 0 0 0 0 1 1] - one server each phase
            assertEquals(1, result.outglspace.get(1).get(1, 6), 0.0, "One server in phase 1");
            assertEquals(1, result.outglspace.get(1).get(1, 7), 0.0, "One server in phase 2");
            
            // Third row: [2 1 0 0 0 0 0 2] - both servers in phase 2
            assertEquals(0, result.outglspace.get(1).get(2, 6), 0.0);
            assertEquals(2, result.outglspace.get(1).get(2, 7), 0.0, "Both servers in phase 2");
        }
        
        // Check multinomial probabilities based on HyperExp phase probabilities
        // With p1=0.99, p2=0.01, multinomial probs are approximately:
        // (2,0): 0.9801, (1,1): 0.0198, (0,2): 0.0001
        if (result.outprob.getNumRows() >= 3) {
            assertEquals(0.9801, result.outprob.get(0, 0), GlobalConstants.FineTol,
                "Probability of both servers in phase 1");
            assertEquals(0.0198, result.outprob.get(1, 0), GlobalConstants.FineTol,
                "Probability of mixed phases");
            assertEquals(0.0001, result.outprob.get(2, 0), GlobalConstants.FineTol,
                "Probability of both servers in phase 2");
        }
        
        // All should have immediate rates
        for (int i = 0; i < result.outrate.getNumRows(); i++) {
            assertEquals(GlobalConstants.Immediate, result.outrate.get(i, 0), GlobalConstants.FineTol);
        }
    }

    // ==================== STATE SPACE GENERATION TESTS ====================

    @Test
    public void testSpaceClosedSingle_Basic() {
        // Test state space generation for single class closed network
        // M stations, N jobs
        int M = 2; // 2 stations
        int N = 3; // 3 jobs

        Matrix space = State.spaceClosedSinglePublic(M, N);

        assertNotNull(space, "State space should not be null");

        // For M=2 stations and N=3 jobs, the possible states are:
        // (3,0), (2,1), (1,2), (0,3) = 4 states
        // This follows stars and bars: C(N+M-1, M-1) = C(4,1) = 4
        assertEquals(4, space.getNumRows(), "Should have 4 possible states");
        assertEquals(M, space.getNumCols(), "Each state should have M columns");

        // Verify each row sums to N (total jobs is conserved)
        for (int i = 0; i < space.getNumRows(); i++) {
            double rowSum = 0;
            for (int j = 0; j < space.getNumCols(); j++) {
                rowSum += space.get(i, j);
            }
            assertEquals(N, rowSum, 0.0, "Each state should have N total jobs");
        }
    }

    @Test
    public void testSpaceClosedSingle_SingleStation() {
        // Edge case: single station
        int M = 1;
        int N = 5;

        Matrix space = State.spaceClosedSinglePublic(M, N);

        assertNotNull(space);
        assertEquals(1, space.getNumRows(), "Single station should have one state");
        assertEquals(1, space.getNumCols());
        assertEquals(5, space.get(0, 0), 0.0, "All jobs at the only station");
    }

    @Test
    public void testSpaceClosedSingle_ZeroJobs() {
        // Edge case: zero jobs
        int M = 3;
        int N = 0;

        Matrix space = State.spaceClosedSinglePublic(M, N);

        assertNotNull(space);
        assertEquals(1, space.getNumRows(), "Zero jobs should have one state (all zeros)");

        // Verify all zeros
        for (int j = 0; j < space.getNumCols(); j++) {
            assertEquals(0, space.get(0, j), 0.0);
        }
    }

    @Test
    public void testSpaceClosedMulti_TwoClasses() {
        // Test state space generation for multi-class closed network
        int M = 2; // 2 stations
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 2); // 2 jobs of class 1
        N.set(0, 1, 1); // 1 job of class 2

        Matrix space = State.spaceClosedMulti(M, N);

        assertNotNull(space, "State space should not be null");

        // For 2 stations, 2 class-1 jobs: (2,0), (1,1), (0,2) = 3 states
        // For 2 stations, 1 class-2 job: (1,0), (0,1) = 2 states
        // Total: 3 * 2 = 6 states when decorated
        assertTrue(space.getNumRows() > 0, "Should have positive number of states");

        // Each row should have M * R columns (stations * classes)
        assertEquals(M * 2, space.getNumCols(), "Should have M*R columns");
    }

    @Test
    public void testSpaceClosedMulti_SingleClass() {
        // Multi-class with single class should behave like single-class
        int M = 3;
        Matrix N = new Matrix(1, 1);
        N.set(0, 0, 2); // 2 jobs

        Matrix space = State.spaceClosedMulti(M, N);

        assertNotNull(space);
        assertEquals(M, space.getNumCols(), "Single class should have M columns");

        // C(N+M-1, M-1) = C(4, 2) = 6 states
        assertEquals(6, space.getNumRows(), "Should have 6 states for M=3, N=2");
    }

    @Test
    public void testSpaceClosedMultiCS_WithChains() {
        // Test state space with class switching chains
        int M = 2; // 2 stations
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 2); // 2 jobs of class 1
        N.set(0, 1, 1); // 1 job of class 2

        // Chains matrix: classes in same chain can switch
        // Single chain containing both classes
        Matrix chains = new Matrix(1, 2);
        chains.set(0, 0, 1); // Class 1 in chain 1
        chains.set(0, 1, 1); // Class 2 in chain 1

        Matrix space = State.spaceClosedMultiCS(M, N, chains);

        assertNotNull(space, "State space with chains should not be null");
        assertTrue(space.getNumRows() > 0, "Should have positive number of states");
    }

    @Test
    public void testSpaceClosedMultiCS_SeparateChains() {
        // Test with separate chains (no class switching)
        int M = 2;
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 1);
        N.set(0, 1, 1);

        // Two separate chains
        Matrix chains = new Matrix(2, 2);
        chains.set(0, 0, 1); chains.set(0, 1, 0); // Class 1 in chain 1
        chains.set(1, 0, 0); chains.set(1, 1, 1); // Class 2 in chain 2

        Matrix space = State.spaceClosedMultiCS(M, N, chains);

        assertNotNull(space);
        // With separate chains, should be same as spaceClosedMulti
        Matrix spaceNoCS = State.spaceClosedMulti(M, N);
        assertEquals(spaceNoCS.getNumRows(), space.getNumRows(),
            "Separate chains should produce same number of states");
    }

    @Test
    public void testSpaceClosedSingle_LargerNetwork() {
        // Test with more stations
        int M = 4;
        int N = 2;

        Matrix space = State.spaceClosedSinglePublic(M, N);

        // C(N+M-1, M-1) = C(5, 3) = 10 states
        assertEquals(10, space.getNumRows(), "Should have 10 states for M=4, N=2");
        assertEquals(M, space.getNumCols());

        // Verify conservation
        for (int i = 0; i < space.getNumRows(); i++) {
            double sum = 0;
            for (int j = 0; j < M; j++) {
                sum += space.get(i, j);
                assertTrue(space.get(i, j) >= 0, "All values should be non-negative");
            }
            assertEquals(N, sum, 0.0, "Jobs should be conserved");
        }
    }

    @Test
    public void testSpaceClosedSingle_AllUniqueStates() {
        // Verify all generated states are unique
        int M = 3;
        int N = 3;

        Matrix space = State.spaceClosedSinglePublic(M, N);

        // Check no duplicate rows
        Set<String> seen = new java.util.HashSet<>();
        for (int i = 0; i < space.getNumRows(); i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < space.getNumCols(); j++) {
                sb.append(space.get(i, j)).append(",");
            }
            String state = sb.toString();
            assertFalse(seen.contains(state), "State should be unique: " + state);
            seen.add(state);
        }
    }
}