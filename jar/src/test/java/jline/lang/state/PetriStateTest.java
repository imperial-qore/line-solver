package jline.lang.state;

import jline.io.Ret;
import jline.lang.GlobalSync;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.Sync;
import jline.lang.constant.EventType;
import jline.GlobalConstants;
import jline.lang.nodes.Transition;
import jline.lang.state.AfterGlobalEvent.AfterGlobalEventResult;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for State.afterEvent and State.afterGlobalEvent methods,
 * replicating every check in MATLAB test_SPN_events.m.
 *
 * Model: spn_basic_closed (P1 -> T1 -> P1, 3 modes: Exp, Erlang, HyperExp)
 *
 * Sync ordering (from model, 0-based):
 *   sync[0]: DEP P1 -> ARV T1
 *   sync[1]: PHASE T1 mode0 (Exp)
 *   sync[2]: PHASE T1 mode1 (Erlang)
 *   sync[3]: PHASE T1 mode2 (HyperExp)
 *   sync[4]: DEP T1 -> ARV P1
 *
 * GSync ordering (from model, 0-based, enables first then fires):
 *   gsync[0]: ENABLE mode0 (Exp)
 *   gsync[1]: ENABLE mode1 (Erlang)
 *   gsync[2]: ENABLE mode2 (HyperExp)
 *   gsync[3]: FIRE mode0 (Exp)
 *   gsync[4]: FIRE mode1 (Erlang)
 *   gsync[5]: FIRE mode2 (HyperExp)
 */
public class PetriStateTest {

    // For afterEvent tests (checks 1-5): all modes 1 server
    private Network model;
    private NetworkStruct sn;
    private int transitionIndex;

    // For afterGlobalEvent tests (checks 6-12): mode1=2 servers, mode2=1, mode3=1
    private Network mode1TwoServerModel;
    private NetworkStruct mode1TwoServerSn;
    private int mode1TwoServerTransitionIndex;

    // For afterGlobalEvent check 13: mode1=2, mode2=1, mode3=2
    private Network mode1Mode3TwoServerModel;
    private NetworkStruct mode1Mode3TwoServerSn;
    private int mode1Mode3TwoServerTransitionIndex;

    @BeforeEach
    public void setUp() {
        // Model for afterEvent tests (checks 1-5)
        model = PetriStateTestFixtures.createTestModelWithThreeModes();
        sn = model.getStruct();
        Transition T1 = (Transition) model.getNodeByName("T1");
        transitionIndex = T1.getNodeIndex();

        // Model for afterGlobalEvent checks 6-12 (mode1=2 servers)
        mode1TwoServerModel = PetriStateTestFixtures.createTestModelMode1TwoServers();
        mode1TwoServerSn = mode1TwoServerModel.getStruct();
        Transition mode1T1 = (Transition) mode1TwoServerModel.getNodeByName("T1");
        mode1TwoServerTransitionIndex = mode1T1.getNodeIndex();

        // Model for afterGlobalEvent check 13 (mode1=2, mode3=2 servers)
        mode1Mode3TwoServerModel = PetriStateTestFixtures.createTestModelMode1Mode3TwoServers();
        mode1Mode3TwoServerSn = mode1Mode3TwoServerModel.getStruct();
        Transition mode13T1 = (Transition) mode1Mode3TwoServerModel.getNodeByName("T1");
        mode1Mode3TwoServerTransitionIndex = mode13T1.getNodeIndex();
    }

    // ==================== STATE AFTER EVENT TESTS (checks 1-5) ====================

    @Test
    public void testPhaseTransition_ExponentialMode_AllDisabled() {
        // Check 1 (MATLAB lines 19-25): s=2, all modes disabled -> empty outspace
        Sync sync = sn.sync.get(1); // PHASE T1 mode0 (Exp)
        EventType event = sync.active.get(0).getEvent();
        int jobClass = sync.active.get(0).getJobClass();

        // State: all modes disabled [1 1 1 0 0 0 0 0]
        Matrix inspace = new Matrix(1, 8);
        inspace.set(0, 0, 1);
        inspace.set(0, 1, 1);
        inspace.set(0, 2, 1);

        Ret.EventResult result = State.afterEvent(sn, transitionIndex, inspace, event, jobClass, false);

        assertTrue(result.outspace.isEmpty(), "Output space should be empty when all modes are disabled");
    }

    @Test
    public void testPhaseTransition_ExponentialMode_CannotProgress() {
        // Check 2 (MATLAB lines 27-31): s=2, Exp enabled in phase 1 -> empty (only 1 phase)
        Sync sync = sn.sync.get(1);
        EventType event = sync.active.get(0).getEvent();
        int jobClass = sync.active.get(0).getJobClass();

        // State: mode 1 enabled [0 1 1 1 0 0 0 0]
        Matrix inspace = new Matrix(1, 8);
        inspace.set(0, 1, 1);
        inspace.set(0, 2, 1);
        inspace.set(0, 3, 1);

        Ret.EventResult result = State.afterEvent(sn, transitionIndex, inspace, event, jobClass, false);

        assertTrue(result.outspace.isEmpty(), "Exponential mode should not have phase transitions");
    }

    @Test
    public void testPhaseTransition_ErlangMode_Phase1ToPhase2() {
        // Check 3 (MATLAB lines 33-40): s=3, Erlang phase 1 -> phase 2, rate=2
        Sync sync = sn.sync.get(2); // PHASE T1 mode1 (Erlang)
        EventType event = sync.active.get(0).getEvent();
        int jobClass = sync.active.get(0).getJobClass();

        // State: mode 2 phase 1 enabled [1 0 1 0 1 0 0 0]
        Matrix inspace = new Matrix(1, 8);
        inspace.set(0, 0, 1);
        inspace.set(0, 2, 1);
        inspace.set(0, 4, 1);

        Ret.EventResult result = State.afterEvent(sn, transitionIndex, inspace, event, jobClass, false);

        assertFalse(result.outspace.isEmpty(), "Should return output space for phase transition");

        // Expected: [1 0 1 0 0 1 0 0]
        double[] expected = {1, 0, 1, 0, 0, 1, 0, 0};
        for (int i = 0; i < 8; i++) {
            assertEquals(expected[i], result.outspace.get(0, i), 0.0,
                    "outspace[" + i + "] mismatch");
        }

        assertEquals(2.0, result.outrate.get(0, 0), GlobalConstants.FineTol,
                "Erlang rate should be 2.0");
    }

    @Test
    public void testPhaseTransition_ErlangMode_Phase2CannotProgress() {
        // Check 4 (MATLAB lines 42-48): s=3, Erlang in phase 2 -> rate=0
        Sync sync = sn.sync.get(2);
        EventType event = sync.active.get(0).getEvent();
        int jobClass = sync.active.get(0).getJobClass();

        // State: mode 2 phase 2 enabled [1 0 1 0 0 1 0 0]
        Matrix inspace = new Matrix(1, 8);
        inspace.set(0, 0, 1);
        inspace.set(0, 2, 1);
        inspace.set(0, 5, 1);

        Ret.EventResult result = State.afterEvent(sn, transitionIndex, inspace, event, jobClass, false);

        assertNotNull(result.outrate, "Should return outrate");
        assertEquals(0.0, result.outrate.get(0, 0), GlobalConstants.FineTol,
                "Rate should be 0 for Erlang phase 2 (final phase)");
    }

    @Test
    public void testPhaseTransition_HyperExpMode_Phase1RateZero() {
        // Check 5 (MATLAB lines 50-56): s=4, HyperExp phase 1 -> rate=0
        Sync sync = sn.sync.get(3); // PHASE T1 mode2 (HyperExp)
        EventType event = sync.active.get(0).getEvent();
        int jobClass = sync.active.get(0).getJobClass();

        // State: mode 3 phase 1 enabled [1 0 0 0 1 0 1 0]
        Matrix inspace = new Matrix(1, 8);
        inspace.set(0, 0, 1);
        inspace.set(0, 4, 1);
        inspace.set(0, 6, 1);

        Ret.EventResult result = State.afterEvent(sn, transitionIndex, inspace, event, jobClass, false);

        assertNotNull(result.outrate, "Should return rate for HyperExp phase transition");
        assertEquals(0.0, result.outrate.get(0, 0), GlobalConstants.FineTol,
                "HyperExp phase transition rate should be 0");
    }

    // ==================== STATE AFTER GLOBAL EVENT TESTS (checks 6-13) ====================

    @Test
    public void testEnablingEvent_Mode1Exp_FullyEnabled() {
        // Check 6 (MATLAB lines 67-75): gsync{1} Enable Mode1 Exp, 2 jobs -> 2 servers enabled
        GlobalSync glevent = mode1TwoServerSn.gsync.get(0); // ENABLE mode0

        List<Matrix> glspace = new ArrayList<Matrix>();
        Matrix buffer = new Matrix(1, 2);
        buffer.set(0, 0, 1);
        buffer.set(0, 1, 1);
        glspace.add(buffer);

        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        glspace.add(transState);

        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(
                mode1TwoServerSn, mode1TwoServerTransitionIndex, glspace, glevent, false);

        assertNotNull(result.outglspace, "Should return output global space");
        assertEquals(2, result.outglspace.size(), "Should have 2 elements in output");

        // Expected transition state: [0 1 1 2 0 0 0 0]
        Matrix outTransState = result.outglspace.get(1);
        double[] expectedTrans = {0, 1, 1, 2, 0, 0, 0, 0};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedTrans[i], outTransState.get(0, i), 0.0,
                    "outglspace{2}[" + i + "] mismatch");
        }

        assertEquals(GlobalConstants.Immediate, result.outrate.get(0, 0), GlobalConstants.FineTol,
                "Enabling events should have immediate rate");
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol,
                "Probability should be 1");
    }

    @Test
    public void testEnablingEvent_Mode1Exp_PartiallyEnabled() {
        // Check 7 (MATLAB lines 77-85): gsync{1} Enable Mode1 Exp, 1 job -> 1 server enabled
        GlobalSync glevent = mode1TwoServerSn.gsync.get(0);

        List<Matrix> glspace = new ArrayList<Matrix>();
        Matrix buffer = new Matrix(1, 1);
        buffer.set(0, 0, 1);
        glspace.add(buffer);

        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        glspace.add(transState);

        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(
                mode1TwoServerSn, mode1TwoServerTransitionIndex, glspace, glevent, false);

        assertNotNull(result.outglspace, "Should return output global space");

        // Expected transition state: [1 1 1 1 0 0 0 0]
        Matrix outTransState = result.outglspace.get(1);
        double[] expectedTrans = {1, 1, 1, 1, 0, 0, 0, 0};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedTrans[i], outTransState.get(0, i), 0.0,
                    "outglspace{2}[" + i + "] mismatch");
        }

        assertEquals(GlobalConstants.Immediate, result.outrate.get(0, 0), GlobalConstants.FineTol);
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol);
    }

    @Test
    public void testEnablingEvent_Mode1Exp_NotEnabled() {
        // Check 8 (MATLAB lines 87-95): gsync{1} Enable Mode1 Exp, 0 jobs -> no change
        GlobalSync glevent = mode1TwoServerSn.gsync.get(0);

        List<Matrix> glspace = new ArrayList<Matrix>();
        Matrix buffer = new Matrix(1, 1);
        buffer.set(0, 0, 0);
        glspace.add(buffer);

        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        glspace.add(transState);

        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(
                mode1TwoServerSn, mode1TwoServerTransitionIndex, glspace, glevent, false);

        assertNotNull(result.outglspace, "Should return output global space");

        // Expected transition state: [2 1 1 0 0 0 0 0] - no change
        Matrix outTransState = result.outglspace.get(1);
        double[] expectedTrans = {2, 1, 1, 0, 0, 0, 0, 0};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedTrans[i], outTransState.get(0, i), 0.0,
                    "outglspace{2}[" + i + "] mismatch");
        }

        assertEquals(GlobalConstants.Immediate, result.outrate.get(0, 0), GlobalConstants.FineTol);
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol);
    }

    @Test
    public void testEnablingEvent_Mode2Erlang() {
        // Check 9 (MATLAB lines 97-106): gsync{2} Enable Mode2 Erlang, 2 phase outcomes
        GlobalSync glevent = mode1TwoServerSn.gsync.get(1); // ENABLE mode1 (Erlang)

        List<Matrix> glspace = new ArrayList<Matrix>();
        Matrix buffer = new Matrix(1, 2);
        buffer.set(0, 0, 1);
        buffer.set(0, 1, 1);
        glspace.add(buffer);

        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        glspace.add(transState);

        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(
                mode1TwoServerSn, mode1TwoServerTransitionIndex, glspace, glevent, false);

        assertNotNull(result.outglspace, "Should return output global space");

        // Expected: 2 rows in transition state
        Matrix outTransState = result.outglspace.get(1);
        assertEquals(2, outTransState.getNumRows(),
                "Should have 2 possible phase outcomes for Erlang");

        // Row 0: [2 0 1 0 1 0 0 0] - phase 1 enabled
        double[] expectedRow0 = {2, 0, 1, 0, 1, 0, 0, 0};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedRow0[i], outTransState.get(0, i), 0.0,
                    "Row 0 col " + i + " mismatch");
        }

        // Row 1: [2 0 1 0 0 1 0 0] - phase 2 enabled
        double[] expectedRow1 = {2, 0, 1, 0, 0, 1, 0, 0};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedRow1[i], outTransState.get(1, i), 0.0,
                    "Row 1 col " + i + " mismatch");
        }

        // Rates: both Immediate
        assertEquals(2, result.outrate.getNumRows(), "Should have 2 rates");
        assertEquals(GlobalConstants.Immediate, result.outrate.get(0, 0), GlobalConstants.FineTol);
        assertEquals(GlobalConstants.Immediate, result.outrate.get(1, 0), GlobalConstants.FineTol);

        // Probabilities: [1; 0] (Erlang deterministically starts in phase 1)
        assertEquals(2, result.outprob.getNumRows(), "Should have 2 probabilities");
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol, "Phase 1 probability");
        assertEquals(0.0, result.outprob.get(1, 0), GlobalConstants.FineTol, "Phase 2 probability");
    }

    @Test
    public void testEnablingEvent_Mode3HyperExp() {
        // Check 10 (MATLAB lines 108-117): gsync{3} Enable Mode3 HyperExp, prob=map_pie
        GlobalSync glevent = mode1TwoServerSn.gsync.get(2); // ENABLE mode2 (HyperExp)

        List<Matrix> glspace = new ArrayList<Matrix>();
        Matrix buffer = new Matrix(1, 2);
        buffer.set(0, 0, 1);
        buffer.set(0, 1, 1);
        glspace.add(buffer);

        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        glspace.add(transState);

        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(
                mode1TwoServerSn, mode1TwoServerTransitionIndex, glspace, glevent, false);

        assertNotNull(result.outglspace, "Should return output global space");

        // Expected: 2 rows in transition state
        Matrix outTransState = result.outglspace.get(1);
        assertEquals(2, outTransState.getNumRows(),
                "Should have 2 possible phase outcomes for HyperExp");

        // Row 0: [2 1 0 0 0 0 1 0] - phase 1 selected
        double[] expectedRow0 = {2, 1, 0, 0, 0, 0, 1, 0};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedRow0[i], outTransState.get(0, i), 0.0,
                    "Row 0 col " + i + " mismatch");
        }

        // Row 1: [2 1 0 0 0 0 0 1] - phase 2 selected
        double[] expectedRow1 = {2, 1, 0, 0, 0, 0, 0, 1};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedRow1[i], outTransState.get(1, i), 0.0,
                    "Row 1 col " + i + " mismatch");
        }

        // Rates: both Immediate
        assertEquals(2, result.outrate.getNumRows(), "Should have 2 rates");
        assertEquals(GlobalConstants.Immediate, result.outrate.get(0, 0), GlobalConstants.FineTol);
        assertEquals(GlobalConstants.Immediate, result.outrate.get(1, 0), GlobalConstants.FineTol);

        // Probabilities: map_pie of HyperExp(1,4) = [0.99, 0.01]
        assertEquals(2, result.outprob.getNumRows(), "Should have 2 probabilities");
        assertEquals(0.99, result.outprob.get(0, 0), GlobalConstants.FineTol,
                "HyperExp phase 1 probability");
        assertEquals(0.01, result.outprob.get(1, 0), GlobalConstants.FineTol,
                "HyperExp phase 2 probability");
    }

    @Test
    public void testFiringEvent_ConsumeAndProduceTokens() {
        // Check 11 (MATLAB lines 119-130): gsync{4} Fire Mode1, consume 2 produce 1
        GlobalSync glevent = mode1TwoServerSn.gsync.get(3); // FIRE mode0
        int ind = glevent.getActive().get(0).getNode();

        List<Matrix> glspace = new ArrayList<Matrix>();
        // FCFS buffer: [2,1,2,1,1]
        Matrix buffer = new Matrix(1, 5);
        buffer.set(0, 0, 2);
        buffer.set(0, 1, 1);
        buffer.set(0, 2, 2);
        buffer.set(0, 3, 1);
        buffer.set(0, 4, 1);
        glspace.add(buffer);

        // Transition state: [0 1 1 2 0 0 0 0] - 2 mode1 servers active
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        transState.set(0, 3, 2);
        glspace.add(transState);

        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(
                mode1TwoServerSn, ind, glspace, glevent, false);

        assertNotNull(result.outglspace, "Should return output global space");

        // Expected buffer: [1 2 1 2]
        Matrix outBuffer = result.outglspace.get(0);
        assertEquals(4, outBuffer.getNumCols(), "Buffer should have 4 tokens after firing");
        double[] expectedBuffer = {1, 2, 1, 2};
        for (int i = 0; i < 4; i++) {
            assertEquals(expectedBuffer[i], outBuffer.get(0, i), 0.0,
                    "Buffer[" + i + "] mismatch");
        }

        // Expected transition state: [1 1 1 1 0 0 0 0]
        Matrix outTransState = result.outglspace.get(1);
        double[] expectedTrans = {1, 1, 1, 1, 0, 0, 0, 0};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedTrans[i], outTransState.get(0, i), 0.0,
                    "outglspace{2}[" + i + "] mismatch");
        }

        assertEquals(2.0, result.outrate.get(0, 0), GlobalConstants.FineTol,
                "Two Mode 1 servers firing with rate 1 each");
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol);
    }

    @Test
    public void testFiringEvent_SingleToken() {
        // Check 12 (MATLAB lines 132-142): gsync{4} Fire Mode1, consume 1 produce 1
        GlobalSync glevent = mode1TwoServerSn.gsync.get(3); // FIRE mode0
        int ind = glevent.getActive().get(0).getNode();

        List<Matrix> glspace = new ArrayList<Matrix>();
        // FCFS buffer: [2,2,1]
        Matrix buffer = new Matrix(1, 3);
        buffer.set(0, 0, 2);
        buffer.set(0, 1, 2);
        buffer.set(0, 2, 1);
        glspace.add(buffer);

        // Transition state: [1 1 1 1 0 0 0 0] - 1 mode1 server active
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 1);
        transState.set(0, 1, 1);
        transState.set(0, 2, 1);
        transState.set(0, 3, 1);
        glspace.add(transState);

        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(
                mode1TwoServerSn, ind, glspace, glevent, false);

        assertNotNull(result.outglspace, "Should return output global space");

        // Expected buffer: [1 2 2]
        Matrix outBuffer = result.outglspace.get(0);
        assertEquals(3, outBuffer.getNumCols(), "Buffer should have 3 tokens after firing");
        double[] expectedBuffer = {1, 2, 2};
        for (int i = 0; i < 3; i++) {
            assertEquals(expectedBuffer[i], outBuffer.get(0, i), 0.0,
                    "Buffer[" + i + "] mismatch");
        }

        // Expected transition state: [2 1 1 0 0 0 0 0]
        Matrix outTransState = result.outglspace.get(1);
        double[] expectedTrans = {2, 1, 1, 0, 0, 0, 0, 0};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedTrans[i], outTransState.get(0, i), 0.0,
                    "outglspace{2}[" + i + "] mismatch");
        }

        assertEquals(1.0, result.outrate.get(0, 0), GlobalConstants.FineTol,
                "One Mode 1 server firing with rate 1");
        assertEquals(1.0, result.outprob.get(0, 0), GlobalConstants.FineTol);
    }

    @Test
    public void testEnablingEvent_Mode3HyperExp_MultipleServers() {
        // Check 13 (MATLAB lines 144-161): gsync{3} Enable Mode3 HyperExp with 2 servers
        // Uses model with mode1=2, mode3=2 servers
        GlobalSync glevent = mode1Mode3TwoServerSn.gsync.get(2); // ENABLE mode2 (HyperExp)

        List<Matrix> glspace = new ArrayList<Matrix>();
        Matrix buffer = new Matrix(1, 3);
        buffer.set(0, 0, 1);
        buffer.set(0, 1, 1);
        buffer.set(0, 2, 1);
        glspace.add(buffer);

        // Transition state: [2 1 2 0 0 0 0 0]
        Matrix transState = new Matrix(1, 8);
        transState.set(0, 0, 2);
        transState.set(0, 1, 1);
        transState.set(0, 2, 2);
        glspace.add(transState);

        AfterGlobalEventResult result = AfterGlobalEvent.afterGlobalEvent(
                mode1Mode3TwoServerSn, mode1Mode3TwoServerTransitionIndex, glspace, glevent, false);

        assertNotNull(result.outglspace, "Should return output global space");

        // Expected: 3 rows (multinomial with 2 servers, 2 phases)
        Matrix outTransState = result.outglspace.get(1);
        assertEquals(3, outTransState.getNumRows(),
                "Should have 3 possible outcomes for 2 HyperExp servers: (2,0), (1,1), (0,2)");

        // Row 0: [2 1 0 0 0 0 2 0] - both servers in phase 1
        double[] expectedRow0 = {2, 1, 0, 0, 0, 0, 2, 0};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedRow0[i], outTransState.get(0, i), 0.0,
                    "Row 0 col " + i + " mismatch");
        }

        // Row 1: [2 1 0 0 0 0 1 1] - one server each phase
        double[] expectedRow1 = {2, 1, 0, 0, 0, 0, 1, 1};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedRow1[i], outTransState.get(1, i), 0.0,
                    "Row 1 col " + i + " mismatch");
        }

        // Row 2: [2 1 0 0 0 0 0 2] - both servers in phase 2
        double[] expectedRow2 = {2, 1, 0, 0, 0, 0, 0, 2};
        for (int i = 0; i < 8; i++) {
            assertEquals(expectedRow2[i], outTransState.get(2, i), 0.0,
                    "Row 2 col " + i + " mismatch");
        }

        // Rates: all Immediate
        assertEquals(3, result.outrate.getNumRows(), "Should have 3 rates");
        for (int i = 0; i < 3; i++) {
            assertEquals(GlobalConstants.Immediate, result.outrate.get(i, 0), GlobalConstants.FineTol);
        }

        // Multinomial probabilities with p1=0.99, p2=0.01:
        // (2,0): 0.9801, (1,1): 0.0198, (0,2): 0.0001
        assertEquals(3, result.outprob.getNumRows(), "Should have 3 probabilities");
        assertEquals(0.9801, result.outprob.get(0, 0), GlobalConstants.FineTol,
                "Probability of both servers in phase 1");
        assertEquals(0.0198, result.outprob.get(1, 0), GlobalConstants.FineTol,
                "Probability of mixed phases");
        assertEquals(0.0001, result.outprob.get(2, 0), GlobalConstants.FineTol,
                "Probability of both servers in phase 2");
    }
}
