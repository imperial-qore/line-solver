package jline.lang.state;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.lang.*;
import jline.lang.constant.EventType;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.DiscreteDistribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.processes.Zipf;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for Cache afterEvent with Random Replacement (RR) policy.
 *
 * Ground truth from MATLAB:
 *
 * n = 5, m=[3]
 * READ item 1 RR: [1,2,3,4,5] -> outstate=[3,2,1,4,5; 4,2,3,1,5; 5,2,3,4,1], outprob=[1/3; 1/3; 1/3]
 * READ item 3 RR: [1,2,3,4,5] -> outstate=[1,2,3,4,5], outprob = 1.0
 *
 * n = 5, m=[2,1]
 * READ item 1 RR: [1,2,3,4,5] -> outstate=[3,2,1,4,5; 4,2,3,1,5], outprob=[1/2; 1/2]
 * READ item 3 RR: [1,2,3,4,5] -> outstate=[1,2,5,4,3], outprob = 1.0
 * READ item 5 RR: [1,2,3,4,5] -> outstate=[1,2,3,4,5], outprob = 1.0
 *
 * n = 5, m=[1,2]
 * READ item 1 RR: [1,2,3,4,5] -> outstate=[3,2,1,4,5], outprob = 1.0
 * READ item 3 RR: [1,2,3,4,5] -> outstate=[1,2,4,3,5; 1,2,5,3,4], outprob=[1/2; 1/2]
 * READ item 5 RR: [1,2,3,4,5] -> outstate=[1,2,3,4,5], outprob = 1.0
 *
 * n = 5, m=[1,1,1]
 * READ item 1 RR: [1,2,3,4,5] -> outstate=[3,2,1,4,5], outprob = 1.0
 * READ item 3 RR: [1,2,3,4,5] -> outstate=[1,2,4,3,5], outprob = 1.0
 * READ item 4 RR: [1,2,3,4,5] -> outstate=[1,2,3,5,4], outprob = 1.0
 * READ item 5 RR: [1,2,3,4,5] -> outstate=[1,2,3,4,5], outprob = 1.0
 *
 * Java State Representation:
 * - spaceVar has sum(m) elements containing only cached item IDs
 * - Java may generate outcomes for all lists, filtered by accessCost (zero rate for inaccessible)
 */
public class CacheAfterEventTest {

    private static final double TOL = GlobalConstants.FineTol;

    // ==================== n=5, m=[3] ====================

    @Test
    public void testRR_n5_m3_readItem1_miss() {
        // Java: [3,4,5] -> [1,4,5; 3,1,5; 3,4,1], prob=[1/3; 1/3; 1/3]
        TestCacheModel model = createCacheModel(5, new int[]{3}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output states for cache miss");
        assertEquals(3, result.outspace.getNumRows(), "Should have 3 possible outcomes for m=[3]");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{1, 4, 5});
        assertCacheStateRow(result.outspace, 1, varStart, new int[]{3, 1, 5});
        assertCacheStateRow(result.outspace, 2, varStart, new int[]{3, 4, 1});

        // Verify equal probabilities
        double totalRate = result.outrate.elementSum();
        assertTrue(totalRate > 0, "Total rate should be positive");
        for (int i = 0; i < 3; i++) {
            assertEquals(1.0/3.0, result.outrate.get(i, 0) / totalRate, TOL,
                "Each outcome should have probability 1/3");
        }
    }

    @Test
    public void testRR_n5_m3_readItem3_hit() {
        // Java: [3,4,5] -> [3,4,5], prob=1.0
        TestCacheModel model = createCacheModel(5, new int[]{3}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 2); // item 3

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 outcome for hit in highest list");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5});
    }

    @Test
    public void testRR_n5_m3_readItem5_hit() {
        TestCacheModel model = createCacheModel(5, new int[]{3}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 4); // item 5

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 outcome for hit in highest list");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5});
    }

    // ==================== n=5, m=[2,1] ====================

    @Test
    public void testRR_n5_m21_readItem1_miss() {
        // Java: [3,4,5] -> [1,4,5; 3,1,5; ...], 2 non-zero rate outcomes
        TestCacheModel model = createCacheModel(5, new int[]{2, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output states for cache miss");
        assertEquals(3, result.outspace.getNumRows(), "Should have 3 total outcomes");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{1, 4, 5});
        assertCacheStateRow(result.outspace, 1, varStart, new int[]{3, 1, 5});

        int nonZeroCount = countNonZeroRates(result.outrate);
        assertEquals(2, nonZeroCount, "Should have 2 outcomes with non-zero rate");
    }

    @Test
    public void testRR_n5_m21_readItem3_hitList0() {
        // Java: [3,4,5] -> [5,4,3], 1 non-zero rate outcome
        TestCacheModel model = createCacheModel(5, new int[]{2, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 2); // item 3

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertTrue(result.outspace.getNumRows() >= 1, "Should have at least 1 outcome");

        int nonZeroCount = countNonZeroRates(result.outrate);
        assertEquals(1, nonZeroCount, "Should have 1 outcome with non-zero rate");

        int varStart = model.numClasses;
        for (int i = 0; i < result.outspace.getNumRows(); i++) {
            if (result.outrate.get(i, 0) > 0) {
                assertCacheStateRow(result.outspace, i, varStart, new int[]{5, 4, 3});
            }
        }
    }

    @Test
    public void testRR_n5_m21_readItem5_hitList1() {
        // Java: [3,4,5] -> [3,4,5], prob=1.0
        TestCacheModel model = createCacheModel(5, new int[]{2, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 4); // item 5

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 outcome for hit in highest list");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5});
    }

    // ==================== n=5, m=[1,2] ====================

    @Test
    public void testRR_n5_m12_readItem1_miss() {
        // Java: [3,4,5] -> [1,4,5; ...], 1 non-zero rate outcome
        TestCacheModel model = createCacheModel(5, new int[]{1, 2}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache miss");
        assertEquals(3, result.outspace.getNumRows(), "Should have 3 total outcomes");

        int nonZeroCount = countNonZeroRates(result.outrate);
        assertEquals(1, nonZeroCount, "Should have 1 outcome with non-zero rate for m[0]=1");

        int varStart = model.numClasses;
        assertTrue(result.outrate.get(0, 0) > 0, "First outcome should have non-zero rate");
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{1, 4, 5});
    }

    @Test
    public void testRR_n5_m12_readItem3_hitList0() {
        // Java: [3,4,5] -> [4,3,5; 5,4,3], 2 non-zero rate outcomes
        TestCacheModel model = createCacheModel(5, new int[]{1, 2}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 2); // item 3

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output states for cache hit with promotion");
        assertTrue(result.outspace.getNumRows() >= 2, "Should have at least 2 outcomes");

        int nonZeroCount = countNonZeroRates(result.outrate);
        assertEquals(2, nonZeroCount, "Should have 2 outcomes with non-zero rate for m[1]=2");

        int varStart = model.numClasses;
        int foundCount = 0;
        for (int i = 0; i < result.outspace.getNumRows(); i++) {
            if (result.outrate.get(i, 0) > 0) {
                if (foundCount == 0) {
                    assertCacheStateRow(result.outspace, i, varStart, new int[]{4, 3, 5});
                } else {
                    assertCacheStateRow(result.outspace, i, varStart, new int[]{5, 4, 3});
                }
                foundCount++;
            }
        }
    }

    @Test
    public void testRR_n5_m12_readItem5_hitList1() {
        // Java: [3,4,5] -> [3,4,5], prob=1.0
        TestCacheModel model = createCacheModel(5, new int[]{1, 2}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 4); // item 5

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 outcome for hit in highest list");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5});
    }

    // ==================== n=5, m=[1,1,1] ====================

    @Test
    public void testRR_n5_m111_readItem1_miss() {
        // Java: [3,4,5] -> [1,4,5; ...], 1 non-zero rate outcome
        TestCacheModel model = createCacheModel(5, new int[]{1, 1, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache miss");
        assertEquals(3, result.outspace.getNumRows(), "Should have 3 total outcomes");

        int nonZeroCount = countNonZeroRates(result.outrate);
        assertEquals(1, nonZeroCount, "Should have 1 outcome with non-zero rate for m[0]=1");

        int varStart = model.numClasses;
        assertTrue(result.outrate.get(0, 0) > 0, "First outcome should have non-zero rate");
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{1, 4, 5});
    }

    @Test
    public void testRR_n5_m111_readItem3_hitList0() {
        // Java: [3,4,5] -> [4,3,5; ...], 1 non-zero rate outcome
        TestCacheModel model = createCacheModel(5, new int[]{1, 1, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 2); // item 3

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertTrue(result.outspace.getNumRows() >= 1, "Should have at least 1 outcome");

        int nonZeroCount = countNonZeroRates(result.outrate);
        assertEquals(1, nonZeroCount, "Should have 1 outcome with non-zero rate for m[1]=1");

        int varStart = model.numClasses;
        for (int i = 0; i < result.outspace.getNumRows(); i++) {
            if (result.outrate.get(i, 0) > 0) {
                assertCacheStateRow(result.outspace, i, varStart, new int[]{4, 3, 5});
            }
        }
    }

    @Test
    public void testRR_n5_m111_readItem4_hitList1() {
        // Java: [3,4,5] -> item 4 is at position 1 (list 1, middle list)
        // MATLAB expects promotion to list 2: [1,2,3,5,4] -> Java [3,5,4]
        TestCacheModel model = createCacheModel(5, new int[]{1, 1, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 3); // item 4

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertTrue(result.outspace.getNumRows() >= 1, "Should have at least 1 outcome");

        int nonZeroCount = countNonZeroRates(result.outrate);
        assertEquals(1, nonZeroCount, "Should have 1 outcome with non-zero rate for m[2]=1");

        // Verify promotion from list 1 to list 2: [3,4,5] -> [3,5,4]
        int varStart = model.numClasses;
        for (int i = 0; i < result.outspace.getNumRows(); i++) {
            if (result.outrate.get(i, 0) > 0) {
                assertCacheStateRow(result.outspace, i, varStart, new int[]{3, 5, 4});
            }
        }
    }

    @Test
    public void testRR_n5_m111_readItem5_hitList2() {
        // Java: [3,4,5] -> [3,4,5], prob=1.0
        TestCacheModel model = createCacheModel(5, new int[]{1, 1, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheState(model, new int[]{3, 4, 5});

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 4); // item 5

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 outcome for hit in highest list");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5});
    }

    // ============= Cache with Retrieval System ==============

    // =========== n=5, m=[3] with Retrieval System ===========
    @Test
    public void testRS_n5_m3_readItem1_retrieval() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{3}, ReplacementStrategy.RR);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{}, 0);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output states for cache miss");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 possible outcome");

        // READ switches the job to the retrieval (pending) class and marks item 1 as being
        // fetched in block A, so that a concurrent request for it is a delayed hit.
        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5, 1, 0, 0, 0, 0});

        // Should switch to the retrieval class for item 1
        assertJobClass(result.outspace, 0, varStart, (int) model.retrievalPendingClasses.get(0, 0));

        assertTrue(result.outrate.elementSum() > 0, "Total rate should be positive");
    }

    @Test
    public void testRS_n5_m3_readItem1_retrieval_non_empty_system() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{3}, ReplacementStrategy.RR);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{2}, 0);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output states for cache miss");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 possible outcome");

        // READ switches the job to the retrieval (pending) class and marks item 1 as being
        // fetched in block A, so that a concurrent request for it is a delayed hit.
        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5, 1, 1, 0, 0, 0});

        // Should switch to the retrieval class for item 1
        assertJobClass(result.outspace, 0, varStart, (int) model.retrievalPendingClasses.get(0, 0));

        assertTrue(result.outrate.elementSum() > 0, "Total rate should be positive");
    }

    @Test
    public void testRS_n5_m3_readItem1_miss() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{3}, ReplacementStrategy.RR);
        int retrievalCompleteClass = (int) model.retrievalCompleteClasses.get(0, 0);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{1}, retrievalCompleteClass);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ,
                retrievalCompleteClass, false);

        assertNotNull(result.outspace, "Should return output states for cache miss");
        assertEquals(3, result.outspace.getNumRows(), "Should have 3 possible outcomes for m=[3]");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{1, 4, 5});
        assertCacheStateRow(result.outspace, 1, varStart, new int[]{3, 1, 5});
        assertCacheStateRow(result.outspace, 2, varStart, new int[]{3, 4, 1});

        // Verify equal probabilities
        double totalRate = result.outrate.elementSum();
        assertTrue(totalRate > 0, "Total rate should be positive");
        for (int i = 0; i < 3; i++) {
            assertEquals(1.0/3.0, result.outrate.get(i, 0) / totalRate, TOL,
                    "Each outcome should have probability 1/3");
        }

        // Verify switch to miss class
        assertJobClass(result.outspace, 0, varStart, (int) model.missClass.get(0, retrievalCompleteClass));
    }

    @Test
    public void testRS_n5_m3_readItem1_miss_first_retrieval() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{3}, ReplacementStrategy.RR);
        int retrievalCompleteClass = (int) model.retrievalCompleteClasses.get(0, 0);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{1, 2}, retrievalCompleteClass);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ,
                retrievalCompleteClass, false);

        assertNotNull(result.outspace, "Should return output states for cache miss");
        assertEquals(3, result.outspace.getNumRows(), "Should have 3 possible outcomes for m=[3]");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{1, 4, 5, 0, 1, 0, 0, 0});
        assertCacheStateRow(result.outspace, 1, varStart, new int[]{3, 1, 5, 0, 1, 0, 0, 0});
        assertCacheStateRow(result.outspace, 2, varStart, new int[]{3, 4, 1, 0, 1, 0, 0, 0});

        // Verify equal probabilities
        double totalRate = result.outrate.elementSum();
        assertTrue(totalRate > 0, "Total rate should be positive");
        for (int i = 0; i < 3; i++) {
            assertEquals(1.0/3.0, result.outrate.get(i, 0) / totalRate, TOL,
                    "Each outcome should have probability 1/3");
        }

        // Verify switch to miss class
        assertJobClass(result.outspace, 0, varStart, (int) model.missClass.get(0, retrievalCompleteClass));
    }

    @Test
    public void testRS_n5_m3_readItem1_delayed_hit() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{3}, ReplacementStrategy.RR);
        seedNodeSpace(model, 1);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{1}, 0);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output states for cache miss");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 possible outcome");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5, 1, 0, 0, 0, 0});

        int numActiveClasses = 0;
        for (int i = 0; i < varStart; i++) {
            numActiveClasses += (int) result.outspace.get(i);
        }

        assertEquals(0, numActiveClasses, "Should be no outgoing class for delayed hit");

        assertTrue(result.outrate.elementSum() > 0, "Total rate should be positive");
    }

    @Test
    public void testRS_n5_m3_readItem3_hit() {
        // Java: [3,4,5] -> [3,4,5], prob=1.0
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{3}, ReplacementStrategy.RR);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{}, 0);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 2); // item 3

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 outcome for hit in highest list");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5, 0, 0, 0, 0, 0});

        double accessProb = model.pAccess.get(2); // Access probability for item 3
        assertEquals(accessProb * GlobalConstants.Immediate, result.outrate.get(0, 0), TOL,
                "Outcome should have an instant rate scaled by access probability");

        // Verify switch to hit class
        assertJobClass(result.outspace, 0, varStart, (int) model.hitClass.get(0, 0));
    }

    @Test
    public void testRS_n5_m3_readItem5_hit() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{3}, ReplacementStrategy.RR);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{}, 0);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 4); // item 5

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 outcome for hit in highest list");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5, 0, 0, 0, 0, 0});

        double accessProb = model.pAccess.get(4); // Access probability for item 5
        assertEquals(accessProb * GlobalConstants.Immediate, result.outrate.get(0, 0), TOL,
                "Outcome should have an instant rate scaled by access probability");

        // Verify switch to hit class
        assertJobClass(result.outspace, 0, varStart, (int) model.hitClass.get(0, 0));
    }

    // ==================== n=5, m=[2,1] with retrieval system ====================

    @Test
    public void testRS_n5_m21_readItem1_retrieval() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{2, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{}, 0);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output states for cache miss");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 possible outcome");

        // READ marks item 1 as being fetched in block A.
        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5, 1, 0, 0, 0, 0});

        assertJobClass(result.outspace, 0, varStart, (int) model.retrievalPendingClasses.get(0, 0));

        assertTrue(result.outrate.elementSum() > 0, "Total rate should be positive");
    }

    @Test
    public void testRS_n5_m21_readItem1_retrieval_non_empty_system() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{2, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{2}, 0);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output states for cache miss");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 possible outcome");

        // READ marks item 1 as being fetched in block A.
        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5, 1, 1, 0, 0, 0});

        assertJobClass(result.outspace, 0, varStart, (int) model.retrievalPendingClasses.get(0, 0));

        assertTrue(result.outrate.elementSum() > 0, "Total rate should be positive");
    }

    @Test
    public void testRS_n5_m21_readItem1_delayed_hit() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{2, 1}, ReplacementStrategy.RR);
        seedNodeSpace(model, 1);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{1}, 0);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output states for delayed hit");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 possible outcome");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5, 1, 0, 0, 0, 0});

        int numActiveClasses = 0;
        for (int i = 0; i < varStart; i++) {
            numActiveClasses += (int) result.outspace.get(i);
        }
        assertEquals(0, numActiveClasses, "Should be no outgoing class for delayed hit");

        assertTrue(result.outrate.elementSum() > 0, "Total rate should be positive");
    }

    @Test
    public void testRS_n5_m21_readItem1_miss() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{2, 1}, ReplacementStrategy.RR);
        int retrievalCompleteClass = (int) model.retrievalCompleteClasses.get(0, 0);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{1}, retrievalCompleteClass);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ,
                retrievalCompleteClass, false);

        assertNotNull(result.outspace, "Should return output states for retrieval complete");
        assertEquals(3, result.outspace.getNumRows(), "Should have 3 possible outcomes for m=[2,1]");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{1, 4, 5, 0, 0, 0, 0, 0});
        assertCacheStateRow(result.outspace, 1, varStart, new int[]{3, 1, 5, 0, 0, 0, 0, 0});

        int nonZeroCount = countNonZeroRates(result.outrate);
        assertEquals(2, nonZeroCount, "Should have 2 outcomes with non-zero rate (list 0 has 2 slots)");

        double totalRate = result.outrate.elementSum();
        assertTrue(totalRate > 0, "Total rate should be positive");
        for (int i = 0; i < result.outspace.getNumRows(); i++) {
            if (result.outrate.get(i, 0) > 0) {
                assertEquals(1.0 / 2.0, result.outrate.get(i, 0) / totalRate, TOL,
                        "Each non-zero outcome should have probability 1/2");
            }
        }

        assertJobClass(result.outspace, 0, varStart, (int) model.missClass.get(0, retrievalCompleteClass));
    }

    @Test
    public void testRS_n5_m21_readItem1_miss_first_retrieval() {
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{2, 1}, ReplacementStrategy.RR);
        int retrievalCompleteClass = (int) model.retrievalCompleteClasses.get(0, 0);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{1, 2}, retrievalCompleteClass);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 0); // item 1

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ,
                retrievalCompleteClass, false);

        assertNotNull(result.outspace, "Should return output states for retrieval complete");
        assertEquals(3, result.outspace.getNumRows(), "Should have 3 possible outcomes for m=[2,1]");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{1, 4, 5, 0, 1, 0, 0, 0});
        assertCacheStateRow(result.outspace, 1, varStart, new int[]{3, 1, 5, 0, 1, 0, 0, 0});

        int nonZeroCount = countNonZeroRates(result.outrate);
        assertEquals(2, nonZeroCount, "Should have 2 outcomes with non-zero rate (list 0 has 2 slots)");

        double totalRate = result.outrate.elementSum();
        assertTrue(totalRate > 0, "Total rate should be positive");
        for (int i = 0; i < result.outspace.getNumRows(); i++) {
            if (result.outrate.get(i, 0) > 0) {
                assertEquals(1.0 / 2.0, result.outrate.get(i, 0) / totalRate, TOL,
                        "Each non-zero outcome should have probability 1/2");
            }
        }

        assertJobClass(result.outspace, 0, varStart, (int) model.missClass.get(0, retrievalCompleteClass));
    }

    @Test
    public void testRS_n5_m21_readItem3_hitList0() {
        // Item 3 is in list 0, promoted to list 1: [3,4,5] -> [5,4,3]
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{2, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{}, 0);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 2); // item 3

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertTrue(result.outspace.getNumRows() >= 1, "Should have at least 1 outcome");

        int nonZeroCount = countNonZeroRates(result.outrate);
        assertEquals(1, nonZeroCount, "Should have 1 outcome with non-zero rate");

        int varStart = model.numClasses;
        for (int i = 0; i < result.outspace.getNumRows(); i++) {
            if (result.outrate.get(i, 0) > 0) {
                assertCacheStateRow(result.outspace, i, varStart, new int[]{5, 4, 3, 0, 0, 0, 0, 0});
                assertJobClass(result.outspace, i, varStart, (int) model.hitClass.get(0, 0));
                double accessProb = model.pAccess.get(2); // item 3, 0-indexed
                assertEquals(accessProb * GlobalConstants.Immediate, result.outrate.get(i, 0), TOL,
                        "Rate should be access probability scaled by immediate rate");
            }
        }
    }

    @Test
    public void testRS_n5_m21_readItem5_hitList1() {
        // Item 5 is in list 1 (top list): [3,4,5] -> [3,4,5]
        TestCacheModel model = createCacheModelWithRetrievalSystem(5, new int[]{2, 1}, ReplacementStrategy.RR);
        Matrix inspace = createCacheStateWithRetrievalSystem(model, new int[]{3, 4, 5}, new int[]{}, 0);

        model.sn.varsparam.set(model.cacheNodeIndex, 0, 4); // item 5

        Ret.EventResult result = State.afterEvent(model.sn, model.cacheNodeIndex, inspace, EventType.READ, 0, false);

        assertNotNull(result.outspace, "Should return output state for cache hit");
        assertEquals(1, result.outspace.getNumRows(), "Should have 1 outcome for hit in highest list");

        int varStart = model.numClasses;
        assertCacheStateRow(result.outspace, 0, varStart, new int[]{3, 4, 5, 0, 0, 0, 0, 0});

        double accessProb = model.pAccess.get(4); // item 5, 0-indexed
        assertEquals(accessProb * GlobalConstants.Immediate, result.outrate.get(0, 0), TOL,
                "Outcome should have an instant rate scaled by access probability");

        assertJobClass(result.outspace, 0, varStart, (int) model.hitClass.get(0, 0));
    }

    // ==================== Helper methods ====================

    private static class TestCacheModel {
        Network network;
        NetworkStruct sn;
        int cacheNodeIndex;
        int numClasses;
        int numItems;
        int[] capacity;
        int totalCapacity;
        Matrix hitClass;
        Matrix missClass;
        int retrievalSystemCapacity;
        Matrix retrievalPendingClasses;
        Matrix retrievalCompleteClasses;
        Matrix pAccess;
    }

    private TestCacheModel createCacheModel(int nItems, int[] capacity, ReplacementStrategy strategy) {
        Network model = new Network("CacheAfterEventTest");

        Delay clientDelay = new Delay(model, "Client");

        Matrix itemLevelCap = new Matrix(1, capacity.length);
        int totalCap = 0;
        for (int i = 0; i < capacity.length; i++) {
            itemLevelCap.set(0, i, capacity[i]);
            totalCap += capacity[i];
        }

        Cache cacheNode = new Cache(model, "Cache", nItems, itemLevelCap, strategy);
        Delay cacheDelay = new Delay(model, "CacheDelay");

        ClosedClass clientClass = new ClosedClass(model, "ClientClass", 1, clientDelay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, clientDelay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, clientDelay, 0);

        clientDelay.setService(clientClass, new Immediate());
        cacheDelay.setService(hitClass, Exp.fitMean(0.2));
        cacheDelay.setService(missClass, Exp.fitMean(1.0));

        cacheNode.setRead(clientClass, new Zipf(1.0, nItems));
        cacheNode.setHitClass(clientClass, hitClass);
        cacheNode.setMissClass(clientClass, missClass);

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(clientClass, clientClass, clientDelay, cacheNode, 1.0);
        P.set(hitClass, hitClass, cacheNode, cacheDelay, 1.0);
        P.set(missClass, missClass, cacheNode, cacheDelay, 1.0);
        P.set(hitClass, clientClass, cacheDelay, clientDelay, 1.0);
        P.set(missClass, clientClass, cacheDelay, clientDelay, 1.0);

        model.link(P);

        NetworkStruct sn = model.getStruct();

        int cacheNodeIndex = -1;
        for (int i = 0; i < sn.nnodes; i++) {
            if (sn.nodes.get(i) == cacheNode) {
                cacheNodeIndex = i;
                break;
            }
        }

        TestCacheModel result = new TestCacheModel();
        result.network = model;
        result.sn = sn;
        result.cacheNodeIndex = cacheNodeIndex;
        result.numClasses = sn.nclasses;
        result.numItems = nItems;
        result.capacity = capacity;
        result.totalCapacity = totalCap;

        return result;
    }

    private TestCacheModel createCacheModelWithRetrievalSystem(int nItems, int[] capacity, ReplacementStrategy strategy) {
        Network model = new Network("CacheAfterEventTest");

        Delay clientDelay = new Delay(model, "Client");

        Matrix itemLevelCap = new Matrix(1, capacity.length);
        int totalCap = 0;
        for (int i = 0; i < capacity.length; i++) {
            itemLevelCap.set(0, i, capacity[i]);
            totalCap += capacity[i];
        }

        Cache cacheNode = new Cache(model, "Cache", nItems, itemLevelCap, strategy);
        Delay cacheDelay = new Delay(model, "CacheDelay");
        Queue queue = new Queue(model, "Queue");

        ClosedClass clientClass = new ClosedClass(model, "ClientClass", 1, clientDelay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, clientDelay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, clientDelay, 0);

        clientDelay.setService(clientClass, new Immediate());
        cacheDelay.setService(hitClass, Exp.fitMean(0.2));
        cacheDelay.setService(missClass, Exp.fitMean(1.0));

        Zipf pAccess = new Zipf(1.0, nItems);

        cacheNode.setRead(clientClass, pAccess);
        cacheNode.setHitClass(clientClass, hitClass);
        cacheNode.setMissClass(clientClass, missClass);

        // Single-queue retrieval system: read-class service (rate 1.0) and cache <-> queue
        // routing for every item, wired via the per-item methods.
        queue.setService(clientClass, new Exp(1.0));
        cacheNode.setRetrievalSystem(clientClass, missClass, new Queue[] {queue});
        for (int item = 0; item < nItems; item++) {
            cacheNode.setItemRoutingProb(clientClass, item, cacheNode, queue, 1.0); // cache -> queue
            cacheNode.setItemRoutingProb(clientClass, item, queue, cacheNode, 1.0);  // queue -> cache
        }

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(clientClass, clientClass, clientDelay, cacheNode, 1.0);
        P.set(hitClass, hitClass, cacheNode, cacheDelay, 1.0);
        P.set(missClass, missClass, cacheNode, cacheDelay, 1.0);
        P.set(hitClass, clientClass, cacheDelay, clientDelay, 1.0);
        P.set(missClass, clientClass, cacheDelay, clientDelay, 1.0);

        model.link(P);

        NetworkStruct sn = model.getStruct();

        int cacheNodeIndex = -1;
        for (int i = 0; i < sn.nnodes; i++) {
            if (sn.nodes.get(i) == cacheNode) {
                cacheNodeIndex = i;
                break;
            }
        }

        TestCacheModel result = new TestCacheModel();
        result.network = model;
        result.sn = sn;
        result.cacheNodeIndex = cacheNodeIndex;
        result.numClasses = sn.nclasses;
        result.numItems = nItems;
        result.capacity = capacity;
        result.totalCapacity = totalCap;
        result.hitClass = cacheNode.getHitClass();
        result.missClass = cacheNode.getMissClass();
        result.retrievalSystemCapacity = cacheNode.getRetrievalSystemCapacity();
        result.retrievalPendingClasses = cacheNode.getRetrievalClasses();
        result.retrievalCompleteClasses = cacheNode.getRetrievalClasses();
        result.pAccess = pAccess.evalPMF();

        return result;
    }

    private Matrix createCacheState(TestCacheModel model, int[] cachedItems) {
        assertEquals(model.totalCapacity, cachedItems.length,
            "cachedItems length must match total cache capacity");

        int numCols = model.numClasses + model.totalCapacity;
        Matrix state = new Matrix(1, numCols);

        state.set(0, 0, 1);
        for (int r = 1; r < model.numClasses; r++) {
            state.set(0, r, 0);
        }

        for (int i = 0; i < model.totalCapacity; i++) {
            state.set(0, model.numClasses + i, cachedItems[i]);
        }

        return state;
    }

    /**
     * Store the enumerated local-variable space of the cache node, as the state-space
     * generator would. afterEventCache reads the delayed-hit truncation level off it,
     * so a merge onto an in-flight fetch is only representable once it is present.
     */
    private void seedNodeSpace(TestCacheModel model, int maxPending) {
        model.sn.space.put(model.sn.stateful.get((int) model.sn.nodeToStateful.get(model.cacheNodeIndex)),
                State.spaceLocalVarsPublic(model.sn, model.cacheNodeIndex, maxPending));
    }

    private Matrix createCacheStateWithRetrievalSystem(TestCacheModel model, int[] cachedItems, int[] itemsInRetrieval,
                                                       int jobClassIndex) {
        assertTrue(itemsInRetrieval.length <= model.retrievalSystemCapacity,
                "itemsInRetrieval length must not exceed retrievalSystemCapacity");

        Matrix cacheStateNoRetrieval = createCacheState(model, cachedItems);
        int noRetrievalCols = cacheStateNoRetrieval.getNumCols();

        // Local-variable layout of a retrieval cache: cache contents, then block A (a
        // per-item occupancy bitmap), then block B (one count per retrieval class)
        int widthB = State.cacheRetrievalClassMap(model.sn, model.cacheNodeIndex)[0].length;
        int numCols = noRetrievalCols + model.numItems + widthB;
        Matrix state = new Matrix(1, numCols);

        for (int i = 0; i < noRetrievalCols; i++) {
            state.set(0, i, cacheStateNoRetrieval.get(i));
        }
        state.set(0, 0, 0);
        state.set(0, jobClassIndex, 1);

        for (int i = 0; i < itemsInRetrieval.length; i++) {
            // itemsInRetrieval holds 1-indexed item identifiers; set the corresponding occupancy bit
            state.set(0, noRetrievalCols + itemsInRetrieval[i] - 1, 1);
        }

        return state;
    }

    private int countNonZeroRates(Matrix outrate) {
        int count = 0;
        for (int i = 0; i < outrate.getNumRows(); i++) {
            if (outrate.get(i, 0) > 0) {
                count++;
            }
        }
        return count;
    }

    private void assertCacheStateRow(Matrix outspace, int row, int varStart, int[] expected) {
        StringBuilder actual = new StringBuilder("[");
        StringBuilder exp = new StringBuilder("[");

        for (int i = 0; i < expected.length; i++) {
            double actualVal = outspace.get(row, varStart + i);
            actual.append((int) actualVal);
            exp.append(expected[i]);
            if (i < expected.length - 1) {
                actual.append(", ");
                exp.append(", ");
            }
        }
        actual.append("]");
        exp.append("]");

        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], (int) outspace.get(row, varStart + i),
                String.format("Row %d, position %d mismatch. Expected %s, got %s",
                    row, i, exp.toString(), actual.toString()));
        }
    }

    private void assertJobClass(Matrix outspace, int row, int varStart, int expectedClassIndex) {
        StringBuilder actual = new StringBuilder("[");
        StringBuilder exp = new StringBuilder("[");
        int[] expected = new int[varStart];
        expected[expectedClassIndex] = 1;

        for (int i = 0; i < varStart; i++) {
            double actualVal = outspace.get(row, i);
            actual.append((int) actualVal);
            exp.append(expected[i]);
            if (i < expected.length - 1) {
                actual.append(", ");
                exp.append(", ");
            }
        }

        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], (int) outspace.get(row, i),
                    String.format("Row %d, position %d mismatch. Expected %s, got %s",
                            row, i, exp.toString(), actual.toString()));
        }
    }
}
