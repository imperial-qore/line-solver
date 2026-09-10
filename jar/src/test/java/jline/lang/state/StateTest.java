package jline.lang.state;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for State utility methods (state space generation).
 * SPN event tests (afterEvent, afterGlobalEvent) are in PetriStateTest.
 */
public class StateTest {

    // ==================== STATE SPACE GENERATION TESTS ====================

    @Test
    public void testSpaceClosedSingle_Basic() {
        int M = 2;
        int N = 3;

        Matrix space = State.spaceClosedSinglePublic(M, N);

        assertNotNull(space, "State space should not be null");
        // C(N+M-1, M-1) = C(4,1) = 4
        assertEquals(4, space.getNumRows(), "Should have 4 possible states");
        assertEquals(M, space.getNumCols(), "Each state should have M columns");

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
        int M = 3;
        int N = 0;

        Matrix space = State.spaceClosedSinglePublic(M, N);

        assertNotNull(space);
        assertEquals(1, space.getNumRows(), "Zero jobs should have one state (all zeros)");

        for (int j = 0; j < space.getNumCols(); j++) {
            assertEquals(0, space.get(0, j), 0.0);
        }
    }

    @Test
    public void testSpaceClosedMulti_TwoClasses() {
        int M = 2;
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 2);
        N.set(0, 1, 1);

        Matrix space = State.spaceClosedMulti(M, N);

        assertNotNull(space, "State space should not be null");
        assertTrue(space.getNumRows() > 0, "Should have positive number of states");
        assertEquals(M * 2, space.getNumCols(), "Should have M*R columns");
    }

    @Test
    public void testSpaceClosedMulti_SingleClass() {
        int M = 3;
        Matrix N = new Matrix(1, 1);
        N.set(0, 0, 2);

        Matrix space = State.spaceClosedMulti(M, N);

        assertNotNull(space);
        assertEquals(M, space.getNumCols(), "Single class should have M columns");
        // C(N+M-1, M-1) = C(4, 2) = 6
        assertEquals(6, space.getNumRows(), "Should have 6 states for M=3, N=2");
    }

    @Test
    public void testSpaceClosedMultiCS_WithChains() {
        int M = 2;
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 2);
        N.set(0, 1, 1);

        Matrix chains = new Matrix(1, 2);
        chains.set(0, 0, 1);
        chains.set(0, 1, 1);

        Matrix space = State.spaceClosedMultiCS(M, N, chains);

        assertNotNull(space, "State space with chains should not be null");
        assertTrue(space.getNumRows() > 0, "Should have positive number of states");
    }

    @Test
    public void testSpaceClosedMultiCS_SeparateChains() {
        int M = 2;
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 1);
        N.set(0, 1, 1);

        Matrix chains = new Matrix(2, 2);
        chains.set(0, 0, 1); chains.set(0, 1, 0);
        chains.set(1, 0, 0); chains.set(1, 1, 1);

        Matrix space = State.spaceClosedMultiCS(M, N, chains);

        assertNotNull(space);
        Matrix spaceNoCS = State.spaceClosedMulti(M, N);
        assertEquals(spaceNoCS.getNumRows(), space.getNumRows(),
            "Separate chains should produce same number of states");
    }

    @Test
    public void testSpaceClosedSingle_LargerNetwork() {
        int M = 4;
        int N = 2;

        Matrix space = State.spaceClosedSinglePublic(M, N);

        // C(N+M-1, M-1) = C(5, 3) = 10
        assertEquals(10, space.getNumRows(), "Should have 10 states for M=4, N=2");
        assertEquals(M, space.getNumCols());

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
        int M = 3;
        int N = 3;

        Matrix space = State.spaceClosedSinglePublic(M, N);

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

    @Test
    public void testSpaceCache_NoRetrievalSystem() {
        int n = 1;
        Matrix m = Matrix.singleton(1);
        int retrievalSystemCapacity = 0;

        Matrix space = State.spaceCachePublic(n, m, retrievalSystemCapacity);

        assertNotNull(space, "State space should not be null");
        assertEquals(1, space.getNumRows(), "Should have 1 state");
        assertEquals(1, space.getNumCols(), "Should have 1 column");

        assertEquals(1, space.get(0, 0), 0.0, "Only possible cache configuration is 1");
    }

    @Test
    public void testSpaceCache_MultipleItemsNoRetrievalSystem() {
        int n = 4;
        Matrix m = Matrix.singleton(2);
        int totalCacheCapacity = 0;

        Matrix space = State.spaceCachePublic(n, m, totalCacheCapacity);

        assertNotNull(space, "State space should not be null");
        assertEquals(12, space.getNumRows(), "Should have 12 states");
        assertEquals(2, space.getNumCols(), "Should have 2 columns");

        double[][] expected = {
                {1.0, 2.0},
                {1.0, 3.0},
                {1.0, 4.0},
                {2.0, 1.0},
                {2.0, 3.0},
                {2.0, 4.0},
                {3.0, 1.0},
                {3.0, 2.0},
                {3.0, 4.0},
                {4.0, 1.0},
                {4.0, 2.0},
                {4.0, 3.0},
        };
        Matrix expectedStateSpace = new Matrix(expected);

        for (int rowIdx = 0; rowIdx < expectedStateSpace.getNumRows(); rowIdx++) {
            Matrix row = expectedStateSpace.getRow(rowIdx);

            List<Integer> matchingRows = Matrix.findRows(space, row);
            assertEquals(1, matchingRows.size(), "Cache state space missing state: " + row +
                    " expected state space: \n" + expectedStateSpace + "\nActual space:\n" + space);
        }
    }

    @Test
    public void testSpaceCache_WithIndividualRetrievalSystem() {
        int n = 2;
        Matrix m = Matrix.singleton(1);
        int retrievalSystemCapacity = 1;

        Matrix space = State.spaceCachePublic(n, m, retrievalSystemCapacity);

        assertNotNull(space, "State space should not be null");
        assertEquals(4, space.getNumRows(), "Should have 4 states");
        assertEquals(3, space.getNumCols(), "Should have 3 columns");

        // Layout: [cache slot | retrieval bitmap item1 | retrieval bitmap item2]
        double[][] expected = {
                {1.0, 0.0, 0.0},
                {1.0, 0.0, 1.0},
                {2.0, 0.0, 0.0},
                {2.0, 1.0, 0.0}
        };
        Matrix expectedStateSpace = new Matrix(expected);

        for (int rowIdx = 0; rowIdx < expectedStateSpace.getNumRows(); rowIdx++) {
            Matrix row = expectedStateSpace.getRow(rowIdx);

            List<Integer> matchingRows = Matrix.findRows(space, row);
            assertEquals(1, matchingRows.size(), "Cache state space missing state: " + row +
                    " expected state space: \n" + expectedStateSpace + "\nActual space:\n" + space);
        }
    }

    @Test
    public void testSpaceCache_WithMultipleItemRetrievalSystem() {
        int n = 3;
        Matrix m = Matrix.singleton(1);
        int retrievalSystemCapacity = 2;

        Matrix space = State.spaceCachePublic(n, m, retrievalSystemCapacity);

        assertNotNull(space, "State space should not be null");
        assertEquals(12, space.getNumRows(), "Should have 12 states");
        assertEquals(4, space.getNumCols(), "Should have 4 columns");

        // Layout: [cache slot | retrieval bitmap item1 | item2 | item3]
        double[][] expected = {
                {1.0, 0.0, 0.0, 0.0},
                {1.0, 0.0, 1.0, 0.0},
                {1.0, 0.0, 0.0, 1.0},
                {1.0, 0.0, 1.0, 1.0},
                {2.0, 0.0, 0.0, 0.0},
                {2.0, 1.0, 0.0, 0.0},
                {2.0, 0.0, 0.0, 1.0},
                {2.0, 1.0, 0.0, 1.0},
                {3.0, 0.0, 0.0, 0.0},
                {3.0, 1.0, 0.0, 0.0},
                {3.0, 0.0, 1.0, 0.0},
                {3.0, 1.0, 1.0, 0.0},
        };
        Matrix expectedStateSpace = new Matrix(expected);

        for (int rowIdx = 0; rowIdx < expectedStateSpace.getNumRows(); rowIdx++) {
            Matrix row = expectedStateSpace.getRow(rowIdx);

            List<Integer> matchingRows = Matrix.findRows(space, row);
            assertEquals(1, matchingRows.size(), "Cache state space missing state: " + row);
        }
    }

    @Test
    public void testSpaceCache_WithMultipleItemCacheAndRetrievalSystem() {
        int n = 4;
        Matrix m = Matrix.singleton(2);
        int retrievalSystemCapacity = 2;

        Matrix space = State.spaceCachePublic(n, m, retrievalSystemCapacity);

        assertNotNull(space, "State space should not be null");
        assertEquals(48, space.getNumRows(), "Should have 48 states");
        assertEquals(6, space.getNumCols(), "Should have 6 columns");

        // Layout: [cache slot 0 | cache slot 1 | retrieval bitmap item1 | item2 | item3 | item4]
        double[][] expected = {
                {1.0, 2.0, 0.0, 0.0, 0.0, 0.0},
                {1.0, 2.0, 0.0, 0.0, 1.0, 0.0},
                {1.0, 2.0, 0.0, 0.0, 0.0, 1.0},
                {1.0, 2.0, 0.0, 0.0, 1.0, 1.0},
                {2.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                {2.0, 1.0, 0.0, 0.0, 1.0, 0.0},
                {2.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                {2.0, 1.0, 0.0, 0.0, 1.0, 1.0},
                {1.0, 3.0, 0.0, 0.0, 0.0, 0.0},
                {1.0, 3.0, 0.0, 1.0, 0.0, 0.0},
                {1.0, 3.0, 0.0, 0.0, 0.0, 1.0},
                {1.0, 3.0, 0.0, 1.0, 0.0, 1.0},
                {3.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                {3.0, 1.0, 0.0, 1.0, 0.0, 0.0},
                {3.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                {3.0, 1.0, 0.0, 1.0, 0.0, 1.0},
                {1.0, 4.0, 0.0, 0.0, 0.0, 0.0},
                {1.0, 4.0, 0.0, 1.0, 0.0, 0.0},
                {1.0, 4.0, 0.0, 0.0, 1.0, 0.0},
                {1.0, 4.0, 0.0, 1.0, 1.0, 0.0},
                {4.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                {4.0, 1.0, 0.0, 1.0, 0.0, 0.0},
                {4.0, 1.0, 0.0, 0.0, 1.0, 0.0},
                {4.0, 1.0, 0.0, 1.0, 1.0, 0.0},
                {2.0, 3.0, 0.0, 0.0, 0.0, 0.0},
                {2.0, 3.0, 1.0, 0.0, 0.0, 0.0},
                {2.0, 3.0, 0.0, 0.0, 0.0, 1.0},
                {2.0, 3.0, 1.0, 0.0, 0.0, 1.0},
                {3.0, 2.0, 0.0, 0.0, 0.0, 0.0},
                {3.0, 2.0, 1.0, 0.0, 0.0, 0.0},
                {3.0, 2.0, 0.0, 0.0, 0.0, 1.0},
                {3.0, 2.0, 1.0, 0.0, 0.0, 1.0},
                {2.0, 4.0, 0.0, 0.0, 0.0, 0.0},
                {2.0, 4.0, 1.0, 0.0, 0.0, 0.0},
                {2.0, 4.0, 0.0, 0.0, 1.0, 0.0},
                {2.0, 4.0, 1.0, 0.0, 1.0, 0.0},
                {4.0, 2.0, 0.0, 0.0, 0.0, 0.0},
                {4.0, 2.0, 1.0, 0.0, 0.0, 0.0},
                {4.0, 2.0, 0.0, 0.0, 1.0, 0.0},
                {4.0, 2.0, 1.0, 0.0, 1.0, 0.0},
                {3.0, 4.0, 0.0, 0.0, 0.0, 0.0},
                {3.0, 4.0, 1.0, 0.0, 0.0, 0.0},
                {3.0, 4.0, 0.0, 1.0, 0.0, 0.0},
                {3.0, 4.0, 1.0, 1.0, 0.0, 0.0},
                {4.0, 3.0, 0.0, 0.0, 0.0, 0.0},
                {4.0, 3.0, 1.0, 0.0, 0.0, 0.0},
                {4.0, 3.0, 0.0, 1.0, 0.0, 0.0},
                {4.0, 3.0, 1.0, 1.0, 0.0, 0.0},
        };

        Matrix expectedStateSpace = new Matrix(expected);

        for (int rowIdx = 0; rowIdx < expectedStateSpace.getNumRows(); rowIdx++) {
            Matrix row = expectedStateSpace.getRow(rowIdx);

            List<Integer> matchingRows = Matrix.findRows(space, row);
            assertEquals(1, matchingRows.size(), "Cache state space missing state: " + row);
        }
    }

    @Test
    public void testSpaceCache_WithMatrixCapacity() {
        int n = 4;
        Matrix m = new Matrix(new double[] {1, 1});
        int retrievalSystemCapacity = 2;

        Matrix space = State.spaceCachePublic(n, m, retrievalSystemCapacity);

        assertNotNull(space, "State space should not be null");
        assertEquals(48, space.getNumRows(), "Should have 48 states");
        assertEquals(6, space.getNumCols(), "Should have 6 columns");

        // Layout: [cache slot 0 | cache slot 1 | retrieval bitmap item1 | item2 | item3 | item4]
        double[][] expected = {
                {1.0, 2.0, 0.0, 0.0, 0.0, 0.0},
                {1.0, 2.0, 0.0, 0.0, 1.0, 0.0},
                {1.0, 2.0, 0.0, 0.0, 0.0, 1.0},
                {1.0, 2.0, 0.0, 0.0, 1.0, 1.0},
                {2.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                {2.0, 1.0, 0.0, 0.0, 1.0, 0.0},
                {2.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                {2.0, 1.0, 0.0, 0.0, 1.0, 1.0},
                {1.0, 3.0, 0.0, 0.0, 0.0, 0.0},
                {1.0, 3.0, 0.0, 1.0, 0.0, 0.0},
                {1.0, 3.0, 0.0, 0.0, 0.0, 1.0},
                {1.0, 3.0, 0.0, 1.0, 0.0, 1.0},
                {3.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                {3.0, 1.0, 0.0, 1.0, 0.0, 0.0},
                {3.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                {3.0, 1.0, 0.0, 1.0, 0.0, 1.0},
                {1.0, 4.0, 0.0, 0.0, 0.0, 0.0},
                {1.0, 4.0, 0.0, 1.0, 0.0, 0.0},
                {1.0, 4.0, 0.0, 0.0, 1.0, 0.0},
                {1.0, 4.0, 0.0, 1.0, 1.0, 0.0},
                {4.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                {4.0, 1.0, 0.0, 1.0, 0.0, 0.0},
                {4.0, 1.0, 0.0, 0.0, 1.0, 0.0},
                {4.0, 1.0, 0.0, 1.0, 1.0, 0.0},
                {2.0, 3.0, 0.0, 0.0, 0.0, 0.0},
                {2.0, 3.0, 1.0, 0.0, 0.0, 0.0},
                {2.0, 3.0, 0.0, 0.0, 0.0, 1.0},
                {2.0, 3.0, 1.0, 0.0, 0.0, 1.0},
                {3.0, 2.0, 0.0, 0.0, 0.0, 0.0},
                {3.0, 2.0, 1.0, 0.0, 0.0, 0.0},
                {3.0, 2.0, 0.0, 0.0, 0.0, 1.0},
                {3.0, 2.0, 1.0, 0.0, 0.0, 1.0},
                {2.0, 4.0, 0.0, 0.0, 0.0, 0.0},
                {2.0, 4.0, 1.0, 0.0, 0.0, 0.0},
                {2.0, 4.0, 0.0, 0.0, 1.0, 0.0},
                {2.0, 4.0, 1.0, 0.0, 1.0, 0.0},
                {4.0, 2.0, 0.0, 0.0, 0.0, 0.0},
                {4.0, 2.0, 1.0, 0.0, 0.0, 0.0},
                {4.0, 2.0, 0.0, 0.0, 1.0, 0.0},
                {4.0, 2.0, 1.0, 0.0, 1.0, 0.0},
                {3.0, 4.0, 0.0, 0.0, 0.0, 0.0},
                {3.0, 4.0, 1.0, 0.0, 0.0, 0.0},
                {3.0, 4.0, 0.0, 1.0, 0.0, 0.0},
                {3.0, 4.0, 1.0, 1.0, 0.0, 0.0},
                {4.0, 3.0, 0.0, 0.0, 0.0, 0.0},
                {4.0, 3.0, 1.0, 0.0, 0.0, 0.0},
                {4.0, 3.0, 0.0, 1.0, 0.0, 0.0},
                {4.0, 3.0, 1.0, 1.0, 0.0, 0.0},
        };

        Matrix expectedStateSpace = new Matrix(expected);

        for (int rowIdx = 0; rowIdx < expectedStateSpace.getNumRows(); rowIdx++) {
            Matrix row = expectedStateSpace.getRow(rowIdx);

            List<Integer> matchingRows = Matrix.findRows(space, row);
            assertEquals(1, matchingRows.size(), "Cache state space missing state: " + row);
        }
    }
}
