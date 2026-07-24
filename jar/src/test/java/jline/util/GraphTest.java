package jline.util;

import jline.util.graph.DirectedGraph;
import jline.util.graph.DirectedGraph.SCCResult;
import jline.util.graph.UndirectedGraph;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.Set;

import static org.ejml.UtilEjml.assertTrue;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class GraphTest {

    @Test
    public void testWeaklyConnectedComponents1() {
        // Initialize the Matrix using a string representation for an undirected graph
        Matrix A = new Matrix("[0.2,0.3,0,0.5,0,0; 0,0.5,0.5,0,0,0; 0,0.5,0.5,0,0,0; 0,0,0,0,1,0; 0,0,0,0,0,1; 0,0,0,1,0,0]");

        // Optional: Set of columns to ignore (none in this case)
        Set<Integer> colsToIgnore = new HashSet<>();

        // Create the UndirectedGraph instance using the adjacency matrix
        UndirectedGraph undirectedGraph = new UndirectedGraph(A, colsToIgnore);

        // Compute the weakly connected components (WCC)
        undirectedGraph.computeWeaklyConnectedComponents();

        // Get the weakly connected components (WCC)
        Set<Set<Integer>> wcc = undirectedGraph.getWCC();

        // Expected WCC output: nodes grouped into connected components
        // For this adjacency matrix, the expected WCCs could be:
        Set<Set<Integer>> expectedWCC = new HashSet<>();

        Set<Integer> component1 = new HashSet<>();
        component1.add(0);
        component1.add(1);
        component1.add(2);
        component1.add(3);
        component1.add(4);
        component1.add(5);
        expectedWCC.add(component1);  // First component

        // Validate that the number of components matches
        assertEquals(expectedWCC.size(), wcc.size(), "Number of weakly connected components does not match.");

        // Validate that the WCC sets are as expected
        for (Set<Integer> expectedComponent : expectedWCC) {
            assertTrue(wcc.contains(expectedComponent), "Expected component " + expectedComponent + " is missing.");
        }
    }

    @Test
    public void testStronglyConnectedComponents() {
        Matrix adjacencyMatrix = new Matrix("[0.2,0.3,0,0.5,0,0; 0,0.5,0.5,0,0,0; 0,0.5,0.5,0,0,0; 0,0,0,0,1,0; 0,0,0,0,0,1; 0,0,0,1,0,0]");

        Set<Integer> colsToIgnore = new HashSet<>();

        DirectedGraph directedGraph = new DirectedGraph(adjacencyMatrix, colsToIgnore);

        SCCResult result = directedGraph.stronglyconncomp();

        int[] expectedI = {3, 2, 2, 1, 1, 1}; // SCC assignments for each node

        boolean[] expectedRecurrent = {true, true, false};

        assertArrayEquals(expectedI, result.I, "SCC group assignments do not match the expected values");
        assertArrayEquals(expectedRecurrent, result.recurrent, "Recurrent status of SCCs does not match the expected values");
    }
}
