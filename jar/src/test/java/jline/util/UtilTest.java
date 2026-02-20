package jline.util;

import jline.util.graph.DirectedGraph;
import jline.util.graph.DirectedGraph.SCCResult;
import jline.util.graph.UndirectedGraph;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.AfterEach;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.ejml.UtilEjml.assertTrue;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Merged test class for utilities in the jline.util package.
 * Contains tests for Lattice, Graph, Maths, and MatFile utilities.
 * These tests were determined by using release 2.0.27 of the LINE
 * software as a baseline. These tests use as inputs and outputs the values observed through the MATLAB
 * workspace whenever the corresponding functions were called/executed.
 */
public class UtilTest {

    // ==================== LATTICE TESTS ====================

    @Test
    void testHashpop1() {
        Matrix n = new Matrix(1, 1);
        Matrix N = new Matrix(1, 1);
        int R = 1;
        Matrix prods = new Matrix(1, 1);
        N.set(0, 0, 16);
        prods.set(0, 0, 1);
        assertEquals(0, PopulationLattice.hashpop(n, N, R, prods));
    }

    @Test
    void testHashpop2() {
        Matrix n = new Matrix(1, 1);
        Matrix N = new Matrix(1, 1);
        int R = 1;
        Matrix prods = new Matrix(1, 1);
        n.set(0, 0, 1);
        N.set(0, 0, 16);
        prods.set(0, 0, 1);
        assertEquals(1, PopulationLattice.hashpop(n, N, R, prods));
    }

    @Test
    void testHashpop3() {
        Matrix n = new Matrix(1, 1);
        Matrix N = new Matrix(1, 1);
        int R = 1;
        Matrix prods = new Matrix(1, 1);
        n.set(0, 0, 2);
        N.set(0, 0, 16);
        prods.set(0, 0, 1);
        assertEquals(2, PopulationLattice.hashpop(n, N, R, prods));
    }

    @Test
    void testHashpop4() {
        Matrix n = new Matrix(1, 1);
        Matrix N = new Matrix(1, 1);
        int R = 1;
        Matrix prods = new Matrix(1, 1);
        N.set(0, 0, 3);
        prods.set(0, 0, 1);
        assertEquals(0, PopulationLattice.hashpop(n, N, R, prods));
    }

    // ==================== GRAPH TESTS ====================

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

    // ==================== MATHS TESTS ====================

    @Test
    void maxpos_shouldReturnIndexOfMaximumValue() {
        Matrix v = new Matrix(1, 5);
        v.set(0, 0, 1.0);
        v.set(0, 1, 5.0);
        v.set(0, 2, 3.0);
        v.set(0, 3, 7.0);
        v.set(0, 4, 2.0);
        
        int maxIndex = Maths.maxpos(v);
        assertEquals(3, maxIndex);
    }

    @Test
    void maxpos_shouldReturnFirstIndexWhenMultipleMaxValues() {
        Matrix v = new Matrix(1, 4);
        v.set(0, 0, 3.0);
        v.set(0, 1, 7.0);
        v.set(0, 2, 7.0);
        v.set(0, 3, 5.0);
        
        int maxIndex = Maths.maxpos(v);
        assertEquals(1, maxIndex);
    }

    @Test
    void maxpos_shouldHandleSingleElement() {
        Matrix v = new Matrix(1, 1);
        v.set(0, 0, 42.0);
        
        int maxIndex = Maths.maxpos(v);
        assertEquals(0, maxIndex);
    }

    @Test
    void maxpos_shouldHandleNegativeValues() {
        Matrix v = new Matrix(1, 3);
        v.set(0, 0, -5.0);
        v.set(0, 1, -2.0);
        v.set(0, 2, -10.0);
        
        int maxIndex = Maths.maxpos(v);
        assertEquals(1, maxIndex);
    }

    @Test
    void maxpos_shouldHandleColumnVector() {
        Matrix v = new Matrix(4, 1);
        v.set(0, 0, 1.0);
        v.set(1, 0, 8.0);
        v.set(2, 0, 3.0);
        v.set(3, 0, 5.0);
        
        int maxIndex = Maths.maxpos(v);
        assertEquals(1, maxIndex);
    }

    @Test
    void maxpos_shouldThrowExceptionForNullMatrix() {
        assertThrows(IllegalArgumentException.class, () -> Maths.maxpos(null));
    }

    @Test
    void maxpos_shouldThrowExceptionForEmptyMatrix() {
        Matrix v = new Matrix(0, 0);
        assertThrows(IllegalArgumentException.class, () -> Maths.maxpos(v));
    }

    @Test
    void maxpos_withN_shouldReturnTopNIndices() {
        Matrix v = new Matrix(1, 5);
        v.set(0, 0, 1.0);
        v.set(0, 1, 5.0);
        v.set(0, 2, 3.0);
        v.set(0, 3, 7.0);
        v.set(0, 4, 2.0);
        
        int[] topIndices = Maths.maxpos(v, 3);
        assertEquals(3, topIndices.length);
        assertEquals(3, topIndices[0]); // 7.0
        assertEquals(1, topIndices[1]); // 5.0
        assertEquals(2, topIndices[2]); // 3.0
    }

    @Test
    void maxpos_withN_shouldReturnAllIndicesWhenNEqualsLength() {
        Matrix v = new Matrix(1, 3);
        v.set(0, 0, 1.0);
        v.set(0, 1, 3.0);
        v.set(0, 2, 2.0);
        
        int[] topIndices = Maths.maxpos(v, 3);
        assertEquals(3, topIndices.length);
        assertEquals(1, topIndices[0]); // 3.0
        assertEquals(2, topIndices[1]); // 2.0
        assertEquals(0, topIndices[2]); // 1.0
    }

    @Test
    void maxpos_withN_shouldThrowExceptionWhenNIsZero() {
        Matrix v = new Matrix(1, 3);
        v.set(0, 0, 1.0);
        v.set(0, 1, 2.0);
        v.set(0, 2, 3.0);
        
        assertThrows(IllegalArgumentException.class, () -> Maths.maxpos(v, 0));
    }

    @Test
    void maxpos_withN_shouldThrowExceptionWhenNIsNegative() {
        Matrix v = new Matrix(1, 3);
        v.set(0, 0, 1.0);
        v.set(0, 1, 2.0);
        v.set(0, 2, 3.0);
        
        assertThrows(IllegalArgumentException.class, () -> Maths.maxpos(v, -1));
    }

    @Test
    void maxpos_withN_shouldThrowExceptionWhenNIsGreaterThanLength() {
        Matrix v = new Matrix(1, 3);
        v.set(0, 0, 1.0);
        v.set(0, 1, 2.0);
        v.set(0, 2, 3.0);
        
        assertThrows(IllegalArgumentException.class, () -> Maths.maxpos(v, 4));
    }

    @Test
    void probchoose_shouldReturnIndexBasedOnProbability() {
        // Set deterministic seed for reproducible test
        Maths.setMatlabRandomSeed(12345);
        
        List<Double> probs = Arrays.asList(0.1, 0.3, 0.4, 0.2);
        
        // Test multiple times to ensure it's working correctly
        int[] counts = new int[4];
        int numTests = 1000;
        
        for (int i = 0; i < numTests; i++) {
            int choice = Maths.probchoose(probs);
            assertTrue(choice >= 0 && choice < 4);
            counts[choice]++;
        }
        
        // Check that we got some selections for each category
        // (with 1000 tests, even low probability should get some hits)
        for (int i = 0; i < 4; i++) {
            assertTrue(counts[i] > 0, "Category " + i + " should have been selected at least once");
        }
    }

    @Test
    void probchoose_shouldReturnZeroForSingleElement() {
        List<Double> probs = Arrays.asList(1.0);
        
        int choice = Maths.probchoose(probs);
        assertEquals(0, choice);
    }

    @Test
    void probchoose_shouldReturnLastElementWhenAllProbabilitiesAreZero() {
        List<Double> probs = Arrays.asList(0.0, 0.0, 0.0, 0.0);
        
        int choice = Maths.probchoose(probs);
        assertEquals(3, choice);
    }

    @Test
    void probchoose_shouldHandleExtremelySmallProbabilities() {
        List<Double> probs = Arrays.asList(1e-10, 1e-10, 1.0 - 2e-10);
        
        int choice = Maths.probchoose(probs);
        assertTrue(choice >= 0 && choice < 3);
    }

    @Test
    void probchoose_shouldThrowExceptionForNullList() {
        assertThrows(IllegalArgumentException.class, () -> Maths.probchoose((List<Double>) null));
    }

    @Test
    void probchoose_shouldThrowExceptionForEmptyList() {
        List<Double> probs = new ArrayList<Double>();
        assertThrows(IllegalArgumentException.class, () -> Maths.probchoose(probs));
    }

    @Test
    void probchoose_matrix_shouldReturnIndexBasedOnProbability() {
        // Set deterministic seed for reproducible test
        Maths.setMatlabRandomSeed(12345);
        
        Matrix probs = new Matrix(1, 4);
        probs.set(0, 0, 0.1);
        probs.set(0, 1, 0.3);
        probs.set(0, 2, 0.4);
        probs.set(0, 3, 0.2);
        
        // Test multiple times to ensure it's working correctly
        int[] counts = new int[4];
        int numTests = 1000;
        
        for (int i = 0; i < numTests; i++) {
            int choice = Maths.probchoose(probs);
            assertTrue(choice >= 0 && choice < 4);
            counts[choice]++;
        }
        
        // Check that we got some selections for each category
        for (int i = 0; i < 4; i++) {
            assertTrue(counts[i] > 0, "Category " + i + " should have been selected at least once");
        }
    }

    @Test
    void probchoose_matrix_shouldHandleColumnVector() {
        Matrix probs = new Matrix(3, 1);
        probs.set(0, 0, 0.3);
        probs.set(1, 0, 0.4);
        probs.set(2, 0, 0.3);
        
        int choice = Maths.probchoose(probs);
        assertTrue(choice >= 0 && choice < 3);
    }

    @Test
    void probchoose_matrix_shouldThrowExceptionForNullMatrix() {
        assertThrows(IllegalArgumentException.class, () -> Maths.probchoose((Matrix) null));
    }

    @Test
    void probchoose_matrix_shouldThrowExceptionForEmptyMatrix() {
        Matrix probs = new Matrix(0, 0);
        assertThrows(IllegalArgumentException.class, () -> Maths.probchoose(probs));
    }

    @Test
    void probchoose_shouldHandleEdgeCaseWithProbabilitiesThatDontSumToOne() {
        // This tests the robustness of the implementation
        List<Double> probs = Arrays.asList(0.1, 0.2, 0.3); // sum = 0.6
        
        int choice = Maths.probchoose(probs);
        assertTrue(choice >= 0 && choice < 3);
    }

    // ==================== MATFILE TESTS ====================

    @BeforeEach
    public void setUp() {
        // Clear any existing system properties that might affect workspace detection
    }

    @AfterEach
    public void tearDown() {
        // Clean up any created directories
        cleanupWorkspaceDirectories();
    }

    @Test
    public void testGetWorkspaceDirectoryDefault() {
        // Test default behavior (should return ../jar/workspace relative to JAR location)
        String workspaceDir = MatFileUtils.getWorkspaceDirectory();
        assertTrue(workspaceDir.endsWith("jar/workspace"));
    }

    @Test
    public void testEnsureWorkspaceDirectoryExists() throws IOException {
        // Test that workspace directory is created
        MatFileUtils.ensureWorkspaceDirectoryExists();
        
        String workspaceDir = MatFileUtils.getWorkspaceDirectory();
        File dir = new File(workspaceDir);
        assertTrue(dir.exists());
        assertTrue(dir.isDirectory());
    }

    @Test
    public void testGenerateFilenameIncludesWorkspace() {
        String filename = MatFileUtils.genFilename("workspace");
        String workspaceDir = MatFileUtils.getWorkspaceDirectory();
        
        assertTrue(filename.startsWith(workspaceDir + "/"));
        assertTrue(filename.contains("workspace"));
        assertTrue(filename.endsWith(".mat"));
    }

    @Test
    public void testSaveWorkspaceToCorrectDirectory() throws IOException {
        // Create a test matrix
        Matrix testMatrix = new Matrix(2, 2);
        testMatrix.set(0, 0, 1.0);
        testMatrix.set(0, 1, 2.0);
        testMatrix.set(1, 0, 3.0);
        testMatrix.set(1, 1, 4.0);
        
        Map<String, Matrix> workspace = new HashMap<String, Matrix>();
        workspace.put("testMatrix", testMatrix);
        
        // Generate filename and save workspace
        String filename = MatFileUtils.genFilename("workspace");
        MatFileUtils.ensureWorkspaceDirectoryExists();
        MatFileUtils.saveWorkspace(workspace, filename);
        
        // Verify file was created in correct location
        File savedFile = new File(filename);
        assertTrue(savedFile.exists());
        assertEquals("workspace", savedFile.getParentFile().getName());
        
        // Clean up
        savedFile.delete();
    }

    private void cleanupWorkspaceDirectories() {
        // Clean up test directories
        String workspaceDir = MatFileUtils.getWorkspaceDirectory();
        deleteDirectory(new File(workspaceDir));
    }

    private void deleteDirectory(File directory) {
        if (directory.exists()) {
            File[] files = directory.listFiles();
            if (files != null) {
                for (File file : files) {
                    if (file.isDirectory()) {
                        deleteDirectory(file);
                    } else {
                        file.delete();
                    }
                }
            }
            directory.delete();
        }
    }
}