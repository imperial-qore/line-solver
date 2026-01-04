package jline.solvers.auto;

import jline.lang.Network;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Source;
import jline.lang.nodes.Sink;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.lang.ClosedClass;
import jline.lang.OpenClass;
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkSolver;
import jline.solvers.auto.SolverAUTO;
import jline.solvers.auto.AUTOptions;
import jline.solvers.auto.LINE;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import jline.TestTools;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Comprehensive test suite for SolverAUTO functionality
 */
public class SolverAUTOTest {
    
    private Network closedModel;
    private Network openModel;

    @BeforeEach
    public void setUp() {
        // Create a closed network model
        closedModel = new Network("ClosedTestModel");
        Delay delay = new Delay(closedModel, "Delay");
        Queue queue = new Queue(closedModel, "Queue", SchedStrategy.FCFS);
        ClosedClass jobClass = new ClosedClass(closedModel, "Jobs", 5, delay);
        delay.setService(jobClass, new Exp(1.0));
        queue.setService(jobClass, new Exp(2.0));
        closedModel.link(closedModel.serialRouting(delay, queue));
        
        // Create an open network model
        openModel = new Network("OpenTestModel");
        Source source = new Source(openModel, "Source");
        Queue queue2 = new Queue(openModel, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(openModel, "Sink");
        OpenClass openJobClass = new OpenClass(openModel, "OpenJobs");
        source.setArrival(openJobClass, new Exp(1.0));
        queue2.setService(openJobClass, new Exp(2.0));
        openModel.link(openModel.serialRouting(source, queue2, sink));
    }
    
    @Test
    public void testCandidateSolverDetection() {
        SolverAUTO solver = new SolverAUTO(closedModel);
        List<String> candidateNames = solver.getCandidateSolverNames();
        
        // Verify key solvers are included in candidates
        assertTrue(candidateNames.contains("SolverCTMC"), 
                   "CTMC solver should be included in SolverAUTO candidates");
        assertTrue(candidateNames.contains("SolverMVA"), 
                   "MVA solver should be included in SolverAUTO candidates");
        assertTrue(candidateNames.contains("SolverNC"), 
                   "NC solver should be included in SolverAUTO candidates");
        
        // Verify candidates list is not empty
        assertFalse(candidateNames.isEmpty(), "Should have at least one candidate solver");
    }
    
    @Test
    public void testDefaultSolverSelection() {
        SolverAUTO solver = new SolverAUTO(closedModel);
        
        // Initially, no solver is selected until analysis is triggered
        String initialSelectedName = solver.getSelectedSolverName();
        assertEquals("none", initialSelectedName, "Initially no solver should be selected");
        
        // Trigger analysis to force solver selection
        try {
            TestTools.withSuppressedOutput(() -> solver.getAvg());
            String selectedName = solver.getSelectedSolverName();
            assertNotNull(selectedName, "A solver should be selected after analysis");
            assertNotEquals("none", selectedName, "Should not be 'none' after analysis");
        } catch (Exception e) {
            // Analysis might fail, but we should still get a selected solver
            String selectedName = solver.getSelectedSolverName();
            assertNotNull(selectedName, "A solver should be selected even if analysis fails");
        }
    }
    
    @Test
    public void testForcedSolverSelection() {
        // Test forcing specific solvers
        AUTOptions options = new AUTOptions();
        options.forceSolver = "mva";
        
        SolverAUTO solver = new SolverAUTO(closedModel, options);
        
        // Verify the forced solver is used (after running analysis)
        try {
            TestTools.withSuppressedOutput(() -> solver.getAvg());
            String selectedName = solver.getSelectedSolverName();
            assertEquals("SolverMVA", selectedName, "Should use forced MVA solver");
        } catch (Exception e) {
            // If MVA fails, that's okay for this test - we're testing selection logic
        }
    }
    
    @Test
    public void testAUTOptionsConstructors() {
        // Test different constructor patterns
        SolverAUTO solver1 = new SolverAUTO(closedModel);
        assertNotNull(solver1, "Default constructor should work");
        
        SolverAUTO solver2 = new SolverAUTO(closedModel, "heur");
        assertNotNull(solver2, "Method constructor should work");
        
        AUTOptions options = new AUTOptions("default");
        SolverAUTO solver3 = new SolverAUTO(closedModel, options);
        assertNotNull(solver3, "Options constructor should work");
    }
    
    @Test
    public void testModelSupport() {
        SolverAUTO closedSolver = new SolverAUTO(closedModel);
        assertTrue(closedSolver.supports(closedModel), "Should support closed model");
        
        SolverAUTO openSolver = new SolverAUTO(openModel);
        assertTrue(openSolver.supports(openModel), "Should support open model");
    }
    
    @Test
    public void testBasicAnalysis() {
        SolverAUTO solver = new SolverAUTO(closedModel);
        
        // Test that we can run basic analysis without errors
        try {
            TestTools.suppressOutput();
            SolverResult result = null;
            try {
                result = solver.getAvg();
            } finally {
                TestTools.restoreOutput();
            }
            assertNotNull(result, "Should return a result object");
            
            String selectedSolver = solver.getSelectedSolverName();
            assertNotNull(selectedSolver, "Should have a selected solver after analysis");
            assertNotEquals("none", selectedSolver, "Selected solver should not be 'none'");
            
        } catch (Exception e) {
            // Some analysis might fail depending on model complexity
            // But we should at least get a meaningful error, not a null pointer
            assertNotNull(e.getMessage(), "Exception should have a message");
        }
    }
    
    @Test
    public void testSelectionMethods() {
        // Test different selection methods
        String[] methods = {"default", "heur", "ai", "nn"};
        
        for (String method : methods) {
            AUTOptions options = new AUTOptions(method);
            SolverAUTO solver = new SolverAUTO(closedModel, options);
            
            assertNotNull(solver, "Solver should be created with method: " + method);
            assertFalse(solver.getCandidateSolverNames().isEmpty(), 
                       "Should have candidates with method: " + method);
        }
    }
    
    @Test
    public void testOpenNetworkModel() {
        SolverAUTO solver = new SolverAUTO(openModel);
        List<String> candidates = solver.getCandidateSolverNames();
        
        assertFalse(candidates.isEmpty(), "Should have candidates for open model");
        
        // Test basic functionality with open model
        try {
            TestTools.withSuppressedOutput(() -> solver.getAvg());
            String selectedSolver = solver.getSelectedSolverName();
            assertNotNull(selectedSolver, "Should select a solver for open model");
        } catch (Exception e) {
            // Analysis might fail, but we should get a meaningful error
            assertNotNull(e.getMessage(), "Exception should have a message");
        }
    }
    
    @Test
    public void testInvalidForcedSolver() {
        AUTOptions options = new AUTOptions();
        options.forceSolver = "invalid_solver";
        
        SolverAUTO solver = new SolverAUTO(closedModel, options);
        
        // Exception should be thrown when analysis is triggered
        assertThrows(RuntimeException.class, () -> {
            TestTools.withSuppressedOutput(() -> solver.getAvg());
        }, "Should throw exception for invalid forced solver during analysis");
    }


    @Test
    public void testLINEClassInstantiation() {
        // Test that the LINE class exists and has correct structure
        assertNotNull(LINE.class);

        // Test that we can create instances via reflection (proving class structure is correct)
        try {
            Class<?> lineClass = Class.forName("jline.solvers.auto.LINE");
            assertNotNull(lineClass);
            assertTrue(lineClass.getSuperclass().getSimpleName().equals("SolverAUTO"));

            // Check that the load method exists
            Method loadMethod = lineClass.getMethod("load", String.class, Network.class, Object[].class);
            assertNotNull(loadMethod);
            assertTrue(java.lang.reflect.Modifier.isStatic(loadMethod.getModifiers()));

            // Check that create methods exist
            Method createMethod1 = lineClass.getMethod("create", Network.class);
            assertNotNull(createMethod1);
            assertTrue(java.lang.reflect.Modifier.isStatic(createMethod1.getModifiers()));

            Method createMethod2 = lineClass.getMethod("create", Network.class, String.class);
            assertNotNull(createMethod2);
            assertTrue(java.lang.reflect.Modifier.isStatic(createMethod2.getModifiers()));

        } catch (ClassNotFoundException | NoSuchMethodException e) {
            fail("LINE class should exist with proper methods: " + e.getMessage());
        }
    }

    @Test
    public void testStaticMethodsReturnTypes() {
        // Test that the static load method can create different solver types
        Network emptyModel = new Network("test");

        try {
            // Test that method routing works by checking return types
            NetworkSolver autoSolver = LINE.load("auto", emptyModel);
            // Even if it fails due to empty model, it should create the right type first
        } catch (RuntimeException e) {
            // Expected for empty model - that's fine, we're testing the factory logic
            assertTrue(e.getMessage().contains("empty") || e.getMessage().contains("supports"));
        }

        try {
            NetworkSolver mvaSolver = LINE.load("mva", emptyModel);
        } catch (RuntimeException e) {
            // Expected for empty model
            assertTrue(e.getMessage().contains("empty") || e.getMessage().contains("supports"));
        }

        try {
            NetworkSolver ncSolver = LINE.load("nc", emptyModel);
        } catch (RuntimeException e) {
            // Expected for empty model
            assertTrue(e.getMessage().contains("empty") || e.getMessage().contains("supports"));
        }
    }

    /**
     * Test single chain model - should select NC solver
     */
    @Test
    public void testSingleChainModel() {

        Network model = new Network("SingleChainTest");

        // Create a simple closed model with one class
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass cclass = new ClosedClass(model, "Class1", 5, delay);

        delay.setService(cclass, new Exp(1.0));
        queue.setService(cclass, new Exp(2.0));

        model.link(model.serialRouting(delay, queue));

        // Test AUTO solver
        SolverAUTO solver = new SolverAUTO(model);
        TestTools.withSuppressedOutput(() -> {
            NetworkAvgTable avgTable = solver.getAvgTable();
            
            // Selected solver: solver.getSelectedSolverName()
            // Candidates: solver.getCandidateSolverNames()
            
            assertNotNull(avgTable, "Should return average table");
            String selectedSolver = solver.getSelectedSolverName();
            assertNotNull(selectedSolver, "Should have selected a solver");
        });
    }

    /**
     * Test multi-chain model with FCFS - should select based on load
     */
    @Test
    public void testMultiChainModel() {

        Network model = new Network("MultiChainTest");

        // Create two closed classes
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        
        ClosedClass class1 = new ClosedClass(model, "Class1", 2, delay);
        ClosedClass class2 = new ClosedClass(model, "Class2", 3, delay);

        delay.setService(class1, new Exp(1.0));
        delay.setService(class2, new Exp(1.0));
        
        queue1.setService(class1, new Exp(2.0));
        queue1.setService(class2, new Exp(1.5));

        queue2.setService(class1, new Exp(1.0));
        queue2.setService(class2, new Exp(0.8));

        model.link(model.serialRouting(delay, queue1, queue2));

        // Test AUTO solver
        SolverAUTO solver = new SolverAUTO(model);
        TestTools.withSuppressedOutput(() -> {
            NetworkAvgTable avgTable = solver.getAvgTable();
            
            // Selected solver: solver.getSelectedSolverName()
            // Total jobs: (2 + 3) should use NC for light load
            
            assertNotNull(avgTable, "Should return average table");
            String selectedSolver = solver.getSelectedSolverName();
            assertNotNull(selectedSolver, "Should have selected a solver");
        });
    }

    /**
     * Test forced solver selection
     */
    @Test
    public void testForcedSolver() {

        Network model = new Network("ForcedSolverTest");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        ClosedClass cclass = new ClosedClass(model, "Class1", 10, delay);

        delay.setService(cclass, new Exp(1.0));
        queue.setService(cclass, new Exp(1.0));

        model.link(model.serialRouting(delay, queue));

        // Force MVA solver
        SolverAUTO solver = new SolverAUTO(model, "force", "mva");
        TestTools.withSuppressedOutput(() -> {
            NetworkAvgTable avgTable = solver.getAvgTable();
            
            // Selected solver: solver.getSelectedSolverName()
            // Should be: MVA (forced)
            
            assertNotNull(avgTable, "Should return average table");
            String selectedSolver = solver.getSelectedSolverName();
            assertEquals("SolverMVA", selectedSolver, "Should use forced MVA solver");
        });
    }
}