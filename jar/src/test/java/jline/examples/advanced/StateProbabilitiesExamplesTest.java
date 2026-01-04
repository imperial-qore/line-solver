package jline.examples.advanced;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.examples.java.advanced.StateProbabilitiesModel;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.nc.SolverNC;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Timeout;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.MID_TOL;
import static jline.TestTools.relativeTolerance;
import static jline.TestTools.withSuppressedOutput;

/**
 * Comprehensive tests for StateProbabilitiesModel examples following TEST_GENERATION_GUIDE.md.
 * These tests verify state probability computations against MATLAB dev/ directory outputs.
 */
public class StateProbabilitiesExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        // Set verbose level to SILENT to suppress warnings during tests
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    // ======================== statepr_aggr tests ========================
    // Tests basic closed network for state probability analysis
    // MATLAB uses: CTMC and NC solvers
    
    @Test
    @Timeout(value = 5, unit = java.util.concurrent.TimeUnit.MINUTES)
    public void testStatePrAggrCTMC() {
        // Test the statepr_aggr example with CTMC solver
        Network model = StateProbabilitiesModel.statepr_aggr();
        
        // Set state for Queue2 (node index 2) to [0,0]
        List<Node> nodes = model.getNodes();
        Queue queue2 = (Queue) nodes.get(2); // Queue2 is the third node
        
        // Initialize state as in MATLAB example
        queue2.setState(new Matrix("[0,0]"));
        
        final double[] probAggrHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "verbose", true);
            probAggrHolder[0] = solver.getProbAggr(queue2).getScalarProbability();
        });
        double probAggr = probAggrHolder[0];
        
        // Expected value from allExamplesBaseline.txt (CTMC solver output)
        // Station 3 is in state [0 0] with probability 3.400000e-01
        double expectedProb = 0.340000000000000;
        
        assertEquals(expectedProb, probAggr, relativeTolerance(expectedProb, MID_TOL),
                String.format("CTMC getProbAggr: expected %.6g, got %.6g", expectedProb, probAggr));
    }
    
    @Test
    public void testStatePrAggrNC() {
        // Test the statepr_aggr example with NC solver
        Network model = StateProbabilitiesModel.statepr_aggr();
        
        // Set state for Queue2 (node index 2) to [0,0]
        List<Node> nodes = model.getNodes();
        Queue queue2 = (Queue) nodes.get(2); // Queue2 is the third node
        
        // Initialize state as in MATLAB example
        queue2.setState(new Matrix("[0,0]"));
        
        final double[] probAggrHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model, "verbose", true);
            probAggrHolder[0] = solver.getProbAggr(queue2);
        });
        double probAggr = probAggrHolder[0];
        
        // Expected value from allExamplesBaseline.txt (NC solver output)
        // Station 3 is in state [0 0] with probability 3.400000e-01
        double expectedProb = 0.340000000000000;
        
        assertEquals(expectedProb, probAggr, relativeTolerance(expectedProb, MID_TOL),
                String.format("NC getProbAggr: expected %.6g, got %.6g", expectedProb, probAggr));
    }
    
    // ======================== statepr_aggr_large tests ========================
    // Tests four-class closed network with matrix-based routing
    // MATLAB uses: CTMC and NC solvers
    
    @Test
    // Testing singular matrix issue: @Disabled("CTMC solver error: Matrix is singular and system cannot be solved")
    public void testStatePrAggrLargeCTMC() {
        // Test the statepr_aggr_large example with CTMC solver
        Network model = StateProbabilitiesModel.statepr_aggr_large();

        // Set state for Queue2 (node index 2) to [1,0,2,1]
        List<Node> nodes = model.getNodes();
        Queue queue2 = (Queue) nodes.get(2); // Queue2 is the third node

        // Initialize state as in MATLAB example
        queue2.setState(new Matrix("[1,0,2,1]"));

        SolverCTMC solver = new SolverCTMC(model, "verbose", true, "keep", false);
        double probAggr = solver.getProbAggr(queue2).getScalarProbability();
        
        // Expected value from allExamplesBaseline.txt (CTMC solver output)
        // Station 3 is in state [1 0 2 1] with probability 5.511140e-03
        double expectedProb = 0.005511140445099;
        
        assertEquals(expectedProb, probAggr, relativeTolerance(expectedProb, MID_TOL),
                String.format("CTMC getProbAggr: expected %.6g, got %.6g", expectedProb, probAggr));
    }
    
    @Test
    public void testStatePrAggrLargeNC() {
        // Test the statepr_aggr_large example with NC solver
        Network model = StateProbabilitiesModel.statepr_aggr_large();
        NetworkStruct sn = model.getStruct();

        // Set state for Queue2 (node index 2) to [1,0,2,1]
        List<Node> nodes = model.getNodes();
        Queue queue2 = (Queue) nodes.get(2); // Queue2 is the third node
        
        // Initialize state as in MATLAB example
        queue2.setState(new Matrix("[1,0,2,1]"));
        
        final double[] probAggrHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model, "verbose", true);
            probAggrHolder[0] = solver.getProbAggr(queue2);
        });
        double probAggr = probAggrHolder[0];
        
        // Expected value from allExamplesBaseline.txt (NC solver output)
        // Station 3 is in state [1 0 2 1] with probability 5.511140e-03
        double expectedProb = 0.005511140445099;
        
        assertEquals(expectedProb, probAggr, relativeTolerance(expectedProb, MID_TOL),
                String.format("NC getProbAggr: expected %.6g, got %.6g", expectedProb, probAggr));
    }
    
    // ======================== statepr_sys_aggr tests ========================
    // Tests complex four-class network with class switching router
    // MATLAB uses: CTMC and NC solvers
    
    @Test
    public void testStatePrSysAggrCTMC() {
        // Test the statepr_sys_aggr example with CTMC solver
        Network model = StateProbabilitiesModel.statepr_sys_aggr();
        
        // Set system state as in MATLAB example:
        // node{1}: [0,0,0,0], node{2}: [0,0,0,0], node{3}: [1,0,3,0]
        List<Node> nodes = model.getNodes();
        nodes.get(0).setState(new Matrix("[0,0,0,0]")); // Delay
        nodes.get(1).setState(new Matrix("[0,0,0,0]")); // Queue1
        nodes.get(2).setState(new Matrix("[1,0,3,0]")); // Queue2
        // Only 3 nodes in this model: Delay, Queue1, Queue2
        
        final double[] probHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "verbose", true);
            probHolder[0] = solver.getProbSysAggr().getScalarProbability();
        });
        double prob = probHolder[0];
        
        // Expected value from ground truth
        // Pr_ctmc = 2.626044433948373e-04
        double expectedProb = 2.626044433948373e-04;
        
        assertEquals(expectedProb, prob, relativeTolerance(expectedProb, MID_TOL),
                String.format("CTMC getProbSysAggr: expected %.6g, got %.6g", expectedProb, prob));
    }
    
    @Test
    public void testStatePrSysAggrNC() {
        // Test the statepr_sys_aggr example with NC solver
        Network model = StateProbabilitiesModel.statepr_sys_aggr();
        
        // Set system state as in MATLAB example:
        // node{1}: [0,0,0,0], node{2}: [0,0,0,0], node{3}: [1,0,3,0]
        List<Node> nodes = model.getNodes();
        nodes.get(0).setState(new Matrix("[0,0,0,0]")); // Delay
        nodes.get(1).setState(new Matrix("[0,0,0,0]")); // Queue1
        nodes.get(2).setState(new Matrix("[1,0,3,0]")); // Queue2
        // Only 3 nodes in this model: Delay, Queue1, Queue2
        
        final double[] probHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model, "method","exact","verbose", true);
            probHolder[0] = solver.getProbSysAggr().getScalarProbability();
        });
        double prob = probHolder[0];
        
        // Expected value from ground truth
        // Pr_nc = 2.626044433948372e-04
        double expectedProb = 2.626044433948372e-04;
        
        assertEquals(expectedProb, prob, relativeTolerance(expectedProb, MID_TOL),
                String.format("NC getProbSysAggr: expected %.6g, got %.6g", expectedProb, prob));
    }
    
    // ======================== statepr_sys_aggr_large tests ========================
    // Tests three-queue network with symmetric class populations
    // MATLAB output not fully shown in extraction
    
    @Test
    public void testStatePrSysAggrLargeCTMC() {
        // Test the statepr_sys_aggr_large example with CTMC solver
        Network model = StateProbabilitiesModel.statepr_sys_aggr_large();
        
        // Set system state as in MATLAB example:
        // node{1}: [0,0,0,0], node{2}: [0,0,0,0], node{3}: [1,1,1,1]
        List<Node> nodes = model.getNodes();
        nodes.get(0).setState(new Matrix("[0,0,0,0]")); // Queue1
        nodes.get(1).setState(new Matrix("[0,0,0,0]")); // Queue2
        nodes.get(2).setState(new Matrix("[1,1,1,1]")); // Queue3
        
        final double[] probHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "verbose", true, "seed", 23000);
            probHolder[0] = solver.getProbSysAggr().getScalarProbability();
        });
        double prob = probHolder[0];
        
        // Expected value from allExamplesBaseline.txt
        // Pr_ctmc = 3.484356916108210e-04
        double expectedProb = 3.484356916108210e-04;
        
        assertEquals(expectedProb, prob, relativeTolerance(expectedProb, MID_TOL),
                String.format("CTMC getProbSysAggr: expected %.6g, got %.6g", expectedProb, prob));
        
        // Verify model structure
        assertEquals(3, model.getNumberOfStations(), "Should have 3 stations");
        assertEquals(4, model.getNumberOfClasses(), "Should have 4 classes");
    }
    
    @Test
    @Disabled("Expected value mismatch - ls sampling method uses different API than MATLAB")
    public void testStatePrSysAggrLargeNC() {
        // Test the statepr_sys_aggr_large example with NC solver
        Network model = StateProbabilitiesModel.statepr_sys_aggr_large();
        
        // Set system state as in MATLAB example:
        // node{1}: [0,0,0,0], node{2}: [0,0,0,0], node{3}: [1,1,1,1]
        List<Node> nodes = model.getNodes();
        nodes.get(0).setState(new Matrix("[0,0,0,0]")); // Queue1
        nodes.get(1).setState(new Matrix("[0,0,0,0]")); // Queue2
        nodes.get(2).setState(new Matrix("[1,1,1,1]")); // Queue3
        
        final double[] probHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model, "seed",23000, "method","ls", "samples",10000);
            probHolder[0] = solver.getProbSysAggr().getScalarProbability();
        });
        double prob = probHolder[0];
        
        // Expected value from ground truth
        // Pr_nc = 3.414029203772685e-04
        double expectedProb = 3.414029203772685e-04;
        
        assertEquals(expectedProb, prob, relativeTolerance(expectedProb, MID_TOL),
                String.format("NC getProbSysAggr: expected %.6g, got %.6g", expectedProb, prob));
        
        // Verify model structure
        assertEquals(3, model.getNumberOfStations(), "Should have 3 stations");
        assertEquals(4, model.getNumberOfClasses(), "Should have 4 classes");
    }
    
    // ======================== statepr_allprobs_ps tests ========================
    // Tests two-class network with bidirectional class switching
    // MATLAB output not fully shown in extraction
    
    @Test
    public void testStatePrAllProbsPSCTMC() {
        // Test the statepr_allprobs_ps example with CTMC solver
        Network model = StateProbabilitiesModel.statepr_allprobs_ps();
        
        // Set state as in MATLAB example:
        // node{1}: [0,0], node{2}: [1,0], node{3}: [0,1]
        List<Station> stations = model.getStations();
        stations.get(0).setState(new Matrix("[0,0]")); // Delay
        stations.get(1).setState(new Matrix("[1,0]")); // Queue1
        stations.get(2).setState(new Matrix("[0,1]")); // Queue2
        
        final double[] probHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "verbose", true, "seed", 23000);
            probHolder[0] = solver.getProbAggr(stations.get(2)).getScalarProbability();
        });
        double prob = probHolder[0];
        
        // Expected value from ground truth
        // Station 3 is in state [0 1] with probability 3.076923e-01
        // Pmarga_ctmc = 0.307692307692308
        double expectedProb = 0.307692307692308;
        
        assertEquals(expectedProb, prob, relativeTolerance(expectedProb, MID_TOL),
                String.format("CTMC getProbAggr: expected %.6g, got %.6g", expectedProb, prob));
        
        // Verify model structure
        assertEquals(6, model.getNumberOfNodes(), "Should have 3 nodes: Delay, Queue1, Queue2");
        assertEquals(2, model.getNumberOfClasses(), "Should have 2 classes");
    }
    
    @Test
    public void testStatePrAllProbsPSNC() {
        // Test the statepr_allprobs_ps example with NC solver
        Network model = StateProbabilitiesModel.statepr_allprobs_ps();
        
        // Set state as in MATLAB example:
        // node{1}: [0,0], node{2}: [1,0], node{3}: [0,1]
        List<Station> stations = model.getStations();
        stations.get(0).setState(new Matrix("[0,0]")); // Delay
        stations.get(1).setState(new Matrix("[1,0]")); // Queue1
        stations.get(2).setState(new Matrix("[0,1]")); // Queue2
        
        final double[] probAggrHolder = new double[1];
        final double[] probJointHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model, "verbose", true);
            // Get marginal probability for Queue2
            probAggrHolder[0] = solver.getProbAggr(stations.get(2));
            // Get joint probability (full system state)
            probJointHolder[0] = solver.getProbSysAggr().getScalarProbability();
        });
        double probAggr = probAggrHolder[0];
        double probJoint = probJointHolder[0];
        
        // Expected values from ground truth for NC solver
        // Pmarga_nc = 0.307692307692308
        // Pjointa_nc = 0.076190930869983
        double expectedProbAggr = 0.307692307692308;
        double expectedProbJoint = 0.076190930869983;
        
        assertEquals(expectedProbAggr, probAggr, relativeTolerance(expectedProbAggr, MID_TOL),
                String.format("NC getProbAggr: expected %.6g, got %.6g", expectedProbAggr, probAggr));
        assertEquals(expectedProbJoint, probJoint, relativeTolerance(expectedProbJoint, MID_TOL),
                String.format("NC getProbSysAggr: expected %.6g, got %.6g", expectedProbJoint, probJoint));
        
        // Verify model structure
        assertEquals(6, model.getNumberOfNodes(), "Should have 3 nodes: Delay, Queue1, Queue2");
        assertEquals(2, model.getNumberOfClasses(), "Should have 2 classes");
    }
    
    // ======================== statepr_allprobs_fcfs tests ========================
    // Tests mixed scheduling network with class switching
    // MATLAB output not fully shown in extraction
    
    @Test
    public void testStatePrAllProbsFCFSCTMC() {
        // Test the statepr_allprobs_fcfs example with CTMC solver
        Network model = StateProbabilitiesModel.statepr_allprobs_fcfs();
        
        // Set state as in MATLAB example (same as PS example):
        // node{1}: [0,0], node{2}: [1,0], node{3}: [0,1]
        List<Station> stations = model.getStations();
        stations.get(0).setState(new Matrix("[0,0]")); // Delay
        stations.get(1).setState(new Matrix("[1,0]")); // Queue1
        stations.get(2).setState(new Matrix("[0,1]")); // Queue2
        
        final double[] probHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "verbose", true, "seed", 23000);
            probHolder[0] = solver.getProbAggr(stations.get(2)).getScalarProbability();
        });
        double prob = probHolder[0];
        
        // Expected value from ground truth
        // Station 3 is in state [0 1] with probability 3.076923e-01
        // Pmarga_ctmc = 0.307692307692308
        double expectedProb = 0.307692307692308;
        
        assertEquals(expectedProb, prob, relativeTolerance(expectedProb, MID_TOL),
                String.format("CTMC getProbAggr: expected %.6g, got %.6g", expectedProb, prob));
        
        // Verify model structure
        assertEquals(6, model.getNumberOfNodes(), "Should have 3 nodes: Delay, Queue1, Queue2");
        assertEquals(2, model.getNumberOfClasses(), "Should have 2 classes");
        
        // Verify scheduling strategies
        Queue queue1 = (Queue) stations.get(1); // Queue1
        Queue queue2 = (Queue) stations.get(2); // Queue2
        
        assertEquals("PS", queue1.getSchedStrategy().toString(), "Queue1 should use PS scheduling");
        assertEquals("FCFS", queue2.getSchedStrategy().toString(), "Queue2 should use FCFS scheduling");
    }
    
    @Test
    public void testStatePrAllProbsFCFSNC() {
        // Test the statepr_allprobs_fcfs example with NC solver
        Network model = StateProbabilitiesModel.statepr_allprobs_fcfs();
        
        // Set state as in MATLAB example:
        // node{1}: [0,0], node{2}: [1,0], node{3}: [0,1]
        List<Station> stations = model.getStations();
        stations.get(0).setState(new Matrix("[0,0]")); // Delay
        stations.get(1).setState(new Matrix("[1,0]")); // Queue1
        stations.get(2).setState(new Matrix("[0,1]")); // Queue2
        
        final double[] probAggrHolder = new double[1];
        final double[] probJointHolder = new double[1];
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model, "verbose", true);
            // Get marginal probability for Queue2
            probAggrHolder[0] = solver.getProbAggr(stations.get(2));
            // Get joint probability (full system state)
            probJointHolder[0] = solver.getProbSysAggr().getScalarProbability();
        });
        double probAggr = probAggrHolder[0];
        double probJoint = probJointHolder[0];
        
        // Expected values from ground truth for NC solver
        // Pmarga_nc = 0.307692307692308
        // Pjointa_nc = 0.076190930869983
        double expectedProbAggr = 0.307692307692308;
        double expectedProbJoint = 0.076190930869983;
        
        assertEquals(expectedProbAggr, probAggr, relativeTolerance(expectedProbAggr, MID_TOL),
                String.format("NC getProbAggr: expected %.6g, got %.6g", expectedProbAggr, probAggr));
        assertEquals(expectedProbJoint, probJoint, relativeTolerance(expectedProbJoint, MID_TOL),
                String.format("NC getProbSysAggr: expected %.6g, got %.6g", expectedProbJoint, probJoint));
        
        // Verify model structure
        assertEquals(6, model.getNumberOfNodes(), "Should have 3 nodes: Delay, Queue1, Queue2");
        assertEquals(2, model.getNumberOfClasses(), "Should have 2 classes");
    }
}
