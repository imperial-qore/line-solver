package jline.examples.advanced;

import jline.examples.java.advanced.RandomEnvironmentModel;
import jline.lang.ClosedClass;
import jline.lang.Environment;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.env.SolverENV;
import jline.solvers.fluid.SolverFLD;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.lang.constant.SolverType;
import jline.GlobalConstants;
import jline.VerboseLevel;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import java.util.List;

import static jline.TestTools.*;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for Random Environment examples.
 *
 * Expected values come from allExamplesBaseline.txt MATLAB execution.
 * Note: allExamplesBaseline.txt only shows QLen, Util, RespT, and Tput for renv examples.
 * ResidT and ArvR are not provided, so these are set to 0.0 in the expected arrays.
 */
public class RandomEnvExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        // Set verbose level to SILENT to suppress warnings during tests
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /**
     * Test for renv_twostages_repairmen example using SolverENV.
     * 
     * This example models a system with two environmental stages (UP and DOWN).
     * The queueing network has 2 stations (Queue1 and Queue2) with service rates that
     * change based on the environmental stage.
     * 
     */
    @Test
    public void testRenvTwoStagesRepairmenEnv() {
        // Create the model
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            Environment envModel = RandomEnvironmentModel.renv_twostages_repairmen();
            int E = envModel.getEnsemble().size();
            
            // Create solver options
            SolverOptions options = new SolverOptions(SolverType.ENV);
            options.iter_tol = 0.01;
            options.timespan[0] = 0;
            options.verbose = VerboseLevel.SILENT;
            
            SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
            // Same fluid configuration as renv_twostages_repairmen.m: FLD defaults
            // (stiff, no step cap) over [0, 1e3]
            fluidOptions.timespan[1] = 1000;
            
            // Create solvers for each stage
            NetworkSolver[] solvers = new NetworkSolver[E];
            for (int e = 0; e < E; e++) {
                solvers[e] = new SolverFLD(envModel.getModel(e));
                solvers[e].options = fluidOptions;
            }
            
            SolverENV solver = new SolverENV(envModel, solvers, options);
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Re-recorded 2026-08-11, after the ENV exit average became a converged
        // quadrature on the stage ODE grid rather than a right-endpoint sum
        // (f08f6b8d7); every renv row moved in all four codebases. Current
        // MATLAB ground truth, re-measured 2026-08-11 on this tree, is
        // QLen=[0.55549349, 0.44450651], Util the same,
        // RespT=[0.80361862, 0.64288334], Tput=[0.69124019, 0.69142640].
        // Re-recorded against it 2026-08-11, after the stage handoff started
        // weighting the exit average by the stage SOJOURN rather than the e->h
        // arc clock: the two throughputs, which a closed cycle admits only one
        // of, closed from 1.4% apart to 1.3e-4, and QLen[0] from 0.59% off
        // MATLAB to 0.0057%.
        // Order: Queue1(Class1), Queue2(Class1)
        double[] expectedQLen = {0.55552531791372340, 0.44447468208627583};
        double[] expectedUtil = {0.55552531791372340, 0.44447468208627583};
        double[] expectedRespT = {0.80359532011849380, 0.64287063057591070};
        double[] expectedResidT = {0.0, 0.0};
        double[] expectedArvR = {0.0, 0.0};
        double[] expectedTput = {0.69129984210436750, 0.69139055502986140};
        
        // Verify table size
        assertEquals(2, avgTable.getQLen().size(), 
            "Expected 2 entries (2 stations × 1 class)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    /**
     * Test for renv_threestages_repairmen example using SolverENV.
     * 
     * This example models a system with three environmental stages (UP, DOWN, FAST).
     * The stages have circular transitions using Erlang distributions.
     * The queueing network has 2 stations with stage-dependent service rates.
     * 
     * This test focuses on the infinitesimal generator computation to match the MATLAB example.
     */
    @Test
    public void testRenvThreeStagesRepairmenEnv() {
        // Create the model with CTMC solver configuration
        final SolverENV.EnvGeneratorResult[] generatorResultHolder = new SolverENV.EnvGeneratorResult[1];
        withSuppressedOutput(() -> {
            // Get the pre-configured model from RandomEnvironmentModel
            Environment envModel = RandomEnvironmentModel.renv_threestages_repairmen();

            // Get the number of stages from the model
            int E = envModel.getEnsemble().size();
            
            // Configure environment solver options  
            SolverOptions envOptions = new SolverOptions(SolverType.ENV);
            envOptions.iter_max = 100;
            envOptions.iter_tol = 0.05;
            envOptions.method = "default";
            envOptions.timespan[0] = 0;
            envOptions.timespan[1] = Double.POSITIVE_INFINITY;
            
            // Configure CTMC solver options
            SolverOptions ctmcOptions = new SolverOptions(SolverType.CTMC);
            ctmcOptions.stiff = false;
            ctmcOptions.verbose = VerboseLevel.SILENT;
            ctmcOptions.timespan[0] = 0;
            ctmcOptions.timespan[1] = Double.POSITIVE_INFINITY;
            
            // Use CTMC solvers instead of Fluid solvers for this test
            
            // Create CTMC solvers for each stage
            NetworkSolver[] solvers = new NetworkSolver[E];
            for (int e = 0; e < E; e++) {
                solvers[e] = new SolverCTMC(envModel.getModel(e));
                solvers[e].options = ctmcOptions;
            }
            
            // Create SolverENV with CTMC solvers
            SolverENV solver = new SolverENV(envModel, solvers, envOptions);
            generatorResultHolder[0] = solver.getGenerator();
        });
        SolverENV.EnvGeneratorResult generatorResult = generatorResultHolder[0];
        
        // Check if results are computed
        assertNotNull(generatorResult);
        assertNotNull(generatorResult.renvInfGen);
        assertNotNull(generatorResult.stageInfGen);

        // Verify that we have 3 stages
        assertEquals(3, generatorResult.stageInfGen.length, 
            "Expected 3 environmental stages");
            
        // Verify each stage generator is 3x3
        for (int i = 0; i < 3; i++) {
            assertEquals(3, generatorResult.stageInfGen[i].getNumRows(),
                "Stage " + i + " generator should have 3 rows");
            assertEquals(3, generatorResult.stageInfGen[i].getNumCols(),
                "Stage " + i + " generator should have 3 columns");
        }
        
        // Check stageInfGen{1} (Stage 1) values
        Matrix stage1 = generatorResult.stageInfGen[0];
        assertEquals(-1.0, stage1.get(0, 0), FINE_TOL, "Stage 1: (1,1) should be -1");
        assertEquals(3.0, stage1.get(1, 0), FINE_TOL, "Stage 1: (2,1) should be 3");
        assertEquals(1.0, stage1.get(0, 1), FINE_TOL, "Stage 1: (1,2) should be 1");
        assertEquals(-4.0, stage1.get(1, 1), FINE_TOL, "Stage 1: (2,2) should be -4");
        assertEquals(6.0, stage1.get(2, 1), FINE_TOL, "Stage 1: (3,2) should be 6");
        assertEquals(1.0, stage1.get(1, 2), FINE_TOL, "Stage 1: (2,3) should be 1");
        assertEquals(-6.0, stage1.get(2, 2), FINE_TOL, "Stage 1: (3,3) should be -6");
        
        // Check stageInfGen{2} (Stage 2) values
        Matrix stage2 = generatorResult.stageInfGen[1];
        assertEquals(-2.0, stage2.get(0, 0), FINE_TOL, "Stage 2: (1,1) should be -2");
        assertEquals(2.0, stage2.get(1, 0), FINE_TOL, "Stage 2: (2,1) should be 2");
        assertEquals(2.0, stage2.get(0, 1), FINE_TOL, "Stage 2: (1,2) should be 2");
        assertEquals(-4.0, stage2.get(1, 1), FINE_TOL, "Stage 2: (2,2) should be -4");
        assertEquals(4.0, stage2.get(2, 1), FINE_TOL, "Stage 2: (3,2) should be 4");
        assertEquals(2.0, stage2.get(1, 2), FINE_TOL, "Stage 2: (2,3) should be 2");
        assertEquals(-4.0, stage2.get(2, 2), FINE_TOL, "Stage 2: (3,3) should be -4");
        
        // Check stageInfGen{3} (Stage 3) values
        Matrix stage3 = generatorResult.stageInfGen[2];
        assertEquals(-3.0, stage3.get(0, 0), FINE_TOL, "Stage 3: (1,1) should be -3");
        assertEquals(1.0, stage3.get(1, 0), FINE_TOL, "Stage 3: (2,1) should be 1");
        assertEquals(3.0, stage3.get(0, 1), FINE_TOL, "Stage 3: (1,2) should be 3");
        assertEquals(-4.0, stage3.get(1, 1), FINE_TOL, "Stage 3: (2,2) should be -4");
        assertEquals(2.0, stage3.get(2, 1), FINE_TOL, "Stage 3: (3,2) should be 2");
        assertEquals(3.0, stage3.get(1, 2), FINE_TOL, "Stage 3: (2,3) should be 3");
        assertEquals(-2.0, stage3.get(2, 2), FINE_TOL, "Stage 3: (3,3) should be -2");


        // Verify the infinitesimal generator dimensions match expected (36x36 sparse)
        assertEquals(36, generatorResult.renvInfGen.getNumRows(),
                "Expected generator matrix to have 36 rows");
        assertEquals(36, generatorResult.renvInfGen.getNumCols(),
                "Expected generator matrix to have 36 columns");

        // Check key values in the main renvInfGen matrix (36x36)
        Matrix renv = generatorResult.renvInfGen;
        
        // Check diagonal values and some key off-diagonal values
        // Note: MATLAB uses 1-based indexing, Java uses 0-based
        assertEquals(-4.0, renv.get(0, 0), FINE_TOL, "renvInfGen: (1,1) should be -4");
        assertEquals(3.0, renv.get(3, 0), FINE_TOL, "renvInfGen: (4,1) should be 3");
        assertEquals(4.0, renv.get(27, 0), FINE_TOL, "renvInfGen: (28,1) should be 4");
        assertEquals(3.0, renv.get(0, 1), FINE_TOL, "renvInfGen: (1,2) should be 3");
        assertEquals(-4.0, renv.get(1, 1), FINE_TOL, "renvInfGen: (2,2) should be -4");
        assertEquals(3.0, renv.get(4, 1), FINE_TOL, "renvInfGen: (5,2) should be 3");
        
        // Check more diagonal values
        assertEquals(-7.0, renv.get(3, 3), FINE_TOL, "renvInfGen: (4,4) should be -7");
        assertEquals(-9.0, renv.get(6, 6), FINE_TOL, "renvInfGen: (7,7) should be -9");
        assertEquals(-7.0, renv.get(9, 9), FINE_TOL, "renvInfGen: (10,10) should be -7");
        assertEquals(-9.0, renv.get(14, 14), FINE_TOL, "renvInfGen: (15,15) should be -9");
        assertEquals(-6.0, renv.get(35, 35), FINE_TOL, "renvInfGen: (36,36) should be -6");
        
        // Count non-zeros
        int nonZeros = 0;
        for (int i = 0; i < renv.getNumRows(); i++) {
            for (int j = 0; j < renv.getNumCols(); j++) {
                if (Math.abs(renv.get(i, j)) > FINE_TOL) {
                    nonZeros++;
                }
            }
        }
        assertEquals(120, nonZeros, "Expected 120 non-zero entries in renvInfGen");

    }
    
    /**
     * Test for renv_fourstages_repairmen example using SolverENV.
     * 
     * This example models a complex system with four environmental stages (UP, DOWN, FAST, SLOW).
     * The queueing network has 3 stations but only 2 are shown in the output (Queue1 and Queue2).
     * The model uses Coxian transitions between stages with SCV=0.5.
     * 
     */
    @Test
    public void testRenvFourStagesRepairmenEnv() {
        // Create the model
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            Environment envModel = RandomEnvironmentModel.renv_fourstages_repairmen();
            int E = envModel.getEnsemble().size();
            
            // Create solver options
            SolverOptions options = new SolverOptions(SolverType.ENV);
            options.iter_tol = 0.05;
            options.timespan[0] = 0;
            options.verbose = VerboseLevel.SILENT;
            
            // Same fluid configuration as renv_fourstages_repairmen.m: FLD defaults
            // (stiff, no step cap) over [0, Inf]
            SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
            fluidOptions.verbose = VerboseLevel.SILENT;
            
            // Create solvers for each stage
            NetworkSolver[] solvers = new NetworkSolver[E];
            for (int e = 0; e < E; e++) {
                solvers[e] = new SolverFLD(envModel.getModel(e));
                solvers[e].options = fluidOptions;
            }
            
            SolverENV solver = new SolverENV(envModel, solvers, options);
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Re-recorded 2026-08-11, after the ENV exit average became a converged
        // quadrature on the stage ODE grid rather than a right-endpoint sum
        // (f08f6b8d7); every renv row moved in all four codebases. Current
        // MATLAB ground truth, re-measured 2026-08-11 on this tree, is
        // QLen=[0.43985144, 29.56014856], Util=[0.43985144, 1],
        // RespT=[0.45267701, 29.56014856], Tput=[0.97166727, 1].
        // Re-recorded against it 2026-08-11, with the sojourn-weighted stage
        // handoff: QLen[0] is now 0.0076% off MATLAB, where the golden it
        // replaces was 0.28% off.
        // Order: Queue1(Class1), Queue2(Class1)
        double[] expectedQLen = {0.43981813634160316, 29.560181863658520};
        double[] expectedUtil = {0.43981813634160316, 1.0};
        double[] expectedRespT = {0.45265641920570404, 29.560181863658520};
        double[] expectedResidT = {0.0, 0.0};
        double[] expectedArvR = {0.0, 0.0};
        double[] expectedTput = {0.97163790831326600, 1.0};
        
        // Verify table size
        assertEquals(2, avgTable.getQLen().size(), 
            "Expected 2 entries (2 visible stations × 1 class)");

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
}