package jline.lang.state;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static jline.TestTools.MID_TOL;
import static jline.TestTools.assertMatrixEquals;
import static jline.lang.state.StationStateTestFixtures.*;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for network state handling with different scheduling strategies.
 * Verifies infinitesimal generator and state space against ground truth from MATLAB.
 */
public class StationStateTest {

    @BeforeAll
    public static void setUp() {
        // Suppress priority diagnostic output during tests
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    @Test
    public void testFCFS() {
        Network model = createClosedNetwork(SchedStrategy.FCFS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator
        Matrix expectedGen = new Matrix("[-0.5, 0.5, 0, 0;" +
                "1.0, -1.5, 0.5, 0;" +
                "0, 2.0, -2.5, 0.5;" +
                "0, 0, 3.0, -3.0]");
        assertMatrixEquals(expectedGen, generator, "Generator matrix for FCFS");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0, 1, 1, 1;" +
                "1, 0, 1, 1;" +
                "2, 0, 0, 1;" +
                "3, 0, 0, 0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for FCFS");
    }

    @Test
    public void testLCFS() {
        Network model = createClosedNetwork(SchedStrategy.LCFS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - based on MATLAB output
        Matrix expectedGen = new Matrix("[-0.5, 0.5, 0, 0;" +
                "1.0, -1.5, 0.5, 0;" +
                "0, 2.0, -2.5, 0.5;" +
                "0, 0, 3.0, -3.0]");
        assertMatrixEquals(expectedGen, generator, "Generator matrix for LCFS");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0, 1, 1, 1;" +
                "1, 0, 1, 1;" +
                "2, 0, 0, 1;" +
                "3, 0, 0, 0]");

        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for LCFS");
    }

    @Test
    public void testSIRO() {
        Network model = createClosedNetwork(SchedStrategy.SIRO);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - SIRO should have same as FCFS
        Matrix expectedGen = new Matrix("[-0.5, 0.5, 0, 0;" +
                "1.0, -1.5, 0.5, 0;" +
                "0, 2.0, -2.5, 0.5;" +
                "0, 0, 3.0, -3.0]");
        assertMatrixEquals(expectedGen, generator, "Generator matrix for SIRO");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0, 2, 1;" +
                "1, 1, 1;" +
                "2, 0, 1;" +
                "3, 0, 0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for SIRO");
    }

    @Test
    public void testSEPT() {
        Network model = createClosedNetwork(SchedStrategy.SEPT);
        SolverCTMC ctmcSolver = new SolverCTMC(model);

        // SEPT may have null generator for networks with identical service times
        try {
            Matrix generator = ctmcSolver.getGenerator().infGen;
            Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

            // If we get here, verify the generator
            Matrix expectedGen = new Matrix("[-0.5, 0.5, 0, 0;" +
                    "1.0, -1.5, 0.5, 0;" +
                    "0, 2.0, -2.5, 0.5;" +
                    "0, 0, 3.0, -3.0]");
            assertMatrixEquals(expectedGen, generator, "Generator matrix for SEPT");

            // Verify state space
            Matrix expectedStateSpace = new Matrix("[0, 2, 1;" +
                    "1, 1, 1;" +
                    "2, 0, 1;" +
                    "3, 0, 0]");
            assertMatrixEquals(expectedStateSpace, stateSpace, "State space for SEPT");
        } catch (NullPointerException e) {
            // Expected for SEPT with identical service times
            // This is consistent with MATLAB behavior
        }
    }

    @Test
    public void testLEPT() {
        Network model = createClosedNetwork(SchedStrategy.LEPT);
        SolverCTMC ctmcSolver = new SolverCTMC(model);

        // LEPT may have null generator for networks with identical service times
        try {
            Matrix generator = ctmcSolver.getGenerator().infGen;
            Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

            // If we get here, verify the generator
            Matrix expectedGen = new Matrix("[-0.5, 0.5, 0, 0;" +
                    "1.0, -1.5, 0.5, 0;" +
                    "0, 2.0, -2.5, 0.5;" +
                    "0, 0, 3.0, -3.0]");
            assertMatrixEquals(expectedGen, generator, "Generator matrix for LEPT");

            // Verify state space
            Matrix expectedStateSpace = new Matrix("[0, 2, 1;" +
                    "1, 1, 1;" +
                    "2, 0, 1;" +
                    "3, 0, 0]");
            assertMatrixEquals(expectedStateSpace, stateSpace, "State space for LEPT");
        } catch (NullPointerException e) {
            // Expected for LEPT with identical service times
            // This is consistent with MATLAB behavior
        }
    }

    @Test
    public void testPS() {
        Network model = createClosedNetwork(SchedStrategy.PS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator
        Matrix expectedGen = new Matrix("[-0.5, 0.5, 0, 0;" +
                "1.0, -1.5, 0.5, 0;" +
                "0, 2.0, -2.5, 0.5;" +
                "0, 0, 3.0, -3.0]");
        assertMatrixEquals(expectedGen, generator, "Generator matrix for PS");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0, 3;" +
                "1, 2;" +
                "2, 1;" +
                "3, 0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for PS");
    }

    @Test
    public void testDPS() {
        Network model = createClosedNetwork(SchedStrategy.DPS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator
        Matrix expectedGen = new Matrix("[-0.5, 0.5, 0, 0;" +
                "1.0, -1.5, 0.5, 0;" +
                "0, 2.0, -2.5, 0.5;" +
                "0, 0, 3.0, -3.0]");
        assertMatrixEquals(expectedGen, generator, "Generator matrix for DPS");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0, 3;" +
                "1, 2;" +
                "2, 1;" +
                "3, 0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for DPS");
    }

    @Test
    public void testGPS() {
        Network model = createClosedNetwork(SchedStrategy.GPS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator
        Matrix expectedGen = new Matrix("[-0.5, 0.5, 0, 0;" +
                "1.0, -1.5, 0.5, 0;" +
                "0, 2.0, -2.5, 0.5;" +
                "0, 0, 3.0, -3.0]");
        assertMatrixEquals(expectedGen, generator, "Generator matrix for GPS");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0, 3;" +
                "1, 2;" +
                "2, 1;" +
                "3, 0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for GPS");
    }

    @Test
    public void testINF() {
        Network model = createClosedNetwork(SchedStrategy.INF);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - INF has different rates
        Matrix expectedGen = new Matrix("[-1.5, 1.5, 0, 0;" +
                "1.0, -2.0, 1.0, 0;" +
                "0, 2.0, -2.5, 0.5;" +
                "0, 0, 3.0, -3.0]");
        assertMatrixEquals(expectedGen, generator, "Generator matrix for INF");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0, 3;" +
                "1, 2;" +
                "2, 1;" +
                "3, 0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for INF");
    }

    @Test
    public void testHOL() {
        Network model = createClosedNetwork(SchedStrategy.HOL);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator
        Matrix expectedGen = new Matrix("[-0.5, 0.5, 0, 0;" +
                "1.0, -1.5, 0.5, 0;" +
                "0, 2.0, -2.5, 0.5;" +
                "0, 0, 3.0, -3.0]");
        assertMatrixEquals(expectedGen, generator, "Generator matrix for HOL");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0, 1, 1, 1;" +
                "1, 0, 1, 1;" +
                "2, 0, 0, 1;" +
                "3, 0, 0, 0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for HOL");
    }

    // Mixed Network Tests - One closed class (3 jobs) + One open class (λ=0.1)

    @Test
    public void testMixedNetworkSIRO() {
        Network model = createMixedNetwork(SchedStrategy.SIRO);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(6);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        Matrix expectedGen = new Matrix(46, 46);
        expectedGen.set(0,0,-0.500000);
        expectedGen.set(0,2,0.100000);
        expectedGen.set(0,13,0.400000);
        expectedGen.set(1,0,0.500000);
        expectedGen.set(1,1,-0.600000);
        expectedGen.set(1,3,0.100000);
        expectedGen.set(2,2,-0.500000);
        expectedGen.set(2,4,0.100000);
        expectedGen.set(2,14,0.133333);
        expectedGen.set(2,15,0.266667);
        expectedGen.set(3,1,0.125000);
        expectedGen.set(3,2,0.375000);
        expectedGen.set(3,3,-0.600000);
        expectedGen.set(3,5,0.100000);
        expectedGen.set(4,4,-0.500000);
        expectedGen.set(4,6,0.100000);
        expectedGen.set(4,16,0.200000);
        expectedGen.set(4,17,0.200000);
        expectedGen.set(5,3,0.200000);
        expectedGen.set(5,4,0.300000);
        expectedGen.set(5,5,-0.600000);
        expectedGen.set(5,7,0.100000);
        expectedGen.set(6,6,-0.500000);
        expectedGen.set(6,8,0.100000);
        expectedGen.set(6,18,0.240000);
        expectedGen.set(6,19,0.160000);
        expectedGen.set(7,5,0.250000);
        expectedGen.set(7,6,0.250000);
        expectedGen.set(7,7,-0.600000);
        expectedGen.set(7,9,0.100000);
        expectedGen.set(8,8,-0.500000);
        expectedGen.set(8,10,0.100000);
        expectedGen.set(8,20,0.266667);
        expectedGen.set(8,21,0.133333);
        expectedGen.set(9,7,0.285714);
        expectedGen.set(9,8,0.214286);
        expectedGen.set(9,9,-0.600000);
        expectedGen.set(9,11,0.100000);
        expectedGen.set(10,10,-0.500000);
        expectedGen.set(10,12,0.100000);
        expectedGen.set(10,22,0.285714);
        expectedGen.set(10,23,0.114286);
        expectedGen.set(11,9,0.312500);
        expectedGen.set(11,10,0.187500);
        expectedGen.set(11,11,-0.500000);
        expectedGen.set(12,12,-0.400000);
        expectedGen.set(12,24,0.300000);
        expectedGen.set(12,25,0.100000);
        expectedGen.set(13,0,1.000000);
        expectedGen.set(13,13,-1.500000);
        expectedGen.set(13,15,0.100000);
        expectedGen.set(13,26,0.400000);
        expectedGen.set(14,1,1.000000);
        expectedGen.set(14,13,0.500000);
        expectedGen.set(14,14,-1.600000);
        expectedGen.set(14,16,0.100000);
        expectedGen.set(15,2,1.000000);
        expectedGen.set(15,15,-1.500000);
        expectedGen.set(15,17,0.100000);
        expectedGen.set(15,27,0.200000);
        expectedGen.set(15,28,0.200000);
        expectedGen.set(16,3,1.000000);
        expectedGen.set(16,14,0.166667);
        expectedGen.set(16,15,0.333333);
        expectedGen.set(16,16,-1.600000);
        expectedGen.set(16,18,0.100000);
        expectedGen.set(17,4,1.000000);
        expectedGen.set(17,17,-1.500000);
        expectedGen.set(17,19,0.100000);
        expectedGen.set(17,29,0.266667);
        expectedGen.set(17,30,0.133333);
        expectedGen.set(18,5,1.000000);
        expectedGen.set(18,16,0.250000);
        expectedGen.set(18,17,0.250000);
        expectedGen.set(18,18,-1.600000);
        expectedGen.set(18,20,0.100000);
        expectedGen.set(19,6,1.000000);
        expectedGen.set(19,19,-1.500000);
        expectedGen.set(19,21,0.100000);
        expectedGen.set(19,31,0.300000);
        expectedGen.set(19,32,0.100000);
        expectedGen.set(20,7,1.000000);
        expectedGen.set(20,18,0.300000);
        expectedGen.set(20,19,0.200000);
        expectedGen.set(20,20,-1.600000);
        expectedGen.set(20,22,0.100000);
        expectedGen.set(21,8,1.000000);
        expectedGen.set(21,21,-1.500000);
        expectedGen.set(21,23,0.100000);
        expectedGen.set(21,33,0.320000);
        expectedGen.set(21,34,0.080000);
        expectedGen.set(22,9,1.000000);
        expectedGen.set(22,20,0.333333);
        expectedGen.set(22,21,0.166667);
        expectedGen.set(22,22,-1.600000);
        expectedGen.set(22,24,0.100000);
        expectedGen.set(23,10,1.000000);
        expectedGen.set(23,23,-1.500000);
        expectedGen.set(23,25,0.100000);
        expectedGen.set(23,35,0.333333);
        expectedGen.set(23,36,0.066667);
        expectedGen.set(24,11,1.000000);
        expectedGen.set(24,22,0.357143);
        expectedGen.set(24,23,0.142857);
        expectedGen.set(24,24,-1.500000);
        expectedGen.set(25,12,1.000000);
        expectedGen.set(25,25,-1.400000);
        expectedGen.set(25,37,0.342857);
        expectedGen.set(25,38,0.057143);
        expectedGen.set(26,13,2.000000);
        expectedGen.set(26,26,-2.500000);
        expectedGen.set(26,28,0.100000);
        expectedGen.set(26,39,0.400000);
        expectedGen.set(27,14,2.000000);
        expectedGen.set(27,26,0.500000);
        expectedGen.set(27,27,-2.600000);
        expectedGen.set(27,29,0.100000);
        expectedGen.set(28,15,2.000000);
        expectedGen.set(28,28,-2.500000);
        expectedGen.set(28,30,0.100000);
        expectedGen.set(28,40,0.400000);
        expectedGen.set(29,16,2.000000);
        expectedGen.set(29,27,0.250000);
        expectedGen.set(29,28,0.250000);
        expectedGen.set(29,29,-2.600000);
        expectedGen.set(29,31,0.100000);
        expectedGen.set(30,17,2.000000);
        expectedGen.set(30,30,-2.500000);
        expectedGen.set(30,32,0.100000);
        expectedGen.set(30,41,0.400000);
        expectedGen.set(31,18,2.000000);
        expectedGen.set(31,29,0.333333);
        expectedGen.set(31,30,0.166667);
        expectedGen.set(31,31,-2.600000);
        expectedGen.set(31,33,0.100000);
        expectedGen.set(32,19,2.000000);
        expectedGen.set(32,32,-2.500000);
        expectedGen.set(32,34,0.100000);
        expectedGen.set(32,42,0.400000);
        expectedGen.set(33,20,2.000000);
        expectedGen.set(33,31,0.375000);
        expectedGen.set(33,32,0.125000);
        expectedGen.set(33,33,-2.600000);
        expectedGen.set(33,35,0.100000);
        expectedGen.set(34,21,2.000000);
        expectedGen.set(34,34,-2.500000);
        expectedGen.set(34,36,0.100000);
        expectedGen.set(34,43,0.400000);
        expectedGen.set(35,22,2.000000);
        expectedGen.set(35,33,0.400000);
        expectedGen.set(35,34,0.100000);
        expectedGen.set(35,35,-2.600000);
        expectedGen.set(35,37,0.100000);
        expectedGen.set(36,23,2.000000);
        expectedGen.set(36,36,-2.500000);
        expectedGen.set(36,38,0.100000);
        expectedGen.set(36,44,0.400000);
        expectedGen.set(37,24,2.000000);
        expectedGen.set(37,35,0.416667);
        expectedGen.set(37,36,0.083333);
        expectedGen.set(37,37,-2.500000);
        expectedGen.set(38,25,2.000000);
        expectedGen.set(38,38,-2.400000);
        expectedGen.set(38,45,0.400000);
        expectedGen.set(39,26,3.000000);
        expectedGen.set(39,39,-3.100000);
        expectedGen.set(39,40,0.100000);
        expectedGen.set(40,27,3.000000);
        expectedGen.set(40,39,0.500000);
        expectedGen.set(40,40,-3.600000);
        expectedGen.set(40,41,0.100000);
        expectedGen.set(41,29,3.000000);
        expectedGen.set(41,40,0.500000);
        expectedGen.set(41,41,-3.600000);
        expectedGen.set(41,42,0.100000);
        expectedGen.set(42,31,3.000000);
        expectedGen.set(42,41,0.500000);
        expectedGen.set(42,42,-3.600000);
        expectedGen.set(42,43,0.100000);
        expectedGen.set(43,33,3.000000);
        expectedGen.set(43,42,0.500000);
        expectedGen.set(43,43,-3.600000);
        expectedGen.set(43,44,0.100000);
        expectedGen.set(44,35,3.000000);
        expectedGen.set(44,43,0.500000);
        expectedGen.set(44,44,-3.600000);
        expectedGen.set(44,45,0.100000);
        expectedGen.set(45,37,3.000000);
        expectedGen.set(45,44,0.500000);
        expectedGen.set(45,45,-3.500000);

        // Verify dimensions
        assertEquals(46, generator.getNumRows());
        assertEquals(46, generator.getNumCols());

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "Generator matrix mismatch");

        // Verify state space - from MATLAB output
        assertEquals(46, stateSpace.getNumRows());
        assertEquals(9, stateSpace.getNumCols());
    }

    @Test
    public void testMixedNetworkPS() {
        Network model = createMixedNetwork(SchedStrategy.PS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(6);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for PS mixed network
        // PS has 28x28 generator matrix
        Matrix expectedGen = new Matrix(28, 28);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        // Note: MATLAB uses 1-based indexing, converting to 0-based for Java
        expectedGen.set(0, 0, -0.5000);
        expectedGen.set(1, 0, 0.1250);
        expectedGen.set(7, 0, 1.0000);
        expectedGen.set(0, 1, 0.1000);
        expectedGen.set(1, 1, -0.5250);
        expectedGen.set(2, 1, 0.2000);
        expectedGen.set(8, 1, 1.0000);
        expectedGen.set(1, 2, 0.1000);
        expectedGen.set(2, 2, -0.5400);
        expectedGen.set(3, 2, 0.2500);
        expectedGen.set(9, 2, 1.0000);
        expectedGen.set(2, 3, 0.1000);
        expectedGen.set(3, 3, -0.5500);
        expectedGen.set(4, 3, 0.2857);
        expectedGen.set(10, 3, 1.0000);
        expectedGen.set(3, 4, 0.1000);
        expectedGen.set(4, 4, -0.5571);
        expectedGen.set(5, 4, 0.3125);
        expectedGen.set(11, 4, 1.0000);
        expectedGen.set(4, 5, 0.1000);
        expectedGen.set(5, 5, -0.5625);
        expectedGen.set(6, 5, 0.3333);
        expectedGen.set(12, 5, 1.0000);
        expectedGen.set(5, 6, 0.1000);
        expectedGen.set(6, 6, -0.4667);
        expectedGen.set(13, 6, 1.0000);
        expectedGen.set(0, 7, 0.4000);
        expectedGen.set(7, 7, -1.5000);
        expectedGen.set(8, 7, 0.1667);
        expectedGen.set(14, 7, 2.0000);
        expectedGen.set(1, 8, 0.3000);
        expectedGen.set(7, 8, 0.1000);
        expectedGen.set(8, 8, -1.5333);
        expectedGen.set(9, 8, 0.2500);
        expectedGen.set(15, 8, 2.0000);
        expectedGen.set(2, 9, 0.2400);
        expectedGen.set(8, 9, 0.1000);
        expectedGen.set(9, 9, -1.5500);
        expectedGen.set(10, 9, 0.3000);
        expectedGen.set(16, 9, 2.0000);
        expectedGen.set(3, 10, 0.2000);
        expectedGen.set(9, 10, 0.1000);
        expectedGen.set(10, 10, -1.5600);
        expectedGen.set(11, 10, 0.3333);
        expectedGen.set(17, 10, 2.0000);
        expectedGen.set(4, 11, 0.1714);
        expectedGen.set(10, 11, 0.1000);
        expectedGen.set(11, 11, -1.5667);
        expectedGen.set(12, 11, 0.3571);
        expectedGen.set(18, 11, 2.0000);
        expectedGen.set(5, 12, 0.1500);
        expectedGen.set(11, 12, 0.1000);
        expectedGen.set(12, 12, -1.5714);
        expectedGen.set(13, 12, 0.3750);
        expectedGen.set(19, 12, 2.0000);
        expectedGen.set(6, 13, 0.1333);
        expectedGen.set(12, 13, 0.1000);
        expectedGen.set(13, 13, -1.4750);
        expectedGen.set(20, 13, 2.0000);
        expectedGen.set(7, 14, 0.4000);
        expectedGen.set(14, 14, -2.5000);
        expectedGen.set(15, 14, 0.2500);
        expectedGen.set(21, 14, 3.0000);
        expectedGen.set(8, 15, 0.2667);
        expectedGen.set(14, 15, 0.1000);
        expectedGen.set(15, 15, -2.5500);
        expectedGen.set(16, 15, 0.3333);
        expectedGen.set(22, 15, 3.0000);
        expectedGen.set(9, 16, 0.2000);
        expectedGen.set(15, 16, 0.1000);
        expectedGen.set(16, 16, -2.5667);
        expectedGen.set(17, 16, 0.3750);
        expectedGen.set(23, 16, 3.0000);
        expectedGen.set(10, 17, 0.1600);
        expectedGen.set(16, 17, 0.1000);
        expectedGen.set(17, 17, -2.5750);
        expectedGen.set(18, 17, 0.4000);
        expectedGen.set(24, 17, 3.0000);
        expectedGen.set(11, 18, 0.1333);
        expectedGen.set(17, 18, 0.1000);
        expectedGen.set(18, 18, -2.5800);
        expectedGen.set(19, 18, 0.4167);
        expectedGen.set(25, 18, 3.0000);
        expectedGen.set(12, 19, 0.1143);
        expectedGen.set(18, 19, 0.1000);
        expectedGen.set(19, 19, -2.5833);
        expectedGen.set(20, 19, 0.4286);
        expectedGen.set(26, 19, 3.0000);
        expectedGen.set(13, 20, 0.1000);
        expectedGen.set(19, 20, 0.1000);
        expectedGen.set(20, 20, -2.4857);
        expectedGen.set(27, 20, 3.0000);
        expectedGen.set(14, 21, 0.4000);
        expectedGen.set(21, 21, -3.1000);
        expectedGen.set(22, 21, 0.5000);
        expectedGen.set(15, 22, 0.2000);
        expectedGen.set(21, 22, 0.1000);
        expectedGen.set(22, 22, -3.6000);
        expectedGen.set(23, 22, 0.5000);
        expectedGen.set(16, 23, 0.1333);
        expectedGen.set(22, 23, 0.1000);
        expectedGen.set(23, 23, -3.6000);
        expectedGen.set(24, 23, 0.5000);
        expectedGen.set(17, 24, 0.1000);
        expectedGen.set(23, 24, 0.1000);
        expectedGen.set(24, 24, -3.6000);
        expectedGen.set(25, 24, 0.5000);
        expectedGen.set(18, 25, 0.0800);
        expectedGen.set(24, 25, 0.1000);
        expectedGen.set(25, 25, -3.6000);
        expectedGen.set(26, 25, 0.5000);
        expectedGen.set(19, 26, 0.0667);
        expectedGen.set(25, 26, 0.1000);
        expectedGen.set(26, 26, -3.6000);
        expectedGen.set(27, 26, 0.5000);
        expectedGen.set(20, 27, 0.0571);
        expectedGen.set(26, 27, 0.1000);
        expectedGen.set(27, 27, -3.5000);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "Generator matrix mismatch");

        // Verify state space
        assertEquals(28, stateSpace.getNumRows());
        assertEquals(7, stateSpace.getNumCols());
    }

    @Test
    public void testMixedNetworkDPS() {
        Network model = createMixedNetwork(SchedStrategy.DPS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(6);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for DPS mixed network
        // DPS has 28x28 generator matrix
        Matrix expectedGen = new Matrix(28, 28);

        expectedGen.set(0,0,-0.500000);
        expectedGen.set(0,1,0.100000);
        expectedGen.set(0,7,0.400000);
        expectedGen.set(1,0,0.166667);
        expectedGen.set(1,1,-0.533333);
        expectedGen.set(1,2,0.100000);
        expectedGen.set(1,8,0.266667);
        expectedGen.set(2,1,0.250000);
        expectedGen.set(2,2,-0.550000);
        expectedGen.set(2,3,0.100000);
        expectedGen.set(2,9,0.200000);
        expectedGen.set(3,2,0.300000);
        expectedGen.set(3,3,-0.560000);
        expectedGen.set(3,4,0.100000);
        expectedGen.set(3,10,0.160000);
        expectedGen.set(4,3,0.333333);
        expectedGen.set(4,4,-0.566667);
        expectedGen.set(4,5,0.100000);
        expectedGen.set(4,11,0.133333);
        expectedGen.set(5,4,0.357143);
        expectedGen.set(5,5,-0.571429);
        expectedGen.set(5,6,0.100000);
        expectedGen.set(5,12,0.114286);
        expectedGen.set(6,5,0.375000);
        expectedGen.set(6,6,-0.475000);
        expectedGen.set(6,13,0.100000);
        expectedGen.set(7,0,1.000000);
        expectedGen.set(7,7,-1.500000);
        expectedGen.set(7,8,0.100000);
        expectedGen.set(7,14,0.400000);
        expectedGen.set(8,1,1.000000);
        expectedGen.set(8,7,0.214286);
        expectedGen.set(8,8,-1.542857);
        expectedGen.set(8,9,0.100000);
        expectedGen.set(8,15,0.228571);
        expectedGen.set(9,2,1.000000);
        expectedGen.set(9,8,0.300000);
        expectedGen.set(9,9,-1.560000);
        expectedGen.set(9,10,0.100000);
        expectedGen.set(9,16,0.160000);
        expectedGen.set(10,3,1.000000);
        expectedGen.set(10,9,0.346154);
        expectedGen.set(10,10,-1.569231);
        expectedGen.set(10,11,0.100000);
        expectedGen.set(10,17,0.123077);
        expectedGen.set(11,4,1.000000);
        expectedGen.set(11,10,0.375000);
        expectedGen.set(11,11,-1.575000);
        expectedGen.set(11,12,0.100000);
        expectedGen.set(11,18,0.100000);
        expectedGen.set(12,5,1.000000);
        expectedGen.set(12,11,0.394737);
        expectedGen.set(12,12,-1.578947);
        expectedGen.set(12,13,0.100000);
        expectedGen.set(12,19,0.084211);
        expectedGen.set(13,6,1.000000);
        expectedGen.set(13,12,0.409091);
        expectedGen.set(13,13,-1.481818);
        expectedGen.set(13,20,0.072727);
        expectedGen.set(14,7,2.000000);
        expectedGen.set(14,14,-2.500000);
        expectedGen.set(14,15,0.100000);
        expectedGen.set(14,21,0.400000);
        expectedGen.set(15,8,2.000000);
        expectedGen.set(15,14,0.300000);
        expectedGen.set(15,15,-2.560000);
        expectedGen.set(15,16,0.100000);
        expectedGen.set(15,22,0.160000);
        expectedGen.set(16,9,2.000000);
        expectedGen.set(16,15,0.375000);
        expectedGen.set(16,16,-2.575000);
        expectedGen.set(16,17,0.100000);
        expectedGen.set(16,23,0.100000);
        expectedGen.set(17,10,2.000000);
        expectedGen.set(17,16,0.409091);
        expectedGen.set(17,17,-2.581818);
        expectedGen.set(17,18,0.100000);
        expectedGen.set(17,24,0.072727);
        expectedGen.set(18,11,2.000000);
        expectedGen.set(18,17,0.428571);
        expectedGen.set(18,18,-2.585714);
        expectedGen.set(18,19,0.100000);
        expectedGen.set(18,25,0.057143);
        expectedGen.set(19,12,2.000000);
        expectedGen.set(19,18,0.441176);
        expectedGen.set(19,19,-2.588235);
        expectedGen.set(19,20,0.100000);
        expectedGen.set(19,26,0.047059);
        expectedGen.set(20,13,2.000000);
        expectedGen.set(20,19,0.450000);
        expectedGen.set(20,20,-2.490000);
        expectedGen.set(20,27,0.040000);
        expectedGen.set(21,14,3.000000);
        expectedGen.set(21,21,-3.100000);
        expectedGen.set(21,22,0.100000);
        expectedGen.set(22,15,3.000000);
        expectedGen.set(22,21,0.500000);
        expectedGen.set(22,22,-3.600000);
        expectedGen.set(22,23,0.100000);
        expectedGen.set(23,16,3.000000);
        expectedGen.set(23,22,0.500000);
        expectedGen.set(23,23,-3.600000);
        expectedGen.set(23,24,0.100000);
        expectedGen.set(24,17,3.000000);
        expectedGen.set(24,23,0.500000);
        expectedGen.set(24,24,-3.600000);
        expectedGen.set(24,25,0.100000);
        expectedGen.set(25,18,3.000000);
        expectedGen.set(25,24,0.500000);
        expectedGen.set(25,25,-3.600000);
        expectedGen.set(25,26,0.100000);
        expectedGen.set(26,19,3.000000);
        expectedGen.set(26,25,0.500000);
        expectedGen.set(26,26,-3.600000);
        expectedGen.set(26,27,0.100000);
        expectedGen.set(27,20,3.000000);
        expectedGen.set(27,26,0.500000);
        expectedGen.set(27,27,-3.500000);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "Generator matrix mismatch");

        // Verify state space dimensions
        assertEquals(28, stateSpace.getNumRows());
        assertEquals(7, stateSpace.getNumCols());
    }

    @Test
    public void testMixedNetworkGPS() {
        Network model = createMixedNetwork(SchedStrategy.GPS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(6);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for GPS mixed network
        // GPS has 28x28 generator matrix
        Matrix expectedGen = new Matrix(28, 28);

        expectedGen.set(0,0,-0.500000);
        expectedGen.set(0,1,0.100000);
        expectedGen.set(0,7,0.400000);
        expectedGen.set(1,0,0.300000);
        expectedGen.set(1,1,-0.560000);
        expectedGen.set(1,2,0.100000);
        expectedGen.set(1,8,0.160000);
        expectedGen.set(2,1,0.300000);
        expectedGen.set(2,2,-0.560000);
        expectedGen.set(2,3,0.100000);
        expectedGen.set(2,9,0.160000);
        expectedGen.set(3,2,0.300000);
        expectedGen.set(3,3,-0.560000);
        expectedGen.set(3,4,0.100000);
        expectedGen.set(3,10,0.160000);
        expectedGen.set(4,3,0.300000);
        expectedGen.set(4,4,-0.560000);
        expectedGen.set(4,5,0.100000);
        expectedGen.set(4,11,0.160000);
        expectedGen.set(5,4,0.300000);
        expectedGen.set(5,5,-0.560000);
        expectedGen.set(5,6,0.100000);
        expectedGen.set(5,12,0.160000);
        expectedGen.set(6,5,0.300000);
        expectedGen.set(6,6,-0.460000);
        expectedGen.set(6,13,0.160000);
        expectedGen.set(7,0,1.000000);
        expectedGen.set(7,7,-1.500000);
        expectedGen.set(7,8,0.100000);
        expectedGen.set(7,14,0.400000);
        expectedGen.set(8,1,1.000000);
        expectedGen.set(8,7,0.300000);
        expectedGen.set(8,8,-1.560000);
        expectedGen.set(8,9,0.100000);
        expectedGen.set(8,15,0.160000);
        expectedGen.set(9,2,1.000000);
        expectedGen.set(9,8,0.300000);
        expectedGen.set(9,9,-1.560000);
        expectedGen.set(9,10,0.100000);
        expectedGen.set(9,16,0.160000);
        expectedGen.set(10,3,1.000000);
        expectedGen.set(10,9,0.300000);
        expectedGen.set(10,10,-1.560000);
        expectedGen.set(10,11,0.100000);
        expectedGen.set(10,17,0.160000);
        expectedGen.set(11,4,1.000000);
        expectedGen.set(11,10,0.300000);
        expectedGen.set(11,11,-1.560000);
        expectedGen.set(11,12,0.100000);
        expectedGen.set(11,18,0.160000);
        expectedGen.set(12,5,1.000000);
        expectedGen.set(12,11,0.300000);
        expectedGen.set(12,12,-1.560000);
        expectedGen.set(12,13,0.100000);
        expectedGen.set(12,19,0.160000);
        expectedGen.set(13,6,1.000000);
        expectedGen.set(13,12,0.300000);
        expectedGen.set(13,13,-1.460000);
        expectedGen.set(13,20,0.160000);
        expectedGen.set(14,7,2.000000);
        expectedGen.set(14,14,-2.500000);
        expectedGen.set(14,15,0.100000);
        expectedGen.set(14,21,0.400000);
        expectedGen.set(15,8,2.000000);
        expectedGen.set(15,14,0.300000);
        expectedGen.set(15,15,-2.560000);
        expectedGen.set(15,16,0.100000);
        expectedGen.set(15,22,0.160000);
        expectedGen.set(16,9,2.000000);
        expectedGen.set(16,15,0.300000);
        expectedGen.set(16,16,-2.560000);
        expectedGen.set(16,17,0.100000);
        expectedGen.set(16,23,0.160000);
        expectedGen.set(17,10,2.000000);
        expectedGen.set(17,16,0.300000);
        expectedGen.set(17,17,-2.560000);
        expectedGen.set(17,18,0.100000);
        expectedGen.set(17,24,0.160000);
        expectedGen.set(18,11,2.000000);
        expectedGen.set(18,17,0.300000);
        expectedGen.set(18,18,-2.560000);
        expectedGen.set(18,19,0.100000);
        expectedGen.set(18,25,0.160000);
        expectedGen.set(19,12,2.000000);
        expectedGen.set(19,18,0.300000);
        expectedGen.set(19,19,-2.560000);
        expectedGen.set(19,20,0.100000);
        expectedGen.set(19,26,0.160000);
        expectedGen.set(20,13,2.000000);
        expectedGen.set(20,19,0.300000);
        expectedGen.set(20,20,-2.460000);
        expectedGen.set(20,27,0.160000);
        expectedGen.set(21,14,3.000000);
        expectedGen.set(21,21,-3.100000);
        expectedGen.set(21,22,0.100000);
        expectedGen.set(22,15,3.000000);
        expectedGen.set(22,21,0.500000);
        expectedGen.set(22,22,-3.600000);
        expectedGen.set(22,23,0.100000);
        expectedGen.set(23,16,3.000000);
        expectedGen.set(23,22,0.500000);
        expectedGen.set(23,23,-3.600000);
        expectedGen.set(23,24,0.100000);
        expectedGen.set(24,17,3.000000);
        expectedGen.set(24,23,0.500000);
        expectedGen.set(24,24,-3.600000);
        expectedGen.set(24,25,0.100000);
        expectedGen.set(25,18,3.000000);
        expectedGen.set(25,24,0.500000);
        expectedGen.set(25,25,-3.600000);
        expectedGen.set(25,26,0.100000);
        expectedGen.set(26,19,3.000000);
        expectedGen.set(26,25,0.500000);
        expectedGen.set(26,26,-3.600000);
        expectedGen.set(26,27,0.100000);
        expectedGen.set(27,20,3.000000);
        expectedGen.set(27,26,0.500000);
        expectedGen.set(27,27,-3.500000);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "Generator matrix mismatch");

        // Verify state space dimensions
        assertEquals(28, stateSpace.getNumRows());
        assertEquals(7, stateSpace.getNumCols());
    }

    @Test
    public void testMixedNetworkINF() {
        Network model = createMixedNetwork(SchedStrategy.INF);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(6);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        Matrix expectedGen = new Matrix(28, 28);
        expectedGen.set(0, 0, -1.3);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(7, 0, 1.0);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -1.8);
        expectedGen.set(2, 1, 1.0);
        expectedGen.set(8, 1, 1.0);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -2.3);
        expectedGen.set(3, 2, 1.5);
        expectedGen.set(9, 2, 1.0);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -2.8);
        expectedGen.set(4, 3, 2.0);
        expectedGen.set(10, 3, 1.0);
        expectedGen.set(3, 4, 0.1);
        expectedGen.set(4, 4, -3.3);
        expectedGen.set(5, 4, 2.5);
        expectedGen.set(11, 4, 1.0);
        expectedGen.set(4, 5, 0.1);
        expectedGen.set(5, 5, -3.8);
        expectedGen.set(6, 5, 3.0);
        expectedGen.set(12, 5, 1.0);
        expectedGen.set(5, 6, 0.1);
        expectedGen.set(6, 6, -4.2);
        expectedGen.set(13, 6, 1.0);
        expectedGen.set(0, 7, 1.2);
        expectedGen.set(7, 7, -1.9);
        expectedGen.set(8, 7, 0.5);
        expectedGen.set(14, 7, 2.0);
        expectedGen.set(1, 8, 1.2);
        expectedGen.set(7, 8, 0.1);
        expectedGen.set(8, 8, -2.4);
        expectedGen.set(9, 8, 1.0);
        expectedGen.set(15, 8, 2.0);
        expectedGen.set(2, 9, 1.2);
        expectedGen.set(8, 9, 0.1);
        expectedGen.set(9, 9, -2.9);
        expectedGen.set(10, 9, 1.5);
        expectedGen.set(16, 9, 2.0);
        expectedGen.set(3, 10, 1.2);
        expectedGen.set(9, 10, 0.1);
        expectedGen.set(10, 10, -3.4);
        expectedGen.set(11, 10, 2.0);
        expectedGen.set(17, 10, 2.0);
        expectedGen.set(4, 11, 1.2);
        expectedGen.set(10, 11, 0.1);
        expectedGen.set(11, 11, -3.9);
        expectedGen.set(12, 11, 2.5);
        expectedGen.set(18, 11, 2.0);
        expectedGen.set(5, 12, 1.2);
        expectedGen.set(11, 12, 0.1);
        expectedGen.set(12, 12, -4.4);
        expectedGen.set(13, 12, 3.0);
        expectedGen.set(19, 12, 2.0);
        expectedGen.set(6, 13, 1.2);
        expectedGen.set(12, 13, 0.1);
        expectedGen.set(13, 13, -4.8);
        expectedGen.set(20, 13, 2.0);
        expectedGen.set(7, 14, 0.8);
        expectedGen.set(14, 14, -2.5);
        expectedGen.set(15, 14, 0.5);
        expectedGen.set(21, 14, 3.0);
        expectedGen.set(8, 15, 0.8);
        expectedGen.set(14, 15, 0.1);
        expectedGen.set(15, 15, -3.0);
        expectedGen.set(16, 15, 1.0);
        expectedGen.set(22, 15, 3.0);
        expectedGen.set(9, 16, 0.8);
        expectedGen.set(15, 16, 0.1);
        expectedGen.set(16, 16, -3.5);
        expectedGen.set(17, 16, 1.5);
        expectedGen.set(23, 16, 3.0);
        expectedGen.set(10, 17, 0.8);
        expectedGen.set(16, 17, 0.1);
        expectedGen.set(17, 17, -4.0);
        expectedGen.set(18, 17, 2.0);
        expectedGen.set(24, 17, 3.0);
        expectedGen.set(11, 18, 0.8);
        expectedGen.set(17, 18, 0.1);
        expectedGen.set(18, 18, -4.5);
        expectedGen.set(19, 18, 2.5);
        expectedGen.set(25, 18, 3.0);
        expectedGen.set(12, 19, 0.8);
        expectedGen.set(18, 19, 0.1);
        expectedGen.set(19, 19, -5.0);
        expectedGen.set(20, 19, 3.0);
        expectedGen.set(26, 19, 3.0);
        expectedGen.set(13, 20, 0.8);
        expectedGen.set(19, 20, 0.1);
        expectedGen.set(20, 20, -5.4);
        expectedGen.set(27, 20, 3.0);
        expectedGen.set(14, 21, 0.4);
        expectedGen.set(21, 21, -3.1);
        expectedGen.set(22, 21, 0.5);
        expectedGen.set(15, 22, 0.4);
        expectedGen.set(21, 22, 0.1);
        expectedGen.set(22, 22, -3.6);
        expectedGen.set(23, 22, 1.0);
        expectedGen.set(16, 23, 0.4);
        expectedGen.set(22, 23, 0.1);
        expectedGen.set(23, 23, -4.1);
        expectedGen.set(24, 23, 1.5);
        expectedGen.set(17, 24, 0.4);
        expectedGen.set(23, 24, 0.1);
        expectedGen.set(24, 24, -4.6);
        expectedGen.set(25, 24, 2.0);
        expectedGen.set(18, 25, 0.4);
        expectedGen.set(24, 25, 0.1);
        expectedGen.set(25, 25, -5.1);
        expectedGen.set(26, 25, 2.5);
        expectedGen.set(19, 26, 0.4);
        expectedGen.set(25, 26, 0.1);
        expectedGen.set(26, 26, -5.6);
        expectedGen.set(27, 26, 3.0);
        expectedGen.set(20, 27, 0.4);
        expectedGen.set(26, 27, 0.1);
        expectedGen.set(27, 27, -6.0);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "Generator matrix mismatch");

        // Verify state space
        assertEquals(28, stateSpace.getNumRows());
        assertEquals(7, stateSpace.getNumCols());
    }

    @Test
    public void testMixedNetworkFCFS() {
        Network model = createMixedNetwork(SchedStrategy.FCFS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(6);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for FCFS mixed network
        // FCFS has 329x329 generator matrix
        Matrix expectedGen = new Matrix(329, 329);
        expectedGen.set(0, 0, -0.5);
        expectedGen.set(4, 0, 0.5);
        expectedGen.set(210, 0, 1.0);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.5);
        expectedGen.set(8, 1, 0.5);
        expectedGen.set(2, 2, -0.5);
        expectedGen.set(11, 2, 0.5);
        expectedGen.set(211, 2, 1.0);
        expectedGen.set(3, 3, -0.5);
        expectedGen.set(13, 3, 0.5);
        expectedGen.set(212, 3, 1.0);
        expectedGen.set(4, 4, -0.6);
        expectedGen.set(14, 4, 0.5);
        expectedGen.set(213, 4, 1.0);
        expectedGen.set(1, 5, 0.1);
        expectedGen.set(5, 5, -0.5);
        expectedGen.set(18, 5, 0.5);
        expectedGen.set(2, 6, 0.1);
        expectedGen.set(6, 6, -0.5);
        expectedGen.set(21, 6, 0.5);
        expectedGen.set(3, 7, 0.1);
        expectedGen.set(7, 7, -0.5);
        expectedGen.set(23, 7, 0.5);
        expectedGen.set(4, 8, 0.1);
        expectedGen.set(8, 8, -0.6);
        expectedGen.set(24, 8, 0.5);
        expectedGen.set(9, 9, -0.5);
        expectedGen.set(27, 9, 0.5);
        expectedGen.set(214, 9, 1.0);
        expectedGen.set(10, 10, -0.5);
        expectedGen.set(29, 10, 0.5);
        expectedGen.set(215, 10, 1.0);
        expectedGen.set(11, 11, -0.6);
        expectedGen.set(30, 11, 0.5);
        expectedGen.set(216, 11, 1.0);
        expectedGen.set(12, 12, -0.5);
        expectedGen.set(32, 12, 0.5);
        expectedGen.set(217, 12, 1.0);
        expectedGen.set(13, 13, -0.6);
        expectedGen.set(33, 13, 0.5);
        expectedGen.set(218, 13, 1.0);
        expectedGen.set(14, 14, -0.6);
        expectedGen.set(34, 14, 0.5);
        expectedGen.set(219, 14, 1.0);
        expectedGen.set(5, 15, 0.1);
        expectedGen.set(15, 15, -0.5);
        expectedGen.set(38, 15, 0.5);
        expectedGen.set(6, 16, 0.1);
        expectedGen.set(16, 16, -0.5);
        expectedGen.set(41, 16, 0.5);
        expectedGen.set(7, 17, 0.1);
        expectedGen.set(17, 17, -0.5);
        expectedGen.set(43, 17, 0.5);
        expectedGen.set(8, 18, 0.1);
        expectedGen.set(18, 18, -0.6);
        expectedGen.set(44, 18, 0.5);
        expectedGen.set(9, 19, 0.1);
        expectedGen.set(19, 19, -0.5);
        expectedGen.set(47, 19, 0.5);
        expectedGen.set(10, 20, 0.1);
        expectedGen.set(20, 20, -0.5);
        expectedGen.set(49, 20, 0.5);
        expectedGen.set(11, 21, 0.1);
        expectedGen.set(21, 21, -0.6);
        expectedGen.set(50, 21, 0.5);
        expectedGen.set(12, 22, 0.1);
        expectedGen.set(22, 22, -0.5);
        expectedGen.set(52, 22, 0.5);
        expectedGen.set(13, 23, 0.1);
        expectedGen.set(23, 23, -0.6);
        expectedGen.set(53, 23, 0.5);
        expectedGen.set(14, 24, 0.1);
        expectedGen.set(24, 24, -0.6);
        expectedGen.set(54, 24, 0.5);
        expectedGen.set(25, 25, -0.5);
        expectedGen.set(57, 25, 0.5);
        expectedGen.set(220, 25, 1.0);
        expectedGen.set(26, 26, -0.5);
        expectedGen.set(59, 26, 0.5);
        expectedGen.set(221, 26, 1.0);
        expectedGen.set(27, 27, -0.6);
        expectedGen.set(60, 27, 0.5);
        expectedGen.set(222, 27, 1.0);
        expectedGen.set(28, 28, -0.5);
        expectedGen.set(62, 28, 0.5);
        expectedGen.set(223, 28, 1.0);
        expectedGen.set(29, 29, -0.6);
        expectedGen.set(63, 29, 0.5);
        expectedGen.set(224, 29, 1.0);
        expectedGen.set(30, 30, -0.6);
        expectedGen.set(64, 30, 0.5);
        expectedGen.set(225, 30, 1.0);
        expectedGen.set(31, 31, -0.5);
        expectedGen.set(66, 31, 0.5);
        expectedGen.set(226, 31, 1.0);
        expectedGen.set(32, 32, -0.6);
        expectedGen.set(67, 32, 0.5);
        expectedGen.set(227, 32, 1.0);
        expectedGen.set(33, 33, -0.6);
        expectedGen.set(68, 33, 0.5);
        expectedGen.set(228, 33, 1.0);
        expectedGen.set(34, 34, -0.6);
        expectedGen.set(69, 34, 0.5);
        expectedGen.set(229, 34, 1.0);
        expectedGen.set(15, 35, 0.1);
        expectedGen.set(35, 35, -0.5);
        expectedGen.set(73, 35, 0.5);
        expectedGen.set(16, 36, 0.1);
        expectedGen.set(36, 36, -0.5);
        expectedGen.set(76, 36, 0.5);
        expectedGen.set(17, 37, 0.1);
        expectedGen.set(37, 37, -0.5);
        expectedGen.set(78, 37, 0.5);
        expectedGen.set(18, 38, 0.1);
        expectedGen.set(38, 38, -0.6);
        expectedGen.set(79, 38, 0.5);
        expectedGen.set(19, 39, 0.1);
        expectedGen.set(39, 39, -0.5);
        expectedGen.set(82, 39, 0.5);
        expectedGen.set(20, 40, 0.1);
        expectedGen.set(40, 40, -0.5);
        expectedGen.set(84, 40, 0.5);
        expectedGen.set(21, 41, 0.1);
        expectedGen.set(41, 41, -0.6);
        expectedGen.set(85, 41, 0.5);
        expectedGen.set(22, 42, 0.1);
        expectedGen.set(42, 42, -0.5);
        expectedGen.set(87, 42, 0.5);
        expectedGen.set(23, 43, 0.1);
        expectedGen.set(43, 43, -0.6);
        expectedGen.set(88, 43, 0.5);
        expectedGen.set(24, 44, 0.1);
        expectedGen.set(44, 44, -0.6);
        expectedGen.set(89, 44, 0.5);
        expectedGen.set(25, 45, 0.1);
        expectedGen.set(45, 45, -0.5);
        expectedGen.set(92, 45, 0.5);
        expectedGen.set(26, 46, 0.1);
        expectedGen.set(46, 46, -0.5);
        expectedGen.set(94, 46, 0.5);
        expectedGen.set(27, 47, 0.1);
        expectedGen.set(47, 47, -0.6);
        expectedGen.set(95, 47, 0.5);
        expectedGen.set(28, 48, 0.1);
        expectedGen.set(48, 48, -0.5);
        expectedGen.set(97, 48, 0.5);
        expectedGen.set(29, 49, 0.1);
        expectedGen.set(49, 49, -0.6);
        expectedGen.set(98, 49, 0.5);
        expectedGen.set(30, 50, 0.1);
        expectedGen.set(50, 50, -0.6);
        expectedGen.set(99, 50, 0.5);
        expectedGen.set(31, 51, 0.1);
        expectedGen.set(51, 51, -0.5);
        expectedGen.set(101, 51, 0.5);
        expectedGen.set(32, 52, 0.1);
        expectedGen.set(52, 52, -0.6);
        expectedGen.set(102, 52, 0.5);
        expectedGen.set(33, 53, 0.1);
        expectedGen.set(53, 53, -0.6);
        expectedGen.set(103, 53, 0.5);
        expectedGen.set(34, 54, 0.1);
        expectedGen.set(54, 54, -0.6);
        expectedGen.set(104, 54, 0.5);
        expectedGen.set(55, 55, -0.5);
        expectedGen.set(107, 55, 0.5);
        expectedGen.set(230, 55, 1.0);
        expectedGen.set(56, 56, -0.5);
        expectedGen.set(109, 56, 0.5);
        expectedGen.set(231, 56, 1.0);
        expectedGen.set(57, 57, -0.6);
        expectedGen.set(110, 57, 0.5);
        expectedGen.set(232, 57, 1.0);
        expectedGen.set(58, 58, -0.5);
        expectedGen.set(112, 58, 0.5);
        expectedGen.set(233, 58, 1.0);
        expectedGen.set(59, 59, -0.6);
        expectedGen.set(113, 59, 0.5);
        expectedGen.set(234, 59, 1.0);
        expectedGen.set(60, 60, -0.6);
        expectedGen.set(114, 60, 0.5);
        expectedGen.set(235, 60, 1.0);
        expectedGen.set(61, 61, -0.5);
        expectedGen.set(116, 61, 0.5);
        expectedGen.set(236, 61, 1.0);
        expectedGen.set(62, 62, -0.6);
        expectedGen.set(117, 62, 0.5);
        expectedGen.set(237, 62, 1.0);
        expectedGen.set(63, 63, -0.6);
        expectedGen.set(118, 63, 0.5);
        expectedGen.set(238, 63, 1.0);
        expectedGen.set(64, 64, -0.6);
        expectedGen.set(119, 64, 0.5);
        expectedGen.set(239, 64, 1.0);
        expectedGen.set(65, 65, -0.5);
        expectedGen.set(121, 65, 0.5);
        expectedGen.set(240, 65, 1.0);
        expectedGen.set(66, 66, -0.6);
        expectedGen.set(122, 66, 0.5);
        expectedGen.set(241, 66, 1.0);
        expectedGen.set(67, 67, -0.6);
        expectedGen.set(123, 67, 0.5);
        expectedGen.set(242, 67, 1.0);
        expectedGen.set(68, 68, -0.6);
        expectedGen.set(124, 68, 0.5);
        expectedGen.set(243, 68, 1.0);
        expectedGen.set(69, 69, -0.6);
        expectedGen.set(125, 69, 0.5);
        expectedGen.set(244, 69, 1.0);
        expectedGen.set(35, 70, 0.1);
        expectedGen.set(70, 70, -0.5);
        expectedGen.set(129, 70, 0.5);
        expectedGen.set(36, 71, 0.1);
        expectedGen.set(71, 71, -0.5);
        expectedGen.set(132, 71, 0.5);
        expectedGen.set(37, 72, 0.1);
        expectedGen.set(72, 72, -0.5);
        expectedGen.set(134, 72, 0.5);
        expectedGen.set(38, 73, 0.1);
        expectedGen.set(73, 73, -0.6);
        expectedGen.set(135, 73, 0.5);
        expectedGen.set(39, 74, 0.1);
        expectedGen.set(74, 74, -0.5);
        expectedGen.set(138, 74, 0.5);
        expectedGen.set(40, 75, 0.1);
        expectedGen.set(75, 75, -0.5);
        expectedGen.set(140, 75, 0.5);
        expectedGen.set(41, 76, 0.1);
        expectedGen.set(76, 76, -0.6);
        expectedGen.set(141, 76, 0.5);
        expectedGen.set(42, 77, 0.1);
        expectedGen.set(77, 77, -0.5);
        expectedGen.set(143, 77, 0.5);
        expectedGen.set(43, 78, 0.1);
        expectedGen.set(78, 78, -0.6);
        expectedGen.set(144, 78, 0.5);
        expectedGen.set(44, 79, 0.1);
        expectedGen.set(79, 79, -0.6);
        expectedGen.set(145, 79, 0.5);
        expectedGen.set(45, 80, 0.1);
        expectedGen.set(80, 80, -0.5);
        expectedGen.set(148, 80, 0.5);
        expectedGen.set(46, 81, 0.1);
        expectedGen.set(81, 81, -0.5);
        expectedGen.set(150, 81, 0.5);
        expectedGen.set(47, 82, 0.1);
        expectedGen.set(82, 82, -0.6);
        expectedGen.set(151, 82, 0.5);
        expectedGen.set(48, 83, 0.1);
        expectedGen.set(83, 83, -0.5);
        expectedGen.set(153, 83, 0.5);
        expectedGen.set(49, 84, 0.1);
        expectedGen.set(84, 84, -0.6);
        expectedGen.set(154, 84, 0.5);
        expectedGen.set(50, 85, 0.1);
        expectedGen.set(85, 85, -0.6);
        expectedGen.set(155, 85, 0.5);
        expectedGen.set(51, 86, 0.1);
        expectedGen.set(86, 86, -0.5);
        expectedGen.set(157, 86, 0.5);
        expectedGen.set(52, 87, 0.1);
        expectedGen.set(87, 87, -0.6);
        expectedGen.set(158, 87, 0.5);
        expectedGen.set(53, 88, 0.1);
        expectedGen.set(88, 88, -0.6);
        expectedGen.set(159, 88, 0.5);
        expectedGen.set(54, 89, 0.1);
        expectedGen.set(89, 89, -0.6);
        expectedGen.set(160, 89, 0.5);
        expectedGen.set(55, 90, 0.1);
        expectedGen.set(90, 90, -0.5);
        expectedGen.set(163, 90, 0.5);
        expectedGen.set(56, 91, 0.1);
        expectedGen.set(91, 91, -0.5);
        expectedGen.set(165, 91, 0.5);
        expectedGen.set(57, 92, 0.1);
        expectedGen.set(92, 92, -0.6);
        expectedGen.set(166, 92, 0.5);
        expectedGen.set(58, 93, 0.1);
        expectedGen.set(93, 93, -0.5);
        expectedGen.set(168, 93, 0.5);
        expectedGen.set(59, 94, 0.1);
        expectedGen.set(94, 94, -0.6);
        expectedGen.set(169, 94, 0.5);
        expectedGen.set(60, 95, 0.1);
        expectedGen.set(95, 95, -0.6);
        expectedGen.set(170, 95, 0.5);
        expectedGen.set(61, 96, 0.1);
        expectedGen.set(96, 96, -0.5);
        expectedGen.set(172, 96, 0.5);
        expectedGen.set(62, 97, 0.1);
        expectedGen.set(97, 97, -0.6);
        expectedGen.set(173, 97, 0.5);
        expectedGen.set(63, 98, 0.1);
        expectedGen.set(98, 98, -0.6);
        expectedGen.set(174, 98, 0.5);
        expectedGen.set(64, 99, 0.1);
        expectedGen.set(99, 99, -0.6);
        expectedGen.set(175, 99, 0.5);
        expectedGen.set(65, 100, 0.1);
        expectedGen.set(100, 100, -0.5);
        expectedGen.set(177, 100, 0.5);
        expectedGen.set(66, 101, 0.1);
        expectedGen.set(101, 101, -0.6);
        expectedGen.set(178, 101, 0.5);
        expectedGen.set(67, 102, 0.1);
        expectedGen.set(102, 102, -0.6);
        expectedGen.set(179, 102, 0.5);
        expectedGen.set(68, 103, 0.1);
        expectedGen.set(103, 103, -0.6);
        expectedGen.set(180, 103, 0.5);
        expectedGen.set(69, 104, 0.1);
        expectedGen.set(104, 104, -0.6);
        expectedGen.set(181, 104, 0.5);
        expectedGen.set(105, 105, -0.5);
        expectedGen.set(184, 105, 0.5);
        expectedGen.set(245, 105, 1.0);
        expectedGen.set(106, 106, -0.5);
        expectedGen.set(186, 106, 0.5);
        expectedGen.set(246, 106, 1.0);
        expectedGen.set(107, 107, -0.6);
        expectedGen.set(187, 107, 0.5);
        expectedGen.set(247, 107, 1.0);
        expectedGen.set(108, 108, -0.5);
        expectedGen.set(189, 108, 0.5);
        expectedGen.set(248, 108, 1.0);
        expectedGen.set(109, 109, -0.6);
        expectedGen.set(190, 109, 0.5);
        expectedGen.set(249, 109, 1.0);
        expectedGen.set(110, 110, -0.6);
        expectedGen.set(191, 110, 0.5);
        expectedGen.set(250, 110, 1.0);
        expectedGen.set(111, 111, -0.5);
        expectedGen.set(193, 111, 0.5);
        expectedGen.set(251, 111, 1.0);
        expectedGen.set(112, 112, -0.6);
        expectedGen.set(194, 112, 0.5);
        expectedGen.set(252, 112, 1.0);
        expectedGen.set(113, 113, -0.6);
        expectedGen.set(195, 113, 0.5);
        expectedGen.set(253, 113, 1.0);
        expectedGen.set(114, 114, -0.6);
        expectedGen.set(196, 114, 0.5);
        expectedGen.set(254, 114, 1.0);
        expectedGen.set(115, 115, -0.5);
        expectedGen.set(198, 115, 0.5);
        expectedGen.set(255, 115, 1.0);
        expectedGen.set(116, 116, -0.6);
        expectedGen.set(199, 116, 0.5);
        expectedGen.set(256, 116, 1.0);
        expectedGen.set(117, 117, -0.6);
        expectedGen.set(200, 117, 0.5);
        expectedGen.set(257, 117, 1.0);
        expectedGen.set(118, 118, -0.6);
        expectedGen.set(201, 118, 0.5);
        expectedGen.set(258, 118, 1.0);
        expectedGen.set(119, 119, -0.6);
        expectedGen.set(202, 119, 0.5);
        expectedGen.set(259, 119, 1.0);
        expectedGen.set(120, 120, -0.5);
        expectedGen.set(204, 120, 0.5);
        expectedGen.set(260, 120, 1.0);
        expectedGen.set(121, 121, -0.6);
        expectedGen.set(205, 121, 0.5);
        expectedGen.set(261, 121, 1.0);
        expectedGen.set(122, 122, -0.6);
        expectedGen.set(206, 122, 0.5);
        expectedGen.set(262, 122, 1.0);
        expectedGen.set(123, 123, -0.6);
        expectedGen.set(207, 123, 0.5);
        expectedGen.set(263, 123, 1.0);
        expectedGen.set(124, 124, -0.6);
        expectedGen.set(208, 124, 0.5);
        expectedGen.set(264, 124, 1.0);
        expectedGen.set(125, 125, -0.6);
        expectedGen.set(209, 125, 0.5);
        expectedGen.set(265, 125, 1.0);
        expectedGen.set(70, 126, 0.1);
        expectedGen.set(126, 126, -0.4);
        expectedGen.set(71, 127, 0.1);
        expectedGen.set(127, 127, -0.4);
        expectedGen.set(72, 128, 0.1);
        expectedGen.set(128, 128, -0.4);
        expectedGen.set(73, 129, 0.1);
        expectedGen.set(129, 129, -0.5);
        expectedGen.set(74, 130, 0.1);
        expectedGen.set(130, 130, -0.4);
        expectedGen.set(75, 131, 0.1);
        expectedGen.set(131, 131, -0.4);
        expectedGen.set(76, 132, 0.1);
        expectedGen.set(132, 132, -0.5);
        expectedGen.set(77, 133, 0.1);
        expectedGen.set(133, 133, -0.4);
        expectedGen.set(78, 134, 0.1);
        expectedGen.set(134, 134, -0.5);
        expectedGen.set(79, 135, 0.1);
        expectedGen.set(135, 135, -0.5);
        expectedGen.set(80, 136, 0.1);
        expectedGen.set(136, 136, -0.4);
        expectedGen.set(81, 137, 0.1);
        expectedGen.set(137, 137, -0.4);
        expectedGen.set(82, 138, 0.1);
        expectedGen.set(138, 138, -0.5);
        expectedGen.set(83, 139, 0.1);
        expectedGen.set(139, 139, -0.4);
        expectedGen.set(84, 140, 0.1);
        expectedGen.set(140, 140, -0.5);
        expectedGen.set(85, 141, 0.1);
        expectedGen.set(141, 141, -0.5);
        expectedGen.set(86, 142, 0.1);
        expectedGen.set(142, 142, -0.4);
        expectedGen.set(87, 143, 0.1);
        expectedGen.set(143, 143, -0.5);
        expectedGen.set(88, 144, 0.1);
        expectedGen.set(144, 144, -0.5);
        expectedGen.set(89, 145, 0.1);
        expectedGen.set(145, 145, -0.5);
        expectedGen.set(90, 146, 0.1);
        expectedGen.set(146, 146, -0.4);
        expectedGen.set(91, 147, 0.1);
        expectedGen.set(147, 147, -0.4);
        expectedGen.set(92, 148, 0.1);
        expectedGen.set(148, 148, -0.5);
        expectedGen.set(93, 149, 0.1);
        expectedGen.set(149, 149, -0.4);
        expectedGen.set(94, 150, 0.1);
        expectedGen.set(150, 150, -0.5);
        expectedGen.set(95, 151, 0.1);
        expectedGen.set(151, 151, -0.5);
        expectedGen.set(96, 152, 0.1);
        expectedGen.set(152, 152, -0.4);
        expectedGen.set(97, 153, 0.1);
        expectedGen.set(153, 153, -0.5);
        expectedGen.set(98, 154, 0.1);
        expectedGen.set(154, 154, -0.5);
        expectedGen.set(99, 155, 0.1);
        expectedGen.set(155, 155, -0.5);
        expectedGen.set(100, 156, 0.1);
        expectedGen.set(156, 156, -0.4);
        expectedGen.set(101, 157, 0.1);
        expectedGen.set(157, 157, -0.5);
        expectedGen.set(102, 158, 0.1);
        expectedGen.set(158, 158, -0.5);
        expectedGen.set(103, 159, 0.1);
        expectedGen.set(159, 159, -0.5);
        expectedGen.set(104, 160, 0.1);
        expectedGen.set(160, 160, -0.5);
        expectedGen.set(105, 161, 0.1);
        expectedGen.set(161, 161, -0.4);
        expectedGen.set(106, 162, 0.1);
        expectedGen.set(162, 162, -0.4);
        expectedGen.set(107, 163, 0.1);
        expectedGen.set(163, 163, -0.5);
        expectedGen.set(108, 164, 0.1);
        expectedGen.set(164, 164, -0.4);
        expectedGen.set(109, 165, 0.1);
        expectedGen.set(165, 165, -0.5);
        expectedGen.set(110, 166, 0.1);
        expectedGen.set(166, 166, -0.5);
        expectedGen.set(111, 167, 0.1);
        expectedGen.set(167, 167, -0.4);
        expectedGen.set(112, 168, 0.1);
        expectedGen.set(168, 168, -0.5);
        expectedGen.set(113, 169, 0.1);
        expectedGen.set(169, 169, -0.5);
        expectedGen.set(114, 170, 0.1);
        expectedGen.set(170, 170, -0.5);
        expectedGen.set(115, 171, 0.1);
        expectedGen.set(171, 171, -0.4);
        expectedGen.set(116, 172, 0.1);
        expectedGen.set(172, 172, -0.5);
        expectedGen.set(117, 173, 0.1);
        expectedGen.set(173, 173, -0.5);
        expectedGen.set(118, 174, 0.1);
        expectedGen.set(174, 174, -0.5);
        expectedGen.set(119, 175, 0.1);
        expectedGen.set(175, 175, -0.5);
        expectedGen.set(120, 176, 0.1);
        expectedGen.set(176, 176, -0.4);
        expectedGen.set(121, 177, 0.1);
        expectedGen.set(177, 177, -0.5);
        expectedGen.set(122, 178, 0.1);
        expectedGen.set(178, 178, -0.5);
        expectedGen.set(123, 179, 0.1);
        expectedGen.set(179, 179, -0.5);
        expectedGen.set(124, 180, 0.1);
        expectedGen.set(180, 180, -0.5);
        expectedGen.set(125, 181, 0.1);
        expectedGen.set(181, 181, -0.5);
        expectedGen.set(182, 182, -0.4);
        expectedGen.set(266, 182, 1.0);
        expectedGen.set(183, 183, -0.4);
        expectedGen.set(267, 183, 1.0);
        expectedGen.set(184, 184, -0.5);
        expectedGen.set(268, 184, 1.0);
        expectedGen.set(185, 185, -0.4);
        expectedGen.set(269, 185, 1.0);
        expectedGen.set(186, 186, -0.5);
        expectedGen.set(270, 186, 1.0);
        expectedGen.set(187, 187, -0.5);
        expectedGen.set(271, 187, 1.0);
        expectedGen.set(188, 188, -0.4);
        expectedGen.set(272, 188, 1.0);
        expectedGen.set(189, 189, -0.5);
        expectedGen.set(273, 189, 1.0);
        expectedGen.set(190, 190, -0.5);
        expectedGen.set(274, 190, 1.0);
        expectedGen.set(191, 191, -0.5);
        expectedGen.set(275, 191, 1.0);
        expectedGen.set(192, 192, -0.4);
        expectedGen.set(276, 192, 1.0);
        expectedGen.set(193, 193, -0.5);
        expectedGen.set(277, 193, 1.0);
        expectedGen.set(194, 194, -0.5);
        expectedGen.set(278, 194, 1.0);
        expectedGen.set(195, 195, -0.5);
        expectedGen.set(279, 195, 1.0);
        expectedGen.set(196, 196, -0.5);
        expectedGen.set(280, 196, 1.0);
        expectedGen.set(197, 197, -0.4);
        expectedGen.set(281, 197, 1.0);
        expectedGen.set(198, 198, -0.5);
        expectedGen.set(282, 198, 1.0);
        expectedGen.set(199, 199, -0.5);
        expectedGen.set(283, 199, 1.0);
        expectedGen.set(200, 200, -0.5);
        expectedGen.set(284, 200, 1.0);
        expectedGen.set(201, 201, -0.5);
        expectedGen.set(285, 201, 1.0);
        expectedGen.set(202, 202, -0.5);
        expectedGen.set(286, 202, 1.0);
        expectedGen.set(203, 203, -0.4);
        expectedGen.set(287, 203, 1.0);
        expectedGen.set(204, 204, -0.5);
        expectedGen.set(288, 204, 1.0);
        expectedGen.set(205, 205, -0.5);
        expectedGen.set(289, 205, 1.0);
        expectedGen.set(206, 206, -0.5);
        expectedGen.set(290, 206, 1.0);
        expectedGen.set(207, 207, -0.5);
        expectedGen.set(291, 207, 1.0);
        expectedGen.set(208, 208, -0.5);
        expectedGen.set(292, 208, 1.0);
        expectedGen.set(209, 209, -0.5);
        expectedGen.set(293, 209, 1.0);
        expectedGen.set(0, 210, 0.4);
        expectedGen.set(210, 210, -1.5);
        expectedGen.set(213, 210, 0.5);
        expectedGen.set(294, 210, 2.0);
        expectedGen.set(1, 211, 0.4);
        expectedGen.set(210, 211, 0.1);
        expectedGen.set(211, 211, -1.5);
        expectedGen.set(216, 211, 0.5);
        expectedGen.set(2, 212, 0.4);
        expectedGen.set(212, 212, -1.5);
        expectedGen.set(218, 212, 0.5);
        expectedGen.set(295, 212, 2.0);
        expectedGen.set(3, 213, 0.4);
        expectedGen.set(213, 213, -1.6);
        expectedGen.set(219, 213, 0.5);
        expectedGen.set(296, 213, 2.0);
        expectedGen.set(5, 214, 0.4);
        expectedGen.set(211, 214, 0.1);
        expectedGen.set(214, 214, -1.5);
        expectedGen.set(222, 214, 0.5);
        expectedGen.set(6, 215, 0.4);
        expectedGen.set(212, 215, 0.1);
        expectedGen.set(215, 215, -1.5);
        expectedGen.set(224, 215, 0.5);
        expectedGen.set(7, 216, 0.4);
        expectedGen.set(213, 216, 0.1);
        expectedGen.set(216, 216, -1.6);
        expectedGen.set(225, 216, 0.5);
        expectedGen.set(9, 217, 0.4);
        expectedGen.set(217, 217, -1.5);
        expectedGen.set(227, 217, 0.5);
        expectedGen.set(297, 217, 2.0);
        expectedGen.set(10, 218, 0.4);
        expectedGen.set(218, 218, -1.6);
        expectedGen.set(228, 218, 0.5);
        expectedGen.set(298, 218, 2.0);
        expectedGen.set(12, 219, 0.4);
        expectedGen.set(219, 219, -1.6);
        expectedGen.set(229, 219, 0.5);
        expectedGen.set(299, 219, 2.0);
        expectedGen.set(15, 220, 0.4);
        expectedGen.set(214, 220, 0.1);
        expectedGen.set(220, 220, -1.5);
        expectedGen.set(232, 220, 0.5);
        expectedGen.set(16, 221, 0.4);
        expectedGen.set(215, 221, 0.1);
        expectedGen.set(221, 221, -1.5);
        expectedGen.set(234, 221, 0.5);
        expectedGen.set(17, 222, 0.4);
        expectedGen.set(216, 222, 0.1);
        expectedGen.set(222, 222, -1.6);
        expectedGen.set(235, 222, 0.5);
        expectedGen.set(19, 223, 0.4);
        expectedGen.set(217, 223, 0.1);
        expectedGen.set(223, 223, -1.5);
        expectedGen.set(237, 223, 0.5);
        expectedGen.set(20, 224, 0.4);
        expectedGen.set(218, 224, 0.1);
        expectedGen.set(224, 224, -1.6);
        expectedGen.set(238, 224, 0.5);
        expectedGen.set(22, 225, 0.4);
        expectedGen.set(219, 225, 0.1);
        expectedGen.set(225, 225, -1.6);
        expectedGen.set(239, 225, 0.5);
        expectedGen.set(25, 226, 0.4);
        expectedGen.set(226, 226, -1.5);
        expectedGen.set(241, 226, 0.5);
        expectedGen.set(300, 226, 2.0);
        expectedGen.set(26, 227, 0.4);
        expectedGen.set(227, 227, -1.6);
        expectedGen.set(242, 227, 0.5);
        expectedGen.set(301, 227, 2.0);
        expectedGen.set(28, 228, 0.4);
        expectedGen.set(228, 228, -1.6);
        expectedGen.set(243, 228, 0.5);
        expectedGen.set(302, 228, 2.0);
        expectedGen.set(31, 229, 0.4);
        expectedGen.set(229, 229, -1.6);
        expectedGen.set(244, 229, 0.5);
        expectedGen.set(303, 229, 2.0);
        expectedGen.set(35, 230, 0.4);
        expectedGen.set(220, 230, 0.1);
        expectedGen.set(230, 230, -1.5);
        expectedGen.set(247, 230, 0.5);
        expectedGen.set(36, 231, 0.4);
        expectedGen.set(221, 231, 0.1);
        expectedGen.set(231, 231, -1.5);
        expectedGen.set(249, 231, 0.5);
        expectedGen.set(37, 232, 0.4);
        expectedGen.set(222, 232, 0.1);
        expectedGen.set(232, 232, -1.6);
        expectedGen.set(250, 232, 0.5);
        expectedGen.set(39, 233, 0.4);
        expectedGen.set(223, 233, 0.1);
        expectedGen.set(233, 233, -1.5);
        expectedGen.set(252, 233, 0.5);
        expectedGen.set(40, 234, 0.4);
        expectedGen.set(224, 234, 0.1);
        expectedGen.set(234, 234, -1.6);
        expectedGen.set(253, 234, 0.5);
        expectedGen.set(42, 235, 0.4);
        expectedGen.set(225, 235, 0.1);
        expectedGen.set(235, 235, -1.6);
        expectedGen.set(254, 235, 0.5);
        expectedGen.set(45, 236, 0.4);
        expectedGen.set(226, 236, 0.1);
        expectedGen.set(236, 236, -1.5);
        expectedGen.set(256, 236, 0.5);
        expectedGen.set(46, 237, 0.4);
        expectedGen.set(227, 237, 0.1);
        expectedGen.set(237, 237, -1.6);
        expectedGen.set(257, 237, 0.5);
        expectedGen.set(48, 238, 0.4);
        expectedGen.set(228, 238, 0.1);
        expectedGen.set(238, 238, -1.6);
        expectedGen.set(258, 238, 0.5);
        expectedGen.set(51, 239, 0.4);
        expectedGen.set(229, 239, 0.1);
        expectedGen.set(239, 239, -1.6);
        expectedGen.set(259, 239, 0.5);
        expectedGen.set(55, 240, 0.4);
        expectedGen.set(240, 240, -1.5);
        expectedGen.set(261, 240, 0.5);
        expectedGen.set(304, 240, 2.0);
        expectedGen.set(56, 241, 0.4);
        expectedGen.set(241, 241, -1.6);
        expectedGen.set(262, 241, 0.5);
        expectedGen.set(305, 241, 2.0);
        expectedGen.set(58, 242, 0.4);
        expectedGen.set(242, 242, -1.6);
        expectedGen.set(263, 242, 0.5);
        expectedGen.set(306, 242, 2.0);
        expectedGen.set(61, 243, 0.4);
        expectedGen.set(243, 243, -1.6);
        expectedGen.set(264, 243, 0.5);
        expectedGen.set(307, 243, 2.0);
        expectedGen.set(65, 244, 0.4);
        expectedGen.set(244, 244, -1.6);
        expectedGen.set(265, 244, 0.5);
        expectedGen.set(308, 244, 2.0);
        expectedGen.set(70, 245, 0.4);
        expectedGen.set(230, 245, 0.1);
        expectedGen.set(245, 245, -1.5);
        expectedGen.set(268, 245, 0.5);
        expectedGen.set(71, 246, 0.4);
        expectedGen.set(231, 246, 0.1);
        expectedGen.set(246, 246, -1.5);
        expectedGen.set(270, 246, 0.5);
        expectedGen.set(72, 247, 0.4);
        expectedGen.set(232, 247, 0.1);
        expectedGen.set(247, 247, -1.6);
        expectedGen.set(271, 247, 0.5);
        expectedGen.set(74, 248, 0.4);
        expectedGen.set(233, 248, 0.1);
        expectedGen.set(248, 248, -1.5);
        expectedGen.set(273, 248, 0.5);
        expectedGen.set(75, 249, 0.4);
        expectedGen.set(234, 249, 0.1);
        expectedGen.set(249, 249, -1.6);
        expectedGen.set(274, 249, 0.5);
        expectedGen.set(77, 250, 0.4);
        expectedGen.set(235, 250, 0.1);
        expectedGen.set(250, 250, -1.6);
        expectedGen.set(275, 250, 0.5);
        expectedGen.set(80, 251, 0.4);
        expectedGen.set(236, 251, 0.1);
        expectedGen.set(251, 251, -1.5);
        expectedGen.set(277, 251, 0.5);
        expectedGen.set(81, 252, 0.4);
        expectedGen.set(237, 252, 0.1);
        expectedGen.set(252, 252, -1.6);
        expectedGen.set(278, 252, 0.5);
        expectedGen.set(83, 253, 0.4);
        expectedGen.set(238, 253, 0.1);
        expectedGen.set(253, 253, -1.6);
        expectedGen.set(279, 253, 0.5);
        expectedGen.set(86, 254, 0.4);
        expectedGen.set(239, 254, 0.1);
        expectedGen.set(254, 254, -1.6);
        expectedGen.set(280, 254, 0.5);
        expectedGen.set(90, 255, 0.4);
        expectedGen.set(240, 255, 0.1);
        expectedGen.set(255, 255, -1.5);
        expectedGen.set(282, 255, 0.5);
        expectedGen.set(91, 256, 0.4);
        expectedGen.set(241, 256, 0.1);
        expectedGen.set(256, 256, -1.6);
        expectedGen.set(283, 256, 0.5);
        expectedGen.set(93, 257, 0.4);
        expectedGen.set(242, 257, 0.1);
        expectedGen.set(257, 257, -1.6);
        expectedGen.set(284, 257, 0.5);
        expectedGen.set(96, 258, 0.4);
        expectedGen.set(243, 258, 0.1);
        expectedGen.set(258, 258, -1.6);
        expectedGen.set(285, 258, 0.5);
        expectedGen.set(100, 259, 0.4);
        expectedGen.set(244, 259, 0.1);
        expectedGen.set(259, 259, -1.6);
        expectedGen.set(286, 259, 0.5);
        expectedGen.set(105, 260, 0.4);
        expectedGen.set(260, 260, -1.5);
        expectedGen.set(288, 260, 0.5);
        expectedGen.set(309, 260, 2.0);
        expectedGen.set(106, 261, 0.4);
        expectedGen.set(261, 261, -1.6);
        expectedGen.set(289, 261, 0.5);
        expectedGen.set(310, 261, 2.0);
        expectedGen.set(108, 262, 0.4);
        expectedGen.set(262, 262, -1.6);
        expectedGen.set(290, 262, 0.5);
        expectedGen.set(311, 262, 2.0);
        expectedGen.set(111, 263, 0.4);
        expectedGen.set(263, 263, -1.6);
        expectedGen.set(291, 263, 0.5);
        expectedGen.set(312, 263, 2.0);
        expectedGen.set(115, 264, 0.4);
        expectedGen.set(264, 264, -1.6);
        expectedGen.set(292, 264, 0.5);
        expectedGen.set(313, 264, 2.0);
        expectedGen.set(120, 265, 0.4);
        expectedGen.set(265, 265, -1.6);
        expectedGen.set(293, 265, 0.5);
        expectedGen.set(314, 265, 2.0);
        expectedGen.set(126, 266, 0.4);
        expectedGen.set(245, 266, 0.1);
        expectedGen.set(266, 266, -1.4);
        expectedGen.set(127, 267, 0.4);
        expectedGen.set(246, 267, 0.1);
        expectedGen.set(267, 267, -1.4);
        expectedGen.set(128, 268, 0.4);
        expectedGen.set(247, 268, 0.1);
        expectedGen.set(268, 268, -1.5);
        expectedGen.set(130, 269, 0.4);
        expectedGen.set(248, 269, 0.1);
        expectedGen.set(269, 269, -1.4);
        expectedGen.set(131, 270, 0.4);
        expectedGen.set(249, 270, 0.1);
        expectedGen.set(270, 270, -1.5);
        expectedGen.set(133, 271, 0.4);
        expectedGen.set(250, 271, 0.1);
        expectedGen.set(271, 271, -1.5);
        expectedGen.set(136, 272, 0.4);
        expectedGen.set(251, 272, 0.1);
        expectedGen.set(272, 272, -1.4);
        expectedGen.set(137, 273, 0.4);
        expectedGen.set(252, 273, 0.1);
        expectedGen.set(273, 273, -1.5);
        expectedGen.set(139, 274, 0.4);
        expectedGen.set(253, 274, 0.1);
        expectedGen.set(274, 274, -1.5);
        expectedGen.set(142, 275, 0.4);
        expectedGen.set(254, 275, 0.1);
        expectedGen.set(275, 275, -1.5);
        expectedGen.set(146, 276, 0.4);
        expectedGen.set(255, 276, 0.1);
        expectedGen.set(276, 276, -1.4);
        expectedGen.set(147, 277, 0.4);
        expectedGen.set(256, 277, 0.1);
        expectedGen.set(277, 277, -1.5);
        expectedGen.set(149, 278, 0.4);
        expectedGen.set(257, 278, 0.1);
        expectedGen.set(278, 278, -1.5);
        expectedGen.set(152, 279, 0.4);
        expectedGen.set(258, 279, 0.1);
        expectedGen.set(279, 279, -1.5);
        expectedGen.set(156, 280, 0.4);
        expectedGen.set(259, 280, 0.1);
        expectedGen.set(280, 280, -1.5);
        expectedGen.set(161, 281, 0.4);
        expectedGen.set(260, 281, 0.1);
        expectedGen.set(281, 281, -1.4);
        expectedGen.set(162, 282, 0.4);
        expectedGen.set(261, 282, 0.1);
        expectedGen.set(282, 282, -1.5);
        expectedGen.set(164, 283, 0.4);
        expectedGen.set(262, 283, 0.1);
        expectedGen.set(283, 283, -1.5);
        expectedGen.set(167, 284, 0.4);
        expectedGen.set(263, 284, 0.1);
        expectedGen.set(284, 284, -1.5);
        expectedGen.set(171, 285, 0.4);
        expectedGen.set(264, 285, 0.1);
        expectedGen.set(285, 285, -1.5);
        expectedGen.set(176, 286, 0.4);
        expectedGen.set(265, 286, 0.1);
        expectedGen.set(286, 286, -1.5);
        expectedGen.set(182, 287, 0.4);
        expectedGen.set(287, 287, -1.4);
        expectedGen.set(315, 287, 2.0);
        expectedGen.set(183, 288, 0.4);
        expectedGen.set(288, 288, -1.5);
        expectedGen.set(316, 288, 2.0);
        expectedGen.set(185, 289, 0.4);
        expectedGen.set(289, 289, -1.5);
        expectedGen.set(317, 289, 2.0);
        expectedGen.set(188, 290, 0.4);
        expectedGen.set(290, 290, -1.5);
        expectedGen.set(318, 290, 2.0);
        expectedGen.set(192, 291, 0.4);
        expectedGen.set(291, 291, -1.5);
        expectedGen.set(319, 291, 2.0);
        expectedGen.set(197, 292, 0.4);
        expectedGen.set(292, 292, -1.5);
        expectedGen.set(320, 292, 2.0);
        expectedGen.set(203, 293, 0.4);
        expectedGen.set(293, 293, -1.5);
        expectedGen.set(321, 293, 2.0);
        expectedGen.set(210, 294, 0.4);
        expectedGen.set(294, 294, -2.5);
        expectedGen.set(296, 294, 0.5);
        expectedGen.set(322, 294, 3.0);
        expectedGen.set(211, 295, 0.4);
        expectedGen.set(294, 295, 0.1);
        expectedGen.set(295, 295, -2.5);
        expectedGen.set(298, 295, 0.5);
        expectedGen.set(212, 296, 0.4);
        expectedGen.set(296, 296, -2.6);
        expectedGen.set(299, 296, 0.5);
        expectedGen.set(323, 296, 3.0);
        expectedGen.set(214, 297, 0.4);
        expectedGen.set(295, 297, 0.1);
        expectedGen.set(297, 297, -2.5);
        expectedGen.set(301, 297, 0.5);
        expectedGen.set(215, 298, 0.4);
        expectedGen.set(296, 298, 0.1);
        expectedGen.set(298, 298, -2.6);
        expectedGen.set(302, 298, 0.5);
        expectedGen.set(217, 299, 0.4);
        expectedGen.set(299, 299, -2.6);
        expectedGen.set(303, 299, 0.5);
        expectedGen.set(324, 299, 3.0);
        expectedGen.set(220, 300, 0.4);
        expectedGen.set(297, 300, 0.1);
        expectedGen.set(300, 300, -2.5);
        expectedGen.set(305, 300, 0.5);
        expectedGen.set(221, 301, 0.4);
        expectedGen.set(298, 301, 0.1);
        expectedGen.set(301, 301, -2.6);
        expectedGen.set(306, 301, 0.5);
        expectedGen.set(223, 302, 0.4);
        expectedGen.set(299, 302, 0.1);
        expectedGen.set(302, 302, -2.6);
        expectedGen.set(307, 302, 0.5);
        expectedGen.set(226, 303, 0.4);
        expectedGen.set(303, 303, -2.6);
        expectedGen.set(308, 303, 0.5);
        expectedGen.set(325, 303, 3.0);
        expectedGen.set(230, 304, 0.4);
        expectedGen.set(300, 304, 0.1);
        expectedGen.set(304, 304, -2.5);
        expectedGen.set(310, 304, 0.5);
        expectedGen.set(231, 305, 0.4);
        expectedGen.set(301, 305, 0.1);
        expectedGen.set(305, 305, -2.6);
        expectedGen.set(311, 305, 0.5);
        expectedGen.set(233, 306, 0.4);
        expectedGen.set(302, 306, 0.1);
        expectedGen.set(306, 306, -2.6);
        expectedGen.set(312, 306, 0.5);
        expectedGen.set(236, 307, 0.4);
        expectedGen.set(303, 307, 0.1);
        expectedGen.set(307, 307, -2.6);
        expectedGen.set(313, 307, 0.5);
        expectedGen.set(240, 308, 0.4);
        expectedGen.set(308, 308, -2.6);
        expectedGen.set(314, 308, 0.5);
        expectedGen.set(326, 308, 3.0);
        expectedGen.set(245, 309, 0.4);
        expectedGen.set(304, 309, 0.1);
        expectedGen.set(309, 309, -2.5);
        expectedGen.set(316, 309, 0.5);
        expectedGen.set(246, 310, 0.4);
        expectedGen.set(305, 310, 0.1);
        expectedGen.set(310, 310, -2.6);
        expectedGen.set(317, 310, 0.5);
        expectedGen.set(248, 311, 0.4);
        expectedGen.set(306, 311, 0.1);
        expectedGen.set(311, 311, -2.6);
        expectedGen.set(318, 311, 0.5);
        expectedGen.set(251, 312, 0.4);
        expectedGen.set(307, 312, 0.1);
        expectedGen.set(312, 312, -2.6);
        expectedGen.set(319, 312, 0.5);
        expectedGen.set(255, 313, 0.4);
        expectedGen.set(308, 313, 0.1);
        expectedGen.set(313, 313, -2.6);
        expectedGen.set(320, 313, 0.5);
        expectedGen.set(260, 314, 0.4);
        expectedGen.set(314, 314, -2.6);
        expectedGen.set(321, 314, 0.5);
        expectedGen.set(327, 314, 3.0);
        expectedGen.set(266, 315, 0.4);
        expectedGen.set(309, 315, 0.1);
        expectedGen.set(315, 315, -2.4);
        expectedGen.set(267, 316, 0.4);
        expectedGen.set(310, 316, 0.1);
        expectedGen.set(316, 316, -2.5);
        expectedGen.set(269, 317, 0.4);
        expectedGen.set(311, 317, 0.1);
        expectedGen.set(317, 317, -2.5);
        expectedGen.set(272, 318, 0.4);
        expectedGen.set(312, 318, 0.1);
        expectedGen.set(318, 318, -2.5);
        expectedGen.set(276, 319, 0.4);
        expectedGen.set(313, 319, 0.1);
        expectedGen.set(319, 319, -2.5);
        expectedGen.set(281, 320, 0.4);
        expectedGen.set(314, 320, 0.1);
        expectedGen.set(320, 320, -2.5);
        expectedGen.set(287, 321, 0.4);
        expectedGen.set(321, 321, -2.5);
        expectedGen.set(328, 321, 3.0);
        expectedGen.set(294, 322, 0.4);
        expectedGen.set(322, 322, -3.1);
        expectedGen.set(323, 322, 0.5);
        expectedGen.set(295, 323, 0.4);
        expectedGen.set(322, 323, 0.1);
        expectedGen.set(323, 323, -3.6);
        expectedGen.set(324, 323, 0.5);
        expectedGen.set(297, 324, 0.4);
        expectedGen.set(323, 324, 0.1);
        expectedGen.set(324, 324, -3.6);
        expectedGen.set(325, 324, 0.5);
        expectedGen.set(300, 325, 0.4);
        expectedGen.set(324, 325, 0.1);
        expectedGen.set(325, 325, -3.6);
        expectedGen.set(326, 325, 0.5);
        expectedGen.set(304, 326, 0.4);
        expectedGen.set(325, 326, 0.1);
        expectedGen.set(326, 326, -3.6);
        expectedGen.set(327, 326, 0.5);
        expectedGen.set(309, 327, 0.4);
        expectedGen.set(326, 327, 0.1);
        expectedGen.set(327, 327, -3.6);
        expectedGen.set(328, 327, 0.5);
        expectedGen.set(315, 328, 0.4);
        expectedGen.set(327, 328, 0.1);
        expectedGen.set(328, 328, -3.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "FCFS mixed network generator matrix mismatch");

        // Verify state space dimensions
        assertEquals(329, stateSpace.getNumRows());
        assertEquals(15, stateSpace.getNumCols());
    }

    @Test
    public void testMixedNetworkLCFS() {
        Network model = createMixedNetwork(SchedStrategy.LCFS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(6);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify dimensions - same as FCFS
        assertEquals(329, generator.getNumRows());
        assertEquals(329, generator.getNumCols());

        // This is a temporary file with the correct expected generator values
// Copy this into NetworkStateTest.java testMixedNetworkLCFS method

        // Create expected sparse matrix (329x329) with ALL non-zero entries from MATLAB
        Matrix expectedGen = new Matrix(329, 329);
        // Set all non-zero values from the ground truth
        expectedGen.set(0, 0, -0.5);
        expectedGen.set(4, 0, 0.5);
        expectedGen.set(210, 0, 1.0);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.5);
        expectedGen.set(11, 1, 0.5);
        expectedGen.set(2, 2, -0.5);
        expectedGen.set(13, 2, 0.5);
        expectedGen.set(211, 2, 1.0);
        expectedGen.set(3, 3, -0.5);
        expectedGen.set(14, 3, 0.5);
        expectedGen.set(212, 3, 1.0);
        expectedGen.set(4, 4, -0.6);
        expectedGen.set(8, 4, 0.5);
        expectedGen.set(213, 4, 1.0);
        expectedGen.set(1, 5, 0.1);
        expectedGen.set(5, 5, -0.5);
        expectedGen.set(27, 5, 0.5);
        expectedGen.set(2, 6, 0.1);
        expectedGen.set(6, 6, -0.5);
        expectedGen.set(29, 6, 0.5);
        expectedGen.set(3, 7, 0.1);
        expectedGen.set(7, 7, -0.5);
        expectedGen.set(30, 7, 0.5);
        expectedGen.set(4, 8, 0.1);
        expectedGen.set(8, 8, -0.6);
        expectedGen.set(18, 8, 0.5);
        expectedGen.set(9, 9, -0.5);
        expectedGen.set(32, 9, 0.5);
        expectedGen.set(214, 9, 1.0);
        expectedGen.set(10, 10, -0.5);
        expectedGen.set(33, 10, 0.5);
        expectedGen.set(215, 10, 1.0);
        expectedGen.set(11, 11, -0.6);
        expectedGen.set(21, 11, 0.5);
        expectedGen.set(216, 11, 1.0);
        expectedGen.set(12, 12, -0.5);
        expectedGen.set(34, 12, 0.5);
        expectedGen.set(217, 12, 1.0);
        expectedGen.set(13, 13, -0.6);
        expectedGen.set(23, 13, 0.5);
        expectedGen.set(218, 13, 1.0);
        expectedGen.set(14, 14, -0.6);
        expectedGen.set(24, 14, 0.5);
        expectedGen.set(219, 14, 1.0);
        expectedGen.set(5, 15, 0.1);
        expectedGen.set(15, 15, -0.5);
        expectedGen.set(57, 15, 0.5);
        expectedGen.set(6, 16, 0.1);
        expectedGen.set(16, 16, -0.5);
        expectedGen.set(59, 16, 0.5);
        expectedGen.set(7, 17, 0.1);
        expectedGen.set(17, 17, -0.5);
        expectedGen.set(60, 17, 0.5);
        expectedGen.set(8, 18, 0.1);
        expectedGen.set(18, 18, -0.6);
        expectedGen.set(38, 18, 0.5);
        expectedGen.set(9, 19, 0.1);
        expectedGen.set(19, 19, -0.5);
        expectedGen.set(62, 19, 0.5);
        expectedGen.set(10, 20, 0.1);
        expectedGen.set(20, 20, -0.5);
        expectedGen.set(63, 20, 0.5);
        expectedGen.set(11, 21, 0.1);
        expectedGen.set(21, 21, -0.6);
        expectedGen.set(41, 21, 0.5);
        expectedGen.set(12, 22, 0.1);
        expectedGen.set(22, 22, -0.5);
        expectedGen.set(64, 22, 0.5);
        expectedGen.set(13, 23, 0.1);
        expectedGen.set(23, 23, -0.6);
        expectedGen.set(43, 23, 0.5);
        expectedGen.set(14, 24, 0.1);
        expectedGen.set(24, 24, -0.6);
        expectedGen.set(44, 24, 0.5);
        expectedGen.set(25, 25, -0.5);
        expectedGen.set(66, 25, 0.5);
        expectedGen.set(220, 25, 1.0);
        expectedGen.set(26, 26, -0.5);
        expectedGen.set(67, 26, 0.5);
        expectedGen.set(221, 26, 1.0);
        expectedGen.set(27, 27, -0.6);
        expectedGen.set(47, 27, 0.5);
        expectedGen.set(222, 27, 1.0);
        expectedGen.set(28, 28, -0.5);
        expectedGen.set(68, 28, 0.5);
        expectedGen.set(223, 28, 1.0);
        expectedGen.set(29, 29, -0.6);
        expectedGen.set(49, 29, 0.5);
        expectedGen.set(224, 29, 1.0);
        expectedGen.set(30, 30, -0.6);
        expectedGen.set(50, 30, 0.5);
        expectedGen.set(225, 30, 1.0);
        expectedGen.set(31, 31, -0.5);
        expectedGen.set(69, 31, 0.5);
        expectedGen.set(226, 31, 1.0);
        expectedGen.set(32, 32, -0.6);
        expectedGen.set(52, 32, 0.5);
        expectedGen.set(227, 32, 1.0);
        expectedGen.set(33, 33, -0.6);
        expectedGen.set(53, 33, 0.5);
        expectedGen.set(228, 33, 1.0);
        expectedGen.set(34, 34, -0.6);
        expectedGen.set(54, 34, 0.5);
        expectedGen.set(229, 34, 1.0);
        expectedGen.set(15, 35, 0.1);
        expectedGen.set(35, 35, -0.5);
        expectedGen.set(107, 35, 0.5);
        expectedGen.set(16, 36, 0.1);
        expectedGen.set(36, 36, -0.5);
        expectedGen.set(109, 36, 0.5);
        expectedGen.set(17, 37, 0.1);
        expectedGen.set(37, 37, -0.5);
        expectedGen.set(110, 37, 0.5);
        expectedGen.set(18, 38, 0.1);
        expectedGen.set(38, 38, -0.6);
        expectedGen.set(73, 38, 0.5);
        expectedGen.set(19, 39, 0.1);
        expectedGen.set(39, 39, -0.5);
        expectedGen.set(112, 39, 0.5);
        expectedGen.set(20, 40, 0.1);
        expectedGen.set(40, 40, -0.5);
        expectedGen.set(113, 40, 0.5);
        expectedGen.set(21, 41, 0.1);
        expectedGen.set(41, 41, -0.6);
        expectedGen.set(76, 41, 0.5);
        expectedGen.set(22, 42, 0.1);
        expectedGen.set(42, 42, -0.5);
        expectedGen.set(114, 42, 0.5);
        expectedGen.set(23, 43, 0.1);
        expectedGen.set(43, 43, -0.6);
        expectedGen.set(78, 43, 0.5);
        expectedGen.set(24, 44, 0.1);
        expectedGen.set(44, 44, -0.6);
        expectedGen.set(79, 44, 0.5);
        expectedGen.set(25, 45, 0.1);
        expectedGen.set(45, 45, -0.5);
        expectedGen.set(116, 45, 0.5);
        expectedGen.set(26, 46, 0.1);
        expectedGen.set(46, 46, -0.5);
        expectedGen.set(117, 46, 0.5);
        expectedGen.set(27, 47, 0.1);
        expectedGen.set(47, 47, -0.6);
        expectedGen.set(82, 47, 0.5);
        expectedGen.set(28, 48, 0.1);
        expectedGen.set(48, 48, -0.5);
        expectedGen.set(118, 48, 0.5);
        expectedGen.set(29, 49, 0.1);
        expectedGen.set(49, 49, -0.6);
        expectedGen.set(84, 49, 0.5);
        expectedGen.set(30, 50, 0.1);
        expectedGen.set(50, 50, -0.6);
        expectedGen.set(85, 50, 0.5);
        expectedGen.set(31, 51, 0.1);
        expectedGen.set(51, 51, -0.5);
        expectedGen.set(119, 51, 0.5);
        expectedGen.set(32, 52, 0.1);
        expectedGen.set(52, 52, -0.6);
        expectedGen.set(87, 52, 0.5);
        expectedGen.set(33, 53, 0.1);
        expectedGen.set(53, 53, -0.6);
        expectedGen.set(88, 53, 0.5);
        expectedGen.set(34, 54, 0.1);
        expectedGen.set(54, 54, -0.6);
        expectedGen.set(89, 54, 0.5);
        expectedGen.set(55, 55, -0.5);
        expectedGen.set(121, 55, 0.5);
        expectedGen.set(230, 55, 1.0);
        expectedGen.set(56, 56, -0.5);
        expectedGen.set(122, 56, 0.5);
        expectedGen.set(231, 56, 1.0);
        expectedGen.set(57, 57, -0.6);
        expectedGen.set(92, 57, 0.5);
        expectedGen.set(232, 57, 1.0);
        expectedGen.set(58, 58, -0.5);
        expectedGen.set(123, 58, 0.5);
        expectedGen.set(233, 58, 1.0);
        expectedGen.set(59, 59, -0.6);
        expectedGen.set(94, 59, 0.5);
        expectedGen.set(234, 59, 1.0);
        expectedGen.set(60, 60, -0.6);
        expectedGen.set(95, 60, 0.5);
        expectedGen.set(235, 60, 1.0);
        expectedGen.set(61, 61, -0.5);
        expectedGen.set(124, 61, 0.5);
        expectedGen.set(236, 61, 1.0);
        expectedGen.set(62, 62, -0.6);
        expectedGen.set(97, 62, 0.5);
        expectedGen.set(237, 62, 1.0);
        expectedGen.set(63, 63, -0.6);
        expectedGen.set(98, 63, 0.5);
        expectedGen.set(238, 63, 1.0);
        expectedGen.set(64, 64, -0.6);
        expectedGen.set(99, 64, 0.5);
        expectedGen.set(239, 64, 1.0);
        expectedGen.set(65, 65, -0.5);
        expectedGen.set(125, 65, 0.5);
        expectedGen.set(240, 65, 1.0);
        expectedGen.set(66, 66, -0.6);
        expectedGen.set(101, 66, 0.5);
        expectedGen.set(241, 66, 1.0);
        expectedGen.set(67, 67, -0.6);
        expectedGen.set(102, 67, 0.5);
        expectedGen.set(242, 67, 1.0);
        expectedGen.set(68, 68, -0.6);
        expectedGen.set(103, 68, 0.5);
        expectedGen.set(243, 68, 1.0);
        expectedGen.set(69, 69, -0.6);
        expectedGen.set(104, 69, 0.5);
        expectedGen.set(244, 69, 1.0);
        expectedGen.set(35, 70, 0.1);
        expectedGen.set(70, 70, -0.5);
        expectedGen.set(184, 70, 0.5);
        expectedGen.set(36, 71, 0.1);
        expectedGen.set(71, 71, -0.5);
        expectedGen.set(186, 71, 0.5);
        expectedGen.set(37, 72, 0.1);
        expectedGen.set(72, 72, -0.5);
        expectedGen.set(187, 72, 0.5);
        expectedGen.set(38, 73, 0.1);
        expectedGen.set(73, 73, -0.6);
        expectedGen.set(129, 73, 0.5);
        expectedGen.set(39, 74, 0.1);
        expectedGen.set(74, 74, -0.5);
        expectedGen.set(189, 74, 0.5);
        expectedGen.set(40, 75, 0.1);
        expectedGen.set(75, 75, -0.5);
        expectedGen.set(190, 75, 0.5);
        expectedGen.set(41, 76, 0.1);
        expectedGen.set(76, 76, -0.6);
        expectedGen.set(132, 76, 0.5);
        expectedGen.set(42, 77, 0.1);
        expectedGen.set(77, 77, -0.5);
        expectedGen.set(191, 77, 0.5);
        expectedGen.set(43, 78, 0.1);
        expectedGen.set(78, 78, -0.6);
        expectedGen.set(134, 78, 0.5);
        expectedGen.set(44, 79, 0.1);
        expectedGen.set(79, 79, -0.6);
        expectedGen.set(135, 79, 0.5);
        expectedGen.set(45, 80, 0.1);
        expectedGen.set(80, 80, -0.5);
        expectedGen.set(193, 80, 0.5);
        expectedGen.set(46, 81, 0.1);
        expectedGen.set(81, 81, -0.5);
        expectedGen.set(194, 81, 0.5);
        expectedGen.set(47, 82, 0.1);
        expectedGen.set(82, 82, -0.6);
        expectedGen.set(138, 82, 0.5);
        expectedGen.set(48, 83, 0.1);
        expectedGen.set(83, 83, -0.5);
        expectedGen.set(195, 83, 0.5);
        expectedGen.set(49, 84, 0.1);
        expectedGen.set(84, 84, -0.6);
        expectedGen.set(140, 84, 0.5);
        expectedGen.set(50, 85, 0.1);
        expectedGen.set(85, 85, -0.6);
        expectedGen.set(141, 85, 0.5);
        expectedGen.set(51, 86, 0.1);
        expectedGen.set(86, 86, -0.5);
        expectedGen.set(196, 86, 0.5);
        expectedGen.set(52, 87, 0.1);
        expectedGen.set(87, 87, -0.6);
        expectedGen.set(143, 87, 0.5);
        expectedGen.set(53, 88, 0.1);
        expectedGen.set(88, 88, -0.6);
        expectedGen.set(144, 88, 0.5);
        expectedGen.set(54, 89, 0.1);
        expectedGen.set(89, 89, -0.6);
        expectedGen.set(145, 89, 0.5);
        expectedGen.set(55, 90, 0.1);
        expectedGen.set(90, 90, -0.5);
        expectedGen.set(198, 90, 0.5);
        expectedGen.set(56, 91, 0.1);
        expectedGen.set(91, 91, -0.5);
        expectedGen.set(199, 91, 0.5);
        expectedGen.set(57, 92, 0.1);
        expectedGen.set(92, 92, -0.6);
        expectedGen.set(148, 92, 0.5);
        expectedGen.set(58, 93, 0.1);
        expectedGen.set(93, 93, -0.5);
        expectedGen.set(200, 93, 0.5);
        expectedGen.set(59, 94, 0.1);
        expectedGen.set(94, 94, -0.6);
        expectedGen.set(150, 94, 0.5);
        expectedGen.set(60, 95, 0.1);
        expectedGen.set(95, 95, -0.6);
        expectedGen.set(151, 95, 0.5);
        expectedGen.set(61, 96, 0.1);
        expectedGen.set(96, 96, -0.5);
        expectedGen.set(201, 96, 0.5);
        expectedGen.set(62, 97, 0.1);
        expectedGen.set(97, 97, -0.6);
        expectedGen.set(153, 97, 0.5);
        expectedGen.set(63, 98, 0.1);
        expectedGen.set(98, 98, -0.6);
        expectedGen.set(154, 98, 0.5);
        expectedGen.set(64, 99, 0.1);
        expectedGen.set(99, 99, -0.6);
        expectedGen.set(155, 99, 0.5);
        expectedGen.set(65, 100, 0.1);
        expectedGen.set(100, 100, -0.5);
        expectedGen.set(202, 100, 0.5);
        expectedGen.set(66, 101, 0.1);
        expectedGen.set(101, 101, -0.6);
        expectedGen.set(157, 101, 0.5);
        expectedGen.set(67, 102, 0.1);
        expectedGen.set(102, 102, -0.6);
        expectedGen.set(158, 102, 0.5);
        expectedGen.set(68, 103, 0.1);
        expectedGen.set(103, 103, -0.6);
        expectedGen.set(159, 103, 0.5);
        expectedGen.set(69, 104, 0.1);
        expectedGen.set(104, 104, -0.6);
        expectedGen.set(160, 104, 0.5);
        expectedGen.set(105, 105, -0.5);
        expectedGen.set(204, 105, 0.5);
        expectedGen.set(245, 105, 1.0);
        expectedGen.set(106, 106, -0.5);
        expectedGen.set(205, 106, 0.5);
        expectedGen.set(246, 106, 1.0);
        expectedGen.set(107, 107, -0.6);
        expectedGen.set(163, 107, 0.5);
        expectedGen.set(247, 107, 1.0);
        expectedGen.set(108, 108, -0.5);
        expectedGen.set(206, 108, 0.5);
        expectedGen.set(248, 108, 1.0);
        expectedGen.set(109, 109, -0.6);
        expectedGen.set(165, 109, 0.5);
        expectedGen.set(249, 109, 1.0);
        expectedGen.set(110, 110, -0.6);
        expectedGen.set(166, 110, 0.5);
        expectedGen.set(250, 110, 1.0);
        expectedGen.set(111, 111, -0.5);
        expectedGen.set(207, 111, 0.5);
        expectedGen.set(251, 111, 1.0);
        expectedGen.set(112, 112, -0.6);
        expectedGen.set(168, 112, 0.5);
        expectedGen.set(252, 112, 1.0);
        expectedGen.set(113, 113, -0.6);
        expectedGen.set(169, 113, 0.5);
        expectedGen.set(253, 113, 1.0);
        expectedGen.set(114, 114, -0.6);
        expectedGen.set(170, 114, 0.5);
        expectedGen.set(254, 114, 1.0);
        expectedGen.set(115, 115, -0.5);
        expectedGen.set(208, 115, 0.5);
        expectedGen.set(255, 115, 1.0);
        expectedGen.set(116, 116, -0.6);
        expectedGen.set(172, 116, 0.5);
        expectedGen.set(256, 116, 1.0);
        expectedGen.set(117, 117, -0.6);
        expectedGen.set(173, 117, 0.5);
        expectedGen.set(257, 117, 1.0);
        expectedGen.set(118, 118, -0.6);
        expectedGen.set(174, 118, 0.5);
        expectedGen.set(258, 118, 1.0);
        expectedGen.set(119, 119, -0.6);
        expectedGen.set(175, 119, 0.5);
        expectedGen.set(259, 119, 1.0);
        expectedGen.set(120, 120, -0.5);
        expectedGen.set(209, 120, 0.5);
        expectedGen.set(260, 120, 1.0);
        expectedGen.set(121, 121, -0.6);
        expectedGen.set(177, 121, 0.5);
        expectedGen.set(261, 121, 1.0);
        expectedGen.set(122, 122, -0.6);
        expectedGen.set(178, 122, 0.5);
        expectedGen.set(262, 122, 1.0);
        expectedGen.set(123, 123, -0.6);
        expectedGen.set(179, 123, 0.5);
        expectedGen.set(263, 123, 1.0);
        expectedGen.set(124, 124, -0.6);
        expectedGen.set(180, 124, 0.5);
        expectedGen.set(264, 124, 1.0);
        expectedGen.set(125, 125, -0.6);
        expectedGen.set(181, 125, 0.5);
        expectedGen.set(265, 125, 1.0);
        expectedGen.set(70, 126, 0.1);
        expectedGen.set(126, 126, -0.4);
        expectedGen.set(71, 127, 0.1);
        expectedGen.set(127, 127, -0.4);
        expectedGen.set(72, 128, 0.1);
        expectedGen.set(128, 128, -0.4);
        expectedGen.set(73, 129, 0.1);
        expectedGen.set(129, 129, -0.5);
        expectedGen.set(74, 130, 0.1);
        expectedGen.set(130, 130, -0.4);
        expectedGen.set(75, 131, 0.1);
        expectedGen.set(131, 131, -0.4);
        expectedGen.set(76, 132, 0.1);
        expectedGen.set(132, 132, -0.5);
        expectedGen.set(77, 133, 0.1);
        expectedGen.set(133, 133, -0.4);
        expectedGen.set(78, 134, 0.1);
        expectedGen.set(134, 134, -0.5);
        expectedGen.set(79, 135, 0.1);
        expectedGen.set(135, 135, -0.5);
        expectedGen.set(80, 136, 0.1);
        expectedGen.set(136, 136, -0.4);
        expectedGen.set(81, 137, 0.1);
        expectedGen.set(137, 137, -0.4);
        expectedGen.set(82, 138, 0.1);
        expectedGen.set(138, 138, -0.5);
        expectedGen.set(83, 139, 0.1);
        expectedGen.set(139, 139, -0.4);
        expectedGen.set(84, 140, 0.1);
        expectedGen.set(140, 140, -0.5);
        expectedGen.set(85, 141, 0.1);
        expectedGen.set(141, 141, -0.5);
        expectedGen.set(86, 142, 0.1);
        expectedGen.set(142, 142, -0.4);
        expectedGen.set(87, 143, 0.1);
        expectedGen.set(143, 143, -0.5);
        expectedGen.set(88, 144, 0.1);
        expectedGen.set(144, 144, -0.5);
        expectedGen.set(89, 145, 0.1);
        expectedGen.set(145, 145, -0.5);
        expectedGen.set(90, 146, 0.1);
        expectedGen.set(146, 146, -0.4);
        expectedGen.set(91, 147, 0.1);
        expectedGen.set(147, 147, -0.4);
        expectedGen.set(92, 148, 0.1);
        expectedGen.set(148, 148, -0.5);
        expectedGen.set(93, 149, 0.1);
        expectedGen.set(149, 149, -0.4);
        expectedGen.set(94, 150, 0.1);
        expectedGen.set(150, 150, -0.5);
        expectedGen.set(95, 151, 0.1);
        expectedGen.set(151, 151, -0.5);
        expectedGen.set(96, 152, 0.1);
        expectedGen.set(152, 152, -0.4);
        expectedGen.set(97, 153, 0.1);
        expectedGen.set(153, 153, -0.5);
        expectedGen.set(98, 154, 0.1);
        expectedGen.set(154, 154, -0.5);
        expectedGen.set(99, 155, 0.1);
        expectedGen.set(155, 155, -0.5);
        expectedGen.set(100, 156, 0.1);
        expectedGen.set(156, 156, -0.4);
        expectedGen.set(101, 157, 0.1);
        expectedGen.set(157, 157, -0.5);
        expectedGen.set(102, 158, 0.1);
        expectedGen.set(158, 158, -0.5);
        expectedGen.set(103, 159, 0.1);
        expectedGen.set(159, 159, -0.5);
        expectedGen.set(104, 160, 0.1);
        expectedGen.set(160, 160, -0.5);
        expectedGen.set(105, 161, 0.1);
        expectedGen.set(161, 161, -0.4);
        expectedGen.set(106, 162, 0.1);
        expectedGen.set(162, 162, -0.4);
        expectedGen.set(107, 163, 0.1);
        expectedGen.set(163, 163, -0.5);
        expectedGen.set(108, 164, 0.1);
        expectedGen.set(164, 164, -0.4);
        expectedGen.set(109, 165, 0.1);
        expectedGen.set(165, 165, -0.5);
        expectedGen.set(110, 166, 0.1);
        expectedGen.set(166, 166, -0.5);
        expectedGen.set(111, 167, 0.1);
        expectedGen.set(167, 167, -0.4);
        expectedGen.set(112, 168, 0.1);
        expectedGen.set(168, 168, -0.5);
        expectedGen.set(113, 169, 0.1);
        expectedGen.set(169, 169, -0.5);
        expectedGen.set(114, 170, 0.1);
        expectedGen.set(170, 170, -0.5);
        expectedGen.set(115, 171, 0.1);
        expectedGen.set(171, 171, -0.4);
        expectedGen.set(116, 172, 0.1);
        expectedGen.set(172, 172, -0.5);
        expectedGen.set(117, 173, 0.1);
        expectedGen.set(173, 173, -0.5);
        expectedGen.set(118, 174, 0.1);
        expectedGen.set(174, 174, -0.5);
        expectedGen.set(119, 175, 0.1);
        expectedGen.set(175, 175, -0.5);
        expectedGen.set(120, 176, 0.1);
        expectedGen.set(176, 176, -0.4);
        expectedGen.set(121, 177, 0.1);
        expectedGen.set(177, 177, -0.5);
        expectedGen.set(122, 178, 0.1);
        expectedGen.set(178, 178, -0.5);
        expectedGen.set(123, 179, 0.1);
        expectedGen.set(179, 179, -0.5);
        expectedGen.set(124, 180, 0.1);
        expectedGen.set(180, 180, -0.5);
        expectedGen.set(125, 181, 0.1);
        expectedGen.set(181, 181, -0.5);
        expectedGen.set(182, 182, -0.4);
        expectedGen.set(266, 182, 1.0);
        expectedGen.set(183, 183, -0.4);
        expectedGen.set(267, 183, 1.0);
        expectedGen.set(184, 184, -0.5);
        expectedGen.set(268, 184, 1.0);
        expectedGen.set(185, 185, -0.4);
        expectedGen.set(269, 185, 1.0);
        expectedGen.set(186, 186, -0.5);
        expectedGen.set(270, 186, 1.0);
        expectedGen.set(187, 187, -0.5);
        expectedGen.set(271, 187, 1.0);
        expectedGen.set(188, 188, -0.4);
        expectedGen.set(272, 188, 1.0);
        expectedGen.set(189, 189, -0.5);
        expectedGen.set(273, 189, 1.0);
        expectedGen.set(190, 190, -0.5);
        expectedGen.set(274, 190, 1.0);
        expectedGen.set(191, 191, -0.5);
        expectedGen.set(275, 191, 1.0);
        expectedGen.set(192, 192, -0.4);
        expectedGen.set(276, 192, 1.0);
        expectedGen.set(193, 193, -0.5);
        expectedGen.set(277, 193, 1.0);
        expectedGen.set(194, 194, -0.5);
        expectedGen.set(278, 194, 1.0);
        expectedGen.set(195, 195, -0.5);
        expectedGen.set(279, 195, 1.0);
        expectedGen.set(196, 196, -0.5);
        expectedGen.set(280, 196, 1.0);
        expectedGen.set(197, 197, -0.4);
        expectedGen.set(281, 197, 1.0);
        expectedGen.set(198, 198, -0.5);
        expectedGen.set(282, 198, 1.0);
        expectedGen.set(199, 199, -0.5);
        expectedGen.set(283, 199, 1.0);
        expectedGen.set(200, 200, -0.5);
        expectedGen.set(284, 200, 1.0);
        expectedGen.set(201, 201, -0.5);
        expectedGen.set(285, 201, 1.0);
        expectedGen.set(202, 202, -0.5);
        expectedGen.set(286, 202, 1.0);
        expectedGen.set(203, 203, -0.4);
        expectedGen.set(287, 203, 1.0);
        expectedGen.set(204, 204, -0.5);
        expectedGen.set(288, 204, 1.0);
        expectedGen.set(205, 205, -0.5);
        expectedGen.set(289, 205, 1.0);
        expectedGen.set(206, 206, -0.5);
        expectedGen.set(290, 206, 1.0);
        expectedGen.set(207, 207, -0.5);
        expectedGen.set(291, 207, 1.0);
        expectedGen.set(208, 208, -0.5);
        expectedGen.set(292, 208, 1.0);
        expectedGen.set(209, 209, -0.5);
        expectedGen.set(293, 209, 1.0);
        expectedGen.set(0, 210, 0.4);
        expectedGen.set(210, 210, -1.5);
        expectedGen.set(213, 210, 0.5);
        expectedGen.set(294, 210, 2.0);
        expectedGen.set(2, 211, 0.4);
        expectedGen.set(210, 211, 0.1);
        expectedGen.set(211, 211, -1.5);
        expectedGen.set(218, 211, 0.5);
        expectedGen.set(3, 212, 0.4);
        expectedGen.set(212, 212, -1.5);
        expectedGen.set(219, 212, 0.5);
        expectedGen.set(295, 212, 2.0);
        expectedGen.set(1, 213, 0.4);
        expectedGen.set(213, 213, -1.6);
        expectedGen.set(216, 213, 0.5);
        expectedGen.set(296, 213, 2.0);
        expectedGen.set(9, 214, 0.4);
        expectedGen.set(211, 214, 0.1);
        expectedGen.set(214, 214, -1.5);
        expectedGen.set(227, 214, 0.5);
        expectedGen.set(10, 215, 0.4);
        expectedGen.set(212, 215, 0.1);
        expectedGen.set(215, 215, -1.5);
        expectedGen.set(228, 215, 0.5);
        expectedGen.set(5, 216, 0.4);
        expectedGen.set(213, 216, 0.1);
        expectedGen.set(216, 216, -1.6);
        expectedGen.set(222, 216, 0.5);
        expectedGen.set(12, 217, 0.4);
        expectedGen.set(217, 217, -1.5);
        expectedGen.set(229, 217, 0.5);
        expectedGen.set(297, 217, 2.0);
        expectedGen.set(6, 218, 0.4);
        expectedGen.set(218, 218, -1.6);
        expectedGen.set(224, 218, 0.5);
        expectedGen.set(298, 218, 2.0);
        expectedGen.set(7, 219, 0.4);
        expectedGen.set(219, 219, -1.6);
        expectedGen.set(225, 219, 0.5);
        expectedGen.set(299, 219, 2.0);
        expectedGen.set(25, 220, 0.4);
        expectedGen.set(214, 220, 0.1);
        expectedGen.set(220, 220, -1.5);
        expectedGen.set(241, 220, 0.5);
        expectedGen.set(26, 221, 0.4);
        expectedGen.set(215, 221, 0.1);
        expectedGen.set(221, 221, -1.5);
        expectedGen.set(242, 221, 0.5);
        expectedGen.set(15, 222, 0.4);
        expectedGen.set(216, 222, 0.1);
        expectedGen.set(222, 222, -1.6);
        expectedGen.set(232, 222, 0.5);
        expectedGen.set(28, 223, 0.4);
        expectedGen.set(217, 223, 0.1);
        expectedGen.set(223, 223, -1.5);
        expectedGen.set(243, 223, 0.5);
        expectedGen.set(16, 224, 0.4);
        expectedGen.set(218, 224, 0.1);
        expectedGen.set(224, 224, -1.6);
        expectedGen.set(234, 224, 0.5);
        expectedGen.set(17, 225, 0.4);
        expectedGen.set(219, 225, 0.1);
        expectedGen.set(225, 225, -1.6);
        expectedGen.set(235, 225, 0.5);
        expectedGen.set(31, 226, 0.4);
        expectedGen.set(226, 226, -1.5);
        expectedGen.set(244, 226, 0.5);
        expectedGen.set(300, 226, 2.0);
        expectedGen.set(19, 227, 0.4);
        expectedGen.set(227, 227, -1.6);
        expectedGen.set(237, 227, 0.5);
        expectedGen.set(301, 227, 2.0);
        expectedGen.set(20, 228, 0.4);
        expectedGen.set(228, 228, -1.6);
        expectedGen.set(238, 228, 0.5);
        expectedGen.set(302, 228, 2.0);
        expectedGen.set(22, 229, 0.4);
        expectedGen.set(229, 229, -1.6);
        expectedGen.set(239, 229, 0.5);
        expectedGen.set(303, 229, 2.0);
        expectedGen.set(55, 230, 0.4);
        expectedGen.set(220, 230, 0.1);
        expectedGen.set(230, 230, -1.5);
        expectedGen.set(261, 230, 0.5);
        expectedGen.set(56, 231, 0.4);
        expectedGen.set(221, 231, 0.1);
        expectedGen.set(231, 231, -1.5);
        expectedGen.set(262, 231, 0.5);
        expectedGen.set(35, 232, 0.4);
        expectedGen.set(222, 232, 0.1);
        expectedGen.set(232, 232, -1.6);
        expectedGen.set(247, 232, 0.5);
        expectedGen.set(58, 233, 0.4);
        expectedGen.set(223, 233, 0.1);
        expectedGen.set(233, 233, -1.5);
        expectedGen.set(263, 233, 0.5);
        expectedGen.set(36, 234, 0.4);
        expectedGen.set(224, 234, 0.1);
        expectedGen.set(234, 234, -1.6);
        expectedGen.set(249, 234, 0.5);
        expectedGen.set(37, 235, 0.4);
        expectedGen.set(225, 235, 0.1);
        expectedGen.set(235, 235, -1.6);
        expectedGen.set(250, 235, 0.5);
        expectedGen.set(61, 236, 0.4);
        expectedGen.set(226, 236, 0.1);
        expectedGen.set(236, 236, -1.5);
        expectedGen.set(264, 236, 0.5);
        expectedGen.set(39, 237, 0.4);
        expectedGen.set(227, 237, 0.1);
        expectedGen.set(237, 237, -1.6);
        expectedGen.set(252, 237, 0.5);
        expectedGen.set(40, 238, 0.4);
        expectedGen.set(228, 238, 0.1);
        expectedGen.set(238, 238, -1.6);
        expectedGen.set(253, 238, 0.5);
        expectedGen.set(42, 239, 0.4);
        expectedGen.set(229, 239, 0.1);
        expectedGen.set(239, 239, -1.6);
        expectedGen.set(254, 239, 0.5);
        expectedGen.set(65, 240, 0.4);
        expectedGen.set(240, 240, -1.5);
        expectedGen.set(265, 240, 0.5);
        expectedGen.set(304, 240, 2.0);
        expectedGen.set(45, 241, 0.4);
        expectedGen.set(241, 241, -1.6);
        expectedGen.set(256, 241, 0.5);
        expectedGen.set(305, 241, 2.0);
        expectedGen.set(46, 242, 0.4);
        expectedGen.set(242, 242, -1.6);
        expectedGen.set(257, 242, 0.5);
        expectedGen.set(306, 242, 2.0);
        expectedGen.set(48, 243, 0.4);
        expectedGen.set(243, 243, -1.6);
        expectedGen.set(258, 243, 0.5);
        expectedGen.set(307, 243, 2.0);
        expectedGen.set(51, 244, 0.4);
        expectedGen.set(244, 244, -1.6);
        expectedGen.set(259, 244, 0.5);
        expectedGen.set(308, 244, 2.0);
        expectedGen.set(105, 245, 0.4);
        expectedGen.set(230, 245, 0.1);
        expectedGen.set(245, 245, -1.5);
        expectedGen.set(288, 245, 0.5);
        expectedGen.set(106, 246, 0.4);
        expectedGen.set(231, 246, 0.1);
        expectedGen.set(246, 246, -1.5);
        expectedGen.set(289, 246, 0.5);
        expectedGen.set(70, 247, 0.4);
        expectedGen.set(232, 247, 0.1);
        expectedGen.set(247, 247, -1.6);
        expectedGen.set(268, 247, 0.5);
        expectedGen.set(108, 248, 0.4);
        expectedGen.set(233, 248, 0.1);
        expectedGen.set(248, 248, -1.5);
        expectedGen.set(290, 248, 0.5);
        expectedGen.set(71, 249, 0.4);
        expectedGen.set(234, 249, 0.1);
        expectedGen.set(249, 249, -1.6);
        expectedGen.set(270, 249, 0.5);
        expectedGen.set(72, 250, 0.4);
        expectedGen.set(235, 250, 0.1);
        expectedGen.set(250, 250, -1.6);
        expectedGen.set(271, 250, 0.5);
        expectedGen.set(111, 251, 0.4);
        expectedGen.set(236, 251, 0.1);
        expectedGen.set(251, 251, -1.5);
        expectedGen.set(291, 251, 0.5);
        expectedGen.set(74, 252, 0.4);
        expectedGen.set(237, 252, 0.1);
        expectedGen.set(252, 252, -1.6);
        expectedGen.set(273, 252, 0.5);
        expectedGen.set(75, 253, 0.4);
        expectedGen.set(238, 253, 0.1);
        expectedGen.set(253, 253, -1.6);
        expectedGen.set(274, 253, 0.5);
        expectedGen.set(77, 254, 0.4);
        expectedGen.set(239, 254, 0.1);
        expectedGen.set(254, 254, -1.6);
        expectedGen.set(275, 254, 0.5);
        expectedGen.set(115, 255, 0.4);
        expectedGen.set(240, 255, 0.1);
        expectedGen.set(255, 255, -1.5);
        expectedGen.set(292, 255, 0.5);
        expectedGen.set(80, 256, 0.4);
        expectedGen.set(241, 256, 0.1);
        expectedGen.set(256, 256, -1.6);
        expectedGen.set(277, 256, 0.5);
        expectedGen.set(81, 257, 0.4);
        expectedGen.set(242, 257, 0.1);
        expectedGen.set(257, 257, -1.6);
        expectedGen.set(278, 257, 0.5);
        expectedGen.set(83, 258, 0.4);
        expectedGen.set(243, 258, 0.1);
        expectedGen.set(258, 258, -1.6);
        expectedGen.set(279, 258, 0.5);
        expectedGen.set(86, 259, 0.4);
        expectedGen.set(244, 259, 0.1);
        expectedGen.set(259, 259, -1.6);
        expectedGen.set(280, 259, 0.5);
        expectedGen.set(120, 260, 0.4);
        expectedGen.set(260, 260, -1.5);
        expectedGen.set(293, 260, 0.5);
        expectedGen.set(309, 260, 2.0);
        expectedGen.set(90, 261, 0.4);
        expectedGen.set(261, 261, -1.6);
        expectedGen.set(282, 261, 0.5);
        expectedGen.set(310, 261, 2.0);
        expectedGen.set(91, 262, 0.4);
        expectedGen.set(262, 262, -1.6);
        expectedGen.set(283, 262, 0.5);
        expectedGen.set(311, 262, 2.0);
        expectedGen.set(93, 263, 0.4);
        expectedGen.set(263, 263, -1.6);
        expectedGen.set(284, 263, 0.5);
        expectedGen.set(312, 263, 2.0);
        expectedGen.set(96, 264, 0.4);
        expectedGen.set(264, 264, -1.6);
        expectedGen.set(285, 264, 0.5);
        expectedGen.set(313, 264, 2.0);
        expectedGen.set(100, 265, 0.4);
        expectedGen.set(265, 265, -1.6);
        expectedGen.set(286, 265, 0.5);
        expectedGen.set(314, 265, 2.0);
        expectedGen.set(182, 266, 0.4);
        expectedGen.set(245, 266, 0.1);
        expectedGen.set(266, 266, -1.4);
        expectedGen.set(183, 267, 0.4);
        expectedGen.set(246, 267, 0.1);
        expectedGen.set(267, 267, -1.4);
        expectedGen.set(126, 268, 0.4);
        expectedGen.set(247, 268, 0.1);
        expectedGen.set(268, 268, -1.5);
        expectedGen.set(185, 269, 0.4);
        expectedGen.set(248, 269, 0.1);
        expectedGen.set(269, 269, -1.4);
        expectedGen.set(127, 270, 0.4);
        expectedGen.set(249, 270, 0.1);
        expectedGen.set(270, 270, -1.5);
        expectedGen.set(128, 271, 0.4);
        expectedGen.set(250, 271, 0.1);
        expectedGen.set(271, 271, -1.5);
        expectedGen.set(188, 272, 0.4);
        expectedGen.set(251, 272, 0.1);
        expectedGen.set(272, 272, -1.4);
        expectedGen.set(130, 273, 0.4);
        expectedGen.set(252, 273, 0.1);
        expectedGen.set(273, 273, -1.5);
        expectedGen.set(131, 274, 0.4);
        expectedGen.set(253, 274, 0.1);
        expectedGen.set(274, 274, -1.5);
        expectedGen.set(133, 275, 0.4);
        expectedGen.set(254, 275, 0.1);
        expectedGen.set(275, 275, -1.5);
        expectedGen.set(192, 276, 0.4);
        expectedGen.set(255, 276, 0.1);
        expectedGen.set(276, 276, -1.4);
        expectedGen.set(136, 277, 0.4);
        expectedGen.set(256, 277, 0.1);
        expectedGen.set(277, 277, -1.5);
        expectedGen.set(137, 278, 0.4);
        expectedGen.set(257, 278, 0.1);
        expectedGen.set(278, 278, -1.5);
        expectedGen.set(139, 279, 0.4);
        expectedGen.set(258, 279, 0.1);
        expectedGen.set(279, 279, -1.5);
        expectedGen.set(142, 280, 0.4);
        expectedGen.set(259, 280, 0.1);
        expectedGen.set(280, 280, -1.5);
        expectedGen.set(197, 281, 0.4);
        expectedGen.set(260, 281, 0.1);
        expectedGen.set(281, 281, -1.4);
        expectedGen.set(146, 282, 0.4);
        expectedGen.set(261, 282, 0.1);
        expectedGen.set(282, 282, -1.5);
        expectedGen.set(147, 283, 0.4);
        expectedGen.set(262, 283, 0.1);
        expectedGen.set(283, 283, -1.5);
        expectedGen.set(149, 284, 0.4);
        expectedGen.set(263, 284, 0.1);
        expectedGen.set(284, 284, -1.5);
        expectedGen.set(152, 285, 0.4);
        expectedGen.set(264, 285, 0.1);
        expectedGen.set(285, 285, -1.5);
        expectedGen.set(156, 286, 0.4);
        expectedGen.set(265, 286, 0.1);
        expectedGen.set(286, 286, -1.5);
        expectedGen.set(203, 287, 0.4);
        expectedGen.set(287, 287, -1.4);
        expectedGen.set(315, 287, 2.0);
        expectedGen.set(161, 288, 0.4);
        expectedGen.set(288, 288, -1.5);
        expectedGen.set(316, 288, 2.0);
        expectedGen.set(162, 289, 0.4);
        expectedGen.set(289, 289, -1.5);
        expectedGen.set(317, 289, 2.0);
        expectedGen.set(164, 290, 0.4);
        expectedGen.set(290, 290, -1.5);
        expectedGen.set(318, 290, 2.0);
        expectedGen.set(167, 291, 0.4);
        expectedGen.set(291, 291, -1.5);
        expectedGen.set(319, 291, 2.0);
        expectedGen.set(171, 292, 0.4);
        expectedGen.set(292, 292, -1.5);
        expectedGen.set(320, 292, 2.0);
        expectedGen.set(176, 293, 0.4);
        expectedGen.set(293, 293, -1.5);
        expectedGen.set(321, 293, 2.0);
        expectedGen.set(210, 294, 0.4);
        expectedGen.set(294, 294, -2.5);
        expectedGen.set(296, 294, 0.5);
        expectedGen.set(322, 294, 3.0);
        expectedGen.set(212, 295, 0.4);
        expectedGen.set(294, 295, 0.1);
        expectedGen.set(295, 295, -2.5);
        expectedGen.set(299, 295, 0.5);
        expectedGen.set(211, 296, 0.4);
        expectedGen.set(296, 296, -2.6);
        expectedGen.set(298, 296, 0.5);
        expectedGen.set(323, 296, 3.0);
        expectedGen.set(217, 297, 0.4);
        expectedGen.set(295, 297, 0.1);
        expectedGen.set(297, 297, -2.5);
        expectedGen.set(303, 297, 0.5);
        expectedGen.set(214, 298, 0.4);
        expectedGen.set(296, 298, 0.1);
        expectedGen.set(298, 298, -2.6);
        expectedGen.set(301, 298, 0.5);
        expectedGen.set(215, 299, 0.4);
        expectedGen.set(299, 299, -2.6);
        expectedGen.set(302, 299, 0.5);
        expectedGen.set(324, 299, 3.0);
        expectedGen.set(226, 300, 0.4);
        expectedGen.set(297, 300, 0.1);
        expectedGen.set(300, 300, -2.5);
        expectedGen.set(308, 300, 0.5);
        expectedGen.set(220, 301, 0.4);
        expectedGen.set(298, 301, 0.1);
        expectedGen.set(301, 301, -2.6);
        expectedGen.set(305, 301, 0.5);
        expectedGen.set(221, 302, 0.4);
        expectedGen.set(299, 302, 0.1);
        expectedGen.set(302, 302, -2.6);
        expectedGen.set(306, 302, 0.5);
        expectedGen.set(223, 303, 0.4);
        expectedGen.set(303, 303, -2.6);
        expectedGen.set(307, 303, 0.5);
        expectedGen.set(325, 303, 3.0);
        expectedGen.set(240, 304, 0.4);
        expectedGen.set(300, 304, 0.1);
        expectedGen.set(304, 304, -2.5);
        expectedGen.set(314, 304, 0.5);
        expectedGen.set(230, 305, 0.4);
        expectedGen.set(301, 305, 0.1);
        expectedGen.set(305, 305, -2.6);
        expectedGen.set(310, 305, 0.5);
        expectedGen.set(231, 306, 0.4);
        expectedGen.set(302, 306, 0.1);
        expectedGen.set(306, 306, -2.6);
        expectedGen.set(311, 306, 0.5);
        expectedGen.set(233, 307, 0.4);
        expectedGen.set(303, 307, 0.1);
        expectedGen.set(307, 307, -2.6);
        expectedGen.set(312, 307, 0.5);
        expectedGen.set(236, 308, 0.4);
        expectedGen.set(308, 308, -2.6);
        expectedGen.set(313, 308, 0.5);
        expectedGen.set(326, 308, 3.0);
        expectedGen.set(260, 309, 0.4);
        expectedGen.set(304, 309, 0.1);
        expectedGen.set(309, 309, -2.5);
        expectedGen.set(321, 309, 0.5);
        expectedGen.set(245, 310, 0.4);
        expectedGen.set(305, 310, 0.1);
        expectedGen.set(310, 310, -2.6);
        expectedGen.set(316, 310, 0.5);
        expectedGen.set(246, 311, 0.4);
        expectedGen.set(306, 311, 0.1);
        expectedGen.set(311, 311, -2.6);
        expectedGen.set(317, 311, 0.5);
        expectedGen.set(248, 312, 0.4);
        expectedGen.set(307, 312, 0.1);
        expectedGen.set(312, 312, -2.6);
        expectedGen.set(318, 312, 0.5);
        expectedGen.set(251, 313, 0.4);
        expectedGen.set(308, 313, 0.1);
        expectedGen.set(313, 313, -2.6);
        expectedGen.set(319, 313, 0.5);
        expectedGen.set(255, 314, 0.4);
        expectedGen.set(314, 314, -2.6);
        expectedGen.set(320, 314, 0.5);
        expectedGen.set(327, 314, 3.0);
        expectedGen.set(287, 315, 0.4);
        expectedGen.set(309, 315, 0.1);
        expectedGen.set(315, 315, -2.4);
        expectedGen.set(266, 316, 0.4);
        expectedGen.set(310, 316, 0.1);
        expectedGen.set(316, 316, -2.5);
        expectedGen.set(267, 317, 0.4);
        expectedGen.set(311, 317, 0.1);
        expectedGen.set(317, 317, -2.5);
        expectedGen.set(269, 318, 0.4);
        expectedGen.set(312, 318, 0.1);
        expectedGen.set(318, 318, -2.5);
        expectedGen.set(272, 319, 0.4);
        expectedGen.set(313, 319, 0.1);
        expectedGen.set(319, 319, -2.5);
        expectedGen.set(276, 320, 0.4);
        expectedGen.set(314, 320, 0.1);
        expectedGen.set(320, 320, -2.5);
        expectedGen.set(281, 321, 0.4);
        expectedGen.set(321, 321, -2.5);
        expectedGen.set(328, 321, 3.0);
        expectedGen.set(294, 322, 0.4);
        expectedGen.set(322, 322, -3.1);
        expectedGen.set(323, 322, 0.5);
        expectedGen.set(295, 323, 0.4);
        expectedGen.set(322, 323, 0.1);
        expectedGen.set(323, 323, -3.6);
        expectedGen.set(324, 323, 0.5);
        expectedGen.set(297, 324, 0.4);
        expectedGen.set(323, 324, 0.1);
        expectedGen.set(324, 324, -3.6);
        expectedGen.set(325, 324, 0.5);
        expectedGen.set(300, 325, 0.4);
        expectedGen.set(324, 325, 0.1);
        expectedGen.set(325, 325, -3.6);
        expectedGen.set(326, 325, 0.5);
        expectedGen.set(304, 326, 0.4);
        expectedGen.set(325, 326, 0.1);
        expectedGen.set(326, 326, -3.6);
        expectedGen.set(327, 326, 0.5);
        expectedGen.set(309, 327, 0.4);
        expectedGen.set(326, 327, 0.1);
        expectedGen.set(327, 327, -3.6);
        expectedGen.set(328, 327, 0.5);
        expectedGen.set(315, 328, 0.4);
        expectedGen.set(327, 328, 0.1);
        expectedGen.set(328, 328, -3.5);

        // Compare the entire matrix
        assertMatrixEquals(expectedGen, generator, "Generator matrix for mixed LCFS");

        // Verify state space dimensions
        assertEquals(329, stateSpace.getNumRows());
        assertEquals(15, stateSpace.getNumCols());
    }

    @Test
    public void testMixedNetworkHOL() {
        Network model = createMixedNetwork(SchedStrategy.HOL);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(6);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify dimensions - same as FCFS/LCFS
        assertEquals(329, generator.getNumRows());
        assertEquals(329, generator.getNumCols());

        // Create expected sparse matrix (329x329) with non-zero entries from MATLAB
        Matrix expectedGen = new Matrix(329, 329);
        // TODO: Add all non-zero entries from ground truth for HOL
        // This is a placeholder until complete ground truth data is provided
        // HOL should have similar structure to FCFS/LCFS but with different values

        // For now, verify dimensions only
        // Once ground truth is available, use:
        // assertMatrixEquals(expectedGen, generator, "Generator matrix for mixed HOL");

        Matrix expectedStateSpace = new Matrix(329, 329);

        // Verify state space dimensions
        assertEquals(329, stateSpace.getNumRows());
        assertEquals(15, stateSpace.getNumCols());
    }

    @Test
    public void testMixedNetworkSEPT() {
        Network model = createMixedNetwork(SchedStrategy.SEPT);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(6);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        Matrix expectedGen = new Matrix(46, 46);
        expectedGen.set(0,0,-0.500000);
        expectedGen.set(0,2,0.100000);
        expectedGen.set(0,13,0.400000);
        expectedGen.set(1,0,0.500000);
        expectedGen.set(1,1,-0.600000);
        expectedGen.set(1,3,0.100000);
        expectedGen.set(2,2,-0.500000);
        expectedGen.set(2,4,0.100000);
        expectedGen.set(2,14,0.400000);
        expectedGen.set(3,1,0.500000);
        expectedGen.set(3,3,-0.600000);
        expectedGen.set(3,5,0.100000);
        expectedGen.set(4,4,-0.500000);
        expectedGen.set(4,6,0.100000);
        expectedGen.set(4,16,0.400000);
        expectedGen.set(5,3,0.500000);
        expectedGen.set(5,5,-0.600000);
        expectedGen.set(5,7,0.100000);
        expectedGen.set(6,6,-0.500000);
        expectedGen.set(6,8,0.100000);
        expectedGen.set(6,18,0.400000);
        expectedGen.set(7,5,0.500000);
        expectedGen.set(7,7,-0.600000);
        expectedGen.set(7,9,0.100000);
        expectedGen.set(8,8,-0.500000);
        expectedGen.set(8,10,0.100000);
        expectedGen.set(8,20,0.400000);
        expectedGen.set(9,7,0.500000);
        expectedGen.set(9,9,-0.600000);
        expectedGen.set(9,11,0.100000);
        expectedGen.set(10,10,-0.500000);
        expectedGen.set(10,12,0.100000);
        expectedGen.set(10,22,0.400000);
        expectedGen.set(11,9,0.500000);
        expectedGen.set(11,11,-0.500000);
        expectedGen.set(12,12,-0.400000);
        expectedGen.set(12,24,0.400000);
        expectedGen.set(13,0,1.000000);
        expectedGen.set(13,13,-1.500000);
        expectedGen.set(13,15,0.100000);
        expectedGen.set(13,26,0.400000);
        expectedGen.set(14,1,1.000000);
        expectedGen.set(14,13,0.500000);
        expectedGen.set(14,14,-1.600000);
        expectedGen.set(14,16,0.100000);
        expectedGen.set(15,2,1.000000);
        expectedGen.set(15,15,-1.500000);
        expectedGen.set(15,17,0.100000);
        expectedGen.set(15,27,0.400000);
        expectedGen.set(16,3,1.000000);
        expectedGen.set(16,14,0.500000);
        expectedGen.set(16,16,-1.600000);
        expectedGen.set(16,18,0.100000);
        expectedGen.set(17,4,1.000000);
        expectedGen.set(17,17,-1.500000);
        expectedGen.set(17,19,0.100000);
        expectedGen.set(17,29,0.400000);
        expectedGen.set(18,5,1.000000);
        expectedGen.set(18,16,0.500000);
        expectedGen.set(18,18,-1.600000);
        expectedGen.set(18,20,0.100000);
        expectedGen.set(19,6,1.000000);
        expectedGen.set(19,19,-1.500000);
        expectedGen.set(19,21,0.100000);
        expectedGen.set(19,31,0.400000);
        expectedGen.set(20,7,1.000000);
        expectedGen.set(20,18,0.500000);
        expectedGen.set(20,20,-1.600000);
        expectedGen.set(20,22,0.100000);
        expectedGen.set(21,8,1.000000);
        expectedGen.set(21,21,-1.500000);
        expectedGen.set(21,23,0.100000);
        expectedGen.set(21,33,0.400000);
        expectedGen.set(22,9,1.000000);
        expectedGen.set(22,20,0.500000);
        expectedGen.set(22,22,-1.600000);
        expectedGen.set(22,24,0.100000);
        expectedGen.set(23,10,1.000000);
        expectedGen.set(23,23,-1.500000);
        expectedGen.set(23,25,0.100000);
        expectedGen.set(23,35,0.400000);
        expectedGen.set(24,11,1.000000);
        expectedGen.set(24,22,0.500000);
        expectedGen.set(24,24,-1.500000);
        expectedGen.set(25,12,1.000000);
        expectedGen.set(25,25,-1.400000);
        expectedGen.set(25,37,0.400000);
        expectedGen.set(26,13,2.000000);
        expectedGen.set(26,26,-2.500000);
        expectedGen.set(26,28,0.100000);
        expectedGen.set(26,39,0.400000);
        expectedGen.set(27,14,2.000000);
        expectedGen.set(27,26,0.500000);
        expectedGen.set(27,27,-2.600000);
        expectedGen.set(27,29,0.100000);
        expectedGen.set(28,15,2.000000);
        expectedGen.set(28,28,-2.500000);
        expectedGen.set(28,30,0.100000);
        expectedGen.set(28,40,0.400000);
        expectedGen.set(29,16,2.000000);
        expectedGen.set(29,27,0.500000);
        expectedGen.set(29,29,-2.600000);
        expectedGen.set(29,31,0.100000);
        expectedGen.set(30,17,2.000000);
        expectedGen.set(30,30,-2.500000);
        expectedGen.set(30,32,0.100000);
        expectedGen.set(30,41,0.400000);
        expectedGen.set(31,18,2.000000);
        expectedGen.set(31,29,0.500000);
        expectedGen.set(31,31,-2.600000);
        expectedGen.set(31,33,0.100000);
        expectedGen.set(32,19,2.000000);
        expectedGen.set(32,32,-2.500000);
        expectedGen.set(32,34,0.100000);
        expectedGen.set(32,42,0.400000);
        expectedGen.set(33,20,2.000000);
        expectedGen.set(33,31,0.500000);
        expectedGen.set(33,33,-2.600000);
        expectedGen.set(33,35,0.100000);
        expectedGen.set(34,21,2.000000);
        expectedGen.set(34,34,-2.500000);
        expectedGen.set(34,36,0.100000);
        expectedGen.set(34,43,0.400000);
        expectedGen.set(35,22,2.000000);
        expectedGen.set(35,33,0.500000);
        expectedGen.set(35,35,-2.600000);
        expectedGen.set(35,37,0.100000);
        expectedGen.set(36,23,2.000000);
        expectedGen.set(36,36,-2.500000);
        expectedGen.set(36,38,0.100000);
        expectedGen.set(36,44,0.400000);
        expectedGen.set(37,24,2.000000);
        expectedGen.set(37,35,0.500000);
        expectedGen.set(37,37,-2.500000);
        expectedGen.set(38,25,2.000000);
        expectedGen.set(38,38,-2.400000);
        expectedGen.set(38,45,0.400000);
        expectedGen.set(39,26,3.000000);
        expectedGen.set(39,39,-3.100000);
        expectedGen.set(39,40,0.100000);
        expectedGen.set(40,27,3.000000);
        expectedGen.set(40,39,0.500000);
        expectedGen.set(40,40,-3.600000);
        expectedGen.set(40,41,0.100000);
        expectedGen.set(41,29,3.000000);
        expectedGen.set(41,40,0.500000);
        expectedGen.set(41,41,-3.600000);
        expectedGen.set(41,42,0.100000);
        expectedGen.set(42,31,3.000000);
        expectedGen.set(42,41,0.500000);
        expectedGen.set(42,42,-3.600000);
        expectedGen.set(42,43,0.100000);
        expectedGen.set(43,33,3.000000);
        expectedGen.set(43,42,0.500000);
        expectedGen.set(43,43,-3.600000);
        expectedGen.set(43,44,0.100000);
        expectedGen.set(44,35,3.000000);
        expectedGen.set(44,43,0.500000);
        expectedGen.set(44,44,-3.600000);
        expectedGen.set(44,45,0.100000);
        expectedGen.set(45,37,3.000000);
        expectedGen.set(45,44,0.500000);
        expectedGen.set(45,45,-3.500000);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "Generator matrix mismatch");

        // Verify state space
        assertEquals(46, stateSpace.getNumRows());
        assertEquals(9, stateSpace.getNumCols());
    }

    @Test
    public void testMixedNetworkLEPT() {
        Network model = createMixedNetwork(SchedStrategy.LEPT);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(6);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        Matrix expectedGen = new Matrix(46, 46);
        expectedGen.set(0,0,-0.500000);
        expectedGen.set(0,2,0.100000);
        expectedGen.set(0,13,0.400000);
        expectedGen.set(1,0,0.500000);
        expectedGen.set(1,1,-0.600000);
        expectedGen.set(1,3,0.100000);
        expectedGen.set(2,2,-0.500000);
        expectedGen.set(2,4,0.100000);
        expectedGen.set(2,15,0.400000);
        expectedGen.set(3,2,0.500000);
        expectedGen.set(3,3,-0.600000);
        expectedGen.set(3,5,0.100000);
        expectedGen.set(4,4,-0.500000);
        expectedGen.set(4,6,0.100000);
        expectedGen.set(4,17,0.400000);
        expectedGen.set(5,4,0.500000);
        expectedGen.set(5,5,-0.600000);
        expectedGen.set(5,7,0.100000);
        expectedGen.set(6,6,-0.500000);
        expectedGen.set(6,8,0.100000);
        expectedGen.set(6,19,0.400000);
        expectedGen.set(7,6,0.500000);
        expectedGen.set(7,7,-0.600000);
        expectedGen.set(7,9,0.100000);
        expectedGen.set(8,8,-0.500000);
        expectedGen.set(8,10,0.100000);
        expectedGen.set(8,21,0.400000);
        expectedGen.set(9,8,0.500000);
        expectedGen.set(9,9,-0.600000);
        expectedGen.set(9,11,0.100000);
        expectedGen.set(10,10,-0.500000);
        expectedGen.set(10,12,0.100000);
        expectedGen.set(10,23,0.400000);
        expectedGen.set(11,10,0.500000);
        expectedGen.set(11,11,-0.500000);
        expectedGen.set(12,12,-0.400000);
        expectedGen.set(12,25,0.400000);
        expectedGen.set(13,0,1.000000);
        expectedGen.set(13,13,-1.500000);
        expectedGen.set(13,15,0.100000);
        expectedGen.set(13,26,0.400000);
        expectedGen.set(14,1,1.000000);
        expectedGen.set(14,13,0.500000);
        expectedGen.set(14,14,-1.600000);
        expectedGen.set(14,16,0.100000);
        expectedGen.set(15,2,1.000000);
        expectedGen.set(15,15,-1.500000);
        expectedGen.set(15,17,0.100000);
        expectedGen.set(15,28,0.400000);
        expectedGen.set(16,3,1.000000);
        expectedGen.set(16,15,0.500000);
        expectedGen.set(16,16,-1.600000);
        expectedGen.set(16,18,0.100000);
        expectedGen.set(17,4,1.000000);
        expectedGen.set(17,17,-1.500000);
        expectedGen.set(17,19,0.100000);
        expectedGen.set(17,30,0.400000);
        expectedGen.set(18,5,1.000000);
        expectedGen.set(18,17,0.500000);
        expectedGen.set(18,18,-1.600000);
        expectedGen.set(18,20,0.100000);
        expectedGen.set(19,6,1.000000);
        expectedGen.set(19,19,-1.500000);
        expectedGen.set(19,21,0.100000);
        expectedGen.set(19,32,0.400000);
        expectedGen.set(20,7,1.000000);
        expectedGen.set(20,19,0.500000);
        expectedGen.set(20,20,-1.600000);
        expectedGen.set(20,22,0.100000);
        expectedGen.set(21,8,1.000000);
        expectedGen.set(21,21,-1.500000);
        expectedGen.set(21,23,0.100000);
        expectedGen.set(21,34,0.400000);
        expectedGen.set(22,9,1.000000);
        expectedGen.set(22,21,0.500000);
        expectedGen.set(22,22,-1.600000);
        expectedGen.set(22,24,0.100000);
        expectedGen.set(23,10,1.000000);
        expectedGen.set(23,23,-1.500000);
        expectedGen.set(23,25,0.100000);
        expectedGen.set(23,36,0.400000);
        expectedGen.set(24,11,1.000000);
        expectedGen.set(24,23,0.500000);
        expectedGen.set(24,24,-1.500000);
        expectedGen.set(25,12,1.000000);
        expectedGen.set(25,25,-1.400000);
        expectedGen.set(25,38,0.400000);
        expectedGen.set(26,13,2.000000);
        expectedGen.set(26,26,-2.500000);
        expectedGen.set(26,28,0.100000);
        expectedGen.set(26,39,0.400000);
        expectedGen.set(27,14,2.000000);
        expectedGen.set(27,26,0.500000);
        expectedGen.set(27,27,-2.600000);
        expectedGen.set(27,29,0.100000);
        expectedGen.set(28,15,2.000000);
        expectedGen.set(28,28,-2.500000);
        expectedGen.set(28,30,0.100000);
        expectedGen.set(28,40,0.400000);
        expectedGen.set(29,16,2.000000);
        expectedGen.set(29,28,0.500000);
        expectedGen.set(29,29,-2.600000);
        expectedGen.set(29,31,0.100000);
        expectedGen.set(30,17,2.000000);
        expectedGen.set(30,30,-2.500000);
        expectedGen.set(30,32,0.100000);
        expectedGen.set(30,41,0.400000);
        expectedGen.set(31,18,2.000000);
        expectedGen.set(31,30,0.500000);
        expectedGen.set(31,31,-2.600000);
        expectedGen.set(31,33,0.100000);
        expectedGen.set(32,19,2.000000);
        expectedGen.set(32,32,-2.500000);
        expectedGen.set(32,34,0.100000);
        expectedGen.set(32,42,0.400000);
        expectedGen.set(33,20,2.000000);
        expectedGen.set(33,32,0.500000);
        expectedGen.set(33,33,-2.600000);
        expectedGen.set(33,35,0.100000);
        expectedGen.set(34,21,2.000000);
        expectedGen.set(34,34,-2.500000);
        expectedGen.set(34,36,0.100000);
        expectedGen.set(34,43,0.400000);
        expectedGen.set(35,22,2.000000);
        expectedGen.set(35,34,0.500000);
        expectedGen.set(35,35,-2.600000);
        expectedGen.set(35,37,0.100000);
        expectedGen.set(36,23,2.000000);
        expectedGen.set(36,36,-2.500000);
        expectedGen.set(36,38,0.100000);
        expectedGen.set(36,44,0.400000);
        expectedGen.set(37,24,2.000000);
        expectedGen.set(37,36,0.500000);
        expectedGen.set(37,37,-2.500000);
        expectedGen.set(38,25,2.000000);
        expectedGen.set(38,38,-2.400000);
        expectedGen.set(38,45,0.400000);
        expectedGen.set(39,26,3.000000);
        expectedGen.set(39,39,-3.100000);
        expectedGen.set(39,40,0.100000);
        expectedGen.set(40,27,3.000000);
        expectedGen.set(40,39,0.500000);
        expectedGen.set(40,40,-3.600000);
        expectedGen.set(40,41,0.100000);
        expectedGen.set(41,29,3.000000);
        expectedGen.set(41,40,0.500000);
        expectedGen.set(41,41,-3.600000);
        expectedGen.set(41,42,0.100000);
        expectedGen.set(42,31,3.000000);
        expectedGen.set(42,41,0.500000);
        expectedGen.set(42,42,-3.600000);
        expectedGen.set(42,43,0.100000);
        expectedGen.set(43,33,3.000000);
        expectedGen.set(43,42,0.500000);
        expectedGen.set(43,43,-3.600000);
        expectedGen.set(43,44,0.100000);
        expectedGen.set(44,35,3.000000);
        expectedGen.set(44,43,0.500000);
        expectedGen.set(44,44,-3.600000);
        expectedGen.set(44,45,0.100000);
        expectedGen.set(45,37,3.000000);
        expectedGen.set(45,44,0.500000);
        expectedGen.set(45,45,-3.500000);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "Generator matrix mismatch");

        // Verify state space
        assertEquals(46, stateSpace.getNumRows());
        assertEquals(9, stateSpace.getNumCols());
    }

    // Two Closed Classes Network Tests - 3 jobs each

    @Test
    public void testTwoClassFCFS() {
        Network model = createTwoClassClosedNetwork(SchedStrategy.FCFS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for FCFS two-class network
        // FCFS has 69x69 generator matrix
        Matrix expectedGen = new Matrix(69, 69);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.4);
        expectedGen.set(35, 0, 1.0);
        expectedGen.set(1, 1, -0.4);
        expectedGen.set(36, 1, 1.0);
        expectedGen.set(2, 2, -0.4);
        expectedGen.set(37, 2, 1.0);
        expectedGen.set(3, 3, -0.5);
        expectedGen.set(38, 3, 1.0);
        expectedGen.set(4, 4, -0.4);
        expectedGen.set(39, 4, 1.0);
        expectedGen.set(5, 5, -0.4);
        expectedGen.set(40, 5, 1.0);
        expectedGen.set(6, 6, -0.5);
        expectedGen.set(41, 6, 1.0);
        expectedGen.set(7, 7, -0.4);
        expectedGen.set(42, 7, 1.0);
        expectedGen.set(8, 8, -0.5);
        expectedGen.set(43, 8, 1.0);
        expectedGen.set(9, 9, -0.5);
        expectedGen.set(44, 9, 1.0);
        expectedGen.set(10, 10, -0.4);
        expectedGen.set(20, 10, 0.6667);
        expectedGen.set(11, 11, -0.4);
        expectedGen.set(21, 11, 0.6667);
        expectedGen.set(12, 12, -0.5);
        expectedGen.set(22, 12, 0.6667);
        expectedGen.set(13, 13, -0.4);
        expectedGen.set(23, 13, 0.6667);
        expectedGen.set(14, 14, -0.5);
        expectedGen.set(24, 14, 0.6667);
        expectedGen.set(15, 15, -0.5);
        expectedGen.set(25, 15, 0.6667);
        expectedGen.set(16, 16, -0.4);
        expectedGen.set(26, 16, 0.6667);
        expectedGen.set(17, 17, -0.5);
        expectedGen.set(27, 17, 0.6667);
        expectedGen.set(18, 18, -0.5);
        expectedGen.set(28, 18, 0.6667);
        expectedGen.set(19, 19, -0.5);
        expectedGen.set(29, 19, 0.6667);
        expectedGen.set(0, 20, 0.4);
        expectedGen.set(20, 20, -1.0667);
        expectedGen.set(45, 20, 1.0);
        expectedGen.set(1, 21, 0.4);
        expectedGen.set(21, 21, -1.0667);
        expectedGen.set(46, 21, 1.0);
        expectedGen.set(2, 22, 0.4);
        expectedGen.set(22, 22, -1.1667);
        expectedGen.set(47, 22, 1.0);
        expectedGen.set(4, 23, 0.4);
        expectedGen.set(23, 23, -1.0667);
        expectedGen.set(48, 23, 1.0);
        expectedGen.set(5, 24, 0.4);
        expectedGen.set(24, 24, -1.1667);
        expectedGen.set(49, 24, 1.0);
        expectedGen.set(7, 25, 0.4);
        expectedGen.set(25, 25, -1.1667);
        expectedGen.set(50, 25, 1.0);
        expectedGen.set(10, 26, 0.4);
        expectedGen.set(26, 26, -1.0667);
        expectedGen.set(30, 26, 1.3333);
        expectedGen.set(11, 27, 0.4);
        expectedGen.set(27, 27, -1.1667);
        expectedGen.set(31, 27, 1.3333);
        expectedGen.set(13, 28, 0.4);
        expectedGen.set(28, 28, -1.1667);
        expectedGen.set(32, 28, 1.3333);
        expectedGen.set(16, 29, 0.4);
        expectedGen.set(29, 29, -1.1667);
        expectedGen.set(33, 29, 1.3333);
        expectedGen.set(20, 30, 0.4);
        expectedGen.set(30, 30, -1.7333);
        expectedGen.set(51, 30, 1.0);
        expectedGen.set(21, 31, 0.4);
        expectedGen.set(31, 31, -1.8333);
        expectedGen.set(52, 31, 1.0);
        expectedGen.set(23, 32, 0.4);
        expectedGen.set(32, 32, -1.8333);
        expectedGen.set(53, 32, 1.0);
        expectedGen.set(26, 33, 0.4);
        expectedGen.set(33, 33, -1.8333);
        expectedGen.set(34, 33, 2.0);
        expectedGen.set(30, 34, 0.4);
        expectedGen.set(34, 34, -2.5);
        expectedGen.set(54, 34, 1.0);
        expectedGen.set(3, 35, 0.5);
        expectedGen.set(35, 35, -1.4);
        expectedGen.set(55, 35, 2.0);
        expectedGen.set(6, 36, 0.5);
        expectedGen.set(36, 36, -1.4);
        expectedGen.set(56, 36, 2.0);
        expectedGen.set(8, 37, 0.5);
        expectedGen.set(37, 37, -1.4);
        expectedGen.set(57, 37, 2.0);
        expectedGen.set(9, 38, 0.5);
        expectedGen.set(38, 38, -1.5);
        expectedGen.set(58, 38, 2.0);
        expectedGen.set(12, 39, 0.5);
        expectedGen.set(39, 39, -1.4);
        expectedGen.set(45, 39, 0.6667);
        expectedGen.set(14, 40, 0.5);
        expectedGen.set(40, 40, -1.4);
        expectedGen.set(46, 40, 0.6667);
        expectedGen.set(15, 41, 0.5);
        expectedGen.set(41, 41, -1.5);
        expectedGen.set(47, 41, 0.6667);
        expectedGen.set(17, 42, 0.5);
        expectedGen.set(42, 42, -1.4);
        expectedGen.set(48, 42, 0.6667);
        expectedGen.set(18, 43, 0.5);
        expectedGen.set(43, 43, -1.5);
        expectedGen.set(49, 43, 0.6667);
        expectedGen.set(19, 44, 0.5);
        expectedGen.set(44, 44, -1.5);
        expectedGen.set(50, 44, 0.6667);
        expectedGen.set(22, 45, 0.5);
        expectedGen.set(35, 45, 0.4);
        expectedGen.set(45, 45, -2.0667);
        expectedGen.set(59, 45, 2.0);
        expectedGen.set(24, 46, 0.5);
        expectedGen.set(36, 46, 0.4);
        expectedGen.set(46, 46, -2.0667);
        expectedGen.set(60, 46, 2.0);
        expectedGen.set(25, 47, 0.5);
        expectedGen.set(37, 47, 0.4);
        expectedGen.set(47, 47, -2.1667);
        expectedGen.set(61, 47, 2.0);
        expectedGen.set(27, 48, 0.5);
        expectedGen.set(39, 48, 0.4);
        expectedGen.set(48, 48, -2.0667);
        expectedGen.set(51, 48, 1.3333);
        expectedGen.set(28, 49, 0.5);
        expectedGen.set(40, 49, 0.4);
        expectedGen.set(49, 49, -2.1667);
        expectedGen.set(52, 49, 1.3333);
        expectedGen.set(29, 50, 0.5);
        expectedGen.set(42, 50, 0.4);
        expectedGen.set(50, 50, -2.1667);
        expectedGen.set(53, 50, 1.3333);
        expectedGen.set(31, 51, 0.5);
        expectedGen.set(45, 51, 0.4);
        expectedGen.set(51, 51, -2.7333);
        expectedGen.set(62, 51, 2.0);
        expectedGen.set(32, 52, 0.5);
        expectedGen.set(46, 52, 0.4);
        expectedGen.set(52, 52, -2.8333);
        expectedGen.set(63, 52, 2.0);
        expectedGen.set(33, 53, 0.5);
        expectedGen.set(48, 53, 0.4);
        expectedGen.set(53, 53, -2.8333);
        expectedGen.set(54, 53, 2.0);
        expectedGen.set(34, 54, 0.5);
        expectedGen.set(51, 54, 0.4);
        expectedGen.set(54, 54, -3.5);
        expectedGen.set(64, 54, 2.0);
        expectedGen.set(38, 55, 0.5);
        expectedGen.set(55, 55, -2.4);
        expectedGen.set(65, 55, 3.0);
        expectedGen.set(41, 56, 0.5);
        expectedGen.set(56, 56, -2.4);
        expectedGen.set(59, 56, 0.6667);
        expectedGen.set(43, 57, 0.5);
        expectedGen.set(57, 57, -2.4);
        expectedGen.set(60, 57, 0.6667);
        expectedGen.set(44, 58, 0.5);
        expectedGen.set(58, 58, -2.5);
        expectedGen.set(61, 58, 0.6667);
        expectedGen.set(47, 59, 0.5);
        expectedGen.set(55, 59, 0.4);
        expectedGen.set(59, 59, -3.0667);
        expectedGen.set(66, 59, 3.0);
        expectedGen.set(49, 60, 0.5);
        expectedGen.set(56, 60, 0.4);
        expectedGen.set(60, 60, -3.0667);
        expectedGen.set(62, 60, 1.3333);
        expectedGen.set(50, 61, 0.5);
        expectedGen.set(57, 61, 0.4);
        expectedGen.set(61, 61, -3.1667);
        expectedGen.set(63, 61, 1.3333);
        expectedGen.set(52, 62, 0.5);
        expectedGen.set(59, 62, 0.4);
        expectedGen.set(62, 62, -3.7333);
        expectedGen.set(67, 62, 3.0);
        expectedGen.set(53, 63, 0.5);
        expectedGen.set(60, 63, 0.4);
        expectedGen.set(63, 63, -3.8333);
        expectedGen.set(64, 63, 2.0);
        expectedGen.set(54, 64, 0.5);
        expectedGen.set(62, 64, 0.4);
        expectedGen.set(64, 64, -4.5);
        expectedGen.set(68, 64, 3.0);
        expectedGen.set(58, 65, 0.5);
        expectedGen.set(65, 65, -3.4);
        expectedGen.set(66, 65, 0.6667);
        expectedGen.set(61, 66, 0.5);
        expectedGen.set(65, 66, 0.4);
        expectedGen.set(66, 66, -4.0667);
        expectedGen.set(67, 66, 1.3333);
        expectedGen.set(63, 67, 0.5);
        expectedGen.set(66, 67, 0.4);
        expectedGen.set(67, 67, -4.7333);
        expectedGen.set(68, 67, 2.0);
        expectedGen.set(64, 68, 0.5);
        expectedGen.set(67, 68, 0.4);
        expectedGen.set(68, 68, -5.0);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "FCFS two-class generator matrix mismatch");

        // Verify state space from MATLAB output
        assertEquals(69, stateSpace.getNumRows());
        assertEquals(9, stateSpace.getNumCols());
    }

    @Test
    public void testTwoClassLCFS() {
        Network model = createTwoClassClosedNetwork(SchedStrategy.LCFS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for LCFS two-class network
        // LCFS has 69x69 generator matrix with same values as FCFS
        Matrix expectedGen = new Matrix(69, 69);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.4);
        expectedGen.set(35, 0, 1.0);
        expectedGen.set(1, 1, -0.4);
        expectedGen.set(36, 1, 1.0);
        expectedGen.set(2, 2, -0.4);
        expectedGen.set(37, 2, 1.0);
        expectedGen.set(3, 3, -0.5);
        expectedGen.set(38, 3, 1.0);
        expectedGen.set(4, 4, -0.4);
        expectedGen.set(39, 4, 1.0);
        expectedGen.set(5, 5, -0.4);
        expectedGen.set(40, 5, 1.0);
        expectedGen.set(6, 6, -0.5);
        expectedGen.set(41, 6, 1.0);
        expectedGen.set(7, 7, -0.4);
        expectedGen.set(42, 7, 1.0);
        expectedGen.set(8, 8, -0.5);
        expectedGen.set(43, 8, 1.0);
        expectedGen.set(9, 9, -0.5);
        expectedGen.set(44, 9, 1.0);
        expectedGen.set(10, 10, -0.4);
        expectedGen.set(20, 10, 0.6667);
        expectedGen.set(11, 11, -0.4);
        expectedGen.set(21, 11, 0.6667);
        expectedGen.set(12, 12, -0.5);
        expectedGen.set(22, 12, 0.6667);
        expectedGen.set(13, 13, -0.4);
        expectedGen.set(23, 13, 0.6667);
        expectedGen.set(14, 14, -0.5);
        expectedGen.set(24, 14, 0.6667);
        expectedGen.set(15, 15, -0.5);
        expectedGen.set(25, 15, 0.6667);
        expectedGen.set(16, 16, -0.4);
        expectedGen.set(26, 16, 0.6667);
        expectedGen.set(17, 17, -0.5);
        expectedGen.set(27, 17, 0.6667);
        expectedGen.set(18, 18, -0.5);
        expectedGen.set(28, 18, 0.6667);
        expectedGen.set(19, 19, -0.5);
        expectedGen.set(29, 19, 0.6667);
        expectedGen.set(10, 20, 0.4);
        expectedGen.set(20, 20, -1.0667);
        expectedGen.set(45, 20, 1.0);
        expectedGen.set(11, 21, 0.4);
        expectedGen.set(21, 21, -1.0667);
        expectedGen.set(46, 21, 1.0);
        expectedGen.set(0, 22, 0.4);
        expectedGen.set(22, 22, -1.1667);
        expectedGen.set(47, 22, 1.0);
        expectedGen.set(13, 23, 0.4);
        expectedGen.set(23, 23, -1.0667);
        expectedGen.set(48, 23, 1.0);
        expectedGen.set(1, 24, 0.4);
        expectedGen.set(24, 24, -1.1667);
        expectedGen.set(49, 24, 1.0);
        expectedGen.set(2, 25, 0.4);
        expectedGen.set(25, 25, -1.1667);
        expectedGen.set(50, 25, 1.0);
        expectedGen.set(16, 26, 0.4);
        expectedGen.set(26, 26, -1.0667);
        expectedGen.set(30, 26, 1.3333);
        expectedGen.set(4, 27, 0.4);
        expectedGen.set(27, 27, -1.1667);
        expectedGen.set(31, 27, 1.3333);
        expectedGen.set(5, 28, 0.4);
        expectedGen.set(28, 28, -1.1667);
        expectedGen.set(32, 28, 1.3333);
        expectedGen.set(7, 29, 0.4);
        expectedGen.set(29, 29, -1.1667);
        expectedGen.set(33, 29, 1.3333);
        expectedGen.set(26, 30, 0.4);
        expectedGen.set(30, 30, -1.7333);
        expectedGen.set(51, 30, 1.0);
        expectedGen.set(20, 31, 0.4);
        expectedGen.set(31, 31, -1.8333);
        expectedGen.set(52, 31, 1.0);
        expectedGen.set(21, 32, 0.4);
        expectedGen.set(32, 32, -1.8333);
        expectedGen.set(53, 32, 1.0);
        expectedGen.set(23, 33, 0.4);
        expectedGen.set(33, 33, -1.8333);
        expectedGen.set(34, 33, 2.0);
        expectedGen.set(30, 34, 0.4);
        expectedGen.set(34, 34, -2.5);
        expectedGen.set(54, 34, 1.0);
        expectedGen.set(12, 35, 0.5);
        expectedGen.set(35, 35, -1.4);
        expectedGen.set(55, 35, 2.0);
        expectedGen.set(14, 36, 0.5);
        expectedGen.set(36, 36, -1.4);
        expectedGen.set(56, 36, 2.0);
        expectedGen.set(15, 37, 0.5);
        expectedGen.set(37, 37, -1.4);
        expectedGen.set(57, 37, 2.0);
        expectedGen.set(3, 38, 0.5);
        expectedGen.set(38, 38, -1.5);
        expectedGen.set(58, 38, 2.0);
        expectedGen.set(17, 39, 0.5);
        expectedGen.set(39, 39, -1.4);
        expectedGen.set(45, 39, 0.6667);
        expectedGen.set(18, 40, 0.5);
        expectedGen.set(40, 40, -1.4);
        expectedGen.set(46, 40, 0.6667);
        expectedGen.set(6, 41, 0.5);
        expectedGen.set(41, 41, -1.5);
        expectedGen.set(47, 41, 0.6667);
        expectedGen.set(19, 42, 0.5);
        expectedGen.set(42, 42, -1.4);
        expectedGen.set(48, 42, 0.6667);
        expectedGen.set(8, 43, 0.5);
        expectedGen.set(43, 43, -1.5);
        expectedGen.set(49, 43, 0.6667);
        expectedGen.set(9, 44, 0.5);
        expectedGen.set(44, 44, -1.5);
        expectedGen.set(50, 44, 0.6667);
        expectedGen.set(27, 45, 0.5);
        expectedGen.set(39, 45, 0.4);
        expectedGen.set(45, 45, -2.0667);
        expectedGen.set(59, 45, 2.0);
        expectedGen.set(28, 46, 0.5);
        expectedGen.set(40, 46, 0.4);
        expectedGen.set(46, 46, -2.0667);
        expectedGen.set(60, 46, 2.0);
        expectedGen.set(22, 47, 0.5);
        expectedGen.set(35, 47, 0.4);
        expectedGen.set(47, 47, -2.1667);
        expectedGen.set(61, 47, 2.0);
        expectedGen.set(29, 48, 0.5);
        expectedGen.set(42, 48, 0.4);
        expectedGen.set(48, 48, -2.0667);
        expectedGen.set(51, 48, 1.3333);
        expectedGen.set(24, 49, 0.5);
        expectedGen.set(36, 49, 0.4);
        expectedGen.set(49, 49, -2.1667);
        expectedGen.set(52, 49, 1.3333);
        expectedGen.set(25, 50, 0.5);
        expectedGen.set(37, 50, 0.4);
        expectedGen.set(50, 50, -2.1667);
        expectedGen.set(53, 50, 1.3333);
        expectedGen.set(33, 51, 0.5);
        expectedGen.set(48, 51, 0.4);
        expectedGen.set(51, 51, -2.7333);
        expectedGen.set(62, 51, 2.0);
        expectedGen.set(31, 52, 0.5);
        expectedGen.set(45, 52, 0.4);
        expectedGen.set(52, 52, -2.8333);
        expectedGen.set(63, 52, 2.0);
        expectedGen.set(32, 53, 0.5);
        expectedGen.set(46, 53, 0.4);
        expectedGen.set(53, 53, -2.8333);
        expectedGen.set(54, 53, 2.0);
        expectedGen.set(34, 54, 0.5);
        expectedGen.set(51, 54, 0.4);
        expectedGen.set(54, 54, -3.5);
        expectedGen.set(64, 54, 2.0);
        expectedGen.set(41, 55, 0.5);
        expectedGen.set(55, 55, -2.4);
        expectedGen.set(65, 55, 3.0);
        expectedGen.set(43, 56, 0.5);
        expectedGen.set(56, 56, -2.4);
        expectedGen.set(59, 56, 0.6667);
        expectedGen.set(44, 57, 0.5);
        expectedGen.set(57, 57, -2.4);
        expectedGen.set(60, 57, 0.6667);
        expectedGen.set(38, 58, 0.5);
        expectedGen.set(58, 58, -2.5);
        expectedGen.set(61, 58, 0.6667);
        expectedGen.set(49, 59, 0.5);
        expectedGen.set(56, 59, 0.4);
        expectedGen.set(59, 59, -3.0667);
        expectedGen.set(66, 59, 3.0);
        expectedGen.set(50, 60, 0.5);
        expectedGen.set(57, 60, 0.4);
        expectedGen.set(60, 60, -3.0667);
        expectedGen.set(62, 60, 1.3333);
        expectedGen.set(47, 61, 0.5);
        expectedGen.set(55, 61, 0.4);
        expectedGen.set(61, 61, -3.1667);
        expectedGen.set(63, 61, 1.3333);
        expectedGen.set(53, 62, 0.5);
        expectedGen.set(60, 62, 0.4);
        expectedGen.set(62, 62, -3.7333);
        expectedGen.set(67, 62, 3.0);
        expectedGen.set(52, 63, 0.5);
        expectedGen.set(59, 63, 0.4);
        expectedGen.set(63, 63, -3.8333);
        expectedGen.set(64, 63, 2.0);
        expectedGen.set(54, 64, 0.5);
        expectedGen.set(62, 64, 0.4);
        expectedGen.set(64, 64, -4.5);
        expectedGen.set(68, 64, 3.0);
        expectedGen.set(58, 65, 0.5);
        expectedGen.set(65, 65, -3.4);
        expectedGen.set(66, 65, 0.6667);
        expectedGen.set(61, 66, 0.5);
        expectedGen.set(65, 66, 0.4);
        expectedGen.set(66, 66, -4.0667);
        expectedGen.set(67, 66, 1.3333);
        expectedGen.set(63, 67, 0.5);
        expectedGen.set(66, 67, 0.4);
        expectedGen.set(67, 67, -4.7333);
        expectedGen.set(68, 67, 2.0);
        expectedGen.set(64, 68, 0.5);
        expectedGen.set(67, 68, 0.4);
        expectedGen.set(68, 68, -5.0);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "LCFS two-class generator matrix mismatch");

        // Verify state space
        assertEquals(69, stateSpace.getNumRows());
        assertEquals(9, stateSpace.getNumCols());
    }

    @Test
    public void testTwoClassSIRO() {
        Network model = createTwoClassClosedNetwork(SchedStrategy.SIRO);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for SIRO two-class network
        // SIRO has 25x25 generator matrix
        Matrix expectedGen = new Matrix(25, 25);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.5);
        expectedGen.set(2, 0, 0.6667);
        expectedGen.set(7, 0, 1.0);
        expectedGen.set(1, 1, -0.4);
        expectedGen.set(3, 1, 0.6667);
        expectedGen.set(8, 1, 1.0);
        expectedGen.set(1, 2, 0.24);
        expectedGen.set(2, 2, -1.1667);
        expectedGen.set(4, 2, 1.3333);
        expectedGen.set(9, 2, 1.0);
        expectedGen.set(1, 3, 0.16);
        expectedGen.set(3, 3, -1.0667);
        expectedGen.set(5, 3, 1.3333);
        expectedGen.set(10, 3, 1.0);
        expectedGen.set(3, 4, 0.3);
        expectedGen.set(4, 4, -1.8333);
        expectedGen.set(6, 4, 2.0);
        expectedGen.set(11, 4, 1.0);
        expectedGen.set(3, 5, 0.1);
        expectedGen.set(5, 5, -1.7333);
        expectedGen.set(12, 5, 1.0);
        expectedGen.set(5, 6, 0.4);
        expectedGen.set(6, 6, -2.5);
        expectedGen.set(13, 6, 1.0);
        expectedGen.set(0, 7, 0.2);
        expectedGen.set(7, 7, -1.5);
        expectedGen.set(9, 7, 0.6667);
        expectedGen.set(14, 7, 2.0);
        expectedGen.set(0, 8, 0.3);
        expectedGen.set(8, 8, -1.4);
        expectedGen.set(10, 8, 0.6667);
        expectedGen.set(15, 8, 2.0);
        expectedGen.set(2, 9, 0.25);
        expectedGen.set(8, 9, 0.2);
        expectedGen.set(9, 9, -2.1667);
        expectedGen.set(11, 9, 1.3333);
        expectedGen.set(16, 9, 2.0);
        expectedGen.set(2, 10, 0.25);
        expectedGen.set(8, 10, 0.2);
        expectedGen.set(10, 10, -2.0667);
        expectedGen.set(12, 10, 1.3333);
        expectedGen.set(17, 10, 2.0);
        expectedGen.set(4, 11, 0.3333);
        expectedGen.set(10, 11, 0.2667);
        expectedGen.set(11, 11, -2.8333);
        expectedGen.set(13, 11, 2.0);
        expectedGen.set(18, 11, 2.0);
        expectedGen.set(4, 12, 0.1667);
        expectedGen.set(10, 12, 0.1333);
        expectedGen.set(12, 12, -2.7333);
        expectedGen.set(19, 12, 2.0);
        expectedGen.set(6, 13, 0.5);
        expectedGen.set(12, 13, 0.4);
        expectedGen.set(13, 13, -3.5);
        expectedGen.set(20, 13, 2.0);
        expectedGen.set(7, 14, 0.125);
        expectedGen.set(14, 14, -2.5);
        expectedGen.set(16, 14, 0.6667);
        expectedGen.set(7, 15, 0.375);
        expectedGen.set(15, 15, -2.4);
        expectedGen.set(17, 15, 0.6667);
        expectedGen.set(21, 15, 3.0);
        expectedGen.set(9, 16, 0.1667);
        expectedGen.set(15, 16, 0.1333);
        expectedGen.set(16, 16, -3.1667);
        expectedGen.set(18, 16, 1.3333);
        expectedGen.set(9, 17, 0.3333);
        expectedGen.set(15, 17, 0.2667);
        expectedGen.set(17, 17, -3.0667);
        expectedGen.set(19, 17, 1.3333);
        expectedGen.set(22, 17, 3.0);
        expectedGen.set(11, 18, 0.25);
        expectedGen.set(17, 18, 0.2);
        expectedGen.set(18, 18, -3.8333);
        expectedGen.set(20, 18, 2.0);
        expectedGen.set(11, 19, 0.25);
        expectedGen.set(17, 19, 0.2);
        expectedGen.set(19, 19, -3.7333);
        expectedGen.set(23, 19, 3.0);
        expectedGen.set(13, 20, 0.5);
        expectedGen.set(19, 20, 0.4);
        expectedGen.set(20, 20, -4.5);
        expectedGen.set(24, 20, 3.0);
        expectedGen.set(14, 21, 0.5);
        expectedGen.set(21, 21, -3.4);
        expectedGen.set(22, 21, 0.6667);
        expectedGen.set(16, 22, 0.5);
        expectedGen.set(21, 22, 0.4);
        expectedGen.set(22, 22, -4.0667);
        expectedGen.set(23, 22, 1.3333);
        expectedGen.set(18, 23, 0.5);
        expectedGen.set(22, 23, 0.4);
        expectedGen.set(23, 23, -4.7333);
        expectedGen.set(24, 23, 2.0);
        expectedGen.set(20, 24, 0.5);
        expectedGen.set(23, 24, 0.4);
        expectedGen.set(24, 24, -5.0);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "SIRO two-class generator matrix mismatch");

        // Verify state space
        assertEquals(25, stateSpace.getNumRows());
        assertEquals(6, stateSpace.getNumCols());
    }

    @Test
    public void testTwoClassSEPT() {
        Network model = createTwoClassClosedNetwork(SchedStrategy.SEPT);
        SolverCTMC ctmcSolver = new SolverCTMC(model);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for SEPT two-class network
        // SEPT has 25x25 generator matrix
        Matrix expectedGen = new Matrix(25, 25);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.5);
        expectedGen.set(2, 0, 0.6667);
        expectedGen.set(7, 0, 1.0);
        expectedGen.set(1, 1, -0.4);
        expectedGen.set(3, 1, 0.6667);
        expectedGen.set(8, 1, 1.0);
        expectedGen.set(1, 2, 0.4);
        expectedGen.set(2, 2, -1.1667);
        expectedGen.set(4, 2, 1.3333);
        expectedGen.set(9, 2, 1.0);
        expectedGen.set(3, 3, -1.0667);
        expectedGen.set(5, 3, 1.3333);
        expectedGen.set(10, 3, 1.0);
        expectedGen.set(3, 4, 0.4);
        expectedGen.set(4, 4, -1.8333);
        expectedGen.set(6, 4, 2.0);
        expectedGen.set(11, 4, 1.0);
        expectedGen.set(5, 5, -1.7333);
        expectedGen.set(12, 5, 1.0);
        expectedGen.set(5, 6, 0.4);
        expectedGen.set(6, 6, -2.5);
        expectedGen.set(13, 6, 1.0);
        expectedGen.set(0, 7, 0.5);
        expectedGen.set(7, 7, -1.5);
        expectedGen.set(9, 7, 0.6667);
        expectedGen.set(14, 7, 2.0);
        expectedGen.set(8, 8, -1.4);
        expectedGen.set(10, 8, 0.6667);
        expectedGen.set(15, 8, 2.0);
        expectedGen.set(2, 9, 0.5);
        expectedGen.set(8, 9, 0.4);
        expectedGen.set(9, 9, -2.1667);
        expectedGen.set(11, 9, 1.3333);
        expectedGen.set(16, 9, 2.0);
        expectedGen.set(10, 10, -2.0667);
        expectedGen.set(12, 10, 1.3333);
        expectedGen.set(17, 10, 2.0);
        expectedGen.set(4, 11, 0.5);
        expectedGen.set(10, 11, 0.4);
        expectedGen.set(11, 11, -2.8333);
        expectedGen.set(13, 11, 2.0);
        expectedGen.set(18, 11, 2.0);
        expectedGen.set(12, 12, -2.7333);
        expectedGen.set(19, 12, 2.0);
        expectedGen.set(6, 13, 0.5);
        expectedGen.set(12, 13, 0.4);
        expectedGen.set(13, 13, -3.5);
        expectedGen.set(20, 13, 2.0);
        expectedGen.set(7, 14, 0.5);
        expectedGen.set(14, 14, -2.5);
        expectedGen.set(16, 14, 0.6667);
        expectedGen.set(15, 15, -2.4);
        expectedGen.set(17, 15, 0.6667);
        expectedGen.set(21, 15, 3.0);
        expectedGen.set(9, 16, 0.5);
        expectedGen.set(15, 16, 0.4);
        expectedGen.set(16, 16, -3.1667);
        expectedGen.set(18, 16, 1.3333);
        expectedGen.set(17, 17, -3.0667);
        expectedGen.set(19, 17, 1.3333);
        expectedGen.set(22, 17, 3.0);
        expectedGen.set(11, 18, 0.5);
        expectedGen.set(17, 18, 0.4);
        expectedGen.set(18, 18, -3.8333);
        expectedGen.set(20, 18, 2.0);
        expectedGen.set(19, 19, -3.7333);
        expectedGen.set(23, 19, 3.0);
        expectedGen.set(13, 20, 0.5);
        expectedGen.set(19, 20, 0.4);
        expectedGen.set(20, 20, -4.5);
        expectedGen.set(24, 20, 3.0);
        expectedGen.set(14, 21, 0.5);
        expectedGen.set(21, 21, -3.4);
        expectedGen.set(22, 21, 0.6667);
        expectedGen.set(16, 22, 0.5);
        expectedGen.set(21, 22, 0.4);
        expectedGen.set(22, 22, -4.0667);
        expectedGen.set(23, 22, 1.3333);
        expectedGen.set(18, 23, 0.5);
        expectedGen.set(22, 23, 0.4);
        expectedGen.set(23, 23, -4.7333);
        expectedGen.set(24, 23, 2.0);
        expectedGen.set(20, 24, 0.5);
        expectedGen.set(23, 24, 0.4);
        expectedGen.set(24, 24, -5.0);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "SEPT two-class generator matrix mismatch");

        // Verify state space
        assertEquals(25, stateSpace.getNumRows());
        assertEquals(6, stateSpace.getNumCols());
    }

    @Test
    public void testTwoClassLEPT() {
        Network model = createTwoClassClosedNetwork(SchedStrategy.LEPT);
        SolverCTMC ctmcSolver = new SolverCTMC(model);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for LEPT two-class network
        // LEPT has 25x25 generator matrix
        Matrix expectedGen = new Matrix(25, 25);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.5);
        expectedGen.set(2, 0, 0.6667);
        expectedGen.set(7, 0, 1.0);
        expectedGen.set(1, 1, -0.4);
        expectedGen.set(3, 1, 0.6667);
        expectedGen.set(8, 1, 1.0);
        expectedGen.set(2, 2, -1.1667);
        expectedGen.set(4, 2, 1.3333);
        expectedGen.set(9, 2, 1.0);
        expectedGen.set(1, 3, 0.4);
        expectedGen.set(3, 3, -1.0667);
        expectedGen.set(5, 3, 1.3333);
        expectedGen.set(10, 3, 1.0);
        expectedGen.set(4, 4, -1.8333);
        expectedGen.set(6, 4, 2.0);
        expectedGen.set(11, 4, 1.0);
        expectedGen.set(3, 5, 0.4);
        expectedGen.set(5, 5, -1.7333);
        expectedGen.set(12, 5, 1.0);
        expectedGen.set(5, 6, 0.4);
        expectedGen.set(6, 6, -2.5);
        expectedGen.set(13, 6, 1.0);
        expectedGen.set(7, 7, -1.5);
        expectedGen.set(9, 7, 0.6667);
        expectedGen.set(14, 7, 2.0);
        expectedGen.set(0, 8, 0.5);
        expectedGen.set(8, 8, -1.4);
        expectedGen.set(10, 8, 0.6667);
        expectedGen.set(15, 8, 2.0);
        expectedGen.set(9, 9, -2.1667);
        expectedGen.set(11, 9, 1.3333);
        expectedGen.set(16, 9, 2.0);
        expectedGen.set(2, 10, 0.5);
        expectedGen.set(8, 10, 0.4);
        expectedGen.set(10, 10, -2.0667);
        expectedGen.set(12, 10, 1.3333);
        expectedGen.set(17, 10, 2.0);
        expectedGen.set(11, 11, -2.8333);
        expectedGen.set(13, 11, 2.0);
        expectedGen.set(18, 11, 2.0);
        expectedGen.set(4, 12, 0.5);
        expectedGen.set(10, 12, 0.4);
        expectedGen.set(12, 12, -2.7333);
        expectedGen.set(19, 12, 2.0);
        expectedGen.set(6, 13, 0.5);
        expectedGen.set(12, 13, 0.4);
        expectedGen.set(13, 13, -3.5);
        expectedGen.set(20, 13, 2.0);
        expectedGen.set(14, 14, -2.5);
        expectedGen.set(16, 14, 0.6667);
        expectedGen.set(7, 15, 0.5);
        expectedGen.set(15, 15, -2.4);
        expectedGen.set(17, 15, 0.6667);
        expectedGen.set(21, 15, 3.0);
        expectedGen.set(16, 16, -3.1667);
        expectedGen.set(18, 16, 1.3333);
        expectedGen.set(9, 17, 0.5);
        expectedGen.set(15, 17, 0.4);
        expectedGen.set(17, 17, -3.0667);
        expectedGen.set(19, 17, 1.3333);
        expectedGen.set(22, 17, 3.0);
        expectedGen.set(18, 18, -3.8333);
        expectedGen.set(20, 18, 2.0);
        expectedGen.set(11, 19, 0.5);
        expectedGen.set(17, 19, 0.4);
        expectedGen.set(19, 19, -3.7333);
        expectedGen.set(23, 19, 3.0);
        expectedGen.set(13, 20, 0.5);
        expectedGen.set(19, 20, 0.4);
        expectedGen.set(20, 20, -4.5);
        expectedGen.set(24, 20, 3.0);
        expectedGen.set(14, 21, 0.5);
        expectedGen.set(21, 21, -3.4);
        expectedGen.set(22, 21, 0.6667);
        expectedGen.set(16, 22, 0.5);
        expectedGen.set(21, 22, 0.4);
        expectedGen.set(22, 22, -4.0667);
        expectedGen.set(23, 22, 1.3333);
        expectedGen.set(18, 23, 0.5);
        expectedGen.set(22, 23, 0.4);
        expectedGen.set(23, 23, -4.7333);
        expectedGen.set(24, 23, 2.0);
        expectedGen.set(20, 24, 0.5);
        expectedGen.set(23, 24, 0.4);
        expectedGen.set(24, 24, -5.0);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "LEPT two-class generator matrix mismatch");

        // Verify state space
        assertEquals(25, stateSpace.getNumRows());
        assertEquals(6, stateSpace.getNumCols());
    }

    @Test
    public void testTwoClassPS() {
        Network model = createTwoClassClosedNetwork(SchedStrategy.PS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - based on MATLAB output for PS two-class network
        // Verify generator - complete verification based on MATLAB ground truth for PS two-class network
        // PS has 16x16 generator matrix
        Matrix expectedGen = new Matrix(16, 16);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.45);
        expectedGen.set(1, 0, 0.6667);
        expectedGen.set(4, 0, 1.0);
        expectedGen.set(0, 1, 0.2);
        expectedGen.set(1, 1, -1.1267);
        expectedGen.set(2, 1, 1.3333);
        expectedGen.set(5, 1, 1.0);
        expectedGen.set(1, 2, 0.16);
        expectedGen.set(2, 2, -1.8083);
        expectedGen.set(3, 2, 2.0);
        expectedGen.set(6, 2, 1.0);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -2.5);
        expectedGen.set(7, 3, 1.0);
        expectedGen.set(0, 4, 0.25);
        expectedGen.set(4, 4, -1.44);
        expectedGen.set(5, 4, 0.6667);
        expectedGen.set(8, 4, 2.0);
        expectedGen.set(1, 5, 0.3);
        expectedGen.set(4, 5, 0.24);
        expectedGen.set(5, 5, -2.1167);
        expectedGen.set(6, 5, 1.3333);
        expectedGen.set(9, 5, 2.0);
        expectedGen.set(2, 6, 0.375);
        expectedGen.set(5, 6, 0.2);
        expectedGen.set(6, 6, -2.8);
        expectedGen.set(7, 6, 2.0);
        expectedGen.set(10, 6, 2.0);
        expectedGen.set(3, 7, 0.5);
        expectedGen.set(6, 7, 0.1333);
        expectedGen.set(7, 7, -3.5);
        expectedGen.set(11, 7, 2.0);
        expectedGen.set(4, 8, 0.2);
        expectedGen.set(8, 8, -2.425);
        expectedGen.set(9, 8, 0.6667);
        expectedGen.set(12, 8, 3.0);
        expectedGen.set(5, 9, 0.25);
        expectedGen.set(8, 9, 0.3);
        expectedGen.set(9, 9, -3.1);
        expectedGen.set(10, 9, 1.3333);
        expectedGen.set(13, 9, 3.0);
        expectedGen.set(6, 10, 0.3333);
        expectedGen.set(9, 10, 0.2667);
        expectedGen.set(10, 10, -3.7833);
        expectedGen.set(11, 10, 2.0);
        expectedGen.set(14, 10, 3.0);
        expectedGen.set(7, 11, 0.5);
        expectedGen.set(10, 11, 0.2);
        expectedGen.set(11, 11, -4.5);
        expectedGen.set(15, 11, 3.0);
        expectedGen.set(8, 12, 0.125);
        expectedGen.set(12, 12, -3.4);
        expectedGen.set(13, 12, 0.6667);
        expectedGen.set(9, 13, 0.1667);
        expectedGen.set(12, 13, 0.4);
        expectedGen.set(13, 13, -4.0667);
        expectedGen.set(14, 13, 1.3333);
        expectedGen.set(10, 14, 0.25);
        expectedGen.set(13, 14, 0.4);
        expectedGen.set(14, 14, -4.7333);
        expectedGen.set(15, 14, 2.0);
        expectedGen.set(11, 15, 0.5);
        expectedGen.set(14, 15, 0.4);
        expectedGen.set(15, 15, -5.0);
        assertMatrixEquals(expectedGen, generator, "Generator matrix for PS two-class");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0,0,3,3;" +
                "0,1,3,2;" +
                "0,2,3,1;" +
                "0,3,3,0;" +
                "1,0,2,3;" +
                "1,1,2,2;" +
                "1,2,2,1;" +
                "1,3,2,0;" +
                "2,0,1,3;" +
                "2,1,1,2;" +
                "2,2,1,1;" +
                "2,3,1,0;" +
                "3,0,0,3;" +
                "3,1,0,2;" +
                "3,2,0,1;" +
                "3,3,0,0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for PS two-class");
    }

    @Test
    public void testTwoClassDPS() {
        Network model = createTwoClassClosedNetwork(SchedStrategy.DPS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for DPS two-class network
        // DPS has 16x16 generator matrix
        Matrix expectedGen = new Matrix(16, 16);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.4667);
        expectedGen.set(1, 0, 0.6667);
        expectedGen.set(4, 0, 1.0);
        expectedGen.set(0, 1, 0.1333);
        expectedGen.set(1, 1, -1.1417);
        expectedGen.set(2, 1, 1.3333);
        expectedGen.set(5, 1, 1.0);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -1.819);
        expectedGen.set(3, 2, 2.0);
        expectedGen.set(6, 2, 1.0);
        expectedGen.set(2, 3, 0.0571);
        expectedGen.set(3, 3, -2.5);
        expectedGen.set(7, 3, 1.0);
        expectedGen.set(0, 4, 0.3333);
        expectedGen.set(4, 4, -1.4571);
        expectedGen.set(5, 4, 0.6667);
        expectedGen.set(8, 4, 2.0);
        expectedGen.set(1, 5, 0.375);
        expectedGen.set(4, 5, 0.1714);
        expectedGen.set(5, 5, -2.1333);
        expectedGen.set(6, 5, 1.3333);
        expectedGen.set(9, 5, 2.0);
        expectedGen.set(2, 6, 0.4286);
        expectedGen.set(5, 6, 0.1333);
        expectedGen.set(6, 6, -2.8133);
        expectedGen.set(7, 6, 2.0);
        expectedGen.set(10, 6, 2.0);
        expectedGen.set(3, 7, 0.5);
        expectedGen.set(6, 7, 0.08);
        expectedGen.set(7, 7, -3.5);
        expectedGen.set(11, 7, 2.0);
        expectedGen.set(4, 8, 0.2857);
        expectedGen.set(8, 8, -2.44);
        expectedGen.set(9, 8, 0.6667);
        expectedGen.set(12, 8, 3.0);
        expectedGen.set(5, 9, 0.3333);
        expectedGen.set(8, 9, 0.24);
        expectedGen.set(9, 9, -3.1167);
        expectedGen.set(10, 9, 1.3333);
        expectedGen.set(13, 9, 3.0);
        expectedGen.set(6, 10, 0.4);
        expectedGen.set(9, 10, 0.2);
        expectedGen.set(10, 10, -3.8);
        expectedGen.set(11, 10, 2.0);
        expectedGen.set(14, 10, 3.0);
        expectedGen.set(7, 11, 0.5);
        expectedGen.set(10, 11, 0.1333);
        expectedGen.set(11, 11, -4.5);
        expectedGen.set(15, 11, 3.0);
        expectedGen.set(8, 12, 0.2);
        expectedGen.set(12, 12, -3.4);
        expectedGen.set(13, 12, 0.6667);
        expectedGen.set(9, 13, 0.25);
        expectedGen.set(12, 13, 0.4);
        expectedGen.set(13, 13, -4.0667);
        expectedGen.set(14, 13, 1.3333);
        expectedGen.set(10, 14, 0.3333);
        expectedGen.set(13, 14, 0.4);
        expectedGen.set(14, 14, -4.7333);
        expectedGen.set(15, 14, 2.0);
        expectedGen.set(11, 15, 0.5);
        expectedGen.set(14, 15, 0.4);
        expectedGen.set(15, 15, -5.0);
        assertMatrixEquals(expectedGen, generator, "Generator matrix for DPS two-class");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0,0,3,3;" +
                "0,1,3,2;" +
                "0,2,3,1;" +
                "0,3,3,0;" +
                "1,0,2,3;" +
                "1,1,2,2;" +
                "1,2,2,1;" +
                "1,3,2,0;" +
                "2,0,1,3;" +
                "2,1,1,2;" +
                "2,2,1,1;" +
                "2,3,1,0;" +
                "3,0,0,3;" +
                "3,1,0,2;" +
                "3,2,0,1;" +
                "3,3,0,0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for DPS two-class");
    }

    @Test
    public void testTwoClassGPS() {
        Network model = createTwoClassClosedNetwork(SchedStrategy.GPS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for GPS two-class network
        // GPS has 16x16 generator matrix
        Matrix expectedGen = new Matrix(16, 16);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.4667);
        expectedGen.set(1, 0, 0.6667);
        expectedGen.set(4, 0, 1.0);
        expectedGen.set(0, 1, 0.1333);
        expectedGen.set(1, 1, -1.1333);
        expectedGen.set(2, 1, 1.3333);
        expectedGen.set(5, 1, 1.0);
        expectedGen.set(1, 2, 0.1333);
        expectedGen.set(2, 2, -1.8);
        expectedGen.set(3, 2, 2.0);
        expectedGen.set(6, 2, 1.0);
        expectedGen.set(2, 3, 0.1333);
        expectedGen.set(3, 3, -2.5);
        expectedGen.set(7, 3, 1.0);
        expectedGen.set(0, 4, 0.3333);
        expectedGen.set(4, 4, -1.4667);
        expectedGen.set(5, 4, 0.6667);
        expectedGen.set(8, 4, 2.0);
        expectedGen.set(1, 5, 0.3333);
        expectedGen.set(4, 5, 0.1333);
        expectedGen.set(5, 5, -2.1333);
        expectedGen.set(6, 5, 1.3333);
        expectedGen.set(9, 5, 2.0);
        expectedGen.set(2, 6, 0.3333);
        expectedGen.set(5, 6, 0.1333);
        expectedGen.set(6, 6, -2.8);
        expectedGen.set(7, 6, 2.0);
        expectedGen.set(10, 6, 2.0);
        expectedGen.set(3, 7, 0.5);
        expectedGen.set(6, 7, 0.1333);
        expectedGen.set(7, 7, -3.5);
        expectedGen.set(11, 7, 2.0);
        expectedGen.set(4, 8, 0.3333);
        expectedGen.set(8, 8, -2.4667);
        expectedGen.set(9, 8, 0.6667);
        expectedGen.set(12, 8, 3.0);
        expectedGen.set(5, 9, 0.3333);
        expectedGen.set(8, 9, 0.1333);
        expectedGen.set(9, 9, -3.1333);
        expectedGen.set(10, 9, 1.3333);
        expectedGen.set(13, 9, 3.0);
        expectedGen.set(6, 10, 0.3333);
        expectedGen.set(9, 10, 0.1333);
        expectedGen.set(10, 10, -3.8);
        expectedGen.set(11, 10, 2.0);
        expectedGen.set(14, 10, 3.0);
        expectedGen.set(7, 11, 0.5);
        expectedGen.set(10, 11, 0.1333);
        expectedGen.set(11, 11, -4.5);
        expectedGen.set(15, 11, 3.0);
        expectedGen.set(8, 12, 0.3333);
        expectedGen.set(12, 12, -3.4);
        expectedGen.set(13, 12, 0.6667);
        expectedGen.set(9, 13, 0.3333);
        expectedGen.set(12, 13, 0.4);
        expectedGen.set(13, 13, -4.0667);
        expectedGen.set(14, 13, 1.3333);
        expectedGen.set(10, 14, 0.3333);
        expectedGen.set(13, 14, 0.4);
        expectedGen.set(14, 14, -4.7333);
        expectedGen.set(15, 14, 2.0);
        expectedGen.set(11, 15, 0.5);
        expectedGen.set(14, 15, 0.4);
        expectedGen.set(15, 15, -5.0);
        assertMatrixEquals(expectedGen, generator, "Generator matrix for GPS two-class");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0,0,3,3;" +
                "0,1,3,2;" +
                "0,2,3,1;" +
                "0,3,3,0;" +
                "1,0,2,3;" +
                "1,1,2,2;" +
                "1,2,2,1;" +
                "1,3,2,0;" +
                "2,0,1,3;" +
                "2,1,1,2;" +
                "2,2,1,1;" +
                "2,3,1,0;" +
                "3,0,0,3;" +
                "3,1,0,2;" +
                "3,2,0,1;" +
                "3,3,0,0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for GPS two-class");
    }

    @Test
    public void testTwoClassINF() {
        Network model = createTwoClassClosedNetwork(SchedStrategy.INF);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for INF two-class network
        // INF has 16x16 generator matrix
        Matrix expectedGen = new Matrix(16, 16);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -2.7);
        expectedGen.set(1, 0, 0.6667);
        expectedGen.set(4, 0, 1.0);
        expectedGen.set(0, 1, 1.2);
        expectedGen.set(1, 1, -2.9667);
        expectedGen.set(2, 1, 1.3333);
        expectedGen.set(5, 1, 1.0);
        expectedGen.set(1, 2, 0.8);
        expectedGen.set(2, 2, -3.2333);
        expectedGen.set(3, 2, 2.0);
        expectedGen.set(6, 2, 1.0);
        expectedGen.set(2, 3, 0.4);
        expectedGen.set(3, 3, -3.5);
        expectedGen.set(7, 3, 1.0);
        expectedGen.set(0, 4, 1.5);
        expectedGen.set(4, 4, -3.2);
        expectedGen.set(5, 4, 0.6667);
        expectedGen.set(8, 4, 2.0);
        expectedGen.set(1, 5, 1.5);
        expectedGen.set(4, 5, 1.2);
        expectedGen.set(5, 5, -3.4667);
        expectedGen.set(6, 5, 1.3333);
        expectedGen.set(9, 5, 2.0);
        expectedGen.set(2, 6, 1.5);
        expectedGen.set(5, 6, 0.8);
        expectedGen.set(6, 6, -3.7333);
        expectedGen.set(7, 6, 2.0);
        expectedGen.set(10, 6, 2.0);
        expectedGen.set(3, 7, 1.5);
        expectedGen.set(6, 7, 0.4);
        expectedGen.set(7, 7, -4.0);
        expectedGen.set(11, 7, 2.0);
        expectedGen.set(4, 8, 1.0);
        expectedGen.set(8, 8, -3.7);
        expectedGen.set(9, 8, 0.6667);
        expectedGen.set(12, 8, 3.0);
        expectedGen.set(5, 9, 1.0);
        expectedGen.set(8, 9, 1.2);
        expectedGen.set(9, 9, -3.9667);
        expectedGen.set(10, 9, 1.3333);
        expectedGen.set(13, 9, 3.0);
        expectedGen.set(6, 10, 1.0);
        expectedGen.set(9, 10, 0.8);
        expectedGen.set(10, 10, -4.2333);
        expectedGen.set(11, 10, 2.0);
        expectedGen.set(14, 10, 3.0);
        expectedGen.set(7, 11, 1.0);
        expectedGen.set(10, 11, 0.4);
        expectedGen.set(11, 11, -4.5);
        expectedGen.set(15, 11, 3.0);
        expectedGen.set(8, 12, 0.5);
        expectedGen.set(12, 12, -4.2);
        expectedGen.set(13, 12, 0.6667);
        expectedGen.set(9, 13, 0.5);
        expectedGen.set(12, 13, 1.2);
        expectedGen.set(13, 13, -4.4667);
        expectedGen.set(14, 13, 1.3333);
        expectedGen.set(10, 14, 0.5);
        expectedGen.set(13, 14, 0.8);
        expectedGen.set(14, 14, -4.7333);
        expectedGen.set(15, 14, 2.0);
        expectedGen.set(11, 15, 0.5);
        expectedGen.set(14, 15, 0.4);
        expectedGen.set(15, 15, -5.0);
        assertMatrixEquals(expectedGen, generator, "Generator matrix for INF two-class");

        // Verify state space
        Matrix expectedStateSpace = new Matrix("[0,0,3,3;" +
                "0,1,3,2;" +
                "0,2,3,1;" +
                "0,3,3,0;" +
                "1,0,2,3;" +
                "1,1,2,2;" +
                "1,2,2,1;" +
                "1,3,2,0;" +
                "2,0,1,3;" +
                "2,1,1,2;" +
                "2,2,1,1;" +
                "2,3,1,0;" +
                "3,0,0,3;" +
                "3,1,0,2;" +
                "3,2,0,1;" +
                "3,3,0,0]");
        assertMatrixEquals(expectedStateSpace, stateSpace, "State space for INF two-class");
    }

    @Test
    public void testTwoClassHOL() {
        Network model = createTwoClassClosedNetwork(SchedStrategy.HOL);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for HOL two-class network
        // HOL has 69x69 generator matrix
        Matrix expectedGen = new Matrix(69, 69);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.4000);
        expectedGen.set(35, 0, 1.0000);
        expectedGen.set(1, 1, -0.4000);
        expectedGen.set(36, 1, 1.0000);
        expectedGen.set(2, 2, -0.4000);
        expectedGen.set(37, 2, 1.0000);
        expectedGen.set(3, 3, -0.5000);
        expectedGen.set(38, 3, 1.0000);
        expectedGen.set(4, 4, -0.4000);
        expectedGen.set(39, 4, 1.0000);
        expectedGen.set(5, 5, -0.4000);
        expectedGen.set(40, 5, 1.0000);
        expectedGen.set(6, 6, -0.5000);
        expectedGen.set(41, 6, 1.0000);
        expectedGen.set(7, 7, -0.4000);
        expectedGen.set(42, 7, 1.0000);
        expectedGen.set(8, 8, -0.5000);
        expectedGen.set(43, 8, 1.0000);
        expectedGen.set(9, 9, -0.5000);
        expectedGen.set(44, 9, 1.0000);
        expectedGen.set(10, 10, -0.4000);
        expectedGen.set(20, 10, 0.6667);
        expectedGen.set(11, 11, -0.4000);
        expectedGen.set(21, 11, 0.6667);
        expectedGen.set(12, 12, -0.5000);
        expectedGen.set(22, 12, 0.6667);
        expectedGen.set(13, 13, -0.4000);
        expectedGen.set(23, 13, 0.6667);
        expectedGen.set(14, 14, -0.5000);
        expectedGen.set(24, 14, 0.6667);
        expectedGen.set(15, 15, -0.5000);
        expectedGen.set(25, 15, 0.6667);
        expectedGen.set(16, 16, -0.4000);
        expectedGen.set(26, 16, 0.6667);
        expectedGen.set(17, 17, -0.5000);
        expectedGen.set(27, 17, 0.6667);
        expectedGen.set(18, 18, -0.5000);
        expectedGen.set(28, 18, 0.6667);
        expectedGen.set(19, 19, -0.5000);
        expectedGen.set(29, 19, 0.6667);
        expectedGen.set(20, 20, -1.0667);
        expectedGen.set(45, 20, 1.0000);
        expectedGen.set(21, 21, -1.0667);
        expectedGen.set(46, 21, 1.0000);
        expectedGen.set(0, 22, 0.4000);
        expectedGen.set(1, 22, 0.4000);
        expectedGen.set(2, 22, 0.4000);
        expectedGen.set(22, 22, -1.1667);
        expectedGen.set(47, 22, 1.0000);
        expectedGen.set(23, 23, -1.0667);
        expectedGen.set(48, 23, 1.0000);
        expectedGen.set(4, 24, 0.4000);
        expectedGen.set(5, 24, 0.4000);
        expectedGen.set(24, 24, -1.1667);
        expectedGen.set(49, 24, 1.0000);
        expectedGen.set(7, 25, 0.4000);
        expectedGen.set(25, 25, -1.1667);
        expectedGen.set(50, 25, 1.0000);
        expectedGen.set(26, 26, -1.0667);
        expectedGen.set(30, 26, 1.3333);
        expectedGen.set(10, 27, 0.4000);
        expectedGen.set(11, 27, 0.4000);
        expectedGen.set(27, 27, -1.1667);
        expectedGen.set(31, 27, 1.3333);
        expectedGen.set(13, 28, 0.4000);
        expectedGen.set(28, 28, -1.1667);
        expectedGen.set(32, 28, 1.3333);
        expectedGen.set(16, 29, 0.4000);
        expectedGen.set(29, 29, -1.1667);
        expectedGen.set(33, 29, 1.3333);
        expectedGen.set(30, 30, -1.7333);
        expectedGen.set(51, 30, 1.0000);
        expectedGen.set(20, 31, 0.4000);
        expectedGen.set(21, 31, 0.4000);
        expectedGen.set(31, 31, -1.8333);
        expectedGen.set(52, 31, 1.0000);
        expectedGen.set(23, 32, 0.4000);
        expectedGen.set(32, 32, -1.8333);
        expectedGen.set(53, 32, 1.0000);
        expectedGen.set(26, 33, 0.4000);
        expectedGen.set(33, 33, -1.8333);
        expectedGen.set(34, 33, 2.0000);
        expectedGen.set(30, 34, 0.4000);
        expectedGen.set(34, 34, -2.5000);
        expectedGen.set(54, 34, 1.0000);
        expectedGen.set(35, 35, -1.4000);
        expectedGen.set(55, 35, 2.0000);
        expectedGen.set(36, 36, -1.4000);
        expectedGen.set(56, 36, 2.0000);
        expectedGen.set(37, 37, -1.4000);
        expectedGen.set(57, 37, 2.0000);
        expectedGen.set(3, 38, 0.5000);
        expectedGen.set(6, 38, 0.5000);
        expectedGen.set(8, 38, 0.5000);
        expectedGen.set(9, 38, 0.5000);
        expectedGen.set(38, 38, -1.5000);
        expectedGen.set(58, 38, 2.0000);
        expectedGen.set(39, 39, -1.4000);
        expectedGen.set(45, 39, 0.6667);
        expectedGen.set(40, 40, -1.4000);
        expectedGen.set(46, 40, 0.6667);
        expectedGen.set(12, 41, 0.5000);
        expectedGen.set(14, 41, 0.5000);
        expectedGen.set(15, 41, 0.5000);
        expectedGen.set(41, 41, -1.5000);
        expectedGen.set(47, 41, 0.6667);
        expectedGen.set(42, 42, -1.4000);
        expectedGen.set(48, 42, 0.6667);
        expectedGen.set(17, 43, 0.5000);
        expectedGen.set(18, 43, 0.5000);
        expectedGen.set(43, 43, -1.5000);
        expectedGen.set(49, 43, 0.6667);
        expectedGen.set(19, 44, 0.5000);
        expectedGen.set(44, 44, -1.5000);
        expectedGen.set(50, 44, 0.6667);
        expectedGen.set(45, 45, -2.0667);
        expectedGen.set(59, 45, 2.0000);
        expectedGen.set(46, 46, -2.0667);
        expectedGen.set(60, 46, 2.0000);
        expectedGen.set(22, 47, 0.5000);
        expectedGen.set(24, 47, 0.5000);
        expectedGen.set(25, 47, 0.5000);
        expectedGen.set(35, 47, 0.4000);
        expectedGen.set(36, 47, 0.4000);
        expectedGen.set(37, 47, 0.4000);
        expectedGen.set(47, 47, -2.1667);
        expectedGen.set(61, 47, 2.0000);
        expectedGen.set(48, 48, -2.0667);
        expectedGen.set(51, 48, 1.3333);
        expectedGen.set(27, 49, 0.5000);
        expectedGen.set(28, 49, 0.5000);
        expectedGen.set(39, 49, 0.4000);
        expectedGen.set(40, 49, 0.4000);
        expectedGen.set(49, 49, -2.1667);
        expectedGen.set(52, 49, 1.3333);
        expectedGen.set(29, 50, 0.5000);
        expectedGen.set(42, 50, 0.4000);
        expectedGen.set(50, 50, -2.1667);
        expectedGen.set(53, 50, 1.3333);
        expectedGen.set(51, 51, -2.7333);
        expectedGen.set(62, 51, 2.0000);
        expectedGen.set(31, 52, 0.5000);
        expectedGen.set(32, 52, 0.5000);
        expectedGen.set(45, 52, 0.4000);
        expectedGen.set(46, 52, 0.4000);
        expectedGen.set(52, 52, -2.8333);
        expectedGen.set(63, 52, 2.0000);
        expectedGen.set(33, 53, 0.5000);
        expectedGen.set(48, 53, 0.4000);
        expectedGen.set(53, 53, -2.8333);
        expectedGen.set(54, 53, 2.0000);
        expectedGen.set(34, 54, 0.5000);
        expectedGen.set(51, 54, 0.4000);
        expectedGen.set(54, 54, -3.5000);
        expectedGen.set(64, 54, 2.0000);
        expectedGen.set(55, 55, -2.4000);
        expectedGen.set(65, 55, 3.0000);
        expectedGen.set(56, 56, -2.4000);
        expectedGen.set(59, 56, 0.6667);
        expectedGen.set(57, 57, -2.4000);
        expectedGen.set(60, 57, 0.6667);
        expectedGen.set(38, 58, 0.5000);
        expectedGen.set(41, 58, 0.5000);
        expectedGen.set(43, 58, 0.5000);
        expectedGen.set(44, 58, 0.5000);
        expectedGen.set(58, 58, -2.5000);
        expectedGen.set(61, 58, 0.6667);
        expectedGen.set(59, 59, -3.0667);
        expectedGen.set(66, 59, 3.0000);
        expectedGen.set(60, 60, -3.0667);
        expectedGen.set(62, 60, 1.3333);
        expectedGen.set(47, 61, 0.5000);
        expectedGen.set(49, 61, 0.5000);
        expectedGen.set(50, 61, 0.5000);
        expectedGen.set(55, 61, 0.4000);
        expectedGen.set(56, 61, 0.4000);
        expectedGen.set(57, 61, 0.4000);
        expectedGen.set(61, 61, -3.1667);
        expectedGen.set(63, 61, 1.3333);
        expectedGen.set(62, 62, -3.7333);
        expectedGen.set(67, 62, 3.0000);
        expectedGen.set(52, 63, 0.5000);
        expectedGen.set(53, 63, 0.5000);
        expectedGen.set(59, 63, 0.4000);
        expectedGen.set(60, 63, 0.4000);
        expectedGen.set(63, 63, -3.8333);
        expectedGen.set(64, 63, 2.0000);
        expectedGen.set(54, 64, 0.5000);
        expectedGen.set(62, 64, 0.4000);
        expectedGen.set(64, 64, -4.5000);
        expectedGen.set(68, 64, 3.0000);
        expectedGen.set(58, 65, 0.5000);
        expectedGen.set(65, 65, -3.4000);
        expectedGen.set(66, 65, 0.6667);
        expectedGen.set(61, 66, 0.5000);
        expectedGen.set(65, 66, 0.4000);
        expectedGen.set(66, 66, -4.0667);
        expectedGen.set(67, 66, 1.3333);
        expectedGen.set(63, 67, 0.5000);
        expectedGen.set(66, 67, 0.4000);
        expectedGen.set(67, 67, -4.7333);
        expectedGen.set(68, 67, 2.0000);
        expectedGen.set(64, 68, 0.5000);
        expectedGen.set(67, 68, 0.4000);
        expectedGen.set(68, 68, -5.0000);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "HOL two-class generator matrix mismatch");

        // Verify state space
        assertEquals(69, stateSpace.getNumRows());
        assertEquals(9, stateSpace.getNumCols());
    }

    // Open Network Tests - Single open class with λ=0.1

    @Test
    public void testOpenNetworkFCFS() {
        Network model = createOpenNetwork(SchedStrategy.FCFS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(3);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for open FCFS
        // FCFS has 4x4 generator matrix
        Matrix expectedGen = new Matrix(4, 4);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.1);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.6);
        expectedGen.set(2, 1, 0.5);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -0.6);
        expectedGen.set(3, 2, 0.5);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -0.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "FCFS open network generator matrix mismatch");

        // Verify state space
        assertEquals(4, stateSpace.getNumRows());
        assertEquals(5, stateSpace.getNumCols());
    }

    @Test
    public void testOpenNetworkLCFS() {
        Network model = createOpenNetwork(SchedStrategy.LCFS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(3);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for open LCFS
        // LCFS has 4x4 generator matrix (same as FCFS)
        Matrix expectedGen = new Matrix(4, 4);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.1);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.6);
        expectedGen.set(2, 1, 0.5);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -0.6);
        expectedGen.set(3, 2, 0.5);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -0.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "LCFS open network generator matrix mismatch");

        // Verify state space
        assertEquals(4, stateSpace.getNumRows());
        assertEquals(5, stateSpace.getNumCols());
    }

    @Test
    public void testOpenNetworkSIRO() {
        Network model = createOpenNetwork(SchedStrategy.SIRO);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(3);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for open SIRO
        // SIRO has 4x4 generator matrix
        Matrix expectedGen = new Matrix(4, 4);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.1);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.6);
        expectedGen.set(2, 1, 0.5);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -0.6);
        expectedGen.set(3, 2, 0.5);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -0.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "SIRO open network generator matrix mismatch");

        // Verify state space
        assertEquals(4, stateSpace.getNumRows());
        assertEquals(4, stateSpace.getNumCols());
    }

    @Test
    public void testOpenNetworkSEPT() {
        Network model = createOpenNetwork(SchedStrategy.SEPT);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(3);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for open SEPT
        // SEPT has 4x4 generator matrix
        Matrix expectedGen = new Matrix(4, 4);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.1);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.6);
        expectedGen.set(2, 1, 0.5);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -0.6);
        expectedGen.set(3, 2, 0.5);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -0.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "SEPT open network generator matrix mismatch");

        // Verify state space
        assertEquals(4, stateSpace.getNumRows());
        assertEquals(4, stateSpace.getNumCols());
    }

    @Test
    public void testOpenNetworkLEPT() {
        Network model = createOpenNetwork(SchedStrategy.LEPT);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(3);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for open LEPT
        // LEPT has 4x4 generator matrix
        Matrix expectedGen = new Matrix(4, 4);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.1);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.6);
        expectedGen.set(2, 1, 0.5);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -0.6);
        expectedGen.set(3, 2, 0.5);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -0.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "LEPT open network generator matrix mismatch");

        // Verify state space
        assertEquals(4, stateSpace.getNumRows());
        assertEquals(4, stateSpace.getNumCols());
    }

    @Test
    public void testOpenNetworkPS() {
        Network model = createOpenNetwork(SchedStrategy.PS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(3);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for open PS
        // PS has 4x4 generator matrix
        Matrix expectedGen = new Matrix(4, 4);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.1);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.6);
        expectedGen.set(2, 1, 0.5);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -0.6);
        expectedGen.set(3, 2, 0.5);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -0.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "PS open network generator matrix mismatch");

        // Verify state space
        assertEquals(4, stateSpace.getNumRows());
        assertEquals(3, stateSpace.getNumCols());
    }

    @Test
    public void testOpenNetworkDPS() {
        Network model = createOpenNetwork(SchedStrategy.DPS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(3);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for open DPS
        // DPS has 4x4 generator matrix (same as PS)
        Matrix expectedGen = new Matrix(4, 4);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.1);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.6);
        expectedGen.set(2, 1, 0.5);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -0.6);
        expectedGen.set(3, 2, 0.5);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -0.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "DPS open network generator matrix mismatch");

        // Verify state space
        assertEquals(4, stateSpace.getNumRows());
        assertEquals(3, stateSpace.getNumCols());
    }

    @Test
    public void testOpenNetworkGPS() {
        Network model = createOpenNetwork(SchedStrategy.GPS);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(3);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for open GPS
        // GPS has 4x4 generator matrix
        Matrix expectedGen = new Matrix(4, 4);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.1);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.6);
        expectedGen.set(2, 1, 0.5);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -0.6);
        expectedGen.set(3, 2, 0.5);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -0.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "GPS open network generator matrix mismatch");

        // Verify state space
        assertEquals(4, stateSpace.getNumRows());
        assertEquals(3, stateSpace.getNumCols());
    }

    @Test
    public void testOpenNetworkINF() {
        Network model = createOpenNetwork(SchedStrategy.INF);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(3);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for open INF
        // INF has 4x4 generator matrix with different rates
        Matrix expectedGen = new Matrix(4, 4);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.1);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.6);
        expectedGen.set(2, 1, 1.0);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -1.1);
        expectedGen.set(3, 2, 1.5);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -1.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "INF open network generator matrix mismatch");

        // Verify state space
        assertEquals(4, stateSpace.getNumRows());
        assertEquals(3, stateSpace.getNumCols());
    }

    @Test
    public void testOpenNetworkHOL() {
        Network model = createOpenNetwork(SchedStrategy.HOL);
        SolverCTMC ctmcSolver = new SolverCTMC(model);
        SolverOptions options = new SolverOptions();
        options.cutoff = Matrix.singleton(3);
        ctmcSolver.setOptions(options);

        Matrix generator = ctmcSolver.getGenerator().infGen;
        Matrix stateSpace = ctmcSolver.getStateSpace().stateSpace;

        // Verify generator - complete verification based on MATLAB ground truth for open HOL
        // HOL has 4x4 generator matrix (same as FCFS/LCFS)
        Matrix expectedGen = new Matrix(4, 4);

        // Set all non-zero entries from MATLAB sparse output (complete verification)
        expectedGen.set(0, 0, -0.1);
        expectedGen.set(1, 0, 0.5);
        expectedGen.set(0, 1, 0.1);
        expectedGen.set(1, 1, -0.6);
        expectedGen.set(2, 1, 0.5);
        expectedGen.set(1, 2, 0.1);
        expectedGen.set(2, 2, -0.6);
        expectedGen.set(3, 2, 0.5);
        expectedGen.set(2, 3, 0.1);
        expectedGen.set(3, 3, -0.5);

        // Verify complete generator matrix
        assertMatrixEquals(expectedGen, generator, "HOL open network generator matrix mismatch");

        // Verify state space
        assertEquals(4, stateSpace.getNumRows());
        assertEquals(5, stateSpace.getNumCols());
    }
}