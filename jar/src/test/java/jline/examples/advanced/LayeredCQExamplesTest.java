package jline.examples.advanced;

import jline.examples.java.advanced.LayeredCQModel;
import jline.lang.Network;
import jline.lang.constant.SolverType;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.ln.SolverLN;
import jline.solvers.AvgTable;
import jline.solvers.LayeredNetworkAvgTable;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.*;

/**
 * Test class for LayeredCQExamples examples.
 *
 * Expected values come from allExamplesBaseline.txt MATLAB execution.
 */
public class LayeredCQExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
    }


    @Test
    public void testLcqSinglehost() {
        // Create the model
        LayeredNetwork model = LayeredCQModel.lcq_singlehost();

        // Create and run the solver
        SolverLN solver = new SolverLN(model, SolverType.MVA);
        AvgTable avgTable = solver.getAvgTable();

        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);

        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;

        // Expected values from allExamplesBaseline.txt
        // Order: P1, PC, T1, C2, E1, I2, A1, AC2, AC2h, AC2m
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.333333322222223, 0.666666644444445};
        double[] expectedUtil = {0.0, 1.0, 0.0, 1.0, Double.NaN, Double.NaN, 0.0, 0.0, 0.333333322222223, 0.666666644444445};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.5, 1.5, 1.5, 0.0, 1.0, 2.0};
        double[] expectedResidT = {Double.NaN, Double.NaN, 0.0, 1.5, Double.NaN, Double.NaN, 0.0, 0.0, 0.5, 1.0};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.666666653333334, 0.666666644444445, 0.666666653333334, 0.666666644444445, 0.666666653333334, 0.0, 0.333333322222223, 0.333333322222223};
        
        // Verify table size
        assertEquals(10, lnAvgTable.getQLen().size(), 
            "Expected 10 entries as shown in allExamplesBaseline.txt");
        
        // Check all metrics against expected values
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testLcqThreehosts() {
        // Create the model
        LayeredNetwork model = LayeredCQModel.lcq_threehosts();

        // Create and run the solver
        SolverLN solver = new SolverLN(model, SolverType.MVA);
        AvgTable avgTable = solver.getAvgTable();

        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);

        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;

        // Expected values for lcq_threehosts
        // Order: P1, Pc, P2, T1, CT, T2, E1, IE, E2, A2, A1, Ac, Ac_hit, Ac_miss (14 entries)
        // 3 processors + 3 tasks + 3 entries + 5 activities = 14

        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN,
                                1.0, 1.0, 0.0625,
                                1.0, 1.0, 0.0625,
                                0.0625, 1.0, 0.0, 0.3125, 0.6875};
        double[] expectedUtil = {0.0, 0.9375, 0.0625,
                                0.0, 0.9375, 0.0625,
                                Double.NaN, Double.NaN, Double.NaN,
                                0.0625, 0.0, 0.0, 0.3125, 0.625};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN,
                                 Double.NaN, Double.NaN, Double.NaN,
                                 1.6, 1.6, 0.2,
                                 0.2, 1.6, 0.0, 1.0, 2.2};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN,
                                  0.0, 1.5, 0.2,
                                  Double.NaN, Double.NaN, Double.NaN,
                                  0.2, 0.0, 0.0, 0.5, 1.0};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN,
                                Double.NaN, Double.NaN, Double.NaN,
                                Double.NaN, Double.NaN, Double.NaN,
                                Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN,
                                0.625, 0.625, 0.3125,
                                0.625, 0.625, 0.3125,
                                0.3125, 0.625, 0.0, 0.3125, 0.3125};

        // Verify table size
        assertEquals(14, lnAvgTable.getQLen().size(),
            "Expected 14 entries (3 processors + 3 tasks + 3 entries + 5 activities)");
        
        // Check all metrics against expected values
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
}
