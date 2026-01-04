package jline.examples.basic;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.examples.java.basic.ClassSwitchingModel;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.solvers.NetworkAvgChainTable;
import jline.solvers.mva.SolverMVA;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.assertTableMetrics;

/**
 * Unit tests for ClassSwitch examples.
 *
 * ANNOTATION: ClassSwitch examples are not present in allExamplesBaseline.txt.
 * Current expected values are based on Java implementation baseline.
 * The Java implementation creates implicit nodes (Sink, ClassSwitch) that may not exist in MATLAB output.
 */
public class ClassSwitchExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        // Set verbose level to SILENT to suppress warnings during tests
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }
    

    // ===== cs_implicit tests =====

    //@Disabled
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testCsImplicitMVA() {
        // Test the cs_implicit example with MVA solver
        Network model = ClassSwitchingModel.cs_implicit();

        SolverMVA solver = new SolverMVA(model, "verbose", VerboseLevel.STD);
        NetworkAvgChainTable avgTable = solver.getAvgChainTable();

        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB output (MVA solver)
        // MATLAB shows chain-based output: Queue 1 (chain) QLen=1.49984049264051, Util=0.6, Tput=0.6
        // Since Java filters zeros and the test expects specific structure, we match the filtered output
        // Only Queue 1 has non-zero values in MATLAB
        double[] expectedQLen = {0.0,1.5};
        double[] expectedUtil = {0.0,0.6};
        double[] expectedRespT = {0.0,2.5};
        double[] expectedResidT = {0.0,2.5};
        double[] expectedArvR = {0.0,0.6};
        double[] expectedTput = {0.6,0.6};
        
        // Verify table size - MATLAB shows only Queue 1 with non-zero values
        assertEquals(2, avgTable.getQLen().size(), "Expected 2 entry matching MATLAB output");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    // ===== cs_multi_diamond tests =====
    @Test
    public void testCsMultiDiamondMVA() {
        // Test the cs_multi_diamond example with MVA solver
        Network model = ClassSwitchingModel.cs_multi_diamond();
        
        SolverMVA solver = new SolverMVA(model, "verbose", VerboseLevel.STD);
        NetworkAvgChainTable avgTable = solver.getAvgChainTable();
        //avgTable.print();

        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        double[] expectedQLen = {0.0, 0.142856887163414, 0.0191082767626744, 0.0212765916692491};
        double[] expectedUtil = {0.0, 0.125, 0.01875, 0.0208333333333333};
        double[] expectedRespT = {0.0, 0.114285509730732, 0.0509554047004651, 0.0340425466707986};
        double[] expectedResidT = {0.0, 0.142856887163414, 0.0191082767626744, 0.0212765916692491};
        double[] expectedArvR = {0.0, 1.25, 0.375, 0.625};
        double[] expectedTput = {1.0, 1.25, 0.375, 0.625};
        
        // Verify table size - MATLAB shows 4 stations
        assertEquals(4, avgTable.getQLen().size(), "Expected 4 entries matching MATLAB output");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    

    // ===== cs_single_diamond tests =====

    @Test
    public void testCsSingleDiamondMVA() {
        // Test the cs_single_diamond example with MVA solver
        Network model = ClassSwitchingModel.cs_single_diamond();
        
        SolverMVA solver = new SolverMVA(model, "verbose", VerboseLevel.STD);
        NetworkAvgChainTable avgTable = solver.getAvgChainTable();

        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/exact", solver.result.method, "MVA solver should use default/exact method");
        //avgTable.print();

        double[] expectedQLen = {0.32258064516129, 0.193548387096774, 0.483870967741935};
        double[] expectedUtil = {0.32258064516129, 0.193548387096774, 0.483870967741935};
        double[] expectedRespT = {1.0, 2.0, 3.0};
        double[] expectedResidT = {1.0, 0.6, 1.5};
        double[] expectedArvR = {0.32258064516129, 0.0967741935483871, 0.161290322580645};
        double[] expectedTput = {0.32258064516129, 0.0967741935483871, 0.161290322580645};
        
        // Verify table size - MATLAB shows 3 queues
        assertEquals(3, avgTable.getQLen().size(), "Expected 3 entries matching MATLAB output (Queue 0, 1, 2)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    // ===== cs_transient_class tests =====
    // This model has a transient class (Class1 -> Class2/Class3) with reducible routing chain.
    // The dtmc_solve_reducible function properly handles transient states.

    @Disabled
    @Test
    public void testCsTransientClassMVA() {
        // Test the cs_transient_class example with MVA solver
        Network model = ClassSwitchingModel.cs_transient_class();

        SolverMVA solver = new SolverMVA(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgChainTable avgTable = solver.getAvgChainTable();

        // Check if results are computed
        assertNotNull(avgTable);

        // Verify solver method - uses exact method for closed networks
        assertEquals("default/exact", solver.result.method, "MVA solver should use default/exact method");

        // Expected values for transient class model:
        // Class1 is transient (jobs eventually move to Class2 or Class3)
        // In steady state: ~45.83% in Class2 cycle, ~54.17% in Class3 cycle
        // Queue 0: visited by all classes, Queue 1: Class2 only, Queue 2: Class3 only
        double[] expectedQLen = {0.5, 0.229166666666667, 0.270833333333333};
        double[] expectedUtil = {0.5, 0.229166666666667, 0.270833333333333};
        double[] expectedRespT = {1.0, 1.0, 1.0};
        double[] expectedResidT = {1.0, 0.458333333333333, 0.541666666666667};
        double[] expectedArvR = {0.5, 0.229166666666667, 0.270833333333333};
        double[] expectedTput = {0.5, 0.229166666666667, 0.270833333333333};

        // Verify table size - 3 queues
        assertEquals(3, avgTable.getQLen().size(), "Expected 3 entries (Queue 0, 1, 2)");

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }
}