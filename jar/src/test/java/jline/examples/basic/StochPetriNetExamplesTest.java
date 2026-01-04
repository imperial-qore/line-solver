package jline.examples.basic;

import jline.examples.java.basic.StochPetriNetModel;
import jline.lang.Network;
import jline.solvers.NetworkAvgTable;
import jline.solvers.jmt.SolverJMT;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.assertTableMetrics;
import static jline.TestTools.withSuppressedOutput;

//GC

/**
 * Tests for StochPetriNet examples comparing solver outputs between MATLAB and Java implementations.
 *
 */
public class StochPetriNetExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
    }
    

    @Test
    public void testSpnBasicClosedJMT() {
        Network model = StochPetriNetModel.spn_basic_closed();
        
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000, "keep", 2, "verbose", 1, "cutoff", 10);
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        assertNotNull(avgTable);
        
        // Expected values from allExamplesBaseline.txt (JMT solver)
        // Order: P1(Class1)
        double[] expectedQLen = {1.0};
        double[] expectedUtil = {0.0};
        double[] expectedRespT = {0.353082760022072};
        double[] expectedResidT = {0.353082760022072};
        double[] expectedArvR = {2.83497450383065};
        double[] expectedTput = {2.83219718781367};
        
        // Verify table size
        assertEquals(1, avgTable.getQLen().size(),
            "Expected 1 entry (1 station × 1 class)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testSpnBasicOpenJMT() {
        Network model = StochPetriNetModel.spn_basic_open();
        
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000, "keep", 2, "verbose", 1, "cutoff", 10);
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        assertNotNull(avgTable);
        
        // Expected values from allExamplesBaseline.txt (JMT solver)
        // Order: Source(Class1), P1(Class1)
        double[] expectedQLen = {0.0, 0.25965595467563};
        double[] expectedUtil = {0.0, 0.0};
        double[] expectedRespT = {0.0, 0.250296313993243};
        double[] expectedResidT = {0.0, 0.250296313993243};
        double[] expectedArvR = {0.0, 1.03117795734301};
        double[] expectedTput = {1.03117795734301, 1.04270293942227};
        
        // Verify table size
        assertEquals(2, avgTable.getQLen().size(),
            "Expected 2 entries (2 stations × 1 class)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testSpnTwoModesJMT() {
        Network model = StochPetriNetModel.spn_twomodes();

        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000, "keep", 2, "verbose", 1, "cutoff", 10);
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];

        assertNotNull(avgTable);
        
        // Expected values from allExamplesBaseline.txt (JMT solver)
        // Order: P1(Class1), P2(Class1), T1(Class1), T2(Class1)
        double[] expectedQLen = {4.66323597507685, 5.33676402492315};
        double[] expectedUtil = {0.0, 0.0};
        double[] expectedRespT = {0.897636905509251, 1.04742169058208};
        double[] expectedResidT = {0.897636905509251, 1.04742169058208};
        double[] expectedArvR = {5.03479972449197, 4.99912582099244};
        double[] expectedTput = {4.99912582099244, 5.03566478091449};
        
        // Verify table size
        assertEquals(2, avgTable.getQLen().size(),
            "Expected 2 entries (2 stations × 1 class)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testSpnFourModesJMT() {
        Network model = StochPetriNetModel.spn_fourmodes();

        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000, "keep", 2, "verbose", 1, "cutoff", 10);
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];

        assertNotNull(avgTable);
        
        // Expected values from allExamplesBaseline.txt (JMT solver) - only 3 entries shown
        double[] expectedQLen = {3.43438201097283, 2.55889058194375, 1.9935067407945};
        double[] expectedUtil = {0.0, 0.0, 0.0};
        double[] expectedRespT = {0.720317832148919, 0.866494597491938, 1.14047266825276};
        double[] expectedResidT = {0.720317832148919, 0.866494597491938, 1.14047266825276};
        double[] expectedArvR = {4.79603919236671, 2.99838910920584, 1.80443712939947};
        double[] expectedTput = {4.78663240627569, 3.00063031856146, 1.74204760405139};
        
        // Verify table size
        assertEquals(3, avgTable.getQLen().size(),
            "Expected 3 entries (3 stations × 1 class)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testSpnInhibitingJMT() {
        Network model = StochPetriNetModel.spn_inhibiting();

        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000, "keep", 2, "verbose", 1, "cutoff", 10);
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];

        assertNotNull(avgTable);
        
        // Expected values from allExamplesBaseline.txt (JMT solver)
        double[] expectedQLen = {1.36667057262832, 0.744957542362553, 1.8402729234593};
        double[] expectedUtil = {0.0, 0.0, 0.0};
        double[] expectedRespT = {0.568015067355994, 0.437092905729172, 2.63513187817813};
        double[] expectedResidT = {0.568015067355994, 0.109273226432293, 1.31756593908906};
        double[] expectedArvR = {2.47866990307651, 1.73626408780108, 0.71579876920703};
        double[] expectedTput = {2.47781363103459, 1.71340935053271, 0.71633948532544};
        
        // Verify table size
        assertEquals(3, avgTable.getQLen().size(),
            "Expected 3 entries (3 stations × 1 class)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testSpnClosedTwoPlacesJMT() {
        Network model = StochPetriNetModel.spn_closed_twoplaces();

        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000, "keep", 2, "verbose", 1, "cutoff", 10);
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];

        assertNotNull(avgTable);
        
        // Expected values from allExamplesBaseline.txt (JMT solver)
        // Order: P1(Class1,Class2), P2(Class1,Class2), T1(Class1,Class2), T2(Class1,Class2), T3(Class1,Class2)
        double[] expectedQLen = {0.881093703530199, 1.70796452117085, 9.1189062964698, 5.29203547882915};
        double[] expectedUtil = {0.0, 0.0, 0.0, 0.0};
        double[] expectedRespT = {1.16878177655751, 1.11616282064011, 12.100452180648, 3.34226767446581};
        double[] expectedResidT = {1.16878177655751, 1.11616282064011, 12.100452180648, 3.34226767446581};
        double[] expectedArvR = {0.753295516771666, 1.59055256845185, 0.75487566926272, 1.60422257322269};
        double[] expectedTput = {0.75487566926272, 1.60422257322269, 0.754954494671509, 1.58207453336418};
        
        // Verify table size
        assertEquals(4, avgTable.getQLen().size(),
            "Expected 4 entries (2 stations × 2 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testSpnClosedFourPlacesJMT() {
        Network model = StochPetriNetModel.spn_closed_fourplaces();
        
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000, "keep", 2, "verbose", 1, "cutoff", 10);
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        assertNotNull(avgTable);
        
        // Expected values from allExamplesBaseline.txt
        // Order: P1(Class1), P2(Class1), P3(Class1), P4(Class1)
        double[] expectedQLen = {0.277478994982809, 0.776764015466386, 0.251642290175409, 0.693864139928483};
        double[] expectedUtil = {0.0, 0.0, 0.0, 0.0};
        double[] expectedRespT = {0.482170727708879, 1.34410606292658, 0.430353189746762, 1.2300544147064};
        double[] expectedResidT = {0.482170727708879, 1.34410606292658, 0.430353189746762, 1.2300544147064};
        double[] expectedArvR = {0.574608661902354, 0.574324747516291, 0.572727382379681, 0.572936886131921};
        double[] expectedTput = {0.574324747516291, 0.572727382379681, 0.572936886131921, 0.574400657885676};
        
        // Verify table size
        assertEquals(4, avgTable.getQLen().size(), 
            "Expected 4 entries (4 stations × 1 class)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testSpnOpenSevenPlacesJMT() {
        Network model = StochPetriNetModel.spn_open_sevenplaces();
        
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000, "keep", 2, "verbose", 1, "cutoff", 10);
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        assertNotNull(avgTable);
        
        // Expected values from MATLAB JMT solver reference implementation
        // Complex open model with 8 entries (Source, P1-P7)
        double[] expectedQLen = {0.0, 0.7350410943010985, 0.0, 0.14403778288135946, 0.5036622428109216, 0.7904727711121464, 1.3831613557669105, 0.21019992817868244};
        double[] expectedUtil = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double[] expectedRespT = {0.0, 0.2533740304227246, 0.0, 0.10278794849508507, 0.3591654216997956, 0.4329319389151261, 0.9896536517663469, 0.5007662819992489};
        double[] expectedResidT = {0.0, 0.487979614147469, 0.0, 0.098980987439711, 0.172931499336938, 0.288621292610084, 0.403192228497400, 0.778969771998831};
        double[] expectedArvR = {0.0, 2.8389053217941407, 2.8375993575458270, 1.4098436448743297, 1.4223165421477546, 1.8031136508818084, 1.4058170835170205, 0.4175884206234671};
        double[] expectedTput = {1.0152620184538280, 2.8375993575458270, 2.8375993575458270, 1.4058170835170205, 1.4232209203641903, 1.8065763625473912, 1.4005668952229813, 0.4193193304556456};
        
        // Verify table size
        assertEquals(8, avgTable.getQLen().size(), 
            "Expected 8 entries (8 stations × 1 class)");
        
        // Check all metrics against expected values (structural validation for complex model)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
}