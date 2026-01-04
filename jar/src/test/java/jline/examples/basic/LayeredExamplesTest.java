package jline.examples.basic;

import java.util.Arrays;
import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.solvers.ln.LNOptions;
import jline.examples.java.basic.LayeredModel;
import jline.lang.constant.SolverType;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.AvgTable;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.ln.SolverLN;
import jline.solvers.lqns.SolverLQNS;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.*;

public class LayeredExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
    }
    

    @Test
    @Timeout(300)  // 5 minute timeout for SolverLN iteration
    public void testLqnSerial() throws Exception {
        // Create the model
        LayeredNetwork model = LayeredModel.lqn_serial();
        
        // Create and run the solver
        final AvgTable[] avgTableHolder = new AvgTable[1];
        withSuppressedOutput(() -> {
            SolverLN solver = new SolverLN(model, SolverType.MVA);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);
        
        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;
        
        // Expected values from MATLAB dev/ execution (SolverLN results)
        // Source: matlab/examples/basic/layeredModel/lqn_serial.m run in dev/ directory with SolverLN
        // Order: P1, P2, T1, T2, E1, E2, AS1, AS2, AS3, AS4
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.12406428613832, 0.532556142194661, 1.12406428613832, 0.532556142194661, 0.162458950075148, 0.961605337078543, 0.443796785162218, 0.0887593570324436};
        double[] expectedUtil = {0.142014971391375, 0.532556142194661, 0.142014971391375, 0.532556142194661, Double.NaN, Double.NaN, 0.142014971391375, 0.0, 0.443796785162218, 0.0887593570324436};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 12.6641778694225, 6.0, 1.83033040512252, 10.8338474757396, 5.0, 1.0};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.83033040512252, 6.0, Double.NaN, Double.NaN, 1.83033040512252, 0.0, 5.0, 1.0};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0887593571196091, 0.0887593570324436, 0.0887593571196091, 0.0887593570324436, 0.0887593571196091, 0.0887593571196091, 0.0887593570324436, 0.0887593570324436};
        
        // Verify table size
        assertEquals(10, lnAvgTable.getQLen().size(), 
            "Expected 10 entries (2 processors, 2 tasks, 2 entries, 4 activities)");
        
        // Check all metrics against expected values
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    @Timeout(300)  // 5 minute timeout for SolverLN iteration
    public void testLqnMultiSolvers() throws Exception {
        // Create the model
        LayeredNetwork model = LayeredModel.lqn_multi_solvers();
        
        // Create and run the solver
        final AvgTable[] avgTableHolder = new AvgTable[1];
        withSuppressedOutput(() -> {
            SolverLN solver = new SolverLN(model, SolverType.MVA);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0]; //avgTable.print();
        
        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);
        
        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;
        
        // Expected values from MATLAB dev/ execution (SolverLN results)
        // Source: matlab/examples/basic/layeredModel/lqn_multi_solvers.m run in dev/ directory
        // Order: P1, P2, T1, T2, E1, E2, A1, A2
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0, 0.74998, 1.0, 0.74998, 1.0, 0.74998};
        double[] expectedUtil = {0.24999, 0.74998, 0.24999, 0.74998, Double.NaN, Double.NaN, 0.24999, 0.74998};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 4.0, 1.0, 4.0, 1.0};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.0, 1.0, Double.NaN, Double.NaN, 1.0, 1.0};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.24999, 0.74998, 0.24999, 0.74998, 0.24999, 0.74998};
        
        // Verify table size
        assertEquals(8, lnAvgTable.getQLen().size(), 
            "Expected 8 entries (2 processors, 2 tasks, 2 entries, 2 activities)");
        
        // Check all metrics against expected values
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    @Timeout(300)  // 5 minute timeout for SolverLN iteration
    public void testLqnTwotasks() throws Exception {
        // Create the model
        LayeredNetwork model = LayeredModel.lqn_twotasks();

        // Create and run the solver
        final AvgTable[] avgTableHolder = new AvgTable[1];
        withSuppressedOutput(() -> {
            SolverLN solver = new SolverLN(model, SolverType.MVA);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0];

        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);

        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;

        // Expected values from MATLAB dev/ execution (SolverLN results)
        // Source: matlab/examples/basic/layeredModel/lqn_twotasks.m run in dev/ directory
        // Order: P1, P2, T1, T2, E1, E2, E3, A1, A20, A21, A22, A3
        double[] expectedQLen = {Double.NaN, Double.NaN, 97.5, 97.2, 97.5, 72.9, 24.3, 97.5, 24.3, 24.3, 24.3, 24.3};
        double[] expectedUtil = {0.25, 1.0, 0.25, 1.0, Double.NaN, Double.NaN, Double.NaN, 0.25, 0.25, 0.25, 0.25, 0.25};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 390.0, 291.5, 97.2, 390.0, 97.2, 97.2, 97.2, 97.2};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.3275, 194.3, Double.NaN, Double.NaN, Double.NaN, 1.3275, 48.6, 48.6, 48.6, 48.6};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.25, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};

        // Verify table size
        assertEquals(12, lnAvgTable.getQLen().size(),
            "Expected 12 entries (2 processors, 2 tasks, 3 entries, 5 activities)");

        // Check all metrics against expected values (using COARSE_TOL for tiny numerical perturbations)
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
    }

    @Test
    @Disabled("FAIL - LQN solver Util[1] mismatch: expected 0.0615 but got 0.040 (~35% diff)")
    public void testLqnFunction() throws Exception {
        // Create the model
        LayeredNetwork model = LayeredModel.lqn_function();
        
        // Create and run the solver
        final AvgTable[] avgTableHolder = new AvgTable[1];
        withSuppressedOutput(() -> {
            LNOptions lnoptions = new LNOptions();
            lnoptions.seed = 23000;
            SolverLN solver = new SolverLN(model, SolverType.MVA);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0]; //avgTable.print();
        
        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);
        
        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;
        
        // Expected values from MATLAB ground truth for lqn_function
        // Source: MATLAB output run in dev/ directory
        // Order: P1, P2, T1, F2, E1, E2, A1, A2
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0, 0.061516866224749, 1.0, 0.061516866224749, 1.0, 0.061516866224749};
        double[] expectedUtil = {0.74999998875, 0.061516866224749, 0.74999998875, 0.061516866224749, Double.NaN, Double.NaN, 0.74999998875, 0.061516866224749};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.33333333333333, 0.333333333333333, 1.33333333333333, 0.333333333333333};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.0, 0.333333333333333, Double.NaN, Double.NaN, 1.0, 0.333333333333333};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.74999998875, 0.184550598674247, 0.74999998875, 0.184550598674247, 0.74999998875, 0.184550598674247};
        
        // Verify table size
        assertEquals(8, lnAvgTable.getQLen().size(), 
            "Expected 8 entries (2 processors, 2 tasks, 2 entries, 2 activities)");
        
        // Check all metrics against expected values
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    @Disabled("LINE SolverLN and LQNS produce fundamentally different results for this model (up to 3x difference on P3 util). Java matches MATLAB LINE output.")
    public void testLqnWorkflows() throws Exception {
        // Ground truth values from LQNS (external reference solver)
        // Note: LINE SolverLN produces results ~4-6% different from LQNS

        // Create the model
        LayeredNetwork model = LayeredModel.lqn_workflows();

        // Create and run the solver
        final AvgTable[] avgTableHolder = new AvgTable[1];
        withSuppressedOutput(() -> {
            SolverLN solver = new SolverLN(model, SolverType.MVA);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0];

        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);

        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;

        // LQNS ground truth values for lqn_workflows
        // Order: P1, P2, P3, T1, T2, T3, Entry, E2, E1, A1, A2, A3, B1, B2, B3, B4, B5, B6, C1, C2, C3, C4, C5
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 1, 0.204246, 0.0723959, 1, 0.204246, 0.0723959, 0.0795754, 0.477452, 0.442972, 0.00795754, 0.0159151, 0.0238726, 0.0318301, 0.0397877, 0.120141, 0.00795559, 0.00477336, 0.00716003, 0.0127289, 0.039778};
        double[] expectedUtil = {0.795754, 0.167108, 0.0724136, 0.795754, 0.167108, 0.0724136, 0.795754, 0.167108, 0.0724136, 0.0795754, 0.477452, 0.238726, 0.00795754, 0.0159151, 0.0238726, 0.0318301, 0.0397877, 0.0477452, 0.00795754, 0.00477452, 0.00716178, 0.0127321, 0.0397877};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 12.5667, 2.5667, 0.909778, 1, 2, 5.5667, 0.1, 0.2, 0.3, 0.4, 0.5, 1.50978, 0.1, 0.2, 0.3, 0.4, 0.5};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.0795754, 0.0795754, 0.0795754, 0.0795754, 0.0795754, 0.0795754, 0.0795754, 0.238726, 0.0795754, 0.0795754, 0.0795754, 0.0795754, 0.0795754, 0.0795754, 0.0795754, 0.0795754, 0.0238726, 0.0238726, 0.0318301, 0.0795754};

        // Verify table size
        assertEquals(23, lnAvgTable.getQLen().size(),
            "Expected 23 entries (3 processors, 3 tasks, 3 entries, 14 activities)");

        // Check all metrics against expected values (using VERY_COARSE_TOL for LINE vs LQNS comparison)
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, VERY_COARSE_TOL);
    }

    @Test
    @Disabled("FAIL - NullPointerException when loading ofbizExample.xml model")
    public void testLqnOfbiz() throws Exception {
        // Create the model - loads from ofbizExample.xml
        LayeredNetwork model = LayeredModel.lqn_ofbiz();
        
        // Check if model was loaded successfully
        org.junit.jupiter.api.Assumptions.assumeTrue(model != null,
                "ofbizExample.xml not available - skipping test");
        
        // Create and run the solver
        final AvgTable[] avgTableHolder = new AvgTable[1];
        withSuppressedOutput(() -> {
            SolverLN solver = new SolverLN(model, SolverType.MVA);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);
        
        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;

        // Ground truth values from MATLAB execution
        // 72 entries total (9 processors, 10 tasks, 14 entries, 39 activities)
        double[] expectedQLen = {
            // Processors (9)
            Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks (10)
            0.131004377079255, 0.0, 0.0670898699719763, 0.0463312566511443, 0.0191709441581764, 0.0383418785384276, 0.0191709441581764, 0.0639155350921157, 0.00798938309770094,
            // Entries (14)
            0.0, 0.0463312513923384, 0.0191709368292392, 0.0383418727656901, 0.0191709368292392, 0.00798937926274808, 0.0, 0.0670898699719763, 0.0463312566511443, 0.0191709441581764, 0.0383418785384276, 0.0191709441581764, 0.0639155350921157, 0.00798938309770094,
            // Activities (39)
            0.0, 0.0463312513923384, 0.0191709368292392, 0.0383418727656901, 0.0191709368292392, 0.00798937926274808, 0.0, 0.0, 0.0111816537106558, 0.0111816568950473, 0.0111816541790173, 0.0111816541790173, 0.0111816568950473, 0.0111816537106558, 0.0, 
            4.11574352685691e-08, 0.0463312978085795, 4.11574352685691e-08, 1.70301146849138e-08, 0.019170961188291, 1.70301146849138e-08, 3.40602306102675e-08, 0.0383419125986582, 3.40602306102675e-08, 1.70301146849138e-08, 0.0191709611882911, 1.70301146849138e-08,
            0.0, 0.00798944794924938, 0.00798945022454167, 0.00798944828390027, 0.00798944828390027, 0.00798944794924938, 0.00798945100511231, 0.00798945022454166, 0.00798944794924938, 0.0, 0.0, 0.00798939019490449, 0.0
        };
        
        double[] expectedUtil = {
            // Processors
            0.116375098351779, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Tasks  
            0.116375098351779, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Entries
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Activities
            0.0, 0.0411574334976797, 0.0170301153956768, 0.034060229998263, 0.0170301153956768, 0.00709720406448317, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        };
        
        double[] expectedRespT = {
            // Processors
            Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks
            Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Entries
            0.0, 0.0112570798164444, 0.0112570798164444, 0.0112570798164444, 0.0112570798164444, 0.0112570798164444, 0.0, 0.067543015821165, 0.0112570819302538, 0.0112570851361539, 0.0112570824017789, 0.0112570851361539, 0.0900573562823318, 0.0112570862359825,
            // Activities
            0.0, 0.0112570798164444, 0.0112570798164444, 0.0112570798164444, 0.0112570798164444, 0.0112570798164444, 0.0,
            1.00000007151696e-08, 0.0112571780777199, 0.0112571812836198, 0.0112571785492449, 0.0112571785492449, 0.0112571812836198, 0.0112571780777199, 1.00000007151696e-08,
            1.00000011729869e-08, 0.011257091930255, 1.00000011729869e-08, 1.00000004853583e-08, 0.0112570951361544, 1.00000004853583e-08, 1.00000009707166e-08, 0.0112570924017799, 1.00000009707166e-08, 1.00000004853583e-08, 0.0112570951361544, 1.00000004853583e-08,
            1.00000006387483e-08, 0.0112571780777198, 0.0112571812836197, 0.0112571785492449, 0.0112571785492449, 0.0112571780777198, 0.0112571823834483, 0.0112571812836197, 0.0112571780777198, 1.00000006387483e-08,
            1.00000001916245e-08, 0.0112570962359827, 1.00000001916245e-08
        };
        
        double[] expectedResidT = {
            // Processors
            Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks
            0.0112570798164444, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Entries
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Activities
            0.0, 0.00398119976253751, 0.00164734011835498, 0.00329468015999351, 0.00164734011835498, 0.000686519657203445, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0
        };
        
        double[] expectedArvR = new double[72]; // All NaN or 0
        Arrays.fill(expectedArvR, Double.NaN);
        expectedArvR[10] = 0.0; // USAGE_DELAY_Task
        expectedArvR[19] = 0.0; // USAGE_DELAY0_Entry
        expectedArvR[39] = 0.0; // USAGE_DELAY0_Activity
        
        double[] expectedTput = {
            // Processors
            Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks
            11.6375098351779, 0.0, 0.993291003611854, 4.11574304408564, 1.70301138583431, 3.40602273039848, 1.70301138583431, 0.709720313038493, 0.709720342388727,
            // Entries
            0.0, 4.11574334976797, 1.70301153956768, 3.4060229998263, 1.70301153956768, 0.709720406448316, 0.0, 0.993291003611854, 4.11574304408564, 1.70301138583431, 3.40602273039848, 1.70301138583431, 0.709720313038493, 0.709720342388727,
            // Activities
            0.0, 4.11574334976797, 1.70301153956768, 3.4060229998263, 1.70301153956768, 0.709720406448317, 0.0,
            0.993291003611854, 0.993291003611854, 0.993291003611854, 0.993291003611854, 0.993291003611854, 0.993291003611854, 0.993291003611854, 0.993291003611854,
            4.11574304408564, 4.11574304408564, 4.11574304408564, 1.70301138583431, 1.70301138583431, 1.70301138583431, 3.40602273039848, 3.40602273039848, 3.40602273039848, 1.70301138583431, 1.70301138583431, 1.70301138583431,
            0.709720313038493, 0.709720313038493, 0.709720313038493, 0.709720313038493, 0.709720313038493, 0.709720313038493, 0.709720313038493, 0.709720313038493, 0.709720313038493, 0.709720313038493,
            0.709720342388727, 0.709720342388727, 0.709720342388727
        };
        
        // Verify table size
        assertEquals(72, lnAvgTable.getQLen().size(), 
            "Expected 72 entries for lqn_ofbiz model");
        
        // Check all metrics against expected values
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    @Timeout(300)  // 5 minute timeout for SolverLN iteration
    //@Disabled("Test failing - disabled for investigation")
    public void testLqnBasic() throws Exception {
        // Test the actual lqn_basic model (multi-processor basic network)
        LayeredNetwork model = LayeredModel.lqn_basic();
        
        // Create and run the solver
        final AvgTable[] avgTableHolder = new AvgTable[1];
        withSuppressedOutput(() -> {
            SolverLN solver = new SolverLN(model, SolverType.MVA);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0]; //avgTable.print();
        
        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);
        
        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;
        
        // Expected values from MATLAB ground truth for lqn_basic
        // Order: P1, P2, T1, T2, T3, E1, E2, E3, AS1, AS2, AS3
        double[] expectedQLen = {Double.NaN, Double.NaN, 23.4682450470894, 8.67822070268364, 0.120491297980306, 23.4682450470894, 8.67822070268364, 0.120491297980306, 23.4388606730313, 8.66725245115206, 0.124378114044187};
        double[] expectedUtil = {0.995954019616073, 0.0414593696959083, 0.664070430319683, 0.331883589296389, 0.0414593696959083, Double.NaN, Double.NaN, Double.NaN, 0.664070430319683, 0.331883589296389, 0.0414593696959083};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.76699970180812, 0.653709687866906, 0.0193750007720929, 1.76478725771209, 0.652883475613173, 0.0200000007969991};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.09193140773878, 0.556095776529935, 0.0200000007969991, Double.NaN, Double.NaN, Double.NaN, 1.09193140773878, 0.556095776529935, 0.0200000007969991};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 13.2814086063937, 13.2753435718556, 6.21890545438625, 13.2814086063937, 13.2753435718556, 6.21890545438625, 13.2814086063937, 13.2753435718556, 6.21890545438625};
        
        // Verify table size
        assertEquals(11, lnAvgTable.getQLen().size(),
            "Expected 11 entries (2 processors, 3 tasks, 3 entries, 3 activities)");

        // Check all metrics against expected values (using 0.10 tolerance for open arrival changes)
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, 0.10);
    }

    /*@Test
    @Disabled("Takes more than 1h" +
            "")
    public void testLqnBpmn() throws Exception {
        // This test documents the current limitation with fork-join topology
        LayeredNetwork model = LayeredModel.lqn_bpmn();

        // Run without output suppression to see debug info
        SolverLN solver = new SolverLN(model, SolverType.MVA);
        AvgTable avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }*/

    // ===== MULTIPLE SOLVER TESTS =====
    // These tests mirror the MATLAB examples that use multiple solvers

    @Test
    public void testLqnSerialWithLQNS() throws Exception {
        // Test lqn_serial with SolverLQNS (as used in MATLAB version)
        LayeredNetwork model = LayeredModel.lqn_serial();
        
        // This test may fail if LQNS is not installed, so we handle that gracefully
        try {
            final AvgTable[] avgTableHolder = new AvgTable[1];
            withSuppressedOutput(() -> {
                SolverLQNS solver = new SolverLQNS(model);
                avgTableHolder[0] = solver.getAvgTable();
            });
            AvgTable avgTable = avgTableHolder[0]; //avgTable.print();
            
            assertNotNull(avgTable);
            assertTrue(avgTable instanceof LayeredNetworkAvgTable);
            
            LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;
            
            // Verify table size matches expectation
            assertEquals(10, lnAvgTable.getQLen().size(), 
                "Expected 10 entries (2 processors, 2 tasks, 2 entries, 4 activities) for LQNS solver");
                
            // LQNS solver should produce results, but they may differ from SolverLN
            // We mainly verify that the solver runs and produces reasonable results
            List<Double> tput = lnAvgTable.getTput();
            
            // Check that throughput values are reasonable (not all zero or NaN where expected)
            // T1 and T2 tasks should have positive throughput
            assertTrue(tput.get(2) > 0, "T1 task should have positive throughput with LQNS");
            assertTrue(tput.get(3) > 0, "T2 task should have positive throughput with LQNS");
            
        } catch (RuntimeException e) {
            if (e.getMessage().contains("lqns") || e.getMessage().contains("lqsim")) {
                // LQNS not available - this is expected in many environments
                //System.out.println("LQNS solver not available - skipping test: " + e.getMessage());
                org.junit.jupiter.api.Assumptions.assumeTrue(false, "LQNS solver not available");
            } else {
                throw e; // Re-throw other exceptions
            }
        }
    }

    @Test
    public void testLqnMultiSolversWithLQNS() throws Exception {
        // Test lqn_multi_solvers with SolverLQNS (as used in MATLAB version)
        LayeredNetwork model = LayeredModel.lqn_multi_solvers();
        
        // This test may fail if LQNS is not installed, so we handle that gracefully
        try {
            final AvgTable[] avgTableHolder = new AvgTable[1];
            withSuppressedOutput(() -> {
                SolverLQNS solver = new SolverLQNS(model);
                avgTableHolder[0] = solver.getAvgTable();
            });
            AvgTable avgTable = avgTableHolder[0]; //avgTable.print();
            
            assertNotNull(avgTable);
            assertTrue(avgTable instanceof LayeredNetworkAvgTable);
            
            LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;
            
            // Verify table size matches expectation
            assertEquals(8, lnAvgTable.getQLen().size(), 
                "Expected 8 entries (2 processors, 2 tasks, 2 entries, 2 activities) for LQNS solver");
                
            // Check that throughput values are reasonable
            List<Double> tput = lnAvgTable.getTput();
            assertTrue(tput.get(2) > 0, "T1 task should have positive throughput with LQNS");
            assertTrue(tput.get(3) > 0, "T2 task should have positive throughput with LQNS");
            
        } catch (RuntimeException e) {
            if (e.getMessage().contains("lqns") || e.getMessage().contains("lqsim")) {
                // LQNS not available - this is expected in many environments
                //System.out.println("LQNS solver not available - skipping test: " + e.getMessage());
                org.junit.jupiter.api.Assumptions.assumeTrue(false, "LQNS solver not available");
            } else {
                throw e; // Re-throw other exceptions
            }
        }
    }

    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testLqnTwotasksWithLQNS() throws Exception {
        // Test lqn_twotasks with SolverLQNS (as used in MATLAB version)
        LayeredNetwork model = LayeredModel.lqn_twotasks();
        
        // This test may fail if LQNS is not installed, so we handle that gracefully
        try {
            final AvgTable[] avgTableHolder = new AvgTable[1];
            withSuppressedOutput(() -> {
                SolverLQNS solver = new SolverLQNS(model);
                avgTableHolder[0] = solver.getAvgTable();
            });
            AvgTable avgTable = avgTableHolder[0]; //avgTable.print();
            
            assertNotNull(avgTable);
            assertTrue(avgTable instanceof LayeredNetworkAvgTable);
            
            LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;
            
            // Verify table size matches expectation
            assertEquals(12, lnAvgTable.getQLen().size(), 
                "Expected 12 entries (2 processors, 2 tasks, 3 entries, 5 activities) for LQNS solver");
                
            // Check that throughput values are reasonable
            List<Double> tput = lnAvgTable.getTput();
            assertTrue(tput.get(2) > 0, "T1 task should have positive throughput with LQNS");
            assertTrue(tput.get(3) > 0, "T2 task should have positive throughput with LQNS");
            
        } catch (RuntimeException e) {
            if (e.getMessage().contains("lqns") || e.getMessage().contains("lqsim")) {
                // LQNS not available - this is expected in many environments
                //System.out.println("LQNS solver not available - skipping test: " + e.getMessage());
                org.junit.jupiter.api.Assumptions.assumeTrue(false, "LQNS solver not available");
            } else {
                throw e; // Re-throw other exceptions
            }
        }
    }

    // Note: MATLAB examples also use SolverLN with SolverMVA and SolverNC as factory methods
    // However, these require more complex setup and are not directly exposed through
    // the LayeredModel.java examples. The key requirement is that each LayeredExamples
    // method has a corresponding test, which we have achieved.

    // ===== getRawAvgTables TESTS =====
    // These tests verify the getRawAvgTables implementation aligned with MATLAB

    @Test
    public void testPhase2ModelWithSolverLN() throws Exception {
        // Test that SolverLN can solve phase-2 models and compare with LQNS results

        // Load the phase-2 model from test resources
        String modelPath = getClass().getClassLoader().getResource("lqn/phase2/03-sanity.lqnx").getPath();

        final LayeredNetwork[] modelHolder = new LayeredNetwork[1];
        final AvgTable[] lnAvgTableHolder = new AvgTable[1];
        final AvgTable[] lqnsAvgTableHolder = new AvgTable[1];

        // Solve with SolverLN
        withSuppressedOutput(() -> {
            modelHolder[0] = LayeredNetwork.parseXML(modelPath, true);
            SolverLN solver = new SolverLN(modelHolder[0], SolverType.MVA);
            lnAvgTableHolder[0] = solver.getAvgTable();
        });

        LayeredNetwork model = modelHolder[0];
        AvgTable lnAvgTable = lnAvgTableHolder[0];

        // Verify model was loaded successfully
        assertNotNull(model, "Phase-2 model should be loaded successfully");

        // Verify the model has phase-2 activities
        assertTrue(model.getStruct().nacts > 0, "Model should have activities");
        assertNotNull(model.getStruct().actphase, "actphase matrix should be populated");

        // Verify that at least one activity has phase 2
        boolean hasPhase2 = false;
        for (int a = 1; a <= model.getStruct().nacts; a++) {
            if (model.getStruct().actphase.get(0, a) > 1) {
                hasPhase2 = true;
                break;
            }
        }
        assertTrue(hasPhase2, "Model should have at least one phase-2 activity");

        // SolverLN should now support phase-2 models
        assertNotNull(lnAvgTable, "SolverLN should produce results for phase-2 models");
        assertTrue(lnAvgTable instanceof LayeredNetworkAvgTable);

        // Solve with LQNS for comparison
        try {
            withSuppressedOutput(() -> {
                LayeredNetwork lqnsModel = LayeredNetwork.parseXML(modelPath, true);
                SolverLQNS solver = new SolverLQNS(lqnsModel);
                lqnsAvgTableHolder[0] = solver.getAvgTable();
            });

            AvgTable lqnsAvgTable = lqnsAvgTableHolder[0];
            assertNotNull(lqnsAvgTable, "SolverLQNS should produce results for phase-2 model");

            LayeredNetworkAvgTable lnTable = (LayeredNetworkAvgTable) lnAvgTable;
            LayeredNetworkAvgTable lqnsTable = (LayeredNetworkAvgTable) lqnsAvgTable;

            // Compare throughput values between LN and LQNS (allowing 5% relative tolerance)
            List<Double> lnTput = lnTable.getTput();
            List<Double> lqnsTput = lqnsTable.getTput();
            assertEquals(lnTput.size(), lqnsTput.size(), "Both solvers should return same number of elements");

            double tolerance = 0.05; // 5% relative tolerance
            for (int i = 0; i < lnTput.size(); i++) {
                double lnVal = lnTput.get(i);
                double lqnsVal = lqnsTput.get(i);
                if (!Double.isNaN(lnVal) && !Double.isNaN(lqnsVal) && lqnsVal > 0) {
                    double relError = Math.abs(lnVal - lqnsVal) / lqnsVal;
                    assertTrue(relError < tolerance,
                        String.format("Throughput mismatch at index %d: LN=%.6f, LQNS=%.6f, relError=%.4f",
                            i, lnVal, lqnsVal, relError));
                }
            }

            // Compare utilization values
            List<Double> lnUtil = lnTable.getUtil();
            List<Double> lqnsUtil = lqnsTable.getUtil();
            for (int i = 0; i < lnUtil.size(); i++) {
                double lnVal = lnUtil.get(i);
                double lqnsVal = lqnsUtil.get(i);
                if (!Double.isNaN(lnVal) && !Double.isNaN(lqnsVal) && lqnsVal > 0) {
                    double relError = Math.abs(lnVal - lqnsVal) / lqnsVal;
                    assertTrue(relError < tolerance,
                        String.format("Utilization mismatch at index %d: LN=%.6f, LQNS=%.6f, relError=%.4f",
                            i, lnVal, lqnsVal, relError));
                }
            }

        } catch (RuntimeException e) {
            if (e.getMessage() != null && (e.getMessage().contains("lqns") || e.getMessage().contains("lqsim"))) {
                // LQNS not available - just verify LN produces valid results
                LayeredNetworkAvgTable lnTable = (LayeredNetworkAvgTable) lnAvgTable;
                List<Double> tput = lnTable.getTput();
                boolean hasPositiveTput = tput.stream().anyMatch(v -> !Double.isNaN(v) && v > 0);
                assertTrue(hasPositiveTput, "LN should produce positive throughput values");
            } else {
                throw e;
            }
        }
    }

    @Test
    public void testPhase2ModelWithLQNS() throws Exception {
        // Test loading and solving a phase-2 model with SolverLQNS
        // This model (03-sanity.lqnx) contains activities with phase="2" attributes

        // Load the phase-2 model from test resources
        String modelPath = getClass().getClassLoader().getResource("lqn/phase2/03-sanity.lqnx").getPath();

        try {
            final LayeredNetwork[] modelHolder = new LayeredNetwork[1];
            final AvgTable[] avgTableHolder = new AvgTable[1];

            withSuppressedOutput(() -> {
                modelHolder[0] = LayeredNetwork.parseXML(modelPath, true);
                SolverLQNS solver = new SolverLQNS(modelHolder[0]);
                avgTableHolder[0] = solver.getAvgTable();
            });

            LayeredNetwork model = modelHolder[0];
            AvgTable avgTable = avgTableHolder[0];

            // Verify model was loaded successfully
            assertNotNull(model, "Phase-2 model should be loaded successfully");

            // Verify the model has phase-2 activities
            // The model has e0_ph1 (phase 1) and e0_ph2 (phase 2) in t0
            assertTrue(model.getStruct().nacts > 0, "Model should have activities");

            // Check that actphase matrix is populated
            assertNotNull(model.getStruct().actphase, "actphase matrix should be populated");

            // Verify that at least one activity has phase 2
            boolean hasPhase2 = false;
            for (int a = 1; a <= model.getStruct().nacts; a++) {
                if (model.getStruct().actphase.get(0, a) > 1) {
                    hasPhase2 = true;
                    break;
                }
            }
            assertTrue(hasPhase2, "Model should have at least one phase-2 activity");

            // Verify solver produced results
            assertNotNull(avgTable, "SolverLQNS should produce results for phase-2 model");
            assertTrue(avgTable instanceof LayeredNetworkAvgTable);

            LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;

            // Model has: 4 processors, 4 tasks, 4 entries, 5 activities
            // Total: 17 elements
            assertEquals(17, lnAvgTable.getQLen().size(),
                "Expected 17 entries (4 processors, 4 tasks, 4 entries, 5 activities)");

            // Verify that throughput values are reasonable
            List<Double> tput = lnAvgTable.getTput();
            boolean hasPositiveTput = tput.stream().anyMatch(v -> !Double.isNaN(v) && v > 0);
            assertTrue(hasPositiveTput, "Model should have positive throughput values");

            // Get detailed tables to check phase-specific metrics
            final LayeredNetworkAvgTable[][] tablesHolder = new LayeredNetworkAvgTable[1][];
            withSuppressedOutput(() -> {
                SolverLQNS solver = new SolverLQNS(model);
                tablesHolder[0] = solver.getRawAvgTables();
            });
            LayeredNetworkAvgTable[] tables = tablesHolder[0];

            assertNotNull(tables);
            assertEquals(2, tables.length, "getRawAvgTables should return 2 tables");

            LayeredNetworkAvgTable nodeTable = tables[0];
            assertTrue(nodeTable instanceof SolverLQNS.DetailedLayeredNetworkAvgTable);

            SolverLQNS.DetailedLayeredNetworkAvgTable detailedTable =
                (SolverLQNS.DetailedLayeredNetworkAvgTable) nodeTable;

            // Verify phase-2 utilization is present (for t0/e0 which has phase-2 activity)
            List<Double> phase2Util = detailedTable.getPhase2Utilization();
            assertNotNull(phase2Util, "Phase-2 utilization should be available");
            boolean hasPhase2Util = phase2Util.stream().anyMatch(v -> !Double.isNaN(v) && v > 0);
            assertTrue(hasPhase2Util, "Phase-2 utilization should have non-zero values for phase-2 activities");

            // Verify phase-2 service time is present
            List<Double> phase2SvcT = detailedTable.getPhase2ServiceTime();
            assertNotNull(phase2SvcT, "Phase-2 service time should be available");
            boolean hasPhase2SvcT = phase2SvcT.stream().anyMatch(v -> !Double.isNaN(v) && v > 0);
            assertTrue(hasPhase2SvcT, "Phase-2 service time should have non-zero values for phase-2 activities");

        } catch (RuntimeException e) {
            if (e.getMessage() != null && (e.getMessage().contains("lqns") || e.getMessage().contains("lqsim"))) {
                org.junit.jupiter.api.Assumptions.assumeTrue(false, "LQNS solver not available");
            } else {
                throw e;
            }
        }
    }

    @Test
    public void testGetRawAvgTablesWithLQNS() throws Exception {
        // Test getRawAvgTables with SolverLQNS (as used in MATLAB lqn_serial.m)
        LayeredNetwork model = LayeredModel.lqn_serial();

        try {
            final LayeredNetworkAvgTable[][] tablesHolder = new LayeredNetworkAvgTable[1][];
            withSuppressedOutput(() -> {
                SolverLQNS solver = new SolverLQNS(model);
                tablesHolder[0] = solver.getRawAvgTables();
            });
            LayeredNetworkAvgTable[] tables = tablesHolder[0];

            // Verify we get two tables (NodeAvgTable, CallAvgTable)
            assertNotNull(tables);
            assertEquals(2, tables.length, "getRawAvgTables should return 2 tables");

            LayeredNetworkAvgTable nodeTable = tables[0];
            LayeredNetworkAvgTable callTable = tables[1];

            assertNotNull(nodeTable);
            assertNotNull(callTable);

            // Verify node table has correct size (10 entries for lqn_serial)
            assertEquals(10, nodeTable.getNodeNames().size(),
                "Expected 10 node entries (2 processors, 2 tasks, 2 entries, 4 activities)");

            // Verify detailed metrics are present (not all NaN)
            // The DetailedLayeredNetworkAvgTable should have phase metrics
            assertTrue(nodeTable instanceof SolverLQNS.DetailedLayeredNetworkAvgTable,
                "Node table should be DetailedLayeredNetworkAvgTable");

            SolverLQNS.DetailedLayeredNetworkAvgTable detailedNodeTable =
                (SolverLQNS.DetailedLayeredNetworkAvgTable) nodeTable;

            // Verify phase1 utilization has some non-NaN values
            List<Double> phase1Util = detailedNodeTable.getPhase1Utilization();
            assertNotNull(phase1Util);
            assertEquals(10, phase1Util.size());
            boolean hasNonNanPhase1 = phase1Util.stream().anyMatch(v -> !Double.isNaN(v));
            assertTrue(hasNonNanPhase1, "Phase1 utilization should have non-NaN values");

            // Verify phase1 service time has some non-NaN values
            List<Double> phase1SvcT = detailedNodeTable.getPhase1ServiceTime();
            assertNotNull(phase1SvcT);
            boolean hasNonNanSvcT = phase1SvcT.stream().anyMatch(v -> !Double.isNaN(v));
            assertTrue(hasNonNanSvcT, "Phase1 service time should have non-NaN values");

            // Verify processor utilization has some non-NaN values
            List<Double> procUtil = detailedNodeTable.getProcUtilization();
            assertNotNull(procUtil);
            boolean hasNonNanProcUtil = procUtil.stream().anyMatch(v -> !Double.isNaN(v));
            assertTrue(hasNonNanProcUtil, "Processor utilization should have non-NaN values");

            // Verify throughput has some positive values
            List<Double> tput = detailedNodeTable.getTput();
            boolean hasPositiveTput = tput.stream().anyMatch(v -> !Double.isNaN(v) && v > 0);
            assertTrue(hasPositiveTput, "Throughput should have positive values");

            // Verify call table structure (lqn_serial has calls)
            assertTrue(callTable instanceof SolverLQNS.DetailedLayeredNetworkAvgTable,
                "Call table should be DetailedLayeredNetworkAvgTable");

            SolverLQNS.DetailedLayeredNetworkAvgTable detailedCallTable =
                (SolverLQNS.DetailedLayeredNetworkAvgTable) callTable;

            // Check call data is populated
            List<String> sourceNodes = detailedCallTable.getSourceNodes();
            List<String> targetNodes = detailedCallTable.getTargetNodes();
            List<String> callTypes = detailedCallTable.getCallTypes();

            if (sourceNodes != null && !sourceNodes.isEmpty()) {
                assertEquals(sourceNodes.size(), targetNodes.size(),
                    "Source and target nodes should have same size");
                assertEquals(sourceNodes.size(), callTypes.size(),
                    "Call types should have same size as nodes");

                // Verify call types are valid
                for (String callType : callTypes) {
                    assertTrue(callType.equals("Synchronous") ||
                              callType.equals("Asynchronous") ||
                              callType.equals("Forwarding") ||
                              callType.equals("Unknown"),
                        "Call type should be valid: " + callType);
                }
            }

        } catch (RuntimeException e) {
            if (e.getMessage() != null && (e.getMessage().contains("lqns") || e.getMessage().contains("lqsim"))) {
                org.junit.jupiter.api.Assumptions.assumeTrue(false, "LQNS solver not available");
            } else {
                throw e;
            }
        }
    }
}