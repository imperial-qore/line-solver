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
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
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
    @Timeout(300)
    public void testLqnWorkflows() throws Exception {
        // Ground truth values from MATLAB SolverLN

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

        // MATLAB SolverLN ground truth values for lqn_workflows (updated after WN dedup fix)
        // Order: P1, P2, P3, T1, T2, T3, Entry, E2, E1, A1, A2, A3, B1, B2, B3, B4, B5, B6, C1, C2, C3, C4, C5
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 1.0, 0.23872, 0.072181, 1.0, 0.23872, 0.072181, 0.076867, 0.4612, 0.46193, 0.0079319, 0.015859, 0.023788, 0.031718, 0.039647, 0.11977, 0.0079319, 0.0047592, 0.0071388, 0.012691, 0.03966};
        double[] expectedUtil = {0.76867, 0.16654, 0.072181, 0.76867, 0.16654, 0.072181, Double.NaN, Double.NaN, Double.NaN, 0.076867, 0.4612, 0.2306, 0.0079319, 0.015859, 0.023788, 0.031718, 0.039647, 0.047592, 0.0079319, 0.0047592, 0.0071388, 0.012691, 0.03966};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 13.01, 3.0096, 0.91, 1, 2, 6.0096, 0.1, 0.2, 0.3, 0.4, 0.5, 1.51, 0.1, 0.2, 0.3, 0.4, 0.5};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, 10, 2.0996, 0.91, Double.NaN, Double.NaN, Double.NaN, 1, 6, 3, 0.1, 0.19994, 0.29991, 0.39988, 0.49984, 0.6, 0.1, 0.06, 0.09, 0.16, 0.5};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.076867, 0.079319, 0.079319, 0.076867, 0.079319, 0.079319, 0.076867, 0.2306, 0.076867, 0.079319, 0.079295, 0.079295, 0.079295, 0.079295, 0.079319, 0.079319, 0.023796, 0.023796, 0.031728, 0.079319};

        // Verify table size
        assertEquals(23, lnAvgTable.getQLen().size(),
            "Expected 23 entries (3 processors, 3 tasks, 3 entries, 14 activities)");

        // Check all metrics against expected values (using COARSE_TOL for MATLAB SolverLN ground truth)
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
    }

    @Test
    @Timeout(600)
    public void testLqnOfbiz() throws Exception {
        // Create the model - loads from ofbizExample.xml
        LayeredNetwork model = LayeredModel.lqn_ofbiz();

        // Check if model was loaded successfully
        org.junit.jupiter.api.Assumptions.assumeTrue(model != null,
                "ofbizExample.xml not available - skipping test");

        // Create and run the solver - use NC inner solver to match MATLAB:
        // solver{2} = LN(model, @(x) NC(x,'verbose',false))
        final AvgTable[] avgTableHolder = new AvgTable[1];
        withSuppressedOutput(() -> {
            SolverLN solver = new SolverLN(model, SolverType.NC);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0];

        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);

        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;

        // Ground truth values from JAR LN(NC) execution, verified against MATLAB LN(NC).
        // 72 entries: 9 processors, 9 tasks (incl USAGE_DELAY), 14 entries, 40 activities
        double[] expectedQLen = {
            // Processors [0-8]
            Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks [9-17]
            0.12014478873612254, 0.0, 0.0658167346853804, 0.04482643353065431, 0.01821000973320567, 0.036894826460213456, 0.018210009733205697, 0.06271765775677658, 0.007591267646677314,
            // Entries [18-31]
            0.0, 0.04285706749983818, 0.017388521618376968, 0.03526202746395451, 0.017388521618376968, 0.007248650535575906, 0.0, 0.0658167346853804, 0.04482643353065431, 0.01821000973320567, 0.036894826460213456, 0.018210009733205697, 0.06271765775677658, 0.007591267646677314,
            // Activities [32-71]
            0.0, 0.04287374723738245, 0.017565341809287194, 0.03538034585691672, 0.017565341809287194, 0.0073204160326338094, 0.0,
            0.0, 0.011018188613418328, 0.010906413389436738, 0.010983765339835111, 0.010983765339835111, 0.01090641338943673, 0.011018188613418326, 0.0,
            0.0, 0.04482643353065429, 0.0, 0.0, 0.018210009733205662, 0.0, 0.0, 0.036894826460213435, 0.0, 0.0, 0.018210009733205686, 0.0,
            0.0, 0.007875504433064236, 0.007795610513758687, 0.007850899604338112, 0.007850899604338112, 0.007875504433064229, 0.007798124221390411, 0.007795610513758683, 0.00787550443306423, 0.0,
            0.0, 0.00759126764667731, 0.0
        };

        double[] expectedUtil = {
            // Processors [0-8]
            0.10784908101591033, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Tasks [9-17]
            0.10784908101591033, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Entries [18-31]
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Activities [32-71]
            0.0, 0.0383073348717404, 0.015694486117561673, 0.031612043358719207, 0.015694486117561673, 0.00654073055032736, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0
        };

        double[] expectedRespT = {
            // Processors [0-8]
            Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks [9-17]
            Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Entries [18-31]
            0.0, 0.011187692290087245, 0.011079382585785247, 0.011154618214272017, 0.011079382585785233, 0.011082325559509176, 0.0, 0.06677408794494746, 0.011183131551230392, 0.01107208464331075, 0.011149084213947735, 0.011072084643310767, 0.089021167503389, 0.011075325341654244,
            // Activities [32-71]
            0.0, 0.011192046479069597, 0.011192046479069597, 0.011192046479069597, 0.011192046479069597, 0.011192046479069597, 0.0,
            0.0, 0.01117845634523462, 0.011065055267653195, 0.011143532359585888, 0.011143532359585887, 0.011065055267653192, 0.01117845634523462, 0.0,
            0.0, 0.011183131551230387, 0.0, 0.0, 0.011072084643310744, 0.0, 0.0, 0.01114908421394773, 0.0, 0.0, 0.011072084643310762, 0.0,
            0.0, 0.011178456345234618, 0.011065055267653195, 0.011143532359585887, 0.011143532359585887, 0.011178456345234618, 0.011068623213207004, 0.011065055267653195, 0.011178456345234618, 0.0,
            0.0, 0.011075325341654238, 0.0
        };

        double[] expectedResidT = {
            // Processors [0-8]
            Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks [9-17]
            0.01114008460752596, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Entries [18-31]
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Activities [32-71]
            0.0, 0.003973799970860531, 0.001612301324645579, 0.0032695714355461844, 0.001612301324645579, 6.721105518280845E-4, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0
        };

        double[] expectedArvR = new double[72]; // All NaN or 0
        Arrays.fill(expectedArvR, Double.NaN);
        expectedArvR[1] = 0.0;  // P2 (USAGE_DELAY processor)
        expectedArvR[10] = 0.0; // USAGE_DELAY_Task
        expectedArvR[24] = 0.0; // USAGE_DELAY0_Entry
        expectedArvR[38] = 0.0; // USAGE_DELAY0_Activity

        double[] expectedTput = {
            // Processors [0-8]
            Double.NaN, 0.0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks [9-17]
            10.784908101591595, 0.0, 0.9856628029070743, 4.008397229819085, 1.6446776121972053, 3.3092248432438303, 1.644677612197205, 0.7045252215366524, 0.6854216388683946,
            // Entries [18-31]
            1.7974846835985996, 3.8307334871742325, 1.5694486117562447, 3.1612043358720916, 1.5694486117562467, 0.6540730550327689, 0.0, 0.9856628029070743, 4.008397229819085, 1.6446776121972053, 3.3092248432438303, 1.644677612197205, 0.7045252215366524, 0.6854216388683946,
            // Activities [32-71]
            1.7974846835985983, 3.83073348717424, 1.5694486117562492, 3.1612043358720854, 1.5694486117562492, 0.6540730550327701, 0.0,
            0.9856628029070743, 0.9856628029070745, 0.9856628029070745, 0.9856628029070743, 0.9856628029070745, 0.9856628029070741, 0.9856628029070743, 0.9856628029070742,
            4.008397229819085, 4.008397229819085, 4.008397229819084, 1.6446776121972053, 1.6446776121972053, 1.644677612197205, 3.3092248432438303, 3.3092248432438303, 3.30922484324383, 1.644677612197205, 1.644677612197205, 1.6446776121972049,
            0.7045252215366539, 0.704525221536654, 0.704525221536654, 0.7045252215366534, 0.7045252215366534, 0.7045252215366534, 0.7045252215366536, 0.7045252215366535, 0.7045252215366536, 0.7045252215366536,
            0.6854216388683946, 0.6854216388683946, 0.6854216388683945
        };

        // Verify table size
        assertEquals(72, lnAvgTable.getQLen().size(),
            "Expected 72 entries for lqn_ofbiz model");

        // Check all metrics against expected values
        // Use COARSE_TOL since LN solver convergence differs slightly between MATLAB and JAR
        assertTableMetrics(lnAvgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
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

    @Test
    @Tag("slow") // ~103s
    @Timeout(value = 30, unit = java.util.concurrent.TimeUnit.MINUTES)
    public void testLqnBpmn() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_bpmn();
        SolverLN solver = new SolverLN(model, SolverType.MVA);
        AvgTable avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }

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