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
import jline.solvers.wrappers.lqns.SolverLQNS;
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
        // Source: matlab/examples/basic/layeredModel/lqn_serial.m run with SolverLN under the
        // default layering method, which the srvn alias resolves to srvn.ph (was srvn.cs)
        // Order: P1, P2, T1, T2, E1, E2, AS1, AS2, AS3, AS4
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.17334899671, 0.529599062451, 1.17334899671, 0.529599062451, 0.161430012703, 1.01191898313, 0.441332552042, 0.0882665104084};
        double[] expectedUtil = {0.141226416036, 0.529599062451, 0.141226416036, 0.529599062451, 0.141226416036, 0.529599062451, 0.141226416036, 0.0, 0.441332552042, 0.0882665104084};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 13.3, 6.0, 1.82889311769, 11.4643592781, 5.0, 1.0};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.82889311769, 6.0, Double.NaN, Double.NaN, 1.82889311769, 0.0, 5.0, 1.0};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0882665100228, 0.0882665104084, 0.0882665100228, 0.0882665104084, 0.0882665100228, 0.0882665100228, 0.0882665104084, 0.0882665104084};
        
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
        double[] expectedUtil = {0.24999, 0.74998, 0.24999, 0.74998, 0.24999, 0.74998, 0.24999, 0.74998};
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
        double[] expectedUtil = {0.25, 1.0, 0.25, 1.0, 0.25, 0.75, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
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
        LayeredNetwork model = LayeredModel.lqn_setup();
        
        // Create and run the solver
        final AvgTable[] avgTableHolder = new AvgTable[1];
        withSuppressedOutput(() -> {
            LNOptions lnoptions = new LNOptions();
            lnoptions.seed = 23000;
            // the routing encoding, which this golden was recorded under; the
            // 'srvn' alias now resolves to 'srvn.ph' on this model
            lnoptions.method = "srvn.cs";
            SolverLN solver = new SolverLN(model, SolverType.MVA, lnoptions);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0]; //avgTable.print();
        
        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);
        
        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;
        
        // Expected values from MATLAB ground truth for lqn_setup
        // Source: MATLAB SolverLN(model, @(m) SolverMVA(m)) with method='srvn.cs',
        // re-recorded 2026-08-11; the JAR reproduces every row to 1e-16.
        // The previous golden (Util 0.43577, F2 QLen 0.564231, F2 ResidT 1.29479)
        // charged F2's declared think time as a per-request delay. Commit
        // a9f6c2cf5 dropped that reading in MATLAB after lqns, lqsim and LDES all
        // contradicted it: a think time on a SERVED task is not a cycle
        // component, only a reference task's is. F2 is served, so its Exp(1/8)
        // leaves the client delay and T1 runs at X=0.5.
        // Order: P1, P2, T1, F2, E1, E2, A1, A2
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0, 0.5, 1.0, 0.5, 1.0, 0.16666676540009195};
        double[] expectedUtil = {0.5, 0.04166666077165185, 0.5, 0.04166666077165185, 0.5, 0.04166666077165185, 0.5, 0.04166666077165185};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 2.0, 1.0, 2.0, 0.3333335779603371};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.0, 0.3333335779603372, Double.NaN, Double.NaN, 1.0, 0.3333335779603372};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        
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
            // the routing encoding, which this golden was recorded under; the
            // 'srvn' alias now resolves to 'srvn.ph' on this model
            LNOptions lnoptions = new LNOptions();
            lnoptions.method = "srvn.cs";
            SolverLN solver = new SolverLN(model, SolverType.MVA, lnoptions);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0];

        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);

        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;

        // MATLAB SolverLN ground truth values for lqn_workflows (regenerated 2026-07-21;
        // MATLAB, JAR and Python-native agree to display precision on every row)
        // Order: P1, P2, P3, T1, T2, T3, Entry, E2, E1, A1, A2, A3, B1, B2, B3, B4, B5, B6, C1, C2, C3, C4, C5
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 1.0, 0.20949, 0.072181, 1.0, 0.20949, 0.072181, 0.079107, 0.47464, 0.44625, 0.0079319, 0.015871, 0.023807, 0.031742, 0.039678, 0.11977, 0.0079319, 0.0047592, 0.0071388, 0.012691, 0.03966};
        double[] expectedUtil = {0.79107, 0.16662, 0.072181, 0.79107, 0.16662, 0.072181, 0.79107, 0.16662, 0.072181, 0.079107, 0.47464, 0.23732, 0.0079319, 0.015871, 0.023807, 0.031742, 0.039678, 0.047592, 0.0079319, 0.0047592, 0.0071388, 0.012691, 0.03966};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 12.641, 2.641, 0.91, 1, 2, 5.641, 0.1, 0.2, 0.3, 0.4, 0.5, 1.51, 0.1, 0.2, 0.3, 0.4, 0.5};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, 10, 2.1, 0.91, Double.NaN, Double.NaN, Double.NaN, 1, 6, 3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.1, 0.06, 0.09, 0.16, 0.5};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.079107, 0.079319, 0.079319, 0.079107, 0.079319, 0.079319, 0.079107, 0.23732, 0.079107, 0.079319, 0.079356, 0.079356, 0.079356, 0.079356, 0.079319, 0.079319, 0.023796, 0.023796, 0.031728, 0.079319};

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
            // The default encoding, NOT the routing one: this golden reproduces
            // under 'srvn.ph' to 1e-9 on every row and under 'srvn.cs' to no
            // better than 17% on Util, so the blanket 'srvn.cs' pin of f81e42237
            // (which re-recorded the OTHER goldens it pinned, but not this one)
            // did not apply here.
            SolverLN solver = new SolverLN(model, SolverType.NC);
            avgTableHolder[0] = solver.getAvgTable();
        });
        AvgTable avgTable = avgTableHolder[0];

        // Check if results are computed
        assertNotNull(avgTable);
        assertTrue(avgTable instanceof LayeredNetworkAvgTable);

        LayeredNetworkAvgTable lnAvgTable = (LayeredNetworkAvgTable) avgTable;

        // Ground truth values from JAR LN(NC) execution (relax_factor=0.5).
        // 72 entries: 9 processors, 9 tasks (incl USAGE_DELAY), 14 entries, 40 activities
        double[] expectedQLen = {
            // Processors [0-8]
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks [9-17]
            0.13148104173001546, 0.0, 0.06733341453388977, 0.04649983203346833, 0.019240689267888116, 0.03848137846513493, 0.019240689267888217, 0.06414762947362763, 0.008018453619040551,
            // Entries [18-31]
            0.0, 0.0464998316044658, 0.01924068918517455, 0.03848137816603592, 0.019240689185174648, 0.008018453589164533, 0.0, 0.06733341453388977, 0.04649983203346833, 0.019240689267888116, 0.03848137846513493, 0.019240689267888217, 0.06414762947362763, 0.008018453619040551,
            // Activities [32-71]
            0.0, 0.04649983161116429, 0.019240689187505915, 0.03848137817127757, 0.01924068918750593, 0.008018453590127497, 0.0,
            0.0, 0.011222235741742827, 0.011222235754285551, 0.011222235746088575, 0.01122223574608858, 0.011222235754285529, 0.011222235741742824, 0.0,
            0.0, 0.046499832007273416, 0.0, 0.0, 0.019240689269877816, 0.0, 0.0, 0.03848137845225767, 0.0, 0.0, 0.019240689269877785, 0.0,
            0.0, 0.008018453674258736, 0.008018453683220699, 0.008018453677363832, 0.008018453677363832, 0.008018453674258727, 0.008018453684642783, 0.008018453683220681, 0.00801845367425873, 0.0,
            0.0, 0.008018453620633878, 0.0
        };

        double[] expectedUtil = {
            // Processors [0-8]
            0.11637229666872696, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Tasks [9-17]
            0.11637229666872696, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Entries [18-31]
            0.0, 0.0411564445137948, 0.0170297037540806, 0.0340594073278384, 0.0170297037540806, 0.00709703731893273, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Activities [32-71]
            0.0, 0.041156444513794764, 0.017029703754080555, 0.03405940732783835, 0.017029703754080562, 0.007097037318932728, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0
        };

        double[] expectedRespT = {
            // Processors [0-8]
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks [9-17]
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Entries [18-31]
            0.0, 0.011298311152431949, 0.01129831115269053, 0.01129831115252055, 0.011298311152690565, 0.011298311152702651, 0.0, 0.06778986691654053, 0.011298311152461465, 0.011298311152837913, 0.011298311152589289, 0.011298311152837972, 0.0903864892221563, 0.011298311152859165,
            // Activities [32-71]
            0.0, 0.011298311154059522, 0.011298311154059522, 0.011298311154059524, 0.011298311154059522, 0.011298311154059524, 0.0,
            0.0, 0.011298311138757018, 0.01129831115138477, 0.01129831114313223, 0.01129831114313223, 0.011298311151384752, 0.011298311138757018, 0.0,
            0.0, 0.01129831114609675, 0.0, 0.0, 0.011298311154006284, 0.0, 0.0, 0.011298311148808467, 0.0, 0.0, 0.011298311154006264, 0.0,
            0.0, 0.011298311138757016, 0.011298311151384768, 0.011298311143132228, 0.011298311143132228, 0.011298311138757015, 0.011298311153388547, 0.011298311151384754, 0.011298311138757015, 0.0,
            0.0, 0.011298311155104223, 0.0
        };

        double[] expectedResidT = {
            // Processors [0-8]
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks [9-17]
            0.01129831115255007, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            // Entries [18-31]
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Activities [32-71]
            0.0, 0.003995781894879102, 0.001653373675134042, 0.0033067473332314564, 0.0016533736751340425, 6.890345741714277E-4, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0
        };

        // Nothing reports an arrival rate on an LQN, the USAGE_DELAY* component
        // included: being unreachable makes its measures idle, not its missing
        // ones present.
        double[] expectedArvR = new double[72];
        Arrays.fill(expectedArvR, Double.NaN);

        double[] expectedTput = {
            // Processors [0-8]
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN,
            // Tasks [9-17]
            11.637229666872786, 0.0, 0.9932666576375977, 4.115644489339259, 1.7029703827067317, 3.405940759236035, 1.7029703827067317, 0.7097037403008594, 0.7097037345277388,
            // Entries [18-31]
            0.0, 4.115644451379511, 1.7029703754080676, 3.4059407327838613, 1.702970375408071, 0.7097037318932796, 0.0, 0.9932666576375977, 4.115644489339259, 1.7029703827067317, 3.405940759236035, 1.7029703827067317, 0.7097037403008594, 0.7097037345277388,
            // Activities [32-71]
            0.0, 4.115644451379509, 1.7029703754080687, 3.405940732783862, 1.70297037540807, 0.7097037318932784, 0.0,
            0.9932666576375977, 0.993266657637598, 0.993266657637598, 0.9932666576375977, 0.993266657637598, 0.9932666576375976, 0.9932666576375977, 0.9932666576375976,
            4.115644489339259, 4.115644489339259, 4.1156444893392585, 1.7029703827067317, 1.7029703827067317, 1.7029703827067313, 3.405940759236035, 3.405940759236035, 3.4059407592360342, 1.7029703827067317, 1.7029703827067317, 1.7029703827067313,
            0.709703740300861, 0.7097037403008611, 0.7097037403008611, 0.7097037403008604, 0.7097037403008604, 0.7097037403008604, 0.7097037403008607, 0.7097037403008605, 0.7097037403008607, 0.7097037403008607,
            0.7097037345277388, 0.7097037345277388, 0.7097037345277387
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
        //
        // Re-recorded 2026-08-10, when a served task's declared think time stopped
        // entering the thread cycle as a per-request delay (SolverLN.refThinkTime).
        // T3 is a served task with think time 4 and demand 0.02, and charging it per
        // request capped E3 at 25/(4+0.02) = 6.2189 completions per second, the rate
        // the rows here replace. AS2 makes five calls per E2 request at X(E2)=13.279,
        // so call-flow balance alone puts X(E3) at 66.4, and LDES (500k samples, seed
        // 23000) reads 66.744 with P2 utilization 0.44518 against the 0.44263 below.
        // MATLAB SolverLN reproduces these rows digit for digit.
        double[] expectedQLen = {Double.NaN, Double.NaN, 23.438871284695853, 8.7, 1.3278875842291844, 23.438871284695853, 8.7, 1.3278875842291844, 23.441553568169066, 8.7, 1.3278875842291264};
        double[] expectedUtil = {0.9959344198398389, 0.44262927744750075, 0.6639624468785424, 0.33197197296129644, 0.44262927744750075, 0.6639624468785424, 0.33197197296129644, 0.44262927744750075, 0.6639624468785424, 0.33197197296129644, 0.44262927744750075};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.765075072761118, 0.6554023155872206, 0.019999996263038014, 1.7652770633620782, 0.6555033822101997, 0.019999996263037143};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.090033703344034, 0.5551306450196893, 0.019999996263038014, Double.NaN, Double.NaN, Double.NaN, 1.090033703344034, 0.5551306450196893, 0.019999996263038014};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 13.279248937570843, 13.278878918451866, 66.4, 13.279248937570843, 13.278878918451866, 66.4, 13.279248937570848, 13.27887891845186, 66.4};
        
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
        {
            modelHolder[0] = LayeredNetwork.parseXML(modelPath, false);
            SolverLN solver = new SolverLN(modelHolder[0], SolverType.MVA);
            lnAvgTableHolder[0] = solver.getAvgTable();
        }

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
                LayeredNetwork lqnsModel = LayeredNetwork.parseXML(modelPath, false);
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

            // Cross-solver tolerance vs the LQNS rolia reference: after the
            // LQNS-parity rework of LN forwarding chains (cascaded, callmean-scaled
            // chain delays and per-call chain classes), SolverLN gives X(t0)=0.9556
            // vs LQNS rolia 0.9619 on 03-sanity, with all rows within ~2%.
            double tolerance = 0.05;
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
                modelHolder[0] = LayeredNetwork.parseXML(modelPath, false);
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

    @Test
    @Timeout(300)  // 5 minute timeout for SolverLN iteration
    public void testLqnSockshop() throws Exception {
        // Create the model
        LayeredNetwork model = LayeredModel.lqn_sockshop();

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

        // Verify table size: 7 processors + 7 tasks + 12 entries + 24 activities = 50
        assertEquals(50, lnAvgTable.getQLen().size(),
            "Expected 50 entries (7 processors, 7 tasks, 12 entries, 24 activities)");

        // Verify key metrics are computed (not all NaN/zero)
        List<Double> tput = lnAvgTable.getTput();
        boolean hasPositiveTput = tput.stream().anyMatch(v -> !Double.isNaN(v) && v > 0);
        assertTrue(hasPositiveTput, "Throughput should have positive values");

        List<Double> util = lnAvgTable.getUtil();
        boolean hasPositiveUtil = util.stream().anyMatch(v -> !Double.isNaN(v) && v > 0);
        assertTrue(hasPositiveUtil, "Utilization should have positive values");
    }
}