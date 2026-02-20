package jline.examples.basic;
import jline.GlobalConstants;

import jline.TestTools;
import jline.examples.java.basic.PrioModel;
import jline.lang.Network;
import jline.solvers.NetworkAvgTable;
import jline.solvers.mva.SolverMVA;
import jline.solvers.ssa.SolverSSA;
import jline.solvers.ssa.SSAOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ctmc.CTMCOptions;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.VerboseLevel;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import static jline.TestTools.*;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for Priority examples, comparing solver outputs between MATLAB and Java implementations.
 * Expected values are obtained by running the corresponding MATLAB examples in the dev/ directory.
 */
public class PrioExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        // Suppress priority diagnostic output during tests
        jline.GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }
    

    @Test
    public void testPrioHolOpenMVA() {
        // Create the model
        Network model = PrioModel.prio_hol_open();
        
        // Create solver options matching MATLAB
        SolverOptions options = new SolverOptions();
        options.seed = 23000;
        options.verbose = VerboseLevel.SILENT;
        
        // Create and run the solver
        SolverMVA solver = new SolverMVA(model, options);
        NetworkAvgTable avgTable = solver.getAvgTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/egflin", solver.result.method, 
            "MVA solver should use default/egflin) method for priority models");
        
        // Expected values from MATLAB output (MVA solver ground truth)
        // Order: Source(Class1,Class2,Class3), WebServer(Class1,Class2,Class3), Storage1(Class1,Class2,Class3),
        //        Storage2(Class1,Class2,Class3), Storage3(Class1,Class2,Class3)
        double[] expectedQLen = {0.0, 0.0, 0.0, 0.37452833642695, 0.454528336426914, 0.494528336426896,
                                0.194426157598126, 0.214426157598117, 0.234426157598108, 0.499946814499407, 0.524944155224377, 0.474949473774437,
                                0.333377782277657, 1.46140608543643, 0.754517533480026};
        double[] expectedUtil = {0.0, 0.0, 0.0, 0.12, 0.2, 0.24,
                                0.11, 0.13, 0.15, 0.2, 0.21, 0.19,
                                0.25, 0.19, 0.43};
        double[] expectedRespT = {0.0, 0.0, 0.0, 0.936320841067375, 1.13632084106728, 1.23632084106724,
                                 1.94426157598126, 2.14426157598117, 2.34426157598108, 4.99946814499407, 5.24944155224377, 4.74949473774437,
                                 3.33377782277657, 14.6140608543643, 7.54517533480026};
        double[] expectedResidT = {0.0, 0.0, 0.0, 3.7452833642695, 4.54528336426914, 4.94528336426896,
                                  1.94426157598126, 2.14426157598117, 2.34426157598108, 4.99946814499407, 5.24944155224377, 4.74949473774437,
                                  3.33377782277657, 14.6140608543643, 7.54517533480026};
        double[] expectedArvR = {0.0, 0.0, 0.0, 0.4, 0.4, 0.4,
                                0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1};
        double[] expectedTput = {0.1, 0.1, 0.1, 0.4, 0.4, 0.4,
                                0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                0.1, 0.1, 0.1};
        
        // Verify table size
        assertEquals(15, avgTable.getQLen().size(), 
            "Expected 15 entries (5 stations × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    // @Disabled - MAPE was 0.0156%, Max APE was 0.0957%
    public void testPrioHolOpenSSA() {
        // Create the model
        Network model = PrioModel.prio_hol_open();
        
        // Create solver options matching MATLAB
        SSAOptions options = new SSAOptions();
        options.seed = 23000;
        options.samples = 10000;
        options.verbose = VerboseLevel.SILENT;
        
        // Create and run the solver
        SolverSSA solver = new SolverSSA(model, options);
        NetworkAvgTable avgTable = solver.getAvgTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/serial", solver.result.method, 
            "SSA solver should use default/serial method");
        
        // Expected values from MATLAB SSA output ground truth (10000 samples)
        // Order: Source(Class1,Class2,Class3), WebServer(Class1,Class2,Class3), Storage1(Class1,Class2,Class3),
        //        Storage2(Class1,Class2,Class3), Storage3(Class1,Class2,Class3)
        // Updated after HOL priority fix (min instead of max)
        double[] expectedQLen = {0.0, 0.0, 0.0, 0.333450898474251, 0.457690647368325, 0.446060053134282,
                                0.208480125467293, 0.222029480342276, 0.235992457066272, 0.648650958922419, 0.748636209389274, 0.570269359850089,
                                1.11194326596433, 4.7361897805337, 1.76709481787512};
        double[] expectedUtil = {0.0, 0.0, 0.0, 0.113570424512756, 0.198728459186634, 0.256533276077918,
                                0.106715460584875, 0.122581179397769, 0.14209071227756, 0.194028110154318, 0.198015751334857, 0.179981568884909,
                                0.242535137692897, 0.179157108350585, 0.407326708529004};
        double[] expectedRespT = {0.0, 0.0, 0.0, 0.859285023724256, 1.2134771514314, 1.17722345988812,
                                 2.20398203434299, 2.23951978742532, 2.16077335302854, 6.34263277136271, 7.0689174895663, 5.78195346566283,
                                 13.6088784259488, 51.2519147611441, 14.7615116648311};
        double[] expectedResidT = {0.0, 0.0, 0.0, 3.43714009489702, 4.8539086057256, 4.70889383955247,
                                  2.20398203434299, 2.23951978742532, 2.16077335302854, 6.34263277136271, 7.0689174895663, 5.78195346566283,
                                  13.6088784259488, 51.2519147611441, 14.7615116648311};
        double[] expectedArvR = {0.0, 0.0, 0.0, 0.378568081709188, 0.397456918373268, 0.427555460129864,
                                0.0970140550771589, 0.0942932149213604, 0.0947271415183731, 0.0970140550771589, 0.0942932149213604, 0.0947271415183731,
                                0.0970140550771589, 0.0942932149213604, 0.0947271415183731};
        double[] expectedTput = {0.1, 0.1, 0.1, 0.388056220308635, 0.377172859685442, 0.378908566073492,
                                0.0945924795296441, 0.0991415577522239, 0.109216663902072, 0.102268408451946, 0.105905354036776, 0.0986291853154365,
                                0.0817071937275981, 0.0924100065842686, 0.119709610912355};
        
        // Verify table size
        assertEquals(15, avgTable.getQLen().size(), 
            "Expected 15 entries (5 stations × 3 classes)");
        
        // Check all metrics against expected values with relaxed tolerance for SSA
        // SSA is stochastic, so use COARSE_TOL (1%) instead of MID_TOL (0.01%)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, TestTools.COARSE_TOL);
    }

    @Test
    //@Disabled
    public void testPrioHolClosedMVA() {
        // Create the model
        Network model = PrioModel.prio_hol_closed();
        
        // Create solver options matching MATLAB
        SolverOptions options = new SolverOptions();
        options.seed = 23000;
        options.verbose = VerboseLevel.SILENT;
        
        // Create and run the solver
        SolverMVA solver = new SolverMVA(model, options);
        NetworkAvgTable avgTable = solver.getAvgTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/egflin", solver.result.method, 
            "MVA solver should use default/egflin) method for priority models");
        
        // Expected values from MATLAB output (MVA solver ground truth)
        // Order: SlowDelay(Class1,Class2,Class3), FCFSQueue(Class1,Class2,Class3), SIROQueue(Class1,Class2,Class3),
        //        PSQueue(Class1,Class2,Class3), HOLQueue(Class1,Class2,Class3), FastDelay(Class1,Class2,Class3)
        double[] expectedQLen = {3.43636107729059, 0.000526312921992351, 2.18586482169469, 6.99976018194001, 0.00113874342457379, 4.68043242970567,
                                1.34409547636036, 0.000220710040261501, 0.943185705572623, 1.32120027060668, 0.000222590477315224, 1.59678617724074,
                                4.55489158652334, 17.9978563455951, 8.37514825093958, 0.343636107729059, 5.26312921992351e-05, 0.218586482169469};
        double[] expectedUtil = {3.43637163449741, 0.000526312415168568, 2.18586435205896, 0.412364596139689, 0.000105262483033714, 0.52460744449415,
                                0.378000879794715, 6.84206139719138e-05, 0.327879652808844, 0.343637163449741, 5.78943656685425e-05, 0.415314226891202,
                                0.859092908624353, 9.99993588820279e-05, 0.939921671385352, 0.343637163449741, 5.26312415168568e-05, 0.218586435205896};
        double[] expectedRespT = {9.99996927804108, 10.0000096297136, 10.0000021485127, 5.09240626921, 5.40906594521942, 5.35306825569593,
                                 3.91137984863777, 4.19351765036383, 4.31493246451545, 3.84475374358021, 4.22924618344663, 7.30505612453332,
                                 13.254944665464, 341961.462942704, 38.3150411097133, 0.999996927804108, 1.00000096297136, 1.00000021485127};
        double[] expectedResidT = {9.99996927804108, 10.0000096297136, 10.0000021485127, 20.36962507684, 21.6362637808777, 21.4122730227837,
                                  3.91137984863777, 4.19351765036383, 4.31493246451545, 3.84475374358021, 4.22924618344663, 7.30505612453332,
                                  13.254944665464, 341961.462942704, 38.3150411097133, 0.999996927804108, 1.00000096297136, 1.00000021485127};
        double[] expectedArvR = {0.343637163449741, 5.26312415168568e-05, 0.218586435205896, 1.37454865379896, 0.000210524966067427, 0.874345740823583,
                                0.343637163449741, 5.26312415168568e-05, 0.218586435205896, 0.343637163449741, 5.26312415168568e-05, 0.218586435205896,
                                0.343637163449741, 5.26312415168568e-05, 0.218586435205896, 0.343637163449741, 5.26312415168568e-05, 0.218586435205896};
        double[] expectedTput = {0.343637163449741, 5.26312415168568e-05, 0.218586435205896, 1.37454865379896, 0.000210524966067427, 0.874345740823583,
                                0.343637163449741, 5.26312415168568e-05, 0.218586435205896, 0.343637163449741, 5.26312415168568e-05, 0.218586435205896,
                                0.343637163449741, 5.26312415168568e-05, 0.218586435205896, 0.343637163449741, 5.26312415168568e-05, 0.218586435205896};
        
        // Verify table size
        assertEquals(18, avgTable.getQLen().size(), 
            "Expected 18 entries (6 stations × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, MID_TOL);
    }
    
    @Test
    // @Disabled - MAPE was 0.0096%, Max APE was 0.0669%
    public void testPrioHolClosedSSA() {
        // Create the model
        Network model = PrioModel.prio_hol_closed();
        
        // Create solver options matching MATLAB
        SSAOptions options = new SSAOptions();
        options.seed = 23000;
        options.samples = 10000;
        options.verbose = VerboseLevel.SILENT;
        
        // Create and run the solver
        SolverSSA solver = new SolverSSA(model, options);
        NetworkAvgTable avgTable = solver.getAvgTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/serial", solver.result.method, 
            "SSA solver should use default/serial method");
        
        // Expected values from MATLAB SSA output ground truth (10000 samples)
        // Order: SlowDelay(Class1,Class2,Class3), FCFSQueue(Class1,Class2,Class3), SIROQueue(Class1,Class2,Class3),
        //        PSQueue(Class1,Class2,Class3), HOLQueue(Class1,Class2,Class3), FastDelay(Class1,Class2,Class3)
        // Updated after HOL priority fix (min instead of max)
        double[] expectedQLen = {1.70708855766482, 0.0804039034390601, 1.62478314583399, 0.936613725956255, 0.244263960951659, 1.04147838119686,
                                0.341219081427859, 0.0181980378169357, 0.365678544981648, 0.381911029752903, 0.0197512627807857, 0.629775202536069,
                                14.4602748002778, 17.6329070136647, 14.182651906317, 0.172892804920345, 0.004475821346873, 0.155632819134412};
        double[] expectedUtil = {1.70708855766482, 0.0804039034390601, 1.62478314583399, 0.195892173467485, 0.0114901296786689, 0.371901352507472,
                                0.17399674412926, 0.00676427061795343, 0.233479636996806, 0.158178858299327, 0.00572361359980675, 0.295740873529287,
                                0.395447145748317, 0.00988624167239347, 0.669308292724176, 0.172892804920345, 0.004475821346873, 0.155632819134412};
        double[] expectedRespT = {10, 10, 10, 1.48030801338803, 11.7360454353477, 1.67275569712388,
                                 2.12150270946362, 2.34111625508587, 2.3592087718123, 2.31499942615152, 3.4004242423306, 3.90148557741801,
                                 92.4249365611837, 12982.8644936051, 100.630834744004, 1, 1, 1};
        double[] expectedResidT = {10, 10, 10, 5.92123205355213, 46.944181741391, 6.69102278849554,
                                  2.12150270946362, 2.34111625508587, 2.3592087718123, 2.31499942615152, 3.4004242423306, 3.90148557741801,
                                  92.4249365611837, 12982.8644936051, 100.630834744004, 1, 1, 1};
        double[] expectedArvR = {0.172892804920345, 0.004475821346873, 0.155632819134412, 0.652973911558282, 0.0229802593573377, 0.619835587512453,
                                0.158178858299327, 0.00520328509073341, 0.155653091331204, 0.158178858299327, 0.00520328509073341, 0.155653091331204,
                                0.158178858299327, 0.00520328509073341, 0.155653091331204, 0.158178858299327, 0.00520328509073341, 0.155653091331204};
        double[] expectedTput = {0.170708855766482, 0.00804039034390601, 0.162478314583399, 0.632715433197308, 0.0208131403629336, 0.622612365324815,
                                0.160838390592548, 0.00777323115731738, 0.155000502435713, 0.164972407957697, 0.0058084701711362, 0.161419333748467,
                                0.156454257241555, 0.00135816768497816, 0.140937436744875, 0.172892804920345, 0.004475821346873, 0.155632819134412};
        
        // Verify table size
        assertEquals(18, avgTable.getQLen().size(), 
            "Expected 18 entries (6 stations × 3 classes)");
        
        // Check all metrics against expected values with relaxed tolerance for SSA
        // SSA is stochastic, so use COARSE_TOL (1%) instead of MID_TOL (0.01%)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, TestTools.COARSE_TOL);
    }

    @Test
    public void testPrioPsprioSSA() {
        // Create the model
        Network model = PrioModel.prio_psprio();

        // Create solver options matching MATLAB
        SSAOptions options = new SSAOptions();
        options.seed = 23000;
        options.samples = 5000;
        options.verbose = VerboseLevel.SILENT;

        // Create and run the SSA solver (MVA doesn't support PSPRIO)
        SolverSSA solver = new SolverSSA(model, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        // Check if results are computed
        assertNotNull(avgTable);

        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/serial", solver.result.method,
            "SSA solver should use default/serial method");

        // Expected values from MATLAB ground truth output (SSA solver with 5000 samples, seed 23000)
        // Order: Delay(Class1,Class2), PSPRIOQueue(Class1,Class2)
        // MATLAB SSA results:
        //   Delay Class1: QLen=1.428, Util=1.428, RespT=0.65673, ArvR=2.1485, Tput=2.1744
        //   Delay Class2: QLen=0.11833, Util=0.11833, RespT=0.20959, ArvR=0.53552, Tput=0.56457
        //   Queue Class1: QLen=0.57202, Util=0.41313, RespT=0.26624, ArvR=2.1744, Tput=2.1485
        //   Queue Class2: QLen=1.8817, Util=0.56457, RespT=3.5137, ArvR=0.56457, Tput=0.53552
        double[] expectedQLen = {1.428, 0.11833, 0.57202, 1.8817};
        double[] expectedUtil = {1.428, 0.11833, 0.41313, 0.56457};
        double[] expectedRespT = {0.65673, 0.20959, 0.26624, 3.5137};
        double[] expectedResidT = {0.65673, 0.20959, 0.26624, 3.5137};
        double[] expectedArvR = {2.1485, 0.53552, 2.1744, 0.56457};
        double[] expectedTput = {2.1744, 0.56457, 2.1485, 0.53552};

        // Verify table size
        assertEquals(4, avgTable.getQLen().size(),
            "Expected 4 entries (2 stations × 2 classes)");

        // Check all metrics against expected values
        // SSA is stochastic, so use COARSE_TOL (1%) for comparison
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, TestTools.COARSE_TOL);
    }
    
    @Test
    public void testPrioIdenticalCTMC() {
        // Create the model
        Network model = PrioModel.prio_identical();

        // Create solver options matching MATLAB
        CTMCOptions options = new CTMCOptions();
        options.verbose = VerboseLevel.SILENT;

        // Create and run the CTMC solver (MVA doesn't support GPSPRIO)
        SolverCTMC solver = new SolverCTMC(model, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        // Check if results are computed
        assertNotNull(avgTable);

        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default", solver.result.method,
            "CTMC solver should use default method");

        // Expected values from Java CTMC solver output
        // Order: Delay(Class1,Class2,Class3,Class4), Queue1(Class1,Class2,Class3,Class4)
        double[] expectedQLen = {5.637201085468526, 0.9138354467020573, 3.0493048042157045, 0.003195443300095893, 0.3627989145314713, 3.0861645532979445, 0.9506951957842928, 0.9968045566999044};
        double[] expectedUtil = {5.637201085468526, 0.9138354467020573, 3.0493048042157045, 0.003195443300095893, 0.2818600542734262, 0.45691772335102865, 0.25410873368464204, 0.006390886600191786};
        double[] expectedRespT = {0.6666666666666669, 1.0, 1.0, 0.5, 0.042905324708827346, 3.377155662363173, 0.3117744065696374, 155.97281239060902};
        double[] expectedResidT = {0.6666666666666669, 1.0, 1.0, 0.5, 0.042905324708827346, 3.377155662363173, 0.3117744065696374, 155.97281239060902};
        double[] expectedArvR = {8.455801628202781, 0.9138354467020313, 3.0493048042156956, 0.006390886600182386, 8.455801628202787, 0.9138354467020573, 3.0493048042157045, 0.006390886600191786};
        double[] expectedTput = {8.455801628202787, 0.9138354467020573, 3.0493048042157045, 0.006390886600191786, 8.455801628202781, 0.9138354467020313, 3.0493048042156956, 0.006390886600182386};
        
        // Verify table size
        assertEquals(8, avgTable.getQLen().size(), 
            "Expected 8 entries (2 stations × 4 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

}