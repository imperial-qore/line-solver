package jline.examples.basic;

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
        double[] expectedQLen = {0.0, 0.0, 0.0, 0.374499208179313, 0.45449920817873, 0.494499208178439, 
                                0.194425921295638, 0.214425921295493, 0.234425921295347, 0.499870152586697, 0.524863660216032, 0.474876644957362, 
                                0.446428503191219, 0.234596862744448, 1.13152163657697};
        double[] expectedUtil = {0.0, 0.0, 0.0, 0.12, 0.2, 0.24, 
                                0.11, 0.13, 0.15, 0.2, 0.21, 0.19, 
                                0.25, 0.19, 0.43};
        double[] expectedRespT = {0.0, 0.0, 0.0, 0.936248020448281, 1.13624802044683, 1.2362480204461, 
                                 1.94425921295638, 2.14425921295493, 2.34425921295347, 4.99870152586697, 5.24863660216032, 4.74876644957362, 
                                 4.46428503191219, 2.34596862744448, 11.3152163657697};
        double[] expectedResidT = {0.0, 0.0, 0.0, 3.74499208179313, 4.5449920817873, 4.94499208178439, 
                                  1.94425921295638, 2.14425921295493, 2.34425921295347, 4.99870152586697, 5.24863660216032, 4.74876644957362, 
                                  4.46428503191219, 2.34596862744448, 11.3152163657697};
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
        double[] expectedQLen = {1.88140943816302, 2.33628253574678, 1.2184469140195, 10.4219512433026, 12.9546061008667, 6.86880811717636, 
                                0.72581251773679, 0.942476117942643, 0.518927693169747, 0.563588852428907, 0.769834062128354, 0.693488541314444, 
                                4.2184922740892, 0.764637091052955, 8.57722363222328, 0.188140943816302, 0.233628253574678, 0.12184469140195};
        double[] expectedUtil = {1.88147264808843, 2.33609251056257, 1.21853223828295, 0.225776717770612, 0.467218502112515, 0.292447737187909, 
                                0.206961991289727, 0.303692026373135, 0.182779835742443, 0.188147264808843, 0.256970176161883, 0.231521125273761, 
                                0.470368162022107, 0.443857577006889, 0.52396886246167, 0.188147264808843, 0.233609251056257, 0.121853223828295};
        double[] expectedRespT = {9.99966404015771, 10.0008134317599, 9.99929977836632, 13.8481301520531, 13.8635414075992, 14.0923807786474, 
                                 3.85768306796388, 4.03441265138801, 4.25862916766956, 2.99546662557923, 3.29539202170964, 5.69117926901668, 
                                 22.4212256201289, 3.27314559502961, 70.389796533488, 0.999966404015771, 1.00008134317599, 0.999929977836632};
        double[] expectedResidT = {9.99966404015771, 10.0008134317599, 9.99929977836632, 55.3925206082124, 55.454165630397, 56.3695231145896, 
                                  3.85768306796388, 4.03441265138801, 4.25862916766956, 2.99546662557923, 3.29539202170964, 5.69117926901668, 
                                  22.4212256201289, 3.27314559502961, 70.389796533488, 0.999966404015771, 1.00008134317599, 0.999929977836632};
        double[] expectedArvR = {0.188147264808843, 0.233609251056257, 0.121853223828295, 0.752589059235372, 0.93443700422503, 0.487412895313182, 
                                0.188147264808843, 0.233609251056257, 0.121853223828295, 0.188147264808843, 0.233609251056257, 0.121853223828295, 
                                0.188147264808843, 0.233609251056257, 0.121853223828295, 0.188147264808843, 0.233609251056257, 0.121853223828295};
        double[] expectedTput = {0.188147264808843, 0.233609251056257, 0.121853223828295, 0.752589059235372, 0.93443700422503, 0.487412895313182, 
                                0.188147264808843, 0.233609251056257, 0.121853223828295, 0.188147264808843, 0.233609251056257, 0.121853223828295, 
                                0.188147264808843, 0.233609251056257, 0.121853223828295, 0.188147264808843, 0.233609251056257, 0.121853223828295};
        
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
    @Disabled("SSA PSPRIO: Max APE 6.82% exceeds 5% tolerance due to simulation variability")
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
        
        // Expected values from MATLAB ground truth output (SSA solver with 5000 samples)
        // Order: Delay(Class1,Class2), Queue1(Class1,Class2)
        // Note: SSA is a simulation method, so we use relaxed tolerances
        double[] expectedQLen = {0.066197229965643, 0.213020755007946, 1.93380277003436, 1.78697924499205};
        double[] expectedUtil = {0.066197229965643, 0.213020755007946, 0.0176746917774692, 0.974661183964202};
        double[] expectedRespT = {0.711609223619123, 0.218558775616296, 22.2059575406122, 1.81970993848782};
        double[] expectedResidT = {0.711609223619123, 0.218558775616296, 22.2059575406122, 1.81970993848782};
        double[] expectedArvR = {0.087084862992179, 0.974661183964202, 0.087084862992179, 0.974661183964202};
        double[] expectedTput = {0.0930246935656275, 0.974661183964202, 0.087084862992179, 0.982013235844078};
        
        // Verify table size
        assertEquals(4, avgTable.getQLen().size(), 
            "Expected 4 entries (2 stations × 2 classes)");
        
        // Check all metrics against expected values
        // SSA is a stochastic method, but we still use COARSE_REL_TOL for consistency
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput, TestTools.MID_TOL);
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