package jline.examples.basic;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.examples.java.basic.ForkJoinModel;
import jline.lang.Network;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.mam.MAMOptions;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.assertTableMetrics;
import static jline.TestTools.withSuppressedOutput;
import static jline.TestTools.COARSE_TOL;
import java.util.List;

/**
 * Unit tests for fork-join network examples, comparing Java results with MATLAB expected values.
 *
 * ANNOTATION: The following fork-join examples are not present in allExamplesBaseline.txt:
 * - fj_asymm (both JMT and MVA solvers)
 * - fj_basic_closed (both JMT and MVA solvers)
 *
 * The examples present in allExamplesBaseline.txt are:
 * - fj_basic_open (lines 853-869)
 * - fj_twoclasses_forked (lines 873-893)
 * - fj_basic_nesting (lines 897-917)
 * - fj_nojoin (lines 921-938)
 * Current expected values for these are based on Java implementation baseline.
 *
 * Examples updated with values from allExamplesBaseline.txt:
 * - fj_basic_open
 * - fj_basic_nesting
 * - fj_twoclasses_forked
 * - fj_nojoin
 *
 * METHOD NAMING (amva vs egflin): the 10 *MVA fork-join tests assert
 * "default/egflin". SolverMVA forces an AMVA sub-solve on fork-join models, and
 * the JAR labels the resulting method by the specific approximation applied
 * (egflin = Extended Gelenbe-Fourneau linearizer), whereas MATLAB reports the
 * generic family name and its banner prints "default/amva". The two names denote
 * the same algorithm, so the numbers agree to ~1e-15 (see below); only the label
 * differs. These tests previously asserted "default/exact" with pre-fix JAR
 * output mislabelled "MATLAB ground truth"; both the label and the values were
 * stale. Expected values are now taken from a MATLAB SolverMVA run on the
 * matching matlab/examples/basic/forkJoin/<model>.m script and verified to agree
 * with the JAR to <= 5.6e-07 relative error (max, on fj_complex_serial RespT;
 * all others ~1e-15), well inside MID_TOL = 1e-4.
 */
public class ForkJoinExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        // Set verbose level to SILENT to suppress warnings during tests
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }


    @Test
    public void testFjAsymmJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_asymm();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];

        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from MATLAB ground truth (JMT solver)
        // Order: Delay1, Queue1, Queue2, Queue3, Join (5 entries total)
        double[] expectedQLen = {1.81080548167309, 6.15130102418248, 0.824031979378134, 5.12694767609361, 4.29300474254624};
        double[] expectedUtil = {1.81080548167309, 0.952809692713114, 0.470643060703271, 0.930831259442108, 0};
        double[] expectedRespT = {1.96139330322507, 6.58921973251044, 0.855700789438968, 5.40486666385993, 2.23960378697351};
        double[] expectedResidT = {1.96139330322507, 6.58921973251044, 0.855700789438968, 5.40486666385993, 2.23960378697351};
        double[] expectedArvR = {0.949137874469015, 0.949824984786469, 0.949824984786469, 0.951921370080833, 1.89539084624116};
        double[] expectedTput = {0.949824984786469, 0.948637654696326, 0.951921370080833, 0.944117963860923, 0.946849658018034};

        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }


    @Test
    public void testFjAsymmMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_asymm();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (SolverMVA, method default/egflin)
        // Order: Delay1, Queue1, Queue2, Queue3, Join (5 entries total)
        double[] expectedQLen = {1.71385178080942, 5.22364332535547, 0.736020657148363, 5.22364332535547, 5.61591953878298};
        double[] expectedUtil = {1.71385165162851, 0.856926317551341, 0.42846315877567, 0.856926317551341, 0.0};
        double[] expectedRespT = {2.00000015074923, 6.09579052290281, 0.858907752129185, 6.09579052290281, 3.27678087588115};
        double[] expectedResidT = {2.00000015074923, 6.09579052290281, 0.858907752129185, 6.09579052290281, 3.27678087588115};
        double[] expectedArvR = {0.856925825814257, 0.856925825814257, 0.856925825814257, 0.856926317551341, 1.71385263510268};
        double[] expectedTput = {0.856925825814257, 0.856926317551341, 0.856926317551341, 0.856926317551341, 0.856925825814257};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }




    @Test
    public void testFjBasicClosedJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_basic_closed();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from MATLAB ground truth (JMT solver)
        // Order: Delay, Queue1, Queue2, Join (4 entries total)
        double[] expectedQLen = {0.866250468016722, 2.64746759356294, 2.58634009650029, 2.96334227612133};
        double[] expectedUtil = {0.866250468016722, 0.888116827362, 0.893364626013386, 0};
        double[] expectedRespT = {0.996113064593301, 3.09452771040838, 2.81503789297669, 1.68503517449312};
        double[] expectedResidT = {0.996113064593301, 3.09452771040838, 2.81503789297669, 1.68503517449312};
        double[] expectedArvR = {0.89441678542412, 0.897680290928009, 0.897680290928009, 1.72396472683137};
        double[] expectedTput = {0.897680290928009, 0.897605937561591, 0.897570691458667, 0.887985227840386};

        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    //
    // HT results look better than MATLAB version, Disable?
    public void testFjBasicClosedMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_basic_closed();
        SolverOptions options = SolverMVA.defaultOptions();
        options.config.fork_join = "ht";
        options.method = "lin";
        SolverMVA solver = new SolverMVA(model, options);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("lin", solver.result.method, "MVA solver should use lin method");
        
        // Expected values from MATLAB ground truth (MVA solver)
        // Order: Delay, Queue1, Queue2, Join (4 entries total)
        double[] expectedQLen = {0.933095134562007, 2.71120038218631, 2.71120038218631, 2.71135226091848};
        double[] expectedUtil = {0.933095134562007, 0.933123478023213, 0.933123478023213, 0};
        double[] expectedRespT = {1, 2.90551084185545, 2.90551084185545, 1.45283680283256};
        double[] expectedResidT = {1, 2.90551084185545, 2.90551084185545, 1.45283680283256};
        double[] expectedArvR = {0.933095134562007, 0.933095134562007, 0.933095134562007, 1.86624695604643};
        double[] expectedTput = {0.933095134562007, 0.933123478023213, 0.933123478023213, 0.933095134562007};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }





    @Test
    public void testFjBasicNestingJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_basic_nesting();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from MATLAB output (JMT solver)
        // Order (9): Delay(c1,c2), Join1(c1,c2), Join1_1(c2), Queue1(c1,c2), Queue2(c1,c2)
        double[] expectedQLen = {2.03678408566221, 0.986416327912589, 2.20264565363805, 1.04793990425473, 0.541063167623458, 1.34469141104993, 0.829687618795429, 2.32992918028417, 1.11491852392952};
        double[] expectedUtil = {2.03678408566221, 0.986416327912589, 0, 0, 0, 0.495092471762105, 0.268125576910165, 0.663103825463966, 0.258694532021825};
        double[] expectedRespT = {4.02829661047477, 3.91380668072405, 2.1777244924427, 1.00107063051576, 1.05553569668762, 2.64853822683008, 1.62968498287804, 4.6827043023585, 2.19470719986523};
        double[] expectedResidT = {4.02829661047477, 3.91380668072405, 2.1777244924427, 1.00107063051576, 1.05553569668762, 2.64853822683008, 1.62968498287804, 4.6827043023585, 2.19470719986523};
        double[] expectedArvR = {0.50636916366692, 0.254265878516369, 0.997305649879641, 1.04448493362023, 0.511834697073124, 0.504868860727642, 0.511780527922063, 0.504868860727642, 0.511780527922063};
        double[] expectedTput = {0.504868860727642, 0.254101099226443, 0.507981831481078, 0.511834697073124, 0.254100490574222, 0.505101824537261, 0.512040780439188, 0.504465290617657, 0.513517879521406};
        

        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    // Retains COARSE_TOL from when the baseline was stale pre-fix JAR output (0.23% off).
    // Against the current MATLAB ground truth the JAR agrees to 3.4e-15, so this could be
    // tightened to the default MID_TOL; left as-is pending approval to change a tolerance.
    public void testFjBasicNestingMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_basic_nesting();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB output (MVA solver)
        // Order (9): Delay(c1,c2), Join1(c1,c2), Join1_1(c2), Queue1(c1,c2), Queue2(c1,c2)
        //
        // 2026-08-14: every row moved when the MMT stopped SCALING the order statistic by
        // tasksPerLink and started taking it over the sibling multiset. Fork1_1 has ONE
        // outgoing link and w = 2, so its join now synchronises on E[max of 2] = 1.5*R
        // rather than on 2*R. The total throughput at the Delay rises from 0.661037 to
        // 0.676917, i.e. TOWARDS SolverJMT 0.7556 and SolverLDES 0.8041. See
        // _kb/05-solvers-overview.md.
        double[] expectedQLen = {1.7618157938612555, 0.9414507504257009, 2.403736104206229, 0.9988789418839202, 0.7176814856530526, 1.2833912068955193, 0.6931871725711469, 2.899609674220357, 1.1786605286479332};
        double[] expectedUtil = {1.7618156897237687, 0.9414507265821999, 0.0, 0.0, 0.0, 0.44045387633869637, 0.23536279063272852, 0.5872718351182618, 0.23536279063272852};
        double[] expectedRespT = {4.000000236432193, 4.00000010130536, 2.728703541205555, 0.8487993950325206, 1.5246286447946946, 2.913792512314339, 1.4725929504567052, 6.5832311394863074, 2.5039228279867998};
        double[] expectedResidT = {4.000000236432193, 4.00000010130536, 2.728703541205555, 0.8487993950325206, 1.5246286447946946, 2.913792512314339, 1.4725929504567052, 6.5832311394863074, 2.5039228279867998};
        double[] expectedArvR = {0.44045392243094217, 0.23536268164554988, 0.8809077526773927, 0.9414511625309141, 0.47072543737343664, 0.44045392243094217, 0.23536268164554996, 0.44045392243094217, 0.23536268164554996};
        double[] expectedTput = {0.44045392243094217, 0.23536268164554996, 0.44045392243094217, 0.47072543737343664, 0.23536268164554988, 0.44045387633869637, 0.47072558126545705, 0.44045387633869637, 0.47072558126545705};

        // Check all metrics against expected values (relaxed 1% tolerance for numerical precision)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
    }

    @Test
    public void testFjBasicOpenJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_basic_open();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from allExamplesBaseline.txt (JMT solver)
        // Order: Source, Queue1, Queue2, Join (4 entries total)
        double[] expectedQLen = {0, 0.0519230302101923, 0.0266905244079356, 0.0444621888114224};
        double[] expectedUtil = {0, 0.0495240984212132, 0.0263285055754063, 0};
        double[] expectedRespT = {0, 1.04775456858692, 0.509727536551953, 0.432372216704833};
        double[] expectedResidT = {0, 1.04775456858692, 0.509727536551953, 0.432372216704833};
        double[] expectedArvR = {0, 0.0504940257909067, 0.0504940257909067, 0.103950920242054};
        double[] expectedTput = {0.0504940257909067, 0.0508116878670037, 0.0508958235892095, 0.0508568672155207};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjBasicOpenMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_basic_open();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB output (MVA solver)
        // Order: Source, Queue1, Queue2, Join (4 entries total)
        double[] expectedQLen = {0, 0.052631082553342, 0.0256408668165232, 0.0437892733328934};
        double[] expectedUtil = {0, 0.0499999880790734, 0.0249999940395367, 0};
        double[] expectedRespT = {0, 1.05262190203141, 0.51281745859565, 0.437892837730701};
        double[] expectedResidT = {0, 1.05262190203141, 0.51281745859565, 0.437892837730701};
        double[] expectedArvR = {0, 0.05, 0.05, 0.0999999761581469};
        double[] expectedTput = {0.05, 0.0499999880790734, 0.0499999880790734, 0.05};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjBasicOpenMAMDecSourceFJ() {
        Network model = ForkJoinModel.fj_basic_open();
        SolverMAM solver = new SolverMAM(model, new MAMOptions().method("dec.source.mmap"));
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];

        assertNotNull(avgTable);
        assertEquals("dec.source.mmap", solver.result.method, "MAM solver should use dec.source.mmap");

        double[] expectedQLen = {0, 0.0526315789473684, 0.0256410256410256, 0.043797123015873};
        double[] expectedUtil = {0, 0.05, 0.025, 0};
        double[] expectedRespT = {0, 1.05263157894737, 0.512820512820513, 0.43797123015873};
        double[] expectedResidT = {0, 1.05263157894737, 0.512820512820513, 0.43797123015873};
        double[] expectedArvR = {0, 0.05, 0.05, 0.1};
        double[] expectedTput = {0.05, 0.05, 0.05, 0.05};

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
    }

    @Test
    public void testFjTwoclassesForkedJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_twoclasses_forked();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from allExamplesBaseline.txt (JMT solver)
        // Order: Source(class1,class2), Queue1(class1,class2), Queue2(class1,class2), Join(class1,class2)
        double[] expectedQLen = {0, 0, 1.42515823093542, 0, 8.79577525314906, 3.4512327314647, 14.4680460376634, 6.75404224340901};
        double[] expectedUtil = {0, 0, 0.492305382371796, 0, 0.655665790777323, 0.25836437326815, 0, 0};
        double[] expectedRespT = {0, 0, 2.94021918610871, 0, 17.4055164446147, 6.73949559625344, 14.6730599499623, 6.30501162418713};
        double[] expectedResidT = {0, 0, 2.94021918610871, 0, 17.4055164446147, 6.73949559625344, 14.6730599499623, 6.30501162418713};
        double[] expectedArvR = {0, 0, 0.502526834989652, 0.505311572388569, 0.502526834989652, 0.505311572388569, 0.986802343233743, 1.0280464825212};
        double[] expectedTput = {0.251277570391642, 0.251840954693361, 0.490039441787383, 0.505311572388569, 0.508626139874509, 0.504633607485049, 0.247629026738526, 0.251541038492206};

        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testFjTwoclassesForkedMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_twoclasses_forked();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB output (MVA solver)
        // Order: Source(class1,class2), Queue1(class1,class2), Queue2(class1,class2), Join(class1,class2)
        //
        // 2026-08-14: the Join rows moved (class1 23.4443944123188 -> 15.082319339601389,
        // class2 8.99998091225116 -> 5.999987271500773) when the MMT stopped SCALING the
        // order statistic by tasksPerLink and started taking it over the sibling multiset,
        // each branch replicated tasksPerLink times. This fork has B = 2 links and w = 2,
        // so the join now synchronises on E[X_(4)] over [2,16,2,16] = 24.0824 rather than
        // on 2*E[X_(2)] over [2,16] = 32.4444; the delay is that minus the mean branch
        // time, 9. See _kb/05-solvers-overview.md: on this OPEN model the old number was
        // closer to SolverJMT/SolverLDES, but only because it compensated the transform's
        // Poisson treatment of what are really BATCHES of w tasks at a branch, which
        // under-states Queue2's queue on its own (11.0 here against LDES 15.47, and 11.0
        // is exactly rho/(1-rho) at that station's 0.917 utilisation).
        double[] expectedQLen = {0, 0, 0.999999983525951, 0, 7.999984025450736, 2.999994009544026, 15.082319339601389, 5.999987271500773};
        double[] expectedUtil = {0, 0, 0.4999999944120648, 0, 0.6666666592160864, 0.2499999972060324, 0, 0};
        double[] expectedRespT = {0, 0, 1.9999999894036427, 1.9999999894036427e-08, 15.999968229715044, 5.999988086143142, 15.082319508159438, 5.999987338555854};
        double[] expectedResidT = {0, 0, 1.9999999894036427, 0, 15.999968229715044, 5.999988086143142, 15.082319508159438, 5.999987338555854};
        double[] expectedArvR = {0, 0, 0.25, 0.25, 0.25, 0.25, 0.9999999888241295, 0.9999999888241295};
        double[] expectedTput = {0.25, 0.25, 0.4999999944120648, 0.4999999944120648, 0.4999999944120648, 0.4999999944120648, 0.25, 0.25};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjNojoinJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_nojoin();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from allExamplesBaseline.txt (JMT solver)
        // Order: Source, Queue1, Queue2, Queue3 (4 entries total)
        double[] expectedQLen = {0, 0.978243762828807, 0.327217087002569, 0.203795885841022};
        double[] expectedUtil = {0, 0.499522680470569, 0.2445235847298, 0.168070695581533};
        double[] expectedRespT = {0, 1.92285725036037, 0.667302743083952, 0.411455155419831};
        double[] expectedResidT = {0, 1.92285725036037, 0.667302743083952, 0.411455155419831};
        double[] expectedArvR = {0, 0.506199639620098, 0.506199639620098, 0.506199639620098};
        double[] expectedTput = {0.506199639620098, 0.501034118121665, 0.503983920778436, 0.506212543708862};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testFjNojoinMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_nojoin();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        //avgTable.print();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB output (MVA solver)
        // Order: Source, Queue1, Queue2, Queue3 (4 entries total)
        double[] expectedQLen = {0, 0.999998334258592, 0.333144466323753, 0.199129124223733};
        double[] expectedUtil = {0, 0.499999920527142, 0.249999960263571, 0.166666473509074};
        double[] expectedRespT = {0, 1.99999698640814, 0.666289038551295, 0.398258311748919};
        double[] expectedResidT = {0, 1.99999698640814, 0.666289038551295, 0.398258311748919};
        double[] expectedArvR = {0, 0.5, 0.5, 0.5};
        double[] expectedTput = {0.5, 0.499999920527142, 0.499999920527142, 0.499999920527142};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjDelaysJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_delays();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from MATLAB ground truth (JMT solver)
        // Order: Delay1, Delay2, Queue1, Queue2, Join (5 entries total)
        double[] expectedQLen = {1.97277663989911, 0.482194017145704, 6.04286454966976, 5.14634145263595, 3.97485244797253};
        double[] expectedUtil = {1.97277663989911, 0.482194017145704, 0.957044097028173, 0.945146738165965, 0};
        double[] expectedRespT = {2.03502847930761, 0.494330080933262, 6.08099768249331, 5.71924985902066, 2.03214940404023};
        double[] expectedResidT = {2.03502847930761, 0.494330080933262, 6.08099768249331, 5.71924985902066, 2.03214940404023};
        double[] expectedArvR = {0.967341427527411, 0.95414852606745, 0.954124124787102, 0.954124124787102, 1.9436324068198};
        double[] expectedTput = {0.95414852606745, 0.954124124787102, 0.964943098987972, 0.966163775525125, 0.967439971339213};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testFjDelaysMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_delays();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (SolverMVA, method default/egflin)
        // Order: Delay1, Delay2, Queue1, Queue2, Join (5 entries total)
        double[] expectedQLen = {1.71838707653906, 0.429596769134765, 5.30909244476314, 5.30909244476314, 5.30913707249604};
        double[] expectedUtil = {1.718386917399, 0.42959672934975, 0.859193874148908, 0.859193874148908, 0.0};
        double[] expectedRespT = {2.00000018522029, 0.500000046305072, 6.17915537400935, 6.17915537400935, 3.0896036577048};
        double[] expectedResidT = {2.00000018522029, 0.500000046305072, 6.17915537400935, 6.17915537400935, 3.0896036577048};
        double[] expectedArvR = {0.8591934586995, 0.8591934586995, 0.8591934586995, 0.8591934586995, 1.71838774829782};
        double[] expectedTput = {0.8591934586995, 0.8591934586995, 0.859193874148908, 0.859193874148908, 0.8591934586995};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjComplexSerialJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_complex_serial();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from MATLAB ground truth (JMT solver)
        // Order: Delay1, Queue1, Queue2, Queue3, Queue4, Queue5, Join (7 entries total)
        double[] expectedQLen = {1.51092356776505, 2.57533998809744, 0.591398625247183, 2.7793315598252, 0.32517425934687, 5.33883626711227, 5.33588343103474};
        double[] expectedUtil = {1.51092356776505, 0.762181810742536, 0.387378879253077, 0.763536260060229, 0.240915295180491, 0.954134610942161, 0};
        double[] expectedRespT = {1.9610778220606, 3.32929445729872, 0.794338897812184, 3.44717014256587, 0.418432269474096, 7.26332042070373, 3.46977199882321};
        double[] expectedResidT = {1.9610778220606, 3.32929445729872, 0.794338897812184, 3.44717014256587, 0.418432269474096, 7.26332042070373, 3.46977199882321};
        double[] expectedArvR = {0.765790622324946, 0.758226075205898, 0.758226075205898, 0.765777437113877, 0.759287624419417, 0.75933035365896, 1.51250860083258};
        double[] expectedTput = {0.758226075205898, 0.759287624419417, 0.765777437113877, 0.75992449185232, 0.75933035365896, 0.765790622324946, 0.758645940348197};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testFjComplexSerialMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_complex_serial();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (SolverMVA, method default/egflin)
        // Order: Delay1, Queue1, Queue2, Queue3, Queue4, Queue5, Join (7 entries total)
        double[] expectedQLen = {1.38670514152894, 2.14201184217135, 0.52374150008887, 2.14201184217135, 0.298357235581938, 5.60878291149622, 6.70973292646954};
        double[] expectedUtil = {1.38670510387262, 0.69334821931348, 0.34667410965674, 0.69334821931348, 0.231116073104493, 0.86668527414185, 0.0};
        double[] expectedRespT = {2.0000000543105, 3.08937382761619, 0.755380176222929, 3.08937382761619, 0.430313697029982, 8.08941705662671, 4.8386458200709};
        double[] expectedResidT = {2.0000000543105, 3.08937382761619, 0.755380176222929, 3.08937382761619, 0.430313697029982, 8.08941705662671, 4.8386458200709};
        double[] expectedArvR = {0.693352551936311, 0.693352551936311, 0.693352551936311, 0.69334821931348, 0.69334821931348, 0.69334821931348, 1.38669643862696};
        double[] expectedTput = {0.693352551936311, 0.69334821931348, 0.69334821931348, 0.69334821931348, 0.69334821931348, 0.69334821931348, 0.693352551936311};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjCsMultiVisitsMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_cs_multi_visits();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        //avgTable.print();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (MVA solver)
        // Order: Source(class1), Queue1(class1,class2), Queue2(class1,class2), Join(class1,class2)
        double[] expectedQLen = {0, 0.124736478395296, 0.124736478395296, 0.124736478395296, 0.124736478395296, 0.124736466947047, 0.124736466947047};
        double[] expectedUtil = {0, 0.0999999761581445, 0.0999999761581446, 0.0999999761581445, 0.0999999761581445, 0, 0};
        double[] expectedRespT = {0, 1.24736508134794, 1.24736508134794, 1.24736508134794, 1.24736508134794, 0.623682483432712, 0.623682483432712};
        double[] expectedResidT = {0, 1.24736508134794, 1.24736508134794, 1.24736508134794, 1.24736508134794, 0.623682483432712, 0.623682483432712};
        double[] expectedArvR = {0, 0.1, 0.1, 0.1, 0.1, 0.199999952316289, 0.199999952316289};
        double[] expectedTput = {0.1, 0.0999999761581445, 0.0999999761581446, 0.0999999761581445, 0.0999999761581445, 0.1, 0.1};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjCsMultiVisitsJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_cs_multi_visits();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        //avgTable.print();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from MATLAB ground truth (JMT solver)
        // Order: Source(class1), Queue1(class1,class2), Queue2(class1,class2), Join(class1,class2)
        double[] expectedQLen = {0, 0.128187709145439, 0.126895293186275, 0.126422095474972, 0.124299459664713, 0.134933191582454, 0.128880770894937};
        double[] expectedUtil = {0, 0.0970098413322245, 0.101098979917457, 0.101807219240264, 0.100702168523447, 0, 0};
        double[] expectedRespT = {0, 1.23082223883547, 1.24796055070012, 1.24402562013115, 1.23179853053884, 0.662883511843912, 0.640097183268929};
        double[] expectedResidT = {0, 1.23082223883547, 1.24796055070012, 1.24402562013115, 1.23179853053884, 0.662883511843912, 0.640097183268929};
        double[] expectedArvR = {0, 0.0999354173308032, 0.0999313110783683, 0.0999354173308032, 0.0999313110783683, 0.198813116156121, 0.198815854237995};
        double[] expectedTput = {0.0999354173308032, 0.099936811348074, 0.0999276613946777, 0.0999313110783683, 0.0999809178450228, 0.0999313110783683, 0.0999235968027199};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjCsPostforkMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_cs_postfork();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (SolverMVA, method default/egflin)
        // Order: Delay(class1), Join1(class1), Queue1(class1), Queue2(class1)
        double[] expectedQLen = {1.63644835900474, 0.246376289494928, 0.246388435090161, 0.246388435090161};
        double[] expectedUtil = {1.63644848775781, 0.0, 0.204522551255609, 0.204522551255609};
        double[] expectedRespT = {3.9999996852866, 0.301160297461538, 0.602350287480594, 0.602350287480594};
        double[] expectedResidT = {3.9999996852866, 0.301160297461538, 0.602350287480594, 0.602350287480594};
        double[] expectedArvR = {0.409112121939452, 0.818090205022437, 0.409112121939452, 0.409112121939452};
        double[] expectedTput = {0.409112121939452, 0.409112121939452, 0.409045102511219, 0.409045102511219};

        // Check all metrics against expected values (relaxed tolerance for numerical precision)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
    }

    @Test
    public void testFjCsPostforkJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_cs_postfork();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        //avgTable.print();

        // Check if results are computed
        assertNotNull(avgTable);

        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");

        // Expected values from MATLAB ground truth (JMT solver, seed=23000)
        // Order: Delay(class1,class2), Join1(class1,class2), Queue1(class1), Queue2(class1)
        // Class2 is dropped at Queue1/Queue2 due to post-fork class switch (class2 -> class1)
        double[] expectedQLen  = {0.839190241238998, 0.828942889175662, 0.213740249793893, 0, 0.212837842704092, 0.221079157084476};
        double[] expectedUtil  = {0.839190241238998, 0.828942889175662, 0, 0, 0.195306612719814, 0.199257505512527};
        double[] expectedRespT = {4.00076305880139, 3.88747366676825, 0.264904093025746, 0, 0.541401830774152, 0.54182474693752};
        double[] expectedResidT= {4.00076305880139, 0, 0.264904093025746, 0, 0.541401830774152, 0.54182474693752};
        double[] expectedArvR  = {0.206970443211937, 0.209478338486398, 0.824102008919696, 0, 0.414712442993219, 0.414712442993219};
        double[] expectedTput  = {0.209772151417075, 0.211918467751353, 0.207741310657977, 0.211915508719549, 0.414633802343402, 0.414695626983187};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjCsPreforkJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_cs_prefork();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from MATLAB ground truth (JMT solver)
        // Order: Delay1(class1), Delay2(class2), Queue1(class1), Queue2(class1), Join(class1)
        double[] expectedQLen = {1.90828256808915, 0.482123253308903, 12.7723873195026, 8.78386579695215, 13.6318001512482};
        double[] expectedUtil = {1.90828256808915, 0.482123253308903, 0.975939063267673, 0.942426828475797, 0};
        double[] expectedRespT = {2.02263207091286, 0.495351963945263, 12.4623458186284, 8.99882025866619, 7.08367429839283};
        double[] expectedResidT = {2.02263207091286, 0.495351963945263, 12.4623458186284, 8.99882025866619, 7.08367429839283};
        double[] expectedArvR = {0.97769851327727, 0.977387837106569, 0.977373629133931, 0.977373629133931, 1.92000407949589};
        double[] expectedTput = {0.977387837106569, 0.977373629133931, 0.978227185254143, 0.973432893602492, 0.97769851327727};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testFjCsPreforkMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_cs_prefork();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (SolverMVA, method default/egflin)
        // Order: Delay1(class1), Delay2(class2), Queue1(class1), Queue2(class1), Join(class1)
        double[] expectedQLen = {1.86584829890864, 0.466462074727161, 11.8695369405169, 11.8695369405169, 11.8701194207029};
        double[] expectedUtil = {1.86584822391925, 0.466462055979812, 0.932928561685477, 0.932928561685477, 0.0};
        double[] expectedRespT = {2.00000008038102, 0.500000020095255, 12.7228787154643, 12.7228787154643, 6.36175153607568};
        double[] expectedResidT = {2.00000008038102, 0.500000020095255, 12.7228787154643, 12.7228787154643, 6.36175153607568};
        double[] expectedArvR = {0.932924111959625, 0.932924111959625, 0.932924111959625, 0.932924111959625, 1.86585712337095};
        double[] expectedTput = {0.932924111959625, 0.932924111959625, 0.932928561685477, 0.932928561685477, 0.932924111959625};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjDeepNestingJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_deep_nesting();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from MATLAB ground truth (JMT solver)
        // Order: Delay1, Queue1, Queue2, Join, Queue3, Queue4, Join2 (7 entries total)
        double[] expectedQLen = {0.503009053607737, 0.247901049205163, 0.243412008475434, 0.309736361652484, 0.124683520044267, 0.127825430618943, 0.126024496115059};
        double[] expectedUtil = {0.503009053607737, 0.247901049205163, 0.243412008475434, 0, 0.124683520044267, 0.127825430618943, 0};
        double[] expectedRespT = {1.9944934497095, 1.01923289959026, 0.957810899502042, 0.626763745777582, 0.505112389067662, 0.504306520112183, 0.246947131178907};
        double[] expectedResidT = {1.9944934497095, 1.01923289959026, 0.957810899502042, 0.626763745777582, 0.505112389067662, 0.504306520112183, 0.246947131178907};
        double[] expectedArvR = {0.251697831929587, 0.249066498366798, 0.249066498366798, 0.505877820516287, 0.249063153894484, 0.249063153894484, 0.507680138946161};
        double[] expectedTput = {0.249066498366798, 0.249063153894484, 0.251709937811002, 0.249064137179797, 0.24844277504574, 0.249066577964212, 0.248448066935248};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    // Retains COARSE_TOL from when the baseline was stale pre-fix JAR output (0.12% off).
    // Against the current MATLAB ground truth the JAR agrees to 1.0e-15, so this could be
    // tightened to the default MID_TOL; left as-is pending approval to change a tolerance.
    public void testFjDeepNestingMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_deep_nesting();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (SolverMVA, method default/egflin)
        // Order: Delay1, Queue1, Queue2, Join, Queue3, Queue4, Join2 (7 entries total)
        double[] expectedQLen = {0.451926889149808, 0.273065938973322, 0.273065938973322, 0.392271337778741, 0.12634342279544, 0.12634342279544, 0.126338199914298};
        double[] expectedUtil = {0.451926547411844, 0.225964020751795, 0.225964020751795, 0.0, 0.112947182559932, 0.112947182559932, 0.0};
        double[] expectedRespT = {2.00000151236065, 1.20844875243774, 1.20844875243774, 0.867995126997723, 0.5593031182004, 0.5593031182004, 0.223698055285634};
        double[] expectedResidT = {2.00000151236065, 1.20844875243774, 1.20844875243774, 0.867995126997723, 0.5593031182004, 0.5593031182004, 0.223698055285634};
        double[] expectedArvR = {0.225963273705922, 0.225963273705922, 0.225963273705922, 0.45192804150359, 0.225964020751795, 0.225964020751795, 0.451788730239729};
        double[] expectedTput = {0.225963273705922, 0.225964020751795, 0.225964020751795, 0.225963273705922, 0.225894365119864, 0.225894365119864, 0.225964020751795};

        // Check all metrics against expected values (relaxed 1% tolerance for numerical precision)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
    }

    @Test
    public void testFjRouteOverlapJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_route_overlap();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Previous MAPE: 0.0008%, Max APE: 0.0317%
        // Expected values from MATLAB ground truth (JMT solver)
        // Mixed class model - Delay1(class1,class2), Queue1(class1,class2), Queue2(class1,class2), Join(class1,class2)
        double[] expectedQLen = {1.04336327249769, 2.24820736548409, 8.95663672750231, 7.75179263451591, 0.528170890708182, 0.4712231435402, 8.4186269515058, 7.23385678636148};
        double[] expectedUtil = {1.04336327249769, 2.24820736548409, 0.538334723283713, 0.458080172394163, 0.272610008241875, 0.226076314018791, 0, 0};
        double[] expectedRespT = {1.95675836848845, 4.92371582868094, 16.5670419560199, 16.7431049346619, 0.997967325760237, 1.0142734608016, 7.83145174892689, 7.72449632592672};
        double[] expectedResidT = {1.95675836848845, 4.92371582868094, 16.5670419560199, 16.7431049346619, 0.997967325760237, 1.0142734608016, 7.83145174892689, 7.72449632592672};
        double[] expectedArvR = {0.535611465239907, 0.46316176919595, 0.533611228036742, 0.46313440519477, 0.533611228036742, 0.46313440519477, 1.09160149698717, 0.917575442340489};
        double[] expectedTput = {0.533611228036742, 0.46313440519477, 0.535611465239907, 0.462601009210036, 0.533104076679514, 0.462901553323331, 0.535611465239907, 0.462601009210036};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testFjRouteOverlapMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_route_overlap();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (SolverMVA, method default/egflin)
        // Mixed class model - Delay1(class1,class2), Queue1(class1,class2), Queue2(class1,class2), Join(class1,class2)
        double[] expectedQLen = {1.02412733913278, 2.21501190185207, 9.00964434084338, 7.80752604994228, 0.484322573742232, 0.419358159481455, 8.57473967791028, 7.43092391666856};
        double[] expectedUtil = {1.02412729090567, 2.21501181161158, 0.512063659285942, 0.443002372665442, 0.256031829642971, 0.221501186332721, 0, 0};
        double[] expectedRespT = {2.00000009418188, 5.00000020370204, 17.5947739650321, 17.6241179092704, 0.945824928130235, 0.946627344134243, 8.37272819737638, 8.38700239003058};
        double[] expectedResidT = {2.00000009418188, 5.00000020370204, 17.5947739650321, 17.6241179092704, 0.945824928130235, 0.946627344134243, 8.37272819737638, 8.38700239003058};
        double[] expectedArvR = {0.512063645452833, 0.443002362322317, 0.512063645452832, 0.443002362322317, 0.512063645452832, 0.443002362322317, 1.02412731857188, 0.886004745330884};
        double[] expectedTput = {0.512063645452832, 0.443002362322317, 0.512063659285942, 0.443002372665442, 0.512063659285942, 0.443002372665442, 0.512063645452833, 0.443002362322317};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjSerialfjsClosedJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_serialfjs_closed();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from MATLAB ground truth (JMT solver)
        // Order: Delay1, Queue1, Queue2, Join, Queue3, Queue4, Join2 (7 entries total)
        double[] expectedQLen = {1.54783690531209, 2.69841819702249, 2.80349706122689, 3.18562234948136, 2.55975111144918, 2.68042013831799, 2.8827839833709};
        double[] expectedUtil = {1.54783690531209, 0.817439819940281, 0.803434078561481, 0, 0.80584195265541, 0.800791571059644, 0};
        double[] expectedRespT = {2.01840310683681, 3.43410471738878, 3.42270833270948, 1.93180181909228, 3.19438458010616, 3.39938878060967, 1.75375855116248};
        double[] expectedResidT = {2.01840310683681, 3.43410471738878, 3.42270833270948, 1.93180181909228, 3.19438458010616, 3.39938878060967, 1.75375855116248};
        double[] expectedArvR = {0.808194254806168, 0.813656520259378, 0.813656520259378, 1.58795842942525, 0.808876084116358, 0.808876084116358, 1.62085151889386};
        double[] expectedTput = {0.813656520259378, 0.808090489094702, 0.812672922993953, 0.808876084116358, 0.808154103094597, 0.809517418807689, 0.808001078018123};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testFjSerialfjsClosedMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_serialfjs_closed();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (SolverMVA, method default/egflin)
        // Order: Delay1, Queue1, Queue2, Join, Queue3, Queue4, Join2 (7 entries total)
        double[] expectedQLen = {1.51005382458859, 2.86537319828304, 2.86537319828304, 2.86536145405452, 2.86537319828304, 2.86537319828304, 2.86536145405452};
        double[] expectedUtil = {1.51005346867591, 0.755026085879607, 0.755026085879607, 0.0, 0.755026085879607, 0.755026085879607, 0.0};
        double[] expectedRespT = {2.00000047139083, 3.79506516645034, 3.79506516645034, 1.8975248058591, 3.79506516645034, 3.79506516645034, 1.8975248058591};
        double[] expectedResidT = {2.00000047139083, 3.79506516645034, 3.79506516645034, 1.8975248058591, 3.79506516645034, 3.79506516645034, 1.8975248058591};
        double[] expectedArvR = {0.755026734337954, 0.755026734337954, 0.755026734337954, 1.51005217175921, 0.755026734337954, 0.755026734337954, 1.51005217175921};
        double[] expectedTput = {0.755026734337954, 0.755026085879607, 0.755026085879607, 0.755026734337954, 0.755026085879607, 0.755026085879607, 0.755026734337954};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjSerialfjsOpenJMT() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_serialfjs_open();
        SolverJMT solver = new SolverJMT(model, "seed", 23000);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default", solver.result.method, "JMT solver should use default method");
        
        // Expected values from MATLAB ground truth (JMT solver)
        // Order: Source, Queue1, Queue2, Join1, Queue3, Queue4, Join2 (7 entries total)
        double[] expectedQLen = {0, 0.6862786234576, 0.695657281295777, 0.591550104471919, 0.689764915244344, 0.673096961784191, 0.613046849660212};
        double[] expectedUtil = {0, 0.401014929148586, 0.40012117796086, 0, 0.412923654913366, 0.413634954374046, 0};
        double[] expectedRespT = {0, 1.66746174839026, 1.76163085622374, 0.7631227484094, 1.68830734687715, 1.64893391695126, 0.725202815868963};
        double[] expectedResidT = {0, 1.66746174839026, 1.76163085622374, 0.7631227484094, 1.68830734687715, 1.64893391695126, 0.725202815868963};
        double[] expectedArvR = {0, 0.404735260928858, 0.404735260928858, 0.832915424364862, 0.408671723723258, 0.408671723723258, 0.833425388647928};
        double[] expectedTput = {0.404735260928858, 0.409198006343229, 0.404709876671848, 0.408671723723258, 0.408723427945088, 0.412360781875233, 0.412250699271869};
        
        // REBASED 2026-08-14: the fork export now writes an OutPathEntry for EVERY
        // outgoing link, not only the last, which leaves the model identical and the
        // RNG consumption different. At 1e4 samples these rows are pinned sample
        // paths -- five seeds span 4.9%% to 8.3%% on these models and 37%% on
        // fj_twoclasses_forked -- and both writers agree once the run is converged.
        // See _kb/08-build-and-test.md.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testFjSerialfjsOpenMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_serialfjs_open();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (MVA solver)
        // Order: Source, Queue1, Queue2, Join1, Queue3, Queue4, Join2 (7 entries total)
        double[] expectedQLen = {0, 0.66666476639387, 0.66666476639387, 0.666664713404591, 0.66666476639387, 0.66666476639387, 0.666664713404591};
        double[] expectedUtil = {0, 0.399999952316286, 0.399999952316286, 0, 0.399999952316286, 0.399999952316286, 0};
        double[] expectedRespT = {0, 1.66666211466628, 1.66666211466628, 0.833330991096531, 1.66666211466628, 1.66666211466628, 0.833330991096531};
        double[] expectedResidT = {0, 1.66666211466628, 1.66666211466628, 0.833330991096531, 1.66666211466628, 1.66666211466628, 0.833330991096531};
        double[] expectedArvR = {0, 0.4, 0.4, 0.799999904632571, 0.4, 0.4, 0.799999904632571};
        double[] expectedTput = {0.4, 0.399999952316286, 0.399999952316286, 0.4, 0.399999952316286, 0.399999952316286, 0.4};

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testFjThreebranchesMVA() {
        // Create the model and run solver with suppressed output
        Network model = ForkJoinModel.fj_threebranches();
        SolverMVA solver = new SolverMVA(model);
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            avgTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify solver method
        assertEquals("default/egflin", solver.result.method, "MVA solver should use default/egflin method");
        
        // Expected values from MATLAB ground truth (SolverMVA, method default/egflin)
        // Order: Delay1(class1,class2), Queue1(class1,class2), Queue2(class1,class2), Queue3(class1,class2), Join(class1,class2)
        double[] expectedQLen = {1.45910542772117, 0.797552416491446, 1.64128470250385, 0.773824908278897, 4.80364849286599, 1.55729408654501, 3.52897187158018, 7.65154531954084, 7.23150933868093, 8.55498524790506};
        double[] expectedUtil = {1.45910545308043, 0.797552390222842, 0.486368491096253, 0.227872116766866, 0.6632297605858, 0.212680642315742, 0.291821094657752, 0.638041926947226, 0.0, 0.0};
        double[] expectedRespT = {1.99999996524, 1.25000004117066, 2.24971358486439, 1.21281200434796, 6.58437457854062, 2.44073942600612, 4.83717172772443, 11.9922296582458, 4.95612515412023, 6.70409332568255};
        double[] expectedResidT = {1.99999996524, 1.25000004117066, 2.24971358486439, 1.21281200434796, 6.58437457854062, 2.44073942600612, 4.83717172772443, 11.9922296582458, 4.95612515412023, 6.70409332568255};
        double[] expectedArvR = {0.729552726540213, 0.638041912178273, 0.729552726540213, 0.638041912178273, 0.729552726540213, 0.638041912178273, 0.72955273664438, 0.638041926947226, 1.45910547328876, 1.27608385389445};
        double[] expectedTput = {0.729552726540213, 0.638041912178273, 0.72955273664438, 0.638041926947226, 0.72955273664438, 0.638041926947226, 0.72955273664438, 0.638041926947226, 0.729552726540213, 0.638041912178273};

        // REBASED 2026-08-14: the recorded vectors predate a fix that moved this
        // model in BOTH engines. MATLAB SolverMVA (the ground truth) and this JAR
        // now agree on it -- QLen to 1e-14 on the fork-join model and to ~1e-7 on
        // the priority one -- so the goldens were the stale party, not the solvers.
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    // ===== Additional fork-join models =====
    // ANNOTATION: The following fork-join models are not present in allExamplesBaseline.txt:
    // - test_forkJoinCS_1
    // Current expected values for these are based on Java implementation baseline.
}
