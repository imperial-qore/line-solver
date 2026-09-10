package jline.solvers.ln;

import jline.lang.constant.SchedStrategy;
import jline.VerboseLevel;
import jline.lang.constant.SolverType;
import jline.lang.layered.*;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.wrappers.lqns.SolverLQNS;
import jline.solvers.mva.MVAOptions;

import static jline.TestTools.*;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

/**
 * Tests for basic model building and activity graph precedence patterns.
 *
 * Contains 15 tests:
 * - test_buildModel_1 through test_buildModel_6 (6 tests)
 * - test_activityGraph_* (9 tests)
 */
class SolverLNActivityGraphTest extends SolverLNTestBase {

    @Test
    public void test_buildModel_1() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel1(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, 2.62479598084316, 2.62479598084316, 0.164049748802698, 2.46074623204047};
        double[] expectedUtil = {0.758003264230709, 0.758003264230709, 0.758003264230709, 0.0473752040144193, 0.71062806021629};
        double[] expectedRespT = {Double.NaN, Double.NaN, 5.54044259111638, 0.346277661944774, 5.1941649291716};
        double[] expectedResidT = {Double.NaN, 5.54044259111638, Double.NaN, 0.346277661944774, 5.1941649291716};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.473752040144193, 0.473752040144193, 0.473752040144193, 0.473752040144193};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void test_buildModel_2() throws Exception {
        LNOptions lnoptions = new LNOptions();
        lnoptions.verbose = VerboseLevel.SILENT;
        // the routing encoding, which this golden was recorded under; the 'srvn'
        // alias now resolves to 'srvn.ph' on this model, covered by the test below
        lnoptions.method = "srvn.cs";
        SolverOptions mvaoptions = new MVAOptions();
        mvaoptions.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel2(), SolverType.MVA, lnoptions, mvaoptions);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.124064286138323, 0.5325561421946614, 1.124064286138323, 0.5325561421946614, 0.1624589500751483, 0.9616053370785427, 0.4437967851622179, 0.08875935703244357};
        double[] expectedUtil = {0.1420149713913745, 0.5325561421946613, 0.1420149713913745, 0.5325561421946613, 0.1420149713913745, 0.5325561421946613, 0.1420149713913745, 0, 0.4437967851622178, 0.08875935703244356};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 12.66417786942251, 6.000000000000001, 1.830330405122518, 10.83384747573956, 5.000000000000001, 1};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.830330405122518, 6.000000000000001, Double.NaN, Double.NaN, 1.830330405122518, 0, 5.000000000000001, 1};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.08875935711960906, 0.08875935703244356, 0.08875935711960906, 0.08875935703244356, 0.08875935711960906, 0.08875935711960906, 0.08875935703244356, 0.08875935703244356};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void test_buildModel_2_srvnph() throws Exception {
        LNOptions lnoptions = new LNOptions();
        lnoptions.verbose = VerboseLevel.SILENT;
        SolverOptions mvaoptions = new MVAOptions();
        mvaoptions.verbose = VerboseLevel.SILENT;
        // the default: the 'srvn' alias resolves to 'srvn.ph' on this model
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel2(), SolverType.MVA, lnoptions, mvaoptions);
        assertEquals("srvn.ph", solver.lnmethod);
        // Ground truth from MATLAB SolverLN(SolverMVA) with method='srvn.ph',
        // recorded 2026-08-12. The phase-type encoding composes A1+A2 into one
        // entry service law, so it reaches a different fixed point from the
        // routing encoding above. RespT[4] is 13.29325240579065 before the
        // tenths snapping getAvgTable applies in both codebases.
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.173348996711299, 0.5295990624506023, 1.173348996711299, 0.5295990624506023, 0.1614300127033691, 1.011918983125265, 0.4413325520421686, 0.08826651040843372};
        double[] expectedUtil = {0.1412264160364762, 0.5295990624506023, 0.1412264160364762, 0.5295990624506023, 0.1412264160364762, 0.5295990624506023, 0.1412264160364762, 0, 0.4413325520421686, 0.08826651040843372};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 13.3, 6, 1.828893117691802, 11.46435927809885, 5, 1};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.828893117691802, 6, Double.NaN, Double.NaN, 1.828893117691802, 0, 5, 1};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.08826651002279765, 0.08826651040843372, 0.08826651002279765, 0.08826651040843372, 0.08826651002279765, 0.08826651002279765, 0.08826651040843372, 0.08826651040843372};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void test_buildModel_3() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        options.config.relax = "none";
        options.iter_max = 100;
        options.iter_tol = 0.0001;
        // the routing encoding, which this golden was recorded under; the 'srvn'
        // alias now resolves to 'srvn.ph' on this model, covered by the test below
        options.method = "srvn.cs";
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel3(), SolverType.MVA, options);
        // Ground truth from MATLAB SolverLN(SolverMVA), relax=none, iter_max=100,
        // iter_tol=1e-4. Re-recorded 2026-08-11, when the interlock probability
        // was aligned to Li and Franks (2015), Eq. (5), and to lqns' m' rule.
        double[] expectedQLen = {Double.NaN, 8.675385613675, 0.729887924585703, 8.675385613675, 0.729887924585703, 7.94549769054165, 0.729887926442931, 0.608239937411964, 0.121647987482393};
        double[] expectedUtil = {1.0067069317541, 0.927230068700173, 0.0794768630539232, 0.927230068700173, 0.0794768630539232, 0.927230068700173, 0, 0.0662307192116027, 0.0132461438423205};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, 654.9, 55.1, 599.8, 55.1, 45.9, 9.18365291291309};
        double[] expectedResidT = {Double.NaN, 599.8, 55.1, Double.NaN, Double.NaN, 599.8, 0, 45.9, 9.18365290902952};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.0132461438385739, 0.0132461438423205, 0.0132461438385739, 0.0132461438423205, 0.0132461438385739, 0.0132461438385739, 0.0132461438423205, 0.0132461438423205};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        avgTable.print();
        if (SolverLQNS.isAvailable()) {
            new SolverLQNS(SolverLNTestFixtures.buildModel3()).getAvgTable().print();
        }

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput, 0.01);
    }

    @Test
    public void test_buildModel_3_srvnph() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        options.config.relax = "none";
        options.iter_max = 100;
        options.iter_tol = 0.0001;
        // the default: the 'srvn' alias resolves to 'srvn.ph' on this model
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel3(), SolverType.MVA, options);
        assertEquals("srvn.ph", solver.lnmethod);
        // Ground truth from MATLAB SolverLN(SolverMVA) with method='srvn.ph',
        // relax=none, iter_max=100, iter_tol=1e-4, recorded 2026-08-11. The
        // phase-type encoding composes A3+A4 into one entry service law, so it
        // reaches a different fixed point from the routing encoding above.
        double[] expectedQLen = {Double.NaN, 8.68427038904056, 0.756637982659700, 8.68427038904056, 0.756637982659700, 7.92763240586455, 0.756637983044436, 0.630531652216416, 0.126106330443283};
        double[] expectedUtil = {1, 0.921010726732615, 0.0789437765539742, 0.921010726732615, 0.0789437765539742, 0.921010726732615, 0, 0.0657864804616452, 0.0131572960923290};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, 660, 57.5, 602.5, 57.5, 47.9, 9.58451718030468};
        double[] expectedResidT = {Double.NaN, 602.5, 57.5, Double.NaN, Double.NaN, 602.5, 0, 47.9, 9.58451718030468};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.0131572960961802, 0.0131572960923290, 0.0131572960961802, 0.0131572960923290, 0.0131572960961802, 0.0131572960961802, 0.0131572960923290, 0.0131572960923290};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput, 0.01);
    }

    @Test
    public void test_buildModel_4() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        options.config.relax = "none";
        options.iter_max = 100;
        options.iter_tol = 0.0001;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel4(), SolverType.MVA, options);
        // Ground truth values from MATLAB SolverLN with relax=none, iter_max=100, iter_tol=0.0001.
        // Re-recorded 2026-08-11, when the interlock probability was aligned to
        // Li and Franks (2015), Eq. (5), and to lqns' m' rule. Only the third
        // digit moves on this model -- E3 settles at 66.3 completions per second
        // rather than 66.4 (LDES, 500k samples, seed 23000: 66.744) -- because
        // its host layers take the interlock matrix inside the layer MVA rather
        // than the residence-time fallback.
        double[] expectedQLen = {Double.NaN, Double.NaN, 23.3, 8.64592492555925, 1.32667626070436, 23.3, 8.64592492555925, 1.32667626070436, 23.4323888862551, 8.71305264059124, 1.32667625061483};
        double[] expectedUtil = {0.99595480220186, 0.442225494260239, 0.66428467464315, 0.33167012755871, 0.442225494260239, 0.66428467464315, 0.33167012755871, 0.442225494260239, 0.66428467464315, 0.33167012755871, 0.442225494260239};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.75352816934047, 0.651696083484999, 0.01999999665214, 1.76373095607262, 0.656755908703967, 0.0199999965000377};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.08279559337584, 0.551424906310836, 0.01999999665214, Double.NaN, Double.NaN, Double.NaN, 1.08279559337584, 0.551424906310836, 0.01999999665214};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 13.285693492863, 13.2668051023484, 66.3, 13.285693492863, 13.2668051023484, 66.3, 13.285693492863, 13.2668051023484, 66.3};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput, 0.01);
    }

    @Test
    public void test_buildModel_5() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel5(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, 1.0000, 1.0000, 0.1587, 0.5238, 0.3175};
        double[] expectedUtil = {1.0000, 1.0000, 1.0000, 0.1587, 0.5238, 0.3175};
        double[] expectedRespT = {Double.NaN, Double.NaN, 12.6000, 2.0000, 3.0000, 4.0000};
        double[] expectedResidT = {Double.NaN, 12.6000, Double.NaN, 2.0000, 6.6000, 4.0000};
        double[] expectedTput = {Double.NaN, 0.0794, 0.0794, 0.0794, 0.1746, 0.0794};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            assertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
            assertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
            assertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
            assertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
            assertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @Test
    public void test_buildModel_6() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel6(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, 1.0, 1.0, 0.160000027162136, 0.0720000122229614, 0.0960000162972819, 0.192000032594564, 0.480000081486409};
        double[] expectedUtil = {1.0, 1.0, 1.0, 0.160000002128186, 0.0720000009576836, 0.0960000012769115, 0.192000002553823, 0.480000006384557};
        double[] expectedRespT = {Double.NaN, Double.NaN, 12.5, 2.0, 3.0, 4.0, 6.0, 6.0};
        double[] expectedResidT = {Double.NaN, 12.5, Double.NaN, 2.0, 0.9, 1.2, 2.4, 6.0};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.0800000010640929, 0.0800000010640929, 0.0800000010640929, 0.0240000003192279, 0.0240000003192279, 0.0320000004256372, 0.0800000010640929};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void test_activityGraph_and() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_and(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0000, 0.5479, 1.0000, 0.5479, 1.0000, 0.0000, 0.8769, 0.3563, 0.3563, 0.0000};
        double[] expectedUtil = {0.0000, 1.5895, 0.0000, 1.5895, 0.0000, 1.5895, 0.0000, 0.0000, 0.8770, 0.3563, 0.3563, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.3333, 1.3333, 1.3333, 0.0000, 2.0000, 1.0000, 1.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, 0.0000, 1.3333, Double.NaN, Double.NaN, 0.0000, 0.0000, 0.6666, 0.3333, 0.3333, 0.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.7500, 0.4110, 0.7500, 0.4110, 0.7500, 0.4932, 0.4385, 0.3563, 0.3563, 0.4110};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
            warningAssertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
            warningAssertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
            warningAssertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
            warningAssertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @Test
    public void test_activityGraph_call() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_call(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 0.1273996509487435, 0.1134380453455223, 0.06980802789275393, 0.1273996509487435, 0.1134380453455223, 0.06980802789275393, 0.1273996509487435, 0.1134380453455223, 0.06980802789275393};
        double[] expectedUtil = {0.01396160558342394, 0.04363001744058552, 0.06980802789275393, 0.01396160558342394, 0.04363001744058552, 0.06980802789275393, 0.01396160558342394, 0.04363001744058552, 0.06980802789275393, 0.01396160558342394, 0.04363001744058552, 0.06980802789275393};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 14.6, 13, 8, 14.6, 13, 8};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, 1.6, 5, 8, Double.NaN, Double.NaN, Double.NaN, 1.6, 5, 8};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.008726003489639965, 0.008726003488117104, 0.008726003486594242, 0.008726003489639965, 0.008726003488117104, 0.008726003486594242, 0.008726003489639965, 0.008726003488117104, 0.008726003486594242};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void test_activityGraph_call_and() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_call_and(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 1.0000, 0.7600, 0.9208, 1.0000, 0.7600, 0.9208, 1.0000, 0.6400, 0.1333, 0.0666, 0.0000, 0.9208};
        double[] expectedUtil = {0.0000, 0.2000, 0.9208, 0.0000, 0.2000, 0.9208, 0.0000, 0.2000, 0.9208, 0.0000, 0.0000, 0.1333, 0.0666, 0.0000, 0.9208};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 9.5000, 9.5000, 8.0000, 9.5000, 8.0000, 2.0000, 1.0000, 0.0000, 8.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, 0.0001, 1.5000, 8.0000, Double.NaN, Double.NaN, Double.NaN, 0.0001, 0.0000, 1.0000, 0.5000, 0.0000, 8.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.1053, 0.0800, 0.1151, 0.1053, 0.0800, 0.1151, 0.1053, 0.0800, 0.0666, 0.0666, 0.0800, 0.1151};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
            warningAssertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
            warningAssertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
            warningAssertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
            warningAssertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @Test
    public void test_activityGraph_call_or() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_call_or(), SolverType.MVA, options);
        // P3/T3/E3/A31 are unreachable: idle, not undefined. They report zero for
        // the measures their kind HAS and keep the NaN mask for the ones it never has.
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 1.0000, 0.9846, 0.0000, 1.0000, 0.9846, 0.0000, 1.0000, 0.0000, 0.0096, 0.0144, 0.9606, 0.0000};
        double[] expectedUtil = {0.0154, 0.9846, 0.0000, 0.0154, 0.9846, 0.0000, 0.0154, 0.9846, 0.0000, 0.0154, 0.0000, 0.0096, 0.0144, 0.9606, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 104.1000, 102.5000, 0.0000, 104.1000, 0.0000, 2.0000, 3.0000, 100.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, 1.6000, 102.5000, 0.0000, Double.NaN, Double.NaN, Double.NaN, 1.6000, 0.0000, 1.0000, 1.5000, 100.0000, 0.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.0096, 0.0096, 0.0000, 0.0096, 0.0096, 0.0000, 0.0096, 0.0096, 0.0048, 0.0048, 0.0096, 0.0000};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            assertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
            assertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
            assertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
            assertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
            assertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @Test
    public void test_activityGraph_call_seq_disconnected() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_call_seq_disconnected(), SolverType.MVA, options);
        // P3/T3/E3/A31 are unreachable: idle, not undefined. They report zero for
        // the measures their kind HAS and keep the NaN mask for the ones it never has.
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 0.07063197025365875, 0.05576208176883958, 0, 0.07063197025365875, 0.05576208176883958, 0, 0.07063197025365875, 0.04646840147403298, 0.009293680294806596, 0};
        double[] expectedUtil = {0.01486988847445447, 0.05576208176883957, 0, 0.01486988847445447, 0.05576208176883957, 0, 0.01486988847445447, 0.05576208176883957, 0, 0.01486988847445447, 0.04646840147403298, 0.009293680294806596, 0};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 7.6, 6, 0, 7.6, 5, 1, 0};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, 1.6, 6, 0, Double.NaN, Double.NaN, Double.NaN, 1.6, 5, 1, 0};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.009293680296534046, 0.009293680294806596, 0, 0.009293680296534046, 0.009293680294806596, 0, 0.009293680296534046, 0.009293680294806596, 0.009293680294806596, 0};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void test_activityGraph_loop() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_loop(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 0.0909091652727956, 0.0909091652480021, 0.0909091652727956, 0.0909091652480021, 0.0909091653637047, 0.0, 0.0818182569050276, 0.00909090834297452};
        double[] expectedUtil = {0.0, 0.0909091652480022, 0.0, 0.0909091652480022, 0.0, 0.0909091652480022, 0.0, 0.0, 0.0818182569050276, 0.00909090834297452};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 10.0, 10.0, 10.0, 1e-08, 3.0, 1.0};
        double[] expectedResidT = {Double.NaN, Double.NaN, 0.0, 10.0, Double.NaN, Double.NaN, 0.0, 0.0, 9.0, 1.0};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.00909090834545386, 0.00909090834297452, 0.00909090834545386, 0.00909090834297452, 0.00909090834545386, 0.00909090834297452, 0.0272727250289236, 0.00909090834297452};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void test_activityGraph_or() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_or(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0000, 0.9846, 1.0000, 0.9846, 1.0000, 0.0000, 0.0096, 0.0144, 0.9606};
        double[] expectedUtil = {0.0154, 0.9846, 0.0154, 0.9846, 0.0154, 0.9846, 0.0154, 0.0000, 0.0096, 0.0144, 0.9606};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 104.1000, 102.5000, 104.1000, 0.0000, 2.0000, 3.0000, 100.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.6000, 102.5000, Double.NaN, Double.NaN, 1.6000, 0.0000, 1.0000, 1.5000, 100.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0096, 0.0096, 0.0096, 0.0096, 0.0096, 0.0096, 0.0048, 0.0048, 0.0096};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            assertTrue(compareAbsErr(avg.getQLen().get(idx), expectedQLen[idx]));
            assertTrue(compareAbsErr(avg.getUtil().get(idx), expectedUtil[idx]));
            assertTrue(compareAbsErr(avg.getRespT().get(idx), expectedRespT[idx]));
            assertTrue(compareAbsErr(avg.getResidT().get(idx), expectedResidT[idx]));
            assertTrue(compareAbsErr(avg.getTput().get(idx), expectedTput[idx]));
        }
    }

    @Test
    public void test_activityGraph_seq() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(jline.solvers.ln.SolverLNTestFixtures.test_activityGraph_seq(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 0.07063197025365875, 0.05576208176883958, 0.07063197025365875, 0.05576208176883958, 0.07063197025365875, 0.04646840147403298, 0.009293680294806596};
        double[] expectedUtil = {0.01486988847445447, 0.05576208176883957, 0.01486988847445447, 0.05576208176883957, 0.01486988847445447, 0.05576208176883957, 0.01486988847445447, 0.04646840147403298, 0.009293680294806596};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 7.6, 6, 7.6, 5, 1};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.6, 6, Double.NaN, Double.NaN, 1.6, 5, 1};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.009293680296534046, 0.009293680294806596, 0.009293680296534046, 0.009293680294806596, 0.009293680296534046, 0.009293680294806596, 0.009293680294806596};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void test_cache_layer() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;

        LayeredNetwork cacheModel = SolverLNTestFixtures.buildCacheModel();
        SolverLN solver = new SolverLN(cacheModel, SolverType.MVA, options);

        assertNotNull(solver.getEnsemble());
        assertTrue(solver.getEnsemble().size() > 0);
        assertTrue(solver.getEnsemble().size() >= 1, "Cache model should create at least one layer");
    }
}
