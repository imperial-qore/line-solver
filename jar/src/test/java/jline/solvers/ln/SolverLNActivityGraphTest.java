package jline.solvers.ln;

import jline.lang.constant.SchedStrategy;
import jline.VerboseLevel;
import jline.lang.constant.SolverType;
import jline.lang.layered.*;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.lqns.SolverLQNS;
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
        double[] expectedUtil = {0.758003264230709, 0.758003264230709, Double.NaN, 0.0473752040144193, 0.71062806021629};
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
        SolverOptions mvaoptions = new MVAOptions();
        mvaoptions.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel2(), SolverType.MVA, lnoptions, mvaoptions);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.124064286138323, 0.5325561421946614, 1.124064286138323, 0.5325561421946614, 0.1624589500751483, 0.9616053370785427, 0.4437967851622179, 0.08875935703244357};
        double[] expectedUtil = {0.1420149713913745, 0.5325561421946613, 0.1420149713913745, 0.5325561421946613, Double.NaN, Double.NaN, 0.1420149713913745, 0, 0.4437967851622178, 0.08875935703244356};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 12.66417786942251, 6.000000000000001, 1.830330405122518, 10.83384747573956, 5.000000000000001, 1};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.830330405122518, 6.000000000000001, Double.NaN, Double.NaN, 1.830330405122518, 0, 5.000000000000001, 1};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.08875935711960906, 0.08875935703244356, 0.08875935711960906, 0.08875935703244356, 0.08875935711960906, 0.08875935711960906, 0.08875935703244356, 0.08875935703244356};
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
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel3(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, 8.684302439121234, 0.7566547562277667, 8.684302439121234, 0.7566547562277667, 7.92765158428176, 0.7566508559719956, 0.6305456301898056, 0.1261091260379611};
        double[] expectedUtil = {0.9999545067535648, 0.9210092588749577, 0.07894524787860714, Double.NaN, Double.NaN, 0.9210092588749577, 0, 0.06578770656550595, 0.01315754131310119};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, 660.0380668062525, 57.50730613130224, 602.5298937576315, 57.50817313469662, 47.92275510941854, 9.584551021883705};
        double[] expectedResidT = {Double.NaN, 602.5298937576315, 57.50730613130224, Double.NaN, Double.NaN, 602.5298937576315, 0, 47.92275510941854, 9.584551021883705};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.01315727512678511, 0.01315754131310119, 0.01315727512678511, 0.01315754131310119, 0.01315727512678511, 0.01315727512678511, 0.01315754131310119, 0.01315754131310119};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        avgTable.print();
        if (SolverLQNS.isAvailable()) {
            new SolverLQNS(SolverLNTestFixtures.buildModel3()).getAvgTable().print();
        }

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void test_buildModel_4() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        options.config.relax = "none";
        options.iter_max = 100;
        options.iter_tol = 0.0001;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel4(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 23.4682450470894, 8.67822070268364, 0.120491297980306, 23.4682450470894, 8.67822070268364, 0.120491297980306, 23.4388606730313, 8.66725245115206, 0.124378114044187};
        double[] expectedUtil = {0.995954019616073, 0.0414593696959083, 0.664070430319683, 0.331883589296389, 0.0414593696959083, Double.NaN, Double.NaN, Double.NaN, 0.664070430319683, 0.331883589296389, 0.0414593696959083};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.76699970180812, 0.653709687866906, 0.0193750007720929, 1.76478725771209, 0.652883475613173, 0.0200000007969991};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.09193140773878, 0.556095776529935, 0.0200000007969991, Double.NaN, Double.NaN, Double.NaN, 1.09193140773878, 0.556095776529935, 0.0200000007969991};
        double[] expectedArvR = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 13.2814086063937, 13.2753435718556, 6.21890545438625, 13.2814086063937, 13.2753435718556, 6.21890545438625, 13.2814086063937, 13.2753435718556, 6.21890545438625};
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void test_buildModel_5() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel5(), SolverType.MVA, options);
        double[] expectedQLen = {Double.NaN, 1.0000, 1.0000, 0.1587, 0.5238, 0.3175};
        double[] expectedUtil = {1.0000, 1.0000, Double.NaN, 0.1587, 0.5238, 0.3175};
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
        double[] expectedUtil = {1.0, 1.0, Double.NaN, 0.160000002128186, 0.0720000009576836, 0.0960000012769115, 0.192000002553823, 0.480000006384557};
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
        double[] expectedUtil = {0.0000, 1.5895, 0.0000, 1.5895, Double.NaN, Double.NaN, 0.0000, 0.0000, 0.8770, 0.3563, 0.3563, 0.0000};
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
        double[] expectedUtil = {0.01396160558342394, 0.04363001744058552, 0.06980802789275393, 0.01396160558342394, 0.04363001744058552, 0.06980802789275393, Double.NaN, Double.NaN, Double.NaN, 0.01396160558342394, 0.04363001744058552, 0.06980802789275393};
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
        double[] expectedUtil = {0.0000, 0.2000, 0.9208, 0.0000, 0.2000, 0.9208, Double.NaN, Double.NaN, Double.NaN, 0.0000, 0.0000, 0.1333, 0.0666, 0.0000, 0.9208};
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
        double[] expectedQLen = {Double.NaN, Double.NaN, 0.0000, 1.0000, 0.9846, 0.0000, 1.0000, 0.9846, 0.0000, 1.0000, 0.0000, 0.0096, 0.0144, 0.9606, 0.0000};
        double[] expectedUtil = {0.0154, 0.9846, 0.0000, 0.0154, 0.9846, 0.0000, Double.NaN, Double.NaN, 0.0000, 0.0154, 0.0000, 0.0096, 0.0144, 0.9606, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, 0.0000, Double.NaN, Double.NaN, 0.0000, 104.1000, 102.5000, 0.0000, 104.1000, 0.0000, 2.0000, 3.0000, 100.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, 0.0000, 1.6000, 102.5000, 0.0000, Double.NaN, Double.NaN, 0.0000, 1.6000, 0.0000, 1.0000, 1.5000, 100.0000, 0.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0000, 0.0096, 0.0096, 0.0000, 0.0096, 0.0096, 0.0000, 0.0096, 0.0096, 0.0048, 0.0048, 0.0096, 0.0000};
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
        double[] expectedQLen = {Double.NaN, Double.NaN, 0, 0.07063197025365875, 0.05576208176883958, 0, 0.07063197025365875, 0.05576208176883958, 0, 0.07063197025365875, 0.04646840147403298, 0.009293680294806596, 0};
        double[] expectedUtil = {0.01486988847445447, 0.05576208176883957, 0, 0.01486988847445447, 0.05576208176883957, 0, Double.NaN, Double.NaN, 0, 0.01486988847445447, 0.04646840147403298, 0.009293680294806596, 0};
        double[] expectedRespT = {Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0, 7.6, 6, 0, 7.6, 5, 1, 0};
        double[] expectedResidT = {Double.NaN, Double.NaN, 0, 1.6, 6, 0, Double.NaN, Double.NaN, 0, 1.6, 5, 1, 0};
        double[] expectedArvR = {Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, 0, Double.NaN, Double.NaN, Double.NaN, 0};
        double[] expectedTput = {Double.NaN, Double.NaN, 0, 0.009293680296534046, 0.009293680294806596, 0, 0.009293680296534046, 0.009293680294806596, 0, 0.009293680296534046, 0.009293680294806596, 0.009293680294806596, 0};
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
        double[] expectedUtil = {0.0, 0.0909091652480022, 0.0, 0.0909091652480022, Double.NaN, Double.NaN, 0.0, 0.0, 0.0818182569050276, 0.00909090834297452};
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
        double[] expectedUtil = {0.0154, 0.9846, 0.0154, 0.9846, Double.NaN, Double.NaN, 0.0154, 0.0000, 0.0096, 0.0144, 0.9606};
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
        double[] expectedUtil = {0.01486988847445447, 0.05576208176883957, 0.01486988847445447, 0.05576208176883957, Double.NaN, Double.NaN, 0.01486988847445447, 0.04646840147403298, 0.009293680294806596};
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
