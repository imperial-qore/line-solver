package jline.solvers.ln;

import jline.examples.java.basic.LayeredModel;
import jline.VerboseLevel;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;

import static jline.TestTools.*;

/**
 * Regression tests comparing SolverLN against LQNS solver baseline.
 * NaN or 0.0 entries in LQNS are ignored.
 */
class SolverLNDiffTest extends SolverLNTestBase {

    private double allowedMaxRelDiffLQNS = 0.10;

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_1() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel1(), options);
        double[] expectedQLen = {Double.NaN, 2.6133, 2.6133, 0.1633, 2.4500};
        double[] expectedUtil = {0.7582, 0.0000, 0.0000, 0.0474, 0.7108};
        double[] expectedRespT = {Double.NaN, Double.NaN, 5.5148, 0.3447, 5.1701};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.4739, 0.4739, 0.4739, 0.4739};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_2() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel2(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.1188, 0.5328, 1.1188, 0.5328, 0.1625, 0.9563, 0.4440, 0.0888};
        double[] expectedUtil = {0.1421, 0.5329, 0.1421, 0.5329, Double.NaN, Double.NaN, 0.1421, 0.0000, 0.4441, 0.0888};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 12.6000, 6.0000, 1.8299, 10.7678, 5.0000, 1.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, 1.8299, 6.0000, Double.NaN, Double.NaN, 1.8299, 0.0000, 5.0000, 1.0000};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_3() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel3(), options);
        double[] expectedQLen = {Double.NaN, 8.5090, 0.8385, 8.5090, 0.8385, 7.6706, 0.8384, 0.6988, 0.1398};
        double[] expectedUtil = {1.1332, 0.0000, 0.0000, 0.0000, 0.0000, 1.0437, 0.0000, 0.0746, 0.0149};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, 570.7000, 56.2000, 514.5000, 56.2000, 46.9000, 9.3732};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.0149, 0.0149, 0.0149, 0.0149, 0.0149, 0.0149, 0.0149, 0.0149};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_4() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel4(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 19.7749, 8.9502, 1.5927, 19.7749, 8.9502, 1.5927, 19.7749, 8.9502, 1.5927};
        double[] expectedUtil = {2.2668, 1.5111, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 1.5113, 0.7555, 1.5111};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 1.3084, 0.5924, 0.0211, 1.3084, 0.5924, 0.0211};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 15.1000, 15.1000, 75.6000, 15.1000, 15.1000, 75.6000, 15.1000, 15.1000, 75.6000};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_5() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel5(), options);
        double[] expectedQLen = {Double.NaN, 1.0000, 1.0000, 0.1587, 0.5238, 0.3175};
        double[] expectedUtil = {1.0000, 0.0000, 0.0000, 0.1587, 0.5238, 0.3175};
        double[] expectedRespT = {Double.NaN, Double.NaN, 12.6000, 2.0000, 3.0000, 4.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.0794, 0.0794, 0.0794, 0.1746, 0.0794};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_buildModel_6() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.buildModel6(), options);
        double[] expectedQLen = {Double.NaN, 1.0000, 1.0000, 0.1653, 0.0744, 0.0992, 0.1653, 0.4959};
        double[] expectedUtil = {1.0323, 0.0000, 0.0000, 0.1706, 0.0768, 0.1024, 0.1706, 0.5119};
        double[] expectedRespT = {Double.NaN, Double.NaN, 11.7219, 1.9375, 2.9062, 3.8750, 4.8438, 5.8125};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, 0.0853, 0.0853, 0.0853, 0.0256, 0.0256, 0.0341, 0.0853};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_lqn_serial() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(LayeredModel.lqn_serial(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.1215, 0.5327, 1.1215, 0.5327, 0.1627, 0.9587, 0.4439, 0.0888};
        double[] expectedUtil = {0.1421, 0.5327, 0.0000, 0.0000, 0.0000, 0.0000, 0.1421, 0.0000, 0.4439, 0.0888};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 12.6310, 6.0000, 1.8327, 10.8000, 5.0000, 1.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888, 0.0888};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_and() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_and(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.0000, 0.7539, 0.3770, 0.3770, 0.0000};
        double[] expectedUtil = {0.0000, 1.5078, 0.0000, 1.5078, 0.0000, 1.5078, 0.0000, 0.0000, 0.7539, 0.3770, 0.3770, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 2.6528, 2.6528, 2.6528, 0.0000, 2.0000, 1.0000, 1.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770, 0.3770};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_call() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_call(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 0.1274, 0.1134, 0.0698, 0.1274, 0.1134, 0.0698, 0.1274, 0.1134, 0.0698};
        double[] expectedUtil = {0.0140, 0.0436, 0.0698, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0140, 0.0436, 0.0698};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 14.6000, 13.0000, 8.0000, 14.6000, 13.0000, 8.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.0087, 0.0087, 0.0087, 0.0087, 0.0087, 0.0087, 0.0087, 0.0087, 0.0087};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_call_and() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_call_and(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 1.0000, 1.0000, 0.7680, 1.0000, 1.0000, 0.7680, 1.0000, 0.7680, 0.1920, 0.0960, 0.0000, 0.7680};
        double[] expectedUtil = {0.0000, 0.2880, 0.7680, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1920, 0.0960, 0.0000, 0.7680};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 10.4167, 10.4167, 8.0000, 10.4167, 8.0000, 2.0000, 1.0000, 0.0000, 8.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960, 0.0960};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_call_or() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_call_or(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 1.0000, 0.9846, 0.0000, 1.0000, 0.9846, 0.0000, 1.0000, 0.0000, 0.0096, 0.0144, 0.9606, 0.0000};
        double[] expectedUtil = {0.0154, 0.9846, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0154, 0.0000, 0.0096, 0.0144, 0.9606, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 104.1000, 102.5000, Double.NaN, 104.1000, 0.0000, 2.0000, 3.0000, 100.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.0096, 0.0096, 0.0000, 0.0096, 0.0096, 0.0000, 0.0096, 0.0096, 0.0048, 0.0048, 0.0096, 0.0000};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_call_seq_disconnected() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_call_seq_disconnected(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, Double.NaN, 0.0706, 0.0558, 0.0000, 0.0706, 0.0558, 0.0000, 0.0706, 0.0465, 0.0093, 0.0000};
        double[] expectedUtil = {0.0149, 0.0558, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0149, 0.0465, 0.0093, 0.0000};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 7.6000, 6.0000, Double.NaN, 7.6000, 5.0000, 1.0000, 0.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, Double.NaN, 0.0093, 0.0093, 0.0000, 0.0093, 0.0093, 0.0000, 0.0093, 0.0093, 0.0093, 0.0000};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_loop() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_loop(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 0.0909, 0.0909, 0.0909, 0.0909, 0.0909, 0.0000, 0.0818, 0.0091};
        double[] expectedUtil = {0.0000, 0.0909, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0818, 0.0091};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 10.0000, 10.0000, 10.0000, 0.0000, 3.0000, 1.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0091, 0.0091, 0.0091, 0.0091, 0.0091, 0.0091, 0.0273, 0.0091};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_or() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_or(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 1.0000, 0.9846, 1.0000, 0.9846, 1.0000, 0.0000, 0.0096, 0.0144, 0.9606};
        double[] expectedUtil = {0.0154, 0.9846, 0.0000, 0.0000, 0.0000, 0.0000, 0.0154, 0.0000, 0.0096, 0.0144, 0.9606};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 104.1000, 102.5000, 104.1000, 0.0000, 2.0000, 3.0000, 100.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0096, 0.0096, 0.0096, 0.0096, 0.0096, 0.0096, 0.0048, 0.0048, 0.0096};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }

    @org.junit.jupiter.api.Test
    public void difftest_activityGraph_seq() throws Exception {
        SolverOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverLN solver = new SolverLN(SolverLNTestFixtures.test_activityGraph_seq(), options);
        double[] expectedQLen = {Double.NaN, Double.NaN, 0.0706, 0.0558, 0.0706, 0.0558, 0.0706, 0.0465, 0.0093};
        double[] expectedUtil = {0.0149, 0.0558, 0.0000, 0.0000, 0.0000, 0.0000, 0.0149, 0.0465, 0.0093};
        double[] expectedRespT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, 7.6000, 6.0000, 7.6000, 5.0000, 1.0000};
        double[] expectedResidT = {Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN};
        double[] expectedTput = {Double.NaN, Double.NaN, 0.0093, 0.0093, 0.0093, 0.0093, 0.0093, 0.0093, 0.0093};
        LayeredNetworkAvgTable avg = (LayeredNetworkAvgTable) solver.getEnsembleAvg();
        for (int idx = 0; idx < avg.getQLen().size(); idx++) {
            warningAssertTrue(compareRelErrPositive(avg.getQLen().get(idx), expectedQLen[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getUtil().get(idx), expectedUtil[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getRespT().get(idx), expectedRespT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getResidT().get(idx), expectedResidT[idx], allowedMaxRelDiffLQNS));
            warningAssertTrue(compareRelErrPositive(avg.getTput().get(idx), expectedTput[idx], allowedMaxRelDiffLQNS));
        }
    }
}
