package jline.solvers.ssa;

import jline.examples.ClosedModel;
import jline.examples.MixedModel;
import jline.examples.OpenModel;
import jline.lang.Network;
import jline.lang.constant.GlobalConstants;
import jline.lang.constant.SolverType;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static jline.solvers.ssa.MatlabRand.*;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class SolverSSAMixedExamplesTest {



    // change random numbers to be Mersenne Twister generated to match matlab engine
    @BeforeEach
    public void setUp() {
        Maths.setRandomNumbersMatlab(true);
        Maths.setMatlabRandomSeed(1);
    }

    @AfterEach
    public void clean() {
        Maths.setRandomNumbersMatlab(false);
    }


    @Test
    public void MixedExample1ReturnsCorrectAvgTable() {

        List<Double> QLen = Arrays.asList(1.4456, 0.018912, 0.5544, 0.17113);
        List<Double> Util = Arrays.asList(1.4456, 0.018912, 0.40243, 0.08683, 0.0);
        List<Double> RespT = Arrays.asList(0.68251, 0.2178, 0.25093, 1.6718, 0.0);
        List<Double> ResidT = Arrays.asList(0.68251, 0.2178, 0.25093, 1.6718, 0.0);
        List<Double> ArvR = Arrays.asList(2.2094, 0.1, 2.1181, 0.08683, 0.0);
        List<Double> Tval = Arrays.asList(2.1181, 0.08683, 2.2094, 0.10236, 0.1);

        Network sn = MixedModel.ex1_line();
        SolverOptions options = new SolverOptions();
        options.keep = true;
        options.cutoff = 3;
        SolverSSA solver = new SolverSSA(sn, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        // Needed as we only consider non-zero rows
        for (int i = 0; i < avgTable.getQLen().size(); i++) {
            if (avgTable.getQLen().get(i) <= GlobalConstants.Zero && avgTable.getUtil().get(i) <= GlobalConstants.Zero &&
                    avgTable.getRespT().get(i) <= GlobalConstants.Zero && avgTable.getResidT().get(i) <= GlobalConstants.Zero
                    && avgTable.getArvR().get(i) <= GlobalConstants.Zero && avgTable.getTput().get(i) <= GlobalConstants.Zero) {
                avgTable.getQLen().remove(i);
                avgTable.getUtil().remove(i);
                avgTable.getRespT().remove(i);
                avgTable.getResidT().remove(i);
                avgTable.getArvR().remove(i);
                avgTable.getTput().remove(i);
                i--;
            }
        }

        assertTrue(closeEnough(QLen, avgTable.getQLen(), allowedDeviation));
        assertTrue(closeEnough(Util, avgTable.getUtil(), allowedDeviation));
        assertTrue(closeEnough(RespT, avgTable.getRespT(), allowedDeviation));
        assertTrue(closeEnough(ResidT, avgTable.getResidT(), allowedDeviation));
        assertTrue(closeEnough(ArvR, avgTable.getArvR(), allowedDeviation));
        assertTrue(closeEnough(Tval, avgTable.getTput(), allowedDeviation));

    }


    @Test
    public void MixedExample2ReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(2.2422, 1.0193, 0.35076, 0.19411, 0.23882, 0.15184, 0.16825, 0.11084, 0.0);
        List<Double> Util = Arrays.asList(0.67299, 0.27014, 0.17964, 0.090441, 0.074645, 0.049973, 0.044622, 0.023411, 0.0);
        List<Double> RespT = Arrays.asList(3.1204, 3.9846, 0.52211, 0.74754, 0.33451, 0.58009, 0.25, 0.44721, 0.0);
        List<Double> ResidT = Arrays.asList(3.1204, 3.9846, 0.52211, 0.74754, 0.33451, 0.58009, 0.25, 0.44721, 0.0);
        List<Double> ArvR = Arrays.asList(0.67299, 0.27014, 0.71856, 0.2558, 0.67181, 0.25967, 0.71395, 0.26175, 0.0);
        List<Double> Tval = Arrays.asList(0.71856, 0.2558, 0.67181, 0.25967, 0.71395, 0.26175, 0.67299, 0.24784, 0.27014);

        Network sn = MixedModel.ex2();
        SolverOptions options = new SolverOptions();
        options.keep = false;
        options.cutoff = 3;
        SolverSSA solver = new SolverSSA(sn, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        // Needed as we only consider non-zero rows
        for (int i = 0; i < avgTable.getQLen().size(); i++) {
            if (avgTable.getQLen().get(i) <= GlobalConstants.Zero && avgTable.getUtil().get(i) <= GlobalConstants.Zero &&
                    avgTable.getRespT().get(i) <= GlobalConstants.Zero && avgTable.getResidT().get(i) <= GlobalConstants.Zero
                    && avgTable.getArvR().get(i) <= GlobalConstants.Zero && avgTable.getTput().get(i) <= GlobalConstants.Zero) {
                avgTable.getQLen().remove(i);
                avgTable.getUtil().remove(i);
                avgTable.getRespT().remove(i);
                avgTable.getResidT().remove(i);
                avgTable.getArvR().remove(i);
                avgTable.getTput().remove(i);
                i--;
            }
        }

        assertTrue(closeEnough(QLen, avgTable.getQLen(), allowedDeviation));
        assertTrue(closeEnough(Util, avgTable.getUtil(), allowedDeviation));
        assertTrue(closeEnough(RespT, avgTable.getRespT(), allowedDeviation));
        assertTrue(closeEnough(ResidT, avgTable.getResidT(), allowedDeviation));
        assertTrue(closeEnough(ArvR, avgTable.getArvR(), allowedDeviation));
        assertTrue(closeEnough(Tval, avgTable.getTput(), allowedDeviation));
    }

    @Test
    public void MixedExample3ReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(0.6978, 0.66768, 0.22706, 0.19602, 0.13169, 0.16302, 0.10723, 0.13357, 0.0);
        List<Double> Util = Arrays.asList(0.42893, 0.3, 0.10885, 0.10102, 0.04941, 0.051973, 0.024691, 0.025206, 0.0);
        List<Double> RespT = Arrays.asList(1.6026, 2.3366, 0.5106, 0.72584, 0.33333, 0.57848, 0.25, 0.44721, 0.0);
        List<Double> ResidT = Arrays.asList(1.6026, 2.3366, 0.5106, 0.72584, 0.33333, 0.57848, 0.25, 0.44721, 0.0);
        List<Double> ArvR = Arrays.asList(0.42893, 0.3, 0.43541, 0.28574, 0.44469, 0.27006, 0.39506, 0.28181, 0.0);
        List<Double> Tval = Arrays.asList(0.43541, 0.28574, 0.44469, 0.27006, 0.39506, 0.28181, 0.42893, 0.29868, 0.3);

        Network sn = MixedModel.ex3();
        SolverOptions options = new SolverOptions();
        options.keep = false;
        options.cutoff = 3;

        SolverSSA solver = new SolverSSA(sn, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        // Needed as we only consider non-zero rows
        for (int i = 0; i < avgTable.getQLen().size(); i++) {
            if (avgTable.getQLen().get(i) <= GlobalConstants.Zero && avgTable.getUtil().get(i) <= GlobalConstants.Zero &&
                    avgTable.getRespT().get(i) <= GlobalConstants.Zero && avgTable.getResidT().get(i) <= GlobalConstants.Zero
                    && avgTable.getArvR().get(i) <= GlobalConstants.Zero && avgTable.getTput().get(i) <= GlobalConstants.Zero) {
                avgTable.getQLen().remove(i);
                avgTable.getUtil().remove(i);
                avgTable.getRespT().remove(i);
                avgTable.getResidT().remove(i);
                avgTable.getArvR().remove(i);
                avgTable.getTput().remove(i);
                i--;
            }
        }

        assertTrue(closeEnough(QLen, avgTable.getQLen(), allowedDeviation));
        assertTrue(closeEnough(Util, avgTable.getUtil(), allowedDeviation));
        assertTrue(closeEnough(RespT, avgTable.getRespT(), allowedDeviation));
        assertTrue(closeEnough(ResidT, avgTable.getResidT(), allowedDeviation));
        assertTrue(closeEnough(ArvR, avgTable.getArvR(), allowedDeviation));
        assertTrue(closeEnough(Tval, avgTable.getTput(), allowedDeviation));
    }

    @Test
    public void MixedExample5ReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(98.15, 2.8948, 1.0288, 0.051146, 0.49159, 0.026332, 0.32952, 0.0);
        List<Double> Util = Arrays.asList(0.98494, 0.028814, 0.48568, 0.020247, 0.32266, 0.01618, 0.24605, 0.0);
        List<Double> RespT = Arrays.asList(101.04, 101.1, 1.0628, 1.8251, 0.49948, 0.92472, 0.33456, 0.0);
        List<Double> ResidT = Arrays.asList(101.04, 101.1, 1.0628, 1.8251, 0.49948, 0.92472, 0.33456, 0.0);
        List<Double> ArvR = Arrays.asList(0.98494, 0.028814, 0.97135, 0.028633, 0.96799, 0.028024, 0.98421, 0.0);
        List<Double> Tval = Arrays.asList(0.97135, 0.028633, 0.96799, 0.028024, 0.98421, 0.028476, 0.98494, 0.028814);

        Network sn = MixedModel.ex5();
        SolverOptions options = new SolverOptions();
        options.keep = false;
        options.cutoff = 3;
        options.samples = 20000;

        SolverSSA solver = new SolverSSA(sn, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        // Needed as we only consider non-zero rows
        for (int i = 0; i < avgTable.getQLen().size(); i++) {
            if (avgTable.getQLen().get(i) <= GlobalConstants.Zero && avgTable.getUtil().get(i) <= GlobalConstants.Zero &&
                    avgTable.getRespT().get(i) <= GlobalConstants.Zero && avgTable.getResidT().get(i) <= GlobalConstants.Zero
                    && avgTable.getArvR().get(i) <= GlobalConstants.Zero && avgTable.getTput().get(i) <= GlobalConstants.Zero) {
                avgTable.getQLen().remove(i);
                avgTable.getUtil().remove(i);
                avgTable.getRespT().remove(i);
                avgTable.getResidT().remove(i);
                avgTable.getArvR().remove(i);
                avgTable.getTput().remove(i);
                i--;
            }
        }

        assertTrue(closeEnough(QLen, avgTable.getQLen(), allowedDeviation));
        assertTrue(closeEnough(Util, avgTable.getUtil(), allowedDeviation));
        assertTrue(closeEnough(RespT, avgTable.getRespT(), allowedDeviation));
        assertTrue(closeEnough(ResidT, avgTable.getResidT(), allowedDeviation));
        assertTrue(closeEnough(ArvR, avgTable.getArvR(), allowedDeviation));
        assertTrue(closeEnough(Tval, avgTable.getTput(), allowedDeviation));
    }

}
