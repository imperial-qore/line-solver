package jline.solvers.ssa;

import jline.examples.ClosedModel;
import jline.examples.MixedModel;
import jline.examples.OpenModel;
import jline.lang.Network;
import jline.lang.constant.GlobalConstants;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static jline.solvers.ssa.MatlabRand.allowedDeviation;
import static jline.solvers.ssa.MatlabRand.closeEnough;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class SolverSSAOpenExamplesTest {

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
    public void OpenExample1ReturnsCorrectAvgTable() {

        List<Double> QLen = Arrays.asList(0.022365, 0.11565, 0.0);
        List<Double> Util = Arrays.asList(0.022365, 0.10338, 0.0);
        List<Double> RespT = Arrays.asList(0.21634, 1.1215, 0.0);
        List<Double> ResidT = Arrays.asList(0.21634, 1.1215, 0.0);
        List<Double> ArvR = Arrays.asList(0.1, 0.10338, 0.0);
        List<Double> Tval = Arrays.asList(0.10338, 0.10312, 0.1);

        Network sn = OpenModel.ex1_line_v();
        SolverOptions options = new SolverOptions();
        options.keep = true;
        options.cutoff = 10;


        SolverSSA solver = new SolverSSA(sn, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        // Needed as we only consider non-zero rows
        List<Double> aQLen = avgTable.getQLen();
        List<Double> aUtil = avgTable.getUtil();
        List<Double> aRespT = avgTable.getRespT();
        List<Double> aResidT = avgTable.getResidT();
        List<Double> aArvR = avgTable.getArvR();
        List<Double> aTput = avgTable.getTput();

        for (int i = 0; i < aQLen.size(); i++) {
            if (aQLen.get(i) <= GlobalConstants.Zero && aUtil.get(i) <= GlobalConstants.Zero &&
                    aRespT.get(i) <= GlobalConstants.Zero && aResidT.get(i) <= GlobalConstants.Zero
                    && aArvR.get(i) <= GlobalConstants.Zero && aTput.get(i) <= GlobalConstants.Zero) {
                aQLen.remove(i);
                aUtil.remove(i);
                aRespT.remove(i);
                aResidT.remove(i);
                aArvR.remove(i);
                aTput.remove(i);
                i--;
            }
        }

        assertTrue(closeEnough(QLen, aQLen, allowedDeviation));
        assertTrue(closeEnough(Util, aUtil, allowedDeviation));
        assertTrue(closeEnough(RespT, aRespT, allowedDeviation));
        assertTrue(closeEnough(ResidT, aResidT, allowedDeviation));
        assertTrue(closeEnough(ArvR, aArvR, allowedDeviation));
        assertTrue(closeEnough(Tval, aTput, allowedDeviation));
    }

    @Test
    public void OpenExample3ReturnsCorrectAvgTable() {

        List<Double> QLen = Arrays.asList(0.0, 0.0, 1.1902, 0.9329, 0.77401);
        List<Double> Util = Arrays.asList(0.0, 0.0, 0.39602, 0.29981, 0.43654);
        List<Double> RespT = Arrays.asList(0.0, 0.0, 0.62834, 0.9182, 0.26459);
        List<Double> ResidT = Arrays.asList(0.0, 0.0, 0.41889, 0.30607, 0.26459);
        List<Double> ArvR = Arrays.asList(0.0, 0.0, 1.9801, 0.99938, 2.9103);
        List<Double> Tval = Arrays.asList(1.9801, 0.99938, 1.8942, 1.016, 2.9253);

        Network sn = OpenModel.ex3_line();
        SolverOptions options = new SolverOptions();
        options.keep = true;
        options.cutoff = 7;
        SolverSSA solver = new SolverSSA(sn, options);
        NetworkAvgTable avgTable = solver.getAvgTable();


        // Needed as we only consider non-zero rows
        List<Double> aQLen = avgTable.getQLen();
        List<Double> aUtil = avgTable.getUtil();
        List<Double> aRespT = avgTable.getRespT();
        List<Double> aResidT = avgTable.getResidT();
        List<Double> aArvR = avgTable.getArvR();
        List<Double> aTput = avgTable.getTput();

        for (int i = 0; i < aQLen.size(); i++) {
            if (aQLen.get(i) <= GlobalConstants.Zero && aUtil.get(i) <= GlobalConstants.Zero &&
                    aRespT.get(i) <= GlobalConstants.Zero && aResidT.get(i) <= GlobalConstants.Zero
                    && aArvR.get(i) <= GlobalConstants.Zero && aTput.get(i) <= GlobalConstants.Zero) {
                aQLen.remove(i);
                aUtil.remove(i);
                aRespT.remove(i);
                aResidT.remove(i);
                aArvR.remove(i);
                aTput.remove(i);
                i--;
            }
        }

        assertTrue(closeEnough(QLen, aQLen, allowedDeviation));
        assertTrue(closeEnough(Util, aUtil, allowedDeviation));
        assertTrue(closeEnough(RespT, aRespT, allowedDeviation));
        assertTrue(closeEnough(ResidT, aResidT, allowedDeviation));
        assertTrue(closeEnough(ArvR, aArvR, allowedDeviation));
        assertTrue(closeEnough(Tval, aTput, allowedDeviation));
    }

    @Test
    public void OpenExample6ReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(0.0, 0.0, 0.0, 0.28266, 0.25107, 0.31411, 0.22384, 0.13395, 0.19571, 0.25946, 0.20363, 0.19654, 0.31339, 0.15148, 0.28859);
        List<Double> Util = Arrays.asList(0.0, 0.0, 0.0, 0.16576, 0.1708, 0.2457, 0.13697, 0.10054, 0.13469, 0.19263, 0.13453, 0.14619, 0.18678, 0.069602, 0.20652);
        List<Double> RespT = Arrays.asList(0.0, 0.0, 0.0, 0.60156, 0.84801, 0.907, 1.7572, 1.8178, 1.9133, 2.7037, 2.9256, 2.7073, 2.4259, 2.0665, 3.1454);
        List<Double> ResidT = Arrays.asList(0.0, 0.0, 0.0, 2.4062, 3.3921, 3.628, 1.7572, 1.8178, 1.9133, 2.7037, 2.9256, 2.7073, 2.4259, 2.0665, 3.1454);
        List<Double> ArvR = Arrays.asList(0.0, 0.0, 0.0, 0.55253, 0.34159, 0.40949, 0.11747, 0.074017, 0.086579, 0.11747, 0.074017, 0.086579, 0.11747, 0.074017, 0.086579);
        List<Double> Tval = Arrays.asList(0.2, 0.125, 0.14286, 0.46988, 0.29607, 0.34632, 0.12738, 0.073688, 0.10229, 0.095963, 0.069604, 0.072595, 0.12919, 0.073301, 0.09175);

        Network sn = OpenModel.ex6();
        SolverOptions options = new SolverOptions();
        options.cutoff = 1;
        SolverSSA solver = new SolverSSA(sn, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        // Needed as we only consider non-zero rows
        List<Double> aQLen = avgTable.getQLen();
        List<Double> aUtil = avgTable.getUtil();
        List<Double> aRespT = avgTable.getRespT();
        List<Double> aResidT = avgTable.getResidT();
        List<Double> aArvR = avgTable.getArvR();
        List<Double> aTput = avgTable.getTput();

        for (int i = 0; i < aQLen.size(); i++) {
            if (aQLen.get(i) <= GlobalConstants.Zero && aUtil.get(i) <= GlobalConstants.Zero &&
                    aRespT.get(i) <= GlobalConstants.Zero && aResidT.get(i) <= GlobalConstants.Zero
                    && aArvR.get(i) <= GlobalConstants.Zero && aTput.get(i) <= GlobalConstants.Zero) {
                aQLen.remove(i);
                aUtil.remove(i);
                aRespT.remove(i);
                aResidT.remove(i);
                aArvR.remove(i);
                aTput.remove(i);
                i--;
            }
        }

        assertTrue(closeEnough(QLen, aQLen, allowedDeviation));
        assertTrue(closeEnough(Util, aUtil, allowedDeviation));
        assertTrue(closeEnough(RespT, aRespT, allowedDeviation));
        assertTrue(closeEnough(ResidT, aResidT, allowedDeviation));
        assertTrue(closeEnough(ArvR, aArvR, allowedDeviation));
        assertTrue(closeEnough(Tval, aTput, allowedDeviation));
    }
}
