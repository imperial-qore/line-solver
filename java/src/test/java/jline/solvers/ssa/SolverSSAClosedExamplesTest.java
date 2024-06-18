package jline.solvers.ssa;

import jline.examples.ClosedModel;
import jline.lang.Network;
import jline.lang.constant.GlobalConstants;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.*;
import java.nio.file.FileSystems;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import static jline.solvers.ssa.MatlabRand.*;
import static org.junit.jupiter.api.Assertions.assertTrue;


public class SolverSSAClosedExamplesTest {

    @BeforeEach
    public void setUp() {
        Maths.setRandomNumbersMatlab(true);
    }

    @AfterEach
    public void clean() {
        Maths.setRandomNumbersMatlab(false);
    }

    @Test
    public void ClosedExample1ReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(2.1593, 7.8407);
        List<Double> Util = Arrays.asList(2.1593, 0.9717);
        List<Double> RespT = Arrays.asList(1.0, 11.764);
        List<Double> ResidT = Arrays.asList(1.0, 3.5292);
        List<Double> ArvR = Arrays.asList(2.178, 0.6478);
        List<Double> Tval = Arrays.asList(2.1593, 0.6665);

        Network sn = ClosedModel.ex1();
        SolverOptions options = new SolverOptions();
        options.seed = 1;
        SolverSSA solver = new SolverSSA(sn, options);
        NetworkAvgTable avgTable = solver.getAvgTable();


        assertTrue(closeEnough(QLen, avgTable.getQLen(), allowedDeviation));
        assertTrue(closeEnough(Util, avgTable.getUtil(), allowedDeviation));
        assertTrue(closeEnough(RespT, avgTable.getRespT(), allowedDeviation));
        assertTrue(closeEnough(ResidT, avgTable.getResidT(), allowedDeviation));
        assertTrue(closeEnough(ArvR, avgTable.getArvR(), allowedDeviation));
        assertTrue(closeEnough(Tval, avgTable.getTput(), allowedDeviation));
    }

    @Test
    public void ClosedExample2ReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(0.90513, 0.20221, 0.07763, 2.815);
        List<Double> Util = Arrays.asList(0.90513, 0.20221, 0.02597, 0.94854);
        List<Double> RespT = Arrays.asList(0.6622, 0.21318, 0.55506, 2.9625);
        List<Double> ResidT = Arrays.asList(0.39417, 0.086286, 0.033039, 1.1991);
        List<Double> ArvR = Arrays.asList(1.3882, 0.93199, 0.13668, 0.94854);
        List<Double> Tval = Arrays.asList(1.3668, 0.94854, 0.13986, 0.95021);

        Network sn = ClosedModel.ex2_line();
        SolverOptions options = new SolverOptions();
        options.seed = 1;
        SolverSSA solver = new SolverSSA(sn, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertTrue(closeEnough(QLen, avgTable.getQLen(), allowedDeviation));
        assertTrue(closeEnough(Util, avgTable.getUtil(), allowedDeviation));
        assertTrue(closeEnough(RespT, avgTable.getRespT(), allowedDeviation));
        assertTrue(closeEnough(ResidT, avgTable.getResidT(), allowedDeviation));
        assertTrue(closeEnough(ArvR, avgTable.getArvR(), allowedDeviation));
        assertTrue(closeEnough(Tval, avgTable.getTput(), allowedDeviation));

    }

    @Test
    public void ClosedExample3ReturnsCorrectAvgTable() {

        List<Double> QLen = Arrays.asList(1.1548, 0.24746, 0.73372, 0.028159, 0.56954, 0.26628);
        List<Double> Util = Arrays.asList(1.1548, 0.24746, 0.73372, 0.016565, 0.28674, 0.12229);
        List<Double> RespT = Arrays.asList(0.6623, 0.21575, 1.0, 0.1712, 0.51676, 0.34617);
        List<Double> ResidT = Arrays.asList(0.39423, 0.087329, 1.0, 0.010191, 0.20917, 0.34617);
        List<Double> ArvR = Arrays.asList(1.6581, 1.1778, 0.76921, 0.17437, 1.147, 0.73372);
        List<Double> Tval = Arrays.asList(1.7437, 1.147, 0.73372, 0.16448, 1.1021, 0.76921);

        Network sn = ClosedModel.ex3_line();
        SolverOptions options = new SolverOptions();
        options.seed = 1;
        SolverSSA solver = new SolverSSA(sn, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertTrue(closeEnough(QLen, avgTable.getQLen(), allowedDeviation));
        assertTrue(closeEnough(Util, avgTable.getUtil(), allowedDeviation));
        assertTrue(closeEnough(RespT, avgTable.getRespT(), allowedDeviation));
        assertTrue(closeEnough(ResidT, avgTable.getResidT(), allowedDeviation));
        assertTrue(closeEnough(ArvR, avgTable.getArvR(), allowedDeviation));
        assertTrue(closeEnough(Tval, avgTable.getTput(), allowedDeviation));

    }

    @Test
    public void ClosedExample4ReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(0.96691, 0.49397, 0.13524, 0.67856, 1.3888, 1.1503, 0.57277, 0.90066, 0.71276);
        List<Double> Util = Arrays.asList(0.96691, 0.49397, 0.13524, 0.67856, 0.32581, 0.3223, 0.033889, 0.22541, 0.22541);
        List<Double> RespT = Arrays.asList(1.0, 1.0, 0.1, 1.0, 1.3823, 3.0741, 0.5815, 1.3239, 1.9679);
        List<Double> ResidT = Arrays.asList(0.66667, 0.33333, 0.066667, 0.33333, 0.92153, 1.0247, 0.29075, 0.44129, 0.32799);
        List<Double> ArvR = Arrays.asList(1.0047, 0.37418, 1.3472, 0.68033, 0.97743, 0.48346, 1.0167, 0.67622, 0.33811);
        List<Double> Tval = Arrays.asList(0.96691, 0.49397, 1.3524, 0.67856, 1.0047, 0.37418, 0.98498, 0.68033, 0.36219);

        Network sn = ClosedModel.ex4_line();
        SolverOptions options = new SolverOptions();
        options.seed = 1;
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
    public void ClosedExample7PSReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(0.3111, 0.11933, 1.6889, 1.8807);
        List<Double> Util = Arrays.asList(0.3111, 0.11933, 0.45312, 0.53967);
        List<Double> RespT = Arrays.asList(0.68657, 0.22111, 3.6475, 3.5062);
        List<Double> ResidT = Arrays.asList(0.68657, 0.22111, 3.6475, 3.5062);
        List<Double> ArvR = Arrays.asList(0.46303, 0.53638, 0.45312, 0.53967);
        List<Double> Tval = Arrays.asList(0.45312, 0.53967, 0.46303, 0.53638);

        Network sn = ClosedModel.ex7_line_ps();
        SolverOptions options = new SolverOptions();
        options.seed = 1;
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
    public void ClosedExample7FCFSReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(0.31392, 0.12086, 1.6861, 1.8791);
        List<Double> Util = Arrays.asList(0.31392, 0.12086, 0.46169, 0.54902);
        List<Double> RespT = Arrays.asList(0.67994, 0.22014, 3.6571, 3.4916);
        List<Double> ResidT = Arrays.asList(0.67994, 0.22014, 3.6571, 3.4916);
        List<Double> ArvR = Arrays.asList(0.46105, 0.53819, 0.46169, 0.54902);
        List<Double> Tval = Arrays.asList(0.46169, 0.54902, 0.46105, 0.53819);

        Network sn = ClosedModel.ex7_line_fcfs();
        SolverOptions options = new SolverOptions();
        options.seed = 1;
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
    public void ClosedExample8ReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(2.0002, 1.648, 1.9998, 0.35204);
        List<Double> Util = Arrays.asList(2.0002, 1.648, 0.66672, 0.054932);
        List<Double> RespT = Arrays.asList(1.0, 1.0, 1.053, 0.21338);
        List<Double> ResidT = Arrays.asList(1.0, 1.0, 1.053, 0.21338);
        List<Double> ArvR = Arrays.asList(1.8991, 1.6499, 2.0002, 1.648);
        List<Double> Tval = Arrays.asList(2.0002, 1.648, 1.8991, 1.6499);

        Network sn = ClosedModel.ex8_line();
        SolverOptions options = new SolverOptions();
        options.seed = 1;
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
    public void ClosedExample9ReturnsCorrectAvgTable() {
        List<Double> QLen = Arrays.asList(0.17004, 0.16877, 0.48488, 0.49773, 9.3451, 9.3335);
        List<Double> Util = Arrays.asList(0.17004, 0.16877, 0.25505, 0.25316, 0.47219, 0.47765);
        List<Double> RespT = Arrays.asList(1.0, 1.0, 3.0806, 3.1261, 55.897, 56.259);
        List<Double> ResidT = Arrays.asList(1.0, 1.0, 3.0806, 3.1261, 55.897, 56.259);
        List<Double> ArvR = Arrays.asList(0.16718, 0.1659, 0.17004, 0.16877, 0.1574, 0.15922);
        List<Double> Tval = Arrays.asList(0.17004, 0.16877, 0.1574, 0.15922, 0.16718, 0.1659);

        Network sn = ClosedModel.ex9_line();
        SolverOptions options = new SolverOptions();
        options.samples = 5000;
        options.seed = 1;
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
