package jline.lang.state;

import jline.examples.ClosedModel;
import jline.examples.TestModels;
import jline.lang.NetworkStruct;
import jline.util.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class FromMarginalAndStartedTest {

    @Test
    public void returnsCorrectStatesClosedModel1() {
        int[][] n = new int[][]{{3}};
        int[][] s = new int[][]{{3}};
        int[][] correct = new int[][]{{3}};

        testFromMarginalAndStarted(
                ClosedModel.example_closedModel_1().getStruct(true),
                0, n, s, correct);
        correct = new int[][]{{1,1,1,0}};
        s = new int[][]{{0}};
        testFromMarginalAndStarted(
                ClosedModel.example_closedModel_1().getStruct(true),
                1, n, s,correct);
    }

    @Test
    public void returnsCorrectStatesClosedModel2() {
        int[][] n1 = new int[][]{{2, 2}};
        int[][] s = new int[][]{{0,1}};
        int[][] correct1 = new int[][]{
                {2, 0, 2, 0}
        };

        testFromMarginalAndStarted(TestModels.test_closedModel_2().getStruct(true), 0, n1, s, correct1);

        int[][] n2 = new int[][]{{2, 2}};
        int[][] correct2 = new int[][]{{2, 0, 2}};
        testFromMarginalAndStarted(TestModels.test_closedModel_2().getStruct(true), 1, n2, s, correct2);
    }

    @Test
    public void returnsCorrectStatesClosedModel3() {
        int[][] m = new int[][]{{1, 2, 1}};
        int[][] s = new int[][]{{1,1,0}};
        int[][] correct = new int[][]{{1, 0, 2, 1}};

        testFromMarginalAndStarted(TestModels.test_closedModel_3().getStruct(true), 1, m, s, correct);
    }

    @Test
    public void returnsCorrectStatesClosedModel4() {
        int[][] m = new int[][]{{2, 1, 1, 1}};
        int[][] s = new int[][]{{2,1,0,0}};
        int[][] correct = new int[][]{
                {4, 3,2,1,0,0,0},
                {3,4,2,1,0,0,0}};
        testFromMarginalAndStarted(TestModels.test_closedModel_4().getStruct(true), 1, m, s, correct);
    }

    @Test
    public void returnsCorrectStatesClosedModel7LCFSPR() {
        int[][] m = new int[][]{{2, 2}};
        int[][] s = new int[][] {{0, 1}};
        int[][] correct = new int[][] {
                {2, 1, 1, 1, 1, 1, 0, 1},
                {1, 1, 2, 1, 1, 1, 0, 1},
                {1, 1, 1, 1, 2, 1, 0, 1}
        };
        testFromMarginalAndStarted(ClosedModel.example_closedModel_7lcfspr().getStruct(true), 1, m, s, correct);
    }


    private void testFromMarginalAndStarted(NetworkStruct sn, int ind, int[][] n, int[][] s, int[][] correct) {
        Matrix matrixn = new Matrix(n.length, n[0].length);
        matrixn.fromArray2D(n);

        Matrix matrixs = new Matrix(s.length, s[0].length);
        matrixs.fromArray2D(s);

        Matrix correctMatrix = new Matrix(correct.length, correct[0].length);
        correctMatrix.fromArray2D(correct);
        Matrix res = State.fromMarginalAndStarted(sn, ind, matrixn, matrixs);
//        System.out.println(res);

        assertTrue(res.isEqualTo(correctMatrix));
    }


}
