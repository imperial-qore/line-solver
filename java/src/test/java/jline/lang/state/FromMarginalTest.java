package jline.lang.state;

import jline.examples.ClosedModel;
import jline.examples.TestModels;
import jline.lang.NetworkStruct;
import jline.util.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class FromMarginalTest {

    @Test
    public void returnsCorrectStatesClosedModel1() {
        int[][] m = new int[][]{{3}};
        int[][] correct = new int[][]{{3}};
        testFromMarginal(
                ClosedModel.example_closedModel_1().getStruct(true),
                0, m, correct);
        correct = new int[][]{{1,1,1}};
        testFromMarginal(
                ClosedModel.example_closedModel_1().getStruct(true),
                1, m, correct);
    }

    @Test
    public void returnsCorrectStatesClosedModel2() {
        int[][] m1 = new int[][]{{2, 2}};
        int[][] correct1 = new int[][]{
                {2, 0, 2, 0},
                {2, 0, 1, 1},
                {2, 0, 0, 2},
                {1, 1, 2, 0},
                {1, 1, 1, 1},
                {1, 1, 0, 2},
                {0, 2, 2, 0},
                {0, 2, 1, 1},
                {0, 2, 0, 2}
        };
        testFromMarginal(TestModels.test_closedModel_2().getStruct(true), 0, m1, correct1);

        int[][] m2 = new int[][]{{2, 2}};
        int[][] correct2 = new int[][]{{2, 0, 2}, {1, 1, 2}, {0, 2, 2}};
        testFromMarginal(TestModels.test_closedModel_2().getStruct(true), 1, m2, correct2);
    }

    @Test
    public void returnsCorrectStatesClosedModel3() {
        int[][] m = new int[][]{{1, 2, 1}};
        int[][] correct = new int[][]{{1, 0, 2, 1}, {0, 1, 2, 1}};
        testFromMarginal(TestModels.test_closedModel_3().getStruct(true), 1, m, correct);
    }

    @Test
    public void returnsCorrectStatesClosedModel4() {
        int[][] m = new int[][]{{2, 1, 1, 1}};
        int[][] correct = new int[][]{
                {4, 3, 2, 1, 0, 0, 0},
                {4, 3, 2, 0, 1, 0, 0},
                {4, 2, 2, 0, 0, 1, 0},
                {4, 1, 1, 1, 0, 1, 0},
                {4, 1, 1, 0, 1, 1, 0},
                {3, 4, 2, 1, 0, 0, 0},
                {3, 4, 2, 0, 1, 0, 0},
                {3, 2, 2, 0, 0, 0, 1},
                {3, 1, 1, 1, 0, 0, 1},
                {3, 1, 1, 0, 1, 0, 1},
                {2, 4, 2, 0, 0, 1, 0},
                {2, 3, 2, 0, 0, 0, 1},
                {2, 1, 1, 0, 0, 1, 1},
                {1, 4, 1, 1, 0, 1, 0},
                {1, 4, 1, 0, 1, 1, 0},
                {1, 3, 1, 1, 0, 0, 1},
                {1, 3, 1, 0, 1, 0, 1},
                {1, 2, 1, 0, 0, 1, 1},
                {1, 1, 0, 1, 0, 1, 1},
                {1, 1, 0, 0, 1, 1, 1}};

        testFromMarginal(TestModels.test_closedModel_4().getStruct(true), 1, m, correct);
    }

    @Test
    public void returnsCorrectStatesClosedModel7LCFSPR() {
        int[][] m = new int[][]{{2, 2}};
        int[][] correct = new int[][]{
                {2, 1, 2, 1, 1, 1, 1, 0},
                {2, 1, 1, 1, 2, 1, 1, 0},
                {2, 1, 1, 1, 1, 1, 0, 1},
                {1, 1, 2, 1, 2, 1, 1, 0},
                {1, 1, 2, 1, 1, 1, 0, 1},
                {1, 1, 1, 1, 2, 1, 0, 1}
        };
        testFromMarginal(ClosedModel.example_closedModel_7lcfspr().getStruct(true), 1, m, correct);
    }


    private void testFromMarginal(NetworkStruct sn, int ind, int[][] n, int[][] correct) {
        Matrix inputMatrix = new Matrix(n.length, n[0].length);
        inputMatrix.fromArray2D(n);

        Matrix correctMatrix = new Matrix(correct.length, correct[0].length);
        correctMatrix.fromArray2D(correct);

        Matrix res = State.fromMarginal(sn, ind, inputMatrix);
        assertTrue(res.isEqualTo(correctMatrix));
    }



}
