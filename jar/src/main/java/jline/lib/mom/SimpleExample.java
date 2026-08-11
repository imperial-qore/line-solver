package jline.lib.mom;

import jline.lib.mom.util.MomUtils;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class SimpleExample {
    private SimpleExample() {}

    public static void main(String[] args) {
        System.out.println("=== MOM Solver Integration Test ===\n");

        // Test 1: Using jline.util.Maths functions
        System.out.println("Test 1: Binomial coefficient (10 choose 3)");
        double binCoeff = Maths.binomialCoeff(10, 3);
        System.out.println("Result: " + binCoeff);

        // Test 2: Using multichoose
        System.out.println("\nTest 2: Multichoose(3, 2)");
        Matrix multiChoose = Maths.multichoose(3.0, 2.0);
        System.out.println("Results (" + multiChoose.getNumRows() + " combinations):");
        for (int i = 0; i < multiChoose.getNumRows(); i++) {
            System.out.print("  [");
            for (int j = 0; j < multiChoose.getNumCols(); j++) {
                System.out.print((int) multiChoose.get(i, j));
                if (j < multiChoose.getNumCols() - 1) System.out.print(", ");
            }
            System.out.println("]");
        }

        // Test 3: Using Matrix operations
        System.out.println("\nTest 3: Matrix operations");
        Matrix matrix1 = new Matrix(2, 2);
        matrix1.set(0, 0, 1.0);
        matrix1.set(0, 1, 2.0);
        matrix1.set(1, 0, 3.0);
        matrix1.set(1, 1, 4.0);

        Matrix matrix2 = new Matrix(2, 2);
        matrix2.set(0, 0, 5.0);
        matrix2.set(0, 1, 6.0);
        matrix2.set(1, 0, 7.0);
        matrix2.set(1, 1, 8.0);

        Matrix result = matrix1.mult(matrix2);
        System.out.println("Matrix multiplication result:");
        for (int i = 0; i < result.getNumRows(); i++) {
            System.out.print("  [");
            for (int j = 0; j < result.getNumCols(); j++) {
                System.out.print(result.get(i, j));
                if (j < result.getNumCols() - 1) System.out.print(", ");
            }
            System.out.println("]");
        }

        // Test 4: MOM-specific utilities
        System.out.println("\nTest 4: MOM utilities");
        int[][] combinations = new int[][] {
                {2, 0, 0},
                {1, 1, 0},
                {0, 2, 0},
                {1, 0, 1},
                {0, 1, 1},
                {0, 0, 2}
        };

        System.out.println("Original combinations:");
        for (int[] c : combinations) {
            System.out.println("  " + java.util.Arrays.toString(c));
        }

        int[][] sorted = MomUtils.sortByNnzPos(combinations);
        System.out.println("\nSorted by non-zeros:");
        for (int[] c : sorted) {
            System.out.println("  " + java.util.Arrays.toString(c));
        }

        System.out.println("\n=== All tests completed successfully ===");
    }
}
