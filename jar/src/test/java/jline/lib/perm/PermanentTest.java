package jline.lib.perm;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.*;

/**
 * Unit tests for matrix permanent computation implementations.
 * Tests verify correctness of various permanent algorithms against known values.
 */
public class PermanentTest {


    /**
     * Test 1: Simple 2x2 matrix with known permanent value.
     * Matrix: [[1, 2], [3, 4]]
     * Expected permanent: 1*4 + 2*3 = 10
     */
    @Test
    public void testPermanent2x2Matrix() {
        Matrix matrix = new Matrix(2, 2);
        matrix.set(0, 0, 1.0); matrix.set(0, 1, 2.0);
        matrix.set(1, 0, 3.0); matrix.set(1, 1, 4.0);

        double expected = 10.0;  // 1*4 + 2*3 = 10

        // Test exact algorithms
        Permanent perm = new Permanent(matrix, true);
        assertEquals(expected, perm.getValue(), LOOSE_FINE_TOL,
            "Permanent (inclusion-exclusion) failed for 2x2 matrix");

        NaivePermanent naive = new NaivePermanent(matrix, true);
        assertEquals(expected, naive.getValue(), LOOSE_FINE_TOL,
            "NaivePermanent failed for 2x2 matrix");

        RyzerPermanent ryzer = new RyzerPermanent(matrix, "graycode", true);
        assertEquals(expected, ryzer.getValue(), LOOSE_FINE_TOL,
            "RyzerPermanent (graycode) failed for 2x2 matrix");

        RyzerPermanent ryzerNaive = new RyzerPermanent(matrix, "naive", true);
        assertEquals(expected, ryzerNaive.getValue(), LOOSE_FINE_TOL,
            "RyzerPermanent (naive) failed for 2x2 matrix");
    }

    /**
     * Test 2: 3x3 matrix with all unique entries.
     * Matrix: [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
     * Manual calculation:
     * perm = 1*5*9 + 1*6*8 + 2*4*9 + 2*6*7 + 3*4*8 + 3*5*7 = 45 + 48 + 72 + 84 + 96 + 105 = 450
     */
    @Test
    public void testPermanent3x3MatrixUnique() {
        Matrix matrix = new Matrix(3, 3);
        matrix.set(0, 0, 1.0); matrix.set(0, 1, 2.0); matrix.set(0, 2, 3.0);
        matrix.set(1, 0, 4.0); matrix.set(1, 1, 5.0); matrix.set(1, 2, 6.0);
        matrix.set(2, 0, 7.0); matrix.set(2, 1, 8.0); matrix.set(2, 2, 9.0);

        double expected = 450.0;

        Permanent perm = new Permanent(matrix, true);
        assertEquals(expected, perm.getValue(), LOOSE_FINE_TOL,
            "Permanent (inclusion-exclusion) failed for 3x3 unique matrix");

        NaivePermanent naive = new NaivePermanent(matrix, true);
        assertEquals(expected, naive.getValue(), LOOSE_FINE_TOL,
            "NaivePermanent failed for 3x3 unique matrix");

        RyzerPermanent ryzer = new RyzerPermanent(matrix, "graycode", true);
        assertEquals(expected, ryzer.getValue(), LOOSE_FINE_TOL,
            "RyzerPermanent (graycode) failed for 3x3 unique matrix");
    }

    /**
     * Test 3: 3x3 matrix with two repeated columns.
     * Matrix: [[1, 1, 2], [2, 2, 3], [3, 3, 4]]
     * Columns 0 and 1 are identical.
     * This tests the multiplicity optimization in the Permanent algorithm.
     */
    @Test
    public void testPermanent3x3WithRepeatedColumns() {
        Matrix matrix = new Matrix(3, 3);
        matrix.set(0, 0, 1.0); matrix.set(0, 1, 1.0); matrix.set(0, 2, 2.0);
        matrix.set(1, 0, 2.0); matrix.set(1, 1, 2.0); matrix.set(1, 2, 3.0);
        matrix.set(2, 0, 3.0); matrix.set(2, 1, 3.0); matrix.set(2, 2, 4.0);

        // Compute using different algorithms - they should all agree
        Permanent perm = new Permanent(matrix, true);
        double permValue = perm.getValue();

        NaivePermanent naive = new NaivePermanent(matrix, true);
        assertEquals(permValue, naive.getValue(), LOOSE_FINE_TOL,
            "NaivePermanent should match Permanent for 3x3 matrix with repeated columns");

        RyzerPermanent ryzer = new RyzerPermanent(matrix, "graycode", true);
        assertEquals(permValue, ryzer.getValue(), LOOSE_FINE_TOL,
            "RyzerPermanent should match Permanent for 3x3 matrix with repeated columns");

        // Verify the value is positive
        assertTrue(permValue > 0, "Permanent should be positive for this matrix");
    }

    /**
     * Test 4: 3x3 matrix with two repeated rows.
     * Matrix: [[1, 2, 3], [1, 2, 3], [4, 5, 6]]
     * Rows 0 and 1 are identical.
     * This tests the row-based multiplicity optimization in the Permanent algorithm.
     */
    @Test
    public void testPermanent3x3WithRepeatedRows() {
        Matrix matrix = new Matrix(3, 3);
        matrix.set(0, 0, 1.0); matrix.set(0, 1, 2.0); matrix.set(0, 2, 3.0);
        matrix.set(1, 0, 1.0); matrix.set(1, 1, 2.0); matrix.set(1, 2, 3.0);
        matrix.set(2, 0, 4.0); matrix.set(2, 1, 5.0); matrix.set(2, 2, 6.0);

        // Compute using different algorithms - they should all agree
        Permanent perm = new Permanent(matrix, true);
        double permValue = perm.getValue();

        NaivePermanent naive = new NaivePermanent(matrix, true);
        assertEquals(permValue, naive.getValue(), LOOSE_FINE_TOL,
            "NaivePermanent should match Permanent for 3x3 matrix with repeated rows");

        RyzerPermanent ryzer = new RyzerPermanent(matrix, "graycode", true);
        assertEquals(permValue, ryzer.getValue(), LOOSE_FINE_TOL,
            "RyzerPermanent should match Permanent for 3x3 matrix with repeated rows");

        // Verify the value is positive
        assertTrue(permValue > 0, "Permanent should be positive for this matrix");
    }

    /**
     * Test 5: 3x3 matrix with both repeated rows and columns.
     * Matrix: [[1, 1, 2], [1, 1, 2], [3, 3, 4]]
     * Rows 0 and 1 are identical, columns 0 and 1 are identical.
     * This tests handling of both types of repetition simultaneously.
     */
    @Test
    public void testPermanent3x3WithRepeatedRowsAndColumns() {
        Matrix matrix = new Matrix(3, 3);
        matrix.set(0, 0, 1.0); matrix.set(0, 1, 1.0); matrix.set(0, 2, 2.0);
        matrix.set(1, 0, 1.0); matrix.set(1, 1, 1.0); matrix.set(1, 2, 2.0);
        matrix.set(2, 0, 3.0); matrix.set(2, 1, 3.0); matrix.set(2, 2, 4.0);

        // Compute using different algorithms - they should all agree
        Permanent perm = new Permanent(matrix, true);
        double permValue = perm.getValue();

        NaivePermanent naive = new NaivePermanent(matrix, true);
        assertEquals(permValue, naive.getValue(), LOOSE_FINE_TOL,
            "NaivePermanent should match Permanent for 3x3 matrix with repeated rows and columns");

        RyzerPermanent ryzer = new RyzerPermanent(matrix, "graycode", true);
        assertEquals(permValue, ryzer.getValue(), LOOSE_FINE_TOL,
            "RyzerPermanent should match Permanent for 3x3 matrix with repeated rows and columns");

        // Verify the value is positive
        assertTrue(permValue > 0, "Permanent should be positive for this matrix");
    }

    /**
     * Test 6: Identity matrix - permanent should equal 1.
     * For any nxn identity matrix, perm(I) = 1.
     */
    @Test
    public void testPermanentIdentityMatrix() {
        int n = 4;
        Matrix identity = Matrix.eye(n);

        double expected = 1.0;

        Permanent perm = new Permanent(identity, true);
        assertEquals(expected, perm.getValue(), LOOSE_FINE_TOL,
            "Permanent of identity matrix should be 1");

        NaivePermanent naive = new NaivePermanent(identity, true);
        assertEquals(expected, naive.getValue(), LOOSE_FINE_TOL,
            "NaivePermanent of identity matrix should be 1");

        RyzerPermanent ryzer = new RyzerPermanent(identity, "graycode", true);
        assertEquals(expected, ryzer.getValue(), LOOSE_FINE_TOL,
            "RyzerPermanent of identity matrix should be 1");
    }

    // ================================================================================
    // Heuristic Permanent Tests
    // ================================================================================

    /**
     * Test 7: HeuristicPermanent on small 2x2 matrix - compare with exact permanent.
     * Matrix: [[1, 2], [3, 4]]
     * Expected permanent: 1*4 + 2*3 = 10
     * The heuristic should give a reasonable approximation.
     */
    @Test
    public void testHeuristicPermanent2x2Matrix() {
        Matrix matrix = new Matrix(2, 2);
        matrix.set(0, 0, 1.0); matrix.set(0, 1, 2.0);
        matrix.set(1, 0, 3.0); matrix.set(1, 1, 4.0);

        double exact = 10.0;  // 1*4 + 2*3 = 10

        HeuristicPermanent heuristic = new HeuristicPermanent(matrix, true);
        double approx = heuristic.getValue();

        // Heuristic should be within 50% of exact value
        double relativeError = Math.abs(exact - approx) / exact;
        assertTrue(relativeError < 0.5,
            "HeuristicPermanent relative error too large for 2x2 matrix: " + relativeError);
    }

    /**
     * Test 8: HeuristicPermanent on 3x3 matrix - compare with exact permanent.
     * Matrix: [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
     * Expected permanent: 450
     */
    @Test
    public void testHeuristicPermanent3x3Matrix() {
        Matrix matrix = new Matrix(3, 3);
        matrix.set(0, 0, 1.0); matrix.set(0, 1, 2.0); matrix.set(0, 2, 3.0);
        matrix.set(1, 0, 4.0); matrix.set(1, 1, 5.0); matrix.set(1, 2, 6.0);
        matrix.set(2, 0, 7.0); matrix.set(2, 1, 8.0); matrix.set(2, 2, 9.0);

        Permanent perm = new Permanent(matrix, true);
        double exact = perm.getValue();

        HeuristicPermanent heuristic = new HeuristicPermanent(matrix, true);
        double approx = heuristic.getValue();

        // Heuristic should be within 50% of exact value
        double relativeError = Math.abs(exact - approx) / exact;
        assertTrue(relativeError < 0.5,
            "HeuristicPermanent relative error too large for 3x3 matrix: " + relativeError);
    }

    /**
     * Test 9: HeuristicPermanent on identity matrix - should be positive.
     * For any nxn identity matrix, perm(I) = 1.
     * Note: The heuristic is designed for dense positive matrices and may not
     * be accurate for sparse matrices like the identity matrix.
     */
    @Test
    public void testHeuristicPermanentIdentityMatrix() {
        // The identity is zero off the diagonal, so it has no full support.
        // The heuristic used to floor those zeros to 1e-15 before Sinkhorn and
        // return a positive number, and this test asserted only that
        // positivity -- which the floor guaranteed by construction and which
        // said nothing about accuracy. The floor is not invertible: it changes
        // the permanent by n! eps, O(1) by n = 18. On the identity the old
        // estimate was the van der Waerden bound n!/n^n, its worst case:
        // 3.63e-4 at n = 10 and 6.75e-17 at n = 40, against an exact 1.
        int n = 4;
        Matrix identity = Matrix.eye(n);

        assertThrows(IllegalArgumentException.class, () -> {
            new HeuristicPermanent(identity, true);
        }, "HeuristicPermanent should refuse a matrix without full support");
    }

    /**
     * Test 10: HeuristicPermanent on matrix of all ones.
     * For a matrix of all ones, permanent = n!
     */
    @Test
    public void testHeuristicPermanentOnesMatrix() {
        int n = 4;
        Matrix ones = Matrix.ones(n, n);

        // Exact permanent is n!
        double exact = 1.0;
        for (int i = 1; i <= n; i++) {
            exact *= i;
        }

        HeuristicPermanent heuristic = new HeuristicPermanent(ones, true);
        double approx = heuristic.getValue();

        // Heuristic should be within 50% of exact value
        double relativeError = Math.abs(exact - approx) / exact;
        assertTrue(relativeError < 0.5,
            "HeuristicPermanent relative error too large for ones matrix: " + relativeError);
    }

    /**
     * Test 11: HeuristicPermanent on diagonal matrix.
     * Permanent of diagonal matrix = product of diagonal elements.
     * Note: The heuristic is designed for dense positive matrices and may not
     * be accurate for sparse matrices like diagonal matrices.
     */
    @Test
    public void testHeuristicPermanentDiagonalMatrix() {
        int n = 5;
        Matrix diagonal = new Matrix(n, n);
        double product = 1.0;
        for (int i = 0; i < n; i++) {
            double val = i + 1.0;  // 1, 2, 3, 4, 5
            diagonal.set(i, i, val);
            product *= val;
        }

        // Zero off the diagonal, so no full support: refused rather than
        // floored. See testHeuristicPermanentIdentityMatrix.
        assertThrows(IllegalArgumentException.class, () -> {
            new HeuristicPermanent(diagonal, true);
        }, "HeuristicPermanent should refuse a matrix without full support");

        // The exact engine remains correct on the same matrix.
        double exact = new RyzerPermanent(diagonal, "graycode", true).getValue();
        assertEquals(product, exact, 1e-9,
            "The exact permanent of a diagonal matrix is the product of its diagonal");
    }

    /**
     * Test 12: HeuristicPermanent positivity test.
     * Permanent of a positive matrix should always be positive.
     */
    @Test
    public void testHeuristicPermanentPositivity() {
        Matrix matrix = new Matrix(6, 6);
        // Fill with random positive values
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                matrix.set(i, j, Math.random() + 0.1);
            }
        }

        HeuristicPermanent heuristic = new HeuristicPermanent(matrix, true);
        double approx = heuristic.getValue();

        assertTrue(approx > 0,
            "HeuristicPermanent should be positive for positive matrix");
    }

    /**
     * Test 13: HeuristicPermanent error handling - negative matrix.
     * Should throw IllegalArgumentException for matrices with negative elements.
     */
    @Test
    public void testHeuristicPermanentNegativeMatrix() {
        Matrix matrix = new Matrix(2, 2);
        matrix.set(0, 0, 1.0); matrix.set(0, 1, 2.0);
        matrix.set(1, 0, -1.0); matrix.set(1, 1, 4.0);  // Contains negative

        assertThrows(IllegalArgumentException.class, () -> {
            new HeuristicPermanent(matrix, true);
        }, "HeuristicPermanent should throw IllegalArgumentException for negative matrix");
    }

    /**
     * Test 14: HeuristicPermanent on doubly stochastic matrix.
     * For doubly stochastic matrices, permanent is bounded: n!/n^n <= perm <= n!
     */
    @Test
    public void testHeuristicPermanentDoublyStochastic() {
        int n = 5;
        Matrix matrix = Matrix.ones(n, n);
        // Make it doubly stochastic
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix.set(i, j, 1.0 / n);
            }
        }

        HeuristicPermanent heuristic = new HeuristicPermanent(matrix, true);
        double approx = heuristic.getValue();

        // Compute factorial
        double factorial = 1.0;
        for (int i = 1; i <= n; i++) {
            factorial *= i;
        }

        double lowerBound = factorial / Math.pow(n, n);
        double upperBound = factorial;

        // Allow 10% margin for numerical errors
        assertTrue(approx >= lowerBound * 0.9 && approx <= upperBound * 1.1,
            "HeuristicPermanent should be within van der Waerden bounds for doubly stochastic matrix. " +
            "Value: " + approx + ", bounds: [" + lowerBound + ", " + upperBound + "]");
    }

    /**
     * Test 15: HeuristicPermanent scaling property.
     * If A is scaled by factor c, permanent should scale by c^n.
     */
    @Test
    public void testHeuristicPermanentScaling() {
        Matrix matrix = new Matrix(4, 4);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                matrix.set(i, j, Math.random() + 0.1);
            }
        }

        double scale = 2.0;
        Matrix scaledMatrix = matrix.copy();
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                scaledMatrix.set(i, j, scale * matrix.get(i, j));
            }
        }

        HeuristicPermanent heuristic1 = new HeuristicPermanent(matrix, true);
        HeuristicPermanent heuristic2 = new HeuristicPermanent(scaledMatrix, true);

        double ratio = heuristic2.getValue() / heuristic1.getValue();
        double expectedRatio = Math.pow(scale, 4);

        // Allow 20% error for heuristic approximation
        double relativeError = Math.abs(expectedRatio - ratio) / expectedRatio;
        assertTrue(relativeError < 0.2,
            "HeuristicPermanent scaling property violated. Expected ratio: " + expectedRatio +
            ", actual ratio: " + ratio + ", relative error: " + relativeError);
    }

    /**
     * Test 16: HeuristicPermanent performance on larger matrix.
     * Should complete in reasonable time for 10x10 matrix.
     */
    @Test
    public void testHeuristicPermanentPerformance() {
        Matrix matrix = new Matrix(10, 10);
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                matrix.set(i, j, Math.random() + 0.1);
            }
        }

        long startTime = System.currentTimeMillis();
        HeuristicPermanent heuristic = new HeuristicPermanent(matrix, true);
        long endTime = System.currentTimeMillis();

        double approx = heuristic.getValue();

        assertTrue(approx > 0,
            "HeuristicPermanent should produce positive result");
        assertTrue(endTime - startTime < 1000,
            "HeuristicPermanent should complete in less than 1 second for 10x10 matrix");
    }

    /**
     * Test 17: SaddlePointPermanent is exact on a single column.
     *
     * At h == 1 no direction survives the homogeneity of the integrand, the
     * Laplace factor is empty, and the expansion returns the exact
     * n! prod_k a(k,0).
     */
    @Test
    public void testSaddlePointPermanentExactOnSingleColumn() {
        for (int n = 1; n <= 6; n++) {
            Matrix col = new Matrix(n, 1);
            double product = 1.0;
            for (int i = 0; i < n; i++) {
                double v = 0.1 + 0.13 * (i + 1);
                col.set(i, 0, v);
                product *= v;
            }
            double factorial = 1.0;
            for (int i = 2; i <= n; i++) {
                factorial *= i;
            }
            double expected = factorial * product;
            SaddlePointPermanent spm = new SaddlePointPermanent(col, new int[]{n}, true);
            assertEquals(expected, spm.getValue(), 1e-12 * expected,
                "The saddle point must be exact for a single column");
        }
    }

    /**
     * Test 18: SaddlePointPermanent closed form on the all-ones matrix.
     *
     * J_n scales to itself (xi = 1), so phi = n log n, the reduced Laplacian
     * I - J/n has determinant 1/n, and the estimate is (2 pi)^(-(n-1)/2)
     * n^(n+1/2). That is the analytic anchor the MATLAB, native Python and C++
     * twins are held to as well.
     */
    @Test
    public void testSaddlePointPermanentClosedFormOnAllOnes() {
        for (int n = 2; n <= 10; n++) {
            Matrix ones = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    ones.set(i, j, 1.0);
                }
            }
            double expected = Math.pow(2.0 * Math.PI, -0.5 * (n - 1)) * Math.pow(n, n + 0.5);
            SaddlePointPermanent spm = new SaddlePointPermanent(ones, true);
            assertEquals(expected, spm.getValue(), 1e-9 * expected,
                "The saddle point on J_" + n + " must match its closed form");
        }
    }

    /**
     * Test 19: SaddlePointPermanent lands nearer than the capacity it corrects.
     *
     * exp(getLogCapacity()) is the Gurvits capacity, an upper bound of the
     * permanent. The Gaussian factor is what turns that e^n-scale bound into a
     * usable estimate, so it must land strictly nearer the exact value.
     */
    @Test
    public void testSaddlePointPermanentBeatsItsCapacity() {
        java.util.Random random = new java.util.Random(4711);
        for (int n = 3; n <= 8; n++) {
            Matrix matrix = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    matrix.set(i, j, random.nextDouble() + 0.05);
                }
            }
            double reference = new RyzerPermanent(matrix, "graycode", true).getValue();
            SaddlePointPermanent spm = new SaddlePointPermanent(matrix, true);
            double capacity = Math.exp(spm.getLogCapacity());
            assertTrue(capacity >= reference * (1.0 - 1e-9),
                "The capacity must be an upper bound of the permanent");
            assertTrue(Math.abs(spm.getValue() - reference) < Math.abs(capacity - reference),
                "The Gaussian factor must bring the capacity nearer the permanent");
            assertEquals(Math.log(spm.getValue()), spm.getLogValue(), 1e-12,
                "getLogValue must be the logarithm of getValue");
            // measured bias at unit multiplicities, near (e/sqrt(2 pi))^n
            double bias = Math.pow(Math.E / Math.sqrt(2.0 * Math.PI), n);
            assertTrue(spm.getValue() >= reference && spm.getValue() <= reference * bias * 1.1,
                "The saddle point overestimates by about (e/sqrt(2 pi))^n at unit multiplicities");
        }
    }

    /**
     * Test 20: SaddlePointPermanent sharpens as the column multiplicities grow.
     *
     * With h fixed this is a genuine asymptotic expansion in min(m). The matrix
     * whose 2k rows all equal (a, b) with m = (k, k) has permanent
     * (2k)! (ab)^k, and the expansion returns exactly 4^k / (C(2k,k) sqrt(pi k))
     * times it: an analytic ratio, free of the matrix, decreasing to 1 like
     * 1 + 1/(8k).
     */
    @Test
    public void testSaddlePointPermanentSharpensWithMultiplicities() {
        final double a = 0.7;
        final double b = 1.3;
        double previous = Double.POSITIVE_INFINITY;
        for (int k = 1; k <= 8; k++) {
            Matrix rows = new Matrix(2 * k, 2);
            for (int i = 0; i < 2 * k; i++) {
                rows.set(i, 0, a);
                rows.set(i, 1, b);
            }
            double exact = 1.0;
            for (int i = 2; i <= 2 * k; i++) {
                exact *= i;
            }
            exact *= Math.pow(a * b, k);
            double estimate = new SaddlePointPermanent(rows, new int[]{k, k}, true).getValue();
            double ratio = estimate / exact;
            double binomial = 1.0;
            for (int i = 1; i <= k; i++) {
                binomial = binomial * (k + i) / i;
            }
            double closed = Math.pow(4.0, k) / (binomial * Math.sqrt(Math.PI * k));
            assertEquals(closed, ratio, 1e-9 * closed,
                "The ratio at m = (k,k) must match its closed form");
            assertTrue(Math.abs(ratio - (1.0 + 1.0 / (8.0 * k))) < 0.01 / k,
                "The ratio must approach 1 + 1/(8k)");
            assertTrue(ratio < previous, "The expansion must sharpen as k grows");
            previous = ratio;
        }
    }

    /**
     * Test 21: SaddlePointPermanent refuses what it cannot expand.
     *
     * A structural zero has no full support and flooring it is not invertible;
     * sum(m) != n makes the coefficient identically zero by homogeneity.
     */
    @Test
    public void testSaddlePointPermanentRefusals() {
        Matrix zeroRow = new Matrix(3, 3);
        zeroRow.set(0, 0, 1.0); zeroRow.set(0, 1, 2.0); zeroRow.set(0, 2, 3.0);
        zeroRow.set(2, 0, 4.0); zeroRow.set(2, 1, 5.0); zeroRow.set(2, 2, 6.0);
        assertThrows(RuntimeException.class, () -> new SaddlePointPermanent(zeroRow, true));

        Matrix negative = new Matrix(2, 2);
        negative.set(0, 0, -1.0); negative.set(0, 1, 2.0);
        negative.set(1, 0, 3.0); negative.set(1, 1, 4.0);
        assertThrows(RuntimeException.class, () -> new SaddlePointPermanent(negative, true));

        Matrix wide = new Matrix(4, 2);
        for (int i = 0; i < 4; i++) {
            wide.set(i, 0, 1.0);
            wide.set(i, 1, 1.0);
        }
        assertThrows(RuntimeException.class, () -> new SaddlePointPermanent(wide, true));
        assertThrows(RuntimeException.class,
            () -> new SaddlePointPermanent(wide, new int[]{1, 1}, true));

        // sum(m) = 4 = rows is accepted, and the columns are interchangeable here
        SaddlePointPermanent ok = new SaddlePointPermanent(wide, new int[]{2, 2}, true);
        assertTrue(ok.getValue() > 0.0, "A feasible multiplicity vector must be accepted");
    }
}
