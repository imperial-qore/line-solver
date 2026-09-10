/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Phase 0 proof-of-concept tests for in-place dense/sparse storage switching
 * (see jar/MATRIX-REFACTOR-PLAN.md).
 *
 * Verifies that toDense()/toSparse() swap the delegate in place and that the
 * dispatched operations (mult, add, sub, transpose, elementMult) and the
 * already-delegated ones (elementSum, sumRows, sumCols, get/set) return
 * identical results for sparse-sparse, dense-dense, and mixed operand formats.
 * Also includes a micro-benchmark comparing CSC and row-major backends at
 * full and 5% density to calibrate auto-switching thresholds.
 */
public class FormatSwitchPoCTest {

    private static final double TOL = 1e-12;

    // benchmark chatter is diagnostic only; keep the test run quiet unless
    // explicitly requested via -Djline.poc.verbose=true
    private static final boolean VERBOSE =
            Boolean.getBoolean("jline.poc.verbose");

    // deterministic sample matrix, mixed sign, roughly 50% fill
    private static Matrix sampleA(int n) {
        Matrix a = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((i + 2 * j) % 2 == 0) {
                    a.set(i, j, Math.sin(1.0 + i + n * j));
                }
            }
        }
        return a;
    }

    // deterministic sample matrix with a different sparsity pattern
    private static Matrix sampleB(int n) {
        Matrix b = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((i * j) % 3 != 1) {
                    b.set(i, j, Math.cos(1.0 + 2 * i + j));
                }
            }
        }
        return b;
    }

    private static Matrix withDensity(int n, double density) {
        Matrix m = new Matrix(n, n);
        int period = Math.max(1, (int) Math.round(1.0 / density));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((i + n * j) % period == 0) {
                    m.set(i, j, Math.sin(1.0 + i + n * j));
                }
            }
        }
        return m;
    }

    private static Matrix dense(Matrix m) {
        return new Matrix(m).toDense();
    }

    private static void assertMatrixEquals(Matrix expected, Matrix actual) {
        assertEquals(expected.getNumRows(), actual.getNumRows(), "row count");
        assertEquals(expected.getNumCols(), actual.getNumCols(), "col count");
        for (int i = 0; i < expected.getNumRows(); i++) {
            for (int j = 0; j < expected.getNumCols(); j++) {
                assertEquals(expected.get(i, j), actual.get(i, j), TOL, "(" + i + "," + j + ")");
            }
        }
    }

    @Test
    public void roundTripPreservesValues() {
        Matrix a = sampleA(17);
        Matrix r = new Matrix(a);
        assertTrue(r.isSparse());
        assertFalse(r.isDense());
        r.toDense();
        assertTrue(r.isDense());
        assertFalse(r.isSparse());
        assertMatrixEquals(a, r);
        r.toSparse();
        assertTrue(r.isSparse());
        assertMatrixEquals(a, r);
        // explicit zeros must be dropped when sparsifying
        assertEquals(a.getNonZeros(), r.getNonZeros());
    }

    @Test
    public void getSetAfterSwitch() {
        Matrix a = sampleA(9).toDense();
        a.set(3, 4, 42.0);
        assertEquals(42.0, a.get(3, 4), 0.0);
        a.set(3, 4, 0.0);
        assertEquals(0.0, a.get(3, 4), 0.0);
        a.toSparse();
        assertEquals(0.0, a.get(3, 4), 0.0);
    }

    @Test
    public void binaryOpsMatchAcrossFormats() {
        int n = 13;
        Matrix a = sampleA(n);
        Matrix b = sampleB(n);
        // reference results computed sparse-sparse
        Matrix mult = a.mult(b);
        Matrix add = a.add(2.5, b);
        Matrix sub = a.sub(0.5, b);
        // dense-dense
        assertMatrixEquals(mult, dense(a).mult(dense(b)));
        assertMatrixEquals(add, dense(a).add(2.5, dense(b)));
        assertMatrixEquals(sub, dense(a).sub(0.5, dense(b)));
        // mixed operands
        assertMatrixEquals(mult, dense(a).mult(b));
        assertMatrixEquals(mult, a.mult(dense(b)));
        assertMatrixEquals(add, dense(a).add(2.5, b));
        assertMatrixEquals(add, a.add(2.5, dense(b)));
        assertMatrixEquals(sub, a.sub(0.5, dense(b)));
        // scalar add on dense
        assertMatrixEquals(a.add(0.75), dense(a).add(0.75));
    }

    @Test
    public void elementWiseAndReductionsMatch() {
        int n = 11;
        Matrix a = sampleA(n);
        Matrix b = sampleB(n);
        Matrix emult = a.elementMult(b, null);
        assertMatrixEquals(emult, dense(a).elementMult(dense(b), null));
        assertMatrixEquals(emult, dense(a).elementMult(b, null));
        assertEquals(a.elementSum(), dense(a).elementSum(), TOL);
        assertMatrixEquals(a.sumRows(), dense(a).sumRows());
        assertMatrixEquals(a.sumCols(), dense(a).sumCols());
        assertMatrixEquals(a.transpose(), dense(a).transpose());
    }

    @Test
    public void benchmarkDenseVsSparse() {
        benchmarkAtDensity(200, 1.0);
        benchmarkAtDensity(500, 0.05);
        benchmarkFill(200);
    }

    // element-wise fill: the CSC worst case, O(nz) memmove per insertion
    private void benchmarkFill(int n) {
        long t0 = System.nanoTime();
        Matrix s = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                s.set(i, j, i + j + 1.0);
            }
        }
        long tSparse = System.nanoTime() - t0;
        t0 = System.nanoTime();
        Matrix d = new Matrix(n, n).toDense();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                d.set(i, j, i + j + 1.0);
            }
        }
        long tDense = System.nanoTime() - t0;
        assertMatrixEquals(s, d);
        if (VERBOSE) {
            System.out.printf("PoC benchmark fill n=%d: sparse %.1f ms, dense %.1f ms, ratio %.2fx%n",
                    n, tSparse / 1e6, tDense / 1e6, (double) tSparse / tDense);
        }
    }

    private void benchmarkAtDensity(int n, double density) {
        Matrix as = withDensity(n, density);
        Matrix bs = withDensity(n, density);
        Matrix ad = dense(as);
        Matrix bd = dense(bs);
        // warm-up and correctness check
        assertMatrixEquals(as.mult(bs), ad.mult(bd));
        int reps = 3;
        long t0 = System.nanoTime();
        for (int r = 0; r < reps; r++) {
            as.mult(bs).add(1.0, bs).transpose();
        }
        long tSparse = System.nanoTime() - t0;
        t0 = System.nanoTime();
        for (int r = 0; r < reps; r++) {
            ad.mult(bd).add(1.0, bd).transpose();
        }
        long tDense = System.nanoTime() - t0;
        if (VERBOSE) {
            System.out.printf("PoC benchmark n=%d density=%.2f: sparse %.1f ms, dense %.1f ms, ratio %.2fx%n",
                    n, density, tSparse / 1e6, tDense / 1e6, (double) tSparse / tDense);
        }
    }
}
