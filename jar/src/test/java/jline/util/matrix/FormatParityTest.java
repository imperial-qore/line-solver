/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import org.junit.jupiter.api.Test;

import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Characterization safety net for the dense/sparse storage refactoring
 * (MATRIX-REFACTOR-PLAN.md, Phase 1).
 *
 * Every operation is executed with sparse-backed, dense-backed, and mixed
 * operands and the results compared elementwise. Cases are added here as each
 * facade section is migrated to format-neutral dispatch in Phase 3. The
 * edge-semantics tests freeze current behavior (explicit zeros, getNonZeros
 * meaning per format) before any change is made.
 */
public class FormatParityTest {

    private static final double TOL = 1e-12;

    // deterministic fixtures

    private static Matrix square(int n) {
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

    private static Matrix square2(int n) {
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

    private static Matrix rect(int m, int n) {
        Matrix a = new Matrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if ((i + j) % 3 != 2) {
                    a.set(i, j, 1.0 + Math.sin(i + m * j));
                }
            }
        }
        return a;
    }

    private static Matrix dense(Matrix m) {
        return new Matrix(m).toDense();
    }

    private static Matrix sparse(Matrix m) {
        return new Matrix(m);
    }

    private static void assertMatrixEquals(Matrix expected, Matrix actual, String ctx) {
        assertEquals(expected.getNumRows(), actual.getNumRows(), ctx + ": row count");
        assertEquals(expected.getNumCols(), actual.getNumCols(), ctx + ": col count");
        for (int i = 0; i < expected.getNumRows(); i++) {
            for (int j = 0; j < expected.getNumCols(); j++) {
                assertEquals(expected.get(i, j), actual.get(i, j), TOL,
                        ctx + " at (" + i + "," + j + ")");
            }
        }
    }

    /**
     * Runs a producer op on sparse-backed and dense-backed copies of the
     * operand and asserts identical results.
     */
    private static void checkUnary(String name, Function<Matrix, Matrix> op, Matrix a) {
        Matrix ref = op.apply(sparse(a));
        assertMatrixEquals(ref, op.apply(dense(a)), name + " [dense]");
    }

    /**
     * Runs a scalar-valued op on sparse-backed and dense-backed copies.
     */
    private static void checkUnaryScalar(String name, Function<Matrix, Double> op, Matrix a) {
        double ref = op.apply(sparse(a));
        assertEquals(ref, op.apply(dense(a)), TOL, name + " [dense]");
    }

    /**
     * Runs a binary producer op with all four operand format combinations
     * and asserts identical results.
     */
    private static void checkBinary(String name, BiFunction<Matrix, Matrix, Matrix> op,
                                    Matrix a, Matrix b) {
        Matrix ref = op.apply(sparse(a), sparse(b));
        assertMatrixEquals(ref, op.apply(dense(a), dense(b)), name + " [dd]");
        assertMatrixEquals(ref, op.apply(dense(a), sparse(b)), name + " [ds]");
        assertMatrixEquals(ref, op.apply(sparse(a), dense(b)), name + " [sd]");
    }

    /**
     * Runs an in-place mutator on sparse-backed and dense-backed copies and
     * asserts the mutated matrices agree.
     */
    private static void checkInPlace(String name, Consumer<Matrix> op, Matrix a) {
        Matrix s = sparse(a);
        Matrix d = dense(a);
        op.accept(s);
        op.accept(d);
        assertMatrixEquals(s, d, name + " [in-place]");
    }

    // Section 3-4: arithmetic and element-wise (dispatched in Phase 0 PoC)

    @Test
    public void arithmeticParity() {
        Matrix a = square(12);
        Matrix b = square2(12);
        checkBinary("mult", (x, y) -> x.mult(y), a, b);
        checkBinary("mult(out=null)", (x, y) -> x.mult(y, null), a, b);
        checkBinary("add", (x, y) -> x.add(2.5, y), a, b);
        checkBinary("sub", (x, y) -> x.sub(0.5, y), a, b);
        checkBinary("elementMult", (x, y) -> x.elementMult(y, null), a, b);
        checkUnary("add(scalar)", x -> x.add(0.75), a);
        checkUnary("transpose", Matrix::transpose, a);
        checkUnary("transpose rect", Matrix::transpose, rect(7, 12));
    }

    @Test
    public void reductionParity() {
        Matrix a = square(11);
        checkUnaryScalar("elementSum", Matrix::elementSum, a);
        checkUnary("sumRows", Matrix::sumRows, a);
        checkUnary("sumCols", Matrix::sumCols, a);
    }

    // Section 2: element access

    @Test
    public void elementAccessParity() {
        Matrix a = square(9);
        checkUnaryScalar("get(i,j)", x -> x.get(4, 5), a);
        checkUnaryScalar("get(idx)", x -> x.get(31), a);
        checkInPlace("set nonzero", x -> x.set(2, 7, 3.25), a);
        checkInPlace("set zero", x -> x.set(0, 0, 0.0), a);
        checkInPlace("set int", x -> x.set(1, 1, 5), a);
    }

    // Section 7: slicing (extractRows/extractColumn had dense paths pre-PoC)

    @Test
    public void slicingParity() {
        Matrix a = rect(9, 13);
        checkUnary("extractRows", x -> Matrix.extractRows(x, 2, 6, null), a);
        checkUnary("extractColumn", x -> Matrix.extractColumn(x, 5, null), a);
        checkUnary("extractColumns", x -> Matrix.extractColumns(x, 3, 8), a);
    }

    // Frozen edge semantics (pre-refactor behavior, must not change silently)

    @Test
    public void explicitZeroSemantics() {
        // sparse: set(i,j,0) removes the entry
        Matrix s = square(6);
        int nzBefore = s.getNonZeros();
        assertTrue(s.get(0, 0) != 0.0);
        s.set(0, 0, 0.0);
        assertEquals(nzBefore - 1, s.getNonZeros(), "sparse set-zero removes entry");
        // dense: getNonZeros counts actual nonzero values
        Matrix d = dense(square(6));
        assertEquals(square(6).getNonZeros(), d.getNonZeros(), "dense counts values");
        d.set(0, 0, 0.0);
        assertEquals(nzBefore - 1, d.getNonZeros(), "dense set-zero reduces count");
    }

    // Section 3-5, 12-13: in-place mutators (migrated in Phase 3)

    @Test
    public void mutatorParity() {
        Matrix a = square(10);
        checkInPlace("apply equal", x -> x.apply(0.0, 7.0, "equal"), a);
        checkInPlace("apply notequal", x -> x.apply(0.0, 3.0, "notequal"), a);
        checkInPlace("apply great", x -> x.apply(0.2, 9.0, "great"), a);
        checkInPlace("apply less", x -> x.apply(0.2, -9.0, "less"), a);
        checkInPlace("addEq scalar", x -> x.addEq(1.5), a);
        checkInPlace("scaleEq", x -> x.scaleEq(2.5), a);
        checkInPlace("mulByMinusOne", Matrix::mulByMinusOne, a);
        checkInPlace("multEq", x -> x.multEq(square2(10)), a);
        checkInPlace("divideEq", x -> x.divideEq(4.0), a);
        checkInPlace("replace", x -> x.replace(0.0, 5.5), a);
        checkInPlace("remove", x -> x.remove(0, 0), a);
        checkInPlace("zero", Matrix::zero, a);
        checkInPlace("ones", Matrix::ones, a);
        checkInPlace("fill", x -> x.fill(3.25), a);
        checkInPlace("expandMatrix", x -> x.expandMatrix(14, 15, 40), a);
        checkInPlace("expandMatrixToSquare", Matrix::expandMatrixToSquare, rect(6, 9));
        checkInPlace("removeRows", x -> x.removeRows(java.util.Arrays.asList(1, 3)), a);
        checkInPlace("removeCols", x -> x.removeCols(java.util.Arrays.asList(0, 2)), a);
        checkInPlace("setTo", x -> x.setTo(square2(10)), a);
        checkInPlace("unsafeSet", x -> x.unsafeSet(2, 3, 1.25), a);
    }

    // Sections 6-8: linear algebra, slicing, concatenation

    @Test
    public void linearAlgebraParity() {
        Matrix a = square(8);
        // make well-conditioned: diagonally dominant
        for (int i = 0; i < 8; i++) {
            a.set(i, i, 10.0 + i);
        }
        Matrix b = square2(8);
        checkUnary("inv", x -> x.inv(), a);
        checkUnaryScalar("det", Matrix::det, a);
        checkUnaryScalar("rank", x -> (double) x.rank(), a);
        checkBinary("kron", (x, y) -> x.kron(y), square(4), square2(3));
        checkBinary("krons", (x, y) -> x.krons(y), square(4), square2(3));
        checkUnary("elementPower", x -> x.elementPower(2.0), square(7));
        checkBinary("solve", (x, y) -> {
            Matrix sol = new Matrix(0, 0);
            Matrix.solve(x, y, sol);
            return sol;
        }, a, b);
    }

    @Test
    public void slicingAndConcatParity() {
        Matrix a = rect(8, 11);
        Matrix b = rect(8, 11);
        checkBinary("concatColumns", (x, y) -> Matrix.concatColumns(x, y, null), a, b);
        checkBinary("concatRows", (x, y) -> Matrix.concatRows(x, y, null), a, b);
        checkUnary("getSlice", x -> x.getSlice(1, 6, 2, 9), a);
        checkUnary("getRow", x -> x.getRow(3), a);
        checkUnary("getColumn", x -> x.getColumn(5), a);
        checkUnary("getRowsFrom", x -> x.getRowsFrom(2), square(9));
        checkUnary("repmat", x -> x.repmat(2, 3), a);
        checkUnary("extract", x -> Matrix.extract(x, 1, 5, 2, 7), a);
        checkUnary("extractDiag", x -> {
            Matrix d = new Matrix(0, 0);
            Matrix.extractDiag(x, d);
            return d;
        }, square(9));
        checkUnary("transposeSquare", Matrix::transpose, square(9));
    }

    // Sections 9-11: statistics, properties, find/filter

    @Test
    public void statsAndFindParity() {
        Matrix a = square(9);
        checkUnaryScalar("elementMax", Matrix::elementMax, a);
        checkUnaryScalar("elementMin", Matrix::elementMin, a);
        checkUnaryScalar("count", x -> (double) x.count(0.0), a);
        checkUnary("countEachRow", x -> x.countEachRow(0.0), a);
        checkUnary("find", Matrix::find, a);
        checkUnary("findNonNegative", Matrix::findNonNegative, a);
        checkUnaryScalar("value", Matrix::value, square(1));
        checkUnaryScalar("get(idx) vec", x -> x.sumRows().get(4), a);
    }

    @Test
    public void formatPromotionRules() {
        Matrix a = square(6);
        Matrix b = square2(6);
        // isolate the promotion rule from the fill-based auto-switching
        Matrix.setAutoFormatSwitch(false);
        try {
            promotionRuleAssertions(a, b);
        } finally {
            Matrix.setAutoFormatSwitch(true);
        }
    }

    private void promotionRuleAssertions(Matrix a, Matrix b) {
        assertTrue(sparse(a).mult(sparse(b)).isSparse(), "sparse op sparse stays sparse");
        assertTrue(dense(a).mult(dense(b)).isDense(), "dense op dense stays dense");
        assertTrue(dense(a).mult(sparse(b)).isDense(), "mixed promotes to dense");
        assertTrue(sparse(a).mult(dense(b)).isDense(), "mixed promotes to dense");
        // in-place writes preserve the destination format
        Matrix out = new Matrix(6, 6);
        dense(a).mult(dense(b), out);
        assertTrue(out.isSparse(), "sparse out preserved on dense mult");
        Matrix outD = dense(new Matrix(6, 6));
        sparse(a).mult(sparse(b), outD);
        assertTrue(outD.isDense(), "dense out preserved on sparse mult");
    }

    // Phase 5: opt-in automatic format selection

    @Test
    public void autoFormatSwitching() {
        // full-fill fixtures
        Matrix full = new Matrix(8, 8);
        full.fill(1.0);
        Matrix nearlyEmptyDense = dense(new Matrix(20, 20));
        nearlyEmptyDense.set(0, 0, 1.0);

        // on by default
        assertTrue(Matrix.isAutoFormatSwitch(), "auto-switching defaults to on");
        assertTrue(full.mult(full).isDense(), "dense result above densify threshold");
        Matrix sparseProduct = nearlyEmptyDense.mult(nearlyEmptyDense);
        assertTrue(sparseProduct.isSparse(), "sparse result below sparsify threshold");
        // values unaffected by switching
        assertMatrixEquals(dense(full).mult(dense(full)), full.mult(full), "autoFormat values");

        // disabled: formats never change spontaneously
        Matrix.setAutoFormatSwitch(false);
        try {
            assertTrue(full.mult(full).isSparse(), "no switching when disabled");
        } finally {
            Matrix.setAutoFormatSwitch(true);
        }
    }

    @Test
    public void copyConstructorPreservesFormat() {
        Matrix s = square(5);
        assertTrue(new Matrix(s).isSparse());
        Matrix d = dense(square(5));
        assertTrue(new Matrix(d).isDense());
        assertMatrixEquals(s, new Matrix(d), "copy of dense");
    }
}
