package jline.util.matrix;

import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Coverage tests for matrix utility classes with 0% coverage.
 */
public class MatrixUtilsCoverageTest {

    private static final double TOL = 1e-6;

    // ========== ComplexMatrix ==========

    @Test
    void testComplexMatrix() {
        // Test dimension constructor
        ComplexMatrix cm = new ComplexMatrix(2, 3);
        assertEquals(2, cm.getNumRows());
        assertEquals(3, cm.getNumCols());
        assertEquals(6, cm.getNumElements());
        assertFalse(cm.isEmpty());

        // Test set/get with Complex values
        cm.set(0, 0, new Complex(1.0, 2.0));
        Complex val = cm.get(0, 0);
        assertEquals(1.0, val.getReal(), TOL);
        assertEquals(2.0, val.getImaginary(), TOL);

        // Test set/get with double (real only)
        cm.set(1, 1, 5.0);
        assertEquals(5.0, cm.get(1, 1).getReal(), TOL);
        assertEquals(0.0, cm.get(1, 1).getImaginary(), TOL);

        // Test copy
        ComplexMatrix copy = cm.copy();
        assertNotNull(copy);
        assertEquals(cm.get(0, 0).getReal(), copy.get(0, 0).getReal(), TOL);

        // Test zero
        cm.zero();
        assertEquals(0.0, cm.get(0, 0).getReal(), TOL);

        // Test scale
        ComplexMatrix cm2 = new ComplexMatrix(2, 2);
        cm2.set(0, 0, new Complex(2.0, 1.0));
        cm2.scale(2.0);
        assertEquals(4.0, cm2.get(0, 0).getReal(), TOL);
        assertEquals(2.0, cm2.get(0, 0).getImaginary(), TOL);

        // Test sumRows
        ComplexMatrix cm3 = new ComplexMatrix(2, 2);
        cm3.set(0, 0, new Complex(1.0, 1.0));
        cm3.set(0, 1, new Complex(2.0, 2.0));
        ComplexMatrix rowSums = cm3.sumRows();
        assertEquals(2, rowSums.getNumRows());
        assertEquals(1, rowSums.getNumCols());

        // Test construction from real Matrix
        Matrix realMatrix = new Matrix(2, 2);
        realMatrix.set(0, 0, 3.0);
        ComplexMatrix fromReal = new ComplexMatrix(realMatrix);
        assertEquals(3.0, fromReal.get(0, 0).getReal(), TOL);
        assertEquals(0.0, fromReal.get(0, 0).getImaginary(), TOL);

        // Test determinant for simple 2x2 matrix
        ComplexMatrix det2x2 = new ComplexMatrix(2, 2);
        det2x2.set(0, 0, new Complex(1.0, 0.0));
        det2x2.set(0, 1, new Complex(2.0, 0.0));
        det2x2.set(1, 0, new Complex(3.0, 0.0));
        det2x2.set(1, 1, new Complex(4.0, 0.0));
        Complex det = det2x2.det();
        // det = 1*4 - 2*3 = -2
        assertEquals(-2.0, det.getReal(), TOL);
    }

    // ========== RowView ==========

    @Test
    void testRowView() {
        // Create a sparse matrix with some non-zero elements
        Matrix m = new Matrix(3, 4);
        m.set(0, 0, 1.0);
        m.set(0, 2, 3.0);
        m.set(1, 1, 2.0);
        m.set(2, 3, 4.0);

        // Get row view for row 0
        RowView row0 = m.getRowView(0);
        assertNotNull(row0);
        assertEquals(0, row0.getRowIndex());
        assertEquals(4, row0.getNumCols());
        assertTrue(row0.hasNonZeros());
        assertEquals(2, row0.getNonZeroCount());

        // Check non-zero values
        assertEquals(1.0, row0.get(0), TOL);
        assertEquals(0.0, row0.get(1), TOL);  // zero element
        assertEquals(3.0, row0.get(2), TOL);

        // Test getNonZeroCol and getNonZeroValue
        assertEquals(0, row0.getNonZeroCol(0));
        assertEquals(1.0, row0.getNonZeroValue(0), TOL);

        // Test dot product
        Matrix colVec = new Matrix(4, 1);
        colVec.set(0, 0, 1.0);
        colVec.set(1, 0, 1.0);
        colVec.set(2, 0, 1.0);
        colVec.set(3, 0, 1.0);
        double dotResult = row0.dotProduct(colVec);
        assertEquals(4.0, dotResult, TOL);  // 1*1 + 0*1 + 3*1 + 0*1 = 4

        // Test empty row
        RowView emptyRow = m.getRowView(1);
        emptyRow = m.getRowView(1);  // row 1 has only one non-zero
        assertTrue(emptyRow.hasNonZeros());
    }

    // ========== ColumnView ==========

    @Test
    void testColumnView() {
        // Create a sparse matrix with some non-zero elements
        Matrix m = new Matrix(4, 3);
        m.set(0, 0, 1.0);
        m.set(2, 0, 3.0);
        m.set(1, 1, 2.0);
        m.set(3, 2, 4.0);

        // Get column view for column 0
        ColumnView col0 = m.getColumnView(0);
        assertNotNull(col0);
        assertEquals(0, col0.getColumnIndex());
        assertEquals(4, col0.getNumRows());
        assertTrue(col0.hasNonZeros());
        assertEquals(2, col0.getNonZeroCount());

        // Check values
        assertEquals(1.0, col0.get(0), TOL);
        assertEquals(0.0, col0.get(1), TOL);  // zero element
        assertEquals(3.0, col0.get(2), TOL);

        // Test getNonZeroRow and getNonZeroValue
        assertEquals(0, col0.getNonZeroRow(0));
        assertEquals(1.0, col0.getNonZeroValue(0), TOL);

        // Test column with single non-zero
        ColumnView col1 = m.getColumnView(1);
        assertTrue(col1.hasNonZeros());
        assertEquals(1, col1.getNonZeroCount());
    }

    // ========== DenseMatrixImpl ==========

    @Test
    void testDenseMatrixImpl() {
        // DenseMatrixImpl is package-private, but we're in the same package so we can access it

        // Test constructor with dimensions
        DenseMatrixImpl dense = new DenseMatrixImpl(3, 3);
        assertNotNull(dense);
        assertEquals(3, dense.getNumRows());
        assertEquals(3, dense.getNumCols());

        // Test set/get
        dense.set(0, 0, 1.0);
        dense.set(1, 1, 2.0);
        dense.set(2, 2, 3.0);
        assertEquals(1.0, dense.get(0, 0), TOL);
        assertEquals(2.0, dense.get(1, 1), TOL);
        assertEquals(3.0, dense.get(2, 2), TOL);
        assertEquals(0.0, dense.get(0, 1), TOL);  // check zero

        // Test copy functionality
        BaseMatrix denseCopy = dense.copy();
        assertNotNull(denseCopy);
        assertEquals(dense.get(0, 0), denseCopy.get(0, 0), TOL);

        // Test toString
        String str = dense.toString();
        assertNotNull(str);
        assertTrue(str.contains("DenseMatrix"));

        // Test equals
        DenseMatrixImpl other = new DenseMatrixImpl(3, 3);
        other.set(0, 0, 1.0);
        other.set(1, 1, 2.0);
        other.set(2, 2, 3.0);
        assertTrue(dense.equals(other));

        // Test not equal (different value)
        other.set(0, 0, 99.0);
        assertFalse(dense.equals(other));

        // Test not equal (different size)
        DenseMatrixImpl different = new DenseMatrixImpl(2, 2);
        assertFalse(dense.equals(different));

        // Test constructor with nzLength (ignored for dense)
        DenseMatrixImpl withNz = new DenseMatrixImpl(2, 2, 100);
        assertEquals(2, withNz.getNumRows());
    }
}
