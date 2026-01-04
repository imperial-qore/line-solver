/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.sparse.csc.CommonOps_DSCC;

/**
 * A matrix class for handling complex-valued matrices using separate real and imaginary components.
 *
 * <p>ComplexMatrix provides a comprehensive implementation for complex matrix operations by storing
 * the real and imaginary parts as separate Matrix objects. This approach allows for efficient
 * manipulation of complex numbers in matrix computations while leveraging the existing Matrix
 * infrastructure.</p>
 *
 * <p>Key features:</p>
 * <ul>
 *   <li>Separate storage of real and imaginary components</li>
 *   <li>Support for all standard complex matrix operations</li>
 *   <li>Integration with Apache Commons Math for complex number operations</li>
 *   <li>Efficient determinant calculation using LU decomposition</li>
 *   <li>Row and column extraction operations</li>
 * </ul>
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */
public class ComplexMatrix {
    /**
     * The real component of the complex matrix
     */
    public Matrix real;
    /**
     * The imaginary component of the complex matrix
     */
    public Matrix im;

    /**
     * Creates a new complex matrix with the specified dimensions.
     * All elements are initialized to zero (both real and imaginary parts).
     *
     * @param i the number of rows
     * @param j the number of columns
     */
    public ComplexMatrix(int i, int j) {
        this.real = new Matrix(i, j);
        this.im = new Matrix(i, j);
    }

    /**
     * Creates a complex matrix from separate real and imaginary component matrices.
     *
     * @param real the real component matrix
     * @param im   the imaginary component matrix
     * @throws AssertionError if the dimensions of real and imaginary matrices don't match
     */
    public ComplexMatrix(Matrix real, Matrix im) {
        assert (real.getNumRows() == im.getNumRows() && real.getNumCols() == im.getNumCols());
        this.real = real;
        this.im = im;
    }

    /**
     * Creates a complex matrix from EJML sparse matrices representing real and imaginary components.
     *
     * @param matrix_real the real component as DMatrixSparseCSC
     * @param matrix_im   the imaginary component as DMatrixSparseCSC
     */
    public ComplexMatrix(DMatrixSparseCSC matrix_real, DMatrixSparseCSC matrix_im) {
        this.real = new Matrix(matrix_real);
        this.im = new Matrix(matrix_im);

    }

    /**
     * Creates a complex matrix from a real matrix with zero imaginary component.
     *
     * @param real the real component matrix (imaginary part will be set to zero)
     */
    public ComplexMatrix(Matrix real) {
        this.real = real;
        this.im = new Matrix(real.getNumRows(), real.getNumCols());
        this.im.zero();
    }

    /**
     * Concatenates two complex matrices vertically (row-wise).
     *
     * @param top    the upper complex matrix
     * @param bottom the lower complex matrix
     * @param out    the output matrix, or null to create a new one
     * @return a complex matrix with the concatenated rows
     */
    public static ComplexMatrix concatRows(ComplexMatrix top, ComplexMatrix bottom, ComplexMatrix out) {
        if (out == null) {
            return new ComplexMatrix(Matrix.concatRows(top.real, bottom.real, null), Matrix.concatRows(top.im, bottom.im, null));
        } else {
            Matrix.concatRows(top.real, bottom.real, out.real);
            Matrix.concatRows(top.im, bottom.im, out.im);
            return out;
        }
    }

    /**
     * Extracts a range of rows from a complex matrix.
     *
     * @param A    the source complex matrix
     * @param row0 the first row to extract (inclusive)
     * @param row1 the last row to extract (exclusive)
     * @param out  the output matrix, or null to create a new one
     * @return a complex matrix containing the extracted rows
     */
    public static ComplexMatrix extractRows(ComplexMatrix A, int row0, int row1, ComplexMatrix out) {
        if (out == null) {
            return new ComplexMatrix(CommonOps_DSCC.extractRows((DMatrixSparseCSC) A.real.getData(), row0, row1, null), CommonOps_DSCC.extractRows((DMatrixSparseCSC) A.im.getData(), row0, row1, null));
        } else {
            CommonOps_DSCC.extractRows((DMatrixSparseCSC) A.real.getData(), row0, row1, (DMatrixSparseCSC) out.real.getData());
            CommonOps_DSCC.extractRows((DMatrixSparseCSC) A.im.getData(), row0, row1, (DMatrixSparseCSC) out.im.getData());
            return out;
        }
    }

    /**
     * Creates a deep copy of this complex matrix.
     *
     * @return a new ComplexMatrix that is a copy of this instance
     */
    public ComplexMatrix copy() {
        return new ComplexMatrix(this.real.copy(), this.im.copy());
    }

    /**
     * Computes the determinant of this complex matrix using LU decomposition.
     *
     * @return the complex determinant of the matrix
     */
    public Complex det() {
        FieldMatrix<Complex> a = MatrixUtils.createFieldMatrix(ComplexField.getInstance(), this.getNumRows(), this.getNumCols());
        for (int i = 0; i < this.getNumRows(); i++) {
            for (int j = 0; j < this.getNumCols(); j++) {
                a.setEntry(i, j, this.get(i, j));
            }
        }
        FieldLUDecomposition<Complex> LU = new FieldLUDecomposition<>(a);
        Complex det = (Complex) LU.getDeterminant();
        return det;
    }

    /**
     * Gets the complex element at the specified position.
     *
     * @param i the row index
     * @param j the column index
     * @return the complex value at the specified position
     */
    public Complex get(int i, int j) {
        return new Complex(this.real.get(i, j), this.im.get(i, j));
    }

    /**
     * Gets the complex element at the specified linear index.
     *
     * @param idx the linear index
     * @return the complex value at the specified index
     */
    public Complex get(int idx) {
        return new Complex(this.real.get(idx), this.im.get(idx));
    }

    /**
     * Returns the number of columns in this complex matrix.
     *
     * @return the number of columns
     */
    public int getNumCols() {
        return this.real.getNumCols();
    }

    /**
     * Returns the total number of elements in this complex matrix.
     *
     * @return the total number of elements (rows * columns)
     */
    public int getNumElements() {
        return this.real.getNumElements();
    }

    /**
     * Returns the number of rows in this complex matrix.
     *
     * @return the number of rows
     */
    public int getNumRows() {
        return this.real.getNumRows();
    }

    /**
     * Checks if this complex matrix is empty (both real and imaginary parts have no elements).
     *
     * @return true if the matrix is empty, false otherwise
     */
    public boolean isEmpty() {
        return this.real.isEmpty() && this.im.isEmpty();
    }

    /**
     * Scales this complex matrix by a real scalar value.
     * Both real and imaginary components are multiplied by the scalar.
     *
     * @param a the scalar value to multiply by
     */
    public void scale(double a) {
        this.real.scaleEq(a);
        this.im.scaleEq(a);
    }

    /**
     * Sets the element at the specified position to a complex value.
     *
     * @param i   the row index
     * @param j   the column index
     * @param val the complex value to set
     */
    public void set(int i, int j, Complex val) {
        this.real.set(i, j, val.getReal());
        this.im.set(i, j, val.getImaginary());
    }

    /**
     * Sets the element at the specified position to an integer value (real part only).
     * Special handling for Integer.MAX_VALUE and Integer.MIN_VALUE as positive and negative infinity.
     *
     * @param row the row index
     * @param col the column index
     * @param val the integer value to set
     */
    public void set(int row, int col, int val) {
        if (val == Integer.MAX_VALUE) {
            this.real.set(row, col, Inf);
        } else if (val == Integer.MIN_VALUE) {
            this.real.set(row, col, NegInf);
        } else {
            this.real.set(row, col, (double) val);
            if (val == 0)
                this.real.remove(row, col); //Remove to ensure that getNonZeroElement not contains the value with 0
        }
    }

    /**
     * Sets the element at the specified position to a real value (imaginary part becomes zero).
     *
     * @param i   the row index
     * @param j   the column index
     * @param val the real value to set
     */
    public void set(int i, int j, double val) {
        this.real.set(i, j, val);
    }

    /**
     * Sets the element at the specified linear index to a complex value.
     *
     * @param idx the linear index
     * @param val the complex value to set
     */
    public void set(int idx, Complex val) {
        this.real.set(idx, val.getReal());
        this.im.set(idx, val.getImaginary());
    }

    /**
     * Computes the sum of each row, returning a complex matrix with one column.
     *
     * @return a complex matrix containing the row sums
     */
    public ComplexMatrix sumRows() {
        return new ComplexMatrix(this.real.sumRows(), this.im.sumRows());
    }

    /**
     * Sets all elements of this complex matrix to zero (both real and imaginary parts).
     */
    public void zero() {
        this.real.zero();
        this.im.zero();
    }


}
