/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.ArrayFieldVector;
import org.apache.commons.math3.linear.FieldDecompositionSolver;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.FieldVector;
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

    /**
     * Creates a complex identity matrix of the specified size.
     *
     * @param n the size of the identity matrix
     * @return an n x n complex identity matrix
     */
    public static ComplexMatrix eye(int n) {
        return new ComplexMatrix(Matrix.eye(n), new Matrix(n, n));
    }

    /**
     * Creates a complex zero matrix of the specified dimensions.
     *
     * @param rows the number of rows
     * @param cols the number of columns
     * @return a rows x cols complex zero matrix
     */
    public static ComplexMatrix zeros(int rows, int cols) {
        return new ComplexMatrix(new Matrix(rows, cols), new Matrix(rows, cols));
    }

    /**
     * Multiplies this complex matrix by another complex matrix.
     * (A+iB)(C+iD) = (AC-BD) + i(AD+BC)
     *
     * @param other the matrix to multiply by
     * @return the product of the two complex matrices
     */
    public ComplexMatrix mult(ComplexMatrix other) {
        Matrix ac = this.real.mult(other.real);
        Matrix bd = this.im.mult(other.im);
        Matrix ad = this.real.mult(other.im);
        Matrix bc = this.im.mult(other.real);
        return new ComplexMatrix(ac.add(-1.0, bd), ad.add(1.0, bc));
    }

    /**
     * Adds another complex matrix to this one, returning a new matrix.
     *
     * @param other the matrix to add
     * @return a new ComplexMatrix that is the sum of this and other
     */
    public ComplexMatrix add(ComplexMatrix other) {
        return new ComplexMatrix(this.real.add(1.0, other.real), this.im.add(1.0, other.im));
    }

    /**
     * Subtracts another complex matrix from this one, returning a new matrix.
     *
     * @param other the matrix to subtract
     * @return a new ComplexMatrix that is this minus other
     */
    public ComplexMatrix sub(ComplexMatrix other) {
        return new ComplexMatrix(this.real.add(-1.0, other.real), this.im.add(-1.0, other.im));
    }

    /**
     * Returns the transpose of this complex matrix.
     *
     * @return the transpose
     */
    public ComplexMatrix transpose() {
        return new ComplexMatrix(this.real.transpose(), this.im.transpose());
    }

    /**
     * Returns the conjugate transpose (Hermitian transpose) of this complex matrix.
     * For (A+iB), the conjugate transpose is (A^T - iB^T).
     *
     * @return the conjugate transpose
     */
    public ComplexMatrix conjugateTranspose() {
        Matrix negIm = this.im.scale(-1.0);
        return new ComplexMatrix(this.real.transpose(), negIm.transpose());
    }

    /**
     * Scales this complex matrix by a complex scalar, returning a new matrix.
     * (a+ib)(C+iD) = (aC-bD) + i(aD+bC)
     *
     * @param z the complex scalar
     * @return a new ComplexMatrix scaled by z
     */
    public ComplexMatrix scaleComplex(Complex z) {
        double a = z.getReal();
        double b = z.getImaginary();
        Matrix newReal = this.real.scale(a).add(-1.0, this.im.scale(b));
        Matrix newIm = this.real.scale(b).add(1.0, this.im.scale(a));
        return new ComplexMatrix(newReal, newIm);
    }

    /**
     * Solves the complex linear system A*x = b. For square systems, uses LU decomposition.
     * For overdetermined systems (more rows than columns), uses the real-embedding approach
     * with SVD-based pseudo-inverse, which correctly handles rank-deficient systems.
     *
     * <p>The complex system (A_r + i*A_i)(x_r + i*x_i) = (b_r + i*b_i) is converted to
     * the equivalent real system:
     * <pre>
     *   [A_r, -A_i] [x_r]   [b_r]
     *   [A_i,  A_r] [x_i] = [b_i]
     * </pre>
     * This 2m x 2n real system is solved via the real Matrix.leftMatrixDivide (SVD-based).</p>
     *
     * @param b the right-hand side complex column vector or matrix
     * @return the solution vector x
     */
    public ComplexMatrix leftMatrixDivide(ComplexMatrix b) {
        int m = this.getNumRows();
        int n = this.getNumCols();
        int nrhs = b.getNumCols();

        if (m == n) {
            // Square system: use LU decomposition directly
            return solveLU(this, b);
        }

        // Convert to real-embedded system
        // [A_r, -A_i; A_i, A_r] * [x_r; x_i] = [b_r; b_i]
        Matrix Ar = this.real;
        Matrix Ai = this.im;

        // Build 2m x 2n real matrix
        Matrix realSys = new Matrix(2 * m, 2 * n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                double rval = Ar.get(i, j);
                double ival = Ai.get(i, j);
                if (rval != 0.0) {
                    realSys.set(i, j, rval);
                    realSys.set(m + i, n + j, rval);
                }
                if (ival != 0.0) {
                    realSys.set(i, n + j, -ival);
                    realSys.set(m + i, j, ival);
                }
            }
        }

        ComplexMatrix result = new ComplexMatrix(n, nrhs);
        for (int col = 0; col < nrhs; col++) {
            // Build 2m x 1 real RHS
            Matrix realRHS = new Matrix(2 * m, 1);
            for (int i = 0; i < m; i++) {
                Complex bval = b.get(i, col);
                if (bval.getReal() != 0.0) realRHS.set(i, 0, bval.getReal());
                if (bval.getImaginary() != 0.0) realRHS.set(m + i, 0, bval.getImaginary());
            }

            // Solve using real SVD-based least squares
            Matrix realSol = realSys.leftMatrixDivide(realRHS);

            // Extract complex solution
            for (int j = 0; j < n; j++) {
                double xr = realSol.get(j, 0);
                double xi = realSol.get(n + j, 0);
                if (xr != 0.0 || xi != 0.0) {
                    result.set(j, col, new Complex(xr, xi));
                }
            }
        }
        return result;
    }

    /**
     * Solves a square complex linear system A*x = b using LU decomposition
     * via Apache Commons Math3 FieldMatrix.
     *
     * @param A the coefficient matrix (must be square)
     * @param b the right-hand side (column vector or matrix)
     * @return the solution x
     */
    private static ComplexMatrix solveLU(ComplexMatrix A, ComplexMatrix b) {
        int n = A.getNumRows();
        int nrhs = b.getNumCols();

        // Convert to Apache Commons Math3 FieldMatrix
        FieldMatrix<Complex> fieldA = MatrixUtils.createFieldMatrix(ComplexField.getInstance(), n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fieldA.setEntry(i, j, A.get(i, j));
            }
        }

        FieldLUDecomposition<Complex> lu = new FieldLUDecomposition<Complex>(fieldA);
        FieldDecompositionSolver<Complex> solver = lu.getSolver();

        ComplexMatrix result = new ComplexMatrix(nrhs > 0 ? n : 0, nrhs);
        for (int col = 0; col < nrhs; col++) {
            Complex[] bCol = new Complex[n];
            for (int i = 0; i < n; i++) {
                bCol[i] = b.get(i, col);
            }
            FieldVector<Complex> bVec = new ArrayFieldVector<Complex>(ComplexField.getInstance(), bCol);
            FieldVector<Complex> xVec = solver.solve(bVec);
            for (int i = 0; i < n; i++) {
                result.set(i, col, xVec.getEntry(i));
            }
        }
        return result;
    }

    /**
     * Converts this ComplexMatrix to an Apache Commons Math3 FieldMatrix.
     *
     * @return the equivalent FieldMatrix
     */
    public FieldMatrix<Complex> toFieldMatrix() {
        int m = getNumRows();
        int n = getNumCols();
        FieldMatrix<Complex> fm = MatrixUtils.createFieldMatrix(ComplexField.getInstance(), m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                fm.setEntry(i, j, get(i, j));
            }
        }
        return fm;
    }

}
