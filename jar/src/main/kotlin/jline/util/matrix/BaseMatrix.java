/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import org.ejml.data.DMatrix;

import java.io.Serializable;

/**
 * Common base class for matrix implementations, providing a unified interface
 * for both dense and sparse matrix operations.
 *
 * <p>This abstract class defines the common contract that all matrix implementations
 * must support, enabling the delegation pattern in Matrix.java to work with either
 * dense or sparse underlying representations.</p>
 *
 * <p>The class provides:</p>
 * <ul>
 *   <li>Basic matrix operations (get, set, dimensions)</li>
 *   <li>Common utility methods for matrix manipulation</li>
 *   <li>Abstract methods that must be implemented by concrete subclasses</li>
 *   <li>Factory methods for creating new instances</li>
 * </ul>
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */
public abstract class BaseMatrix implements Serializable {

    /**
     * Replaces each value in the matrix with its absolute value, in-place.
     */
    public abstract void absEq();

    /**
     * Performs the operation: output = alpha*A + beta*B.
     *
     * @param alpha  scalar coefficient for matrix A
     * @param A      first input matrix
     * @param beta   scalar coefficient for matrix B
     * @param B      second input matrix
     * @param output the result matrix
     */
    protected abstract DMatrix addMatrices(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output);

    /**
     * Checks if the current matrix contains any non-zero elements.
     *
     * @return {@code true} if at least one non-zero element exists; {@code false} otherwise.
     */
    public abstract boolean any();

    /**
     * Applies a transformation to matrix elements based on conditions.
     * Replaces elements matching the condition with the target value.
     */
    protected abstract void applyConditionalTransform(double source, double target, double tol, String operation);

    // Additional mathematical operations
    protected abstract void changeSign();

    /**
     * Creates and returns a deep copy of this matrix.
     *
     * @return a new matrix that is a copy of this instance
     */
    public abstract BaseMatrix copy();

    // Common sparse matrix operations that may need default implementations for dense matrices

    /**
     * Adds a scalar value to each element in the specified column.
     *
     * @param col Column index to update.
     * @param a   Scalar value to add to each element in the column.
     */
    public abstract void colIncrease(int col, double a);

    /**
     * Concatenates matrices column-wise.
     *
     * @param left   the left matrix
     * @param right  the right matrix
     * @param output the result matrix
     */
    protected abstract void concatColumnsInPlace(DMatrix left, DMatrix right, DMatrix output);

    // Abstract methods for sparse-specific operations that dense matrices must implement

    /**
     * Concatenates matrices row-wise.
     *
     * @param top    the top matrix
     * @param bottom the bottom matrix
     * @param output the result matrix
     */
    protected abstract void concatRowsInPlace(DMatrix top, DMatrix bottom, DMatrix output);

    /**
     * Counts occurrences of each unique value in each row.
     */
    public abstract BaseMatrix countEachRow(double val);

    /**
     * Creates a new instance of the concrete matrix implementation.
     * This factory method allows the base class to create instances of the correct subtype.
     *
     * @param rows     the number of rows for the new matrix
     * @param cols     the number of columns for the new matrix
     * @param nzLength the initial capacity for non-zero elements (may be ignored by dense matrices)
     * @return a new instance of the concrete matrix type
     */
    protected abstract BaseMatrix createNewInstance(int rows, int cols, int nzLength);

    protected abstract double determinant();

    protected abstract void divideInPlace(double scalar);

    protected abstract void divideMatrix(double scalar, DMatrix output);

    protected abstract void divideRowsByArray(double[] diag, int offset);

    protected abstract void divideScalarByMatrix(double scalar, DMatrix output);

    protected abstract double elementMax();

    protected abstract double elementMaxAbs();

    protected abstract double elementMin();

    // Element-wise operations
    protected abstract DMatrix elementMult(DMatrix A, DMatrix B, DMatrix output);

    // Data access methods - these may need different implementations for dense vs sparse

    protected abstract double elementSum();

    /**
     * Indicates whether some other object is "equal to" this matrix.
     *
     * @param obj the reference object with which to compare
     * @return {@code true} if this object is the same as the obj argument; {@code false} otherwise
     */
    public abstract boolean equals(Object obj);

    // Extraction methods
    protected abstract void extractMatrix(DMatrix src, int srcX0, int srcX1, int srcY0, int srcY1, DMatrix dst, int dstY0, int dstX0);

    // Advanced matrix operations that both dense and sparse matrices should support

    /**
     * Fills the matrix with the specified value.
     *
     * @param value the value to fill the matrix with
     */
    protected abstract void fillMatrix(double value);

    /**
     * Returns a matrix containing indices of all non-negative elements.
     */
    public abstract BaseMatrix findNonNegative();

    /**
     * Returns the value at the specified matrix position.
     *
     * @param row the row index
     * @param col the column index
     * @return the value at position (row, col)
     */
    public abstract double get(int row, int col);

    protected abstract int getColumnIndex(int col);

    protected abstract int[] getColumnIndicesArray();

    // Matrix operation methods that may need sparse matrix parameters but should work for both types

    /**
     * Returns the underlying data structure.
     * For sparse matrices, returns DMatrix.
     * For dense matrices, may return a wrapped DMatrixRMaj or converted sparse representation.
     *
     * @return the underlying data structure
     */
    protected abstract Object getData();

    /**
     * Sets the underlying data structure.
     * Implementation depends on the concrete matrix type.
     *
     * @param newData the new data structure
     */
    protected abstract void setData(Object newData);

    /**
     * Returns the number of stored non-zero elements.
     * This may be different from actual non-zero count for dense matrices.
     *
     * @return number of stored non-zero elements
     */
    public abstract int getNonZeroLength();

    protected abstract void setNonZeroLength(int length);

    protected abstract int getNonZeroRow(int index);

    // Internal data access methods - these will have very different implementations

    /**
     * Returns array of row indices of non-zero entries.
     * For dense matrices, this may compute row indices dynamically.
     * For sparse matrices, this returns the compressed row index array.
     *
     * @return array of row indices
     */
    public abstract int[] getNonZeroRows();

    protected abstract double getNonZeroValue(int index);

    // Sparse-specific internal methods that dense matrices must provide equivalents for

    /**
     * Returns array of non-zero values.
     * For dense matrices, this may return all values or only actual non-zeros.
     * For sparse matrices, this returns the compressed storage array.
     *
     * @return array of non-zero values
     */
    public abstract double[] getNonZeroValues();

    /**
     * Returns the number of non-zero elements in the matrix.
     * For dense matrices, this counts actual non-zero values.
     * For sparse matrices, this returns the stored non-zero count.
     *
     * @return the number of non-zero elements
     */
    public abstract int getNonZeros();

    /**
     * Returns the number of columns in the matrix.
     *
     * @return the number of columns
     */
    public abstract int getNumCols();

    /**
     * Returns the number of non-zero elements in the matrix.
     * This is an alias for {@link #getNonZeros()}.
     *
     * @return the number of non-zero elements
     */
    public int getNumNonZeros() {
        return getNonZeros();
    }

    /**
     * Returns the number of rows in the matrix.
     *
     * @return the number of rows
     */
    public abstract int getNumRows();

    /**
     * Checks if the matrix contains duplicate values.
     *
     * @return true if duplicates exist
     */
    public abstract boolean hasDuplicates();

    /**
     * Checks if the matrix contains any finite (non-infinite, non-NaN) values.
     *
     * @return true if at least one finite value exists
     */
    public abstract boolean hasFinite();

    /**
     * Checks if the matrix contains any infinite values.
     *
     * @return true if at least one infinite value exists
     */
    public abstract boolean hasInfinite();

    /**
     * Checks if more than one finite value exists.
     *
     * @return true if more than one finite value exists
     */
    public abstract boolean hasMultipleFinite();

    /**
     * Checks if the matrix contains any NaN values.
     *
     * @return true if at least one NaN exists
     */
    public abstract boolean hasNaN();

    /**
     * Performs matrix multiplication with sparse matrices.
     *
     * @param B      the right matrix
     * @param output the result matrix
     */
    protected abstract DMatrix multMatrix(DMatrix B, DMatrix output);

    /**
     * Removes zero-valued elements from the matrix with default tolerance.
     * For dense matrices, this may be a no-op or convert to sparse representation.
     */
    protected abstract void removeZeros();

    /**
     * Removes zero-valued elements from the matrix with specified tolerance.
     *
     * @param tolerance the tolerance for considering values as zero
     */
    protected abstract void removeZerosWithTol(double tolerance);

    /**
     * Reshapes the matrix to the specified dimensions.
     *
     * @param numRows the new number of rows
     * @param numCols the new number of columns
     */
    public abstract void reshape(int numRows, int numCols);

    /**
     * Scales the matrix in-place by the specified factor.
     *
     * @param alpha the scaling factor
     */
    protected abstract void scaleInPlace(double alpha);

    protected abstract void scaleMatrix(double scalar, DMatrix output);

    /**
     * Sets the value at the specified matrix position.
     *
     * @param row   the row index
     * @param col   the column index
     * @param value the value to set
     */
    public abstract void set(int row, int col, double value);

    protected abstract void setColumnIndex(int col, int value);

    protected abstract void setNonZeroRow(int index, int row);

    protected abstract void setNonZeroValue(int index, double value);

    /**
     * Reduce the maximum number of columns by setting the internal column count.
     *
     * @param newmax the new maximum number of columns
     */
    public abstract void shrinkNumCols(int newmax);

    /**
     * Reduce the maximum number of rows by setting the internal row count.
     *
     * @param newmax the new maximum number of rows
     */
    public abstract void shrinkNumRows(int newmax);

    protected abstract DMatrix sumColsRaw();

    protected abstract DMatrix sumRowsRaw();

    /**
     * Returns a string representation of the matrix.
     * Elements are formatted to 4 decimal places in a grid layout.
     *
     * @return a formatted string representation of the matrix
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                sb.append(String.format("%8.4f ", get(i, j)));
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    /**
     * Computes the transpose of this matrix.
     *
     * @param output the transposed matrix
     */
    protected abstract void transposeMatrix(DMatrix output);

    /**
     * Sets all elements in the matrix to zero.
     */
    public abstract void zero();
}