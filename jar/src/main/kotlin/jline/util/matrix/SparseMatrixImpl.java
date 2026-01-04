/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import org.ejml.data.DMatrixSparseCSC;

/**
 * Concrete implementation of SparseMatrix for delegation pattern.
 *
 * <p>This class provides the minimal implementation required for abstract methods
 * in the SparseMatrix base class. It serves as a concrete instantiable class that
 * implements the factory method pattern, allowing SparseMatrix to create instances
 * of itself through the delegation pattern.</p>
 *
 * <p>This implementation is primarily used internally by the Matrix class for
 * its delegation pattern, where Matrix can delegate operations to either a
 * SparseMatrixImpl or DenseMatrix depending on the matrix characteristics.</p>
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */
class SparseMatrixImpl extends SparseMatrix {

    /**
     * Creates a sparse matrix with specified dimensions and initial capacity for non-zero elements.
     *
     * @param numRows     the number of rows
     * @param numCols     the number of columns
     * @param arrayLength the initial capacity for non-zero elements
     */
    public SparseMatrixImpl(int numRows, int numCols, int arrayLength) {
        super(numRows, numCols, arrayLength);
    }

    /**
     * Creates a sparse matrix with specified dimensions and default capacity.
     *
     * @param numRows the number of rows
     * @param numCols the number of columns
     */
    public SparseMatrixImpl(int numRows, int numCols) {
        super(numRows, numCols);
    }

    /**
     * Creates a sparse matrix by wrapping an existing EJML sparse matrix.
     *
     * @param matrix the EJML sparse matrix to wrap
     */
    public SparseMatrixImpl(DMatrixSparseCSC matrix) {
        super(matrix);
    }

    /**
     * Copy constructor that creates a new sparse matrix from an existing one.
     *
     * @param matrix the sparse matrix to copy
     */
    public SparseMatrixImpl(SparseMatrix matrix) {
        super(matrix);
    }

    /**
     * Creates a deep copy of this sparse matrix.
     *
     * @return a new SparseMatrixImpl that is a copy of this instance
     */
    @Override
    public SparseMatrix copy() {
        return new SparseMatrixImpl(this);
    }

    /**
     * Factory method for creating new instances of this sparse matrix implementation.
     * This method is required by the abstract base class for the delegation pattern.
     *
     * @param rows     the number of rows for the new matrix
     * @param cols     the number of columns for the new matrix
     * @param nzLength the initial capacity for non-zero elements
     * @return a new SparseMatrixImpl instance
     */
    @Override
    protected SparseMatrix createNewInstance(int rows, int cols, int nzLength) {
        return new SparseMatrixImpl(rows, cols, nzLength);
    }

    /**
     * Determines whether this sparse matrix is equal to another object.
     * Two sparse matrices are considered equal if they have the same dimensions
     * and all corresponding elements are equal.
     *
     * @param obj the object to compare with this matrix
     * @return true if the matrices are equal, false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;

        SparseMatrixImpl other = (SparseMatrixImpl) obj;

        if (getNumRows() != other.getNumRows() || getNumCols() != other.getNumCols()) {
            return false;
        }

        // Compare all elements
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                if (Double.compare(get(i, j), other.get(i, j)) != 0) {
                    return false;
                }
            }
        }

        return true;
    }
}