/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import org.ejml.data.DMatrixRMaj;

/**
 * Concrete implementation of DenseMatrix for delegation pattern.
 *
 * <p>This class provides the minimal implementation required for abstract methods
 * in the DenseMatrix base class. It serves as a concrete instantiable class that
 * implements the factory method pattern, allowing DenseMatrix to create instances
 * of itself through the delegation pattern.</p>
 *
 * <p>This implementation is primarily used internally by the Matrix class for
 * its delegation pattern, where Matrix can delegate operations to either a
 * SparseMatrixImpl or DenseMatrixImpl depending on the matrix characteristics.</p>
 *
 * <p>Dense matrices store all elements explicitly (including zeros), making them
 * suitable for matrices with a high percentage of non-zero elements or when
 * dense linear algebra operations are required. They offer O(1) element access
 * but use more memory compared to sparse representations.</p>
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */
class DenseMatrixImpl extends DenseMatrix {

    /**
     * Creates a dense matrix with specified dimensions.
     * All elements are initialized to zero.
     *
     * @param numRows the number of rows
     * @param numCols the number of columns
     */
    public DenseMatrixImpl(int numRows, int numCols) {
        super(numRows, numCols);
    }

    /**
     * Creates a dense matrix with specified dimensions and initial capacity.
     * The nzLength parameter is ignored for dense matrices but maintained for API consistency.
     *
     * @param numRows  the number of rows
     * @param numCols  the number of columns
     * @param nzLength ignored for dense matrices (maintained for compatibility)
     */
    public DenseMatrixImpl(int numRows, int numCols, int nzLength) {
        super(numRows, numCols);
    }

    /**
     * Creates a dense matrix by wrapping an existing EJML dense matrix.
     *
     * @param matrix the EJML dense matrix to wrap
     */
    public DenseMatrixImpl(DMatrixRMaj matrix) {
        super(matrix);
    }

    /**
     * Copy constructor that creates a new dense matrix from an existing one.
     *
     * @param matrix the dense matrix to copy
     */
    public DenseMatrixImpl(DenseMatrix matrix) {
        super(matrix);
    }

    /**
     * Creates a deep copy of this dense matrix.
     *
     * @return a new DenseMatrixImpl that is a copy of this instance
     */
    @Override
    public BaseMatrix copy() {
        return new DenseMatrixImpl(this);
    }

    /**
     * Factory method for creating new instances of this dense matrix implementation.
     * This method is required by the abstract base class for the delegation pattern.
     *
     * @param rows     the number of rows for the new matrix
     * @param cols     the number of columns for the new matrix
     * @param nzLength ignored for dense matrices (maintained for compatibility)
     * @return a new DenseMatrixImpl instance
     */
    @Override
    protected BaseMatrix createNewInstance(int rows, int cols, int nzLength) {
        return new DenseMatrixImpl(rows, cols, nzLength);
    }

    /**
     * Determines whether this dense matrix is equal to another object.
     * Two dense matrices are considered equal if their string representations are the same.
     *
     * @param obj the object to compare with this matrix
     * @return true if the matrices are equal, false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;

        DenseMatrixImpl other = (DenseMatrixImpl) obj;

        if (getNumRows() != other.getNumRows() || getNumCols() != other.getNumCols()) {
            return false;
        }

        // Compare all elements with tolerance
        double tolerance = 1e-10;
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                double diff = Math.abs(get(i, j) - other.get(i, j));
                if (diff > tolerance) {
                    return false;
                }
            }
        }

        return true;
    }


    /**
     * Returns a string representation of this dense matrix.
     * For large matrices, only a subset of elements are shown.
     *
     * @return a string representation of the matrix
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("DenseMatrix[").append(getNumRows()).append("x").append(getNumCols()).append("]\n");

        int maxRowsToShow = Math.min(5, getNumRows());
        int maxColsToShow = Math.min(8, getNumCols());

        for (int i = 0; i < maxRowsToShow; i++) {
            sb.append("[");
            for (int j = 0; j < maxColsToShow; j++) {
                if (j > 0) sb.append(", ");
                sb.append(String.format("%8.4f", get(i, j)));
            }
            if (getNumCols() > maxColsToShow) {
                sb.append(", ...");
            }
            sb.append("]\n");
        }

        if (getNumRows() > maxRowsToShow) {
            sb.append("...\n");
        }

        return sb.toString();
    }
}