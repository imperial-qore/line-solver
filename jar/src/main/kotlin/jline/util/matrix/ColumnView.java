/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import org.ejml.data.DMatrixSparseCSC;

/**
 * A lightweight view into a column of a sparse matrix that doesn't copy data.
 * This provides efficient access to matrix column elements without the overhead
 * of creating a new Matrix object and copying all the data.
 *
 * <p>This is particularly useful for large sparse matrices where you need to
 * access column elements frequently, such as in iterative algorithms.</p>
 *
 * <p>Example usage:</p>
 * <pre>
 * Matrix bigMatrix = ...; // 1M x 1000 sparse matrix
 * ColumnView col = bigMatrix.getColumnView(5);
 *
 * // Efficient access - only processes non-zero elements
 * for (int i = 0; i < col.getNonZeroCount(); i++) {
 *     int row = col.getNonZeroRow(i);
 *     double value = col.getNonZeroValue(i);
 *     // process (row, value) pair
 * }
 *
 * // Or get any element (returns 0.0 for non-stored elements)
 * double val = col.get(100);
 * </pre>
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */
public class ColumnView {
    private final DMatrixSparseCSC matrix;
    private final int columnIndex;
    private final int startIdx;
    private final int endIdx;
    private final int numRows;

    /**
     * Creates a view into the specified column of the matrix.
     *
     * @param matrix      the source sparse matrix
     * @param columnIndex the column index to view
     */
    ColumnView(DMatrixSparseCSC matrix, int columnIndex) {
        if (columnIndex < 0 || columnIndex >= matrix.numCols) {
            throw new IndexOutOfBoundsException("Column index " + columnIndex + " out of bounds [0, " + matrix.numCols + ")");
        }

        this.matrix = matrix;
        this.columnIndex = columnIndex;
        this.startIdx = matrix.col_idx[columnIndex];
        this.endIdx = matrix.col_idx[columnIndex + 1];
        this.numRows = matrix.numRows;
    }

    /**
     * Returns the value at the specified row in this column.
     * This method searches through the non-zero elements, so it's O(nnz_in_column).
     * For frequent access to multiple elements, consider using the non-zero iterators.
     *
     * @param row the row index
     * @return the value at (row, column), or 0.0 if not stored
     */
    public double get(int row) {
        if (row < 0 || row >= numRows) {
            throw new IndexOutOfBoundsException("Row index " + row + " out of bounds [0, " + numRows + ")");
        }

        // Binary search through the non-zero rows in this column
        for (int i = startIdx; i < endIdx; i++) {
            if (matrix.nz_rows[i] == row) {
                return matrix.nz_values[i];
            }
            if (matrix.nz_rows[i] > row) {
                break; // rows are sorted, so we won't find it
            }
        }
        return 0.0;
    }

    /**
     * Returns the column index that this view represents.
     *
     * @return the column index
     */
    public int getColumnIndex() {
        return columnIndex;
    }

    /**
     * Returns the number of non-zero elements in this column.
     *
     * @return the count of non-zero elements
     */
    public int getNonZeroCount() {
        return endIdx - startIdx;
    }

    /**
     * Returns the row index of the i-th non-zero element in this column.
     *
     * @param i the index into the non-zero elements (0 to getNonZeroCount()-1)
     * @return the row index
     */
    public int getNonZeroRow(int i) {
        if (i < 0 || i >= getNonZeroCount()) {
            throw new IndexOutOfBoundsException("Non-zero index " + i + " out of bounds [0, " + getNonZeroCount() + ")");
        }
        return matrix.nz_rows[startIdx + i];
    }

    /**
     * Returns the value of the i-th non-zero element in this column.
     *
     * @param i the index into the non-zero elements (0 to getNonZeroCount()-1)
     * @return the value
     */
    public double getNonZeroValue(int i) {
        if (i < 0 || i >= getNonZeroCount()) {
            throw new IndexOutOfBoundsException("Non-zero index " + i + " out of bounds [0, " + getNonZeroCount() + ")");
        }
        return matrix.nz_values[startIdx + i];
    }

    /**
     * Returns the total number of rows in the column (including zeros).
     *
     * @return the number of rows
     */
    public int getNumRows() {
        return numRows;
    }

    /**
     * Checks if this column has any non-zero elements.
     *
     * @return true if the column contains at least one non-zero element
     */
    public boolean hasNonZeros() {
        return getNonZeroCount() > 0;
    }
}