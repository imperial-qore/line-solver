/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import org.ejml.data.DMatrixSparseCSC;

/**
 * A lightweight view into a row of a sparse matrix that doesn't copy data.
 * This provides efficient access to matrix row elements without the overhead
 * of creating a new Matrix object and copying all the data.
 *
 * <p>This is particularly useful for large sparse matrices where you need to
 * access row elements frequently, such as in iterative algorithms.</p>
 *
 * <p>Example usage:</p>
 * <pre>
 * Matrix bigMatrix = ...; // 1000 x 1M sparse matrix
 * RowView row = bigMatrix.getRowView(5);
 *
 * // Efficient access - only processes non-zero elements
 * for (int i = 0; i < row.getNonZeroCount(); i++) {
 *     int col = row.getNonZeroCol(i);
 *     double value = row.getNonZeroValue(i);
 *     // process (col, value) pair
 * }
 *
 * // Or get any element (returns 0.0 for non-stored elements)
 * double val = row.get(100);
 * </pre>
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */
public class RowView {
    private final DMatrixSparseCSC matrix;
    private final int rowIndex;
    private final int numCols;
    private final int[] nonZeroCols;
    private final double[] nonZeroValues;
    private final int nonZeroCount;

    /**
     * Creates a view into the specified row of the matrix.
     *
     * @param matrix   the source sparse matrix
     * @param rowIndex the row index to view
     */
    RowView(DMatrixSparseCSC matrix, int rowIndex) {
        if (rowIndex < 0 || rowIndex >= matrix.numRows) {
            throw new IndexOutOfBoundsException("Row index " + rowIndex + " out of bounds [0, " + matrix.numRows + ")");
        }

        this.matrix = matrix;
        this.rowIndex = rowIndex;
        this.numCols = matrix.numCols;

        // Extract non-zero elements for this row from the CSC format
        // This requires scanning through all columns to find elements in this row
        int[] tempCols = new int[matrix.numCols]; // Worst case: all columns have this row
        double[] tempValues = new double[matrix.numCols];
        int count = 0;

        // Scan through each column to find non-zero elements in this row
        for (int col = 0; col < matrix.numCols; col++) {
            int colStart = matrix.col_idx[col];
            int colEnd = matrix.col_idx[col + 1];

            // Binary search for rowIndex in this column
            for (int idx = colStart; idx < colEnd; idx++) {
                if (matrix.nz_rows[idx] == rowIndex) {
                    tempCols[count] = col;
                    tempValues[count] = matrix.nz_values[idx];
                    count++;
                    break;
                } else if (matrix.nz_rows[idx] > rowIndex) {
                    break; // rows are sorted, so we won't find it
                }
            }
        }

        // Copy to appropriately sized arrays
        this.nonZeroCount = count;
        this.nonZeroCols = new int[count];
        this.nonZeroValues = new double[count];
        System.arraycopy(tempCols, 0, this.nonZeroCols, 0, count);
        System.arraycopy(tempValues, 0, this.nonZeroValues, 0, count);
    }

    /**
     * Computes the dot product of this row with a column vector.
     * This is optimized to only iterate through the non-zero elements of the row.
     *
     * @param columnVector the column vector to multiply with (must have getNumRows() == this.getNumCols())
     * @return the scalar result of the dot product
     * @throws IllegalArgumentException if dimensions don't match
     */
    public double dotProduct(Matrix columnVector) {
        if (columnVector.getNumCols() != 1) {
            throw new IllegalArgumentException("Second argument must be a column vector (1 column)");
        }
        if (columnVector.getNumRows() != numCols) {
            throw new IllegalArgumentException("Dimension mismatch: row has " + numCols +
                    " columns but column vector has " + columnVector.getNumRows() + " rows");
        }

        double result = 0.0;

        // Iterate through non-zero elements in the row
        for (int i = 0; i < nonZeroCount; i++) {
            int col = nonZeroCols[i];
            double rowValue = nonZeroValues[i];
            double colValue = columnVector.get(col, 0);
            result += rowValue * colValue;
        }

        return result;
    }

    /**
     * Returns the value at the specified column in this row.
     * This method searches through the non-zero elements, so it's O(nnz_in_row).
     * For frequent access to multiple elements, consider using the non-zero iterators.
     *
     * @param col the column index
     * @return the value at (row, col), or 0.0 if not stored
     */
    public double get(int col) {
        if (col < 0 || col >= numCols) {
            throw new IndexOutOfBoundsException("Column index " + col + " out of bounds [0, " + numCols + ")");
        }

        // Linear search through the non-zero columns in this row
        for (int i = 0; i < nonZeroCount; i++) {
            if (nonZeroCols[i] == col) {
                return nonZeroValues[i];
            }
            if (nonZeroCols[i] > col) {
                break; // columns are sorted, so we won't find it
            }
        }
        return 0.0;
    }

    /**
     * Returns the column index of the i-th non-zero element in this row.
     *
     * @param i the index into the non-zero elements (0 to getNonZeroCount()-1)
     * @return the column index
     */
    public int getNonZeroCol(int i) {
        if (i < 0 || i >= nonZeroCount) {
            throw new IndexOutOfBoundsException("Non-zero index " + i + " out of bounds [0, " + nonZeroCount + ")");
        }
        return nonZeroCols[i];
    }

    /**
     * Returns the number of non-zero elements in this row.
     *
     * @return the count of non-zero elements
     */
    public int getNonZeroCount() {
        return nonZeroCount;
    }

    /**
     * Returns the value of the i-th non-zero element in this row.
     *
     * @param i the index into the non-zero elements (0 to getNonZeroCount()-1)
     * @return the value
     */
    public double getNonZeroValue(int i) {
        if (i < 0 || i >= nonZeroCount) {
            throw new IndexOutOfBoundsException("Non-zero index " + i + " out of bounds [0, " + nonZeroCount + ")");
        }
        return nonZeroValues[i];
    }

    /**
     * Returns the total number of columns in the row (including zeros).
     *
     * @return the number of columns
     */
    public int getNumCols() {
        return numCols;
    }

    /**
     * Returns the row index that this view represents.
     *
     * @return the row index
     */
    public int getRowIndex() {
        return rowIndex;
    }

    /**
     * Checks if this row has any non-zero elements.
     *
     * @return true if the row contains at least one non-zero element
     */
    public boolean hasNonZeros() {
        return nonZeroCount > 0;
    }
}