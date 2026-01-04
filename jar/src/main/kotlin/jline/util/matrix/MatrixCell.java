/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * An ordered collection of Matrix objects that provides indexing and manipulation operations.
 *
 * <p>MatrixCell serves as a container for multiple Matrix instances, allowing them to be
 * stored, retrieved, and operated upon as a group. This is particularly useful in scenarios
 * where multiple related matrices need to be managed together, such as in matrix decompositions
 * or iterative algorithms.</p>
 *
 * <p>Key features:</p>
 * <ul>
 *   <li>Integer-indexed access to matrices</li>
 *   <li>Support for matrix operations across all contained matrices</li>
 *   <li>Automatic cloning for deep copy operations</li>
 *   <li>Null-safe operations and cleanup utilities</li>
 * </ul>
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */
public class MatrixCell implements Serializable {
    private final Map<Integer, Matrix> data;

    /**
     * Creates an empty MatrixCell with no initial capacity constraints.
     */
    public MatrixCell() {
        data = new HashMap<>();
    }

    /**
     * Creates a deep copy of another MatrixCell.
     * All matrices in the source MatrixCell are cloned to ensure independence.
     *
     * @param x the MatrixCell to copy
     */
    public MatrixCell(MatrixCell x) {
        data = new HashMap<>(x.size());
        for (int i = 0; i < x.size(); i++) {
            data.put(i, x.get(i).copy());
        }
    }

    /**
     * Creates a MatrixCell from an array of matrices.
     * The matrices are stored with indices 0, 1, 2, ... corresponding to their array positions.
     *
     * @param x the array of matrices to store
     */
    public MatrixCell(Matrix[] x) {
        data = new HashMap<>(x.length);
        for (int i = 0; i < x.length; i++) {
            data.put(i, x[i]);
        }
    }

    /**
     * Creates an empty MatrixCell with a specified initial capacity.
     *
     * @param K the initial capacity for the underlying storage
     */
    public MatrixCell(int K) {
        data = new HashMap<>(K);
    }

    /**
     * Creates a MatrixCell containing exactly two matrices.
     * This is a convenience constructor for the common case of matrix pairs.
     *
     * @param D0 the first matrix (stored at index 0)
     * @param D1 the second matrix (stored at index 1)
     */
    public MatrixCell(Matrix D0, Matrix D1) {
        data = new HashMap<>(2);
        data.put(0, D0);
        data.put(1, D1);
    }

    /**
     * Computes the sum of elements at a specific position across all matrices in the cell.
     * This operation is useful for aggregating corresponding elements from multiple matrices.
     *
     * @param row the row index
     * @param col the column index
     * @return the sum of elements at position (row, col) across all matrices
     */
    public double cellsum(int row, int col) {
        double sum = 0.0;
        for (int i = 0; i < data.size(); i++) {
            sum += data.get(i).get(row, col);
        }
        return sum;
    }

    /**
     * Computes the element-wise sum of all matrices in this cell.
     * Each matrix must have the same dimensions. The matrices are summed in the order
     * of their integer keys (0 to N-1). The result is accumulated in-place for efficiency.
     *
     * @return A Matrix representing the element-wise sum of all matrices in the cell
     * @throws RuntimeException if the cell is empty or if any matrices have incompatible dimensions
     */
    public Matrix cellsum() {
        if (data.isEmpty()) {
            throw new RuntimeException("MatrixCell is empty");
        }

        Matrix result = data.get(0).copy();

        for (int i = 1; i < data.size(); i++) {
            Matrix m = data.get(i);
            if (m.getNumRows() != result.getNumRows() || m.getNumCols() != result.getNumCols()) {
                throw new RuntimeException("Matrices have incompatible sizes for add operation");
            }

            result = result.add(1, m);
        }

        return result;
    }

    /**
     * Retrieves the matrix stored at the specified index.
     *
     * @param i the index of the matrix to retrieve
     * @return the matrix at the specified index, or null if no matrix exists at that index
     */
    public Matrix get(int i) {
        return data.get(i);
    }

    /**
     * Checks whether this MatrixCell contains any matrices.
     *
     * @return true if the MatrixCell is empty, false otherwise
     */
    public boolean isEmpty() {
        return data.isEmpty();
    }

    /**
     * Prints all matrices in the collection to the standard output.
     * Each matrix is printed using its own print() method.
     */
    public void print() {
        for (int i = 0; i < data.size(); i++) {
            data.get(i).print();
        }
    }

    /**
     * Removes the matrix at the specified index.
     *
     * @param i the index of the matrix to remove
     */
    public void remove(int i) {
        data.remove(i);
    }

    /**
     * Removes all null matrices from the collection.
     * This method iterates backwards through the indices to safely remove entries
     * without affecting the iteration process.
     */
    public void removeNull() {
        for (int i = data.size() - 1; i >= 0; i--) {
            if (data.get(i) == null) {
                data.remove(i);
            }
        }
    }

    /**
     * Stores a matrix at the specified index.
     * If a matrix already exists at that index, it is replaced.
     *
     * @param i the index at which to store the matrix
     * @param D the matrix to store
     * @return the previous matrix at this index, or null if none existed
     */
    public Matrix set(int i, Matrix D) {
        return data.put(i, D);
    }

    /**
     * Returns the number of matrices currently stored in this MatrixCell.
     *
     * @return the number of matrices in the collection
     */
    public int size() {
        return data.size();
    }

    /**
     * Creates a copy of the internal mapping from indices to matrices.
     * This method provides access to the underlying data structure while maintaining encapsulation.
     *
     * @return a new Map containing copies of the index-to-matrix mappings
     */
    public Map<Integer, Matrix> toMap() {
        Map<Integer, Matrix> map = new HashMap<>();
        for (int i = 0; i < data.size(); i++) {
            map.put(i, data.get(i));
        }
        return map;
    }
}
