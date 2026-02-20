/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;

import jline.io.Ret;
import jline.GlobalConstants;
import jline.util.*;
import jline.util.graph.UndirectedGraph;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.ejml.interfaces.decomposition.QRSparseDecomposition;
import org.ejml.ops.DConvertMatrixStruct;
import org.ejml.simple.SimpleMatrix;
import org.ejml.sparse.csc.decomposition.qr.QrLeftLookingDecomposition_DSCC;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

import static jline.io.InputOutputKt.mfilename;
import static jline.util.Utils.isInf;
import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.line_warning;

/**
 * A sparse matrix data structure supporting linear algebra functions similar to those available in
 * MATLAB.
 * 
 * TABLE OF CONTENTS:
 * 1. CONSTRUCTORS AND INITIALIZATION
 * 2. CORE ELEMENT ACCESS AND MODIFICATION
 * 3. BASIC ARITHMETIC OPERATIONS
 * 4. ELEMENT-WISE OPERATIONS
 * 5. MATRIX TRANSFORMATIONS
 * 6. LINEAR ALGEBRA OPERATIONS
 * 7. MATRIX SLICING AND EXTRACTION
 * 8. MATRIX CONCATENATION AND ASSEMBLY
 * 9. STATISTICAL OPERATIONS
 * 10. MATRIX PROPERTIES AND COMPARISONS
 * 11. FINDING AND FILTERING
 * 12. MATRIX MANIPULATION
 * 13. SPECIAL MATRIX OPERATIONS
 * 14. I/O OPERATIONS
 * 15. FACTORY METHODS
 * 16. UTILITY METHODS
 */
public final class Matrix implements Serializable {

    /**
     * The delegate matrix that handles the underlying operations
     */
    private final BaseMatrix delegate;

    // ================================================================================
    // SECTION 1: CONSTRUCTORS AND INITIALIZATION
    // ================================================================================
    // Methods for creating and initializing matrices
    
    /**
     * Creates a sparse matrix with the specified dimensions and capacity.
     *
     * @param numRows     number of rows in the matrix
     * @param numCols     number of columns in the matrix
     * @param arrayLength initial capacity for non-zero elements
     */
    public Matrix(int numRows, int numCols, int arrayLength) {
        this.delegate = new SparseMatrixImpl(numRows, numCols, arrayLength);
    }

    /**
     * Creates a sparse matrix with the specified dimensions.
     *
     * @param numRows number of rows in the matrix
     * @param numCols number of columns in the matrix
     */
    public Matrix(int numRows, int numCols) {
        this.delegate = new SparseMatrixImpl(numRows, numCols);
    }

    /**
     * Creates a copy of the specified matrix.
     *
     * @param matrix the matrix to copy
     * @throws UnsupportedOperationException if the matrix type is not supported
     */
    public Matrix(Matrix matrix) {
        // Create a new sparse matrix implementation by copying the delegate's data
        if (matrix.delegate instanceof SparseMatrix) {
            this.delegate = new SparseMatrixImpl((SparseMatrix) matrix.delegate);
        } else if (matrix.delegate instanceof DenseMatrix) {
            this.delegate = new DenseMatrixImpl((DenseMatrix) matrix.delegate);
        } else {
            throw new UnsupportedOperationException("Copying from non-sparse matrix delegate not yet supported");
        }
    }

    /**
     * Creates a matrix from an EJML DMatrix (sparse).
     *
     * @param matrix the EJML sparse matrix to wrap
     */
    public Matrix(DMatrix matrix) {
        this.delegate = new SparseMatrixImpl((DMatrixSparseCSC) matrix);
    }

    /**
     * Creates a matrix from an EJML DMatrixRMaj (dense).
     *
     * @param matrix the EJML dense matrix to wrap
     */
    public Matrix(DMatrixRMaj matrix) {
        this.delegate = new DenseMatrixImpl(matrix);
    }

    /**
     * Creates a matrix from a 2D double array.
     * Uses EJML's triplet format for efficient single-pass batch insertion.
     *
     * @param arrays 2D array containing the matrix elements
     */
    public Matrix(double[][] arrays) {
        int numRows = arrays.length;
        if (numRows == 0) {
            this.delegate = new SparseMatrixImpl(new DMatrixSparseCSC(0, 0, 0));
            return;
        }
        int numCols = arrays[0].length;

        // Single-pass: estimate initial capacity and let triplet grow if needed
        // Use 10% density estimate as initial capacity, minimum of 16
        int estimatedNz = Math.max(16, (numRows * numCols) / 10);
        DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(numRows, numCols, estimatedNz);

        for (int i = 0; i < numRows; i++) {
            double[] row = arrays[i];
            int rowLen = row.length;
            for (int j = 0; j < rowLen; j++) {
                double val = row[j];
                if (val != 0.0) {
                    triplet.addItem(i, j, val);
                }
            }
        }

        // Convert triplet to CSC format
        DMatrixSparseCSC csc = DConvertMatrixStruct.convert(triplet, (DMatrixSparseCSC) null);
        this.delegate = new SparseMatrixImpl(csc);
    }

    /**
     * Creates a matrix from a Kotlin Array<DoubleArray> type.
     * This constructor enables direct initialization from Kotlin's arrayOf(doubleArrayOf(...)) syntax.
     * Uses EJML's triplet format for efficient single-pass batch insertion.
     *
     * @param arrays Kotlin Array containing DoubleArray rows
     */
    public Matrix(Double[][] arrays) {
        int numRows = arrays.length;
        int numCols = arrays[0].length;

        // Single-pass: estimate initial capacity and let triplet grow if needed
        int estimatedNz = Math.max(16, (numRows * numCols) / 10);
        DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(numRows, numCols, estimatedNz);

        for (int i = 0; i < numRows; i++) {
            Double[] row = arrays[i];
            int rowLen = row.length;
            for (int j = 0; j < rowLen; j++) {
                double val = row[j];
                if (val != 0.0) {
                    triplet.addItem(i, j, val);
                }
            }
        }

        // Convert triplet to CSC format
        DMatrixSparseCSC csc = DConvertMatrixStruct.convert(triplet, (DMatrixSparseCSC) null);
        this.delegate = new SparseMatrixImpl(csc);
    }

    /**
     * Creates a column vector from an integer array.
     *
     * @param array integer array to convert to column vector
     */
    public Matrix(int[] array) {
        this.delegate = new SparseMatrixImpl(array.length, 1, array.length);
        for (int i = 0; i < array.length; i++) {
            this.set(i, 0, (double) array[i]);
        }
    }

    /**
     * Creates a column vector from a double array.
     *
     * @param array double array to convert to column vector
     */
    public Matrix(double[] array) {
        this.delegate = new SparseMatrixImpl(array.length, 1, array.length);
        for (int i = 0; i < array.length; i++) {
            this.set(i, 0, array[i]);
        }
    }

    // construct column vector
    public Matrix(List<Double> array) {
        this.delegate = new SparseMatrixImpl(array.size(), 1, array.size());
        for (int i = 0; i < array.size(); i++) this.set(i, 0, array.get(i));
    }

    public Matrix(ArrayList<List<Double>> arrays) {
        // in this case we concatenate the arrays as columns
        // Use maximum size among all arrays to determine row count (handles empty first array)
        int maxSize = 0;
        for (List<Double> arr : arrays) {
            if (arr.size() > maxSize) {
                maxSize = arr.size();
            }
        }
        this.delegate = new SparseMatrixImpl(maxSize, arrays.size(), maxSize * arrays.size());
        for (int j = 0; j < arrays.size(); j++)
            for (int i = 0; i < arrays.get(j).size(); i++) this.set(i, j, arrays.get(j).get(i));
    }

    public Matrix(SimpleMatrix matrix) {
        this.delegate = new SparseMatrixImpl(matrix.numRows(), matrix.numCols());
        for (int i = 0; i < matrix.numRows(); i++) {
            for (int j = 0; j < matrix.numCols(); j++) {
                this.set(i, j, matrix.get(i, j));
            }
        }
    }

    /**
     * Parse matrix from string in MATLAB or Python formats MATLAB format assumes commas between row
     * elements and semicolon between rows
     *
     * @param matrixString input string, e.g., "[10,4;5,9]" or "[[10,4],[4,9]]"
     */
    public Matrix(String matrixString) {
        if (matrixString.isEmpty()) {
            this.delegate = new SparseMatrixImpl(0, 0, 0);
            return;
        }

        // Remove spaces
        matrixString = matrixString.replaceAll("\\s", "");

        // Determine if the string is in NumPy (Python) style based on the presence of "[["
        boolean isNumPyStyle = matrixString.startsWith("[[");

        if (isNumPyStyle) {
            // NumPy style: [[1,2],[3,4]]
            matrixString =
                    matrixString.substring(1, matrixString.length() - 1); // Remove outer square brackets
            String[] rows = matrixString.split("\\],\\[");

            int numRows = rows.length;
            int numCols = rows[0].split(",").length;

            this.delegate = new SparseMatrixImpl(numRows, numCols);
            for (int i = 0; i < numRows; i++) {
                String row = rows[i].replaceAll("\\[|\\]", ""); // Remove inner square brackets
                String[] elements = row.split(",");
                for (int j = 0; j < elements.length; j++) {
                    this.set(i, j, Double.parseDouble(elements[j]));
                }
            }
        } else {
            // MATLAB style: [1,2;3,4]
            matrixString = matrixString.replaceAll("\\[|\\]", "");
            String[] rows = matrixString.split(";");

            int numRows = rows.length;
            int numCols = rows[0].split(",").length;

            this.delegate = new SparseMatrixImpl(numRows, numCols);
            for (int i = 0; i < numRows; i++) {
                String[] elements = rows[i].split(",");
                for (int j = 0; j < elements.length; j++) {
                    this.set(i, j, Double.parseDouble(elements[j]));
                }
            }
        }
    }

    // ================================================================================
    // SECTION 2: CORE ELEMENT ACCESS AND MODIFICATION
    // ================================================================================
    // Basic operations for accessing and modifying matrix elements
    
    /**
     * Returns all elements in a matrix except the first ones
     *
     * @param y - a matrix
     * @param k - a position of an element in the matrix
     * @return - A row vector with all elements of y but ypos
     */
    public static Matrix allbut(Matrix y, int k) {
        int total = y.getNumElements();
        if (k < 0 || k >= total) {
            return new Matrix(1, total);
        }

        Matrix ret = new Matrix(1, total - 1);
        int limit = Math.min(k, total);
        for (int j = 0; j < limit; j++) {
            ret.set(j, y.get(j));
        }
        for (int j = k + 1; j < total; j++) {
            ret.set(j - 1, y.get(j));
        }

        return ret;
    }

    /**
     * Computes the element-wise sum of a column vector and a row vector,
     * producing a matrix where each entry (i, j) is the sum of colVector[i] and rowVector[j].
     * This is equivalent to broadcasting the column vector along columns and the row vector along rows.
     *
     * @param colVector A column vector of size (m x 1)
     * @param rowVector A row vector of size (1 x n)
     * @return A matrix of size (m x n) where each element is colVector[i] + rowVector[j]
     */
    public static Matrix broadcastColPlusRow(Matrix colVector, Matrix rowVector) {
        int m = colVector.getNumRows();
        int n = rowVector.getNumCols();

        Matrix result = new Matrix(m, n, m * n);

        for (int i = 0; i < m; i++) {
            double colVal = colVector.get(i);
            for (int j = 0; j < n; j++) {
                double sum = colVal + rowVector.get(j);
                if (sum != 0.0) result.set(i, j, sum);  // preserve sparsity
            }
        }

        return result;
    }

    /**
     * Cartesian product of two matrices. It replicates elements of the first input matrix and pairs
     * them with each row of the second input matrix.
     *
     * @param matrixA first input matrix
     * @param matrixB second input matrix
     */
    public static Matrix cartesian(Matrix matrixA, Matrix matrixB) {
        if (matrixA.isEmpty()) return matrixB;
        if (matrixB.isEmpty()) return matrixA;

        int n1 = matrixA.getNumRows();
        int m1 = matrixA.getNumCols();
        int n2 = matrixB.getNumRows();
        int m2 = matrixB.getNumCols();

        int totalRows = n1 * n2;
        int totalCols = m1 + m2;

        Matrix result = new Matrix(totalRows, totalCols);

        // Fill left block: matrixA repeated n2 times
        for (int i = 0; i < n2; i++) {
            for (int r = 0; r < n1; r++) {
                int destRow = i * n1 + r;
                for (int c = 0; c < m1; c++) {
                    double val = matrixA.get(r, c);
                    if (val != 0.0) result.set(destRow, c, val);  // skip zero to preserve sparsity
                }
            }
        }

        // Fill right block: each row of matrixB repeated n1 times
        for (int i = 0; i < n2; i++) {
            for (int r = 0; r < n1; r++) {
                int destRow = i * n1 + r;
                for (int c = 0; c < m2; c++) {
                    double val = matrixB.get(i, c);
                    if (val != 0.0) result.set(destRow, m1 + c, val);  // preserve sparsity
                }
            }
        }

        return result;
    }

    /**
     * Concatenates a collection of matrices stored in a map into a single row vector.
     * Each matrix is flattened in column-major order and concatenated sequentially based on the map's value order.
     * Only non-zero values are inserted to preserve sparse structure.
     *
     * @param cellArray A map from integer indices to Matrix objects to be concatenated
     * @return A single-row Matrix containing all elements from the input matrices in sequence
     */
    public static Matrix cell2mat(Map<Integer, Matrix> cellArray) {
        int totalElements = Matrix.getColIndexSum(cellArray);
        Matrix newMatrix = new Matrix(1, totalElements);

        int currentColIndex = 0;
        for (Matrix matrix : cellArray.values()) {
            DMatrixSparseCSC data = (DMatrixSparseCSC) matrix.getData();
            int numCols = data.getNumCols();
            int[] colIndices = matrix.getColumnIndicesArray();
            double[] values = matrix.getNonZeroValues();
            int[] rowIndices = matrix.getNonZeroRows();

            for (int col = 0; col < numCols; col++) {
                int idx0 = colIndices[col];
                int idx1 = colIndices[col + 1];
                for (int i = idx0; i < idx1; i++) {
                    newMatrix.set(0, currentColIndex + rowIndices[i], values[i]);
                }
                currentColIndex += data.numRows;
            }
        }

        return newMatrix;
    }

    /**
     * Computes the element-wise sum of all matrices stored in a map.
     * Each matrix must have the same dimensions. The matrices are summed in the order
     * of their integer keys (0 to N-1). The result is accumulated in-place for efficiency.
     *
     * @param cellArray A map from integer indices to Matrix objects to be summed
     * @return A Matrix representing the element-wise sum of all matrices in the map
     * @throws RuntimeException if the map is empty or if any matrices have incompatible dimensions
     */
    public static Matrix cellsum(Map<Integer, Matrix> cellArray) {
        // Create a MatrixCell and delegate to its cellsum method
        MatrixCell cell = new MatrixCell();
        for (Map.Entry<Integer, Matrix> entry : cellArray.entrySet()) {
            cell.set(entry.getKey(), entry.getValue());
        }
        return cell.cellsum();
    }

    // ================================================================================
    // ================================================================================
    // SECTION 16: UTILITY METHODS
    // ================================================================================
    // Helper methods and internal utilities
    // ================================================================================
    // Helper and utility methods for various matrix operations
    
    /**
     * Converts a column matrix (of size m x 1) to a dense double array of length m.
     *
     * @param columnMatrix A Matrix with exactly one column
     * @return A double array containing the values from the column matrix
     * @throws IllegalArgumentException if the input matrix does not have exactly one column
     */
    public static double[] columnMatrixToDoubleArray(Matrix columnMatrix) {
        if (columnMatrix.getNumCols() != 1) {
            throw new IllegalArgumentException("Input matrix must be a column matrix with only one column.");
        }

        DMatrixSparseCSC data = (DMatrixSparseCSC) columnMatrix.getData();
        double[] resultArray = new double[data.getNumRows()];

        int[] colIndices = columnMatrix.getColumnIndicesArray();
        int[] rowIndices = columnMatrix.getNonZeroRows();
        double[] values = columnMatrix.getNonZeroValues();

        int start = colIndices[0];
        int end = colIndices[1];
        for (int i = start; i < end; i++) {
            resultArray[rowIndices[i]] = values[i];
        }

        return resultArray;
    }

    /**
     * Compares two matrices element-wise using a specified comparison operator.
     * The matrices must have the same dimensions. Supported operators are:
     * "eq", "equal", "lt", "lessthan", "lte", "lessthanequal",
     * "gt", "greater", "gte", "greaterthanequal".
     *
     * @param A  The first matrix
     * @param B  The second matrix
     * @param op The comparison operator as a string
     * @return true if the comparison holds for all elements, false otherwise
     */
    public static boolean compare(Matrix A, Matrix B, String op) {
        if (A.getNumRows() != B.getNumRows() || A.getNumCols() != B.getNumCols()) {
            throw new IllegalArgumentException("Matrix dimensions must match for comparison.");
        }

        int rows = A.getNumRows();
        int cols = A.getNumCols();

        switch (op.toLowerCase()) {
            case "eq":
            case "equal":
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        if (A.get(i, j) != B.get(i, j)) return false;
                    }
                }
                return true;
            case "lt":
            case "lessthan":
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        if (A.get(i, j) >= B.get(i, j)) return false;
                    }
                }
                return true;
            case "lte":
            case "lessthanequal":
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        if (A.get(i, j) > B.get(i, j)) return false;
                    }
                }
                return true;
            case "gt":
            case "greater":
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        if (A.get(i, j) <= B.get(i, j)) return false;
                    }
                }
                return true;
            case "gte":
            case "greaterthanequal":
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        if (A.get(i, j) < B.get(i, j)) return false;
                    }
                }
                return true;
            default:
                throw new IllegalArgumentException("Invalid comparison operator: " + op);
        }
    }

    /**
     * Concatenates two matrices horizontally (column-wise).
     * If either matrix is empty, returns the non-empty matrix.
     * If an output matrix is provided, writes the result into it; otherwise, returns a newly allocated matrix.
     *
     * @param left  The matrix to appear on the left
     * @param right The matrix to appear on the right
     * @param out   (Optional) Output matrix to store the result; can be null
     * @return A matrix representing [left | right]
     */
    public static Matrix concatColumns(Matrix left, Matrix right, Matrix out) {
        if (right.getNumRows() == 0 && right.getNumCols() == 0) {
            return left;
        }

        if (left.getNumRows() == 0 && left.getNumCols() == 0) {
            return right;
        }

        if (out == null) {
            return new Matrix(SparseMatrix.concatColumnsStatic(left.getData(), right.getData(), null));
        } else {
            SparseMatrix.concatColumnsInPlaceStatic(left.getData(), right.getData(), out.getData());
            return out;
        }
    }

    /**
     * Concatenates two matrices vertically (row-wise).
     * If either matrix is empty, returns the non-empty matrix.
     * If an output matrix is provided, writes the result into it; otherwise, returns a newly allocated matrix.
     *
     * @param top    The matrix to appear on top
     * @param bottom The matrix to appear below
     * @param out    (Optional) Output matrix to store the result; can be null
     * @return A matrix representing [top; bottom]
     */
    public static Matrix concatRows(Matrix top, Matrix bottom, Matrix out) {
        if (top.getNumCols() == 0 && top.getNumRows() == 0) {
            return bottom;
        }

        if (bottom.getNumCols() == 0 && bottom.getNumRows() == 0) {
            return top;
        }

        if (out == null) {
            return new Matrix(SparseMatrix.concatRowsStatic(top.getData(), bottom.getData(), null));
        } else {
            SparseMatrix.concatRowsInPlaceStatic(top.getData(), bottom.getData(), out.getData());
            return out;
        }
    }

    /**
     * Creates a new empty matrix with the same shape and internal structure as the given matrix.
     * The contents are uninitialized and all values are implicitly zero (sparse).
     *
     * @param B The matrix to copy the structure from
     * @return A new matrix with the same dimensions and sparsity structure as {@code B}
     */
    public static Matrix createLike(Matrix B) {
        return new Matrix((DMatrixSparseCSC) B.getData().createLike());
    }

    public static Matrix decorate(Matrix inspace1, Matrix inspace2) {
        if (inspace2 == null || inspace2.getNumCols() == 0 || inspace2.getNumRows() == 0) {
            Matrix C = inspace1;
            inspace1 = new Matrix(0, 0);
            for (int c = 0; c < C.length(); c++) {
                Matrix tempMatrix = new Matrix(1, 1);
                tempMatrix.set(0, 0, C.get(c));
                inspace1 = decorate(inspace1, tempMatrix);
            }
            return inspace1;
        }
        if (inspace1.isEmpty()) {
            inspace1 = inspace2.copy();
            return inspace1;
        }
        if (inspace2.isEmpty()) {
            return inspace1;
        }
        int n1 = inspace1.getNumRows();
        int m1 = inspace1.getNumCols();
        int n2 = inspace2.getNumRows();
        int m2 = inspace2.getNumCols();
        inspace1 = inspace1.repmat(n2, 1);
        List<Integer> curStates = new ArrayList<>();
        for (int i = 0; i < n1; i++) {
            curStates.add(i);
        }
        for (int s = 0; s < n2; s++) {
            Matrix inspace2_rowS = inspace2.getRow(s);
            Matrix value = inspace2_rowS.repmat(curStates.size(), 1);
            int valueRowIdx = 0;
            inspace1.expandMatrix(inspace1.getNumRows(), m1 + m2, inspace1.getNumNonZeros());
            for (int curState : curStates) {
                int valueColIdx = 0;
                for (int col = m1; col < m1 + m2; col++) {
                    inspace1.set(curState, col, value.get(valueRowIdx, valueColIdx));
                    valueColIdx += 1;
                }
                valueRowIdx += 1;
            }
            for (int i = 0; i < curStates.size(); i++) {
                int updatedValue = curStates.get(i) + n1;
                curStates.set(i, updatedValue);
            }
        }
        return inspace1;
    }

    // ================================================================================
    // ================================================================================
    // SECTION 15: FACTORY METHODS
    // ================================================================================
    // Static methods for creating special matrices: zeros, ones, eye, diag
    // ================================================================================
    // Static methods for creating special matrices: zeros, ones, eye, diag, etc.
    
    /**
     * Creates a square diagonal matrix from the given values.
     * The input values are placed on the main diagonal; all off-diagonal entries are zero.
     *
     * @param values The values to place on the diagonal
     * @return A square sparse matrix with {@code values[i]} at position (i, i)
     */
    public static Matrix diag(double... values) {
        return new Matrix(SparseMatrix.diagStatic(values));
    }

    /**
     * Creates a square diagonal matrix from a given array of values.
     * This is a low-level wrapper that constructs a diagonal matrix using EJML internals.
     * Equivalent to placing {@code values[i]} at position (i, i) for all i.
     *
     * @param values An array of values to place on the main diagonal
     * @return A square sparse matrix with {@code values[i]} at position (i, i)
     */
    public static Matrix diagMatrix(double[] values) {
        return new Matrix(SparseMatrix.diagWithRetStatic(null, values, 0, values.length));
    }

    /**
     * Creates a square diagonal matrix from the elements of a column or row vector matrix.
     * Each element in the input matrix is placed on the main diagonal of the resulting matrix.
     *
     * @param values A matrix interpreted as a vector (either row or column) containing diagonal values
     * @return A square sparse matrix with the input values along the main diagonal
     */
    public static Matrix diagMatrix(Matrix values) {
        return Matrix.diagMatrix(null, values.toArray1D(), 0, values.length());
    }

    /**
     * Creates or fills a square diagonal matrix with the specified values.
     * If the input matrix {@code A} is null, a new diagonal matrix is created.
     * Otherwise, the provided matrix {@code A} is filled with the diagonal values.
     *
     * @param A      An optional matrix to be filled with the diagonal values; can be null
     * @param values The array of values to place on the main diagonal
     * @param offset The starting index in {@code values} to read from
     * @param length The number of values to place on the diagonal
     * @return A matrix with {@code values[offset + i]} placed at position (i, i)
     */
    public static Matrix diagMatrix(Matrix A, double[] values, int offset, int length) {
        if (A == null) {
            return new Matrix(SparseMatrix.diagWithRetStatic(null, values, offset, length));
        } else {
            return new Matrix(SparseMatrix.diagWithMatrixStatic((DMatrixSparseCSC) A.getData(), values, offset, length));
        }
    }

    /**
     * Returns the smallest non-zero positive element in the given matrix.
     * If no such element exists, returns 0.
     *
     * @param matrix The input matrix
     * @return The minimum strictly positive element, or 0 if all elements are zero or negative
     */
    public static double elementMinNonZero(Matrix matrix) {
        double[] values = matrix.getNonZeroValues();
        double min = Double.MAX_VALUE;

        for (double v : values) {
            if (v > 0 && v < min) {
                min = v;
            }
        }

        return min == Double.MAX_VALUE ? 0 : min;
    }

    /**
     * Extracts a rectangular submatrix from a source matrix and stores it in a destination matrix.
     * Equivalent to slicing rows [srcY0:srcY1) and columns [srcX0:srcX1) from {@code src}.
     *
     * @param src   The source matrix
     * @param srcX0 Starting column (inclusive)
     * @param srcX1 Ending column (exclusive)
     * @param srcY0 Starting row (inclusive)
     * @param srcY1 Ending row (exclusive)
     * @param dst   The destination matrix to store the result
     * @param dstY0 Starting row index in {@code dst}
     * @param dstX0 Starting column index in {@code dst}
     */
    public static void extract(
            Matrix src, int srcX0, int srcX1, int srcY0, int srcY1, Matrix dst, int dstY0, int dstX0) {
        SparseMatrix.extractMatrixStatic(src.getData(), srcX0, srcX1, srcY0, srcY1, dst.getData(), dstY0, dstX0);
    }

    /**
     * Extracts a rectangular submatrix from the given source matrix.
     * Equivalent to {@code src[srcY0:srcY1, srcX0:srcX1]}.
     *
     * @param src   The source matrix
     * @param srcX0 Starting column (inclusive)
     * @param srcX1 Ending column (exclusive)
     * @param srcY0 Starting row (inclusive)
     * @param srcY1 Ending row (exclusive)
     * @return A new matrix containing the extracted submatrix
     */
    public static Matrix extract(Matrix src, int srcX0, int srcX1, int srcY0, int srcY1) {
        if (srcX0 >= srcX1 || srcY0 >= srcY1) {
            return new Matrix(0, 0);
        }
        Matrix out = new Matrix(srcX1 - srcX0, srcY1 - srcY0);
        extract(src, srcX0, srcX1, srcY0, srcY1, out, 0, 0);
        return out;
    }

    /**
     * Extracts a single column from the given matrix.
     *
     * @param A      The input matrix
     * @param column The index of the column to extract
     * @param out    (Optional) Matrix to store the output; if null, a new matrix is created
     * @return A matrix containing the extracted column
     */
    public static Matrix extractColumn(Matrix A, int column, Matrix out) {
        if (out == null) {
            return new Matrix(SparseMatrix.extractColumnStatic(A.getData(), column, null));
        } else {
            SparseMatrix.extractColumnInPlaceStatic(A.getData(), column, out.getData());
            return out;
        }
    }

    /**
     * Extracts a range of columns from the matrix [col0:col1).
     *
     * @param A    The source matrix
     * @param col0 Starting column index (inclusive)
     * @param col1 Ending column index (exclusive)
     * @return A new matrix containing the specified columns
     */
    public static Matrix extractColumns(Matrix A, int col0, int col1) {
        return extractColumns(A, col0, col1, null);
    }

    /**
     * Extracts a range of columns from the matrix [col0:col1) into a destination matrix.
     *
     * @param A    The source matrix
     * @param col0 Starting column index (inclusive)
     * @param col1 Ending column index (exclusive)
     * @param out  Output matrix to hold the result
     * @return A matrix with the specified columns
     */
    public static Matrix extractColumns(Matrix A, int col0, int col1, Matrix out) {
        // Extract a range of columns from col0 (inclusive) to col1 (exclusive)
        // Validate bounds
        if (col0 < 0 || col1 > A.getNumCols() || col0 > col1) {
            throw new IllegalArgumentException("Invalid column range: col0=" + col0 + 
                ", col1=" + col1 + ", matrix has " + A.getNumCols() + " columns");
        }
        
        // If out is null, create it with appropriate dimensions
        if (out == null) {
            out = new Matrix(A.getNumRows(), col1 - col0);
        }
        
        // Extract the columns
        for (int row = 0; row < A.getNumRows(); row++) {
            for (int col = col0; col < col1; col++) {
                out.set(row, col - col0, A.get(row, col));
            }
        }
        return out;
    }

    /**
     * Extracts the diagonal elements of a matrix and stores them in a destination matrix.
     *
     * @param A       The source matrix
     * @param outputB The matrix to store the diagonal values
     */
    public static void extractDiag(Matrix A, Matrix outputB) {
        SparseMatrix.extractDiagStatic(A.getData(), outputB.getData());
    }

    /**
     * Extracts a range of rows from the matrix [row0:row1).
     *
     * @param A    The source matrix
     * @param row0 Starting row index (inclusive)
     * @param row1 Ending row index (exclusive)
     * @return A matrix containing the specified rows
     */
    public static Matrix extractRows(Matrix A, int row0, int row1) {
        Matrix out = new Matrix(row1 - row0, A.getNumCols());
        return extractRows(A, row0, row1, out);
    }

    /**
     * Extracts a range of rows from the matrix [row0:row1) into a destination matrix.
     *
     * @param A    The source matrix
     * @param row0 Starting row index (inclusive)
     * @param row1 Ending row index (exclusive)
     * @param out  (Optional) Output matrix to store the result; can be null
     * @return A matrix containing the extracted rows
     */
    public static Matrix extractRows(Matrix A, int row0, int row1, Matrix out) {
        if (out == null) {
            return new Matrix(SparseMatrix.extractRowsStatic(A.getData(), row0, row1, null));
        } else {
            SparseMatrix.extractRowsInPlaceStatic(A.getData(), row0, row1, out.getData());
            return out;
        }
    }

    /**
     * Creates an identity matrix of given size.
     *
     * @param length The number of rows and columns in the identity matrix
     * @return A sparse identity matrix of size {@code length x length}
     */
    public static Matrix eye(int length) {
        if (length < 0) {
            throw new IllegalArgumentException("Matrix dimensions must be non-negative, got: " + length);
        }
        return new Matrix(SparseMatrix.identityStatic(length));
    }

    /**
     * Computes the natural logarithm of the factorial (log(x!)) for each element in the input matrix.
     * The input matrix is assumed to contain non-negative integers or values appropriate for {@code factln()}.
     *
     * @param n Input matrix
     * @return A matrix where each element is {@code log(n_i!)}
     */
    public static Matrix factln(Matrix n) {
        Matrix ret = n.copy();
        for (int i = 0; i < n.getNumElements(); i++) {
            ret.set(i, Maths.factln(ret.get(i)));
        }
        return ret;
    }

    /**
     * Finds the indices of all rows or columns in a matrix that have a sum of zero.
     *
     * @param matrix The matrix to check
     * @param isRow  If true, checks rows; if false, checks columns
     * @return A list of indices corresponding to zero-sum rows or columns
     */
    public static List<Integer> findIndexWithZeroSum(Matrix matrix, boolean isRow) {
        List<Integer> zeroSumIndices = new ArrayList<>();
        int rows = matrix.getNumRows();
        int cols = matrix.getNumCols();

        if (isRow) {
            for (int row = 0; row < rows; row++) {
                double sum = 0.0;
                for (int col = 0; col < cols; col++) {
                    sum += matrix.get(row, col);
                }
                if (sum == 0.0) {
                    zeroSumIndices.add(row);
                }
            }
        } else {
            for (int col = 0; col < cols; col++) {
                double sum = 0.0;
                for (int row = 0; row < rows; row++) {
                    sum += matrix.get(row, col);
                }
                if (sum == 0.0) {
                    zeroSumIndices.add(col);
                }
            }
        }
        return zeroSumIndices;
    }

    /**
     * Finds the indices of all rows in a matrix that exactly match a given row vector.
     *
     * @param matrix The matrix to search
     * @param row    A 1-row matrix to compare against
     * @return A list of row indices in {@code matrix} that match {@code row}
     */
    public static List<Integer> findRows(Matrix matrix, Matrix row) {
        List<Integer> matchingRows = new ArrayList<>();

        // Validate inputs
        if (matrix == null || row == null) {
            return matchingRows;
        }
        if (matrix.getNumRows() == 0 || row.getNumRows() == 0) {
            return matchingRows;
        }
        if (matrix.getNumCols() != row.getNumCols()) {
            return matchingRows;
        }

        int numCols = matrix.getNumCols();

        for (int i = 0; i < matrix.getNumRows(); i++) {
            boolean match = true;
            for (int j = 0; j < numCols; j++) {
                if (matrix.get(i, j) != row.get(0, j)) {
                    match = false;
                    break;
                }
            }
            if (match) {
                matchingRows.add(i);
            }
        }
        return matchingRows;
    }

    /**
     * Computes the 1-norm (maximum absolute column sum) of the matrix.
     *
     * @param a The input matrix
     * @return The 1-norm of the matrix
     */
    public static double firstNorm(Matrix a) {
        Matrix b = a.copy();
        b.absEq();
        double norm = 0;
        for (int i = 0; i < b.getNumCols(); i++) {
            norm = FastMath.max(norm, b.sumCols(i));
        }
        return norm;
    }

    /**
     * Computes the total number of elements across all matrices in a map,
     * summing {@code numRows * numCols} for each matrix.
     *
     * @param cellArray A map of matrices
     * @return The total number of scalar elements across all matrices
     */
    public static int getColIndexSum(Map<Integer, Matrix> cellArray) {
        int sum = 0;
        for (Matrix matrix : cellArray.values()) {
            sum += matrix.getNumCols() * matrix.getNumRows();
        }
        return sum;
    }

    /**
     * Extracts a submatrix from the source matrix using zero-based index ranges.
     * The submatrix spans rows [x0:x1) and columns [y0:y1).
     *
     * @param sourceMatrix The source matrix
     * @param x0           Starting row index (inclusive)
     * @param x1           Ending row index (exclusive)
     * @param y0           Starting column index (inclusive)
     * @param y1           Ending column index (exclusive)
     * @return A new matrix containing the specified submatrix
     */
    public static Matrix getSubMatrix(Matrix sourceMatrix, int x0, int x1, int y0, int y1) {
        Matrix subMatrix = new Matrix(x1 - x0, y1 - y0);
        for (int i = x0; i < x1; i++) {
            for (int j = y0; j < y1; j++) {
                subMatrix.set(i - x0, j - y0, sourceMatrix.get(i, j));
            }
        }
        return subMatrix;
    }

    /**
     * Computes the infinity norm (maximum absolute row sum) of the matrix.
     *
     * @param a The input matrix
     * @return The infinity norm of the matrix
     */
    public static double infNorm(Matrix a) {
        Matrix b = a.copy();
        b.absEq();
        return b.sumRows().elementMax();
    }

    /**
     * Computes the infinity norm of this matrix.
     * The infinity norm is the maximum row sum of absolute values.
     *
     * @return The infinity norm of this matrix
     */
    public double infinityNorm() {
        return infNorm(this);
    }

    /**
     * Extracts a range of columns from this matrix (instance method).
     *
     * @param col0 The starting column index (inclusive)
     * @param col1 The ending column index (exclusive)
     * @return A new matrix containing columns [col0, col1)
     */
    public Matrix extractCols(int col0, int col1) {
        return extractColumns(this, col0, col1);
    }

    /**
     * Extracts rows from this matrix (instance method).
     *
     * @param row0 The starting row index (inclusive)
     * @param row1 The ending row index (exclusive)
     * @return A new matrix containing rows [row0, row1)
     */
    public Matrix extractRows(int row0, int row1) {
        return extractRows(this, row0, row1);
    }

    /**
     * Computes the intersection of scalar values present in two matrices.
     * The result is a list of distinct values that appear in both matrices.
     *
     * @param matrixA First input matrix
     * @param matrixB Second input matrix
     * @return A list of values common to both matrices
     */
    public static List<Double> intersect(Matrix matrixA, Matrix matrixB) {
        Set<Double> matrixAValues = new HashSet<>();
        Set<Double> matrixBValues = new HashSet<>();
        List<Double> outputValues = new LinkedList<>();

        int rows = matrixA.getNumRows();
        int cols = matrixA.getNumCols();
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                matrixAValues.add(matrixA.get(row, col));
            }
        }

        rows = matrixB.getNumRows();
        cols = matrixB.getNumCols();
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                matrixBValues.add(matrixB.get(row, col));
            }
        }

        for (double value : matrixAValues) {
            if (matrixBValues.contains(value)) {
                outputValues.add(value);
            }
        }

        return outputValues;
    }

    /**
     * Computes the inverse of the given matrix.
     * Delegates to the matrix's {@code inv()} method.
     *
     * @param m The input matrix
     * @return The inverse of the matrix
     */
    public static Matrix inv(Matrix m) {
        return m.inv();
    }

    /**
     * Computes the sum of the natural logarithms of all elements in the matrix.
     *
     * @param matrix The input matrix
     * @return The sum of log(x) over all elements x in the matrix
     * @throws ArithmeticException if any element is ≤ 0
     */
    public static double logSum(Matrix matrix) {
        double sum = 0;
        for (int i = 0; i < matrix.getNumRows(); i++) {
            for (int j = 0; j < matrix.getNumCols(); j++) {
                double val = matrix.get(i, j);
                if (val <= 0) throw new ArithmeticException("Log of non-positive number: " + val);
                sum += FastMath.log(val);
            }
        }
        return sum;
    }

    /**
     * Computes log(sum_i exp(x_i)) in a numerically stable way (log-sum-exp trick).
     * This is commonly used in probabilistic computations to avoid underflow.
     *
     * @param x Input matrix treated as a vector of values x₁, x₂, ..., xₙ
     * @return The log-sum-exp of the input values
     */
    public static double logsumexp(Matrix x) {
        int n = x.getNumElements();
        double a = x.elementMax();

        double s = 0;
        for (int i = 0; i < n; i++) {
            if (x.get(i) != a) {
                s += FastMath.exp(x.get(i) - a);
            }
        }
        return a + FastMath.log1p(s);
    }

    /**
     * Computes the solution to the continuous-time Lyapunov equation AX + XAᵀ + Q = 0.
     * Currently redirects to {@link #sylv} assuming it implements the solver.
     *
     * @param A The matrix A
     * @param B Ignored
     * @param C The matrix Q
     * @param D Ignored
     * @return The solution matrix X
     */
    public static Matrix lyap(Matrix A, Matrix B, Matrix C, Matrix D) {
        return sylv(A, B, C);
    }

    /**
     * Returns the index of the row in the matrix that exactly matches the given row vector.
     *
     * @param matrix The matrix to search
     * @param row    A single-row matrix to match
     * @return The index of the matching row, or -1 if no match is found
     */
    public static int matchrow(Matrix matrix, Matrix row) {
        if (matrix.getNumCols() != row.getNumCols()) return -1;
        for (int i = 0; i < matrix.getNumRows(); i++) {
            boolean rowsEqual = true;
            for (int j = 0; j < matrix.getNumCols(); j++) {
                if (matrix.get(i, j) != row.get(j)) {
                    rowsEqual = false;
                    break;
                }
            }
            if (rowsEqual) return i;
        }
        return -1;
    }

    /**
     * Adds a vector to each row or column of the matrix, depending on the vector's orientation.
     *
     * @param matrix The input matrix
     * @param vector A column vector (adds to each row) or a row vector (adds to each column)
     * @return A new matrix with the vector added
     * @throws RuntimeException if the vector shape is incompatible with the matrix
     */
    public static Matrix matrixAddVector(Matrix matrix, Matrix vector) {
        Matrix result = matrix.copy();
        if (vector.getNumCols() == 1 && vector.getNumRows() == matrix.getNumRows()) {
            for (int i = 0; i < vector.getNumRows(); i++) {
                result.rowIncrease(i, vector.get(i));
            }
        } else if (vector.getNumRows() == 1 && vector.getNumCols() == matrix.getNumCols()) {
            for (int i = 0; i < vector.getNumCols(); i++) {
                result.colIncrease(i, vector.get(i));
            }
        } else {
            throw new RuntimeException("The size of matrix and vector are not compatible");
        }
        return result;
    }

    /**
     * Computes the maximum relative absolute difference between corresponding elements
     * of two matrices: max(abs((a - b) / b)).
     *
     * @param a First matrix
     * @param b Second matrix (denominator for relative difference)
     * @return The maximum relative absolute difference
     */
    public static double maxAbsDiff(Matrix a, Matrix b) {
        double maxDiff = 0;
        for (int i = 0; i < a.getNumRows(); i++) {
            for (int j = 0; j < a.getNumCols(); j++) {
                double denom = b.get(i, j);
                double diff = denom != 0.0 ? FastMath.abs((a.get(i, j) - denom) / denom) : 0.0;
                if (diff > maxDiff) {
                    maxDiff = diff;
                }
            }
        }
        return maxDiff;
    }

    /**
     * Negates all elements in the matrix and returns the result.
     * Operates directly on the sparse structure for performance.
     *
     * @param a The input matrix
     * @return A new matrix with all values negated
     */
    public static Matrix negative(Matrix a) {
        Matrix b = a.copy();
        for (int i = 0; i < b.getNonZeros(); i++) {
            b.getNonZeroValues()[i] = -b.getNonZeroValues()[i];
        }
        return b;
    }

    /**
     * Computes the matrix 1 - A, where diagonal entries become 1 - A(i, i)
     * and off-diagonal entries become -A(i, j).
     *
     * @param matrix The input matrix
     * @return A matrix where each element is replaced by (1 - A(i, i)) on the diagonal and -A(i, j) elsewhere
     */
    public static Matrix oneMinusMatrix(Matrix matrix) {
        Matrix ret = matrix.copy();
        for (int i = 0; i < ret.getNumRows(); i++) {
            for (int j = 0; j < ret.getNumCols(); j++) {
                if (i == j) {
                    ret.set(i, j, 1 - ret.get(i, j));
                } else {
                    ret.set(i, j, -ret.get(i, j));
                }
            }
        }
        return ret;
    }

    /**
     * Decreases a single element of an integer vector by one.
     *
     * @param N The input vector
     * @param s The index to decrement
     * @return A new vector with element {@code s} decreased by one, if in bounds
     */
    public static Matrix oner(Matrix N, Integer s) {
        Matrix res = N.copy();
        if (s >= 0 && s < N.length()) {
            res.set(s, res.get(s) - 1);
        }
        return res;
    }

    /**
     * Decreases multiple elements of an integer vector by one.
     *
     * @param N The input vector
     * @param r A list of indices to decrement
     * @return A new vector with elements at positions in {@code r} decreased by one
     */
    public static Matrix oner(Matrix N, List<Integer> r) {
        Matrix res = N.copy();
        for (Integer s : r) {
            if (s >= 0 && s < N.length()) {
                res.set(s, res.get(s) - 1);
            }
        }
        return res;
    }

    /**
     * Creates a matrix of the given shape filled with ones.
     *
     * @param rows Number of rows
     * @param cols Number of columns
     * @return A matrix of shape (rows x cols) filled with ones
     */
    public static Matrix ones(int rows, int cols) {
        if (rows < 0 || cols < 0) {
            throw new IllegalArgumentException("Matrix dimensions must be non-negative, got: " + rows + "x" + cols);
        }
        Matrix e = new Matrix(rows, cols);
        e.ones();
        return e;
    }

    /**
     * Computes the matrix power A^b for a non-negative integer exponent b.
     * This is done by repeated multiplication and assumes b ≥ 0.
     *
     * @param a The base matrix
     * @param b The exponent (must be ≥ 0)
     * @return The matrix {@code a} raised to the power {@code b}
     */
    public static Matrix pow(Matrix a, int b) {
        if (b < 0) {
            throw new IllegalArgumentException("Exponent must be non-negative");
        } else if (b == 0) {
            return eye(a.getNumRows());  // assuming square
        }

        Matrix result = a.copy();
        for (int i = 1; i < b; i++) {
            result = result.mult(a);
        }
        return result;
    }

    /**
     * Reads a CSV-formatted matrix from a file. Each line is interpreted as a row.
     * Supports numeric values as well as "Inf" and "-Inf".
     *
     * @param fileName Path to the CSV file
     * @return A new Matrix containing the data read from the file
     * @throws RuntimeException if the file cannot be read or parsed
     */
    public static Matrix readFromFile(String fileName) {
        ArrayList<double[]> data = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                double[] row = new double[values.length];
                for (int i = 0; i < values.length; i++) {
                    if (values[i].equalsIgnoreCase("Inf")) {
                        row[i] = Inf;
                    } else if (values[i].equalsIgnoreCase("-Inf")) {
                        row[i] = NegInf;
                    } else {
                        row[i] = Double.parseDouble(values[i]);
                    }
                }
                data.add(row);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return new Matrix(data.toArray(new double[0][]));
    }

    public static double[][][] removeRows(double[][][] array, Matrix rowsToRemove) {
        if (rowsToRemove.getNumCols() == 0 || rowsToRemove.getNumRows() == 0) {
            return array;
        }
        List<Integer> rowsList = new ArrayList<>();
        for (int i = 0; i < rowsToRemove.getNumCols(); i++) {
            rowsList.add((int) rowsToRemove.get(i));
        }

        int newSize = array.length - rowsToRemove.getNumCols();
        double[][][] newArray = new double[newSize][array[0].length][array[0][0].length];

        int newIndex = 0;
        for (int i = 0; i < array.length; i++) {
            if (!rowsList.contains(i)) {
                newArray[newIndex++] = array[i];
            }
        }
        return newArray;
    }

    public static String removeTrailingNewLine(String str) {
        if (str.endsWith("\n")) {
            return str.substring(0, str.length() - 1);
        } else if (str.endsWith("\r\n")) {
            return str.substring(0, str.length() - 2);
        }
        return str;
    }

    /**
     * Multiplies all elements of a matrix by a scalar.
     * (Marked for removal — consider using {@code matrix.scaleEq(n)} directly.)
     *
     * @param a The input matrix
     * @param n The scalar multiplier
     * @return A new matrix equal to {@code a * n}
     */
    public static Matrix scaleMult(Matrix a, double n) {
        Matrix b = a.copy();
        b.scaleEq(n);
        return b;
    }

    /**
     * Creates a 1×1 matrix containing a single scalar value.
     *
     * @param value The scalar to store in the matrix
     * @return A 1×1 matrix with the given value
     */
    public static Matrix singleton(double value) {
        return new Matrix(1, 1, 1).fill(value);
    }

    /**
     * Solves the sparse linear system Ax = b.
     *
     * @param a The coefficient matrix A
     * @param b The right-hand side vector or matrix b
     * @param x The solution matrix x (output)
     * @return true if a solution was found, false otherwise
     */
    public static boolean solve(Matrix a, Matrix b, Matrix x) {
        if (a.getNumRows() != a.getNumCols()) {
            throw new IllegalArgumentException("Matrix A must be square for linear system solving");
        }
        if (a.getNumRows() != b.getNumRows()) {
            throw new IllegalArgumentException("Matrix dimensions incompatible: A is " + a.getNumRows() + "x" + a.getNumCols() +
                    " but b is " + b.getNumRows() + "x" + b.getNumCols());
        }

        // Check if matrix is singular
        double det = a.det();
        if (Math.abs(det) < 1e-14) {
            throw new RuntimeException("Matrix is singular and system cannot be solved (determinant = " + det + ")");
        }

        return SparseMatrix.solveStatic(a.getData(), b.getData(), x.getData());
    }

    /**
     * Solves the linear system Ax = b for x, handling singular matrices gracefully.
     * This method does not throw exceptions for singular matrices, instead returning false
     * and filling the result matrix with NaN values to match MATLAB behavior.
     *
     * @param a The coefficient matrix A
     * @param b The right-hand side vector/matrix b
     * @param x The solution matrix x (output)
     * @return true if a solution was found, false if matrix is singular
     */
    public static boolean solveSafe(Matrix a, Matrix b, Matrix x) {
        if (a.getNumRows() != a.getNumCols()) {
            throw new IllegalArgumentException("Matrix A must be square for linear system solving");
        }
        if (a.getNumRows() != b.getNumRows()) {
            throw new IllegalArgumentException("Matrix dimensions incompatible: A is " + a.getNumRows() + "x" + a.getNumCols() +
                    " but b is " + b.getNumRows() + "x" + b.getNumCols());
        }

        // Check if matrix is singular before attempting solve
        // For small matrices use determinant (fast), for larger matrices use rank (numerically stable)
        int n = a.getNumRows();
        if (n <= 3) {
            double det = a.det();
            if (Math.abs(det) < 1e-14) {
                x.reshape(n, b.getNumCols());
                x.fill(Double.NaN);
                return false;
            }
        } else {
            try {
                if (a.rank() < n) {
                    x.reshape(n, b.getNumCols());
                    x.fill(Double.NaN);
                    return false;
                }
            } catch (RuntimeException e) {
                // SVD decomposition can fail on large or ill-conditioned matrices;
                // skip rank check and try solving directly — check result for validity below
            }
        }

        boolean success = SparseMatrix.solveStatic(a.getData(), b.getData(), x.getData());

        // Check if solver failed (result contains NaN or Inf)
        if (!success || x.hasNaN() || x.hasInfinite()) {
            x.fill(Double.NaN);
            return false;
        }

        // Verify solution quality via residual check: ||A*x - b||_inf / ||b||_inf
        // This catches cases where LU returns a wrong non-NaN answer for near-singular systems
        // that slip past the rank check (e.g., when SVD fails and rank check is skipped).
        double bNormInf = 0;
        for (int i = 0; i < n; i++) {
            bNormInf = Math.max(bNormInf, Math.abs(b.get(i, 0)));
        }
        if (bNormInf > 0) {
            Matrix residual = a.mult(x).add(-1.0, b);
            double resNormInf = 0;
            for (int i = 0; i < n; i++) {
                resNormInf = Math.max(resNormInf, Math.abs(residual.get(i, 0)));
            }
            if (resNormInf / bNormInf > 1e-6) {
                x.fill(Double.NaN);
                return false;
            }
        }

        return true;
    }

    /**
     * Solves the linear system A*x = b directly using LU decomposition without
     * singularity checks. Faster than {@link #solveSafe} for cases where the
     * matrix is known to be non-singular or where the caller handles failures.
     *
     * @param a The coefficient matrix A (must be square)
     * @param b The right-hand side matrix b
     * @param x The solution matrix x (output)
     * @return true if a solution was found, false if the LU solver failed
     */
    public static boolean solveDirect(Matrix a, Matrix b, Matrix x) {
        return SparseMatrix.solveStatic(a.getData(), b.getData(), x.getData());
    }

    /**
     * Computes the spectral decomposition of a matrix A using its eigendecomposition.
     * Returns a diagonal matrix of eigenvalues and a cell array of spectral projectors.
     * For diagonalizable matrices: A = V * D * V⁻¹, where D is the spectrum and each projector is v_k * (v_k⁻¹)^T.
     * For defective matrices, falls back to Jordan decomposition or provides approximation.
     *
     * @param A The input matrix
     * @return A {@code SpectralDecomposition} containing the diagonal spectrum and associated projectors
     * @throws RuntimeException if the matrix cannot be decomposed
     */
    public static Ret.SpectralDecomposition spectd(Matrix A) {
        Ret.Eigs ret = A.eigvec();
        Matrix spectrum = Matrix.diagMatrix(ret.values);
        Matrix V = ret.vectors;

        try {
            // Attempt standard diagonalization
            Matrix iV = V.inv();
            MatrixCell projectors = new MatrixCell(ret.values.length());

            for (int k = 0; k < A.getNumRows(); k++) {
                Matrix proj = V.getColumn(k).mult(iV.getRow(k));
                projectors.set(k, proj);
            }

            return new Ret.SpectralDecomposition(spectrum, projectors);

        } catch (Exception e) {
            // Handle defective matrices (non-diagonalizable) - use warning not error since we have a fallback
            line_warning(mfilename(new Object[]{}), "Defective matrix in spectral decomposition, using approximate projectors: %s", e.getMessage());

            // For defective matrices, create approximate projectors
            // Use a fallback approach that creates orthogonal projectors based on eigenvalues
            MatrixCell projectors = new MatrixCell(ret.values.length());

            for (int k = 0; k < Math.min(A.getNumRows(), ret.values.length()); k++) {
                try {
                    // Create a simple projector: outer product of normalized eigenvector with itself
                    Matrix eigenvec = V.getColumn(k);
                    double norm = 0;
                    for (int i = 0; i < eigenvec.getNumRows(); i++) {
                        norm += eigenvec.get(i, 0) * eigenvec.get(i, 0);
                    }
                    norm = Math.sqrt(norm);

                    if (norm > 1e-12) {
                        // Normalize the eigenvector
                        Matrix normalizedVec = Matrix.scaleMult(eigenvec, 1.0 / norm);
                        // Create projector as v * v^T
                        Matrix proj = normalizedVec.mult(normalizedVec.transpose());
                        projectors.set(k, proj);
                    } else {
                        // Create zero projector for degenerate case
                        projectors.set(k, Matrix.zeros(A.getNumRows(), A.getNumCols()));
                    }
                } catch (Exception innerE) {
                    // If calculation fails, create zero projector
                    projectors.set(k, Matrix.zeros(A.getNumRows(), A.getNumCols()));
                }
            }

            return new Ret.SpectralDecomposition(spectrum, projectors);
        }
    }

    /**
     * Computes the sum of the cumulative product along a row vector.
     * For vector [a, b, c], returns a + ab + abc.
     *
     * @param matrix A 1-row matrix (row vector)
     * @return The scalar sum of cumulative products
     */
    public static double sumCumprod(Matrix matrix) {
        double sum = 0;
        double cumprod = 1;
        for (int i = 0; i < matrix.getNumCols(); i++) {
            cumprod *= matrix.get(0, i);
            sum += cumprod;
        }
        return sum;
    }

    /**
     * Robust linear solve for A·X = B that handles singular matrices.
     * Tries standard leftMatrixDivide first, falls back to pseudo-inverse on failure.
     *
     * @param A Coefficient matrix (may be singular)
     * @param B Right-hand side matrix
     * @return Solution matrix X
     */
    private static Matrix robustLeftDivide(Matrix A, Matrix B) {
        try {
            // Try standard solver first
            return A.leftMatrixDivide(B);
        } catch (RuntimeException e) {
            // If singular, use pseudo-inverse via SVD
            if (e.getMessage() != null && e.getMessage().contains("singular")) {
                // Compute pseudo-inverse using SVD
                Ret.SVD svd = A.svd();
                Matrix U = svd.u;
                Matrix S = svd.s;  // Column vector of singular values
                Matrix V = svd.v;

                // Create S+ (pseudo-inverse of S)
                int m = A.getNumRows();
                int n = A.getNumCols();
                Matrix Splus = new Matrix(n, m);

                double tol = 1e-10 * Math.max(m, n) * S.get(0, 0);  // Relative tolerance
                int rank = S.getNumRows();
                for (int i = 0; i < rank; i++) {
                    double sigma = S.get(i, 0);
                    if (Math.abs(sigma) > tol) {
                        Splus.set(i, i, 1.0 / sigma);
                    }
                }

                // A+ = V * S+ * U^T
                Matrix Aplus = V.mult(Splus).mult(U.transpose());
                return Aplus.mult(B);
            }
            throw e;  // Re-throw if not a singularity error
        }
    }

    /**
     * Solves the Sylvester equation A·X + X·B = -C using a Schur decomposition-based method.
     * Handles the case where B = Aᵀ with a specialized backward solve, otherwise proceeds forward.
     *
     * @param A Left coefficient matrix
     * @param B Right coefficient matrix
     * @param C Constant matrix
     * @return Solution matrix X
     */
    public static Matrix sylv(Matrix A, Matrix B, Matrix C) {
        int n = C.getNumCols();
        Map<String, Matrix> schur_decomposition_A = A.schur();
        Matrix ZA = schur_decomposition_A.get("U");
        Matrix TA = schur_decomposition_A.get("T");
        Matrix ZB;
        Matrix TB;
        String solver_direction;
        if (A.transpose().isEqualTo(B)) {
            ZB = ZA;
            TB = TA.transpose();
            solver_direction = "backward";
        } else {
            Map<String, Matrix> schur_decomposition_B = B.schur();
            ZB = schur_decomposition_B.get("U");
            TB = schur_decomposition_B.get("T");
            solver_direction = "forward";
        }

        Matrix F = ZA.transpose().mult(C).mult(ZB);

        //		if(TA.isDiag_withintol() && TB.isDiag_withintol()){
        //			Matrix TA_diag = new Matrix(0,0,0);
        //			Matrix TB_diag = new Matrix(0,0,0);
        //			Matrix.extractDiag(TA,TA_diag);
        //			Matrix.extractDiag(TB,TB_diag);
        //			Matrix L = TA.compatible_sizes_add(TB.transpose());
        //			for (int i=0;i<L.numRows;i++){
        //				for (int j=0;j<L.numCols;j++){
        //					L.set(i,j,-1/L.get(i,j));
        //				}
        //			}
        //			return ZA.mult(L.elementMult(F,null)).mult(ZB.transpose());
        //		}

        Matrix Y = new Matrix(C.getNumRows(), C.getNumCols(), C.getNumRows() * C.getNumCols());
        Matrix P = new Matrix(0, 0, 0);
        Matrix.extractDiag(TA, P);

        if (solver_direction.equals("backward")) {
            for (int k = n - 1; k > 0; k--) {
                Matrix rhs =
                        Matrix.extractColumn(F, k, null).add(1, Y.mult(Matrix.extractColumn(TB, k, null)));
                for (int i = 0; i < TA.getNumRows(); i++) {
                    TA.set(i, i, P.get(i) + TB.get(k, k));
                }
                Y.insertSubMatrix(
                        0, k, Y.getNumRows(), k + 1, robustLeftDivide(TA, Matrix.scaleMult(rhs, -1)));
            }
        } else {
            for (int k = 0; k < n; k++) {
                Matrix rhs =
                        Matrix.extractColumn(F, k, null).add(1, Y.mult(Matrix.extractColumn(TB, k, null)));
                for (int i = 0; i < TA.getNumRows(); i++) {
                    TA.set(i, i, P.get(i) + TB.get(k, k));
                }
                Y.insertSubMatrix(
                        0, k, Y.getNumRows(), k + 1, robustLeftDivide(TA, Matrix.scaleMult(rhs, -1)));
            }
        }

        return ZA.mult(Y).mult(ZB.transpose());
    }

    /**
     * Returns the lower triangular part of a matrix (zeroing elements above the main diagonal).
     *
     * @param matrix The input matrix
     * @return A lower triangular matrix with upper entries set to zero
     */
    public static Matrix tril(Matrix matrix) {
        return tril(matrix, 0);
    }

    /**
     * Returns the elements on and below the kth diagonal of matrix A.
     * For k = 0, returns the main diagonal and below.
     * For k = -1, returns only below the main diagonal (strict lower triangular).
     * For k = 1, returns the main diagonal, below, and first superdiagonal.
     *
     * @param matrix the input matrix
     * @param k the diagonal offset (0 = main diagonal, -1 = below main, 1 = above main)
     * @return matrix with elements on and below the kth diagonal, others set to zero
     */
    public static Matrix tril(Matrix matrix, int k) {
        int numRows = matrix.getNumRows();
        int numCols = matrix.getNumCols();
        Matrix result = matrix.copy();

        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                if (j > i + k) {
                    result.set(i, j, 0);
                }
            }
        }

        return result;
    }

    /**
     * Computes the union of two matrices A and B.
     * For each non-zero element in B, overwrites the corresponding value in A.
     * The resulting matrix is a clone of A with B's non-zero entries inserted.
     *
     * @param A The base matrix
     * @param B The matrix whose non-zero entries override A
     * @return A matrix representing the union of A and B
     */
    public static Matrix union(Matrix A, Matrix B) {
        Matrix unionMatrix = A.copy();
        for (int colB = 0; colB < B.getNumCols(); colB++) {
            for (int rowB = 0; rowB < B.getNumRows(); rowB++) {
                if (B.get(rowB, colB) != 0) {
                    unionMatrix.set(rowB, colB, B.get(rowB, colB));
                }
            }
        }
        return unionMatrix;
    }

    /**
     * Finds the indices of unique rows in a matrix without returning a matrix of the unique rows themselves.
     * Returns:
     * - vi: indices of the first occurrence of each unique row
     * - vj_map: maps unique row indices to a list of positions where that row occurs
     *
     * @param m The input matrix
     * @return A {@code UniqueRowResult} containing index mappings of unique rows
     */
    public static UniqueRowResult uniqueRowIndexes(Matrix m) {
        Map<RowHashKey, List<Integer>> rowToIndex = new HashMap<>();
        for (int i = 0; i < m.getNumRows(); i++) {
            RowHashKey key = new RowHashKey(m, i);
            if (!rowToIndex.containsKey(key)) {
                rowToIndex.put(key, new LinkedList<>());
            }
            rowToIndex.get(key).add(i);
        }

        List<Pair<RowHashKey, List<Integer>>> pairs =
                rowToIndex.entrySet().stream()
                        .map(e -> new Pair<>(e.getKey(), e.getValue()))
                        .sorted(
                                (i1, i2) -> {
                                    RowHashKey row1 = i1.getLeft();
                                    RowHashKey row2 = i2.getLeft();
                                    return row1.compareTo(row2);
                                })
                        .collect(Collectors.toList());

        Matrix vi = new Matrix(pairs.size(), 1);
        for (int i = 0; i < pairs.size(); i++) {
            vi.set(i, 0, pairs.get(i).getRight().get(0));
        }

        Matrix vj = new Matrix(m.getNumRows(), 1);
        for (int i = 0; i < pairs.size(); i++) {
            for (int j : pairs.get(i).getRight()) {
                vj.set(j, i);
            }
        }

        Map<Integer, List<Integer>> vj_map = new HashMap<>();
        for (int i = 0; i < vj.getNumElements(); i++) {
            if (!vj_map.containsKey((int) vj.get(i))) {
                vj_map.put((int) vj.get(i), new ArrayList<>());
            }
            vj_map.get((int) vj.get(i)).add(i);
        }

        return new UniqueRowResult(null, vi, vj_map);
    }

    /**
     * Identifies unique rows in a matrix starting from a specified column index.
     * Does not return the unique row data itself, only index mappings.
     * Uses string-based keys to handle state spaces that exceed integer hash range.
     *
     * @param m        The input matrix
     * @param startCol The column index from which uniqueness is considered
     * @return A {@code UniqueRowResult} containing:
     * - vi: indices of first occurrence of each unique row pattern
     * - vj_map: groupings of row indices by unique pattern
     */
    public static UniqueRowResult uniqueRowIndexesFromColumn(Matrix m, int startCol) {
        int numRows = m.getNumRows();
        int numCols = m.getNumCols();
        int partialLen = numCols - startCol;

        // Pre-compute all row keys in a single pass, storing them for reuse
        String[] rowKeys = new String[numRows];
        // Estimate StringBuilder capacity: ~20 chars per double value
        int estimatedCapacity = partialLen * 20;

        for (int i = 0; i < numRows; i++) {
            StringBuilder sb = new StringBuilder(estimatedCapacity);
            for (int j = startCol; j < numCols; j++) {
                if (j > startCol) sb.append(',');
                sb.append(m.get(i, j));
            }
            rowKeys[i] = sb.toString();
        }

        // Build map from row key to list of row indices
        Map<String, List<Integer>> rowToIndex = new HashMap<String, List<Integer>>();
        List<String> uniqueKeysOrdered = new ArrayList<String>();

        for (int i = 0; i < numRows; i++) {
            String rowKey = rowKeys[i];
            List<Integer> indices = rowToIndex.get(rowKey);
            if (indices == null) {
                indices = new LinkedList<Integer>();
                rowToIndex.put(rowKey, indices);
                uniqueKeysOrdered.add(rowKey);
            }
            indices.add(i);
        }

        // Sort unique keys lexicographically
        Collections.sort(uniqueKeysOrdered);

        // Build vi (first occurrence indices) and vj_map (row groupings)
        Matrix vi = new Matrix(uniqueKeysOrdered.size(), 1);
        Map<Integer, List<Integer>> vj_map = new HashMap<Integer, List<Integer>>();

        for (int i = 0; i < uniqueKeysOrdered.size(); i++) {
            String rowKey = uniqueKeysOrdered.get(i);
            List<Integer> indices = rowToIndex.get(rowKey);
            vi.set(i, 0, indices.get(0));
            vj_map.put(i, new ArrayList<Integer>(indices));
        }

        return new UniqueRowResult(null, vi, vj_map);
    }

    /**
     * Finds all unique rows in a matrix and returns the unique sorted rows,
     * along with mapping indices to/from the original matrix.
     * Uses efficient array-based hashing instead of Matrix objects as keys.
     *
     * @param m The input matrix
     * @return A {@code UniqueRowResult} containing:
     * - the matrix of sorted unique rows,
     * - vi: indices of first occurrences,
     * - vj_map: row groupings in the original matrix
     */
    public static UniqueRowResult uniqueRows(Matrix m) {
        double[][] arr = m.toArray2D();
        int numCols = m.getNumCols();

        // Use DoubleArrayWrapper for efficient hashing (O(n) hashCode vs O(n²) Matrix equality)
        Map<DoubleArrayWrapper, List<Integer>> rowToIndex = new HashMap<DoubleArrayWrapper, List<Integer>>();
        List<DoubleArrayWrapper> uniqueRowsOrdered = new ArrayList<DoubleArrayWrapper>();

        for (int i = 0; i < m.getNumRows(); i++) {
            DoubleArrayWrapper rowWrapper = new DoubleArrayWrapper(arr[i]);
            List<Integer> indices = rowToIndex.get(rowWrapper);
            if (indices == null) {
                indices = new LinkedList<Integer>();
                rowToIndex.put(rowWrapper, indices);
                uniqueRowsOrdered.add(rowWrapper);
            }
            indices.add(i);
        }

        // Sort unique rows lexicographically
        Collections.sort(uniqueRowsOrdered, new Comparator<DoubleArrayWrapper>() {
            @Override
            public int compare(DoubleArrayWrapper w1, DoubleArrayWrapper w2) {
                double[] row1 = w1.getArray();
                double[] row2 = w2.getArray();
                int minLen = Math.min(row1.length, row2.length);
                for (int k = 0; k < minLen; k++) {
                    int cmp = Double.compare(row1[k], row2[k]);
                    if (cmp != 0) return cmp;
                }
                return Integer.compare(row1.length, row2.length);
            }
        });

        // Build vi (first occurrence indices) and vj_map (row groupings)
        Matrix vi = new Matrix(uniqueRowsOrdered.size(), 1);
        Map<Integer, List<Integer>> vj_map = new HashMap<Integer, List<Integer>>();

        for (int i = 0; i < uniqueRowsOrdered.size(); i++) {
            DoubleArrayWrapper rowWrapper = uniqueRowsOrdered.get(i);
            List<Integer> indices = rowToIndex.get(rowWrapper);
            vi.set(i, 0, indices.get(0));
            vj_map.put(i, new ArrayList<Integer>(indices));
        }

        // Create sorted unique rows matrix
        double[][] sortedData = new double[uniqueRowsOrdered.size()][numCols];
        for (int i = 0; i < uniqueRowsOrdered.size(); i++) {
            System.arraycopy(uniqueRowsOrdered.get(i).getArray(), 0, sortedData[i], 0, numCols);
        }

        return new UniqueRowResult(new Matrix(sortedData), vi, vj_map);
    }

    /**
     * Lightweight wrapper for double[] arrays to use as HashMap keys.
     * Provides efficient hashCode() and equals() implementations using Arrays utility methods.
     */
    private static final class DoubleArrayWrapper {
        private final double[] array;
        private final int hashCode;

        public DoubleArrayWrapper(double[] array) {
            this.array = array;
            this.hashCode = Arrays.hashCode(array);
        }

        public double[] getArray() {
            return array;
        }

        @Override
        public int hashCode() {
            return hashCode;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof DoubleArrayWrapper)) return false;
            DoubleArrayWrapper other = (DoubleArrayWrapper) obj;
            return Arrays.equals(this.array, other.array);
        }
    }

    /**
     * Weakly-connected components of a sub-matrix.
     * Constructs an undirected graph from the input matrix, ignoring specified columns,
     * and returns the set of weakly connected components.
     *
     * @param param        Input matrix
     * @param colsToIgnore Indexes to be ignored
     * @return A set of sets, each representing a weakly connected component
     */
    public static Set<Set<Integer>> weaklyConnect(Matrix param, Set<Integer> colsToIgnore) {
        if (param.hasNaN()) {
            param.setNaNTo(1.0);
        }
        UndirectedGraph graph = new UndirectedGraph(param, colsToIgnore);
        graph.computeWeaklyConnectedComponents();
        return graph.getWCC();
    }

    /**
     * Creates a matrix of the specified shape filled with zeros.
     *
     * @param rows Number of rows
     * @param cols Number of columns
     * @return A zero-initialized matrix of size (rows x cols)
     */
    public static Matrix zeros(int rows, int cols) {
        if (rows < 0 || cols < 0) {
            throw new IllegalArgumentException("Matrix dimensions must be non-negative, got: " + rows + "x" + cols);
        }
        return new Matrix(rows, cols);
    }

    // ================================================================================
    // ================================================================================
    // SECTION 3: BASIC ARITHMETIC OPERATIONS
    // ================================================================================
    // Addition, subtraction, multiplication, division operations
    // ================================================================================
    // Methods for basic matrix arithmetic: add, subtract, multiply, divide
    
    public void absEq() {
        delegate.absEq();
    }

    /**
     * Adds another matrix to this matrix: {@code this + matrix}.
     *
     * @param matrix The matrix to add
     * @return A new matrix containing the sum
     */
    public Matrix add(Matrix matrix) {
        return add(1.0, matrix);
    }

    /**
     * Adds a scalar multiple of the all-ones matrix to this matrix: {@code this + alpha * 1}.
     *
     * @param alpha The scalar multiplier
     * @return A new matrix representing the addition
     */
    public Matrix add(double alpha) {
        return new Matrix(
                addMatrices(
                        1,
                        this.getData(),
                        alpha,
                        Matrix.ones(this.getNumRows(), this.getNumCols()).getData(),
                        null));
    }

    /**
     * Adds {@code alpha * matrix} to this matrix.
     *
     * @param alpha  The scalar multiplier
     * @param matrix The matrix to add
     * @return A new matrix containing the result
     */
    public Matrix add(double alpha, Matrix matrix) {
        return new Matrix(addMatrices(1, this.getData(), alpha, matrix.getData(), null));
    }

    /**
     * Adds a scalar value to each element of the matrix in-place.
     *
     * @param alpha The scalar value to add
     */
    public void addEq(double alpha) {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
                this.set(row, col, this.getData().get(row, col) + alpha);
            }
        }
    }

    /**
     * Adds another matrix to this matrix in-place: {@code this += matrix}.
     *
     * @param matrix The matrix to add
     */
    public void addEq(Matrix matrix) {
        addEq(1.0, matrix);
    }

    /**
     * Adds {@code alpha * matrix} to this matrix in-place: {@code this += alpha * matrix}.
     *
     * @param alpha  The scalar multiplier
     * @param matrix The matrix to add
     */
    public void addEq(double alpha, Matrix matrix) {
        this.setTo(this.add(alpha, matrix));
    }

    private DMatrix addMatrices(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output) {
        return delegate.addMatrices(alpha, A, beta, B, output);
    }

    // ================================================================================
    // ================================================================================
    // SECTION 10: MATRIX PROPERTIES AND COMPARISONS
    // ================================================================================
    // Methods for checking matrix properties and comparing matrices
    // ================================================================================
    // Methods for checking matrix properties and comparing matrices
    
    /**
     * Checks whether all elements in the matrix are equal to 1.
     *
     * @return true if all elements are exactly 1.0, false otherwise
     */
    public boolean allEqualToOne() {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
                if (this.get(row, col) != 1) {
                    return false;
                }
            }
        }
        return true;
    }

    public boolean any() {
        return delegate.any();
    }

    /**
     * Applies a conditional element-wise transformation to the matrix.
     * For each element matching a condition based on {@code source} and {@code op}, the value is replaced with {@code target}.
     * <p>
     * Supported operations include: "equal", "notequal", "great", "greatequal", "less", "lessequal".
     * <p>
     * This method safely handles special cases like zero, NaN, and Infinity.
     *
     * @param source The reference value for the comparison
     * @param target The value to assign if the condition holds
     * @param op     The comparison operator
     */
    public void apply(double source, double target, String op) {
        double tol = GlobalConstants.Zero;
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        switch (op) {
            case "equal":
                if (Math.abs(source - 0) < tol) {
                    if (Math.abs(target - 0) < tol) return;
                    for (int i = 0; i < numRows; i++) {
                        for (int j = 0; j < numCols; j++) {
                            if (Math.abs(this.get(i, j) - 0) < tol) getData().set(i, j, target);
                        }
                    }
                } else if (Double.isNaN(source)) {
                    applyConditionalTransform(source, target, tol, "NaN");
                } else if (isInf(source)) {
                    applyConditionalTransform(source, target, tol, "Inf");
                } else {
                    applyConditionalTransform(source, target, tol, "Value");
                }
                break;

            case "notequal":
                if (Math.abs(source - 0) < tol) {
                    if (Math.abs(target - 0) < tol) this.getData().zero();
                    for (int colIdx = 0; colIdx < numCols; colIdx++) {
                        int col1 = getColumnIndex(colIdx);
                        int col2 = getColumnIndex(colIdx + 1);
                        for (int i = col1; i < col2; i++) {
                            if ((Math.abs(this.getNonZeroValue(i) - 0) >= tol)
                                    || (Double.isNaN(this.getNonZeroValue(i)))) {
                                getData().set(this.getNonZeroRow(i), colIdx, target);
                            }
                        }
                    }
                } else if (Double.isNaN(source)) {
                    for (int row = 0; row < numRows; row++) {
                        for (int col = 0; col < numCols; col++) {
                            if (!Double.isNaN(this.get(row, col))) {
                                getData().set(row, col, target);
                            }
                        }
                    }
                } else if (isInf(source)) {
                    for (int row = 0; row < numRows; row++) {
                        for (int col = 0; col < numCols; col++) {
                            if (!isInf(this.get(row, col))) {
                                getData().set(row, col, target);
                            }
                        }
                    }
                } else {
                    for (int row = 0; row < numRows; row++) {
                        for (int col = 0; col < numCols; col++) {
                            if ((Math.abs(this.get(row, col) - source) >= tol)
                                    || (Double.isNaN(this.get(row, col)))) {
                                getData().set(row, col, target);
                            }
                        }
                    }
                }
                break;

            case "great":
                if (Math.abs(source - 0) < tol) {
                    for (int i = 0; i < numRows; i++) {
                        for (int j = 0; j < numCols; j++) {
                            if (Math.abs(this.get(i, j) - 0) >= tol && Double.compare(this.get(i, j), 0) > 0) {
                                getData().set(i, j, target);
                            }
                        }
                    }
                } else if (Double.isNaN(source)) {
                    throw new RuntimeException("Cannot compare with NaN");
                } else if (isInf(source)) {
                    throw new RuntimeException("Cannot compare with Infinite");
                } else {
                    for (int row = 0; row < numRows; row++) {
                        for (int col = 0; col < numCols; col++) {
                            if (Math.abs(this.get(row, col) - source) >= tol
                                    && Double.compare(this.get(row, col), source) > 0) {
                                getData().set(row, col, target);
                            }
                        }
                    }
                }
                break;

            case "greatequal":
                if (Math.abs(source - 0) < tol) {
                    for (int i = 0; i < numRows; i++) {
                        for (int j = 0; j < numCols; j++) {
                            if ((Math.abs(this.get(i, j) - 0) < tol)
                                    || (Math.abs(this.get(i, j) - 0) >= tol
                                    && Double.compare(this.get(i, j), 0) > 0)) {
                                getData().set(i, j, target);
                            }
                        }
                    }
                } else if (Double.isNaN(source)) {
                    throw new RuntimeException("Cannot compare with NaN");
                } else if (isInf(source)) {
                    throw new RuntimeException("Cannot compare with Infinite");
                } else {
                    for (int row = 0; row < numRows; row++) {
                        for (int col = 0; col < numCols; col++) {
                            if ((Math.abs(this.get(row, col) - source) < tol)
                                    || (Math.abs(this.get(row, col) - source) >= tol
                                    && Double.compare(this.get(row, col), source) > 0)) {
                                getData().set(row, col, target);
                            }
                        }
                    }
                }
                break;

            case "less":
                if (Math.abs(source - 0) < tol) {
                    for (int i = 0; i < numRows; i++) {
                        for (int j = 0; j < numCols; j++) {
                            if (Math.abs(this.get(i, j) - 0) >= tol && Double.compare(this.get(i, j), 0) < 0) {
                                getData().set(i, j, target);
                            }
                        }
                    }
                } else if (Double.isNaN(source)) {
                    throw new RuntimeException("Cannot compare with NaN");
                } else if (isInf(source)) {
                    throw new RuntimeException("Cannot compare with Infinite");
                } else {
                    for (int row = 0; row < numRows; row++) {
                        for (int col = 0; col < numCols; col++) {
                            if (Math.abs(this.get(row, col) - source) >= tol
                                    && Double.compare(this.get(row, col), source) < 0) {
                                getData().set(row, col, target);
                            }
                        }
                    }
                }
                break;

            case "lessequal":
                if (Math.abs(source - 0) < tol) {
                    for (int i = 0; i < numRows; i++) {
                        for (int j = 0; j < numCols; j++) {
                            if ((Math.abs(this.get(i, j) - 0) < tol)
                                    || (Math.abs(this.get(i, j) - 0) >= tol
                                    && Double.compare(this.get(i, j), 0) < 0)) {
                                getData().set(i, j, target);
                            }
                        }
                    }
                } else if (Double.isNaN(source)) {
                    throw new RuntimeException("Cannot compare with NaN");
                } else if (isInf(source)) {
                    throw new RuntimeException("Cannot compare with Infinite");
                } else {
                    for (int row = 0; row < numRows; row++) {
                        for (int col = 0; col < numCols; col++) {
                            if ((Math.abs(this.get(row, col) - source) < tol)
                                    || (Math.abs(this.get(row, col) - source) >= tol
                                    && Double.compare(this.get(row, col), source) < 0)) {
                                getData().set(row, col, target);
                            }
                        }
                    }
                }
                break;

            default:
                throw new RuntimeException("Operation is not supported");
        }

        if (target == 0) removeZeros();
    }

    private void applyConditionalTransform(double source, double target, double tol, String operation) {
        delegate.applyConditionalTransform(source, target, tol, operation);
    }

    /**
     * Returns a new matrix where each element is the ceiling of the corresponding element in this matrix.
     *
     * @return A new matrix with the ceiling applied element-wise.
     */
    public Matrix ceil() {
        Matrix output = new Matrix(this);
        int numRows = output.getNumRows();
        int numCols = output.getNumCols();
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                output.set(i, j, FastMath.ceil(output.get(i, j)));
            }
        }
        return output;
    }

    /**
     * Applies the ceiling operation in-place to each element of this matrix.
     *
     * @return This matrix after applying ceiling element-wise.
     */
    public Matrix ceilEq() {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                this.set(i, j, FastMath.ceil(this.get(i, j)));
            }
        }
        return this;
    }

    /**
     * Negates all values in the matrix, in-place.
     */
    public void changeSign() {
        delegate.changeSign();
    }

    /**
     * Returns a deep copy of this matrix.
     *
     * @return A cloned matrix identical to this one.
     */
    public Matrix copy() {
        return new Matrix(this);
    }

    public void colIncrease(int col, double a) {
        delegate.colIncrease(col, a);
    }

    /**
     * Returns the matrix in column-major order.
     *
     * @return The same matrix interpreted in column-major order.
     */
    public Matrix colon() {
        return this.columnMajorOrder();
    }

    /**
     * Equivalent to the colon operator in MATLAB (:).
     * Flattens the matrix in column-major order into a single column vector.
     *
     * @return A column vector representing the matrix in column-major order.
     */
    public Matrix columnMajorOrder() {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        Matrix output = new Matrix(numRows * numCols, 1);
        int idx = 0;
        for (int j = 0; j < numCols; j++) {
            for (int i = 0; i < numRows; i++) {
                output.set(idx++, 0, this.get(i, j));
            }
        }
        return output;
    }

    /**
     * Compares this matrix to another matrix for approximate equality.
     * Tolerance is defined by GlobalConstants.FineTol.
     *
     * @param matrix The matrix to compare against.
     * @return true if all elements are approximately equal, false otherwise.
     */
    public boolean compareMatrix(Matrix matrix) {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        if (numRows != matrix.getNumRows() || numCols != matrix.getNumCols()) {
            return false;
        }

        boolean check = true;
        double epsilon = GlobalConstants.FineTol;
        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
                double diff = Math.abs(this.get(row, col) - matrix.get(row, col));
                if (diff > epsilon) {
                    check = false;
                    break;
                }
            }
            if (!check) break;
        }
        return check;
    }

    /**
     * Performs element-wise addition with shape broadcasting where applicable.
     *
     * @param b Matrix to be added.
     * @return Result of the addition with compatible broadcasting.
     */
    public Matrix compatibleSizesAdd(Matrix b) {
        if (getNumRows() == b.getNumRows() && getNumCols() == b.getNumCols()) {
            return this.add(1, b);
        }
        if (this.length() == 1) return b.elementIncrease(this.get(0));
        if (b.length() == 1) return this.elementIncrease(b.get(0));

        if ((getNumCols() == 1 || getNumRows() == 1) && b.getNumCols() != 1 && b.getNumRows() != 1) {
            return matrixAddVector(b, this);
        }
        if ((b.getNumCols() == 1 || b.getNumRows() == 1) && getNumCols() != 1 && getNumRows() != 1) {
            return matrixAddVector(this, b);
        }

        if (getNumCols() == 1 && b.getNumRows() == 1) {
            return broadcastColPlusRow(this, b);
        }
        if (getNumRows() == 1 && b.getNumCols() == 1) {
            return broadcastColPlusRow(b, this);
        }

        throw new RuntimeException("The size of matrices are not compatible");
    }

    // ================================================================================
    // ================================================================================
    // SECTION 8: MATRIX CONCATENATION AND ASSEMBLY
    // ================================================================================
    // Methods for combining matrices: concatenation, block operations
    // ================================================================================
    // Methods for combining matrices: concat, block diagonal, etc.
    
    /**
     * Concatenates this matrix with another matrix horizontally.
     *
     * @param other Matrix to concatenate.
     * @return A new matrix with columns of both matrices concatenated.
     */
    public Matrix concatCols(Matrix other) {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        int otherNumCols = other.getNumCols();
        if (numRows != other.getNumRows()) {
            throw new IllegalArgumentException("Matrices must have the same number of rows.");
        }

        Matrix res = new Matrix(numRows, numCols + otherNumCols);

        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                res.set(i, j, this.get(i, j));
            }
            for (int j = 0; j < otherNumCols; j++) {
                res.set(i, numCols + j, other.get(i, j));
            }
        }

        return res;
    }

    // Additional methods needed by Java code
    private void concatColumnsInPlace(DMatrix left, DMatrix right, DMatrix output) {
        delegate.concatColumnsInPlace(left, right, output);
    }

    private void concatRowsInPlace(DMatrix top, DMatrix bottom, DMatrix output) {
        delegate.concatRowsInPlace(top, bottom, output);
    }

    /**
     * Counts how many elements in the matrix are equal to a specified value.
     *
     * @param val Value to count.
     * @return Number of occurrences of the value.
     */
    public int count(double val) {
        if (val == 0) {
            return this.getNumCols() * this.getNumRows() - ((DMatrixSparseCSC) getData()).getNonZeroLength();
        }
        int res = 0;
        for (int i = 0; i < this.getNonZeros(); i++) {
            if (this.getNonZeroValues()[i] == val) res++;
        }
        return res;
    }

    /**
     * Counts the number of occurrences of a value in each row of the matrix.
     *
     * @param val Value to count.
     * @return A column vector where each entry is the count for the corresponding row.
     */
    public Matrix countEachRow(double val) {
        BaseMatrix tempRes = delegate.countEachRow(val);
        Matrix res;
        if (tempRes instanceof SparseMatrix) {
            res = new Matrix(((SparseMatrix) tempRes).getData());
        } else {
            throw new UnsupportedOperationException("Non-sparse matrix results not yet supported");
        }

        if (val == 0) {
            int numRows = this.getNumRows();
            int numCols = this.getNumCols();
            for (int i = 0; i < numRows; i++) {
                for (int j = 0; j < numCols; j++) {
                    if (this.get(i, j) == 0) {
                        res.set(i, 0, res.get(i, 0) + 1);
                    }
                }
            }
        }

        return res;
    }

    /**
     * Creates a block diagonal matrix by placing the current matrix in the top-left
     * and another matrix (if provided) in the bottom-right.
     *
     * @param matrix2 The second matrix to place on the block diagonal. Can be null.
     * @return A new block diagonal matrix combining this and matrix2.
     */
    public Matrix createBlockDiagonal(Matrix matrix2) {
        int m1rows = this.getNumRows();
        int m2rows = 0;
        int m1cols = this.getNumCols();
        int m2cols = 0;

        if (matrix2 != null) {
            m2rows = matrix2.getNumRows();
            m2cols = matrix2.getNumCols();
        }

        Matrix output = new Matrix(m1rows + m2rows, m1cols + m2cols);
        for (int i = 0; i < m1rows; i++) {
            for (int j = 0; j < m1cols; j++) {
                output.set(i, j, this.get(i, j));
            }
        }
        for (int i = 0; i < m2rows; i++) {
            for (int j = 0; j < m2cols; j++) {
                output.set(i + m1rows, j + m1cols, matrix2.get(i, j));
            }
        }

        return output;
    }

    /**
     * Computes the cumulative sum of the matrix down each column.
     *
     * @return A matrix where each element is the cumulative sum along its column.
     */
    public Matrix cumsumViaCol() {
        Matrix res = new Matrix(this.getNumRows(), this.getNumCols(), this.getNumRows() * this.getNumCols());
        for (int i = 0; i < this.getNumCols(); i++) res.set(0, i, this.get(0, i));

        for (int i = 0; i < this.getNumCols(); i++) {
            for (int j = 1; j < this.getNumRows(); j++) {
                res.set(j, i, this.get(j, i) + res.get(j - 1, i));
            }
        }
        return res;
    }

    /**
     * Computes the cumulative sum of the matrix across each row.
     *
     * @return A matrix where each element is the cumulative sum along its row.
     */
    public Matrix cumsumViaRow() {
        Matrix res = new Matrix(this.getNumRows(), this.getNumCols(), this.getNumRows() * this.getNumCols());
        if (this.getNumRows() == 0 || this.getNumCols() == 0) {
            return res;
        }
        for (int i = 0; i < this.getNumRows(); i++) res.set(i, 0, this.get(i, 0));

        for (int i = 0; i < this.getNumRows(); i++) {
            for (int j = 1; j < this.getNumCols(); j++) {
                res.set(i, j, this.get(i, j) + res.get(i, j - 1));
            }
        }
        return res;
    }

    // ================================================================================
    // ================================================================================
    // SECTION 6: LINEAR ALGEBRA OPERATIONS
    // ================================================================================
    // Determinant, inverse, solve, decompositions, eigenvalues
    // ================================================================================
    // Methods for linear algebra: solve, inverse, determinant, eigenvalues, etc.
    
    /**
     * Computes the determinant of the matrix.
     *
     * @return Determinant value.
     */
    public double det() {
        if (getNumRows() != getNumCols()) {
            throw new IllegalArgumentException("Determinant can only be computed for square matrices, got: " +
                    getNumRows() + "x" + getNumCols());
        }
        return delegate.determinant();
    }

    /**
     * Performs element-wise division between this matrix and the provided matrix.
     *
     * @param den The denominator matrix.
     * @return A new matrix resulting from element-wise division.
     */
    public Matrix div(Matrix den) {
        Matrix ret = this.copy();
        for (int i = 0; i < den.getNumElements(); i++) {
            ret.set(i, ret.get(i) / den.get(i));
        }
        return ret;
    }

    /**
     * Performs in-place element-wise division between this matrix and the provided matrix.
     *
     * @param den The denominator matrix.
     */
    public void divEq(Matrix den) {
        for (int i = 0; i < den.getNumElements(); i++) {
            this.set(i, this.get(i) / den.get(i));
        }
    }

    /**
     * Divides this matrix by a scalar and stores the result in the output matrix.
     *
     * @param scalar  The scalar divisor.
     * @param outputB The matrix to store the result.
     * @param flag    If true, performs matrix / scalar; if false, scalar / matrix.
     */
    public void divide(double scalar, Matrix outputB, boolean flag) {
        if (flag) divideMatrix(scalar, outputB.getData());
        else divideScalarByMatrix(scalar, outputB.getData());
    }

    /**
     * Performs in-place division of this matrix by a scalar.
     *
     * @param scalar The scalar divisor.
     */
    public void divideEq(double scalar) {
        divideInPlace(scalar);
    }

    private void divideInPlace(double scalar) {
        delegate.divideInPlace(scalar);
    }

    private void divideMatrix(double scalar, DMatrix output) {
        delegate.divideMatrix(scalar, output);
    }

    /**
     * Divides each row of this matrix by the corresponding diagonal element (with an offset).
     *
     * @param diag   Array of diagonal values.
     * @param offset Starting offset in the diagonal array.
     */
    public void divideRows(double[] diag, int offset) {
        divideRowsByArray(diag, offset);
    }

    private void divideRowsByArray(double[] diag, int offset) {
        delegate.divideRowsByArray(diag, offset);
    }

    private void divideScalarByMatrix(double scalar, DMatrix output) {
        delegate.divideScalarByMatrix(scalar, output);
    }

    /**
     * Computes the eigenvalues of this square matrix.
     *
     * @return A {@link Ret.Eigs} object containing the eigenvalues as a column matrix.
     * If the matrix has complex eigenvalues, the result contains NaN values.
     */
    public Ret.Eigs eigval() {

        if (this.getNumCols() != this.getNumRows()) {
            throw new RuntimeException("Only square matrix can be eigen decomposited");
        }
        DMatrixRMaj matrix = new DMatrixRMaj(this.getNumRows(), this.getNumCols());
        DConvertMatrixStruct.convert(this.getData(), matrix);

        EigenDecomposition_F64<DMatrixRMaj> eig = DecompositionFactory_DDRM.eig(matrix.numCols, false);
        eig.decompose(matrix);

        Matrix D = new Matrix(1, getNumCols(), getNumCols());
        for (int i = 0; i < getNumCols(); i++) {
            if (eig.getEigenvalue(i).isReal()) {
                D.set(i, eig.getEigenvalue(i).getReal());
            } else {
                D.set(i, Double.NaN);
            }
        }

        return new Ret.Eigs(D, null);
    }

    /**
     * Computes the eigenvalues of this square matrix, including complex eigenvalues.
     *
     * @return A list of complex eigenvalues (using Apache Commons Math3 Complex class).
     */
    public List<org.apache.commons.math3.complex.Complex> eig() {
        if (this.getNumCols() != this.getNumRows()) {
            throw new RuntimeException("Only square matrix can be eigen decomposed");
        }
        DMatrixRMaj matrix = new DMatrixRMaj(this.getNumRows(), this.getNumCols());
        DConvertMatrixStruct.convert(this.getData(), matrix);

        EigenDecomposition_F64<DMatrixRMaj> eigDecomp = DecompositionFactory_DDRM.eig(matrix.numCols, false);
        eigDecomp.decompose(matrix);

        List<org.apache.commons.math3.complex.Complex> eigenvalues = new ArrayList<>();
        for (int i = 0; i < getNumCols(); i++) {
            org.ejml.data.Complex_F64 ev = eigDecomp.getEigenvalue(i);
            eigenvalues.add(new org.apache.commons.math3.complex.Complex(ev.getReal(), ev.getImaginary()));
        }
        return eigenvalues;
    }

    /**
     * Computes the eigenvalues and eigenvectors of this square matrix.
     *
     * @return A {@link Ret.Eigs} object containing the eigenvalues and eigenvectors.
     */
    public Ret.Eigs eigvec() {
        if (this.getNumCols() != this.getNumRows()) {
            throw new RuntimeException("Only a square matrix can be eigen-decomposed");
        }
        RealMatrix matrix = MatrixUtils.createRealMatrix(this.toArray2D());
        org.apache.commons.math3.linear.EigenDecomposition eigenDecomposition =
                new org.apache.commons.math3.linear.EigenDecomposition(matrix);
        double[] eigenvalues = eigenDecomposition.getRealEigenvalues();
        RealMatrix eigenvectors = eigenDecomposition.getV();

        Matrix V = new Matrix(eigenvalues.length, eigenvalues.length);
        for (int i = 0; i < eigenvectors.getRowDimension(); i++) {
            V.setRow(i, new Matrix(eigenvectors.getRowVector(i).toArray()));
        }
        return new Ret.Eigs(new Matrix(Arrays.stream(eigenvalues).toArray()).transpose(), V);
    }

    /**
     * Computes the Singular Value Decomposition (SVD) of this matrix.
     * The decomposition is A = U * S * V^T where:
     * - U contains the left singular vectors
     * - S is a diagonal matrix of singular values
     * - V contains the right singular vectors
     *
     * @return A {@link Ret.SVD} object containing U, S (as column vector), and V matrices.
     */
    public Ret.SVD svd() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(this.toArray2D());
        org.apache.commons.math3.linear.SingularValueDecomposition svd =
                new org.apache.commons.math3.linear.SingularValueDecomposition(matrix);

        double[] singularValues = svd.getSingularValues();
        RealMatrix uMatrix = svd.getU();
        RealMatrix vMatrix = svd.getV();

        Matrix U = new Matrix(uMatrix.getRowDimension(), uMatrix.getColumnDimension());
        for (int i = 0; i < uMatrix.getRowDimension(); i++) {
            for (int j = 0; j < uMatrix.getColumnDimension(); j++) {
                U.set(i, j, uMatrix.getEntry(i, j));
            }
        }

        Matrix S = new Matrix(singularValues.length, 1);
        for (int i = 0; i < singularValues.length; i++) {
            S.set(i, 0, singularValues[i]);
        }

        Matrix V = new Matrix(vMatrix.getRowDimension(), vMatrix.getColumnDimension());
        for (int i = 0; i < vMatrix.getRowDimension(); i++) {
            for (int j = 0; j < vMatrix.getColumnDimension(); j++) {
                V.set(i, j, vMatrix.getEntry(i, j));
            }
        }

        return new Ret.SVD(U, S, V);
    }

    /**
     * Computes the matrix exponential.
     * This is a convenience method that delegates to expm_higham().
     *
     * @return The matrix exponential e^A
     * @throws IllegalArgumentException if the matrix is not square
     */
    public Matrix expm() {
        return expm_higham();
    }

    /**
     * Computes the matrix exponential using Higham's scaling and squaring method.
     * This is an implementation of the algorithm from:
     * Higham, N.J. (2005). "The scaling and squaring method for the matrix exponential revisited."
     * SIAM Journal on Matrix Analysis and Applications, 26(4), 1179-1193.
     *
     * @return The matrix exponential e^A
     * @throws IllegalArgumentException if the matrix is not square
     */
    public Matrix expm_higham() {
        if (getNumRows() != getNumCols()) {
            throw new IllegalArgumentException("Matrix exponential can only be computed for square matrices");
        }
        
        int n = getNumRows();
        
        // Balance the matrix first
        // For now, we'll implement a simple version without actual balancing
        // A full implementation would use LAPACK's DGEBAL routine
        Matrix A = this.copy();
        double[] scale = new double[n];
        Arrays.fill(scale, 1.0);
        int ilo = 0;
        int ihi = n - 1;
        
        double nrmA = A.norm1();
        
        Matrix expA;
        
        if (nrmA <= 2.1) {
            // For small norm, use Padé approximation directly
            double[] b;
            if (nrmA > 0.95) {
                b = new double[]{17643225600.0, 8821612800.0, 2075673600.0, 302702400.0,
                                30270240.0, 2162160.0, 110880.0, 3960.0, 90.0, 1.0};
            } else if (nrmA > 0.25) {
                b = new double[]{17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0};
            } else if (nrmA > 0.015) {
                b = new double[]{30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0};
            } else {
                b = new double[]{120.0, 60.0, 12.0, 1.0};
            }
            
            Matrix A2 = A.mult(A);
            Matrix Q = new Matrix(n, n);
            Matrix U = Matrix.eye(n).scale(b[1]);
            Matrix V = Matrix.eye(n).scale(b[0]);
            int K = b.length / 2 - 1;
            
            for (int k = 1; k <= K; k++) {
                if (k == 1) {
                    Q = A2.copy();
                } else {
                    Q = Q.mult(A2);
                }
                U = U.add(Q.scale(b[2*k + 1]));
                V = V.add(Q.scale(b[2*k]));
            }
            
            U = A.mult(U);
            Matrix denominator = V.sub(U);
            Matrix numerator = V.add(U);
            
            // Solve denominator * expA = numerator
            expA = denominator.inv().mult(numerator);
            
        } else {
            // For large norm, use scaling and squaring
            int si = (int) Math.ceil(Math.log(nrmA / 5.4) / Math.log(2));
            if (si > 0) {
                A = A.scale(1.0 / Math.pow(2, si));
            }
            
            double[] b = new double[]{64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
                                     1187353796428800.0, 129060195264000.0, 10559470521600.0,
                                     670442572800.0, 33522128640.0, 1323241920.0,
                                     40840800.0, 960960.0, 16380.0, 182.0, 1.0};
            
            Matrix A2 = A.mult(A);
            Matrix A4 = A2.mult(A2);
            Matrix A6 = A2.mult(A4);
            Matrix I = Matrix.eye(n);
            
            // U = A * (A6*(b[13]*A6 + b[11]*A4 + b[9]*A2) + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*I)
            Matrix temp1 = A6.scale(b[13]).add(A4.scale(b[11])).add(A2.scale(b[9]));
            Matrix U = A.mult(
                A6.mult(temp1)
                .add(A6.scale(b[7]))
                .add(A4.scale(b[5]))
                .add(A2.scale(b[3]))
                .add(I.scale(b[1]))
            );
            
            // V = A6*(b[12]*A6 + b[10]*A4 + b[8]*A2) + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*I
            Matrix temp2 = A6.scale(b[12]).add(A4.scale(b[10])).add(A2.scale(b[8]));
            Matrix V = A6.mult(temp2)
                      .add(A6.scale(b[6]))
                      .add(A4.scale(b[4]))
                      .add(A2.scale(b[2]))
                      .add(I.scale(b[0]));
            
            Matrix denominator = V.sub(U);
            Matrix numerator = V.add(U);
            
            // Solve denominator * expA = numerator
            expA = denominator.inv().mult(numerator);
            
            // Square si times
            for (int t = 0; t < si; t++) {
                expA = expA.mult(expA);
            }
        }
        
        // Unbalance the result
        Matrix result = expA.copy();
        for (int j = ilo; j <= ihi; j++) {
            double scj = scale[j];
            for (int i = 0; i < n; i++) {
                if (i == j) continue;
                result.set(j, i, result.get(j, i) * scj);
                result.set(i, j, result.get(i, j) / scj);
            }
        }
        
        return result;
    }
    
    /**
     * Computes the 1-norm of the matrix (maximum absolute column sum).
     *
     * @return The 1-norm of the matrix
     */
    public double norm1() {
        double maxSum = 0;
        for (int j = 0; j < getNumCols(); j++) {
            double colSum = 0;
            for (int i = 0; i < getNumRows(); i++) {
                colSum += Math.abs(get(i, j));
            }
            maxSum = Math.max(maxSum, colSum);
        }
        return maxSum;
    }

    // ================================================================================
    // ================================================================================
    // SECTION 4: ELEMENT-WISE OPERATIONS
    // ================================================================================
    // Element-by-element operations (Hadamard product, etc.)
    // ================================================================================
    // Methods for element-wise operations: multiplication, division, power, etc.
    
    /**
     * Performs element-wise division of this matrix by another matrix.
     *
     * @param B The divisor matrix (must match dimensions).
     * @return A new matrix with the result of element-wise division A ./ B.
     * @throws IllegalArgumentException if matrix dimensions do not match.
     */
    public Matrix elementDiv(Matrix B) {
        if (this.getNumRows() != B.getNumRows() || this.getNumCols() != B.getNumCols()) {
            throw new IllegalArgumentException("Matrix dimensions should match for element division!");
        }
        Matrix res = new Matrix(this.getNumRows(), this.getNumCols());
        for (int i = 0; i < this.getNumRows(); i++) {
            for (int j = 0; j < this.getNumCols(); j++) {
                res.set(i, j, this.get(i, j) / B.get(i, j));
            }
        }
        return res;
    }

    /**
     * Performs element-wise division with another matrix.
     *
     * @param b the matrix to divide by
     * @return the result of element-wise division
     * @throws RuntimeException if matrix dimensions do not match
     */
    public Matrix elementDivide(Matrix b) {
        if (getNumCols() != b.getNumCols() || getNumRows() != b.getNumRows()) {
            throw new RuntimeException("Element divide function requires two matrices of the same size");
        }
        Matrix result = new Matrix(getNumRows(), getNumCols(), getNumRows() * getNumCols());
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                result.set(i, j, get(i, j) / b.get(i, j));
            }
        }
        return result;
    }

    /**
     * Increases each element of the matrix by a constant value.
     *
     * @param val the value to add to each element
     * @return a new matrix with incremented values
     */
    public Matrix elementIncrease(double val) {
        Matrix res = this.copy();
        for (int row = 0; row < this.getNumRows(); row++) {
            for (int col = 0; col < this.getNumCols(); col++) {
                res.set(row, col, res.get(row, col) + val);
            }
        }
        return res;
    }

    /**
     * Returns the maximum value among all elements in the matrix.
     *
     * @return maximum element value
     */
    public double elementMax() {
        return delegate.elementMax();
    }

    /**
     * Returns the maximum absolute value among all elements in the matrix.
     *
     * @return maximum absolute element value
     */
    public double elementMaxAbs() {
        return delegate.elementMaxAbs();
    }

    /**
     * Returns the minimum value among all elements in the matrix.
     *
     * @return minimum element value
     */
    public double elementMin() {
        return delegate.elementMin();
    }

    private DMatrix elementMult(DMatrix A, DMatrix B, DMatrix output) {
        return delegate.elementMult(A, B, output);
    }

    /**
     * Performs in-place element-wise multiplication with another matrix.
     *
     * @param B the matrix to multiply with
     * @return the updated matrix (this)
     */
    public Matrix elementMult(Matrix B) {
        Matrix output = Matrix.createLike(B);
        this.elementMult(B, output);
        this.setTo(output);
        return this;
    }

    /**
     * Performs element-wise multiplication with another matrix, storing the result in the given output matrix.
     *
     * @param B      the matrix to multiply with
     * @param output the matrix to store the result in (if null, a new matrix is created)
     * @return the result of the element-wise multiplication
     */
    public Matrix elementMult(Matrix B, Matrix output) {
        DMatrix m;
        if (output == null) {
            m = elementMult(this.getData(), B.getData(), null);
        } else {
            m = elementMult(this.getData(), B.getData(), output.getData());
        }
        return new Matrix(m);
    }

    /**
     * Computes the product of the elements of a row or column vector.
     *
     * @return the product of all elements
     * @throws IllegalArgumentException if the matrix is not a row or column vector
     */
    public double elementMult() {
        if (this.getNumRows() != 1 && this.getNumCols() != 1) {
            throw new IllegalArgumentException("Argument should be a row/column vector");
        }
        if (this.getNonZeros() < this.getNumRows() * this.getNumCols()) {
            return 0;
        }
        double product = 1;
        for (int i = 0; i < this.getNumRows(); i++) {
            for (int j = 0; j < this.getNumCols(); j++) {
                product *= this.get(i, j);
            }
        }
        return product;
    }

    /**
     * Performs element-wise multiplication between this matrix and a row vector.
     * Each row of the matrix is scaled by the corresponding value in the vector.
     *
     * @param B the row vector to multiply with
     * @return the result of the element-wise multiplication
     */
    public Matrix elementMultWithVector(Matrix B) {
        if (this.getNumCols() != B.getNumCols()) {
            throw new IllegalArgumentException(
                    "Matrix dimensions should match for element multiplication!");
        }
        Matrix res = new Matrix(this.getNumRows(), this.getNumCols());
        for (int i = 0; i < this.getNumRows(); i++) {
            for (int j = 0; j < this.getNumCols(); j++) {
                res.set(i, j, this.get(i, j) * B.get(i));
            }
        }
        return res;
    }

    /**
     * Raises each non-zero element of the matrix to the specified power.
     *
     * @param a the exponent
     * @return a new matrix with powered elements
     */
    public Matrix elementPow(double a) {
        Matrix b = this.copy();
        for (int i = 0; i < this.getNonZeros(); i++) {
            b.getNonZeroValues()[i] = FastMath.pow(getNonZeroValues()[i], a);
        }
        return b;
    }

    /**
     * Raises each element of the matrix to the given power.
     *
     * @param t the exponent
     * @return a new matrix with each element raised to the power of {@code t}
     */
    public Matrix elementPower(double t) {
        Matrix res = this.copy();
        if (t == 0) {
            SparseMatrix.fillMatrixStatic(res.getData(), 1);
        } else if (t != 1) {
            for (int colIdx = 0; colIdx < this.getNumCols(); colIdx++) {
                int col1 = getColumnIndex(colIdx);
                int col2 = getColumnIndex(colIdx + 1);
                for (int i = col1; i < col2; i++) {
                    int rowIdx = this.getNonZeroRow(i);
                    res.set(rowIdx, colIdx, FastMath.pow(res.get(rowIdx, colIdx), t));
                }
            }
        }
        return res;
    }

    /**
     * Computes the sum of all elements in the matrix.
     *
     * @return the sum of all elements
     */
    public double elementSum() {
        return delegate.elementSum();
    }

    /**
     * Checks for matrix equality.
     *
     * @param obj the object to compare with
     * @return true if matrices are equal, false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof Matrix)) {
            return false;
        }
        return this.isEqualTo((Matrix) obj);
    }

    /**
     * Applies the exponential function to each element of the matrix.
     *
     * @return a new matrix with exponentiated elements
     */
    public Matrix exp() {
        Matrix ret = this.copy();
        for (int i = 0; i < ret.getNumRows(); i++) {
            for (int j = 0; j < ret.getNumCols(); j++) {
                ret.set(i, j, FastMath.exp(ret.get(i, j)));
            }
        }
        return ret;
    }

    /**
     * Expands the matrix dimensions to the specified size while preserving data.
     *
     * @param rows      new row count
     * @param cols      new column count
     * @param nz_length estimated non-zero count
     */
    public void expandMatrix(int rows, int cols, int nz_length) {
        if (rows < this.getNumRows() || cols < this.getNumCols()) {
            return;
        }

        DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(rows, cols, nz_length);
       // Handle case where matrix has no columns
        if (this.getNumCols() == 0) {
            getData().setTo(DConvertMatrixStruct.convert(triplet, (DMatrixSparseCSC) null));
            return;
        }
       for (int colIdx = 0; colIdx < this.getNumCols(); colIdx++) {
            int col1 = getColumnIndex(colIdx);
            int col2 = getColumnIndex(colIdx + 1);

            for (int i = col1; i < col2; i++) {
                int rowIdx = this.getNonZeroRow(i);
                double value = this.getNonZeroValue(i);
                triplet.addItem(rowIdx, colIdx, value);
            }
        }
        getData().setTo(DConvertMatrixStruct.convert(triplet, (DMatrixSparseCSC) null));
    }

    /**
     * Expands the matrix to be square, padding with zeros as needed.
     */
    public void expandMatrixToSquare() {
        if (this.getNumRows() == this.getNumCols()) {
            return;
        }
        Matrix expandedMatrix =
                new Matrix(
                        FastMath.max(this.getNumCols(), this.getNumRows()),
                        FastMath.max(this.getNumCols(), this.getNumRows()),
                        0);
        expandedMatrix.zero();
        for (int row = 0; row < this.getNumRows(); row++) {
            for (int col = 0; col < this.getNumCols(); col++) {
                expandedMatrix.set(row, col, this.get(row, col));
            }
        }
        this.setData(expandedMatrix.getData());
    }

    private void extractMatrix(DMatrix src, int srcX0, int srcX1, int srcY0, int srcY1, DMatrix dst, int dstY0, int dstX0) {
        delegate.extractMatrix(src, srcX0, srcX1, srcY0, srcY1, dst, dstY0, dstX0);
    }

    /**
     * Computes the factorial of each element in the matrix.
     * Equivalent to {@code exp(factln())}.
     *
     * @return a matrix with the factorial of each element
     */
    public Matrix fact() {
        return this.factln().exp();
    }

    /**
     * Computes the natural logarithm of the factorial for each element.
     *
     * @return a matrix where each element is the log-factorial of the original element
     */
    public Matrix factln() {
        Matrix ret = this.copy();
        for (int i = 0; i < ret.getNumRows(); i++) {
            for (int j = 0; j < ret.getNumCols(); j++) {
                ret.set(i, j, Maths.factln(ret.get(i, j)));
            }
        }
        return ret;
    }

    /**
     * Fills all entries in the matrix with the given value.
     *
     * @param val the value to fill with
     * @return this matrix after filling
     */
    public Matrix fill(double val) {
        fillMatrix(val);
        return this;
    }

    private void fillMatrix(double value) {
        delegate.fillMatrix(value);
    }

    // ================================================================================
    // ================================================================================
    // SECTION 11: FINDING AND FILTERING
    // ================================================================================
    // Methods for finding elements, rows, and applying filters
    // ================================================================================
    // Methods for finding elements and filtering matrices
    
    /**
     * Returns the linear indices of all non-zero elements in the matrix.
     *
     * @return a column vector containing indices of non-zero entries
     */
    public Matrix find() {
        Matrix res = new Matrix(this.getNonZeros(), 1, this.getNonZeros());
        int count = 0;
        for (int colIdx = 0; colIdx < this.getNumCols(); colIdx++) {
            int col1 = getColumnIndex(colIdx);
            int col2 = getColumnIndex(colIdx + 1);

            for (int i = col1; i < col2; i++) {
                int rowIdx = this.getNonZeroRow(i);
                res.set(count++, 0, colIdx * this.getNumRows() + rowIdx);
            }
        }
        return res;
    }

    /**
     * Returns the linear indices of all elements that are non-negative (≥ 0).
     *
     * @return a column vector of indices for non-negative elements
     */
    public Matrix findNonNegative() {
        BaseMatrix tempRes = delegate.findNonNegative();
        if (tempRes instanceof SparseMatrix) {
            return new Matrix(((SparseMatrix) tempRes).getData());
        } else {
            throw new UnsupportedOperationException("Non-sparse matrix results not yet supported");
        }
    }

    /**
     * Finds all row indices in a given column where the value is non-zero.
     *
     * @param column the column index to check
     * @return a list of row indices with non-zero values in the given column
     */
    public List<Integer> findNonZeroRowsInColumn(int column) {
        List<Integer> nonZeroRows = new ArrayList<>();
        for (int i = 0; i < this.getNumRows(); i++) {
            if (this.get(i, column) != 0) {
                nonZeroRows.add(i);
            }
        }
        return nonZeroRows;
    }

    /**
     * Finds all linear indices where the matrix has a specific value.
     *
     * @param number the value to search for
     * @return a matrix containing the linear indices where the value matches
     */
    public Matrix findNumber(double number) {
        List<Integer> array = new ArrayList<Integer>();
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        for (int colIdx = 0; colIdx < numCols; colIdx++) {
            for (int rowIdx = 0; rowIdx < numRows; rowIdx++) {
                if (this.get(rowIdx, colIdx) == number)
                    array.add(colIdx * numRows + rowIdx);
            }
        }

        int size = array.size();
        Matrix res = new Matrix(size, 1, size);
        for (int i = 0; i < size; i++)
            res.set(i, 0, array.get(i));
        return res;
    }

    /**
     * Finds all linear indices where the matrix value is zero.
     *
     * @return a matrix containing the linear indices of zero-valued elements
     */
    public Matrix findZero() {
        List<Integer> array = new ArrayList<Integer>();
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        for (int colIdx = 0; colIdx < numCols; colIdx++) {
            for (int rowIdx = 0; rowIdx < numRows; rowIdx++) {
                if (this.get(rowIdx, colIdx) == 0)
                    array.add(colIdx * numRows + rowIdx);
            }
        }

        int size = array.size();
        Matrix res = new Matrix(size, 1, size);
        for (int i = 0; i < size; i++)
            res.set(i, 0, array.get(i));
        return res;
    }

    private String formatInt(int value) {
        if (value == Integer.MIN_VALUE) return "NaN";
        if (value == Integer.MAX_VALUE) return "Inf";
        if (value == Integer.MIN_VALUE + 1) return "-Inf";
        return Integer.toString(value);
    }

    private String formatValue(double value) {
        if (value == 0.0) return "0";
        if (Double.isNaN(value)) return "NaN";
        if (Double.isInfinite(value)) return value > 0 ? "Inf" : "-Inf";
        return Double.toString(Double.parseDouble(String.format("%.5g", value)));
    }

    /**
     * Populates the matrix from a 2D integer array.
     *
     * @param matrix 2D array of integers
     * @return this matrix after population
     */
    public Matrix fromArray2D(int[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                set(i, j, matrix[i][j]);
            }
        }
        return this;
    }

    /**
     * Populates the matrix from a 2D double array.
     *
     * @param matrix 2D array of doubles
     * @return this matrix after population
     */
    public Matrix fromArray2D(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                set(i, j, matrix[i][j]);
            }
        }
        return this;
    }

    /**
     * Returns the value at the specified row and column.
     *
     * @param i row index
     * @param j column index
     * @return value at (i, j)
     */
    public double get(int i, int j) {
        return getData().get(i, j);
    }

    /**
     * Returns the value at the specified index.
     * For vectors (either row or column), returns the idx-th element.
     * For matrices, returns the value at the flattened index (column-major).
     *
     * @param idx absolute index
     * @return value at index
     */
    public double get(int idx) {
        // Check if this is a vector (either row or column)
        if (this.getNumRows() == 1 || this.getNumCols() == 1) {
            // This is a vector - treat idx as the element index
            int totalElements = Math.max(this.getNumRows(), this.getNumCols());
            if (idx >= totalElements) {
                throw new RuntimeException("Index out of matrix");
            }
            
            // For row vector (1 x n), access as (0, idx)
            // For column vector (m x 1), access as (idx, 0)
            if (this.getNumRows() == 1) {
                return getData().get(0, idx);
            } else {
                return getData().get(idx, 0);
            }
        } else {
            // This is a matrix - use column-major flattened indexing
            if (idx >= this.getNumCols() * this.getNumRows()) {
                throw new RuntimeException("Index out of matrix");
            }
            
            int row = idx % this.getNumRows();
            int col = idx / this.getNumRows();
            
            return getData().get(row, col);
        }
    }

    /**
     * Returns internal column index array from the sparse matrix structure.
     *
     * @return array of column indexes
     */
    public int[] getColIndexes() {
        return getColumnIndicesArray();
    }

    /**
     * Returns the row index of the maximum value in the specified column.
     *
     * @param col column index
     * @return row index of maximum value
     */
    public int getColMax(int col) {
        double imaxval = Double.MIN_VALUE;
        int imax = 0;
        for (int i = 0; i < this.getNumRows(); i++) {
            if (this.get(i, col) > imaxval) {
                imax = i;
                imaxval = this.get(i, col);
            }
        }
        return imax;
    }

    // ================================================================================
    // ================================================================================
    // SECTION 7: MATRIX SLICING AND EXTRACTION
    // ================================================================================
    // Methods for extracting submatrices, rows, columns, and diagonals
    // ================================================================================
    // Methods for extracting parts of matrices: getRow, getColumn, getSlice, etc.
    
    /**
     * Returns a column of the matrix as a new single-column matrix.
     *
     * @param j column index
     * @return column matrix
     */
    public Matrix getColumn(int j) {
        int numRows = this.getNumRows();
        double[] colData = new double[numRows];
        for (int i = 0; i < numRows; i++) {
            colData[i] = getData().get(i, j);
        }
        return new Matrix(colData);
    }

    // Additional delegation methods for Matrix functionality
    private int getColumnIndex(int col) {
        return delegate.getColumnIndex(col);
    }

    private int[] getColumnIndicesArray() {
        return delegate.getColumnIndicesArray();
    }

    /**
     * Returns a lightweight view into the specified column without copying data.
     * This is much more efficient than getColumn() for large sparse matrices as it
     * doesn't create a new Matrix object or copy any data.
     *
     * <p>Use this when you need to iterate through column elements or perform
     * operations that don't require a full Matrix object.</p>
     *
     * @param j column index
     * @return a ColumnView providing efficient access to the column
     */
    public ColumnView getColumnView(int j) {
        return new ColumnView((DMatrixSparseCSC) getData(), j);
    }

    // Access to underlying data structure
    public DMatrix getData() {
        if (delegate instanceof SparseMatrix) {
            return ((SparseMatrix) delegate).getData();
            //} else if (delegate instanceof DenseMatrix) {
            //return ((DenseMatrix) delegate).getData();
        } else {
            throw new UnsupportedOperationException("Non-sparse matrix delegate does not support getData() as DMatrix");
        }
    }

    public void setData(DMatrix newData) {
        if (delegate instanceof SparseMatrix) {
            ((SparseMatrix) delegate).setData((DMatrixSparseCSC) newData);
        } else {
            throw new UnsupportedOperationException("Non-sparse matrix delegate does not support setData() with DMatrix");
        }
    }

    /**
     * Returns an array of unique column indices containing non-zero elements.
     *
     * @return array of column indices
     */
    public int[] getNonZeroCols() {
        List<Integer> list = new ArrayList<>();
        for (int row = 0; row < getNumRows(); row++) {
            for (int col = 0; col < getNumCols(); col++) {
                if (this.get(row, col) != 0) {
                    list.add(col);
                }
            }
        }
        return list.stream().mapToInt(i -> i).toArray();
    }

    /**
     * Returns the number of non-zero values in the matrix.
     *
     * @return number of non-zero elements
     */
    public int getNonZeroLength() {
        return delegate.getNonZeroLength();
    }

    private void setNonZeroLength(int length) {
        delegate.setNonZeroLength(length);
    }

    private int getNonZeroRow(int index) {
        return delegate.getNonZeroRow(index);
    }

    /**
     * Returns array of row indices of non-zero entries.
     *
     * @return array of row indices
     */
    public int[] getNonZeroRows() {
        return delegate.getNonZeroRows();
    }

    private double getNonZeroValue(int index) {
        return delegate.getNonZeroValue(index);
    }

    /**
     * Returns array of non-zero values in the matrix.
     *
     * @return array of non-zero values
     */
    public double[] getNonZeroValues() {
        return delegate.getNonZeroValues();
    }

    /**
     * Returns the number of non-zero elements in this matrix.
     *
     * @return the number of non-zero elements
     */
    public int getNonZeros() {
        return delegate.getNonZeros();
    }

    /**
     * Returns the number of columns in this matrix.
     *
     * @return the number of columns
     */
    public int getNumCols() {
        return delegate.getNumCols();
    }

    /**
     * Returns total number of elements in the matrix.
     *
     * @return total element count
     */
    public int getNumElements() {
        return getData().getNumElements();
    }

    public int getNumNonZeros() {
        return delegate.getNumNonZeros();
    }

    /**
     * Returns the number of rows in this matrix.
     *
     * @return the number of rows
     */
    public int getNumRows() {
        return delegate.getNumRows();
    }

    /**
     * Returns the specified row as a single-row matrix.
     *
     * @param i row index
     * @return row matrix
     */
    public Matrix getRow(int i) {
        int numCols = this.getNumCols();
        Matrix row = new Matrix(1, numCols, numCols);
        for (int j = 0; j < numCols; j++) {
            row.set(0, j, getData().get(i, j));
        }
        return row;
    }

    /**
     * Returns column index of maximum element in a specified row.
     *
     * @param row the row index
     * @return column index of max value
     */
    public int getRowMax(int row) {
        double jmaxval = Double.MIN_VALUE;
        int jmax = 0;
        int numCols = this.getNumCols();
        for (int j = 0; j < numCols; j++) {
            if (this.get(row, j) > jmaxval) {
                jmax = j;
                jmaxval = this.get(row, j);
            }
        }
        return jmax;
    }

    /**
     * Returns a lightweight view into the specified row without copying data.
     * This is much more efficient than getRow() for large sparse matrices as it
     * doesn't create a new Matrix object or copy any data.
     *
     * <p>Use this when you need to iterate through row elements or perform
     * operations that don't require a full Matrix object.</p>
     *
     * @param i row index
     * @return a RowView providing efficient access to the row
     */
    public RowView getRowView(int i) {
        return new RowView((DMatrixSparseCSC) getData(), i);
    }

    /**
     * Returns a new matrix consisting of rows from the given start index to the end.
     *
     * @param row starting row index
     * @return matrix slice
     */
    public Matrix getRowsFrom(int row) {
        if (this.getNumRows() < 0 || row >= this.getNumRows()) {
            throw new IllegalArgumentException("The number of rows is out of range");
        }
        int newNumRows = this.getNumRows() - row;
        DMatrixSparseCSC newData = new DMatrixSparseCSC(newNumRows, getData().getNumCols());
        for (int r = row; r < getData().getNumRows(); r++) {
            int start = getColumnIndex(r);
            int end = getColumnIndex(r + 1);
            for (int i = start; i < end; i++) {
                int c = this.getNonZeroRow(i);
                double value = this.getNonZeroValue(i);
                newData.set(r - row, c, value);
            }
        }
        return new Matrix(newData);
    }

    /**
     * Extracts a submatrix based on row/column bounds.
     *
     * @param r0 start row
     * @param r1 end row (exclusive)
     * @param c0 start column
     * @param c1 end column (exclusive)
     * @return submatrix
     */
    public Matrix getSlice(int r0, int r1, int c0, int c1) {
        // allocate output CSC with zero initial nz; EJML will grow as needed
        DMatrixSparseCSC sliceCsc = new DMatrixSparseCSC(r1 - r0, c1 - c0, 0);
        extractMatrix(this.getData(), r0, r1, c0, c1, sliceCsc, 0, 0);
        return new Matrix(sliceCsc);
    }

    /**
     * Extracts a submatrix from rows/columns marked as true in input flags.
     *
     * @param rowFlags boolean array for rows
     * @param colFlags boolean array for columns
     * @return submatrix
     */
    public Matrix getSlice(boolean[] rowFlags, boolean[] colFlags) {
        if (rowFlags.length != this.getNumRows())
            throw new IllegalArgumentException("rowFlags length must be " + this.getNumRows());
        if (colFlags.length != this.getNumCols())
            throw new IllegalArgumentException("colFlags length must be " + this.getNumCols());

        // build mappings from old → new indices
        int[] rowMap = new int[this.getNumRows()];
        int newRowCount = 0;
        for (int i = 0; i < this.getNumRows(); i++) {
            if (rowFlags[i]) {
                rowMap[i] = newRowCount++;
            } else {
                rowMap[i] = -1;
            }
        }

        int[] colMap = new int[this.getNumCols()];
        int newColCount = 0;
        for (int j = 0; j < this.getNumCols(); j++) {
            if (colFlags[j]) {
                colMap[j] = newColCount++;
            } else {
                colMap[j] = -1;
            }
        }

        // allocate destination; let EJML grow nz as needed
        DMatrixSparseCSC dst = new DMatrixSparseCSC(newRowCount, newColCount, 0);

        // iterate only non-zeros in the source
        for (int oldCol = 0; oldCol < this.getNumCols(); oldCol++) {
            int dstCol = colMap[oldCol];
            if (dstCol < 0) continue;

            int idx0 = getColumnIndex(oldCol);
            int idx1 = getColumnIndex(oldCol + 1);
            for (int idx = idx0; idx < idx1; idx++) {
                int oldRow = this.getNonZeroRow(idx);
                int dstRow = rowMap[oldRow];
                if (dstRow < 0) continue;

                double v = this.getNonZeroValue(idx);
                dst.set(dstRow, dstCol, v);
            }
        }

        return new Matrix(dst);
    }

    /**
     * Returns a matrix consisting of rows and columns selected by two index matrices.
     *
     * @param rows matrix of row indices
     * @param cols matrix of column indices
     * @return submatrix
     */
    public Matrix getSubMatrix(Matrix rows, Matrix cols) {
        Matrix subMatrix = new Matrix(rows.getNumRows(), cols.getNumRows());
        for (int i = 0; i < rows.getNumRows(); i++) {
            for (int j = 0; j < cols.getNumRows(); j++) {
                subMatrix.set(i, j, this.get((int) (rows.get(i)), (int) (cols.get(j))));
            }
        }
        return subMatrix;
    }

    /**
     * Increases the internal column capacity of the matrix.
     *
     * @param newmax   new max columns
     * @param preserve true if data should be preserved
     */
    public void growMaxColumns(int newmax, boolean preserve) {
        ((DMatrixSparseCSC) getData()).growMaxColumns(newmax, preserve);
    }

    /**
     * Increases the internal non-zero element capacity of the matrix.
     *
     * @param newmax   new max size
     * @param preserve true if data should be preserved
     */
    public void growMaxLength(int newmax, boolean preserve) {
        ((DMatrixSparseCSC) getData()).growMaxLength(newmax, preserve);
    }

    /**
     * Performs the Hadamard (element-wise) product of two matrices.
     *
     * @param B the other matrix
     * @return matrix resulting from element-wise multiplication
     */
    public Matrix hadamard(Matrix B) {
        return elementMult(B);
    }

    public boolean hasDuplicates() {
        return delegate.hasDuplicates();
    }

    public boolean hasFinite() {
        return delegate.hasFinite();
    }

    public boolean hasInfinite() {
        return delegate.hasInfinite();
    }

    public boolean hasMultipleFinite() {
        return delegate.hasMultipleFinite();
    }

    // Additional delegation methods used by Kotlin code
    public boolean hasNaN() {
        return delegate.hasNaN();
    }

    /**
     * Generates a hash code for the matrix based on its values.
     *
     * @return hash code
     */
    @Override
    public int hashCode() {
        int result = 1;
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                result = result + Double.hashCode(get(i, j));
            }
        }
        return result;
    }

    /**
     * Inserts a sub-matrix into the current matrix at the specified location.
     *
     * @param start_row             Starting row index
     * @param start_col             Starting column index
     * @param end_row               Ending row index
     * @param end_col               Ending column index
     * @param matrix_to_be_inserted Matrix to insert
     */
    public void insertSubMatrix(
            int start_row, int start_col, int end_row, int end_col, Matrix matrix_to_be_inserted) {
        if (end_col - start_col != matrix_to_be_inserted.getNumCols()
                || end_row - start_row != matrix_to_be_inserted.getNumRows()) {
            throw new RuntimeException("matrix_to_be_inserted doesn't fit");
        }
        for (int i = 0; i < matrix_to_be_inserted.getNumRows(); i++) {
            for (int j = 0; j < matrix_to_be_inserted.getNumCols(); j++) {
                this.set(start_row + i, start_col + j, matrix_to_be_inserted.get(i, j));
            }
        }
    }

    /**
     * Computes the inverse of the matrix.
     *
     * @return the inverse matrix
     */
    public Matrix inv() {
        if (getNumRows() != getNumCols()) {
            throw new IllegalArgumentException("Matrix must be square to compute inverse");
        }

        // Check if matrix is singular
        double det = this.det();
        if (Math.abs(det) < 1e-14) {
            throw new RuntimeException("Matrix is singular and cannot be inverted (determinant = " + det + ")");
        }

        DMatrixRMaj inverse = new DMatrixRMaj(this.getNumRows(), this.getNumCols());
        Matrix thisinverse = new Matrix(this.getNumRows(), this.getNumCols());
        SparseMatrix.invertStatic(this.getData(), inverse);
        DConvertMatrixStruct.convert(inverse, thisinverse.getData());
        return thisinverse;
    }

    /**
     * Computes the Moore-Penrose pseudo-inverse of the matrix using SVD.
     * Works for singular, rank-deficient, and non-square matrices.
     * For full-rank square matrices, this is equivalent to inv().
     *
     * @return the pseudo-inverse matrix
     */
    public Matrix pinv() {
        Ret.SVD svd = this.svd();
        Matrix U = svd.u;
        Matrix S = svd.s;  // Column vector of singular values
        Matrix V = svd.v;

        int m = getNumRows();
        int n = getNumCols();
        int rank = S.getNumRows();

        // Compute tolerance for rank determination
        double tol = 1e-10 * Math.max(m, n) * (rank > 0 ? S.get(0, 0) : 1.0);

        // Create diagonal matrix S+ from singular values
        Matrix Splus = new Matrix(rank, rank);
        for (int i = 0; i < rank; i++) {
            double sigma = S.get(i, 0);
            if (Math.abs(sigma) > tol) {
                Splus.set(i, i, 1.0 / sigma);
            }
        }

        // Pseudo-inverse: A+ = V * S+ * U^T
        // Need to extract proper dimensions from U and V
        Matrix Ut = U.transpose();
        return V.mult(Splus).mult(Ut);
    }

    /**
     * Checks whether a specific element is assigned in the sparse matrix.
     *
     * @param row Row index
     * @param col Column index
     * @return true if assigned, false otherwise
     */
    public boolean isAssigned(int row, int col) {
        return ((DMatrixSparseCSC) getData()).isAssigned(row, col);
    }

    /**
     * Determines whether the matrix is a diagonal matrix.
     * Zero matrices are not considered diagonal in this implementation.
     *
     * @return true if diagonal, false otherwise
     */
    public boolean isDiag() {
        if (getNumCols() != getNumRows()) return false;

        // Check if matrix is zero matrix - zero matrices are not considered diagonal
        boolean hasNonZeroDiagonal = false;
        for (int i = 0; i < getNumRows(); i++) {
            if (Math.abs(get(i, i)) > 1e-15) {
                hasNonZeroDiagonal = true;
                break;
            }
        }
        if (!hasNonZeroDiagonal) {
            return false;
        }

        // Check all off-diagonal elements are zero
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                if (i != j && Math.abs(get(i, j)) > 1e-15) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Checks if the matrix is empty (has zero rows or columns).
     *
     * @return true if empty, false otherwise
     */
    public boolean isEmpty() {
        return (this.getNumCols() == 0 || this.getNumRows() == 0);
    }

    /**
     * Checks if two matrices are exactly equal.
     *
     * @param m the matrix to compare
     * @return true if equal, false otherwise
     */
    public boolean isEqualTo(Matrix m) {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        if (numRows != m.getNumRows() || numCols != m.getNumCols()) return false;
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                if (this.get(i, j) != m.get(i, j)) return false;
            }
        }
        return true;
    }

    /**
     * Checks if two matrices are equal within a specified tolerance.
     *
     * @param m   matrix to compare
     * @param tol tolerance value
     * @return true if equal within tolerance, false otherwise
     */
    public boolean isEqualToTol(Matrix m, double tol) {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        if (numRows != m.getNumRows() || numCols != m.getNumCols()) return false;
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                if (Math.abs(this.get(i, j) - m.get(i, j)) > tol) return false;
            }
        }
        return true;
    }

    /**
     * Checks whether all matrix elements are finite.
     *
     * @return true if all elements are finite, false otherwise
     */
    public boolean isFinite() {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                if (Utils.isInf(this.get(i, j)) || Double.isNaN(this.get(i, j))) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Checks if all values in the matrix are integers.
     *
     * @return true if all values are integers, false otherwise
     */
    public boolean isInteger() {
        boolean isInt = true;
        for (int i = 0; i < this.getNonZeroValues().length; i++) {
            if (FastMath.abs(this.getNonZeroValue(i) - FastMath.round(this.getNonZeroValue(i))) > GlobalConstants.Zero) {
                isInt = false;
                break;
            }
        }
        return isInt;
    }

    /**
     * Keeps only the specified columns in the matrix.
     * Any column not listed will be removed.
     * For better performance, prefer passing a {@link HashSet} for `cols`.
     *
     * @param cols Collection of column indices to retain
     */
    public void keepCols(Collection<Integer> cols) {
        Set<Integer> complement = new HashSet<>();
        for (int i = 0; i < this.getNumCols(); i++) {
            complement.add(i);
        }
        complement.removeAll(cols);
        removeCols(complement);
    }

    /**
     * Keeps only the specified rows in the matrix.
     * Any row not listed will be removed.
     * For better performance, prefer passing a {@link HashSet} for `rows`.
     *
     * @param rows Collection of row indices to retain
     */
    public void keepRows(Collection<Integer> rows) {
        Set<Integer> complement = new HashSet<>();
        for (int i = 0; i < this.getNumCols(); i++) {
            complement.add(i);
        }
        complement.removeAll(rows);
        removeRows(complement);
    }

    // ================================================================================
    // ================================================================================
    // SECTION 5: MATRIX TRANSFORMATIONS
    // ================================================================================
    // Kronecker product, transpose, reshape, and other transformations
    // ================================================================================
    // Methods for transforming matrices: transpose, reshape, reverse, kron, etc.
    
    /**
     * Computes the Kronecker product of this matrix and another matrix.
     *
     * @param b The matrix to compute the Kronecker product with
     * @return The Kronecker product matrix
     */
    public Matrix kron(Matrix b) {
        DMatrixRMaj A = new DMatrixRMaj(this.getData());
        DMatrixRMaj B = new DMatrixRMaj(b.getData());
        DMatrixRMaj C = DenseMatrix.kronStatic(A, B, null);
        return new Matrix(SimpleMatrix.wrap(C));
    }

    /**
     * Computes the Kronecker sum of two matrices: A ⊕ B = A \otimes I + I \otimes B,
     * where \otimes is the Kronecker product and I is the identity matrix of matching size.
     *
     * @param other The other matrix
     * @return The Kronecker sum of this matrix and the other matrix
     */
    public Matrix krons(Matrix other) {
        DMatrixRMaj A = new DMatrixRMaj(this.getData());
        DMatrixRMaj B = new DMatrixRMaj(other.getData());
        DMatrixRMaj C = DenseMatrix.kronStatic(A, DenseMatrix.identityStatic(B.numRows), null);
        DMatrixRMaj D = DenseMatrix.kronStatic(DenseMatrix.identityStatic(A.numRows), B, null);
        DMatrixRMaj output = DenseMatrix.addStatic(C, D, null);
        return new Matrix(SimpleMatrix.wrap(output));
    }

    /**
     * Solves the equation AX = B for X, where A is this matrix and B is the right-hand side.
     * For square matrices, this solves the exact system. For rectangular matrices (overdetermined),
     * this computes the least squares solution using the normal equation: X = (A^T * A)^-1 * A^T * B.
     *
     * @param b The right-hand side matrix
     * @return The solution matrix X
     */
    public Matrix leftMatrixDivide(Matrix b) {
        if (getNumRows() == getNumCols()) {
            // Square matrix: try exact solve first, fall back to pseudo-inverse
            try {
                Matrix x = new Matrix(0, 0, 0);
                solve(this, b, x);
                return x;
            } catch (RuntimeException e) {
                // If singular, use pseudo-inverse via SVD
                if (e.getMessage() != null && e.getMessage().contains("singular")) {
                    return pseudoInverseSolve(this, b);
                }
                throw e;
            }
        } else if (getNumRows() > getNumCols()) {
            // Overdetermined system: least squares via pseudo-inverse
            return pseudoInverseSolve(this, b);
        } else {
            // Underdetermined system: minimum norm via pseudo-inverse
            return pseudoInverseSolve(this, b);
        }
    }

    /**
     * Solves A·X = B using pseudo-inverse via SVD.
     * Robust for singular, rank-deficient, and ill-conditioned matrices.
     *
     * @param A Coefficient matrix
     * @param B Right-hand side matrix
     * @return Solution matrix X
     */
    private static Matrix pseudoInverseSolve(Matrix A, Matrix B) {
        // Compute pseudo-inverse using SVD
        Ret.SVD svd = A.svd();
        Matrix U = svd.u;
        Matrix S = svd.s;  // Column vector of singular values
        Matrix V = svd.v;

        // Create S+ (pseudo-inverse of S)
        int m = A.getNumRows();
        int n = A.getNumCols();
        int rank = S.getNumRows();
        Matrix Splus = new Matrix(rank, rank);

        // Compute tolerance for rank determination
        double tol = 1e-10 * Math.max(m, n) * (rank > 0 ? S.get(0, 0) : 1.0);
        
        for (int i = 0; i < rank; i++) {
            double sigma = S.get(i, 0);
            if (Math.abs(sigma) > tol) {
                Splus.set(i, i, 1.0 / sigma);
            }
        }

        // A+ = V * S+ * U^T
        Matrix Aplus = V.mult(Splus).mult(U.transpose());
        return Aplus.mult(B);
    }

    /**
     * Returns the length of the matrix, defined as max(rows, cols).
     *
     * @return The maximum dimension of the matrix
     */
    public int length() {
        return FastMath.max(getNumRows(), getNumCols());
    }

    /**
     * Applies the natural logarithm element-wise.
     *
     * @return A new matrix with log applied to each element
     */
    public Matrix log() {
        Matrix ret = this.copy();
        for (int i = 0; i < ret.getNumRows(); i++) {
            for (int j = 0; j < ret.getNumCols(); j++) {
                ret.set(i, j, FastMath.log(ret.get(i, j)));
            }
        }
        return ret;
    }

    // ================================================================================
    // ================================================================================
    // SECTION 9: STATISTICAL OPERATIONS
    // ================================================================================
    // Mean, sum, variance, norms, and other statistical measures
    // ================================================================================
    // Methods for statistical operations: mean, sum, standard deviation, etc.
    
    /**
     * Computes the mean of each column.
     *
     * @return A 1xN matrix containing the mean of each column
     */
    public Matrix meanCol() {
        Matrix res = new Matrix(1, this.getNumCols());
        for (int col = 0; col < this.getNumCols(); col++) {
            res.set(0, col, this.sumCols(col) / this.getNumRows());
        }
        return res;
    }

    /**
     * Computes the mean of each row.
     *
     * @return An Nx1 matrix containing the mean of each row
     */
    public Matrix meanRow() {
        Matrix res = new Matrix(this.getNumRows(), 1);
        for (int row = 0; row < this.getNumRows(); row++) {
            res.set(row, 0, this.sumRows(row) / this.getNumCols());
        }
        return res;
    }

    /**
     * Multiplies all elements in the matrix by -1 in-place.
     */
    public void mulByMinusOne() {
        for (int row = 0; row < this.getData().getNumRows(); row++) {
            for (int col = 0; col < this.getData().getNumCols(); col++) {
                this.getData().set(row, col, -this.getData().get(row, col));
            }
        }
    }

    /**
     * Performs matrix multiplication: this * B
     *
     * @param B The right-hand side matrix
     * @return Resultant matrix
     */
    public Matrix mult(Matrix B) {
        return new Matrix(multMatrix(B.getData(), null));
    }

    /**
     * Performs matrix multiplication: this * B
     *
     * @param B   The right-hand side matrix
     * @param out Output matrix to write result to (can be null)
     * @return Resultant matrix
     */
    public Matrix mult(Matrix B, Matrix out) {
        if (out == null) {
            return new Matrix(multMatrix(B.getData(), null));
        } else {
            return new Matrix(multMatrix(B.getData(), out.getData()));
        }
    }

    /**
     * Efficiently multiplies this row vector with a column view.
     * This is optimized for the common pattern: rowVector.mult(matrix.getColumn(j))
     * by avoiding the creation of intermediate Matrix objects.
     *
     * @param columnView the column view to multiply with
     * @return the scalar result of the multiplication
     * @throws IllegalArgumentException if this matrix is not a row vector or dimensions don't match
     */
    public double multColumnView(ColumnView columnView) {
        if (getNumRows() != 1) {
            throw new IllegalArgumentException("This matrix must be a row vector (1 row)");
        }
        if (getNumCols() != columnView.getNumRows()) {
            throw new IllegalArgumentException("Dimension mismatch: row vector has " + getNumCols() +
                    " columns but column has " + columnView.getNumRows() + " rows");
        }

        double result = 0.0;

        // Iterate through non-zero elements in the column
        for (int i = 0; i < columnView.getNonZeroCount(); i++) {
            int row = columnView.getNonZeroRow(i);
            double colValue = columnView.getNonZeroValue(i);
            double rowValue = this.get(0, row); // Get value from row vector
            result += rowValue * colValue;
        }

        return result;
    }

    /**
     * Replaces this matrix with the result of this * B.
     *
     * @param B The right-hand side matrix
     */
    public void multEq(Matrix B) {
        Matrix output = new Matrix(multMatrix(B.getData(), null));
        this.getData().setTo(output.getData());
    }

    private DMatrix multMatrix(DMatrix B, DMatrix output) {
        return delegate.multMatrix(B, output);
    }

    /**
     * Efficiently multiplies a row view with this column vector.
     * This is optimized for the common pattern: matrix.getRow(i).mult(columnVector)
     * by avoiding the creation of intermediate Matrix objects.
     *
     * @param rowView the row view to multiply with
     * @return the scalar result of the multiplication
     * @throws IllegalArgumentException if this matrix is not a column vector or dimensions don't match
     */
    public double multRowView(RowView rowView) {
        if (getNumCols() != 1) {
            throw new IllegalArgumentException("This matrix must be a column vector (1 column)");
        }
        if (getNumRows() != rowView.getNumCols()) {
            throw new IllegalArgumentException("Dimension mismatch: column vector has " + getNumRows() +
                    " rows but row has " + rowView.getNumCols() + " columns");
        }

        return rowView.dotProduct(this);
    }

    /**
     * Computes the Euclidean norm of a matrix
     *
     * @return - the Euclidean norm of the given matrix
     */
    public double norm() {
        double sum = 0;
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                double num = this.get(i, j);

                // Handle special values
                if (Double.isNaN(num)) {
                    return Double.NaN;
                }
                if (Double.isInfinite(num)) {
                    return Inf;
                }

                sum += num * num;
            }
        }
        return FastMath.sqrt(sum);
    }

    /**
     * Computes the Frobenius norm of a matrix (alias for norm()).
     *
     * @return - the Frobenius norm of the given matrix
     */
    public double normFrobenius() {
        return norm();
    }

    /**
     * Fills the matrix with ones. Matrix must already be initialized with correct size.
     * Throws an exception if matrix cannot be fully filled.
     */
    public void ones() {
        int rows = this.getNumRows();
        int cols = this.getNumCols();

        if (this.getData().getNumElements() != rows * cols) {
            throw new RuntimeException("Matrix is too small to fill with ones");
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                this.set(i, j, 1);
            }
        }
    }

    /**
     * Computes the sum of powers of absolute values in a column: ∑ |A_{ji}|^alpha
     *
     * @param col   The column index
     * @param alpha The exponent to raise each absolute value
     * @return The computed sum
     */
    public double powerSumCols(int col, double alpha) {
        double sum = 0;
        int numRows = this.getNumRows();
        for (int i = 0; i < numRows; i++) {
            sum += FastMath.pow(Math.abs(this.get(i, col)), alpha);
        }
        return sum;
    }

    /**
     * Computes the sum of powers of absolute values in a row: ∑ |A_{ij}|^alpha
     *
     * @param row   The row index
     * @param alpha The exponent to raise each absolute value
     * @return The computed sum
     */
    public double powerSumRows(int row, double alpha) {
        double sum = 0;
        int numCols = this.getNumCols();
        for (int i = 0; i < numCols; i++) {
            sum += FastMath.pow(Math.abs(this.get(row, i)), alpha);
        }
        return sum;
    }

    /**
     * Pretty prints the matrix with automatic formatting and alignment.
     * Supports real, NaN, Inf, and -Inf values.
     */
    public void prettyPrint() {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        
        if (numRows == 0 || numCols == 0) {
            System.out.println("[]");
            return;
        }
        
        int[] columnWidths = new int[numCols];

        for (int j = 0; j < numCols; j++) {
            for (int i = 0; i < numRows; i++) {
                double value = this.get(i, j);
                String formattedValue = formatValue(value);
                columnWidths[j] = Math.max(columnWidths[j], formattedValue.length());
            }
        }

        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                String formattedValue = formatValue(this.get(i, j));
                sb.append(String.format("%" + columnWidths[j] + "s", formattedValue));
                if (j < numCols - 1) sb.append("  ");
            }
            sb.append(i < numRows - 1 ? " \n " : "]\n");
        }
        System.out.println(sb);
    }

    /**
     * Pretty prints the matrix assuming integer values.
     */
    public void prettyPrintInt() {
        int numRows = this.getNumRows();
        int numCols = this.getNumCols();
        
        if (numRows == 0 || numCols == 0) {
            System.out.println("[]");
            return;
        }
        
        int[] columnWidths = new int[numCols];

        for (int j = 0; j < numCols; j++) {
            for (int i = 0; i < numRows; i++) {
                int value = (int) Math.round(this.get(i, j));
                String formatted = formatInt(value);
                columnWidths[j] = Math.max(columnWidths[j], formatted.length());
            }
        }

        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                String formatted = formatInt((int) Math.round(this.get(i, j)));
                sb.append(String.format("%" + columnWidths[j] + "s", formatted));
                if (j < numCols - 1) sb.append("  ");
            }
            sb.append(i < numRows - 1 ? " \n " : "]\n");
        }
        System.out.println(sb);
    }

    // ================================================================================
    // ================================================================================
    // SECTION 14: I/O OPERATIONS
    // ================================================================================
    // Methods for input/output: print, toString, file operations
    // ================================================================================
    // Methods for input/output: print, toString, readFromFile, etc.
    
    /**
     * Prints the matrix with appropriate formatting, either as floats or integers.
     */
    public void print() {
        if (isInteger()) {
            prettyPrintInt();
        } else {
            prettyPrint();
        }
    }

    /**
     * Prints only non-zero values and their positions in the matrix.
     */
    public void printNonZero() {
        ((DMatrixSparseCSC) getData()).printNonZero();
    }

    /**
     * Performs QR decomposition on the matrix.
     * Only square matrices are currently supported.
     *
     * @return A map containing "Q" and "R" matrices from the decomposition.
     */
    public Map<String, Matrix> qr() {
        if (getNumRows() != getNumCols()) {
            throw new RuntimeException("Only square matrix can be decomposed");
        }
        QRSparseDecomposition<DMatrixSparseCSC> decomposition =
                new QrLeftLookingDecomposition_DSCC(null);
        decomposition.decompose(this.toDMatrixSparseCSC());

        Map<String, Matrix> result = new HashMap<>();
        result.put("Q", new Matrix(decomposition.getQ(null, false)));
        result.put("R", new Matrix(decomposition.getR(null, false)));
        return result;
    }

    /**
     * Fills the matrix with random values between 0 and 1.
     *
     * @param length The number of rows and columns in the square matrix.
     */
    public void randMatrix(int length) {
        // Use RandomManager for centralized random number generation
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                this.set(i, j, RandomManager.nextDouble());
            }
        }
    }

    /**
     * Computes the rank of the matrix.
     *
     * @return The rank of the matrix.
     */
    public int rank() {
        DMatrixRMaj matrix = new DMatrixRMaj(this.getNumRows(), this.getNumCols());
        DConvertMatrixStruct.convert(this.getData(), matrix);
        return MatrixFeatures_DDRM.rank(matrix);
    }

    /**
     * Computes the element-wise reciprocal (1/x) of the matrix.
     *
     * @return A new matrix with reciprocal values.
     */
    public Matrix reciprocal() {
        Matrix B = Matrix.createLike(this);
        for (int i = 0; i < this.getNumElements(); i++) {
            B.set(i, 1.0 / this.get(i));
        }
        return B;
    }

    /**
     * Removes a specific element from the matrix at the given row and column.
     *
     * @param row The row index
     * @param col The column index
     */
    public void remove(int row, int col) {
        ((DMatrixSparseCSC) getData()).remove(row, col);
    }

    /**
     * Removes the specified columns from the matrix.
     * It is recommended to use a HashSet for efficiency when checking column indices.
     *
     * @param cols The indices of columns to be removed.
     */
    public void removeCols(Collection<Integer> cols) {
        for (int c : cols) {
            if (c < 0 || c >= this.getNumCols()) {
                throw new IllegalArgumentException(
                        "Cannot remove cols that are outside of the current matrix");
            }
        }
        Matrix newMatrix = new Matrix(this.getNumRows(), this.getNumCols() - cols.size());
        for (int i = 0; i < this.getNumRows(); i++) {
            int factor = 0;
            for (int j = 0; j < this.getNumCols(); j++) {
                if (cols.contains(j)) {
                    factor++;
                } else {
                    newMatrix.set(i, j - factor, this.get(i, j));
                }
            }
        }
        this.getData().setTo(newMatrix.getData());
    }

    /**
     * Replaces all infinite values (positive or negative) in the matrix with 0.
     * Also adjusts internal bookkeeping fields to maintain structural integrity.
     */
    public void removeInfinite() {
        int offset = 0;
        for (int i = 0; i < this.getNumCols(); i++) {
            for (int j = 0; j < this.getNumRows(); j++) {
                if (isInf(get(j, i))) {
                    set(j, i, 0);
                }
            }
            setColumnIndex(i + 1, getColumnIndex(i + 1) - offset);
        }
        setNonZeroLength(getNonZeros() - offset);
    }

    /**
     * Removes all infinite values (positive or negative) from the sparse matrix structure.
     * This shifts remaining values to preserve compactness.
     */
    public void removeInfinity() {
        int offset = 0;
        for (int i = 0; i < this.getNumCols(); i++) {
            int idx0 = getColumnIndex(i) + offset;
            int idx1 = getColumnIndex(i + 1);

            for (int j = idx0; j < idx1; j++) {
                double val = this.getNonZeroValue(j);
                if (!Utils.isInf(val) && !Double.isNaN(val)) {
                    this.setNonZeroRow(j - offset, this.getNonZeroRow(j));
                    this.setNonZeroValue(j - offset, val);
                } else {
                    offset++;
                }
            }
            setColumnIndex(i + 1, getColumnIndex(i + 1) - offset);
        }
        setNonZeroLength(getNonZeros() - offset);
    }

    /**
     * Removes all NaN values from the matrix structure.
     * Internal arrays are compacted and updated accordingly.
     */
    public void removeNaN() {
        if (!hasNaN()) return;

        int offset = 0;
        for (int i = 0; i < this.getNumCols(); i++) {
            int idx0 = getColumnIndex(i) + offset;
            int idx1 = getColumnIndex(i + 1);

            for (int j = idx0; j < idx1; j++) {
                double val = this.getNonZeroValue(j);
                if (!Double.isNaN(val)) {
                    this.setNonZeroRow(j - offset, this.getNonZeroRow(j));
                    this.setNonZeroValue(j - offset, val);
                } else {
                    offset++;
                }
            }
            setColumnIndex(i + 1, getColumnIndex(i + 1) - offset);
        }
        setNonZeroLength(getNonZeros() - offset);
    }

    /**
     * Removes all negative values from the matrix structure.
     * Shifts non-negative values to preserve compactness.
     */
    public void removeNegative() {
        int offset = 0;
        for (int i = 0; i < this.getNumCols(); i++) {
            int idx0 = getColumnIndex(i) + offset;
            int idx1 = getColumnIndex(i + 1);

            for (int j = idx0; j < idx1; j++) {
                double val = this.getNonZeroValue(j);
                if (val > 0) {
                    this.setNonZeroRow(j - offset, this.getNonZeroRow(j));
                    this.setNonZeroValue(j - offset, val);
                } else {
                    offset++;
                }
            }
            setColumnIndex(i + 1, getColumnIndex(i + 1) - offset);
        }
        setNonZeroLength(getNonZeros() - offset);
    }

    /**
     * Removes the specified rows from the matrix. If possible, use a HashSet for better performance.
     *
     * @param rows the indices of the rows to be removed
     */
    public void removeRows(Collection<Integer> rows) {
        for (int r : rows) {
            if (r < 0 || r >= this.getNumRows()) {
                throw new IllegalArgumentException("Cannot remove rows that are outside of the current matrix");
            }
        }
        Matrix newMatrix = new Matrix(this.getNumRows() - rows.size(), this.getNumCols());
        int factor = 0;
        for (int i = 0; i < this.getNumRows(); i++) {
            if (rows.contains(i)) {
                factor++;
            } else {
                for (int j = 0; j < this.getNumCols(); j++) {
                    newMatrix.set(i - factor, j, this.get(i, j));
                }
            }
        }
        this.getData().setTo(newMatrix.getData());
    }

    // Method delegation for operations
    private void removeZeros() {
        delegate.removeZeros();
    }

    /**
     * Removes values from the matrix that are equal to the specified value.
     *
     * @param val the value to be removed (typically 0)
     */
    public void removeZeros(double val) {
        removeZerosWithTol(val);
    }

    private void removeZerosWithTol(double tolerance) {
        delegate.removeZerosWithTol(tolerance);
    }

    /**
     * Replaces all occurrences of a specific value with a new value in the matrix.
     * Uses optimized sparse iteration when replacing non-zero values.
     *
     * @param delete      the value to be replaced
     * @param replacement the value to insert in place of the deleted value
     */
    public void replace(double delete, double replacement) {
        if (delegate instanceof SparseMatrix) {
            if (delete != 0.0) {
                // For non-zero delete values, iterate only over non-zero elements (much faster)
                double[] nzValues = getNonZeroValues();
                for (int i = 0; i < nzValues.length; i++) {
                    if (nzValues[i] == delete) {
                        nzValues[i] = replacement;
                    }
                }
                // If replacement is zero, we should clean up the sparse structure
                if (replacement == 0.0) {
                    removeZeros();
                }
            } else {
                // Replacing zeros with non-zeros requires restructuring - use triplet format
                int numRows = getNumRows();
                int numCols = getNumCols();
                int[] colIndices = getColumnIndicesArray();
                int[] rowIndices = getNonZeroRows();
                double[] values = getNonZeroValues();

                // Build a set of existing (row, col) pairs for quick lookup
                Set<Long> existingEntries = new HashSet<Long>();
                for (int colIdx = 0; colIdx < numCols; colIdx++) {
                    int col1 = colIndices[colIdx];
                    int col2 = colIndices[colIdx + 1];
                    for (int i = col1; i < col2; i++) {
                        long key = ((long) rowIndices[i]) * numCols + colIdx;
                        existingEntries.add(key);
                    }
                }

                // Use triplet format to rebuild with zeros replaced
                int estimatedNz = getNonZeros() + (numRows * numCols - getNonZeros());
                DMatrixSparseTriplet triplet = new DMatrixSparseTriplet(numRows, numCols, Math.min(estimatedNz, numRows * numCols));

                // Add existing non-zero values
                for (int colIdx = 0; colIdx < numCols; colIdx++) {
                    int col1 = colIndices[colIdx];
                    int col2 = colIndices[colIdx + 1];
                    for (int i = col1; i < col2; i++) {
                        triplet.addItem(rowIndices[i], colIdx, values[i]);
                    }
                }

                // Add replacement for zeros
                for (int i = 0; i < numRows; i++) {
                    for (int j = 0; j < numCols; j++) {
                        long key = ((long) i) * numCols + j;
                        if (!existingEntries.contains(key)) {
                            triplet.addItem(i, j, replacement);
                        }
                    }
                }

                getData().setTo(DConvertMatrixStruct.convert(triplet, (DMatrixSparseCSC) null));
            }
        } else {
            // Dense matrix - use simple loop (get/set is O(1) for dense)
            for (int i = 0; i < this.getNumRows(); i++) {
                for (int j = 0; j < this.getNumCols(); j++) {
                    if (this.get(i, j) == delete) {
                        this.set(i, j, replacement);
                    }
                }
            }
        }
    }

    // ================================================================================
    // ================================================================================
    // SECTION 12: MATRIX MANIPULATION
    // ================================================================================
    // Methods for modifying matrix structure: repmat, remove, replace
    // ================================================================================
    // Methods for manipulating matrix structure: repmat, reshape, setRow, setColumn, etc.
    
    /**
     * Repeats the matrix to match the specified row and column replication.
     *
     * @param rows number of times to repeat along the row dimension
     * @param cols number of times to repeat along the column dimension
     * @return a new matrix formed by repeating this matrix
     */
    public Matrix repmat(int rows, int cols) {
        // Handle edge cases where rows or cols is 0
        if (rows == 0 || cols == 0) {
            return new Matrix(rows == 0 ? 0 : this.getNumRows() * rows,
                             cols == 0 ? 0 : this.getNumCols() * cols);
        }
        Matrix res = this.copy();
        for (int i = 1; i < rows; i++) {
            Matrix tmp = new Matrix(res.getNumRows() + this.getNumRows(), res.getNumCols());
            concatRowsInPlace(res.getData(), this.getData(), tmp.getData());
            res = tmp;
        }
        Matrix colBase = res.copy();
        for (int i = 1; i < cols; i++) {
            Matrix tmp = new Matrix(res.getNumRows(), res.getNumCols() + colBase.getNumCols());
            concatColumnsInPlace(res.getData(), colBase.getData(), tmp.getData());
            res = tmp;
        }
        return res;
    }

    /**
     * Reshapes the matrix to the specified number of rows and columns.
     *
     * @param numRows number of rows
     * @param numCols number of columns
     */
    public void reshape(int numRows, int numCols) {
        if (numRows < 0 || numCols < 0) {
            throw new IllegalArgumentException("Cannot reshape matrix: dimensions must be non-negative, got (" +
                    numRows + "x" + numCols + ")");
        }
        int totalElements = this.getNumRows() * this.getNumCols();
        if (numRows * numCols != totalElements) {
            throw new IllegalArgumentException("Cannot reshape matrix: new dimensions (" +
                    numRows + "x" + numCols + ") do not match total number of elements (" +
                    totalElements + ")");
        }
        ((DMatrixSparseCSC) getData()).reshape(numRows, numCols);
    }

    /**
     * Reshapes the matrix with new dimensions and internal storage length.
     *
     * @param numRows     number of rows
     * @param numCols     number of columns
     * @param arrayLength number of non-zero entries the matrix can store
     */
    public void reshape(int numRows, int numCols, int arrayLength) {
        ((DMatrixSparseCSC) getData()).reshape(numRows, numCols, arrayLength);
    }

    /**
     * Reverses the elements of a vector (row or column) matrix.
     * Equivalent to MATLAB's v(end:-1:1).
     *
     * @return a matrix with reversed elements
     */
    public Matrix reverse() {
        Matrix ret = Matrix.createLike(this);
        for (int i = 0; i < ret.getNumElements(); i++) {
            ret.set(i, this.get(ret.getNumElements() - 1 - i));
        }
        return ret;
    }

    /**
     * Reverses the order of rows in the matrix.
     *
     * @return a new matrix with rows in reverse order
     */
    public Matrix reverseRows() {
        int numRows = getNumRows();
        int numCols = getNumCols();
        Matrix ret = Matrix.createLike(this);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                ret.set(i, j, this.get(numRows - 1 - i, j));
            }
        }
        return ret;
    }

    /**
     * Performs right matrix division A / B = A * inv(B)
     *
     * @param b the matrix to divide by (on the right)
     * @return result of A / B
     */
    public Matrix rightMatrixDivide(Matrix b) {
        // Solve X * B = A (MATLAB: X = A / B)
        // Transpose: B^T * X^T = A^T
        // Solve for X^T: X^T = (B^T) \ (A^T)
        // Return X: X = ((B^T) \ (A^T))^T
        return b.transpose().leftMatrixDivide(this.transpose()).transpose();
    }

    /**
     * Increases all elements in the specified row by a scalar value.
     *
     * @param row the row to modify
     * @param a   the scalar to add to each element in the row
     */
    public void rowIncrease(int row, double a) {
        if (row >= 0 && row < getNumRows()) {
            for (int i = 0; i < getNumCols(); i++) {
                set(row, i, get(row, i) + a);
            }
        } else {
            throw new RuntimeException("Row index out of range");
        }
    }

    // REMOVE
    public Matrix safeMult(Matrix B) {
        if (B.length() == 1 && this.length() != 1) {
            return Matrix.scaleMult(this, B.get(0));
        } else if (this.getNumCols() == B.getNumRows()) {
            return this.mult(B);
        } else {
            throw new RuntimeException("matrix product of X and Y failed");
        }
    }

    /**
     * Scales the matrix by the given scalar value.
     *
     * @param scalar the value to scale each element by
     * @return a new matrix scaled by the scalar
     */
    public Matrix scale(double scalar) {
        Matrix output = new Matrix(this);
        scaleMatrix(scalar, output.getData());
        return output;
    }

    /**
     * Returns the negation of this matrix (all elements multiplied by -1).
     *
     * @return a new matrix with all elements negated
     */
    public Matrix neg() {
        return scale(-1.0);
    }

    /**
     * Scales this matrix in-place by the given scalar value.
     *
     * @param scalar the value to scale each element by
     */
    public void scaleEq(double scalar) {
        Matrix output = new Matrix(this);
        scaleMatrix(scalar, output.getData());
        this.getData().setTo(output.getData());
    }

    /**
     * Scales this matrix by the given scalar value and stores the result in the provided output matrix.
     *
     * @param scalar the scalar multiplier
     * @param output the matrix to store the result
     */
    public void scaleEq(double scalar, Matrix output) {
        scaleMatrix(scalar, output.getData());
    }

    private void scaleMatrix(double scalar, DMatrix output) {
        delegate.scaleMatrix(scalar, output);
    }

    // ================================================================================
    // ================================================================================
    // SECTION 13: SPECIAL MATRIX OPERATIONS
    // ================================================================================
    // Specialized operations: Schur complement, factorial, reciprocal
    // ================================================================================
    // Methods for specialized matrix operations: Schur decomposition, etc.
    
    /**
     * Computes the Schur decomposition of the matrix.
     *
     * @param method method name, currently only "default" is supported
     * @param iter   number of iterations for the default method (null for default=1)
     * @return a map containing the keys "T" (upper triangular matrix) and "U" (unitary matrix)
     */
    public Map<String, Matrix> schur(String method, Integer iter) {
        if (getNumRows() != getNumCols()) {
            throw new RuntimeException("Only square matrix can be decomposed");
        }
        if (this.getNonZeroLength() == 0) {
            Map<String, Matrix> result = new HashMap<>();
            result.put("U", Matrix.eye(getNumRows()));
            result.put("T", this.copy());
            return result;
        }
        int n = getNumRows();
        Map<String, Matrix> result = new HashMap<>();

        // Use Apache Commons Math3 SchurTransformer for numerical accuracy.
        // SchurTransformer is package-private, so we access it via reflection.
        try {
            RealMatrix rm = MatrixUtils.createRealMatrix(this.toArray2D());
            Class<?> stClass = Class.forName("org.apache.commons.math3.linear.SchurTransformer");
            java.lang.reflect.Constructor<?> ctor = stClass.getDeclaredConstructor(RealMatrix.class);
            ctor.setAccessible(true);
            Object st = ctor.newInstance(rm);

            java.lang.reflect.Method getT = stClass.getMethod("getT");
            java.lang.reflect.Method getP = stClass.getMethod("getP");

            RealMatrix tMat = (RealMatrix) getT.invoke(st);
            RealMatrix pMat = (RealMatrix) getP.invoke(st);

            Matrix T = new Matrix(n, n);
            Matrix U = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double tv = tMat.getEntry(i, j);
                    if (tv != 0.0) T.set(i, j, tv);
                    double uv = pMat.getEntry(i, j);
                    if (uv != 0.0) U.set(i, j, uv);
                }
            }
            result.put("T", T);
            result.put("U", U);
        } catch (Exception e) {
            // Fallback to QR iteration if reflection fails
            int it = (iter != null) ? iter : 200 * n;
            Matrix A = copy();
            Matrix U = Matrix.eye(n);
            for (int i = 0; i < it; i++) {
                Map<String, Matrix> qrResult = A.qr();
                Matrix Q = qrResult.get("Q");
                Matrix R = qrResult.get("R");
                A = R.mult(Q);
                U = U.mult(Q);
            }
            result.put("T", A);
            result.put("U", U);
        }
        return result;
    }

    /**
     * Computes the Schur decomposition with the default method and iteration count.
     *
     * @return a map with matrices "T" and "U" from the Schur decomposition
     */
    public Map<String, Matrix> schur() {
        return schur("default", null);
    }

    /**
     * Computes the complex Schur decomposition of this matrix, equivalent to MATLAB's
     * {@code schur(A,'complex')}. The result has a truly upper-triangular T matrix with
     * complex eigenvalues on the diagonal.
     *
     * <p>The method first computes the real Schur decomposition (which may have 2x2 blocks
     * on the diagonal for complex eigenvalue pairs), then converts each 2x2 block to
     * upper-triangular form using unitary Givens rotations.</p>
     *
     * @return a map with keys "U" (ComplexMatrix, unitary) and "T" (ComplexMatrix, upper triangular)
     */
    public Map<String, ComplexMatrix> schurComplex() {
        Map<String, Matrix> realSchur = this.schur();
        Matrix Tr = realSchur.get("T");
        Matrix Ur = realSchur.get("U");
        int n = Tr.getNumRows();

        // Start with the real Schur form as complex matrices
        ComplexMatrix T = new ComplexMatrix(Tr.copy(), new Matrix(n, n));
        ComplexMatrix U = new ComplexMatrix(Ur.copy(), new Matrix(n, n));

        // Scan for 2x2 blocks on the diagonal and diagonalize them
        int k = 0;
        while (k < n - 1) {
            double subdiag = Tr.get(k + 1, k);
            if (Math.abs(subdiag) > 1e-14) {
                // Found a 2x2 block at rows/cols [k, k+1] with complex eigenvalue pair
                double a = Tr.get(k, k);
                double b = Tr.get(k, k + 1);
                double c = Tr.get(k + 1, k);
                double d = Tr.get(k + 1, k + 1);

                // Eigenvalues of [a b; c d] are (a+d)/2 +/- sqrt(disc)/2
                // where disc = (a-d)^2 + 4*b*c  (should be < 0 for complex pair)
                double trace = a + d;
                double disc = (a - d) * (a - d) + 4.0 * b * c;

                // Complex eigenvalue: lambda = trace/2 + i*sqrt(-disc)/2
                double realPart = trace / 2.0;
                double imagPart = Math.sqrt(Math.abs(disc)) / 2.0;

                // We need a unitary 2x2 matrix G such that G^H * [a b; c d] * G = upper triangular
                // with eigenvalue lambda1 = realPart + i*imagPart on (0,0)
                // and lambda2 = realPart - i*imagPart on (1,1).
                //
                // Eigenvector for lambda1 = realPart + i*imagPart:
                //   (a - lambda1)*v1 + b*v2 = 0
                //   v2/v1 = -(a - lambda1)/b = -(a - realPart - i*imagPart)/b
                // Let v = [b, lambda1 - a]^T = [b, (realPart - a) + i*imagPart]^T
                // Normalize: G's first column = v / ||v||
                double vRe0 = b;
                double vIm0 = 0.0;
                double vRe1 = realPart - a;
                double vIm1 = imagPart;
                double norm = Math.sqrt(vRe0 * vRe0 + vIm0 * vIm0 + vRe1 * vRe1 + vIm1 * vIm1);
                vRe0 /= norm;
                vIm0 /= norm;
                vRe1 /= norm;
                vIm1 /= norm;

                // G = [v, w] where w is orthogonal to v (conjugate of rotated v)
                // w = [-conj(v2), conj(v1)]^T
                double wRe0 = -vRe1;
                double wIm0 = vIm1;
                double wRe1 = vRe0;
                double wIm1 = -vIm0;  // conj(v1) but v1 is real so this is just v1

                // Build the full n x n Givens rotation: identity except at rows/cols k, k+1
                // G_full[k,k] = v1, G_full[k,k+1] = w1, G_full[k+1,k] = v2, G_full[k+1,k+1] = w2
                // T <- G^H * T * G, U <- U * G

                // Apply T <- G^H * T (affects rows k and k+1)
                for (int j = 0; j < n; j++) {
                    org.apache.commons.math3.complex.Complex t_kj = T.get(k, j);
                    org.apache.commons.math3.complex.Complex t_k1j = T.get(k + 1, j);
                    // G^H row 0 = conj(v)^T: [conj(v0), conj(v1)]
                    // G^H row 1 = conj(w)^T: [conj(w0), conj(w1)]
                    org.apache.commons.math3.complex.Complex cv0 = new org.apache.commons.math3.complex.Complex(vRe0, -vIm0);
                    org.apache.commons.math3.complex.Complex cv1 = new org.apache.commons.math3.complex.Complex(vRe1, -vIm1);
                    org.apache.commons.math3.complex.Complex cw0 = new org.apache.commons.math3.complex.Complex(wRe0, -wIm0);
                    org.apache.commons.math3.complex.Complex cw1 = new org.apache.commons.math3.complex.Complex(wRe1, -wIm1);

                    org.apache.commons.math3.complex.Complex newRow0 = cv0.multiply(t_kj).add(cv1.multiply(t_k1j));
                    org.apache.commons.math3.complex.Complex newRow1 = cw0.multiply(t_kj).add(cw1.multiply(t_k1j));
                    T.set(k, j, newRow0);
                    T.set(k + 1, j, newRow1);
                }

                // Apply T <- T * G (affects cols k and k+1)
                for (int i = 0; i < n; i++) {
                    org.apache.commons.math3.complex.Complex t_ik = T.get(i, k);
                    org.apache.commons.math3.complex.Complex t_ik1 = T.get(i, k + 1);
                    org.apache.commons.math3.complex.Complex gv0 = new org.apache.commons.math3.complex.Complex(vRe0, vIm0);
                    org.apache.commons.math3.complex.Complex gv1 = new org.apache.commons.math3.complex.Complex(vRe1, vIm1);
                    org.apache.commons.math3.complex.Complex gw0 = new org.apache.commons.math3.complex.Complex(wRe0, wIm0);
                    org.apache.commons.math3.complex.Complex gw1 = new org.apache.commons.math3.complex.Complex(wRe1, wIm1);

                    org.apache.commons.math3.complex.Complex newCol0 = t_ik.multiply(gv0).add(t_ik1.multiply(gv1));
                    org.apache.commons.math3.complex.Complex newCol1 = t_ik.multiply(gw0).add(t_ik1.multiply(gw1));
                    T.set(i, k, newCol0);
                    T.set(i, k + 1, newCol1);
                }

                // Apply U <- U * G (affects cols k and k+1)
                for (int i = 0; i < n; i++) {
                    org.apache.commons.math3.complex.Complex u_ik = U.get(i, k);
                    org.apache.commons.math3.complex.Complex u_ik1 = U.get(i, k + 1);
                    org.apache.commons.math3.complex.Complex gv0 = new org.apache.commons.math3.complex.Complex(vRe0, vIm0);
                    org.apache.commons.math3.complex.Complex gv1 = new org.apache.commons.math3.complex.Complex(vRe1, vIm1);
                    org.apache.commons.math3.complex.Complex gw0 = new org.apache.commons.math3.complex.Complex(wRe0, wIm0);
                    org.apache.commons.math3.complex.Complex gw1 = new org.apache.commons.math3.complex.Complex(wRe1, wIm1);

                    org.apache.commons.math3.complex.Complex newCol0 = u_ik.multiply(gv0).add(u_ik1.multiply(gv1));
                    org.apache.commons.math3.complex.Complex newCol1 = u_ik.multiply(gw0).add(u_ik1.multiply(gw1));
                    U.set(i, k, newCol0);
                    U.set(i, k + 1, newCol1);
                }

                // Clean up: force subdiagonal to zero and set diagonal to exact eigenvalues
                T.set(k + 1, k, new org.apache.commons.math3.complex.Complex(0.0, 0.0));
                T.set(k, k, new org.apache.commons.math3.complex.Complex(realPart, imagPart));
                T.set(k + 1, k + 1, new org.apache.commons.math3.complex.Complex(realPart, -imagPart));

                k += 2;  // Skip past the 2x2 block
            } else {
                k += 1;
            }
        }

        Map<String, ComplexMatrix> result = new HashMap<>();
        result.put("U", U);
        result.put("T", T);
        return result;
    }

    /**
     * Sets the element at absolute index position to the given value.
     *
     * @param idx index in flattened format (row + col * numRows)
     * @param val value to set
     */
    public void set(int idx, double val) {
        if (idx >= this.getNumCols() * this.getNumRows())
            throw new RuntimeException("Index out of matrix");

        int row = idx % this.getNumRows();
        int col = idx / this.getNumRows();

        this.set(row, col, val);
    }

    /**
     * Sets the element at specified row and column.
     * Removes the entry if the value is zero to maintain sparse representation.
     *
     * @param row row index
     * @param col column index
     * @param val value to set
     */
    public void set(int row, int col, double val) {
        getData().set(row, col, val);
        if (val == 0)
            ((DMatrixSparseCSC) getData()).remove(row, col); // maintain sparsity
    }

    /**
     * Sets the element at specified row and column using an integer value.
     * Interprets Integer.MAX_VALUE as +Inf, MIN_VALUE as -Inf.
     *
     * @param row row index
     * @param col column index
     * @param val value to set
     */
    public void set(int row, int col, int val) {
        if (val == Integer.MAX_VALUE) {
            getData().set(row, col, Inf);
        } else if (val == Integer.MIN_VALUE) {
            getData().set(row, col, NegInf);
        } else {
            getData().set(row, col, val);
            if (val == 0)
                ((DMatrixSparseCSC) getData()).remove(row, col);
        }
    }

    /**
     * Sets the specified column to the values in the given vector.
     *
     * @param j   the index of the column to set
     * @param col the column vector to copy from
     * @return this matrix after the operation
     */
    public Matrix setColumn(int j, Matrix col) {
        for (int i = 0; i < this.getNumRows(); i++) {
            this.set(i, j, col.get(i));
        }
        return this;
    }

    private void setColumnIndex(int col, int value) {
        delegate.setColumnIndex(col, value);
    }

    /**
     * Sets multiple columns starting from column index j0 to j1 (exclusive) using values from another matrix.
     *
     * @param j0   starting column index (inclusive)
     * @param j1   ending column index (exclusive)
     * @param cols the matrix containing new column values
     * @return this matrix after the operation
     */
    public Matrix setColumns(int j0, int j1, Matrix cols) {
        if (cols.getNumCols() != j1 - j0) {
            line_error(mfilename(new Object[]{}), "Incorrect number of columns in Matrix argument.");
        }
        for (int k = 0; k < j1 - j0; k++) {
            for (int i = 0; i < this.getNumRows(); i++) {
                this.set(i, j0 + k, cols.get(i, k));
            }
        }
        return this;
    }

    /**
     * Replaces all NaN entries in the matrix with the specified value.
     *
     * @param val value to replace NaNs with
     */
    public void setNaNTo(double val) {
        for (int row = 0; row < getNumRows(); row++) {
            for (int col = 0; col < getNumCols(); col++) {
                if (Double.isNaN(this.get(row, col))) {
                    this.set(row, col, val);
                }
            }
        }
    }

    /**
     * Replaces all NaN entries in the matrix with zero.
     */
    public void setNaNToZero() {
        setNaNTo(0.0);
    }

    private void setNonZeroRow(int index, int row) {
        delegate.setNonZeroRow(index, row);
    }

    //    public Ret.Eigs eigvec() {
    //        if (this.data.numCols != this.data.numRows) {
    //            throw new RuntimeException("Only square matrix can be eigen decomposited");
    //        }
    //        DMatrixRMaj matrix = new DMatrixRMaj(this.data.numRows, this.data.numCols);
    //        DConvertMatrixStruct.convert(this.data, matrix);
    //
    //
    //        EigenDecomposition_F64<DMatrixRMaj> eig = DecompositionFactory_DDRM.eig(matrix.numCols,
    // true);
    //        eig.decompose(matrix);
    //
    //        Matrix D = new Matrix(1, data.numCols, data.numCols);
    //        Matrix V = new Matrix(data.numRows, data.numCols, data.numRows);
    //        for (int i = 0; i < data.numCols; i++) {
    //            if (eig.getEigenvalue(i).isReal()) {
    //                D.set(i, eig.getEigenvalue(i).getReal());
    //                DMatrixRMaj vector = eig.getEigenVector(i);
    //                for (int j = 0; j < data.numRows; j++) {
    //                    V.set(j, i, vector.get(j));
    //                }
    //            } else {
    //                for (int j = 0; j < data.numRows; j++) {
    //                    V.set(i, j, Double.NaN);
    //                }
    //            }
    //        }
    //
    //        if (V.hasNaN()) {
    //            System.out.println("Complex eigenvector detected, return NAN");
    //        }
    //        return new Ret.Eigs(D,  V);
    //    }

    private void setNonZeroValue(int index, double value) {
        delegate.setNonZeroValue(index, value);
    }

    /**
     * Sets the specified row to the values in the given vector.
     *
     * @param j   the index of the row to set
     * @param row the row vector to copy from
     * @return this matrix after the operation
     */
    public Matrix setRow(int j, Matrix row) {
        for (int i = 0; i < this.getNumCols(); i++) {
            this.set(j, i, row.get(i));
        }
        return this;
    }

    /**
     * Sets multiple rows starting from row index i0 to i1 (exclusive) using values from another matrix.
     *
     * @param i0   starting row index (inclusive)
     * @param i1   ending row index (exclusive)
     * @param rows the matrix containing new row values
     * @return this matrix after the operation
     */
    public Matrix setRows(int i0, int i1, Matrix rows) {
        if (rows.getNumRows() != i1 - i0) {
            line_error(mfilename(new Object[]{}), "Incorrect number of columns in Matrix argument.");
        }
        rows.print();
        for (int k = 0; k < i1 - i0; k++) {
            for (int i = 0; i < this.getNumRows(); i++) {
                this.set(i0 + k, i, rows.get(k, i));
            }
        }
        return this;
    }

    /**
     * Returns a new matrix representing the updated slice after assigning values from newSlice.
     *
     * @param rowStart starting row index (inclusive)
     * @param rowEnd   ending row index (exclusive)
     * @param colStart starting column index (inclusive)
     * @param colEnd   ending column index (exclusive)
     * @param newSlice the matrix containing the new values
     * @return the updated submatrix
     */
    public Matrix setSlice(int rowStart, int rowEnd, int colStart, int colEnd, Matrix newSlice) {
        Matrix subMatrix = new Matrix(rowEnd - rowStart + 1, colEnd - colStart + 1);
        for (int row = rowStart; row < rowEnd; row++) {
            for (int col = colStart; col < colEnd; col++) {
                subMatrix.set(row, col, newSlice.get(row - rowStart, col - colStart));
            }
        }
        return subMatrix;
    }

    /**
     * Sets the values of a submatrix (in-place) using the specified newSlice matrix.
     *
     * @param rowStart starting row index (inclusive)
     * @param rowEnd   ending row index (exclusive)
     * @param colStart starting column index (inclusive)
     * @param colEnd   ending column index (exclusive)
     * @param newSlice the matrix to copy values from
     */
    public void setSliceEq(int rowStart, int rowEnd, int colStart, int colEnd, Matrix newSlice) {
        for (int row = rowStart; row < rowEnd; row++) {
            for (int col = colStart; col < colEnd; col++) {
                this.set(row, col, newSlice.get(row - rowStart, col - colStart));
            }
        }
    }

    /**
     * Copies the data from another matrix into this matrix.
     *
     * @param m the matrix to copy from
     */
    public void setTo(Matrix m) {
        getData().setTo(m.getData());
    }

    /**
     * Sets all elements of the matrix to NaN.
     */
    public void setToNaN() {
        for (int row = 0; row < getNumRows(); row++) {
            for (int col = 0; col < getNumCols(); col++) {
                this.set(row, col, Double.NaN);
            }
        }
    }

    public void shrinkNumCols(int newmax) {
        delegate.shrinkNumCols(newmax);
    }

    public void shrinkNumRows(int newmax) {
        delegate.shrinkNumRows(newmax);
    }

    /**
     * Returns a new matrix with the non-zero values sorted in ascending order.
     *
     * @return a new matrix with sorted values
     */
    public Matrix sort() {
        Matrix ret = this.copy();
        Arrays.sort(ret.getNonZeroValues());
        return ret;
    }

    /**
     * Sorts the current matrix's non-zero values in place.
     *
     * @return this matrix after sorting
     */
    public Matrix sortEq() {
        Arrays.sort(this.getNonZeroValues());
        return this;
    }

    /**
     * Computes the element-wise square root of the matrix.
     *
     * @return the resulting matrix
     */
    public Matrix sqrt() {
        return this.elementPower(0.5);
    }

    /**
     * Computes the matrix multiplied by itself.
     *
     * @return the squared matrix
     */
    public Matrix square() {
        return this.mult(this.copy());
    }

    /**
     * Subtract a scalar from all elements.
     *
     * @param x the scalar value to subtract
     * @return the resulting matrix
     */
    public Matrix sub(double x) {
        return this.sub(Matrix.ones(this.getNumRows(), this.getNumCols()).scale(x));
    }

    /**
     * Subtracts another matrix.
     *
     * @param matrix the matrix to subtract
     * @return the resulting matrix
     */
    public Matrix sub(Matrix matrix) {
        return sub(1.0, matrix);
    }

    /**
     * Subtracts alpha-scaled version of the provided matrix.
     *
     * @param alpha  the scalar multiplier
     * @param matrix the matrix to subtract
     * @return the resulting matrix
     */
    public Matrix sub(double alpha, Matrix matrix) {
        return new Matrix(addMatrices(1, this.getData(), -alpha, matrix.getData(), null));
    }

    /**
     * Subtracts a scalar from the current matrix in place.
     *
     * @param x the scalar to subtract
     */
    public void subEq(double x) {
        this.subEq(Matrix.ones(this.getNumRows(), this.getNumCols()).scale(x));
    }

    /**
     * Subtracts a matrix from the current matrix in place.
     *
     * @param matrix the matrix to subtract
     */
    public void subEq(Matrix matrix) {
        subEq(1.0, matrix);
    }

    /**
     * Subtracts a scaled matrix from the current matrix in place.
     *
     * @param alpha  the scale factor
     * @param matrix the matrix to subtract
     */
    public void subEq(double alpha, Matrix matrix) {
        this.setTo(this.sub(alpha, matrix));
    }

    /**
     * Sums the absolute values of the given column.
     *
     * @param col the column index
     * @return the sum of absolute values in the column
     */
    public double sumAbsCols(int col) {
        double sum = 0;
        for (int i = 0; i < this.getNumRows(); i++) {
            sum += FastMath.abs(this.get(i, col));
        }
        return sum;
    }

    /**
     * Sums the absolute values of the given row.
     *
     * @param row the row index
     * @return the sum of absolute values in the row
     */
    public double sumAbsRows(int row) {
        double sum = 0;
        for (int i = 0; i < this.getNumRows(); i++) {
            sum += FastMath.abs(this.get(row, i));
        }
        return sum;
    }

    /**
     * Sums the elements in the specified column.
     *
     * @param col the column index
     * @return the sum of the column values
     */
    public double sumCols(int col) {
        double sum = 0;
        int numRows = this.getNumRows();
        for (int i = 0; i < numRows; i++) {
            sum += this.get(i, col);
        }
        return sum;
    }

    /**
     * Sums the values in each column and returns the results as a row vector.
     *
     * @return a row vector with column sums
     */
    public Matrix sumCols() {
        DMatrix sumcols = delegate.sumColsRaw();
        DMatrixSparseCSC tmp = new DMatrixSparseCSC(0, 0);
        DConvertMatrixStruct.convert(sumcols, tmp);
        return new Matrix(tmp);
    }

    /**
     * Computes the sum of a subset of rows for each column.
     *
     * @param startRow the starting row index (inclusive)
     * @param endRow   the ending row index (exclusive)
     * @return a row vector with the sums of each column over the specified row range
     */
    public Matrix sumCols(int startRow, int endRow) {
        int numCols = this.getNumCols();
        Matrix result = new Matrix(1, numCols);
        for (int i = 0; i < numCols; i++) {
            double sum = 0;
            for (int j = startRow; j < endRow; j++) {
                sum += this.get(i, j);
            }
            result.set(0, i, sum);
        }
        return result;
    }

    /**
     * Computes the sum of the values in a specific row.
     *
     * @param row the row index
     * @return the sum of the row values
     */
    public double sumRows(int row) {
        double sum = 0;
        int numCols = this.getNumCols();
        for (int i = 0; i < numCols; i++) {
            sum += this.get(row, i);
        }
        return sum;
    }

    /**
     * Computes the sum of each row and returns the result as a column vector.
     *
     * @return a column vector containing the sum of each row
     */
    public Matrix sumRows() {
        DMatrix sumrows = delegate.sumRowsRaw();
        DMatrixSparseCSC tmp = new DMatrixSparseCSC(0, 0);
        DConvertMatrixStruct.convert(sumrows, tmp);
        return new Matrix(tmp);
    }

    /**
     * Computes the sum over a subrange of columns for each row.
     *
     * @param startCol the starting column (inclusive)
     * @param endCol   the ending column (exclusive)
     * @return a column vector with the row sums over the specified columns
     */
    public Matrix sumRows(int startCol, int endCol) {
        Matrix result = new Matrix(this.getNumRows(), 1);
        for (int i = 0; i < this.getNumRows(); i++) {
            double sum = 0;
            for (int j = startCol; j < endCol; j++) {
                sum += this.get(i, j);
            }
            result.set(i, 0, sum);
        }
        return result;
    }

    /**
     * Computes the sum of the values in a rectangular submatrix.
     *
     * @param startRow the starting row (inclusive)
     * @param endRow   the ending row (exclusive)
     * @param startCol the starting column (inclusive)
     * @param endCol   the ending column (exclusive)
     * @return the sum of the submatrix values
     */
    public double sumSubMatrix(int startRow, int endRow, int startCol, int endCol) {
        double sum = 0;
        for (int i = startRow; i < endRow; i++) {
            for (int j = startCol; j < endCol; j++) {
                sum += this.get(i, j);
            }
        }
        return sum;
    }

    /**
     * Computes the sum of a submatrix defined by specific row and column indices.
     *
     * @param rowIndexes the array of row indices
     * @param colIndexes the array of column indices
     * @return the sum of the selected submatrix
     */
    public double sumSubMatrix(int[] rowIndexes, int[] colIndexes) {
        double sum = 0;
        for (int rowIndex : rowIndexes) {
            for (int colIndex : colIndexes) {
                sum += this.get(rowIndex, colIndex);
            }
        }
        return sum;
    }

    /**
     * Computes the sum of a submatrix defined by boolean selection flags.
     *
     * @param rowIndexes boolean array indicating selected rows
     * @param colIndexes boolean array indicating selected columns
     * @return the sum of the selected submatrix
     */
    public double sumSubMatrix(boolean[] rowIndexes, boolean[] colIndexes) {
        double sum = 0;
        for (int i = 0; i < getNumRows(); i++) {
            if (rowIndexes[i]) {
                for (int j = 0; j < getNumCols(); j++) {
                    if (colIndexes[j]) {
                        sum += this.get(i, j);
                    }
                }
            }
        }
        return sum;
    }

    /**
     * Converts the matrix to a flat 1D array in row-major order.
     *
     * @return the 1D array representation of the matrix
     */
    public double[] toArray1D() {
        double[] array = new double[getNumRows() * getNumCols()];
        int k = 0;
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                array[k] = this.get(i, j);
                k++;
            }
        }
        return array;
    }

    /**
     * Converts the matrix to a 2D array representation.
     *
     * @return the 2D array representation of the matrix
     */
    public double[][] toArray2D() {
        int numRows = getNumRows();
        int numCols = getNumCols();
        double[][] array = new double[numRows][numCols];
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                array[i][j] = this.get(i, j);
            }
        }
        return array;
    }

    /**
     * Converts this matrix to a copy of its underlying {@link DMatrixSparseCSC} structure.
     *
     * @return a deep copy of the internal sparse matrix representation
     */
    public DMatrixSparseCSC toDMatrixSparseCSC() {
        return this.getData().copy();
    }

    /**
     * Converts a specified matrix to a copy of its underlying {@link DMatrixSparseCSC} structure.
     *
     * @param matrix the matrix to convert
     * @return a deep copy of the sparse matrix representation of the given matrix
     */
    public DMatrixSparseCSC toDMatrixSparseCSC(Matrix matrix) {
        return matrix.getData().copy();
    }

    /**
     * Converts this matrix to a single {@code Double} value.
     * Assumes the matrix is scalar (1x1), otherwise behavior is undefined.
     *
     * @return the single value contained in the matrix
     */
    public Double toDouble() {
        return this.value();
    }

    /**
     * Converts the matrix to a nested list of {@code Double} values.
     *
     * @return a list of lists representing rows and columns of the matrix
     */
    public List<List<Double>> toDoubleList() {
        List<List<Double>> array = new ArrayList<>();
        for (int i = 0; i < getNumRows(); i++) {
            List<Double> row = new ArrayList<>();
            for (int j = 0; j < getNumCols(); j++) {
                row.add(this.get(i, j));
            }
            array.add(row);
        }
        return array;
    }

    /**
     * Converts the matrix to a 1D array of {@code int}, flattening in row-major order.
     *
     * @return a 1D int array representing the matrix values cast to integers
     */
    public int[] toIntArray1D() {
        int[] array = new int[getNumRows() * getNumCols()];
        int k = 0;
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                array[k] = (int) this.get(i, j);
                k++;
            }
        }
        return array;
    }

    /**
     * Converts the matrix into a 1D {@link List} of {@code Double} values in row-major order.
     *
     * @return a list containing all matrix elements in row-major order
     */
    public List<Double> toList1D() {
        List<Double> list = new ArrayList<>();
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                list.add(this.get(i, j));
            }
        }
        return list;
    }

    /**
     * Returns a formatted string representation of the matrix.
     * Each row is separated by a semicolon and new line.
     *
     * @return the string representation of the matrix
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < this.getNumRows(); i++) {
            for (int j = 0; j < this.getNumCols(); j++) {
                if (i == 0 && j == 0) {
                    sb.append(" " + this.get(i, j) + " ");
                } else {
                    sb.append("  ").append(this.get(i, j)).append(" ");
                }
            }
            if (i < this.getNumRows() - 1) {
                sb.append(";\n");
            } else {
                sb.append("]\n");
            }
        }
        return sb.toString();
    }

    /**
     * Computes the transpose of the current matrix.
     *
     * @return a new matrix representing the transpose of this matrix
     */
    public Matrix transpose() {
        Matrix res = new Matrix(0, 0);
        transposeMatrix(res.getData());
        return res;
    }

    private void transposeMatrix(DMatrix output) {
        delegate.transposeMatrix(output);
    }

    /**
     * Finds unique integer values in the specified column of the matrix.
     *
     * @param colIdx the index of the column to search
     * @return a column matrix of unique integer values found in the specified column
     */
    public Matrix uniqueInCol(int colIdx) {
        List<Integer> array = new ArrayList<>();
        for (int rowIdx = 0; rowIdx < this.getNumRows(); rowIdx++) {
            array.add((int) this.get(rowIdx, colIdx));
        }

        List<Integer> unique_array = Utils.unique(array);
        Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
        for (int i = 0; i < unique_array.size(); i++) {
            res.set(i, 0, unique_array.get(i));
        }
        return res;
    }

    /**
     * Finds unique integer values in the specified row of the matrix.
     *
     * @param rowIdx the index of the row to search
     * @return a column matrix of unique integer values found in the specified row
     */
    public Matrix uniqueInRow(int rowIdx) {
        List<Integer> array = new ArrayList<>();
        for (int colIdx = 0; colIdx < this.getNumCols(); colIdx++) {
            array.add((int) this.get(rowIdx, colIdx));
        }

        List<Integer> unique_array = Utils.unique(array);
        Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
        for (int i = 0; i < unique_array.size(); i++) {
            res.set(i, 0, unique_array.get(i));
        }
        return res;
    }

    /**
     * Finds unique positive values (strictly greater than 0) in the specified column.
     *
     * @param colIdx the index of the column to search
     * @return a column matrix of unique positive values
     */
    public Matrix uniqueNonNegativeInCol(int colIdx) {
        List<Integer> array = new ArrayList<>();
        for (int rowIdx = 0; rowIdx < this.getNumRows(); rowIdx++) {
            int val = (int) this.get(rowIdx, colIdx);
            if (val > 0) {
                array.add(val);
            }
        }

        List<Integer> unique_array = Utils.unique(array);
        Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
        for (int i = 0; i < unique_array.size(); i++) res.set(i, 0, unique_array.get(i));
        return res;
    }

    /**
     * Finds unique positive values (strictly greater than 0) in the specified row.
     *
     * @param rowIdx the index of the row to search
     * @return a column matrix of unique positive values
     */
    public Matrix uniqueNonNegativeInRow(int rowIdx) {
        List<Integer> array = new ArrayList<>();
        for (int colIdx = 0; colIdx < this.getNumCols(); colIdx++) {
            int val = (int) this.get(rowIdx, colIdx);
            if (val > 0) {
                array.add(val);
            }
        }

        List<Integer> unique_array = Utils.unique(array);
        Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
        for (int i = 0; i < unique_array.size(); i++) res.set(i, 0, unique_array.get(i));
        return res;
    }

    /**
     * Finds unique non-zero values in the specified column.
     *
     * @param colIdx the index of the column to search
     * @return a column matrix of unique non-zero values
     */
    public Matrix uniqueNonZerosInCol(int colIdx) {
        List<Integer> array = new ArrayList<>();
        for (int rowIdx = 0; rowIdx < this.getNumRows(); rowIdx++) {
            int val = (int) this.get(rowIdx, colIdx);
            if (val != 0) {
                array.add(val);
            }
        }

        List<Integer> unique_array = Utils.unique(array);
        Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
        for (int i = 0; i < unique_array.size(); i++) res.set(i, 0, unique_array.get(i));
        return res;
    }

    /**
     * Finds unique non-zero values in the specified row.
     *
     * @param rowIdx the index of the row to search
     * @return a column matrix of unique non-zero values
     */
    public Matrix uniqueNonZerosInRow(int rowIdx) {
        List<Integer> array = new ArrayList<>();
        for (int colIdx = 0; colIdx < this.getNumCols(); colIdx++) {
            int val = (int) this.get(rowIdx, colIdx);
            if (val != 0) {
                array.add(val);
            }
        }

        List<Integer> unique_array = Utils.unique(array);
        Matrix res = new Matrix(unique_array.size(), 1, unique_array.size());
        for (int i = 0; i < unique_array.size(); i++) res.set(i, 0, unique_array.get(i));
        return res;
    }

    /**
     * Sets a value in the matrix without bounds checking.
     * Use this only if you're certain the row and column are valid.
     *
     * @param row the row index
     * @param col the column index
     * @param val the value to set
     */
    public void unsafeSet(int row, int col, double val) {
        getData().unsafe_set(row, col, val);
    }

    /**
     * Returns the value at position (0, 0).
     * Useful for singleton matrices (1x1).
     *
     * @return the scalar value at the top-left of the matrix
     */
    public double value() {
        if (getNumRows() == 0 || getNumCols() == 0) {
            return 0.0;
        }
        return getData().get(0, 0);
    }

    /**
     * Sets all entries in the matrix to zero.
     */
    public void zero() {
        getData().zero();
    }

    /**
     * Efficient hash key for matrix rows used in uniqueRowIndexes operation.
     * Uses numeric hashing instead of string conversion for better performance
     * with large matrices (>10,000 samples).
     */
    private static class RowHashKey implements Comparable<RowHashKey> {
        private final double[] rowData;
        private final int hashCode;

        /**
         * Creates a hash key for a specific row in a matrix
         *
         * @param matrix The matrix containing the row
         * @param rowIndex The index of the row to use as key
         */
        RowHashKey(Matrix matrix, int rowIndex) {
            int numCols = matrix.getNumCols();
            this.rowData = new double[numCols];

            for (int j = 0; j < numCols; j++) {
                this.rowData[j] = matrix.get(rowIndex, j);
            }

            this.hashCode = Arrays.hashCode(this.rowData);
        }

        @Override
        public int hashCode() {
            return hashCode;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof RowHashKey)) {
                return false;
            }
            RowHashKey other = (RowHashKey) obj;

            if (this.hashCode != other.hashCode) {
                return false;
            }

            return Arrays.equals(this.rowData, other.rowData);
        }

        @Override
        public int compareTo(RowHashKey other) {
            int minLen = Math.min(this.rowData.length, other.rowData.length);
            for (int i = 0; i < minLen; i++) {
                int cmp = Double.compare(this.rowData[i], other.rowData[i]);
                if (cmp != 0) {
                    return cmp;
                }
            }
            return Integer.compare(this.rowData.length, other.rowData.length);
        }
    }

}
