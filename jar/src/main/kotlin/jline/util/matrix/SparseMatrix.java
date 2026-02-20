/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import jline.util.Utils;
import org.ejml.data.DMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.sparse.csc.CommonOps_DSCC;

/**
 * Base class for sparse matrix implementations, containing the core data structure
 * and methods that directly manipulate the underlying sparse matrix representation.
 *
 * <p>This class encapsulates all interactions with EJML's {@link DMatrixSparseCSC} and
 * {@link CommonOps_DSCC} classes, providing a clean abstraction layer for sparse matrix
 * operations. It serves as the foundation for concrete sparse matrix implementations like
 * {@code Matrix}.</p>
 *
 * <p>The class provides:</p>
 * <ul>
 *   <li>Basic matrix operations (get, set, dimensions)</li>
 *   <li>Encapsulated CommonOps_DSCC method calls</li>
 *   <li>Memory management for sparse matrix data structures</li>
 *   <li>Utility methods for matrix manipulation</li>
 * </ul>
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */
public abstract class SparseMatrix extends BaseMatrix {

    /**
     * The underlying EJML sparse matrix data structure
     */
    private DMatrixSparseCSC data;

    /**
     * Constructs a sparse matrix with specified dimensions and initial capacity.
     *
     * @param numRows     the number of rows in the matrix
     * @param numCols     the number of columns in the matrix
     * @param arrayLength the initial capacity for non-zero elements
     */
    public SparseMatrix(int numRows, int numCols, int arrayLength) {
        setData(new DMatrixSparseCSC(numRows, numCols, arrayLength));
    }

    /**
     * Constructs an empty sparse matrix with specified dimensions.
     *
     * @param numRows the number of rows in the matrix
     * @param numCols the number of columns in the matrix
     */
    public SparseMatrix(int numRows, int numCols) {
        setData(new DMatrixSparseCSC(numRows, numCols, 0));
    }

    /**
     * Constructs a sparse matrix by copying an existing EJML sparse matrix.
     *
     * @param matrix the EJML sparse matrix to copy
     */
    public SparseMatrix(DMatrix matrix) {
        setData(new DMatrixSparseCSC((DMatrixSparseCSC) matrix));
    }

    /**
     * Copy constructor for creating a new sparse matrix from an existing one.
     *
     * @param matrix the sparse matrix to copy
     */
    protected SparseMatrix(SparseMatrix matrix) {
        setData(new DMatrixSparseCSC(matrix.getData().copy()));
    }

    /**
     * Protected constructor for delayed initialization.
     * Subclasses must ensure that data is properly initialized.
     */
    protected SparseMatrix() {
        // data will be initialized later by subclass
    }

    protected static void addMatricesInPlaceStatic(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output) {
        CommonOps_DSCC.add(alpha, (DMatrixSparseCSC) A, beta, (DMatrixSparseCSC) B, (DMatrixSparseCSC) output, null, null);
    }

    protected static DMatrixSparseCSC addMatricesStatic(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output) {
        return CommonOps_DSCC.add(alpha, (DMatrixSparseCSC) A, beta, (DMatrixSparseCSC) B, (DMatrixSparseCSC) output, null, null);
    }

    protected static void changeSignStatic(DMatrix input, DMatrix output) {
        CommonOps_DSCC.changeSign((DMatrixSparseCSC) input, (DMatrixSparseCSC) output);
    }

    protected static void concatColumnsInPlaceStatic(DMatrix left, DMatrix right, DMatrix output) {
        CommonOps_DSCC.concatColumns((DMatrixSparseCSC) left, (DMatrixSparseCSC) right, (DMatrixSparseCSC) output);
    }

    protected static DMatrixSparseCSC concatColumnsStatic(DMatrix left, DMatrix right, DMatrix output) {
        return CommonOps_DSCC.concatColumns((DMatrixSparseCSC) left, (DMatrixSparseCSC) right, (DMatrixSparseCSC) output);
    }

    protected static void concatRowsInPlaceStatic(DMatrix top, DMatrix bottom, DMatrix output) {
        CommonOps_DSCC.concatRows((DMatrixSparseCSC) top, (DMatrixSparseCSC) bottom, (DMatrixSparseCSC) output);
    }

    protected static DMatrixSparseCSC concatRowsStatic(DMatrix top, DMatrix bottom, DMatrix output) {
        return CommonOps_DSCC.concatRows((DMatrixSparseCSC) top, (DMatrixSparseCSC) bottom, (DMatrixSparseCSC) output);
    }

    protected static double determinantStatic(DMatrix matrix) {
        return CommonOps_DSCC.det((DMatrixSparseCSC) matrix);
    }

    // Diag methods
    protected static DMatrixSparseCSC diagStatic(double[] values) {
        return CommonOps_DSCC.diag(values);
    }

    protected static DMatrixSparseCSC diagWithMatrixStatic(DMatrixSparseCSC A, double[] values, int offset, int length) {
        return CommonOps_DSCC.diag(A, values, offset, length);
    }

    protected static DMatrixSparseCSC diagWithRetStatic(DMatrixSparseCSC ret, double[] values, int offset, int length) {
        return CommonOps_DSCC.diag(ret, values, offset, length);
    }

    protected static void divideInPlaceStatic(DMatrix matrix, double scalar) {
        DMatrixSparseCSC sparseMatrix = (DMatrixSparseCSC) matrix;
        CommonOps_DSCC.divide(sparseMatrix, scalar, sparseMatrix);
    }

    protected static void divideMatrixStatic(DMatrix input, double scalar, DMatrix output) {
        CommonOps_DSCC.divide((DMatrixSparseCSC) input, scalar, (DMatrixSparseCSC) output);
    }

    protected static void divideRowsByArrayStatic(double[] diag, int offset, DMatrix matrix) {
        CommonOps_DSCC.divideRows(diag, offset, (DMatrixSparseCSC) matrix);
    }

    protected static void divideScalarByMatrixStatic(double scalar, DMatrix input, DMatrix output) {
        CommonOps_DSCC.divide(scalar, (DMatrixSparseCSC) input, (DMatrixSparseCSC) output);
    }

    protected static DMatrixSparseCSC elementMultStatic(DMatrix A, DMatrix B, DMatrix output) {
        return CommonOps_DSCC.elementMult((DMatrixSparseCSC) A, (DMatrixSparseCSC) B, (DMatrixSparseCSC) output, null, null);
    }

    protected static void extractColumnInPlaceStatic(DMatrix input, int column, DMatrix output) {
        CommonOps_DSCC.extractColumn((DMatrixSparseCSC) input, column, (DMatrixSparseCSC) output);
    }

    protected static DMatrixSparseCSC extractColumnStatic(DMatrix input, int column, DMatrix output) {
        return CommonOps_DSCC.extractColumn((DMatrixSparseCSC) input, column, (DMatrixSparseCSC) output);
    }

    protected static void extractDiagStatic(DMatrix input, DMatrix output) {
        CommonOps_DSCC.extractDiag((DMatrixSparseCSC) input, (DMatrixSparseCSC) output);
    }

    protected static void extractMatrixStatic(DMatrix src, int srcX0, int srcX1, int srcY0, int srcY1, DMatrix dst, int dstY0, int dstX0) {
        CommonOps_DSCC.extract((DMatrixSparseCSC) src, srcX0, srcX1, srcY0, srcY1, (DMatrixSparseCSC) dst, dstY0, dstX0);
    }

    protected static void extractRowsInPlaceStatic(DMatrix input, int row0, int row1, DMatrix output) {
        CommonOps_DSCC.extractRows((DMatrixSparseCSC) input, row0, row1, (DMatrixSparseCSC) output);
    }

    protected static DMatrixSparseCSC extractRowsStatic(DMatrix input, int row0, int row1, DMatrix output) {
        return CommonOps_DSCC.extractRows((DMatrixSparseCSC) input, row0, row1, (DMatrixSparseCSC) output);
    }

    protected static void fillMatrixStatic(DMatrix matrix, double value) {
        CommonOps_DSCC.fill((DMatrixSparseCSC) matrix, value);
    }

    // Identity method
    protected static DMatrixSparseCSC identityStatic(int length) {
        return CommonOps_DSCC.identity(length);
    }

    // Invert method (static version for Matrix.java)
    protected static void invertStatic(DMatrix input, DMatrixRMaj output) {
        CommonOps_DSCC.invert((DMatrixSparseCSC) input, output);
    }

    protected static DMatrixSparseCSC multMatrixStatic(DMatrix A, DMatrix B, DMatrix output) {
        return CommonOps_DSCC.mult((DMatrixSparseCSC) A, (DMatrixSparseCSC) B, (DMatrixSparseCSC) output);
    }

    protected static void removeZerosWithTolStatic(DMatrix matrix, double tolerance) {
        CommonOps_DSCC.removeZeros((DMatrixSparseCSC) matrix, tolerance);
    }

    protected static void scaleMatrixStatic(double scalar, DMatrix input, DMatrix output) {
        CommonOps_DSCC.scale(scalar, (DMatrixSparseCSC) input, (DMatrixSparseCSC) output);
    }

    // Solve method
    protected static boolean solveStatic(DMatrix a, DMatrix b, DMatrix x) {
        return CommonOps_DSCC.solve((DMatrixSparseCSC) a, (DMatrixSparseCSC) b, (DMatrixSparseCSC) x);
    }

    protected static DMatrixRMaj sumColsStatic(DMatrix matrix) {
        return CommonOps_DSCC.sumCols((DMatrixSparseCSC) matrix, null);
    }

    protected static DMatrixRMaj sumRowsStatic(DMatrix matrix) {
        return CommonOps_DSCC.sumRows((DMatrixSparseCSC) matrix, null);
    }

    protected static void transposeMatrixStatic(DMatrix input, DMatrix output) {
        CommonOps_DSCC.transpose((DMatrixSparseCSC) input, (DMatrixSparseCSC) output, null);
    }

    /**
     * Replaces each value in the matrix with its absolute value, in-place.
     */
    public void absEq() {
        for (int i = 0; i < data.nz_length; i++) {
            this.data.nz_values[i] = Math.abs(this.data.nz_values[i]);
        }
    }

    /**
     * Performs the operation: output = alpha*A + beta*B.
     *
     * @param alpha  scalar coefficient for matrix A
     * @param A      first input matrix
     * @param beta   scalar coefficient for matrix B
     * @param B      second input matrix
     * @param output the result matrix
     */
    protected void add(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output) {
        CommonOps_DSCC.add(alpha, (DMatrixSparseCSC) A, beta, (DMatrixSparseCSC) B, (DMatrixSparseCSC) output, null, null);
    }

    protected DMatrix addMatrices(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output) {
        return CommonOps_DSCC.add(alpha, (DMatrixSparseCSC) A, beta, (DMatrixSparseCSC) B, (DMatrixSparseCSC) output, null, null);
    }

    protected void addMatricesInPlace(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output) {
        CommonOps_DSCC.add(alpha, (DMatrixSparseCSC) A, beta, (DMatrixSparseCSC) B, (DMatrixSparseCSC) output, null, null);
    }

    /**
     * Checks if the current matrix contains any non-zero elements.
     *
     * @return {@code true} if at least one non-zero element exists; {@code false} otherwise.
     */
    public boolean any() {
        return this.data.nz_length > 0;
    }

    /**
     * Applies a transformation to matrix elements based on conditions.
     * Replaces elements matching the condition with the target value.
     */
    protected void applyConditionalTransform(double source, double target, double tol, String operation) {
        int[] colIndices = data.col_idx;
        int[] rowIndices = data.nz_rows;
        double[] values = data.nz_values;

        for (int colIdx = 0; colIdx < data.numCols; colIdx++) {
            switch (operation) {
                case "NaN":
                    int col1 = colIndices[colIdx];
                    int col2 = colIndices[colIdx + 1];
                    for (int i = col1; i < col2; i++) {
                        if (Double.isNaN(values[i])) {
                            this.set(rowIndices[i], colIdx, target);
                        }
                    }
                    break;
                case "Inf":
                    col1 = colIndices[colIdx];
                    col2 = colIndices[colIdx + 1];
                    for (int i = col1; i < col2; i++) {
                        if (Double.isInfinite(values[i])) {
                            this.set(rowIndices[i], colIdx, target);
                        }
                    }
                    break;
                case "Value":
                    col1 = colIndices[colIdx];
                    col2 = colIndices[colIdx + 1];
                    for (int i = col1; i < col2; i++) {
                        if (Math.abs(values[i] - source) < tol) {
                            this.set(rowIndices[i], colIdx, target);
                        }
                    }
                    break;
            }
        }
    }

    // ChangeSign method
    protected void changeSign() {
        CommonOps_DSCC.changeSign(this.getData(), this.getData());
    }

    /**
     * Creates and returns a deep copy of this sparse matrix.
     *
     * @return a new sparse matrix that is a copy of this instance
     */
    public abstract BaseMatrix copy();

    /**
     * Adds a scalar value to each element in the specified column.
     *
     * @param col Column index to update.
     * @param a   Scalar value to add to each element in the column.
     */
    public void colIncrease(int col, double a) {
        int[] colIndices = data.col_idx;
        int[] rowIndices = data.nz_rows;
        double[] values = data.nz_values;

        if (col < 0 || col >= data.numCols) {
            throw new IllegalArgumentException("Column index out of bounds");
        }

        int start = colIndices[col];
        int end = colIndices[col + 1];

        for (int i = start; i < end; i++) {
            values[i] += a;
        }
    }

    protected DMatrixSparseCSC concatColumns(DMatrix left, DMatrix right, DMatrix output) {
        return CommonOps_DSCC.concatColumns((DMatrixSparseCSC) left, (DMatrixSparseCSC) right, (DMatrixSparseCSC) output);
    }

    protected void concatColumnsInPlace(DMatrix left, DMatrix right, DMatrix output) {
        CommonOps_DSCC.concatColumns((DMatrixSparseCSC) left, (DMatrixSparseCSC) right, (DMatrixSparseCSC) output);
    }

    protected DMatrixSparseCSC concatRows(DMatrix top, DMatrix bottom, DMatrix output) {
        return CommonOps_DSCC.concatRows((DMatrixSparseCSC) top, (DMatrixSparseCSC) bottom, (DMatrixSparseCSC) output);
    }

    protected void concatRowsInPlace(DMatrix top, DMatrix bottom, DMatrix output) {
        CommonOps_DSCC.concatRows((DMatrixSparseCSC) top, (DMatrixSparseCSC) bottom, (DMatrixSparseCSC) output);
    }

    /**
     * Counts occurrences of each unique value in each row.
     */
    public BaseMatrix countEachRow(double val) {
        BaseMatrix res = createNewInstance(data.numRows, 1, data.numRows);
        int[] colIndices = data.col_idx;
        int[] rowIndices = data.nz_rows;
        double[] values = data.nz_values;

        for (int colIdx = 0; colIdx < data.numCols; colIdx++) {
            int col1 = colIndices[colIdx];
            int col2 = colIndices[colIdx + 1];

            for (int i = col1; i < col2; i++) {
                int rowIdx = rowIndices[i];
                if (values[i] == val) {
                    res.set(rowIdx, 0, res.get(rowIdx, 0) + 1);
                }
            }
        }
        return res;
    }

    /**
     * Creates a new instance of the concrete sparse matrix implementation.
     * This factory method allows the base class to create instances of the correct subtype.
     *
     * @param rows     the number of rows for the new matrix
     * @param cols     the number of columns for the new matrix
     * @param nzLength the initial capacity for non-zero elements
     * @return a new instance of the concrete sparse matrix type
     */
    protected abstract BaseMatrix createNewInstance(int rows, int cols, int nzLength);

    // Determinant method
    protected double determinant() {
        return CommonOps_DSCC.det(this.getData());
    }

    protected void divideInPlace(double scalar) {
        CommonOps_DSCC.divide(this.getData(), scalar, this.getData());
    }

    // Divide methods
    protected void divideMatrix(double scalar, DMatrix output) {
        CommonOps_DSCC.divide(this.getData(), scalar, (DMatrixSparseCSC) output);
    }

    // DivideRows method
    protected void divideRowsByArray(double[] diag, int offset) {
        CommonOps_DSCC.divideRows(diag, offset, this.getData());
    }

    protected void divideScalarByMatrix(double scalar, DMatrix output) {
        CommonOps_DSCC.divide(scalar, this.getData(), (DMatrixSparseCSC) output);
    }

    // Element methods
    protected double elementMax() {
        return CommonOps_DSCC.elementMax(this.getData());
    }

    protected double elementMaxAbs() {
        return CommonOps_DSCC.elementMaxAbs(this.getData());
    }

    protected double elementMin() {
        return CommonOps_DSCC.elementMin(this.getData());
    }

    protected DMatrix elementMult(DMatrix A, DMatrix B, DMatrix output) {
        return CommonOps_DSCC.elementMult((DMatrixSparseCSC) A, (DMatrixSparseCSC) B, (DMatrixSparseCSC) output, null, null);
    }

    protected double elementSum() {
        return CommonOps_DSCC.elementSum(this.getData());
    }

    /**
     * Indicates whether some other object is "equal to" this sparse matrix.
     *
     * @param obj the reference object with which to compare
     * @return {@code true} if this object is the same as the obj argument; {@code false} otherwise
     */
    public abstract boolean equals(Object obj);

    protected DMatrixSparseCSC extractColumn(DMatrix input, int column, DMatrix output) {
        return CommonOps_DSCC.extractColumn((DMatrixSparseCSC) input, column, (DMatrixSparseCSC) output);
    }

    protected void extractColumnInPlace(DMatrix input, int column, DMatrix output) {
        CommonOps_DSCC.extractColumn((DMatrixSparseCSC) input, column, (DMatrixSparseCSC) output);
    }

    protected void extractDiag(DMatrix input, DMatrix output) {
        CommonOps_DSCC.extractDiag((DMatrixSparseCSC) input, (DMatrixSparseCSC) output);
    }

    // Extract methods
    protected void extractMatrix(DMatrix src, int srcX0, int srcX1, int srcY0, int srcY1, DMatrix dst, int dstY0, int dstX0) {
        CommonOps_DSCC.extract((DMatrixSparseCSC) src, srcX0, srcX1, srcY0, srcY1, (DMatrixSparseCSC) dst, dstY0, dstX0);
    }

    protected DMatrixSparseCSC extractRows(DMatrix input, int row0, int row1, DMatrix output) {
        return CommonOps_DSCC.extractRows((DMatrixSparseCSC) input, row0, row1, (DMatrixSparseCSC) output);
    }

    protected void extractRowsInPlace(DMatrix input, int row0, int row1, DMatrix output) {
        CommonOps_DSCC.extractRows((DMatrixSparseCSC) input, row0, row1, (DMatrixSparseCSC) output);
    }

    // Fill methods
    protected void fillMatrix(double value) {
        CommonOps_DSCC.fill(this.getData(), value);
    }

    /**
     * Returns a matrix containing indices of all non-negative elements.
     */
    public BaseMatrix findNonNegative() {
        double[] values = data.nz_values;
        int[] rowIndices = data.nz_rows;
        int[] colIndices = data.col_idx;

        java.util.List<Integer> resultIndices = new java.util.ArrayList<>();

        for (int colIdx = 0; colIdx < data.numCols; colIdx++) {
            int col1 = colIndices[colIdx];
            int col2 = colIndices[colIdx + 1];

            for (int i = col1; i < col2; i++) {
                if (values[i] >= 0) {
                    resultIndices.add(rowIndices[i] * data.numCols + colIdx);
                }
            }
        }

        BaseMatrix res = createNewInstance(resultIndices.size(), 1, resultIndices.size());
        for (int i = 0; i < resultIndices.size(); i++) {
            res.set(i, 0, resultIndices.get(i));
        }
        return res;
    }

    /**
     * Returns the value at the specified matrix position.
     *
     * @param row the row index
     * @param col the column index
     * @return the value at position (row, col)
     */
    public double get(int row, int col) {
        try {
            return getData().get(row, col);
        } catch (IllegalArgumentException e) {
            // Handle case where position is not explicitly stored in sparse matrix
            // Sparse matrices don't store zero values, so accessing an unset position should return 0
            if (e.getMessage() != null && e.getMessage().contains("Outside of matrix bounds")) {
                return 0.0;
            }
            throw e;
        }
    }

    protected int getColumnIndex(int col) {
        return data.col_idx[col];
    }

    protected int[] getColumnIndices() {
        return data.col_idx;
    }

    protected int[] getColumnIndicesArray() {
        return data.col_idx;
    }

    /**
     * Returns the underlying EJML sparse matrix data structure.
     * This method provides direct access to the internal data for subclasses.
     *
     * @return the underlying {@link DMatrixSparseCSC} instance
     */
    public DMatrixSparseCSC getData() {
        return data;
    }

    /**
     * Sets the underlying EJML sparse matrix data structure.
     * This method allows subclasses to replace the internal data.
     *
     * @param newData the new {@link DMatrixSparseCSC} instance to use
     */
    public void setData(DMatrixSparseCSC newData) {
        this.data = newData;
    }

    /**
     * Sets the underlying data structure for BaseMatrix compatibility.
     *
     * @param newData the new data structure (must be DMatrixSparseCSC)
     */
    @Override
    protected void setData(Object newData) {
        if (newData instanceof DMatrixSparseCSC) {
            this.data = (DMatrixSparseCSC) newData;
        } else {
            throw new IllegalArgumentException("SparseMatrix requires DMatrixSparseCSC data");
        }
    }

    public int getNonZeroLength() {
        return data.nz_length;
    }

    protected void setNonZeroLength(int length) {
        data.nz_length = length;
    }

    protected int getNonZeroRow(int index) {
        return data.nz_rows[index];
    }

    public int[] getNonZeroRows() {
        return data.nz_rows;
    }

    protected double getNonZeroValue(int index) {
        return data.nz_values[index];
    }

    // Access to internal arrays for specialized operations
    public double[] getNonZeroValues() {
        return data.nz_values;
    }

    // Methods that access nz_length directly

    /**
     * Returns the number of non-zero elements in the matrix.
     *
     * @return the number of non-zero elements
     */
    public int getNonZeros() {
        return getData().nz_length;
    }

    /**
     * Returns the number of columns in the matrix.
     *
     * @return the number of columns
     */
    public int getNumCols() {
        return getData().getNumCols();
    }

    protected int getNumColsInternal() {
        return data.numCols;
    }

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
    public int getNumRows() {
        return getData().getNumRows();
    }

    // Additional field access methods
    protected int getNumRowsInternal() {
        return data.numRows;
    }

    /**
     * Checks if the matrix contains duplicate values.
     *
     * @return true if duplicates exist
     */
    public boolean hasDuplicates() {
        double[] values = data.nz_values;
        for (int i = 0; i < data.nz_length - 1; i++) {
            for (int j = i + 1; j < data.nz_length; j++) {
                if (values[i] == values[j]) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Checks if the matrix contains any finite (non-infinite, non-NaN) values.
     *
     * @return true if at least one finite value exists
     */
    public boolean hasFinite() {
        for (int i = 0; i < this.data.nz_length; i++) {
            if (!Utils.isInf(this.data.nz_values[i]) && !Double.isNaN(this.data.nz_values[i])) return true;
        }
        return false;
    }

    /**
     * Checks if the matrix contains any infinite values.
     *
     * @return true if at least one infinite value exists
     */
    public boolean hasInfinite() {
        for (int i = 0; i < this.data.nz_length; i++) {
            if (Double.isInfinite(this.data.nz_values[i])) return true;
        }
        return false;
    }

    /**
     * Checks if more than one finite value exists.
     *
     * @return true if more than one finite value exists
     */
    public boolean hasMultipleFinite() {
        boolean hasOneFinite = false;
        for (int i = 0; i < this.data.nz_length; i++) {
            if (!Utils.isInf(this.data.nz_values[i]) && !Double.isNaN(this.data.nz_values[i])) {
                if (hasOneFinite) return true;
                hasOneFinite = true;
            }
        }
        return false;
    }

    /**
     * Checks if the matrix contains any NaN values.
     *
     * @return true if at least one NaN exists
     */
    public boolean hasNaN() {
        for (int i = 0; i < this.data.nz_length; i++) {
            if (Double.isNaN(this.data.nz_values[i])) return true;
        }
        return false;
    }

    /**
     * Performs matrix multiplication: output = A * B.
     *
     * @param A      the left matrix
     * @param B      the right matrix
     * @param output the result matrix
     */
    protected void mult(DMatrix A, DMatrix B, DMatrix output) {
        CommonOps_DSCC.mult((DMatrixSparseCSC) A, (DMatrixSparseCSC) B, (DMatrixSparseCSC) output);
    }

    // Additional mult methods (complement existing mult method)
    protected DMatrix multMatrix(DMatrix B, DMatrix output) {
        return CommonOps_DSCC.mult(this.getData(), (DMatrixSparseCSC) B, (DMatrixSparseCSC) output);
    }

    /**
     * Removes zero-valued elements from the sparse matrix with default tolerance.
     */
    protected void removeZeros() {
        CommonOps_DSCC.removeZeros(getData(), 0);
    }

    // Additional removeZeros methods (complement existing removeZeros method)
    protected void removeZerosWithTol(double tolerance) {
        CommonOps_DSCC.removeZeros(this.getData(), tolerance);
    }

    /**
     * Reshapes the matrix to the specified dimensions.
     *
     * @param numRows the new number of rows
     * @param numCols the new number of columns
     */
    public void reshape(int numRows, int numCols) {
        getData().reshape(numRows, numCols);
    }

    /**
     * Scales the matrix in-place by the specified factor.
     *
     * @param alpha the scaling factor
     */
    protected void scaleInPlace(double alpha) {
        CommonOps_DSCC.scale(alpha, getData(), getData());
    }

    // Additional scale methods (complement existing scaleInPlace method)
    protected void scaleMatrix(double scalar, DMatrix output) {
        CommonOps_DSCC.scale(scalar, this.getData(), (DMatrixSparseCSC) output);
    }

    /**
     * Sets the value at the specified matrix position.
     *
     * @param row   the row index
     * @param col   the column index
     * @param value the value to set
     */
    public void set(int row, int col, double value) {
        getData().set(row, col, value);
    }

    protected void setColumnIndex(int col, int value) {
        data.col_idx[col] = value;
    }

    // Methods that perform operations on internal sparse structure

    protected void setNonZeroRow(int index, int row) {
        data.nz_rows[index] = row;
    }

    protected void setNonZeroValue(int index, double value) {
        data.nz_values[index] = value;
    }

    /**
     * Reduce the maximum number of columns by setting the internal column count.
     *
     * @param newmax the new maximum number of columns
     */
    public void shrinkNumCols(int newmax) {
        data.numCols = newmax;
    }

    /**
     * Reduce the maximum number of rows by setting the internal row count.
     *
     * @param newmax the new maximum number of rows
     */
    public void shrinkNumRows(int newmax) {
        data.numRows = newmax;
    }

    // SumCols method
    protected DMatrix sumColsRaw() {
        return CommonOps_DSCC.sumCols(this.getData(), null);
    }

    // SumRows method
    protected DMatrix sumRowsRaw() {
        return CommonOps_DSCC.sumRows(this.getData(), null);
    }

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
     * Computes the transpose of the input matrix.
     *
     * @param input  the matrix to transpose
     * @param output the transposed matrix
     */
    protected void transpose(DMatrix input, DMatrix output) {
        CommonOps_DSCC.transpose((DMatrixSparseCSC) input, (DMatrixSparseCSC) output, null);
    }

    // Transpose method (complement existing transpose method)
    protected void transposeMatrix(DMatrix output) {
        CommonOps_DSCC.transpose(this.getData(), (DMatrixSparseCSC) output, null);
    }

    /**
     * Sets all elements in the matrix to zero.
     */
    public void zero() {
        getData().zero();
    }
}