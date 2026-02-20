/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util.matrix;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;

import jline.util.Utils;
import org.ejml.data.DMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;

/**
 * Base class for dense matrix implementations, containing the core data structure
 * and methods that directly manipulate the underlying dense matrix representation.
 *
 * <p>This class encapsulates all interactions with EJML's {@link DMatrixRMaj} and
 * {@link CommonOps_DDRM} classes, providing a clean abstraction layer for dense matrix
 * operations. It serves as the foundation for concrete dense matrix implementations.</p>
 *
 * <p>The class provides:</p>
 * <ul>
 *   <li>Basic matrix operations (get, set, dimensions)</li>
 *   <li>Encapsulated CommonOps_DDRM method calls</li>
 *   <li>Memory management for dense matrix data structures</li>
 *   <li>Utility methods for matrix manipulation</li>
 * </ul>
 *
 * <p>Dense matrices store all elements (including zeros) explicitly, making them
 * suitable for matrices with a high percentage of non-zero elements or when
 * dense linear algebra operations are required.</p>
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */
public abstract class DenseMatrix extends BaseMatrix {

    /**
     * The underlying EJML dense matrix data structure
     */
    private DMatrixRMaj data;

    /**
     * Constructs a dense matrix with specified dimensions.
     * All elements are initialized to zero.
     *
     * @param numRows the number of rows in the matrix
     * @param numCols the number of columns in the matrix
     */
    public DenseMatrix(int numRows, int numCols) {
        setData(new DMatrixRMaj(numRows, numCols));
    }

    /**
     * Constructs a dense matrix by copying an existing EJML dense matrix.
     *
     * @param matrix the EJML dense matrix to copy
     */
    public DenseMatrix(DMatrixRMaj matrix) {
        setData(new DMatrixRMaj(matrix));
    }

    /**
     * Constructs a dense matrix by copying an existing DMatrix.
     *
     * @param matrix the DMatrix to copy
     */
    public DenseMatrix(DMatrix matrix) {
        setData(new DMatrixRMaj((DMatrixRMaj) matrix));
    }

    /**
     * Copy constructor for creating a new dense matrix from an existing one.
     *
     * @param matrix the dense matrix to copy
     */
    protected DenseMatrix(DenseMatrix matrix) {
        setData(new DMatrixRMaj(matrix.getData()));
    }

    /**
     * Protected constructor for delayed initialization.
     * Subclasses must ensure that data is properly initialized.
     */
    protected DenseMatrix() {
        // data will be initialized later by subclass
    }

    protected static void addMatricesInPlaceStatic(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output) {
        CommonOps_DDRM.add(alpha, (DMatrixRMaj) A, beta, (DMatrixRMaj) B, (DMatrixRMaj) output);
    }

    protected static DMatrixRMaj addMatricesStatic(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output) {
        return CommonOps_DDRM.add(alpha, (DMatrixRMaj) A, beta, (DMatrixRMaj) B, (DMatrixRMaj) output);
    }

    /**
     * Performs element-wise addition of two dense matrices: output = A + B.
     *
     * @param A      the first matrix
     * @param B      the second matrix
     * @param output the result matrix (can be null to create a new matrix)
     * @return the sum matrix
     */
    protected static DMatrixRMaj addStatic(DMatrixRMaj A, DMatrixRMaj B, DMatrixRMaj output) {
        return CommonOps_DDRM.add(A, B, output);
    }

    protected static void changeSignStatic(DMatrix input, DMatrix output) {
        CommonOps_DDRM.changeSign((DMatrixRMaj) input, (DMatrixRMaj) output);
    }

    protected static double determinantStatic(DMatrix matrix) {
        return CommonOps_DDRM.det((DMatrixRMaj) matrix);
    }

    protected static DMatrixRMaj elementMultStatic(DMatrix A, DMatrix B, DMatrix output) {
        return CommonOps_DDRM.elementMult((DMatrixRMaj) A, (DMatrixRMaj) B, (DMatrixRMaj) output);
    }

    // Static utility methods
    protected static void extractMatrixStatic(DMatrix src, int srcX0, int srcX1, int srcY0, int srcY1, DMatrix dst, int dstY0, int dstX0) {
        CommonOps_DDRM.extract(src, srcX0, srcX1, srcY0, srcY1, dst, dstY0, dstX0);
    }

    protected static void fillMatrixStatic(DMatrix matrix, double value) {
        CommonOps_DDRM.fill((DMatrixRMaj) matrix, value);
    }

    /**
     * Creates an identity matrix of the specified size.
     * An identity matrix has ones on the main diagonal and zeros elsewhere.
     *
     * @param size the number of rows and columns in the identity matrix
     * @return a new identity matrix
     */
    protected static DMatrixRMaj identityStatic(int size) {
        return CommonOps_DDRM.identity(size);
    }

    protected static void invertStatic(DMatrix input, DMatrix output) {
        CommonOps_DDRM.invert((DMatrixRMaj) input, (DMatrixRMaj) output);
    }

    // CommonOps_DDRM encapsulation methods

    /**
     * Computes the Kronecker product of two dense matrices: A \otimes B.
     * The Kronecker product creates a larger matrix where each element of A
     * is multiplied by the entire matrix B.
     *
     * @param A      the left matrix
     * @param B      the right matrix
     * @param output the result matrix (can be null to create a new matrix)
     * @return the Kronecker product matrix
     */
    protected static DMatrixRMaj kronStatic(DMatrixRMaj A, DMatrixRMaj B, DMatrixRMaj output) {
        return CommonOps_DDRM.kron(A, B, output);
    }

    protected static DMatrixRMaj multMatrixStatic(DMatrix A, DMatrix B, DMatrix output) {
        return CommonOps_DDRM.mult((DMatrixRMaj) A, (DMatrixRMaj) B, (DMatrixRMaj) output);
    }

    protected static void scaleMatrixStatic(double scalar, DMatrix input, DMatrix output) {
        CommonOps_DDRM.scale(scalar, (DMatrixRMaj) input, (DMatrixRMaj) output);
    }

    // BaseMatrix implementation for dense matrices

    protected static void transposeMatrixStatic(DMatrix input, DMatrix output) {
        CommonOps_DDRM.transpose((DMatrixRMaj) input, (DMatrixRMaj) output);
    }

    @Override
    public void absEq() {
        // Direct array access is much faster than get/set method calls
        double[] dataArray = getData().data;
        for (int i = 0; i < dataArray.length; i++) {
            if (dataArray[i] < 0) {
                dataArray[i] = -dataArray[i];
            }
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
        CommonOps_DDRM.add(alpha, (DMatrixRMaj) A, beta, (DMatrixRMaj) B, (DMatrixRMaj) output);
    }

    @Override
    protected DMatrix addMatrices(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output) {
        return CommonOps_DDRM.add(alpha, (DMatrixRMaj) A, beta, (DMatrixRMaj) B, (DMatrixRMaj) output);
    }

    // Additional clean methods to match SparseMatrix interface
    protected void addMatricesInPlace(double alpha, DMatrix A, double beta, DMatrix B, DMatrix output) {
        CommonOps_DDRM.add(alpha, (DMatrixRMaj) A, beta, (DMatrixRMaj) B, (DMatrixRMaj) output);
    }

    @Override
    public boolean any() {
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                if (get(i, j) != 0.0) {
                    return true;
                }
            }
        }
        return false;
    }

    // Stub implementations for sparse-specific methods
    @Override
    protected void applyConditionalTransform(double source, double target, double tol, String operation) {
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                double val = get(i, j);
                boolean replace = false;

                switch (operation) {
                    case "NaN":
                        replace = Double.isNaN(val);
                        break;
                    case "Inf":
                        replace = Double.isInfinite(val);
                        break;
                    case "Value":
                        replace = Math.abs(val - source) < tol;
                        break;
                }

                if (replace) {
                    set(i, j, target);
                }
            }
        }
    }

    @Override
    protected void changeSign() {
        // Direct array access is much faster than get/set method calls
        double[] dataArray = getData().data;
        for (int i = 0; i < dataArray.length; i++) {
            dataArray[i] = -dataArray[i];
        }
    }

    /**
     * Creates and returns a deep copy of this dense matrix.
     *
     * @return a new dense matrix that is a copy of this instance
     */
    public abstract BaseMatrix copy();

    @Override
    public void colIncrease(int col, double a) {
        for (int i = 0; i < getNumRows(); i++) {
            set(i, col, get(i, col) + a);
        }
    }

    @Override
    protected void concatColumnsInPlace(DMatrix left, DMatrix right, DMatrix output) {
        CommonOps_DDRM.concatColumns((DMatrixRMaj) left, (DMatrixRMaj) right, (DMatrixRMaj) output);
    }

    @Override
    protected void concatRowsInPlace(DMatrix top, DMatrix bottom, DMatrix output) {
        CommonOps_DDRM.concatRows((DMatrixRMaj) top, (DMatrixRMaj) bottom, (DMatrixRMaj) output);
    }

    @Override
    public BaseMatrix countEachRow(double val) {
        BaseMatrix res = createNewInstance(getNumRows(), 1, getNumRows());

        for (int i = 0; i < getNumRows(); i++) {
            int count = 0;
            for (int j = 0; j < getNumCols(); j++) {
                if (get(i, j) == val) {
                    count++;
                }
            }
            res.set(i, 0, count);
        }
        return res;
    }

    /**
     * Creates a new instance of the concrete dense matrix implementation.
     * This factory method allows the base class to create instances of the correct subtype.
     *
     * @param rows the number of rows for the new matrix
     * @param cols the number of columns for the new matrix
     * @return a new instance of the concrete dense matrix type
     */
    protected abstract BaseMatrix createNewInstance(int rows, int cols, int nzLength);

    @Override
    protected double determinant() {
        return CommonOps_DDRM.det(this.getData());
    }

    @Override
    protected void divideInPlace(double scalar) {
        // Direct array access is much faster than get/set method calls
        double[] dataArray = getData().data;
        for (int i = 0; i < dataArray.length; i++) {
            dataArray[i] /= scalar;
        }
    }

    @Override
    protected void divideMatrix(double scalar, DMatrix output) {
        CommonOps_DDRM.divide(scalar, this.getData(), (DMatrixRMaj) output);
    }

    @Override
    protected void divideRowsByArray(double[] diag, int offset) {
        for (int i = 0; i < getNumRows() && i + offset < diag.length; i++) {
            for (int j = 0; j < getNumCols(); j++) {
                set(i, j, get(i, j) / diag[i + offset]);
            }
        }
    }

    @Override
    protected void divideScalarByMatrix(double scalar, DMatrix output) {
        CommonOps_DDRM.divide(scalar, this.getData(), (DMatrixRMaj) output);
    }

    @Override
    protected double elementMax() {
        double max = NegInf;
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                double val = get(i, j);
                if (val > max) {
                    max = val;
                }
            }
        }
        return max;
    }

    @Override
    protected double elementMaxAbs() {
        double maxAbs = 0.0;
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                double absVal = Math.abs(get(i, j));
                if (absVal > maxAbs) {
                    maxAbs = absVal;
                }
            }
        }
        return maxAbs;
    }

    @Override
    protected double elementMin() {
        double min = Inf;
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                double val = get(i, j);
                if (val < min) {
                    min = val;
                }
            }
        }
        return min;
    }

    @Override
    protected DMatrix elementMult(DMatrix A, DMatrix B, DMatrix output) {
        CommonOps_DDRM.elementMult((DMatrixRMaj) A, (DMatrixRMaj) B, (DMatrixRMaj) output);
        return output;
    }

    @Override
    protected double elementSum() {
        // Direct array access is much faster than get(i, j) method calls
        double[] dataArray = getData().data;
        double sum = 0.0;
        for (int i = 0; i < dataArray.length; i++) {
            sum += dataArray[i];
        }
        return sum;
    }

    /**
     * Indicates whether some other object is "equal to" this dense matrix.
     *
     * @param obj the reference object with which to compare
     * @return {@code true} if this object is the same as the obj argument; {@code false} otherwise
     */
    public abstract boolean equals(Object obj);

    @Override
    protected void extractMatrix(DMatrix src, int srcX0, int srcX1, int srcY0, int srcY1, DMatrix dst, int dstY0, int dstX0) {
        CommonOps_DDRM.extract(src, srcX0, srcX1, srcY0, srcY1, dst, dstY0, dstX0);
    }

    @Override
    protected void fillMatrix(double value) {
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                set(i, j, value);
            }
        }
    }

    @Override
    public BaseMatrix findNonNegative() {
        java.util.List<Integer> resultIndices = new java.util.ArrayList<>();

        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                if (get(i, j) >= 0) {
                    resultIndices.add(i * getNumCols() + j);
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
        return getData().get(row, col);
    }

    @Override
    protected int getColumnIndex(int col) {
        throw new UnsupportedOperationException("Dense matrix does not use column index arrays");
    }

    @Override
    protected int[] getColumnIndicesArray() {
        throw new UnsupportedOperationException("Dense matrix does not use column index arrays");
    }

    /**
     * Returns the underlying EJML dense matrix data structure.
     * This method provides direct access to the internal data for subclasses.
     *
     * @return the underlying {@link DMatrixRMaj} instance
     */
    protected DMatrixRMaj getData() {
        return data;
    }

    /**
     * Sets the underlying EJML dense matrix data structure.
     * This method allows subclasses to replace the internal data.
     *
     * @param newData the new {@link DMatrixRMaj} instance to use
     */
    protected void setData(DMatrixRMaj newData) {
        this.data = newData;
    }

    /**
     * Sets the underlying EJML dense matrix data structure.
     * This method allows subclasses to replace the internal data.
     *
     * @param newData the new DMatrix instance to use
     */
    public void setData(DMatrix newData) {
        this.data = (DMatrixRMaj) newData;
    }

    @Override
    protected void setData(Object newData) {
        if (newData instanceof DMatrixRMaj) {
            this.data = (DMatrixRMaj) newData;
        } else {
            throw new IllegalArgumentException("DenseMatrix requires DMatrixRMaj data");
        }
    }

    @Override
    public int getNonZeroLength() {
        return getNonZeros();
    }

    @Override
    protected void setNonZeroLength(int length) {
        throw new UnsupportedOperationException("Dense matrix does not use non-zero length");
    }

    @Override
    protected int getNonZeroRow(int index) {
        int[] rows = getNonZeroRows();
        if (index < rows.length) {
            return rows[index];
        }
        return -1;
    }

    @Override
    public int[] getNonZeroRows() {
        // First pass: count non-zeros for exact array allocation
        double[] dataArray = getData().data;
        int numRows = getNumRows();
        int numCols = getNumCols();
        int count = 0;
        for (int i = 0; i < dataArray.length; i++) {
            if (dataArray[i] != 0.0) {
                count++;
            }
        }

        // Second pass: fill result array using direct array access
        int[] result = new int[count];
        int idx = 0;
        for (int i = 0; i < numRows; i++) {
            int rowOffset = i * numCols;
            for (int j = 0; j < numCols; j++) {
                if (dataArray[rowOffset + j] != 0.0) {
                    result[idx++] = i;
                }
            }
        }
        return result;
    }

    // Sparse indexing methods that don't apply to dense matrices
    @Override
    protected double getNonZeroValue(int index) {
        double[] values = getNonZeroValues();
        if (index < values.length) {
            return values[index];
        }
        return 0.0;
    }

    // Implementation of abstract methods from BaseMatrix using DMatrixRMaj operations

    @Override
    public double[] getNonZeroValues() {
        // First pass: count non-zeros for exact array allocation
        double[] dataArray = getData().data;
        int count = 0;
        for (int i = 0; i < dataArray.length; i++) {
            if (dataArray[i] != 0.0) {
                count++;
            }
        }

        // Second pass: fill result array using direct array access
        double[] result = new double[count];
        int idx = 0;
        for (int i = 0; i < dataArray.length; i++) {
            if (dataArray[i] != 0.0) {
                result[idx++] = dataArray[i];
            }
        }
        return result;
    }

    @Override
    public int getNonZeros() {
        // Direct array access is much faster than get(i, j) method calls
        double[] dataArray = getData().data;
        int count = 0;
        for (int i = 0; i < dataArray.length; i++) {
            if (dataArray[i] != 0.0) {
                count++;
            }
        }
        return count;
    }

    /**
     * Returns the number of columns in the matrix.
     *
     * @return the number of columns
     */
    public int getNumCols() {
        return getData().getNumCols();
    }

    /**
     * Returns the number of rows in the matrix.
     *
     * @return the number of rows
     */
    public int getNumRows() {
        return getData().getNumRows();
    }

    @Override
    public boolean hasDuplicates() {
        for (int i1 = 0; i1 < getNumRows(); i1++) {
            for (int j1 = 0; j1 < getNumCols(); j1++) {
                double val1 = get(i1, j1);
                for (int i2 = 0; i2 < getNumRows(); i2++) {
                    for (int j2 = 0; j2 < getNumCols(); j2++) {
                        if (i1 != i2 || j1 != j2) {
                            if (get(i2, j2) == val1) {
                                return true;
                            }
                        }
                    }
                }
            }
        }
        return false;
    }

    @Override
    public boolean hasFinite() {
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                if (!Utils.isInf(get(i, j)) && !Double.isNaN(get(i, j))) {
                    return true;
                }
            }
        }
        return false;
    }

    @Override
    public boolean hasInfinite() {
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                if (Utils.isInf(get(i, j))) {
                    return true;
                }
            }
        }
        return false;
    }

    @Override
    public boolean hasMultipleFinite() {
        boolean hasOneFinite = false;
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                if (!Utils.isInf(get(i, j)) && !Double.isNaN(get(i, j))) {
                    if (hasOneFinite) return true;
                    hasOneFinite = true;
                }
            }
        }
        return false;
    }

    @Override
    public boolean hasNaN() {
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                if (Double.isNaN(get(i, j))) {
                    return true;
                }
            }
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
        CommonOps_DDRM.mult((DMatrixRMaj) A, (DMatrixRMaj) B, (DMatrixRMaj) output);
    }

    // Additional methods to match SparseMatrix interface

    @Override
    protected DMatrix multMatrix(DMatrix B, DMatrix output) {
        CommonOps_DDRM.mult(this.getData(), (DMatrixRMaj) B, (DMatrixRMaj) output);
        return output;
    }

    @Override
    protected void removeZeros() {
        // No-op for dense matrices
    }

    @Override
    protected void removeZerosWithTol(double tolerance) {
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                if (Math.abs(get(i, j)) < tolerance) {
                    set(i, j, 0.0);
                }
            }
        }
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

    @Override
    protected void scaleInPlace(double alpha) {
        // Direct array access is much faster than get/set method calls
        double[] dataArray = getData().data;
        for (int i = 0; i < dataArray.length; i++) {
            dataArray[i] *= alpha;
        }
    }

    @Override
    protected void scaleMatrix(double scalar, DMatrix output) {
        CommonOps_DDRM.scale(scalar, this.getData(), (DMatrixRMaj) output);
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

    @Override
    protected void setColumnIndex(int col, int value) {
        throw new UnsupportedOperationException("Dense matrix does not use column index arrays");
    }

    @Override
    protected void setNonZeroRow(int index, int row) {
        throw new UnsupportedOperationException("Dense matrix does not support sparse-style indexing");
    }

    @Override
    protected void setNonZeroValue(int index, double value) {
        throw new UnsupportedOperationException("Dense matrix does not support sparse-style indexing");
    }

    @Override
    public void shrinkNumCols(int newmax) {
        if (newmax < getNumCols()) {
            reshape(getNumRows(), newmax);
        }
    }

    @Override
    public void shrinkNumRows(int newmax) {
        if (newmax < getNumRows()) {
            reshape(newmax, getNumCols());
        }
    }

    @Override
    protected DMatrix sumColsRaw() {
        DMatrixRMaj result = new DMatrixRMaj(1, getNumCols());
        for (int j = 0; j < getNumCols(); j++) {
            double sum = 0.0;
            for (int i = 0; i < getNumRows(); i++) {
                sum += get(i, j);
            }
            result.set(0, j, sum);
        }
        return result;
    }

    @Override
    protected DMatrix sumRowsRaw() {
        DMatrixRMaj result = new DMatrixRMaj(getNumRows(), 1);
        for (int i = 0; i < getNumRows(); i++) {
            double sum = 0.0;
            for (int j = 0; j < getNumCols(); j++) {
                sum += get(i, j);
            }
            result.set(i, 0, sum);
        }
        return result;
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
        CommonOps_DDRM.transpose((DMatrixRMaj) input, (DMatrixRMaj) output);
    }

    @Override
    protected void transposeMatrix(DMatrix output) {
        CommonOps_DDRM.transpose(this.getData(), (DMatrixRMaj) output);
    }

    /**
     * Sets all elements in the matrix to zero.
     */
    public void zero() {
        getData().zero();
    }
}