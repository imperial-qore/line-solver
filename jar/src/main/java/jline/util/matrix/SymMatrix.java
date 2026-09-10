/**
 * @file A dense matrix of exact symbolic scalars.
 *
 * @since LINE 3.0
 */
package jline.util.matrix;

import jline.util.symbolic.SymContext;
import jline.util.symbolic.SymExpr;

/**
 * A dense matrix whose entries are {@link SymExpr} rational functions.
 *
 * <p>The SYMBOLIC counterpart of {@link ComplexMatrix}, and it exists for the
 * same reason that one does: {@link Matrix} stores doubles, so a scalar that is
 * not a double needs its own container rather than a reinterpretation of that
 * one. The surface here is deliberately only what the pfqn_gld symbolic arm
 * reads -- shape, element access, row and column slicing -- and not a linear
 * algebra API. Nothing in that arm multiplies two matrices.
 *
 * <p>Every entry shares ONE {@link SymContext}, which is what makes the entries
 * addable: Rings fixes the indeterminates when the ring is built.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class SymMatrix {

    private final SymContext ctx;
    private final int rows;
    private final int cols;
    private final SymExpr[] data;

    /**
     * A zero-filled matrix.
     *
     * @param ctx  the shared symbolic context
     * @param rows row count
     * @param cols column count
     */
    public SymMatrix(SymContext ctx, int rows, int cols) {
        if (rows < 0 || cols < 0) {
            throw new RuntimeException("SymMatrix: negative dimensions.");
        }
        this.ctx = ctx;
        this.rows = rows;
        this.cols = cols;
        this.data = new SymExpr[rows * cols];
        SymExpr zero = ctx.zero();
        for (int i = 0; i < this.data.length; i++) {
            this.data[i] = zero;
        }
    }

    /**
     * A matrix of the given numeric values, carried in exactly.
     *
     * @param ctx    the shared symbolic context
     * @param values row-major values
     * @return the matrix
     */
    public static SymMatrix of(SymContext ctx, double[][] values) {
        int r = values.length;
        int c = r == 0 ? 0 : values[0].length;
        SymMatrix out = new SymMatrix(ctx, r, c);
        for (int i = 0; i < r; i++) {
            if (values[i].length != c) {
                throw new RuntimeException("SymMatrix: ragged input.");
            }
            for (int j = 0; j < c; j++) {
                out.set(i, j, ctx.constant(values[i][j]));
            }
        }
        return out;
    }

    /** @return the shared symbolic context */
    public SymContext context() {
        return this.ctx;
    }

    /** @return the row count */
    public int getNumRows() {
        return this.rows;
    }

    /** @return the column count */
    public int getNumCols() {
        return this.cols;
    }

    /** @return the element count */
    public int getNumElements() {
        return this.data.length;
    }

    /** @return whether the matrix holds no elements */
    public boolean isEmpty() {
        return this.data.length == 0;
    }

    /**
     * @param i row
     * @param j column
     * @return the entry
     */
    public SymExpr get(int i, int j) {
        checkIndex(i, j);
        return this.data[i * this.cols + j];
    }

    /**
     * Linear access in row-major order, as {@link ComplexMatrix#get(int)} offers.
     *
     * @param idx linear index
     * @return the entry
     */
    public SymExpr get(int idx) {
        if (idx < 0 || idx >= this.data.length) {
            throw new RuntimeException("SymMatrix: index " + idx + " is outside a matrix of "
                    + this.data.length + " elements.");
        }
        return this.data[idx];
    }

    /**
     * @param i   row
     * @param j   column
     * @param val the entry
     */
    public void set(int i, int j, SymExpr val) {
        checkIndex(i, j);
        if (val.context() != this.ctx) {
            throw new RuntimeException("SymMatrix: the value comes from a different symbolic context.");
        }
        this.data[i * this.cols + j] = val;
    }

    /**
     * Rows {@code [row0,row1)} as a new matrix.
     *
     * @param row0 first row, inclusive
     * @param row1 last row, exclusive
     * @return the slice
     */
    public SymMatrix extractRows(int row0, int row1) {
        if (row0 < 0 || row1 > this.rows || row0 > row1) {
            throw new RuntimeException("SymMatrix: row range [" + row0 + "," + row1 + ") is outside "
                    + this.rows + " rows.");
        }
        SymMatrix out = new SymMatrix(this.ctx, row1 - row0, this.cols);
        for (int i = row0; i < row1; i++) {
            for (int j = 0; j < this.cols; j++) {
                out.set(i - row0, j, get(i, j));
            }
        }
        return out;
    }

    /**
     * Columns {@code [col0,col1)} as a new matrix.
     *
     * @param col0 first column, inclusive
     * @param col1 last column, exclusive
     * @return the slice
     */
    public SymMatrix extractColumns(int col0, int col1) {
        if (col0 < 0 || col1 > this.cols || col0 > col1) {
            throw new RuntimeException("SymMatrix: column range [" + col0 + "," + col1 + ") is outside "
                    + this.cols + " columns.");
        }
        SymMatrix out = new SymMatrix(this.ctx, this.rows, col1 - col0);
        for (int i = 0; i < this.rows; i++) {
            for (int j = col0; j < col1; j++) {
                out.set(i, j - col0, get(i, j));
            }
        }
        return out;
    }

    /** @return an independent copy */
    public SymMatrix copy() {
        SymMatrix out = new SymMatrix(this.ctx, this.rows, this.cols);
        System.arraycopy(this.data, 0, out.data, 0, this.data.length);
        return out;
    }

    private void checkIndex(int i, int j) {
        if (i < 0 || i >= this.rows || j < 0 || j >= this.cols) {
            throw new RuntimeException("SymMatrix: (" + i + "," + j + ") is outside a "
                    + this.rows + "x" + this.cols + " matrix.");
        }
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                if (j > 0) {
                    sb.append("  ");
                }
                sb.append(get(i, j));
            }
            sb.append('\n');
        }
        return sb.toString();
    }
}
