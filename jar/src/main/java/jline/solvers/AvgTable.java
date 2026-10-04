/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.List;

/**
 * Abstract base class for representing average performance metrics tables.
 *
 * <p>AvgTable provides a unified interface for storing and accessing
 * performance metrics computed by various solvers in LINE. The table
 * stores metrics such as throughput, utilization, queue length, and
 * response time organized by stations and job classes.
 *
 * <p>Subclasses implement specific table formats for different types
 * of performance metrics and provide specialized printing methods
 * for displaying results in human-readable formats.
 *
 * <p>The underlying data is stored as a Matrix and can be accessed
 * either column-wise or as a complete matrix. The class includes
 * utilities for sanitizing numerical results to handle small
 * floating-point perturbations.
 *
 * @see NetworkAvgTable
 * @see LayeredNetworkAvgTable
 * @see Matrix
 */
public abstract class AvgTable {

    /**
     * Solver options used when computing this table
     */
    SolverOptions options;

    /**
     * Matrix containing the average performance metrics data
     */
    Matrix T;

    /**
     * Constructs an AvgTable from a Matrix.
     *
     * @param table the matrix containing performance metrics data
     */
    public AvgTable(Matrix table) {
        this.T = new Matrix(table);
    }

    /**
     * Constructs an AvgTable from a list of lists.
     *
     * @param table the data as a list of lists, where each inner list represents a row
     */
    public AvgTable(ArrayList<List<Double>> table) {
        this.T = new Matrix(table);
    }

    /**
     * Retrieves a specific column from the table as a list.
     *
     * @param col the column index to retrieve
     * @return the column data as a List of Double values
     */
    public List<Double> get(int col) {
        return this.T.getColumn(col).toList1D();
    }

    /**
     * Returns the underlying matrix data.
     *
     * @return the Matrix containing the performance metrics
     */
    public Matrix getData() {
        return this.T;
    }

    /**
     * Prints the table contents in a human-readable format.
     * Implementation is provided by subclasses for specific table types.
     */
    public abstract void print();

    /**
     * Sanitizes the table data to fix small numerical perturbations.
     *
     * <p>This method rounds values that are very close to integers to their
     * nearest integer values, and sets very small values to zero to eliminate
     * numerical noise that can accumulate during floating-point computations.
     */
    public void sanitize() {
        double val = 0.0;
        double rounded_val = 0.0;
        for (int i = 0; i < T.getNumRows(); i++)
            for (int j = 0; j < T.getNumCols(); j++) {
                val = T.get(i, j) * 10;
                rounded_val = FastMath.round(val);
                if (FastMath.abs(val - rounded_val) < GlobalConstants.CoarseTol * val) {
                    T.set(i, j, rounded_val / 10);
                }
            }
    }

    /**
     * Sets the solver options associated with this table.
     *
     * @param options the solver options used when computing this table
     */
    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    /**
     * Returns the table data as a Matrix after sanitizing numerical values.
     *
     * @return the sanitized Matrix containing the performance metrics
     */
    public Matrix toMatrix() {
        sanitize();
        return T;
    }

    /**
     * Formats a metric value compactly with the given number of significant
     * digits, matching the MATLAB table display (e.g. 0, 0.5, 1, 0.33333,
     * 1.5e-06). Integers print without a fractional part.
     *
     * @param value  the value to format
     * @param digits number of significant digits
     * @return the formatted value
     */
    protected static String fmtValue(double value, int digits) {
        if (Double.isNaN(value)) {
            return "NaN";
        }
        if (Double.isInfinite(value)) {
            return value > 0 ? "Inf" : "-Inf";
        }
        if (value == java.lang.Math.rint(value) && java.lang.Math.abs(value) < 1e10) {
            return Long.toString((long) value);
        }
        String s = String.format(java.util.Locale.ROOT, "%." + digits + "g", value);
        String mantissa = s;
        String exponent = "";
        int e = s.indexOf('e');
        if (e < 0) {
            e = s.indexOf('E');
        }
        if (e >= 0) {
            mantissa = s.substring(0, e);
            exponent = s.substring(e);
        }
        if (mantissa.indexOf('.') >= 0) {
            int last = mantissa.length();
            while (last > 0 && mantissa.charAt(last - 1) == '0') {
                last--;
            }
            if (last > 0 && mantissa.charAt(last - 1) == '.') {
                last--;
            }
            mantissa = mantissa.substring(0, last);
        }
        return mantissa + exponent;
    }

    /**
     * Prints a header row followed by data rows, sizing every column to its
     * widest cell with a two-space gap and no trailing padding, as in the
     * MATLAB and Python table displays. Name columns are left-aligned and
     * numeric columns (all data cells numeric) are right-aligned, following
     * the MATLAB table display convention.
     *
     * @param headers the column headers
     * @param rows    the data rows, each with one cell per header
     */
    protected static void printFormattedTable(String[] headers, List<String[]> rows) {
        int ncols = headers.length;
        int[] width = new int[ncols];
        boolean[] rightAlign = new boolean[ncols];
        for (int j = 0; j < ncols; j++) {
            width[j] = headers[j].length();
            rightAlign[j] = !rows.isEmpty();
        }
        for (String[] row : rows) {
            for (int j = 0; j < ncols; j++) {
                if (row[j].length() > width[j]) {
                    width[j] = row[j].length();
                }
                if (!isNumericCell(row[j])) {
                    rightAlign[j] = false;
                }
            }
        }
        StringBuilder sb = new StringBuilder();
        appendFormattedRow(sb, headers, width, rightAlign);
        for (String[] row : rows) {
            appendFormattedRow(sb, row, width, rightAlign);
        }
        System.out.print(sb);
        System.out.flush();
    }

    private static boolean isNumericCell(String cell) {
        return cell.matches("-?(\\d+(\\.\\d+)?([eE][+-]?\\d+)?|NaN|Inf)");
    }

    private static void appendFormattedRow(StringBuilder sb, String[] cells, int[] width, boolean[] rightAlign) {
        for (int j = 0; j < cells.length; j++) {
            if (rightAlign[j]) {
                for (int p = cells[j].length(); p < width[j]; p++) {
                    sb.append(' ');
                }
                sb.append(cells[j]);
                if (j < cells.length - 1) {
                    sb.append("  ");
                }
            } else {
                sb.append(cells[j]);
                if (j < cells.length - 1) {
                    for (int p = cells[j].length(); p < width[j] + 2; p++) {
                        sb.append(' ');
                    }
                }
            }
        }
        sb.append('\n');
    }

}
