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

}
