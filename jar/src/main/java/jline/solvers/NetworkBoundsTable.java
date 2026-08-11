/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.VerboseLevel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Table of bound brackets produced by SolverBA, in the layout of
 * {@link NetworkAvgTable}.
 *
 * <p>Each row is a station-class pair carrying the lower and upper side of a
 * bound family: queue length ({@code Qlower}, {@code Qupper}) and throughput
 * ({@code Tlower}, {@code Tupper}). The two sides bracket the exact solution.
 *
 * <p>One-sided families (cub is upper-only, mbjb and ldbcmp are lower-only)
 * carry NaN on the missing side. NaN is preserved rather than replaced by zero,
 * so a reader can tell "no bound on this side" apart from "the bound is zero".
 *
 * @see NetworkAvgTable
 * @see jline.solvers.ba.SolverBA
 */
public class NetworkBoundsTable extends AvgTable {

    /**
     * Names of job classes, one entry per row
     */
    List<String> classNames;

    /**
     * Names of stations, one entry per row
     */
    List<String> stationNames;

    /**
     * Number of decimal digits to display
     */
    int nDigits = 5;

    /**
     * Creates a bounds table from the four bracket columns.
     *
     * @param Qlower list of queue-length lower bounds
     * @param Qupper list of queue-length upper bounds
     * @param Tlower list of throughput lower bounds
     * @param Tupper list of throughput upper bounds
     */
    public NetworkBoundsTable(List<Double> Qlower, List<Double> Qupper,
                              List<Double> Tlower, List<Double> Tupper) {
        super(new ArrayList<List<Double>>(Arrays.asList(Qlower, Qupper, Tlower, Tupper)));
    }

    /**
     * Returns the queue-length lower bounds.
     *
     * @return the Qlower column
     */
    public List<Double> getQlower() {
        return get(0);
    }

    /**
     * Returns the queue-length upper bounds.
     *
     * @return the Qupper column
     */
    public List<Double> getQupper() {
        return get(1);
    }

    /**
     * Returns the throughput lower bounds.
     *
     * @return the Tlower column
     */
    public List<Double> getTlower() {
        return get(2);
    }

    /**
     * Returns the throughput upper bounds.
     *
     * @return the Tupper column
     */
    public List<Double> getTupper() {
        return get(3);
    }

    /**
     * Returns the job class names, one per row.
     *
     * @return the class names
     */
    public List<String> getClassNames() {
        return classNames;
    }

    /**
     * Sets the job class names, one per row.
     *
     * @param classNames the class names
     */
    public void setClassNames(List<String> classNames) {
        this.classNames = classNames;
    }

    /**
     * Returns the station names, one per row.
     *
     * @return the station names
     */
    public List<String> getStationNames() {
        return stationNames;
    }

    /**
     * Sets the station names, one per row.
     *
     * @param stationNames the station names
     */
    public void setStationNames(List<String> stationNames) {
        this.stationNames = stationNames;
    }

    /**
     * Sets the number of decimal digits used when printing.
     *
     * @param nd number of decimal digits
     */
    public void setNumberOfDigits(int nd) {
        this.nDigits = nd;
    }

    /**
     * Sets the solver options associated with this table.
     *
     * @param options the solver options
     */
    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    /**
     * Prints the table using the options attached to it.
     */
    public void print() {
        this.print(this.options);
    }

    /**
     * Prints the table.
     *
     * <p>Unlike {@link NetworkAvgTable#print}, no row is suppressed on the
     * grounds of being all-zero: a bound of zero is a meaningful result, and a
     * one-sided family would otherwise lose rows whose only present side is
     * zero. Row filtering happens once, in SolverBA.getBoundsTable.
     *
     * @param options the solver options controlling verbosity
     */
    public void print(SolverOptions options) {
        if (options == null || options.verbose != VerboseLevel.SILENT) {
            String[] headers = {"Station", "JobClass", "Qlower", "Qupper", "Tlower", "Tupper"};
            List<String[]> rows = new ArrayList<String[]>();
            for (int i = 0; i < stationNames.size(); i++) {
                rows.add(new String[]{
                        stationNames.get(i),
                        classNames.get(i),
                        fmtValue(getQlower().get(i), nDigits),
                        fmtValue(getQupper().get(i), nDigits),
                        fmtValue(getTlower().get(i), nDigits),
                        fmtValue(getTupper().get(i), nDigits)});
            }
            printFormattedTable(headers, rows);
        }
    }

    /**
     * Prints the table using the options attached to it.
     */
    public void printTable() {
        this.print();
    }

    /**
     * Prints the table.
     *
     * @param options the solver options controlling verbosity
     */
    public void printTable(SolverOptions options) {
        this.print(options);
    }
}
