/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.VerboseLevel;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Table of exact higher moments of the per-class performance measures, one row per
 * (Station, JobClass).
 *
 * <p>Which columns are present depends on the moment orders requested of
 * {@link NetworkSolver#getMomentTable(int[])}: order 1 contributes QLen and RespT,
 * order 2 contributes QLenVar, QLenSCV, RespTVar and RespTSCV, order 3 contributes
 * QLenSkew and RespTSkew. They always appear in the order Station, JobClass, QLen,
 * QLenVar, QLenSCV, QLenSkew, RespT, RespTVar, RespTSCV, RespTSkew, filtered to the
 * requested set; see {@link #getVariableNames()}.</p>
 *
 * <p>QLenSkew is a genuinely per-class quantity: the generating parameter need not
 * scale a whole demand column, and scaling L(i,r) alone is Akyildiz-Strelen Theorem
 * 1 with the class subset T = {r}, which generates the moments of n(i,r) itself. It
 * is available only for closed single-server models, the scope of
 * {@code Pfqn_sens_mom}, and is NaN otherwise.</p>
 *
 * <p>The means agree with {@link NetworkAvgTable}, which reaches them by a different
 * code path. RespTVar, RespTSCV and RespTSkew are NaN at any station that is not
 * FCFS: the sojourn-time distribution of a processor-sharing or LCFS center is not
 * known in general, so there is no correct value to report.</p>
 *
 * @see NetworkSolver#getMomentTable()
 * @see NetworkMomentResult
 * @see NetworkMomentStationTable
 */
public class NetworkMomentTable extends AvgTable {

    /**
     * Names of job classes, one per row
     */
    List<String> classNames;

    /**
     * Names of stations, one per row
     */
    List<String> stationNames;

    /**
     * Names of the numeric columns, in table order
     */
    List<String> columnNames;

    /**
     * Raw moment results behind this table
     */
    NetworkMomentResult moments;

    /**
     * Number of decimal digits to display
     */
    int nDigits = 5;

    /**
     * Creates a moment table from the numeric columns that the requested moment
     * orders select.
     *
     * @param columnNames names of the numeric columns, in table order
     * @param columns     the numeric columns, aligned with {@code columnNames}
     */
    public NetworkMomentTable(List<String> columnNames, List<List<Double>> columns) {
        super(new ArrayList<List<Double>>(columns));
        this.columnNames = new ArrayList<String>(columnNames);
    }

    /**
     * The names of every column of the table, the two name columns first: Station,
     * JobClass, then the numeric columns that the requested moment orders selected.
     * Mirrors {@code T.Properties.VariableNames} of the MATLAB table.
     *
     * @return the column names, in table order
     */
    public List<String> getVariableNames() {
        List<String> out = new ArrayList<String>();
        out.add("Station");
        out.add("JobClass");
        out.addAll(columnNames);
        return out;
    }

    /**
     * Whether a numeric column is present, i.e. whether its moment order was
     * requested.
     *
     * @param name the column name
     * @return true if the column is present
     */
    public boolean hasColumn(String name) {
        return columnNames.contains(name);
    }

    /**
     * Returns a numeric column by name.
     *
     * @param name the column name
     * @return the column values
     * @throws RuntimeException if the column was not selected by the requested
     *                          moment orders
     */
    public List<Double> getColumn(String name) {
        int j = columnNames.indexOf(name);
        if (j < 0) {
            throw new RuntimeException("Unrecognized moment table column '" + name
                    + "'; the columns present are " + getVariableNames()
                    + ". Request a higher moment order to obtain it.");
        }
        return this.T.getColumn(j).toList1D();
    }

    /**
     * Returns the values of the given numeric column.
     *
     * @param col column index into {@link #getVariableNames()} minus the two name
     *            columns
     * @return the column values
     */
    public List<Double> get(int col) {
        return this.T.getColumn(col).toList1D();
    }

    /**
     * @return the mean queue length column (order 1)
     */
    public List<Double> getQLen() {
        return getColumn("QLen");
    }

    /**
     * @return the queue-length variance column (order 2)
     */
    public List<Double> getQLenVar() {
        return getColumn("QLenVar");
    }

    /**
     * @return the queue-length SCV column (order 2)
     */
    public List<Double> getQLenSCV() {
        return getColumn("QLenSCV");
    }

    /**
     * @return the mean response time column (order 1)
     */
    public List<Double> getRespT() {
        return getColumn("RespT");
    }

    /**
     * @return the response-time variance column (order 2), NaN away from FCFS
     */
    public List<Double> getRespTVar() {
        return getColumn("RespTVar");
    }

    /**
     * @return the response-time SCV column (order 2), NaN away from FCFS
     */
    public List<Double> getRespTSCV() {
        return getColumn("RespTSCV");
    }

    /**
     * @return the response-time skewness column (order 3), NaN away from FCFS
     */
    public List<Double> getRespTSkew() {
        return getColumn("RespTSkew");
    }

    /**
     * @return the per-class queue-length skewness column (order 3), NaN unless the
     *         model is closed single-server
     */
    public List<Double> getQLenSkew() {
        return getColumn("QLenSkew");
    }

    /**
     * @return the job class names, one per row
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
     * @return the station names, one per row
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
     * Returns the raw moment results, including the full per-station covariance
     * blocks that this table shows only the diagonal of.
     *
     * @return the raw moment results
     */
    public NetworkMomentResult getMoments() {
        return moments;
    }

    /**
     * Sets the raw moment results.
     *
     * @param moments the raw moment results
     */
    public void setMoments(NetworkMomentResult moments) {
        this.moments = moments;
    }

    /**
     * Returns the index of the row for a (station, class) pair, or -1 if absent.
     * A pair is absent when the class does not visit the station.
     *
     * @param stationName the station name
     * @param className   the job class name
     * @return the row index, or -1
     */
    public int findRow(String stationName, String className) {
        for (int i = 0; i < stationNames.size(); i++) {
            if (stationNames.get(i).equals(stationName) && classNames.get(i).equals(className)) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Sets the number of decimal digits to display in formatted output.
     *
     * @param nDigits the number of digits
     */
    public void setNDigits(int nDigits) {
        this.nDigits = nDigits;
    }

    public void print() {
        this.print(this.options);
    }

    /**
     * Prints the table.
     *
     * @param options the solver options controlling verbosity
     */
    public void print(SolverOptions options) {
        if (options != null && options.verbose == VerboseLevel.SILENT) {
            return;
        }
        int maxStationLength = "Station".length();
        int maxJobClassLength = "JobClass".length();
        for (String name : stationNames) {
            if (name.length() > maxStationLength) {
                maxStationLength = name.length();
            }
        }
        for (String name : classNames) {
            if (name.length() > maxJobClassLength) {
                maxJobClassLength = name.length();
            }
        }
        StringBuilder fmt = new StringBuilder();
        fmt.append(String.format("%%-%ds%%-%ds", maxStationLength + 2, maxJobClassLength + 2));
        for (int j = 0; j < columnNames.size(); j++) {
            fmt.append("%-14s");
        }
        String format = fmt.toString();
        Object[] header = new Object[2 + columnNames.size()];
        header[0] = "Station";
        header[1] = "JobClass";
        for (int j = 0; j < columnNames.size(); j++) {
            header[2 + j] = columnNames.get(j);
        }
        System.out.printf(format, header);
        System.out.println();
        DecimalFormat nf = new DecimalFormat("#0.#####");
        nf.setMinimumFractionDigits(nDigits);
        int tableWidth = (maxStationLength + 2) + (maxJobClassLength + 2) + columnNames.size() * 14;
        printRule(tableWidth);
        for (int i = 0; i < stationNames.size(); i++) {
            Object[] row = new Object[2 + columnNames.size()];
            row[0] = stationNames.get(i);
            row[1] = classNames.get(i);
            for (int j = 0; j < columnNames.size(); j++) {
                row[2 + j] = formatValue(this.T.get(i, j), nf);
            }
            System.out.format(format + "\n", row);
        }
        printRule(tableWidth);
    }

    /**
     * Prints the table.
     */
    public void printTable() {
        this.print();
    }

    private static void printRule(int width) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < width; i++) {
            sb.append('-');
        }
        System.out.println(sb.toString());
    }

    private static String formatValue(double value, DecimalFormat nf) {
        if (Double.isNaN(value)) {
            return "NaN";
        } else if (value == 0.0) {
            return "0";
        } else if (Math.abs(value) < 1e-5) {
            return String.format("%.1e", value);
        } else {
            return nf.format(value);
        }
    }
}
