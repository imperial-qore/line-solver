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
 * Table of exact higher moments of the per-chain queue length, one row per
 * (Station, Chain).
 *
 * <p>This is the chain-level analogue of {@link NetworkAvgChainTable}, and it sits
 * between the two other moment tables: {@link NetworkMomentTable} is per class,
 * {@link NetworkMomentStationTable} is per station total, and this one is per chain,
 * i.e. per group of classes that circulate together. Row (i,c) reports the moments
 * of Q_(i,c) = sum_(r in chain c) n(i,r).</p>
 *
 * <p>All three tables are the same recursion under different groupings of the
 * classes: the generating parameter scales the service times of a class subset T at
 * a station, and the moments it produces are those of sum_(r in T) n(i,r). T = {r}
 * gives the per-class table, T = chain gives this one, T = all classes gives the
 * station table. That is Theorem 1 of Akyildiz and Strelen; Strelen's own x_i is the
 * last case. Unlike the per-class table, order 3 is fully available here.</p>
 *
 * <p>Which columns are present depends on the moment orders requested: order 1
 * contributes QLen, order 2 contributes QLenVar and QLenSCV, order 3 contributes
 * QLenM3 and QLenSkew; see {@link #getVariableNames()}.</p>
 *
 * @see NetworkSolver#getMomentChainTable(int[])
 * @see NetworkMomentStationResult
 * @see NetworkMomentTable
 */
public class NetworkMomentChainTable extends AvgTable {

    /**
     * Names of stations, one per row
     */
    List<String> stationNames;

    /**
     * Names of chains, one per row
     */
    List<String> chainNames;

    /**
     * Names of the numeric columns, in table order
     */
    List<String> columnNames;

    /**
     * Raw moment results behind this table
     */
    NetworkMomentStationResult moments;

    /**
     * Number of decimal digits to display
     */
    int nDigits = 5;

    /**
     * Creates a chain moment table from the numeric columns that the requested
     * moment orders select.
     *
     * @param columnNames names of the numeric columns, in table order
     * @param columns     the numeric columns, aligned with {@code columnNames}
     */
    public NetworkMomentChainTable(List<String> columnNames, List<List<Double>> columns) {
        super(new ArrayList<List<Double>>(columns));
        this.columnNames = new ArrayList<String>(columnNames);
    }

    /**
     * The names of every column of the table, the two name columns first: Station,
     * Chain, then the numeric columns that the requested moment orders selected.
     * Mirrors {@code T.Properties.VariableNames} of the MATLAB table.
     *
     * @return the column names, in table order
     */
    public List<String> getVariableNames() {
        List<String> out = new ArrayList<String>();
        out.add("Station");
        out.add("Chain");
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
            throw new RuntimeException("Unrecognized moment chain table column '" + name
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
     * @return the mean chain queue length column (order 1)
     */
    public List<Double> getQLen() {
        return getColumn("QLen");
    }

    /**
     * @return the chain queue-length variance column (order 2)
     */
    public List<Double> getQLenVar() {
        return getColumn("QLenVar");
    }

    /**
     * @return the chain queue-length SCV column (order 2)
     */
    public List<Double> getQLenSCV() {
        return getColumn("QLenSCV");
    }

    /**
     * @return the chain queue-length third moment column (order 3)
     */
    public List<Double> getQLenM3() {
        return getColumn("QLenM3");
    }

    /**
     * @return the chain queue-length skewness column (order 3)
     */
    public List<Double> getQLenSkew() {
        return getColumn("QLenSkew");
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
     * @return the chain names, one per row
     */
    public List<String> getChainNames() {
        return chainNames;
    }

    /**
     * Sets the chain names, one per row.
     *
     * @param chainNames the chain names
     */
    public void setChainNames(List<String> chainNames) {
        this.chainNames = chainNames;
    }

    /**
     * Returns the raw moment results, including the cross-chain and cross-station
     * covariances that this table does not show.
     *
     * @return the raw moment results
     */
    public NetworkMomentStationResult getMoments() {
        return moments;
    }

    /**
     * Sets the raw moment results.
     *
     * @param moments the raw moment results
     */
    public void setMoments(NetworkMomentStationResult moments) {
        this.moments = moments;
    }

    /**
     * Returns the index of the row for a (station, chain) pair, or -1 if absent.
     *
     * @param stationName the station name
     * @param chainName   the chain name
     * @return the row index, or -1
     */
    public int findRow(String stationName, String chainName) {
        for (int i = 0; i < stationNames.size(); i++) {
            if (stationNames.get(i).equals(stationName) && chainNames.get(i).equals(chainName)) {
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
        int maxChainLength = "Chain".length();
        for (String name : stationNames) {
            if (name.length() > maxStationLength) {
                maxStationLength = name.length();
            }
        }
        for (String name : chainNames) {
            if (name.length() > maxChainLength) {
                maxChainLength = name.length();
            }
        }
        StringBuilder fmt = new StringBuilder();
        fmt.append(String.format("%%-%ds%%-%ds", maxStationLength + 2, maxChainLength + 2));
        for (int j = 0; j < columnNames.size(); j++) {
            fmt.append("%-14s");
        }
        String format = fmt.toString();
        Object[] header = new Object[2 + columnNames.size()];
        header[0] = "Station";
        header[1] = "Chain";
        for (int j = 0; j < columnNames.size(); j++) {
            header[2 + j] = columnNames.get(j);
        }
        System.out.printf(format, header);
        System.out.println();
        DecimalFormat nf = new DecimalFormat("#0.#####");
        nf.setMinimumFractionDigits(nDigits);
        int tableWidth = (maxStationLength + 2) + (maxChainLength + 2) + columnNames.size() * 14;
        printRule(tableWidth);
        for (int i = 0; i < stationNames.size(); i++) {
            Object[] row = new Object[2 + columnNames.size()];
            row[0] = stationNames.get(i);
            row[1] = chainNames.get(i);
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
