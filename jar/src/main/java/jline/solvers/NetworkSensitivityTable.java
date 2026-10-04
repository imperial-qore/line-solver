/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.VerboseLevel;
import jline.io.Ret;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Table of exact analytic performance sensitivities, one row per (Station, JobClass).
 *
 * <p>Each row gives the exact analytic derivative of that row's mean performance
 * measures with respect to that station-class service RATE: dTput_dRate,
 * dRespT_dRate, dQLen_dRate and dUtil_dRate. The columns always appear in the order
 * Station, JobClass, dTput_dRate, dRespT_dRate, dQLen_dRate, dUtil_dRate; see
 * {@link #getVariableNames()}.</p>
 *
 * <p>Every derivative is with respect to the row's OWN rate, i.e. the diagonal of the
 * full parameter Jacobian. The off-diagonal cross-station and cross-class terms are
 * not shown here; they are available from the raw {@link #getSens()} result for a
 * closed model, whose {@code dX}, {@code dQ}, {@code dU} and {@code dR} carry the
 * whole Jacobian.</p>
 *
 * <p>Closed product-form networks obtain the derivatives from {@code pfqn_sens}
 * (differentiated MVA); open product-form networks use the exact closed-form BCMP
 * sensitivities, the stations there being decoupled. Rate derivatives follow from
 * d(.)/d(rate) = -(L/rate) d(.)/dL, since L(i,r) = visits(i,r)/rate(i,r).</p>
 *
 * <p>Rows are emitted only for the (station, class) pairs that the model defines: a
 * pair whose rate is non-finite or non-positive, or whose demand D(i,r) is
 * non-positive, is skipped, since the class does not visit that station and there is
 * no rate there to differentiate.</p>
 *
 * @see NetworkSolver#getSensitivityTable()
 * @see NetworkMomentTable
 */
public class NetworkSensitivityTable extends AvgTable {

    /**
     * Names of the numeric columns, in table order. Mirrors the MATLAB table, which
     * always carries all four.
     */
    private static final List<String> COLUMN_NAMES = new ArrayList<String>(Arrays.asList(
            "dTput_dRate", "dRespT_dRate", "dQLen_dRate", "dUtil_dRate"));

    /**
     * Names of job classes, one per row
     */
    List<String> classNames;

    /**
     * Names of stations, one per row
     */
    List<String> stationNames;

    /**
     * Raw differentiated-MVA result behind this table; null on an open model
     */
    Ret.pfqnSens sens;

    /**
     * Branch that produced the table, "exact" or "fd"
     */
    String method;

    /**
     * Number of decimal digits to display
     */
    int nDigits = 5;

    /**
     * Creates a sensitivity table from its four numeric columns.
     *
     * @param dTput  the dTput_dRate column
     * @param dRespT the dRespT_dRate column
     * @param dQLen  the dQLen_dRate column
     * @param dUtil  the dUtil_dRate column
     */
    public NetworkSensitivityTable(List<Double> dTput, List<Double> dRespT,
                                   List<Double> dQLen, List<Double> dUtil) {
        super(columnsOf(dTput, dRespT, dQLen, dUtil));
    }

    private static ArrayList<List<Double>> columnsOf(List<Double> dTput, List<Double> dRespT,
                                                     List<Double> dQLen, List<Double> dUtil) {
        ArrayList<List<Double>> cols = new ArrayList<List<Double>>();
        cols.add(dTput);
        cols.add(dRespT);
        cols.add(dQLen);
        cols.add(dUtil);
        return cols;
    }

    /**
     * The names of every column of the table, the two name columns first: Station,
     * JobClass, then the four numeric columns. Mirrors
     * {@code T.Properties.VariableNames} of the MATLAB table.
     *
     * @return the column names, in table order
     */
    public List<String> getVariableNames() {
        List<String> out = new ArrayList<String>();
        out.add("Station");
        out.add("JobClass");
        out.addAll(COLUMN_NAMES);
        return out;
    }

    /**
     * Returns a numeric column by name.
     *
     * @param name the column name
     * @return the column values
     * @throws RuntimeException if the name is not one of the four numeric columns
     */
    public List<Double> getColumn(String name) {
        int j = COLUMN_NAMES.indexOf(name);
        if (j < 0) {
            throw new RuntimeException("Unrecognized sensitivity table column '" + name
                    + "'; the columns present are " + getVariableNames() + ".");
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
     * @return the dTput_dRate column, the derivative of the class throughput with
     *         respect to the row's own service rate; identically zero on an open
     *         model, whose throughput is fixed by the arrival rate
     */
    public List<Double> getDTput() {
        return getColumn("dTput_dRate");
    }

    /**
     * @return the dRespT_dRate column
     */
    public List<Double> getDRespT() {
        return getColumn("dRespT_dRate");
    }

    /**
     * @return the dQLen_dRate column
     */
    public List<Double> getDQLen() {
        return getColumn("dQLen_dRate");
    }

    /**
     * @return the dUtil_dRate column
     */
    public List<Double> getDUtil() {
        return getColumn("dUtil_dRate");
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
     * Returns the raw differentiated-MVA result behind this table, the second output
     * that MATLAB's {@code getSensitivityTable} returns alongside the table. It
     * carries the FULL parameter Jacobian, of which the table shows only the
     * diagonal, plus the base metrics the derivatives were taken at.
     *
     * @return the {@code pfqn_sens} result, or null on an open model, whose
     *         closed-form branch runs no differentiated MVA (MATLAB returns an empty
     *         {@code sens} there)
     */
    public Ret.pfqnSens getSens() {
        return sens;
    }

    /**
     * Sets the raw differentiated-MVA result.
     *
     * @param sens the {@code pfqn_sens} result
     */
    public void setSens(Ret.pfqnSens sens) {
        this.sens = sens;
    }

    /**
     * Returns the index of the row for a (station, class) pair, or -1 if absent. A
     * pair is absent when the class does not visit the station.
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
     * Returns the branch that produced the table, the counterpart of MATLAB's
     * {@code SensTable.Properties.UserData.method}.
     *
     * @return "exact" for the analytic branch, "fd" for finite differences
     */
    public String getMethod() {
        return method;
    }

    /**
     * Sets the branch that produced the table.
     *
     * @param method "exact" or "fd"
     */
    public void setMethod(String method) {
        this.method = method;
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
        for (int j = 0; j < COLUMN_NAMES.size(); j++) {
            fmt.append("%-16s");
        }
        String format = fmt.toString();
        Object[] header = new Object[2 + COLUMN_NAMES.size()];
        header[0] = "Station";
        header[1] = "JobClass";
        for (int j = 0; j < COLUMN_NAMES.size(); j++) {
            header[2 + j] = COLUMN_NAMES.get(j);
        }
        System.out.printf(format, header);
        System.out.println();
        DecimalFormat nf = new DecimalFormat("#0.#####");
        nf.setMinimumFractionDigits(nDigits);
        int tableWidth = (maxStationLength + 2) + (maxJobClassLength + 2) + COLUMN_NAMES.size() * 16;
        printRule(tableWidth);
        for (int i = 0; i < stationNames.size(); i++) {
            Object[] row = new Object[2 + COLUMN_NAMES.size()];
            row[0] = stationNames.get(i);
            row[1] = classNames.get(i);
            for (int j = 0; j < COLUMN_NAMES.size(); j++) {
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
