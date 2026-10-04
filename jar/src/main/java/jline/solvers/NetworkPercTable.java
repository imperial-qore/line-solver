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
 * Table of tail QUANTILES produced by the SolverBA 'snc' family, in the layout
 * of {@link NetworkAvgTable}.
 *
 * <p>Each row is a station-class pair carrying the response-time quantile
 * ({@code RespTPerc}) and the queue-length quantile ({@code QLenPerc}) at one
 * violation probability: the smallest d and b for which
 * {@code P{D > d} <= eps} and {@code P{Q > b} <= eps} are certified.
 *
 * <p>This is a different object from {@link NetworkBoundsTable}, which brackets
 * a MEAN between a lower and an upper side. A quantile has no bracket: the
 * stochastic network calculus bounds the tail from above only, so each column
 * is one-sided by construction rather than by a family being one-sided.
 *
 * @see NetworkAvgTable
 * @see jline.solvers.ba.SolverBA
 */
public class NetworkPercTable extends AvgTable {

    /** Names of job classes, one entry per row. */
    List<String> classNames;

    /** Names of stations, one entry per row. */
    List<String> stationNames;

    /** Number of decimal digits to display. */
    int nDigits = 5;

    /**
     * Creates a quantile table from the two columns.
     *
     * @param respTPerc list of response-time quantiles
     * @param qLenPerc  list of queue-length quantiles, in jobs
     */
    public NetworkPercTable(List<Double> respTPerc, List<Double> qLenPerc) {
        super(new ArrayList<List<Double>>(Arrays.asList(respTPerc, qLenPerc)));
    }

    /**
     * @return the RespTPerc column
     */
    public List<Double> getRespTPerc() {
        return get(0);
    }

    /**
     * @return the QLenPerc column
     */
    public List<Double> getQLenPerc() {
        return get(1);
    }

    /**
     * @return the job class names, one per row
     */
    public List<String> getClassNames() {
        return classNames;
    }

    /**
     * @param classNames the class names, one per row
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
     * @param stationNames the station names, one per row
     */
    public void setStationNames(List<String> stationNames) {
        this.stationNames = stationNames;
    }

    /**
     * @param nd number of decimal digits used when printing
     */
    public void setNumberOfDigits(int nd) {
        this.nDigits = nd;
    }

    /**
     * @param options the solver options associated with this table
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
     * @param options the solver options controlling verbosity
     */
    public void print(SolverOptions options) {
        if (options == null || options.verbose != VerboseLevel.SILENT) {
            String[] headers = {"Station", "JobClass", "RespTPerc", "QLenPerc"};
            List<String[]> rows = new ArrayList<String[]>();
            for (int i = 0; i < stationNames.size(); i++) {
                rows.add(new String[]{
                        stationNames.get(i),
                        classNames.get(i),
                        fmtValue(getRespTPerc().get(i), nDigits),
                        fmtValue(getQLenPerc().get(i), nDigits)});
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
