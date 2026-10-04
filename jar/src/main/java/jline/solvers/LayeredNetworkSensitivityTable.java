/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.VerboseLevel;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Layer-wise table of performance sensitivities of a layered network, one row per
 * (Layer, Station, JobClass).
 *
 * <p>This is the layered sibling of {@link NetworkSensitivityTable}: it is the
 * concatenation of the per-layer tables that {@code SolverLN} obtains from its layer
 * solvers, with a leading Layer column carrying the name of the layer Network. The
 * columns always appear in the order Layer, Station, JobClass, dTput_dRate,
 * dRespT_dRate, dQLen_dRate, dUtil_dRate; see {@link #getVariableNames()}.</p>
 *
 * <p>IMPORTANT, on what these derivatives mean. Each entry is a derivative WITHIN ITS
 * LAYER, taken with the layer parameters that the fixed point produced held fixed. It
 * is a partial derivative of the layer submodel, not the total derivative of the
 * layered model: perturbing a host demand in one layer moves the think times,
 * populations and service rates of the other layers through the fixed-point map, and
 * that indirect term is not included here. The layer table is the right object for
 * attributing a bottleneck inside a layer, and the wrong one for predicting the
 * effect of a parameter change on the solved layered model. For the latter,
 * finite-difference the LayeredNetwork itself.</p>
 *
 * <p>{@link #getMethod()} reports the branch summary, "exact" or "fd" when every
 * layer took the same branch and "mixed" when they disagreed;
 * {@link #getLayerMethods()} reports the branch of each layer separately. These
 * mirror MATLAB's {@code SensTable.Properties.UserData.method} and
 * {@code .layerMethods}.</p>
 *
 * @see NetworkSensitivityTable
 * @see jline.solvers.ln.SolverLN#getSensitivityTable()
 */
public class LayeredNetworkSensitivityTable extends AvgTable {

    /**
     * Names of the numeric columns, in table order. Mirrors the MATLAB table, which
     * always carries all four.
     */
    private static final List<String> COLUMN_NAMES = new ArrayList<String>(Arrays.asList(
            "dTput_dRate", "dRespT_dRate", "dQLen_dRate", "dUtil_dRate"));

    /**
     * Names of the layers, one per row
     */
    List<String> layerNames;

    /**
     * Names of stations, one per row
     */
    List<String> stationNames;

    /**
     * Names of job classes, one per row
     */
    List<String> classNames;

    /**
     * Branch taken by each layer, indexed by layer; an entry is null where the layer
     * had no solver and contributed no rows
     */
    List<String> layerMethods;

    /**
     * Summary of the layer branches: "exact", "fd", or "mixed"
     */
    String method;

    /**
     * Raw differentiated-MVA result of each layer, indexed by layer; an entry is null
     * where the layer took the finite-difference branch or the open exact branch
     */
    List<jline.io.Ret.pfqnSens> layerSens;

    /**
     * Number of decimal digits to display
     */
    int nDigits = 5;

    /**
     * Creates a layered sensitivity table from its four numeric columns.
     *
     * @param dTput  the dTput_dRate column
     * @param dRespT the dRespT_dRate column
     * @param dQLen  the dQLen_dRate column
     * @param dUtil  the dUtil_dRate column
     */
    public LayeredNetworkSensitivityTable(List<Double> dTput, List<Double> dRespT,
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
     * The names of every column of the table, the three name columns first: Layer,
     * Station, JobClass, then the four numeric columns. Mirrors
     * {@code T.Properties.VariableNames} of the MATLAB table.
     *
     * @return the column names, in table order
     */
    public List<String> getVariableNames() {
        List<String> out = new ArrayList<String>();
        out.add("Layer");
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
     * @return the dTput_dRate column
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
     * @return the layer names, one per row
     */
    public List<String> getLayerNames() {
        return layerNames;
    }

    /**
     * Sets the layer names, one per row.
     *
     * @param layerNames the layer names
     */
    public void setLayerNames(List<String> layerNames) {
        this.layerNames = layerNames;
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
     * Returns the branch summary of the table, the counterpart of MATLAB's
     * {@code SensTable.Properties.UserData.method}.
     *
     * @return "exact" or "fd" when every layer took the same branch, "mixed" when the
     *         layers disagreed, and the empty string when no layer contributed
     */
    public String getMethod() {
        return method;
    }

    /**
     * Sets the branch summary of the table.
     *
     * @param method "exact", "fd" or "mixed"
     */
    public void setMethod(String method) {
        this.method = method;
    }

    /**
     * Returns the branch taken by each layer, the counterpart of MATLAB's
     * {@code SensTable.Properties.UserData.layerMethods}.
     *
     * @return one entry per layer, null where the layer contributed no rows
     */
    public List<String> getLayerMethods() {
        return layerMethods;
    }

    /**
     * Sets the branch taken by each layer.
     *
     * @param layerMethods one entry per layer
     */
    public void setLayerMethods(List<String> layerMethods) {
        this.layerMethods = layerMethods;
    }

    /**
     * Returns the raw differentiated-MVA result of each layer, the counterpart of the
     * cell array that MATLAB returns as the second output alongside the table.
     *
     * @return one entry per layer, null where the layer took the finite-difference
     *         branch or the open exact branch
     */
    public List<jline.io.Ret.pfqnSens> getLayerSens() {
        return layerSens;
    }

    /**
     * Sets the raw differentiated-MVA result of each layer.
     *
     * @param layerSens one entry per layer
     */
    public void setLayerSens(List<jline.io.Ret.pfqnSens> layerSens) {
        this.layerSens = layerSens;
    }

    /**
     * Returns the index of the row for a (layer, station, class) triple, or -1 if
     * absent.
     *
     * @param layerName   the layer name
     * @param stationName the station name
     * @param className   the job class name
     * @return the row index, or -1
     */
    public int findRow(String layerName, String stationName, String className) {
        for (int i = 0; i < layerNames.size(); i++) {
            if (layerNames.get(i).equals(layerName) && stationNames.get(i).equals(stationName)
                    && classNames.get(i).equals(className)) {
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
        int maxLayerLength = "Layer".length();
        int maxStationLength = "Station".length();
        int maxJobClassLength = "JobClass".length();
        for (String name : layerNames) {
            if (name.length() > maxLayerLength) {
                maxLayerLength = name.length();
            }
        }
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
        fmt.append(String.format("%%-%ds%%-%ds%%-%ds", maxLayerLength + 2,
                maxStationLength + 2, maxJobClassLength + 2));
        for (int j = 0; j < COLUMN_NAMES.size(); j++) {
            fmt.append("%-16s");
        }
        String format = fmt.toString();
        Object[] header = new Object[3 + COLUMN_NAMES.size()];
        header[0] = "Layer";
        header[1] = "Station";
        header[2] = "JobClass";
        for (int j = 0; j < COLUMN_NAMES.size(); j++) {
            header[3 + j] = COLUMN_NAMES.get(j);
        }
        System.out.printf(format, header);
        System.out.println();
        DecimalFormat nf = new DecimalFormat("#0.#####");
        nf.setMinimumFractionDigits(nDigits);
        int tableWidth = (maxLayerLength + 2) + (maxStationLength + 2) + (maxJobClassLength + 2)
                + COLUMN_NAMES.size() * 16;
        printRule(tableWidth);
        for (int i = 0; i < layerNames.size(); i++) {
            Object[] row = new Object[3 + COLUMN_NAMES.size()];
            row[0] = layerNames.get(i);
            row[1] = stationNames.get(i);
            row[2] = classNames.get(i);
            for (int j = 0; j < COLUMN_NAMES.size(); j++) {
                row[3 + j] = formatValue(this.T.get(i, j), nf);
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
