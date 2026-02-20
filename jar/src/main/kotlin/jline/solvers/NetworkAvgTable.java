/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.lang.JobClass;
import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.nodes.Station;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Table for displaying network performance metrics organized by station and job class.
 * <p>
 * This class provides a structured way to access and display steady-state performance
 * metrics including queue lengths, utilizations, response times, throughputs, and
 * arrival rates. The table is organized with rows representing station-class pairs
 * and columns representing different performance metrics.
 *
 * @see AvgTable
 * @see SolverResult
 */
public class NetworkAvgTable extends AvgTable {
    /**
     * Names of job classes in the network
     */
    List<String> classNames;

    /**
     * Names of stations in the network
     */
    List<String> stationNames;

    /**
     * Number of decimal digits to display
     */
    int nDigits = 5;

    /**
     * Creates a new NetworkAvgTable with the specified performance metric values.
     *
     * @param Qval     list of queue length values
     * @param Uval     list of utilization values
     * @param Rval     list of response time values
     * @param Residval list of residence time values
     * @param ArvRval  list of arrival rate values
     * @param Tval     list of throughput values
     */
    public NetworkAvgTable(List<Double> Qval, List<Double> Uval, List<Double> Rval, List<Double> Residval, List<Double> ArvRval, List<Double> Tval) {
        super(new ArrayList<>(Arrays.asList(Qval, Uval, Rval, Residval, ArvRval, Tval)));
    }

    /**
     * Creates a new NetworkAvgTable with all performance metric values including tardiness.
     *
     * @param Qval     list of queue length values
     * @param Uval     list of utilization values
     * @param Rval     list of response time values
     * @param Residval list of residence time values
     * @param ArvRval  list of arrival rate values
     * @param Tval     list of throughput values
     * @param Tardval  list of tardiness values
     * @param SysTardval list of system tardiness values
     */
    public NetworkAvgTable(List<Double> Qval, List<Double> Uval, List<Double> Rval, List<Double> Residval, List<Double> ArvRval, List<Double> Tval, List<Double> Tardval, List<Double> SysTardval) {
        super(new ArrayList<>(Arrays.asList(Qval, Uval, Rval, Residval, ArvRval, Tval, Tardval, SysTardval)));
    }

    /**
     * Extracts rows from the table matching the specified name (station or class).
     *
     * @param T       the table to extract from
     * @param anyname name of station or class to match
     * @return filtered table containing matching rows
     */
    public static AvgTable tget(NetworkAvgTable T, String anyname) {
        return T.tget(anyname);
    }

    /**
     * Extracts rows from the table matching both station and class names.
     *
     * @param T           the table to extract from
     * @param stationname name of the station to match
     * @param classname   name of the job class to match
     * @return filtered table containing matching rows
     */
    public static AvgTable tget(NetworkAvgTable T, String stationname, String classname) {
        return T.tget(stationname, classname);
    }

    private String formatValue(double value, DecimalFormat nf) {
        if (value == 0.0) {
            return "0";
        } else if (Math.abs(value) < 1e-5) {
            return String.format("%.1e", value);
        } else {
            return nf.format(value);
        }
    }

    /**
     * Returns the values from the specified column.
     *
     * @param col column index (0=QLen, 1=Util, 2=RespT, 3=ResidT, 4=ArvR, 5=Tput)
     * @return list of values from the column
     */
    public List<Double> get(int col) {
        return this.T.getColumn(col).toList1D();
    }

    public List<Double> getArvR() {
        return this.get(4);
    }

    /**
     * Returns the list of job class names.
     *
     * @return list of class names
     */
    public java.util.List<String> getClassNames() {
        return classNames;
    }

    /**
     * Sets the list of job class names.
     *
     * @param classNames list of class names to set
     */
    public void setClassNames(java.util.List<String> classNames) {
        this.classNames = classNames;
    }

    /**
     * Returns the queue length values from the table.
     *
     * @return list of queue length values
     */
    public List<Double> getQLen() {
        return this.get(0);
    }

    /**
     * Returns the residence time values from the table.
     *
     * @return list of residence time values
     */
    public List<Double> getResidT() {
        return this.get(3);
    }

    /**
     * Returns the response time values from the table.
     *
     * @return list of response time values
     */
    public List<Double> getRespT() {
        return this.get(2);
    }

    /**
     * Returns the list of station names.
     *
     * @return list of station names
     */
    public java.util.List<String> getStationNames() {
        return stationNames;
    }

    /**
     * Sets the list of station names.
     *
     * @param stationNames list of station names to set
     */
    public void setStationNames(java.util.List<String> stationNames) {
        this.stationNames = stationNames;
    }

    public List<Double> getTput() {
        return this.get(5);
    }

    /**
     * Returns the tardiness values from the table if available.
     *
     * @return list of tardiness values, or empty list if not available
     */
    public List<Double> getTard() {
        if (this.T.getNumCols() > 6) {
            return this.get(6);
        }
        return new ArrayList<>();
    }

    /**
     * Returns the system tardiness values from the table if available.
     *
     * @return list of system tardiness values, or empty list if not available
     */
    public List<Double> getSysTard() {
        if (this.T.getNumCols() > 7) {
            return this.get(7);
        }
        return new ArrayList<>();
    }

    //    List<Double> getArvR() {
//        return this.T.get(5);
//    }

    /**
     * Returns the utilization values from the table.
     *
     * @return list of utilization values
     */
    public List<Double> getUtil() {
        return this.get(1);
    }

    public void print() {
        this.print(this.options);
    }

    public void print(SolverOptions options) {
        if (options.verbose != VerboseLevel.SILENT) {
            // Calculate the maximum length of names for formatting
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

            // Print header with dynamic width for Station and JobClass
            String format = String.format("%%-%ds%%-%ds%%-12s%%-12s%%-12s%%-12s%%-12s%%-12s",
                    maxStationLength + 2, maxJobClassLength + 2);
            System.out.printf(format, "Station", "JobClass", "QLen", "Util", "RespT", "ResidT", "ArvR", "Tput");
            System.out.println();

            // Calculate and print separator line with dynamic length
            DecimalFormat nf = new DecimalFormat("#0.#####");
            nf.setMinimumFractionDigits(nDigits);
            int tableWidth = (maxStationLength + 2) + (maxJobClassLength + 2) + 6 * (12 - 5 + nDigits);
            for (int i = 0; i < tableWidth; i++) {
                System.out.print("-");
            }
            System.out.println();

            // Print data rows with the same dynamic width format
            for (int i = 0; i < stationNames.size(); i++) {
                if (getQLen().get(i) > GlobalConstants.Zero ||
                        getUtil().get(i) > GlobalConstants.Zero ||
                        getRespT().get(i) > GlobalConstants.Zero ||
                        getResidT().get(i) > GlobalConstants.Zero ||
                        getArvR().get(i) > GlobalConstants.Zero ||
                        getTput().get(i) > GlobalConstants.Zero) {

                    System.out.format(format + "\n",
                            stationNames.get(i),
                            classNames.get(i),
                            formatValue(getQLen().get(i), nf),
                            formatValue(getUtil().get(i), nf),
                            formatValue(getRespT().get(i), nf),
                            formatValue(getResidT().get(i), nf),
                            formatValue(getArvR().get(i), nf),
                            formatValue(getTput().get(i), nf));
                }
            }

            // Print separator line again at the end
            for (int i = 0; i < tableWidth; i++) {
                System.out.print("-");
            }
            System.out.println();
        }
    }

    public void printTable() {
        this.print();
    }

    public void printTable(SolverOptions options) {
        this.print(options);
    }

    /**
     * Sets the number of decimal digits to display in formatted output.
     *
     * @param nd number of decimal digits
     */
    public void setNumberOfDigits(int nd) {
        this.nDigits = nd;
    }

    /**
     * Sets the solver options for this table.
     *
     * @param options solver options to associate with this table
     */
    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    public NetworkAvgTable tget(Station station, JobClass jobclass) {
        return this.tget(station.getName(), jobclass.getName());
    }

//    public AvgTable tget(String metricname, Station station, JobClass jobclass) {
//        return this.tget(metricname, (String) station.getName(), (String) jobclass.getName());
//    }
//
//    public AvgTable tget(String metricname, String stationname, String classname) {
//        int rowIdx = this.stationNames.indexOf(stationname);
//        int colIdx = this.classNames.indexOf(classname);
//        if (rowIdx < 0 || colIdx < 0) {
//            NetworkAvgTable filteredAvgTable = new NetworkAvgTable(new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
//            filteredAvgTable.setOptions(this.options);
//            filteredAvgTable.setStationNames(new ArrayList<String>());
//            filteredAvgTable.setClassNames(new ArrayList<String>());
//            return filteredAvgTable;
//        }
//        int indexToKeep = rowIdx * this.classNames.size() + colIdx - 1;
//        List<Double> MetricVal;
//        AvgTable filteredAvgTable = null;
//        switch (metricname) {
//            case "QLen":
//                MetricVal = this.getQLen();
//                if (indexToKeep > 0) {
//                    MetricVal.subList(0, indexToKeep).clear();
//                }
//                if (indexToKeep < MetricVal.size()) {
//                    MetricVal.subList(1, MetricVal.size()).clear();
//                }
//                filteredAvgTable = new NetworkAvgQLenTable(MetricVal);
//                filteredAvgTable.setOptions(this.options);
//                ((NetworkAvgQLenTable) filteredAvgTable).setStationNames(Collections.singletonList(this.stationNames.get(rowIdx)));
//                ((NetworkAvgQLenTable) filteredAvgTable).setClassNames(Collections.singletonList(this.classNames.get(colIdx)));
//                break;
//            case "Util":
//                MetricVal = this.getUtil();
//                if (indexToKeep > 0) {
//                    MetricVal.subList(0, indexToKeep).clear();
//                }
//                if (indexToKeep < MetricVal.size()) {
//                    MetricVal.subList(1, MetricVal.size()).clear();
//                }
//                filteredAvgTable = new NetworkAvgUtilTable(MetricVal);
//                filteredAvgTable.setOptions(this.options);
//                ((NetworkAvgUtilTable) filteredAvgTable).setStationNames(Collections.singletonList(this.stationNames.get(rowIdx)));
//                ((NetworkAvgUtilTable) filteredAvgTable).setClassNames(Collections.singletonList(this.classNames.get(colIdx)));
//                break;
//            case "RespT":
//                MetricVal = this.getRespT();
//                if (indexToKeep > 0) {
//                    MetricVal.subList(0, indexToKeep).clear();
//                }
//                if (indexToKeep < MetricVal.size()) {
//                    MetricVal.subList(1, MetricVal.size()).clear();
//                }
//                filteredAvgTable = new NetworkAvgRespTTable(MetricVal);
//                filteredAvgTable.setOptions(this.options);
//                ((NetworkAvgRespTTable) filteredAvgTable).setStationNames(Collections.singletonList(this.stationNames.get(rowIdx)));
//                ((NetworkAvgRespTTable) filteredAvgTable).setClassNames(Collections.singletonList(this.classNames.get(colIdx)));
//                break;
//            case "ResidT":
//                MetricVal = this.getResidT();
//                if (indexToKeep > 0) {
//                    MetricVal.subList(0, indexToKeep).clear();
//                }
//                if (indexToKeep < MetricVal.size()) {
//                    MetricVal.subList(1, MetricVal.size()).clear();
//                }
//                filteredAvgTable = new NetworkAvgResidTTable(MetricVal);
//                filteredAvgTable.setOptions(this.options);
//                ((NetworkAvgResidTTable) filteredAvgTable).setStationNames(Collections.singletonList(this.stationNames.get(rowIdx)));
//                ((NetworkAvgResidTTable) filteredAvgTable).setClassNames(Collections.singletonList(this.classNames.get(colIdx)));
//                break;
//            case "ArvR":
//                MetricVal = this.getArvR();
//                if (indexToKeep > 0) {
//                    MetricVal.subList(0, indexToKeep).clear();
//                }
//                if (indexToKeep < MetricVal.size()) {
//                    MetricVal.subList(1, MetricVal.size()).clear();
//                }
//                filteredAvgTable = new NetworkAvgArvRTable(MetricVal);
//                filteredAvgTable.setOptions(this.options);
//                ((NetworkAvgArvRTable) filteredAvgTable).setStationNames(Collections.singletonList(this.stationNames.get(rowIdx)));
//                ((NetworkAvgArvRTable) filteredAvgTable).setClassNames(Collections.singletonList(this.classNames.get(colIdx)));
//                break;
//            case "Tput":
//                MetricVal = this.getTput();
//                if (indexToKeep > 0) {
//                    MetricVal.subList(0, indexToKeep).clear();
//                }
//                if (indexToKeep < MetricVal.size()) {
//                    MetricVal.subList(1, MetricVal.size()).clear();
//                }
//                filteredAvgTable = new NetworkAvgTputTable(MetricVal);
//                filteredAvgTable.setOptions(this.options);
//                ((NetworkAvgTputTable) filteredAvgTable).setStationNames(Collections.singletonList(this.stationNames.get(rowIdx)));
//                ((NetworkAvgTputTable) filteredAvgTable).setClassNames(Collections.singletonList(this.classNames.get(colIdx)));
//                break;
//        }
//        return filteredAvgTable;
//    }

    public NetworkAvgTable tget(String anyname) {
        int rowIdx = this.stationNames.indexOf(anyname);
        int colIdx = this.classNames.indexOf(anyname);
        if (rowIdx < 0 && colIdx < 0) {
            NetworkAvgTable filteredAvgTable = new NetworkAvgTable(new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
            filteredAvgTable.setOptions(this.options);
            filteredAvgTable.setStationNames(new ArrayList<String>());
            filteredAvgTable.setClassNames(new ArrayList<String>());
            return filteredAvgTable;
        }

        List<Integer> indexToKeep = new ArrayList<>();
        List<String> stationNamesFilt = new ArrayList<>();
        List<String> classNamesFilt = new ArrayList<>();
        for (int i = 0; i < this.stationNames.size(); i++) {
            if (this.stationNames.get(i).compareTo(anyname) == 0 || this.classNames.get(i).compareTo(anyname) == 0) {
                stationNamesFilt.add(this.stationNames.get(i));
                classNamesFilt.add(this.classNames.get(i));
                indexToKeep.add(i);
            }
        }

        List<Double> Qval = this.getQLen();
        List<Double> QvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Qval.size()) {
                QvalFilt.add(Qval.get(index));
            }
        }
        List<Double> Uval = this.getUtil();
        List<Double> UvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Uval.size()) {
                UvalFilt.add(Uval.get(index));
            }
        }
        List<Double> Rval = this.getRespT();
        List<Double> RvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Rval.size()) {
                RvalFilt.add(Rval.get(index));
            }
        }
        List<Double> Wval = this.getResidT();
        List<Double> WvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Wval.size()) {
                WvalFilt.add(Wval.get(index));
            }
        }
        List<Double> Aval = this.getArvR();
        List<Double> AvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Aval.size()) {
                AvalFilt.add(Aval.get(index));
            }
        }
        List<Double> Tval = this.getTput();
        List<Double> TvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Tval.size()) {
                TvalFilt.add(Tval.get(index));
            }
        }
        NetworkAvgTable filteredAvgTable = new NetworkAvgTable(QvalFilt, UvalFilt, RvalFilt, WvalFilt, AvalFilt, TvalFilt);
        filteredAvgTable.setOptions(this.options);
        filteredAvgTable.setStationNames(stationNamesFilt);
        filteredAvgTable.setClassNames(classNamesFilt);
        return filteredAvgTable;
    }

    public NetworkAvgTable tget(String stationname, String classname) {
        int rowIdx = this.stationNames.indexOf(stationname);
        int colIdx = this.classNames.indexOf(classname);
        if (rowIdx < 0 || colIdx < 0) {
            NetworkAvgTable filteredAvgTable = new NetworkAvgTable(new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
            filteredAvgTable.setOptions(this.options);
            filteredAvgTable.setStationNames(new ArrayList<String>());
            filteredAvgTable.setClassNames(new ArrayList<String>());
            return filteredAvgTable;
        }
        int indexToKeep = 0;
        for (int i = 0; i < this.stationNames.size(); i++) {
            if (this.stationNames.get(i).compareTo(stationname) == 0 && this.classNames.get(i).compareTo(classname) == 0) {
                indexToKeep = i;
                break;
            }
        }
        List<Double> Qval = this.getQLen();
        if (indexToKeep > 0) {
            Qval.subList(0, indexToKeep).clear();
        }
        if (indexToKeep < Qval.size()) {
            Qval.subList(1, Qval.size()).clear();
        }
        List<Double> Uval = this.getUtil();
        if (indexToKeep > 0) {
            Uval.subList(0, indexToKeep).clear();
        }
        if (indexToKeep < Uval.size()) {
            Uval.subList(1, Uval.size()).clear();
        }
        List<Double> Rval = this.getRespT();
        if (indexToKeep > 0) {
            Rval.subList(0, indexToKeep).clear();
        }
        if (indexToKeep < Rval.size()) {
            Rval.subList(1, Rval.size()).clear();
        }
        List<Double> Residval = this.getResidT();
        if (indexToKeep > 0) {
            Residval.subList(0, indexToKeep).clear();
        }
        if (indexToKeep < Residval.size()) {
            Residval.subList(1, Residval.size()).clear();
        }
        List<Double> ArvRval = this.getArvR();
        if (indexToKeep > 0) {
            ArvRval.subList(0, indexToKeep).clear();
        }
        if (indexToKeep < ArvRval.size()) {
            ArvRval.subList(1, ArvRval.size()).clear();
        }
        List<Double> Tval = this.getTput();
        if (indexToKeep > 0) {
            Tval.subList(0, indexToKeep).clear();
        }
        if (indexToKeep < Tval.size()) {
            Tval.subList(1, Tval.size()).clear();
        }
        NetworkAvgTable filteredAvgTable = new NetworkAvgTable(Qval, Uval, Rval, Residval, ArvRval, Tval);
        filteredAvgTable.setOptions(this.options);
        filteredAvgTable.setStationNames(Collections.singletonList(this.stationNames.get(rowIdx)));
        filteredAvgTable.setClassNames(Collections.singletonList(this.classNames.get(colIdx)));
        return filteredAvgTable;
    }

//    public static AvgTable tget(NetworkAvgTable T, String metricname, String stationname, String classname) {
//        return T.tget(metricname, stationname, classname);
//    }
//
//    public static AvgTable tget(NetworkAvgTable T, String metricname, Station station, JobClass jobclass) {
//        return T.tget(metricname, station, jobclass);
//    }
}
