package jline.solvers;

import jline.lang.JobClass;
import jline.lang.constant.GlobalConstants;
import jline.lang.constant.VerboseLevel;
import jline.lang.nodes.Station;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class NetworkAvgResidTTable extends AvgTable {
    List<String> classNames;
    List<String> stationNames;

    public NetworkAvgResidTTable(List<Double> Metricval) {
        super(new ArrayList<>(Arrays.asList(Metricval)));
    }

    public java.util.List<String> getClassNames() {
        return classNames;
    }

    public java.util.List<String> getStationNames() {
        return stationNames;
    }

    public void setClassNames(java.util.List<String> classNames) {
        this.classNames = classNames;
    }

    public void setStationNames(java.util.List<String> stationNames) {
        this.stationNames = stationNames;
    }

    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    public List<Double> get(int col) {
        return this.T.get(col);
    }

    public List<Double> getResidT() {
        return this.T.get(0);
    }

    public void printTable() {
        this.print();
    }

    public void printTable(SolverOptions options) {
        this.print(options);
    }

    public void print() {
        this.print(this.options);
    }

    public void print(SolverOptions options) {
        if (options.verbose != VerboseLevel.SILENT) {
            System.out.printf(
                    "%-10s\t%-10s\t%-10s",
                    "Station", "JobClass", "ResidT");
            System.out.println(
                    "\n-------------------------------");
            DecimalFormat nf = new DecimalFormat("#0.#####");
            nf.setMinimumFractionDigits(5);
            for (int i = 0; i < stationNames.size(); i++) {
                if ( getResidT().get(i) > GlobalConstants.Zero) {
                    System.out.format(
                            "%-10s\t%-10s\t%-10s\n",
                            stationNames.get(i),
                            classNames.get(i),
                            nf.format(getResidT().get(i)));
                }
            }
            System.out.println(
                    "-------------------------------");
        }
    }

    public NetworkAvgResidTTable tget(Station station, JobClass jobclass) {
        return this.tget((String)station.getName(), (String)jobclass.getName());
    }

    public NetworkAvgResidTTable tget(String stationname, String classname) {
        int rowIdx = this.stationNames.indexOf(stationname);
        int colIdx = this.classNames.indexOf(classname);
        if (rowIdx < 0 || colIdx <0) {
            NetworkAvgResidTTable filteredAvgTable = new NetworkAvgResidTTable(new ArrayList<>());
            filteredAvgTable.setOptions(this.options);
            filteredAvgTable.setStationNames(new ArrayList<String>());
            filteredAvgTable.setClassNames(new ArrayList<String>());
            return filteredAvgTable;
        }
        int indexToKeep = rowIdx * this.classNames.size() + colIdx - 1;
        List<Double> myMetric = this.getResidT();
        if (indexToKeep > 0) {myMetric.subList(0, indexToKeep).clear();}
        if (indexToKeep < myMetric.size()) {myMetric.subList(1, myMetric.size()).clear();}
        NetworkAvgResidTTable filteredAvgTable = new NetworkAvgResidTTable(myMetric);
        filteredAvgTable.setOptions(this.options);
        filteredAvgTable.setStationNames(Collections.singletonList(this.stationNames.get(rowIdx)));
        filteredAvgTable.setClassNames(Collections.singletonList(this.classNames.get(colIdx)));
        return filteredAvgTable;
    }
}
