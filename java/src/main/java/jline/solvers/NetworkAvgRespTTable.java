package jline.solvers;

import jline.lang.constant.GlobalConstants;
import jline.lang.constant.VerboseLevel;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class NetworkAvgRespTTable extends AvgTable {
    List<String> classNames;
    List<String> stationNames;

    public NetworkAvgRespTTable(List<Double> Metricval) {
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

    public List<Double> getRespT() {
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
                    "Station", "JobClass", "RespT");
            System.out.println(
                    "\n-------------------------------");
            DecimalFormat nf = new DecimalFormat("#0.#####");
            nf.setMinimumFractionDigits(5);
            for (int i = 0; i < stationNames.size(); i++) {
                if ( getRespT().get(i) > GlobalConstants.Zero) {
                    System.out.format(
                            "%-10s\t%-10s\t%-10s\n",
                            stationNames.get(i),
                            classNames.get(i),
                            nf.format(getRespT().get(i)));
                }
            }
            System.out.println(
                    "-------------------------------");
        }
    }
}
