package jline.solvers;

import jline.lang.constant.GlobalConstants;
import jline.lang.constant.VerboseLevel;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class NetworkAvgUtilTable extends AvgTable {
    List<String> classNames;
    List<String> stationNames;

    public NetworkAvgUtilTable(List<Double> Metricval) {
        super(new ArrayList<>(Arrays.asList(Metricval)));
    }

    public List<String> getClassNames() {
        return classNames;
    }

    public List<String> getStationNames() {
        return stationNames;
    }

    public void setClassNames(List<String> classNames) {
        this.classNames = classNames;
    }

    public void setStationNames(List<String> stationNames) {
        this.stationNames = stationNames;
    }

    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    public List<Double> get(int col) {
        return this.T.get(col);
    }

    public List<Double> getUtil() {
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
                    "Station", "JobClass", "Util");
            System.out.println(
                    "\n-------------------------------");
            DecimalFormat nf = new DecimalFormat("#0.#####");
            nf.setMinimumFractionDigits(5);
            for (int i = 0; i < stationNames.size(); i++) {
                if ( getUtil().get(i) > GlobalConstants.Zero) {
                    System.out.format(
                            "%-10s\t%-10s\t%-10s\n",
                            stationNames.get(i),
                            classNames.get(i),
                            nf.format(getUtil().get(i)));
                }
            }
            System.out.println(
                    "-------------------------------");
        }
    }
}
