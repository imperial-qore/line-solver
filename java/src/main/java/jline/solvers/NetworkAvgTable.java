package jline.solvers;

import java.util.*;

import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.nodes.Station;

import java.util.stream.Collectors;
import jline.lang.constant.GlobalConstants;
import jline.lang.constant.VerboseLevel;
import jline.solvers.jmt.SolverJMT;

import javax.xml.parsers.ParserConfigurationException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;

public class NetworkAvgTable extends AvgTable {
    List<String> classNames;
    List<String> stationNames;

    public NetworkAvgTable(List<Double> Qval, List<Double> Uval, List<Double> Rval, List<Double> Residval, List<Double> ArvRval, List<Double> Tval) {
        super(new ArrayList<>(Arrays.asList(Qval, Uval, Rval, Residval, ArvRval, Tval)));
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
        return new ArrayList<>(this.T.get(col));
    }

    public List<Double> getQLen() {
        return this.get(0);
    }

    public List<Double> getUtil() {
        return this.get(1);
    }

    public List<Double> getRespT() {
        return this.get(2);
    }

    public List<Double> getResidT() {
        return this.get(3);
    }

    public List<Double> getArvR() {
        return this.get(4);
    }

    public List<Double> getTput() {
        return this.get(5);
    }

    //    List<Double> getArvR() {
//        return this.T.get(5);
//    }

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
                    "%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s",
                    "Station", "JobClass", "QLen", "Util", "RespT", "ResidT", "ArvR", "Tput");
            System.out.println(
                    "\n--------------------------------------------------------------------------------------------");
            DecimalFormat nf = new DecimalFormat("#0.#####");
            nf.setMinimumFractionDigits(5);
            for (int i = 0; i < stationNames.size(); i++) {
                if (getQLen().get(i) > GlobalConstants.Zero ||
                        getUtil().get(i) > GlobalConstants.Zero ||
                        getRespT().get(i) > GlobalConstants.Zero ||
                        getResidT().get(i) > GlobalConstants.Zero ||
                        getArvR().get(i) > GlobalConstants.Zero ||
                        getTput().get(i) > GlobalConstants.Zero) {
                    System.out.format(
                            "%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\n",
                            stationNames.get(i),
                            classNames.get(i),
                            nf.format(getQLen().get(i)),
                            nf.format(getUtil().get(i)),
                            nf.format(getRespT().get(i)),
                            nf.format(getResidT().get(i)),
                            nf.format(getArvR().get(i)),
                            nf.format(getTput().get(i)));
                }
            }
            System.out.println(
                    "--------------------------------------------------------------------------------------------");
        }
    }

    public NetworkAvgTable tget(Station station, JobClass jobclass) {
        return this.tget((String)station.getName(), (String)jobclass.getName());
    }

    public NetworkAvgTable tget(String stationname, String classname) {
        int rowIdx = this.stationNames.indexOf(stationname);
        int colIdx = this.classNames.indexOf(classname);
        if (rowIdx < 0 || colIdx <0) {
            NetworkAvgTable filteredAvgTable = new NetworkAvgTable(new ArrayList<>(),new ArrayList<>(),new ArrayList<>(),new ArrayList<>(),new ArrayList<>(),new ArrayList<>());
            filteredAvgTable.setOptions(this.options);
            filteredAvgTable.setStationNames(new ArrayList<String>());
            filteredAvgTable.setClassNames(new ArrayList<String>());
            return filteredAvgTable;
        }
        int indexToKeep = rowIdx * this.classNames.size() + colIdx - 1;
        List<Double> Qval = this.getQLen();
        if (indexToKeep > 0) {Qval.subList(0, indexToKeep).clear();}
        if (indexToKeep < Qval.size()) {Qval.subList(1, Qval.size()).clear();}
        List<Double> Uval = this.getUtil();
        if (indexToKeep > 0) {Uval.subList(0, indexToKeep).clear();}
        if (indexToKeep < Uval.size()) {Uval.subList(1, Uval.size()).clear();}
        List<Double> Rval = this.getRespT();
        if (indexToKeep > 0) {Rval.subList(0, indexToKeep).clear();}
        if (indexToKeep < Rval.size()) {Rval.subList(1, Rval.size()).clear();}
        List<Double> Residval = this.getResidT();
        if (indexToKeep > 0) {Residval.subList(0, indexToKeep).clear();}
        if (indexToKeep < Residval.size()) {Residval.subList(1, Residval.size()).clear();}
        List<Double> ArvRval = this.getArvR();
        if (indexToKeep > 0) {ArvRval.subList(0, indexToKeep).clear();}
        if (indexToKeep < ArvRval.size()) {ArvRval.subList(1, ArvRval.size()).clear();}
        List<Double> Tval = this.getTput();
        if (indexToKeep > 0) {Tval.subList(0, indexToKeep).clear();}
        if (indexToKeep < Tval.size()) {Tval.subList(1, Tval.size()).clear();}
        NetworkAvgTable filteredAvgTable = new NetworkAvgTable(Qval, Uval, Rval, Residval, ArvRval, Tval);
        filteredAvgTable.setOptions(this.options);
        filteredAvgTable.setStationNames(Collections.singletonList(this.stationNames.get(rowIdx)));
        filteredAvgTable.setClassNames(Collections.singletonList(this.classNames.get(colIdx)));
        return filteredAvgTable;
    }

}
