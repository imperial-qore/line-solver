/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.lang.JobClass;
import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.nodes.Station;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class NetworkAvgChainTable extends AvgTable {
    List<String> chainNames = new ArrayList<>();
    List<String> stationNames = new ArrayList<>();
    List<String> inChainNames = new ArrayList<>();
    int nDigits = 5;

    public NetworkAvgChainTable(List<Double> Qval, List<Double> Uval, List<Double> Rval, List<Double> Residval, List<Double> ArvRval, List<Double> Tval) {
        super(new ArrayList<>(Arrays.asList(Qval, Uval, Rval, Residval, ArvRval, Tval)));
    }

    public List<Double> get(int col) {
        return this.T.getColumn(col).toList1D();
    }

    public List<Double> getArvR() {
        return this.get(4);
    }

    public List<String> getChainNames() {
        return chainNames;
    }

    public void setChainNames(List<String> chainNames) {
        this.chainNames = chainNames;
    }

    public java.util.List<String> getInChainNames() {
        return inChainNames;
    }

    public void setInChainNames(java.util.List<String> inChainNames) {
        this.inChainNames = inChainNames;
    }

    public List<Double> getQLen() {
        return this.get(0);
    }

    public List<Double> getResidT() {
        return this.get(3);
    }

    public List<Double> getRespT() {
        return this.get(2);
    }

    public List<String> getStationNames() {
        return stationNames;
    }

    public void setStationNames(List<String> stationNames) {
        this.stationNames = stationNames;
    }

    public List<Double> getTput() {
        return this.get(5);
    }

    public List<Double> getUtil() {
        return this.get(1);
    }

    public void print() {
        this.print(this.options);
    }

    public void print(SolverOptions options) {
        if (options.verbose != VerboseLevel.SILENT) {
            String[] headers = {"Station", "Chain", "JobClasses", "QLen", "Util", "RespT", "ResidT", "ArvR", "Tput"};
            List<String[]> rows = new ArrayList<>();
            for (int i = 0; i < stationNames.size(); i++) {
                if (getQLen().get(i) > GlobalConstants.Zero ||
                        getUtil().get(i) > GlobalConstants.Zero ||
                        getRespT().get(i) > GlobalConstants.Zero ||
                        getResidT().get(i) > GlobalConstants.Zero ||
                        getArvR().get(i) > GlobalConstants.Zero ||
                        getTput().get(i) > GlobalConstants.Zero) {
                    rows.add(new String[]{
                            stationNames.get(i),
                            chainNames.get(i),
                            inChainNames.get(i),
                            fmtValue(getQLen().get(i), nDigits),
                            fmtValue(getUtil().get(i), nDigits),
                            fmtValue(getRespT().get(i), nDigits),
                            fmtValue(getResidT().get(i), nDigits),
                            fmtValue(getArvR().get(i), nDigits),
                            fmtValue(getTput().get(i), nDigits)});
                }
            }
            printFormattedTable(headers, rows);
        }
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

    public void setNumberOfDigits(int nd) {
        this.nDigits = nd;
    }

    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    public AvgTable tget(Station station, JobClass jobclass) {
        return this.tget(station.getName(), jobclass.getName());
    }

    public AvgTable tget(String stationname, String classname) {
        int rowIdx = this.stationNames.indexOf(stationname);
        int colIdx = this.chainNames.indexOf(classname);
        if (rowIdx < 0 || colIdx < 0) {
            NetworkAvgChainTable filteredAvgTable = new NetworkAvgChainTable(new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
            filteredAvgTable.setOptions(this.options);
            filteredAvgTable.setStationNames(new ArrayList<String>());
            filteredAvgTable.setChainNames(new ArrayList<String>());
            return filteredAvgTable;
        }
        int indexToKeep = rowIdx * this.chainNames.size() + colIdx - 1;
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
        NetworkAvgChainTable filteredAvgTable = new NetworkAvgChainTable(Qval, Uval, Rval, Residval, ArvRval, Tval);
        filteredAvgTable.setOptions(this.options);
        filteredAvgTable.setStationNames(Collections.singletonList(this.stationNames.get(rowIdx)));
        filteredAvgTable.setChainNames(Collections.singletonList(this.chainNames.get(colIdx)));
        return filteredAvgTable;
    }

}
