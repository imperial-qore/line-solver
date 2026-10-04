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
 * Table of loss (drop) metrics for every station-class pair that receives
 * offered traffic: offered arrival rate (ArvR), carried throughput (Tput), loss
 * rate (ArvR - Tput, the rate of jobs dropped by finite capacity, blocking, or
 * reneging) and loss ratio (LossRate / ArvR).
 *
 * <p>Only pairs with ArvR &gt; 0 are listed, which excludes the Source (whose
 * offered arrival rate is zero); a lossless station has ArvR = Tput and so
 * LossRate = LossRatio = 0. Port of MATLAB getAvgLossTable.m.</p>
 */
public class NetworkLossTable extends AvgTable {
    List<String> stationNames;
    List<String> classNames;
    String firstColumnName = "Station";

    public NetworkLossTable(List<Double> ArvR, List<Double> Tput, List<Double> LossRate, List<Double> LossRatio) {
        super(new ArrayList<>(Arrays.asList(ArvR, Tput, LossRate, LossRatio)));
    }

    public void setFirstColumnName(String name) {
        this.firstColumnName = name;
    }

    public List<Double> get(int col) {
        return this.T.getColumn(col).toList1D();
    }

    public List<Double> getArvR() {
        return get(0);
    }

    public List<Double> getTput() {
        return get(1);
    }

    public List<Double> getLossRate() {
        return get(2);
    }

    public List<Double> getLossRatio() {
        return get(3);
    }

    public List<String> getStationNames() {
        return stationNames;
    }

    public void setStationNames(List<String> stationNames) {
        this.stationNames = stationNames;
    }

    public List<String> getClassNames() {
        return classNames;
    }

    public void setClassNames(List<String> classNames) {
        this.classNames = classNames;
    }

    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    @Override
    public void print() {
        this.print(this.options);
    }

    public void print(SolverOptions options) {
        if (options != null && options.verbose == VerboseLevel.SILENT) return;
        if (stationNames == null || stationNames.isEmpty()) return;
        String[] headers = {firstColumnName, "JobClass", "ArvR", "Tput", "LossRate", "LossRatio"};
        List<String[]> rows = new ArrayList<>();
        List<Double> a = getArvR();
        List<Double> t = getTput();
        List<Double> lr = getLossRate();
        List<Double> lc = getLossRatio();
        for (int i = 0; i < stationNames.size(); i++) {
            rows.add(new String[]{
                    stationNames.get(i),
                    classNames.get(i),
                    fmtValue(a.get(i), 5),
                    fmtValue(t.get(i), 5),
                    fmtValue(lr.get(i), 5),
                    fmtValue(lc.get(i), 5)});
        }
        printFormattedTable(headers, rows);
    }

    public void printTable() {
        this.print();
    }

    public void printTable(SolverOptions options) {
        this.print(options);
    }
}
