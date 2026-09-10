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
 * Table of the mean orbit length of every retrial station-class pair, with the
 * station population and the in-service population it decomposes into.
 *
 * <p>Reported as a separate table rather than as an extra column of getAvgTable
 * so that the average table keeps its shape for models without retrials. Port of
 * MATLAB getAvgOrbitTable.m.</p>
 */
public class NetworkAvgOrbitTable extends AvgTable {
    List<String> stationNames;
    List<String> classNames;

    public NetworkAvgOrbitTable(List<Double> QLen, List<Double> InService, List<Double> Orbit) {
        super(new ArrayList<>(Arrays.asList(QLen, InService, Orbit)));
    }

    public List<Double> get(int col) {
        return this.T.getColumn(col).toList1D();
    }

    public List<Double> getQLen() {
        return get(0);
    }

    public List<Double> getInService() {
        return get(1);
    }

    public List<Double> getOrbit() {
        return get(2);
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
        String[] headers = {"Station", "JobClass", "QLen", "InService", "Orbit"};
        List<String[]> rows = new ArrayList<>();
        List<Double> q = getQLen();
        List<Double> s = getInService();
        List<Double> o = getOrbit();
        for (int i = 0; i < stationNames.size(); i++) {
            rows.add(new String[]{
                    stationNames.get(i),
                    classNames.get(i),
                    fmtValue(q.get(i), 5),
                    fmtValue(s.get(i), 5),
                    fmtValue(o.get(i), 5)});
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
