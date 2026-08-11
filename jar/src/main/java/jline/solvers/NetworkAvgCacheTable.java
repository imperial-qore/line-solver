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
 * Table of detailed per-class cache performance metrics, one total row per cache
 * node and read class (List=0) plus, where available and the cache has more than
 * one list, one row per cache list (level). Columns: List, ListCap, Items,
 * HitProb, DelayedHitProb, MissProb, HitRate, DelayedHitRate, MissRate, ArvR,
 * ResidT. Port of matlab/src/solvers/@NetworkSolver/getAvgCacheTable.m.
 */
public class NetworkAvgCacheTable extends AvgTable {
    List<String> nodeNames;
    List<String> classNames;

    public NetworkAvgCacheTable(List<Double> List_, List<Double> ListCap, List<Double> Items,
                                List<Double> HitProb, List<Double> DelayedHitProb, List<Double> MissProb,
                                List<Double> HitRate, List<Double> DelayedHitRate, List<Double> MissRate,
                                List<Double> ArvR, List<Double> ResidT) {
        super(new ArrayList<>(Arrays.asList(List_, ListCap, Items, HitProb, DelayedHitProb, MissProb,
                HitRate, DelayedHitRate, MissRate, ArvR, ResidT)));
    }

    public List<Double> get(int col) { return this.T.getColumn(col).toList1D(); }
    public List<Double> getList() { return get(0); }
    public List<Double> getListCap() { return get(1); }
    public List<Double> getItems() { return get(2); }
    public List<Double> getHitProb() { return get(3); }
    public List<Double> getDelayedHitProb() { return get(4); }
    public List<Double> getMissProb() { return get(5); }
    public List<Double> getHitRate() { return get(6); }
    public List<Double> getDelayedHitRate() { return get(7); }
    public List<Double> getMissRate() { return get(8); }
    public List<Double> getArvR() { return get(9); }
    public List<Double> getResidT() { return get(10); }

    public List<String> getNodeNames() { return nodeNames; }
    public void setNodeNames(List<String> nodeNames) { this.nodeNames = nodeNames; }
    public List<String> getClassNames() { return classNames; }
    public void setClassNames(List<String> classNames) { this.classNames = classNames; }
    public void setOptions(SolverOptions options) { this.options = options; }

    @Override
    public void print() { this.print(this.options); }

    public void print(SolverOptions options) {
        if (options != null && options.verbose == VerboseLevel.SILENT) return;
        if (nodeNames == null || nodeNames.isEmpty()) return;
        String[] headers = {"Node", "JobClass", "List", "ListCap", "Items", "HitProb",
                "DelayedHitProb", "MissProb", "HitRate", "DelayedHitRate", "MissRate", "ArvR", "ResidT"};
        List<String[]> rows = new ArrayList<>();
        List<Double> li = getList(), lc = getListCap(), it = getItems(), hp = getHitProb(),
                dhp = getDelayedHitProb(), mp = getMissProb(), hr = getHitRate(), dhr = getDelayedHitRate(),
                mr = getMissRate(), ar = getArvR(), lat = getResidT();
        for (int i = 0; i < nodeNames.size(); i++) {
            rows.add(new String[]{
                    nodeNames.get(i), classNames.get(i),
                    Integer.toString((int) Math.round(li.get(i))),
                    fmtValue(lc.get(i), 5),
                    Integer.toString((int) Math.round(it.get(i))),
                    fmtValue(hp.get(i), 5), fmtValue(dhp.get(i), 5), fmtValue(mp.get(i), 5),
                    fmtValue(hr.get(i), 5), fmtValue(dhr.get(i), 5), fmtValue(mr.get(i), 5),
                    fmtValue(ar.get(i), 5), fmtValue(lat.get(i), 5)});
        }
        printFormattedTable(headers, rows);
    }

    public void printTable() { this.print(); }
    public void printTable(SolverOptions options) { this.print(options); }
}
