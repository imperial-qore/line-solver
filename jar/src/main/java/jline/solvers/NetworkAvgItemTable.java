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
 * Table of item-level cache occupancy, one row per Cache node, item and cache
 * list (level). Columns: Item, List, ListCap, Prob, DelayedHitQLen,
 * DelayedHitQLenFull (Node/class-less names held separately). Prob is the
 * steady-state probability that the item resides in that list;
 * DelayedHitQLen is the mean number of secondary requests waiting on the
 * in-flight fetch of the item (exact under SolverCTMC, NaN where not
 * computed) and DelayedHitQLenFull additionally counts the request that
 * triggered the fetch. Port of matlab getAvgItemTable.m.
 */
public class NetworkAvgItemTable extends AvgTable {
    List<String> nodeNames;

    public NetworkAvgItemTable(List<Double> Item, List<Double> List_, List<Double> ListCap, List<Double> Prob,
                               List<Double> DelayedHitQLen, List<Double> DelayedHitQLenFull) {
        super(new ArrayList<>(Arrays.asList(Item, List_, ListCap, Prob, DelayedHitQLen, DelayedHitQLenFull)));
    }

    public List<Double> get(int col) { return this.T.getColumn(col).toList1D(); }
    public List<Double> getItem() { return get(0); }
    public List<Double> getList() { return get(1); }
    public List<Double> getListCap() { return get(2); }
    public List<Double> getProb() { return get(3); }
    public List<Double> getDelayedHitQLen() { return get(4); }
    public List<Double> getDelayedHitQLenFull() { return get(5); }

    public List<String> getNodeNames() { return nodeNames; }
    public void setNodeNames(List<String> nodeNames) { this.nodeNames = nodeNames; }
    public void setOptions(SolverOptions options) { this.options = options; }

    @Override
    public void print() { this.print(this.options); }

    public void print(SolverOptions options) {
        if (options != null && options.verbose == VerboseLevel.SILENT) return;
        if (nodeNames == null || nodeNames.isEmpty()) return;
        String[] headers = {"Node", "Item", "List", "ListCap", "Prob", "DelayedHitQLen", "DelayedHitQLenFull"};
        List<String[]> rows = new ArrayList<>();
        List<Double> it = getItem(), li = getList(), lc = getListCap(), pr = getProb();
        List<Double> dq = getDelayedHitQLen(), dqf = getDelayedHitQLenFull();
        for (int i = 0; i < nodeNames.size(); i++) {
            rows.add(new String[]{
                    nodeNames.get(i),
                    Integer.toString((int) Math.round(it.get(i))),
                    Integer.toString((int) Math.round(li.get(i))),
                    fmtValue(lc.get(i), 5),
                    fmtValue(pr.get(i), 5),
                    fmtValue(dq.get(i), 5),
                    fmtValue(dqf.get(i), 5)});
        }
        printFormattedTable(headers, rows);
    }

    public void printTable() { this.print(); }
    public void printTable(SolverOptions options) { this.print(options); }
}
