/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.GlobalConstants;
import jline.VerboseLevel;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LayeredNetworkAvgTable extends AvgTable {
    List<String> nodeNames;
    List<String> nodeTypes;
    boolean printZerosDefault = true; // by default, print zero results

    public LayeredNetworkAvgTable(List<Double> Qval, List<Double> Uval, List<Double> Rval, List<Double> Residval, List<Double> Aval, List<Double> Tval) {
        super(new ArrayList<>(Arrays.asList(Qval, Uval, Rval, Residval, Aval, Tval)));
    }

    private String formatValue(double value, DecimalFormat nf) {
        if (value == 0.0) {
            return "0";
        } else {
            return nf.format(value);
        }
    }

    public List<Double> get(int col) {
        return this.T.getColumn(col).toList1D();
    }

    public List<Double> getArvR() {
        return this.T.getColumn(4).toList1D();
    }

    public java.util.List<String> getNodeNames() {
        return nodeNames;
    }

    public void setNodeNames(java.util.List<String> nodeNames) {
        this.nodeNames = nodeNames;
    }

    public java.util.List<String> getNodeTypes() {
        return nodeTypes;
    }

    public void setNodeTypes(java.util.List<String> nodeTypes) {
        this.nodeTypes = nodeTypes;
    }

    public List<Double> getQLen() {
        return this.T.getColumn(0).toList1D();
    }

    public List<Double> getResidT() {
        return this.T.getColumn(3).toList1D();
    }

    public List<Double> getRespT() {
        return this.T.getColumn(2).toList1D();
    }

    public List<Double> getTput() {
        return this.T.getColumn(5).toList1D();
    }

    public List<Double> getUtil() {
        return this.T.getColumn(1).toList1D();
    }

    public void print() {
        this.print(this.options, printZerosDefault);
    }

    public void print(boolean printZeros) {
        this.print(this.options, printZeros);
    }

    public void print(SolverOptions options) {
        print(options, printZerosDefault);
    }

    public void print(SolverOptions options, boolean printZeros) {
        if (options.verbose != VerboseLevel.SILENT) {
            System.out.printf(
                    "%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s",
                    "Node", "NodeType", "QLen", "Util", "RespT", "ResidT", "ArvR", "Tput");
            System.out.println(
                    "\n--------------------------------------------------------------------------------------------");
            DecimalFormat nf = new DecimalFormat("#0.#####");
            nf.setMinimumFractionDigits(5);
            for (int i = 0; i < nodeTypes.size(); i++) {
                if (printZeros || (getQLen().get(i) > GlobalConstants.Zero ||
                        getUtil().get(i) > GlobalConstants.Zero ||
                        getRespT().get(i) > GlobalConstants.Zero ||
                        getResidT().get(i) > GlobalConstants.Zero ||
                        getArvR().get(i) > GlobalConstants.Zero ||
                        getTput().get(i) > GlobalConstants.Zero)) {
                    System.out.format(
                            "%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\n",
                            nodeNames.get(i),
                            nodeTypes.get(i),
                            formatValue(getQLen().get(i), nf),
                            formatValue(getUtil().get(i), nf),
                            formatValue(getRespT().get(i), nf),
                            formatValue(getResidT().get(i), nf),
                            formatValue(getArvR().get(i), nf),
                            formatValue(getTput().get(i), nf));
                }
            }
            System.out.println(
                    "--------------------------------------------------------------------------------------------");
        }
    }

    public void printTable() {
        this.print(this.options.verbose(VerboseLevel.STD));
    }

    public void printTable(SolverOptions options) {
        this.print(options);
    }

    public void setOptions(SolverOptions options) {
        this.options = options;
    }
}
