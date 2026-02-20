/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.lang.JobClass;
import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.nodes.Node;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class NetworkAvgNodeChainTable extends AvgTable {
    List<String> chainNames;
    List<String> nodeNames;
    List<String> inChainNames;
    int nDigits = 5;

    public NetworkAvgNodeChainTable(List<Double> Qval, List<Double> Uval, List<Double> Rval, List<Double> Residval, List<Double> ArvRval, List<Double> Tval) {
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

    public List<String> getInChainNames() {
        return inChainNames;
    }

    public void setInChainNames(List<String> inChainNames) {
        this.inChainNames = inChainNames;
    }

    public List<String> getNodeNames() {
        return nodeNames;
    }

    public void setNodeNames(List<String> nodeNames) {
        this.nodeNames = nodeNames;
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
            // Calculate the maximum length of names for formatting

            int maxNodeLength = "Node".length();
            int maxChainLength = "Chain".length();
            int maxInChainLength = "JobClasses".length();
            for (String name : nodeNames) {
                if (name.length() > maxNodeLength) {
                    maxNodeLength = name.length();
                }
            }
            for (String name : chainNames) {
                if (name.length() > maxChainLength) {
                    maxChainLength = name.length();
                }
            }
            for (String name : inChainNames) {
                if (name.length() > maxChainLength) {
                    maxInChainLength = name.length();
                }
            }


            // Print header with dynamic width for Node and JobClass
            String format = String.format("%%-%ds%%-%ds%%-%ds%%-12s%%-12s%%-12s%%-12s%%-12s%%-12s",
                    maxNodeLength + 2, maxChainLength + 2, maxInChainLength + 2);
            System.out.printf(format, "Node", "Chain", "JobClasses", "QLen", "Util", "RespT", "ResidT", "ArvR", "Tput");
            System.out.println();

            // Calculate and print separator line with dynamic length
            DecimalFormat nf = new DecimalFormat("#0.#####");

            nf.setMinimumFractionDigits(nDigits);
            int tableWidth = (maxNodeLength + 2) + (maxChainLength + 2) + (maxInChainLength + 2) + 6 * (12 - 5 + nDigits);
            for (int i = 0; i < tableWidth; i++) {
                System.out.print("-");
            }
            System.out.println();

            // Print data rows with the same dynamic width format
            for (int row = 0; row < nodeNames.size(); row++) {
                if (getQLen().get(row) > GlobalConstants.Zero ||
                        getUtil().get(row) > GlobalConstants.Zero ||
                        getRespT().get(row) > GlobalConstants.Zero ||
                        getResidT().get(row) > GlobalConstants.Zero ||
                        getArvR().get(row) > GlobalConstants.Zero ||
                        getTput().get(row) > GlobalConstants.Zero) {

                    System.out.format(format + "\n",
                            nodeNames.get(row),
                            chainNames.get(row),
                            inChainNames.get(row),
                            nf.format(getQLen().get(row)),
                            nf.format(getUtil().get(row)),
                            nf.format(getRespT().get(row)),
                            nf.format(getResidT().get(row)),
                            nf.format(getArvR().get(row)),
                            nf.format(getTput().get(row)));
                }
            }

            // Print separator line again at the end
            for (int i = 0; i < tableWidth; i++) {
                System.out.print("-");
            }
            System.out.println();
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

    public AvgTable tget(Node node, JobClass jobclass) {
        return this.tget(node.getName(), jobclass.getName());
    }

    public AvgTable tget(String nodename, String classname) {
        int rowIdx = this.nodeNames.indexOf(nodename);
        int colIdx = this.chainNames.indexOf(classname);
        if (rowIdx < 0 || colIdx < 0) {
            NetworkAvgNodeChainTable filteredAvgTable = new NetworkAvgNodeChainTable(new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
            filteredAvgTable.setOptions(this.options);
            filteredAvgTable.setNodeNames(new ArrayList<String>());
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
        NetworkAvgNodeChainTable filteredAvgTable = new NetworkAvgNodeChainTable(Qval, Uval, Rval, Residval, ArvRval, Tval);
        filteredAvgTable.setOptions(this.options);
        filteredAvgTable.setNodeNames(Collections.singletonList(this.nodeNames.get(rowIdx)));
        filteredAvgTable.setChainNames(Collections.singletonList(this.chainNames.get(colIdx)));
        return filteredAvgTable;
    }

}
