/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.lang.Chain;
import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.util.matrix.Matrix;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class NetworkAvgSysTable extends AvgTable {

    SolverOptions options;
    Matrix T;
    List<String> chainNames;
    List<String> inChainNames;
    int nDigits = 5;

    public NetworkAvgSysTable(List<Double> SysRespTval, List<Double> SysTputval, SolverOptions options) {
        super(new ArrayList<>(Arrays.asList(SysRespTval, SysTputval)));
        this.T = new Matrix(new ArrayList<>(Arrays.asList(SysRespTval, SysTputval)));
        this.options = options;
    }

    public static NetworkAvgSysTable tget(NetworkAvgSysTable T, String chainname) {
        return T.tget(chainname);
    }

    public List<Double> get(int col) {
        return this.T.getColumn(col).toList1D();
    }

    public java.util.List<String> getChainNames() {
        return chainNames;
    }

    public void setChainNames(java.util.List<String> chainNames) {
        this.chainNames = chainNames;
    }

    public java.util.List<String> getInChainNames() {
        return inChainNames;
    }

    public void setInChainNames(java.util.List<String> inChainNames) {
        this.inChainNames = inChainNames;
    }

    public List<Double> getSysRespT() {
        return this.T.getColumn(0).toList1D();
    }

    public List<Double> getSysTput() {
        return this.T.getColumn(1).toList1D();
    }

    public void print() {
        this.print(this.options);
    }

    public void print(SolverOptions options) {
        if (options.verbose != VerboseLevel.SILENT) {
            // Calculate the maximum length of names for formatting

            int maxChainLength = "Chain".length();
            int maxInChainLength = "JobClasses".length();
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

            // Print header with dynamic width for Chain and JobClass
            String format = String.format("%%-%ds%%-%ds%%-%ds%%-%ds",
                    maxChainLength + 2, maxInChainLength + 10, 12, 12);

            System.out.printf(format, "Chain", "JobClasses", "SysRespT", "SysTput");
            System.out.println();

            // Calculate and print separator line with dynamic length
            DecimalFormat nf = new DecimalFormat("#0.#####");
            nf.setMinimumFractionDigits(nDigits);
            int tableWidth = (maxChainLength + 2) + (maxInChainLength + 6) + 2 * (12 - 5 + nDigits);
            for (int i = 0; i < tableWidth; i++) {
                System.out.print("-");
            }
            System.out.println();

            // Print data rows with the same dynamic width format
            for (int row = 0; row < chainNames.size(); row++) {
                if (getSysRespT().get(row) > GlobalConstants.Zero ||
                        getSysTput().get(row) > GlobalConstants.Zero) {
                    System.out.format(format + "\n",
                            chainNames.get(row),
                            inChainNames.get(row),
                            nf.format(getSysRespT().get(row)),
                            nf.format(getSysTput().get(row)));
                }
            }

            // Print separator line again at the end
            for (int i = 0; i < tableWidth; i++) {
                System.out.print("-");
            }
            System.out.println();
        }
    }

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

    public NetworkAvgSysTable tget(Chain chain) {
        return this.tget(chain.getName());
    }

    public NetworkAvgSysTable tget(String chainname) {
        int colIdx = this.chainNames.indexOf(chainname);
        if (colIdx < 0) {
            NetworkAvgSysTable filteredAvgTable = new NetworkAvgSysTable(new ArrayList<>(), new ArrayList<>(), this.options);
            filteredAvgTable.setOptions(this.options);
            filteredAvgTable.setChainNames(new ArrayList<String>());
            filteredAvgTable.setInChainNames(new ArrayList<String>());
            return filteredAvgTable;
        }

        List<Integer> indexToKeep = new ArrayList<>();
        List<String> chainNamesFilt = new ArrayList<>();
        List<String> inChainNamesFilt = new ArrayList<>();
        for (int i = 0; i < this.chainNames.size(); i++) {
            if (this.chainNames.get(i).compareTo(chainname) == 0) {
                chainNamesFilt.add(this.chainNames.get(i));
                inChainNamesFilt.add(this.inChainNames.get(i));
                indexToKeep.add(i);
            }
        }

        List<Double> Rval = this.getSysRespT();
        List<Double> RvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Rval.size()) {
                RvalFilt.add(Rval.get(index));
            }
        }

        List<Double> Tval = this.getSysTput();
        List<Double> TvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Tval.size()) {
                TvalFilt.add(Tval.get(index));
            }
        }
        NetworkAvgSysTable filteredAvgTable = new NetworkAvgSysTable(RvalFilt, TvalFilt, this.options);
        filteredAvgTable.setOptions(this.options);
        filteredAvgTable.setChainNames(chainNamesFilt);
        filteredAvgTable.setInChainNames(inChainNamesFilt);
        return filteredAvgTable;
    }
}
