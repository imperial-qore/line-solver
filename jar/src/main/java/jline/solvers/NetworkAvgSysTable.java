/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.lang.Chain;
import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.util.matrix.Matrix;

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
            String[] headers = {"Chain", "JobClasses", "SysRespT", "SysTput"};
            List<String[]> rows = new ArrayList<>();
            for (int row = 0; row < chainNames.size(); row++) {
                if (getSysRespT().get(row) > GlobalConstants.Zero ||
                        getSysTput().get(row) > GlobalConstants.Zero) {
                    rows.add(new String[]{
                            chainNames.get(row),
                            inChainNames.get(row),
                            fmtValue(getSysRespT().get(row), nDigits),
                            fmtValue(getSysTput().get(row), nDigits)});
                }
            }
            printFormattedTable(headers, rows);
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

    /**
     * The rows of the named chain.
     *
     * @param chainname the chain name
     * @return the rows for that chain
     * @throws IllegalArgumentException if the table has no such chain
     */
    public NetworkAvgSysTable tget(String chainname) {
        int colIdx = this.chainNames.indexOf(chainname);
        if (colIdx < 0) {
            // A name matching nothing is a typo, not an empty result. Returning an
            // empty table made the two indistinguishable: the caller printed a bare
            // header and had no way to tell which had happened.
            List<String> chains = new ArrayList<String>();
            for (int i = 0; i < this.chainNames.size(); i++) {
                if (!chains.contains(this.chainNames.get(i))) {
                    chains.add(this.chainNames.get(i));
                }
            }
            throw new IllegalArgumentException("no chain named '" + chainname
                    + "' in this table. Chains: " + chains + ".");
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
