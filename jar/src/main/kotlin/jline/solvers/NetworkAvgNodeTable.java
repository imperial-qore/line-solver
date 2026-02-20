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

public class NetworkAvgNodeTable extends AvgTable {
    List<String> classNames;
    List<String> nodeNames;

    public NetworkAvgNodeTable(List<Double> Qval, List<Double> Uval, List<Double> Rval, List<Double> Residval, List<Double> ArvRval, List<Double> Tval) {
        super(new ArrayList<>(Arrays.asList(Qval, Uval, Rval, Residval, ArvRval, Tval)));
    }

    public static NetworkAvgNodeTable tget(NetworkAvgNodeTable T, String anyname) {
        return T.tget(anyname);
    }

    public static NetworkAvgNodeTable tget(NetworkAvgNodeTable T, String stationname, String classname) {
        return T.tget(stationname, classname);
    }

    public List<Double> get(int col) {
        return this.T.getColumn(col).toList1D();
    }

    public List<Double> getArvR() {
        return this.get(4);
    }

    public List<String> getClassNames() {
        return classNames;
    }

    public void setClassNames(List<String> classNames) {
        this.classNames = classNames;
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

    //    List<Double> getArvR() {
//        return this.T.get(5);
//    }

    public List<Double> getUtil() {
        return this.get(1);
    }

    public void print() {
        this.print(this.options);
    }

    public void print(SolverOptions options) {
        if (options.verbose != VerboseLevel.SILENT) {
            // Calculate the maximum length of names for formatting
            int maxStationLength = "Station".length();
            int maxJobClassLength = "JobClass".length();
            for (String name : nodeNames) {
                if (name.length() > maxStationLength) {
                    maxStationLength = name.length();
                }
            }
            for (String name : classNames) {
                if (name.length() > maxJobClassLength) {
                    maxJobClassLength = name.length();
                }
            }

            // Print header with dynamic width for Station and JobClass
            String format = String.format("%%-%ds%%-%ds%%-12s%%-12s%%-12s%%-12s%%-12s%%-12s",
                    maxStationLength + 2, maxJobClassLength + 2);
            System.out.printf(format, "Station", "JobClass", "QLen", "Util", "RespT", "ResidT", "ArvR", "Tput");
            System.out.println();

            // Calculate and print separator line with dynamic length
            int tableWidth = (maxStationLength + 2) + (maxJobClassLength + 2) + 6 * 12;
            for (int i = 0; i < tableWidth; i++) {
                System.out.print("-");
            }
            System.out.println();

            DecimalFormat nf = new DecimalFormat("#0.#####");
            nf.setMinimumFractionDigits(5);

            // Print data rows with the same dynamic width format
            for (int i = 0; i < nodeNames.size(); i++) {
                if (getQLen().get(i) > GlobalConstants.Zero ||
                        getUtil().get(i) > GlobalConstants.Zero ||
                        getRespT().get(i) > GlobalConstants.Zero ||
                        getResidT().get(i) > GlobalConstants.Zero ||
                        getArvR().get(i) > GlobalConstants.Zero ||
                        getTput().get(i) > GlobalConstants.Zero) {

                    System.out.format(format + "\n",
                            nodeNames.get(i),
                            classNames.get(i),
                            nf.format(getQLen().get(i)),
                            nf.format(getUtil().get(i)),
                            nf.format(getRespT().get(i)),
                            nf.format(getResidT().get(i)),
                            nf.format(getArvR().get(i)),
                            nf.format(getTput().get(i)));
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

    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    public NetworkAvgNodeTable tget(Node node, JobClass jobclass) {
        return this.tget(node.getName(), jobclass.getName());
    }

    public NetworkAvgNodeTable tget(String anyname) {
        int rowIdx = this.nodeNames.indexOf(anyname);
        int colIdx = this.classNames.indexOf(anyname);
        if (rowIdx < 0 && colIdx < 0) {
            NetworkAvgNodeTable filteredAvgTable = new NetworkAvgNodeTable(new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
            filteredAvgTable.setOptions(this.options);
            filteredAvgTable.setNodeNames(new ArrayList<String>());
            filteredAvgTable.setClassNames(new ArrayList<String>());
            return filteredAvgTable;
        }

        List<Integer> indexToKeep = new ArrayList<>();
        List<String> nodeNamesFilt = new ArrayList<>();
        List<String> classNamesFilt = new ArrayList<>();
        for (int i = 0; i < this.nodeNames.size(); i++) {
            if (this.nodeNames.get(i).compareTo(anyname) == 0 || this.classNames.get(i).compareTo(anyname) == 0) {
                nodeNamesFilt.add(this.nodeNames.get(i));
                classNamesFilt.add(this.classNames.get(i));
                indexToKeep.add(i);
            }
        }

        List<Double> Qval = this.getQLen();
        List<Double> QvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Qval.size()) {
                QvalFilt.add(Qval.get(index));
            }
        }
        List<Double> Uval = this.getUtil();
        List<Double> UvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Uval.size()) {
                UvalFilt.add(Uval.get(index));
            }
        }
        List<Double> Rval = this.getRespT();
        List<Double> RvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Rval.size()) {
                RvalFilt.add(Rval.get(index));
            }
        }
        List<Double> Wval = this.getResidT();
        List<Double> WvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Wval.size()) {
                WvalFilt.add(Wval.get(index));
            }
        }
        List<Double> Aval = this.getArvR();
        List<Double> AvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Aval.size()) {
                AvalFilt.add(Aval.get(index));
            }
        }
        List<Double> Tval = this.getTput();
        List<Double> TvalFilt = new ArrayList<>();
        for (Integer index : indexToKeep) {
            if (index >= 0 && index < Tval.size()) {
                TvalFilt.add(Tval.get(index));
            }
        }
        NetworkAvgNodeTable filteredAvgTable = new NetworkAvgNodeTable(QvalFilt, UvalFilt, RvalFilt, WvalFilt, AvalFilt, TvalFilt);
        filteredAvgTable.setOptions(this.options);
        filteredAvgTable.setNodeNames(nodeNamesFilt);
        filteredAvgTable.setClassNames(classNamesFilt);
        return filteredAvgTable;
    }

    public NetworkAvgNodeTable tget(String nodename, String classname) {
        int rowIdx = this.nodeNames.indexOf(nodename);
        int colIdx = this.classNames.indexOf(classname);
        if (rowIdx < 0 || colIdx < 0) {
            NetworkAvgNodeTable filteredAvgTable = new NetworkAvgNodeTable(new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
            filteredAvgTable.setOptions(this.options);
            filteredAvgTable.setNodeNames(new ArrayList<String>());
            filteredAvgTable.setClassNames(new ArrayList<String>());
            return filteredAvgTable;
        }
        int indexToKeep = 0;
        for (int i = 0; i < this.nodeNames.size(); i++) {
            if (this.nodeNames.get(i).compareTo(nodename) == 0 && this.classNames.get(i).compareTo(classname) == 0) {
                indexToKeep = i;
                break;
            }
        }
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
        NetworkAvgNodeTable filteredAvgTable = new NetworkAvgNodeTable(Qval, Uval, Rval, Residval, ArvRval, Tval);
        filteredAvgTable.setOptions(this.options);
        filteredAvgTable.setNodeNames(Collections.singletonList(this.nodeNames.get(rowIdx)));
        filteredAvgTable.setClassNames(Collections.singletonList(this.classNames.get(colIdx)));
        return filteredAvgTable;
    }
}
