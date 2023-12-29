package jline.solvers;

import jline.lang.constant.GlobalConstants;
import jline.lang.constant.VerboseLevel;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LayeredNetworkAvgTable {

    SolverOptions options;
    ArrayList<List<Double>> T;
    List<String> nodeNames;
    List<String> nodeTypes;

    public LayeredNetworkAvgTable(List<Double> Qval, List<Double> Uval, List<Double> Rval, List<Double> Residval, List<Double> Tval) {
        this.T = new ArrayList<>(Arrays.asList(Qval, Uval, Rval, Residval, Tval));
    }

    public java.util.List<String> getNodeNames() {
        return nodeNames;
    }

    public java.util.List<String> getNodeTypes() {
        return nodeTypes;
    }

    public void setNodeNames(java.util.List<String> nodeNames) {
        this.nodeNames = nodeNames;
    }

    public void setNodeTypes(java.util.List<String> nodeTypes) {
        this.nodeTypes = nodeTypes;
    }

    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    public List<Double> get(int col) {
        return this.T.get(col);
    }

    public List<Double> getQLen() {
        return this.T.get(0);
    }

    public List<Double> getUtil() {
        return this.T.get(1);
    }

    public List<Double> getRespT() {
        return this.T.get(2);
    }

    public List<Double> getResidT() {
        return this.T.get(3);
    }

    public List<Double> getTput() {
        return this.T.get(4);
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

    public void print() {
        this.print(this.options);
    }

    public void print(SolverOptions options) {
        if (options.verbose != VerboseLevel.SILENT) {
            System.out.printf(
                    "%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s",
                    "Node", "NodeType", "QLen", "Util", "RespT", "ResidT", "Tput");
            System.out.println(
                    "\n--------------------------------------------------------------------------------------------");
            DecimalFormat nf = new DecimalFormat("#0.#####");
            nf.setMinimumFractionDigits(5);
            for (int i = 0; i < nodeTypes.size(); i++) {
                if (true || getQLen().get(i) > GlobalConstants.Zero ||
                        getUtil().get(i) > GlobalConstants.Zero ||
                        getRespT().get(i) > GlobalConstants.Zero ||
                        getResidT().get(i) > GlobalConstants.Zero ||
                        getTput().get(i) > GlobalConstants.Zero) {
                    System.out.format(
                            "%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\n",
                            nodeNames.get(i),
                            nodeTypes.get(i),
                            nf.format(getQLen().get(i)),
                            nf.format(getUtil().get(i)),
                            nf.format(getRespT().get(i)),
                            nf.format(getResidT().get(i)),
                            nf.format(getTput().get(i)));
                }
            }
            System.out.println(
                    "--------------------------------------------------------------------------------------------");
        }
    }
}
