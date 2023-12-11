package jline.solvers;

import jline.lang.constant.GlobalConstants;
import jline.lang.constant.VerboseLevel;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class NetworkAvgSysTable {

    SolverOptions options;
    ArrayList<List<Double>> T;
    List<String> chainNames;
    List<String> inChainNames;

    public NetworkAvgSysTable(List<Double> SysRespTval, List<Double> SysTputval, SolverOptions options) {
        this.T = new ArrayList<>(Arrays.asList(SysRespTval,SysTputval));
        this.options = options;
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
    public void setOptions(SolverOptions options) {
        this.options = options;
    }
    public List<Double> get(int col) {
        return this.T.get(col);
    }

    public List<Double> getSysRespT() {
        return this.T.get(0);
    }
    public List<Double> getSysTput() {
        return this.T.get(1);
    }

    public void printTable() {
        this.print();
    }
    public void printTable(SolverOptions options){
        this.print(options);
    }
    public void print() {
        this.print(this.options);
    }
    public void print(SolverOptions options){
        if (options.verbose != VerboseLevel.SILENT) {
            System.out.printf(
                    "%-12s\t %-12s\t %-10s\t %-10s",
                    "Chain", "JobClasses", "SysRespT", "SysTput");
            System.out.println(
                    "\n-----------------------------------------------------------------------------------------------------------------------------------");
            NumberFormat nf = NumberFormat.getNumberInstance();
            //nf.setMinimumFractionDigits(5);
            nf.setMinimumFractionDigits(16);
            for (int i = 0; i < T.get(0).size(); i++) {
                if (getSysRespT().get(i) > GlobalConstants.Zero ||
                        getSysTput().get(i) > GlobalConstants.Zero) {
                    System.out.format(
                            "%-12s\t %-12s\t %-10s\t %-10s\n",
                            chainNames.get(i),
                            inChainNames.get(i),
                            nf.format(getSysRespT().get(i)),
                            nf.format(getSysTput().get(i)));
                }
            }
            System.out.println(
                    "-----------------------------------------------------------------------------------------------------------------------------------");
        }
    }
}
