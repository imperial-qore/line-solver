package jline.solvers.ssa;

import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.solvers.SolverResult;
import jline.util.Matrix;

import java.util.Map;

public class SolverSSAResult extends SolverResult {

    public Map<Integer, Matrix> tranSysState;
    public Matrix tranSync;
    public NetworkStruct sn;
    public Map<Station, Matrix> space;

    // empty to allow population of fields from different sources
    public SolverSSAResult() {

    }
    public SolverSSAResult(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix CN, Matrix XN,
                           Map<Integer, Matrix> tranSysState, Matrix tranSync, NetworkStruct sn) {
        this.QN = QN;
        this.UN = UN;
        this.RN = RN;
        this.TN = TN;
        this.CN = CN;
        this.XN = XN;
        this.tranSysState = tranSysState;
        this.tranSync = tranSync;
        this.sn = sn;
    }


}
