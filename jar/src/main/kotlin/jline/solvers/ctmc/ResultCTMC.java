package jline.solvers.ctmc;

import jline.lang.NetworkStruct;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public class ResultCTMC extends SolverResult {

    private final Matrix Q;
    private final Matrix stateSpaceAggr;
    private final MatrixCell Dfilt;
    private final double[][][] arvRates;
    private final double[][][] depRates;
    private final NetworkStruct sn;
    protected Matrix stateSpace;


    public ResultCTMC(
            Matrix q,
            Matrix stateSpace,
            Matrix stateSpaceAggr,
            MatrixCell dfilt,
            double[][][] arvRates,
            double[][][] depRates,
            NetworkStruct sn) {
        this.Q = q;
        this.stateSpace = stateSpace;
        this.stateSpaceAggr = stateSpaceAggr;
        this.Dfilt = dfilt;
        this.arvRates = arvRates;
        this.depRates = depRates;
        this.sn = sn;
    }

    public double[][][] getArvRates() {
        return arvRates;
    }

    public double[][][] getDepRates() {
        return depRates;
    }

    public MatrixCell getDfilt() {
        return Dfilt;
    }

    public Matrix getQ() {
        return Q;
    }

    public NetworkStruct getSn() {
        return sn;
    }

    public Matrix getStateSpace() {
        return stateSpace;
    }

    public Matrix getStateSpaceAggr() {
        return stateSpaceAggr;
    }
}