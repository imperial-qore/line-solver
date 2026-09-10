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
    /** Declared initial state carried through stochastic complementation, or null. */
    private Matrix pi0;
    /**
     * Derived START filtration, one (state x state) matrix per (station,
     * class): the rate at which a transition starts a class-r service at
     * station i. Kept out of Dfilt, which pairs one-to-one with sn.sync and is
     * summed as D1 -- a START rides on an arc Dfilt already carries.
     */
    private Matrix[][] startFilt;
    /** Derived PREEMPT filtration, laid out like startFilt. */
    private Matrix[][] preemptFilt;


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

    public Matrix[][] getStartFilt() {
        return startFilt;
    }

    public Matrix[][] getPreemptFilt() {
        return preemptFilt;
    }

    public void setAuxFilt(Matrix[][] startFilt, Matrix[][] preemptFilt) {
        this.startFilt = startFilt;
        this.preemptFilt = preemptFilt;
    }

    public Matrix getPi0() {
        return pi0;
    }

    public void setPi0(Matrix pi0) {
        this.pi0 = pi0;
    }
}