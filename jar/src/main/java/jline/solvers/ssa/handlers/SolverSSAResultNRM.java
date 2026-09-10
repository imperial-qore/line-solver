package jline.solvers.ssa.handlers;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

/** Result returned by solver_ssa_nrm. */
public final class SolverSSAResultNRM {
    public final Matrix QN;
    public final Matrix UN;
    public final Matrix RN;
    public final Matrix TN;
    public final Matrix CN;
    public final Matrix XN;
    public final NetworkStruct sn;
    /**
     * Derived START rate [stations x classes]: counts of service starts over
     * the simulated time. The NRM fires one reaction at a time and knows
     * exactly which job takes a server, so these are exact event counts.
     */
    public Matrix StartN;
    /** Derived PREEMPT rate [stations x classes]. */
    public Matrix PreemptN;

    public SolverSSAResultNRM(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix CN, Matrix XN, NetworkStruct sn) {
        this.QN = QN;
        this.UN = UN;
        this.RN = RN;
        this.TN = TN;
        this.CN = CN;
        this.XN = XN;
        this.sn = sn;
    }

    public Matrix getQN() { return QN; }
    public Matrix getUN() { return UN; }
    public Matrix getRN() { return RN; }
    public Matrix getTN() { return TN; }
    public Matrix getCN() { return CN; }
    public Matrix getXN() { return XN; }
    public NetworkStruct getSn() { return sn; }

    public void setTagRates(Matrix StartN, Matrix PreemptN) {
        this.StartN = StartN;
        this.PreemptN = PreemptN;
    }
}
