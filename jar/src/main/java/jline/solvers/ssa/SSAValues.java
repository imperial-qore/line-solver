/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ssa;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

import java.util.Map;

public class SSAValues {

    public final Matrix pi;
    public final Matrix SSq;
    public final Map<Integer, Matrix> arvRates;
    public final Map<Integer, Matrix> depRates;
    /** Rate of the cache merge transitions, i.e. of the delayed hits. */
    public final Map<Integer, Matrix> dlyRates;
    public final Map<Integer, Matrix> tranSysState;
    public final Matrix tranSync;
    public final NetworkStruct sn;
    /**
     * Derived START rates per unique state, laid out like arvRates/depRates:
     * how fast the transitions enabled in that state start a class-r service
     * at a stateful node. Null when the engine did not measure them.
     */
    public Map<Integer, Matrix> startRates;
    /** Derived PREEMPT rates, laid out like startRates. */
    public Map<Integer, Matrix> preemptRates;

    public SSAValues(Matrix pi, Matrix SSq, Map<Integer, Matrix> arvRates, Map<Integer, Matrix> depRates,
                     Map<Integer, Matrix> tranSysState, Matrix tranSync, NetworkStruct sn) {
        this(pi, SSq, arvRates, depRates, null, tranSysState, tranSync, sn);
    }

    public SSAValues(Matrix pi, Matrix SSq, Map<Integer, Matrix> arvRates, Map<Integer, Matrix> depRates,
                     Map<Integer, Matrix> dlyRates,
                     Map<Integer, Matrix> tranSysState, Matrix tranSync, NetworkStruct sn) {
        this.pi = pi;
        this.SSq = SSq;
        this.arvRates = arvRates;
        this.depRates = depRates;
        this.dlyRates = dlyRates;
        this.tranSysState = tranSysState;
        this.tranSync = tranSync;
        this.sn = sn;
    }

    /**
     * Attach the derived START/PREEMPT rates measured over the sampled path.
     */
    public void setTagRates(Map<Integer, Matrix> startRates, Map<Integer, Matrix> preemptRates) {
        this.startRates = startRates;
        this.preemptRates = preemptRates;
    }


}
