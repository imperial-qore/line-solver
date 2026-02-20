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
    public final Map<Integer, Matrix> tranSysState;
    public final Matrix tranSync;
    public final NetworkStruct sn;

    public SSAValues(Matrix pi, Matrix SSq, Map<Integer, Matrix> arvRates, Map<Integer, Matrix> depRates,
                     Map<Integer, Matrix> tranSysState, Matrix tranSync, NetworkStruct sn) {
        this.pi = pi;
        this.SSq = SSq;
        this.arvRates = arvRates;
        this.depRates = depRates;
        this.tranSysState = tranSysState;
        this.tranSync = tranSync;
        this.sn = sn;
    }


}
