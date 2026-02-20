/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */


package jline.solvers;

// Class for handles for the transient performance metrics
public class SolverTranHandles {

    public AvgHandle Qt;
    public AvgHandle Ut;
    public AvgHandle Tt;

    public SolverTranHandles(AvgHandle Qt, AvgHandle Ut, AvgHandle Tt) {
        this.Qt = Qt;
        this.Ut = Ut;
        this.Tt = Tt;
    }

    public AvgHandle getTranQLenHandles() {
        return Qt;
    }

    public AvgHandle getTranTputHandles() {
        return Tt;
    }

    public AvgHandle getTranUtilHandles() {
        return Ut;
    }

}
