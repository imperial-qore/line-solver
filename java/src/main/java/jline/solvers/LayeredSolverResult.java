// Copyright (c) 2012-2024, Imperial College London
// All rights reserved.

package jline.solvers;

import jline.util.Matrix;

public class LayeredSolverResult extends SolverResult{

    public Matrix PN; // mean processor utilization  (nidx x 1)
    public Matrix SN; // mean service time  (nidx x 1)

    public LayeredSolverResult deepCopy() {

        LayeredSolverResult copy = new LayeredSolverResult();

        copy.method = this.method;

        copy.QN = this.QN.clone();
        copy.UN = this.UN.clone();
        copy.RN = this.RN.clone();
        copy.WN = this.WN.clone();
        copy.TN = this.TN.clone();
        copy.AN = this.AN.clone();
        copy.PN = this.PN.clone();
        copy.SN = this.SN.clone();
        copy.runtime = this.runtime;

        return copy;
    }
}
