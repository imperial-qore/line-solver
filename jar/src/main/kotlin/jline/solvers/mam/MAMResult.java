/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import java.util.List;

public class MAMResult extends SolverResult {
    public int iter;
    public double lG;
    public Matrix actionRates;
    public List<Matrix> equilibrium;
    public List<Matrix> generators;
}
