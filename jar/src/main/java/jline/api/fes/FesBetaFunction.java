/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.fes;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import jline.api.pfqn.ld.Ljd;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

/**
 * Wraps a flow-equivalent-server (FES) throughput table as a per-class
 * class-dependence function beta_{i,r}(n).
 *
 * <p>The tables are the linearized vectors of class-r throughputs X_r(n) of the
 * aggregated subnetwork, indexed by {@link Ljd#ljd_linearize} over the
 * population clamped to the cutoffs the table was tabulated on. The function
 * takes the per-class population vector n at the station and returns the
 * length-K row vector [X_1(n), ..., X_K(n)], i.e. Sauer's chain-dependent
 * service rates mu_{r,i}(n) (Sauer 1983, "Computational Algorithms for
 * State-Dependent Queueing Networks", eq. (40)). Because the population is
 * clamped, the rate saturates beyond the tabulated range exactly as the
 * underlying table intends.</p>
 *
 * <p>This is the single class-dependence mechanism used across the solvers: the
 * exact convolution (Pfqn_conv) reads mu_{r,i}(n) from it, and AMVA-QD reads the
 * same function through Pfqn_cdfun.</p>
 *
 * <p>Mirrors the MATLAB helper fes_beta_handle.m.</p>
 */
public class FesBetaFunction implements SerializableFunction<Matrix, Matrix>, Serializable {

    private static final long serialVersionUID = 1L;

    private final List<Matrix> scalingTables;
    private final Matrix cutoffs;

    /**
     * Creates a per-class FES rate function.
     *
     * @param scalingTables linearized per-class throughput tables, one per class
     * @param cutoffs       per-class population cutoffs the tables were built on
     */
    public FesBetaFunction(List<Matrix> scalingTables, Matrix cutoffs) {
        this.scalingTables = new ArrayList<Matrix>(scalingTables);
        this.cutoffs = cutoffs;
    }

    @Override
    public Matrix apply(Matrix n) {
        int K = scalingTables.size();
        Matrix v = new Matrix(1, K);
        v.fill(1.0);

        int C = cutoffs.length();
        Matrix nClamped = new Matrix(1, C);
        for (int k = 0; k < C; k++) {
            double nk = (k < n.length()) ? Math.round(n.get(k)) : 0.0;
            nClamped.set(0, k, Math.max(0.0, Math.min(nk, cutoffs.get(k))));
        }
        int idx = Ljd.ljd_linearize(nClamped, cutoffs);

        double tot = 0.0;
        for (int k = 0; k < C; k++) {
            tot += (k < n.length()) ? Math.round(n.get(k)) : 0.0;
        }

        for (int r = 0; r < K; r++) {
            Matrix tbl = scalingTables.get(r);
            double nr = (r < n.length()) ? Math.round(n.get(r)) : 0.0;
            if (tbl != null && idx >= 0 && idx < tbl.length() && nr > 0) {
                // see _kb/03-api-layer.md for rationale
                v.set(0, r, tbl.get(idx) * tot / nr);
            }
        }
        return v;
    }
}
