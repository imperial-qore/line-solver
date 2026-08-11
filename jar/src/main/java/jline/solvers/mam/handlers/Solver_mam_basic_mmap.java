/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import jline.api.sn.SnIsOpenModel;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.util.matrix.Matrix;

/**
 * Top-level dispatcher for the MAM/MMAP fork-join decomposition.
 *
 * <p>Port of matlab/src/solvers/MAM/solver_mam_basic_mmap.m. Open networks call
 * {@link Solver_mam_basic_mmap_inner} directly with arrival rates derived from
 * the source/refstat. Closed networks go through
 * {@link Solver_mam_basic_mmap_closed}, which wraps the inner algorithm in an
 * MNA-style bisection on per-class throughput.</p>
 */
public final class Solver_mam_basic_mmap {
    private Solver_mam_basic_mmap() {}

    public static MAMResult solver_mam_basic_mmap(NetworkStruct sn, SolverOptions options) {
        if (SnIsOpenModel.snIsOpenModel(sn)) {
            int K = sn.nclasses;
            int C = sn.nchains;
            double[] lambda = new double[K];
            for (int c = 0; c < C; c++) {
                Matrix inchain = sn.inchain.get(c);
                if (inchain == null || inchain.length() == 0) {
                    continue;
                }
                int refClass = (int) inchain.get(0);
                int refStation = (int) sn.refstat.get(refClass);
                double chainLambda = 0.0;
                for (int i = 0; i < inchain.length(); i++) {
                    double rate = sn.rates.get(refStation, (int) inchain.get(i));
                    if (Double.isFinite(rate)) {
                        chainLambda += rate;
                    }
                }
                for (int i = 0; i < inchain.length(); i++) {
                    lambda[(int) inchain.get(i)] = chainLambda;
                }
            }
            return Solver_mam_basic_mmap_inner.solver_mam_basic_mmap_inner(sn, options, lambda);
        } else {
            return Solver_mam_basic_mmap_closed.solver_mam_basic_mmap_closed(sn, options);
        }
    }
}
