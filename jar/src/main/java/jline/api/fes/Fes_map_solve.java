/**
 * @file Closed model left by a MAP flow-equivalent server
 *
 * @since LINE 3.0
 */
package jline.api.fes;

import java.util.List;

import jline.api.mam.Map_moment;
import jline.api.mc.Ctmc_solve;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Solves the reduced model made of a delay and a load-dependent MAP flow-equivalent server.
 *
 * Closes the aggregation of Section 5.2.1 of Casale, Mi, Cherkasova and Smirni, IEEE Trans.
 * Soft. Eng. 37(5), 2011. Once a subnetwork has been replaced by the load-dependent MAP of
 * fes_map_aggregate, the model left is a delay holding the think times and one station,
 * which is a finite level-dependent quasi birth-death process: level k is the number of
 * jobs held by the flow-equivalent server and N-k jobs are thinking. The chain is the same
 * block bidiagonal pair used to measure the inter-departure times, now read as a generator
 * rather than as a MAP, so the delay is a station whose process is scaled by the number of
 * jobs it holds and the marked transitions are the arrivals into the flow-equivalent
 * server.
 *
 * The think time may itself be a MAP, which is how Section 5.3.1 models a bounded flash
 * crowd: burstiness in the stream of requests is carried by (Z0,Z1) and the rates scale
 * with the population at the delay.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Fes_map_solve {
    private Fes_map_solve() {}

    /**
     * @param fes      flow-equivalent server, one MAP per level, from fes_map_aggregate
     * @param thinkMap think time process (Z0,Z1); exponential for a plain think time
     * @param n        number of jobs in the closed model
     * @return system throughput, response time, aggregate queue length and level law
     */
    public static FesMapSolveResult fes_map_solve(List<MatrixCell> fes, MatrixCell thinkMap, int n) {
        if (n < 1) {
            throw new IllegalArgumentException("The population N must be at least 1.");
        }
        List<MatrixCell> fesLev = Fes_map_levels.fes_map_levels(fes, n);
        List<MatrixCell> think = Fes_map_levels.fes_map_levels(thinkMap, n, Double.POSITIVE_INFINITY);

        MatrixCell T = Fes_map_interdeparture.fes_map_interdeparture(think, fesLev, n);
        Matrix T0 = T.get(0);
        Matrix T1 = T.get(1);
        Matrix Q = T0.add(1.0, T1);
        int dim = Q.getNumRows();

        Matrix phi = Ctmc_solve.ctmc_solve(Q);
        Matrix e = new Matrix(dim, 1);
        for (int i = 0; i < dim; i++) {
            e.set(i, 0, 1.0);
        }
        double X = phi.mult(T1).mult(e).get(0, 0);

        int blk = dim / (n + 1);
        double[] pk = new double[n + 1];
        double QN = 0;
        for (int k = 0; k <= n; k++) {
            double s = 0;
            for (int j = 0; j < blk; j++) {
                s += phi.get(0, k * blk + j);
            }
            pk[k] = s;
            QN += k * s;
        }

        double RN = n / X - Map_moment.map_moment(thinkMap, 1);
        return new FesMapSolveResult(X, RN, QN, pk);
    }
}
