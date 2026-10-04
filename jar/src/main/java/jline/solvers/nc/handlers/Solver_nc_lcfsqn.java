/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc.handlers;

import java.util.Arrays;
import java.util.List;

import jline.api.pfqn.lcfs.LcfsqnCaResult;
import jline.api.pfqn.lcfs.Pfqn_lcfsqn_ca;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lib.perm.Permanent;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNC;
import jline.util.matrix.Matrix;

/**
 * Specialized NC solver for LCFS + LCFS-PR 2-station networks.
 */
public final class Solver_nc_lcfsqn {
    private Solver_nc_lcfsqn() {}

    public static SolverNC.SolverNCReturn solver_nc_lcfsqn(
            NetworkStruct sn,
            SolverOptions options,
            int lcfsStat,
            int lcfsprStat) {
        long startTime = System.nanoTime();
        int M = sn.nstations;
        int R = sn.nclasses;
        Matrix njobs = sn.njobs;

        Matrix alpha = new Matrix(1, R);
        Matrix beta = new Matrix(1, R);

        Matrix rates = sn.rates;
        for (int r = 0; r < R; r++) {
            if (njobs.get(r) > 0) {
                double mu_lcfs = rates.get(lcfsStat, r);
                double mu_lcfspr = rates.get(lcfsprStat, r);

                if (mu_lcfs <= 0 || Double.isInfinite(mu_lcfs) || Double.isNaN(mu_lcfs)) {
                    throw new RuntimeException("Invalid service rate at LCFS station for class " + r + ".");
                }
                if (mu_lcfspr <= 0 || Double.isInfinite(mu_lcfspr) || Double.isNaN(mu_lcfspr)) {
                    throw new RuntimeException("Invalid service rate at LCFS-PR station for class " + r + ".");
                }

                alpha.set(r, 1.0 / mu_lcfs);
                beta.set(r, 1.0 / mu_lcfspr);
            }
        }

        Matrix N = njobs;
        int K = (int) N.elementSum();

        LcfsqnCaResult caResult = Pfqn_lcfsqn_ca.pfqn_lcfsqn_ca(alpha, beta, N);
        double G = caResult.getG();

        double lG = (G > 0) ? Math.log(G) : Double.NEGATIVE_INFINITY;

        Matrix Q_lcfs = new Matrix(2, R);
        Matrix U_lcfs = new Matrix(2, R);
        Matrix T_lcfs = new Matrix(2, R);

        // The permanent formulas assume one job per class (K single-job
        // classes); classes with multiplicity N(r) > 1 are expanded into N(r)
        // exchangeable single-job copies, with Gexp = G * prod_r N(r)! the
        // distinguishable-jobs normalizing constant, and per-copy metrics are
        // scaled back by N(r).
        Matrix alphaE = new Matrix(1, K);
        Matrix betaE = new Matrix(1, K);
        int[] ecls = new int[R];
        double Gexp = G;
        int pos = 0;
        for (int r = 0; r < R; r++) {
            int nr = (int) N.get(r);
            ecls[r] = pos;
            for (int f = 2; f <= nr; f++) {
                Gexp *= f;
            }
            for (int rep = 0; rep < nr; rep++) {
                alphaE.set(0, pos, alpha.get(r));
                betaE.set(0, pos, beta.get(r));
                pos++;
            }
        }

        for (int r = 0; r < R; r++) {
            if (njobs.get(r) > 0) {
                int e = ecls[r];
                double Tcopy = 0.0;
                double Qcopy = 0.0;
                for (int xt = 1; xt <= K; xt++) {
                    Matrix Tx = makeTx(alphaE, betaE, xt, K, K, e);
                    double permT = new Permanent(Tx, true).value;
                    Tcopy += Math.pow(alphaE.get(e), xt - 1) * permT / Gexp;

                    Matrix Yx = makeYx(alphaE, betaE, xt, K, K, e);
                    double permY = new Permanent(Yx, true).value;
                    Qcopy += permY / Gexp;
                }
                double nr = N.get(r);
                T_lcfs.set(0, r, nr * Tcopy);
                T_lcfs.set(1, r, nr * Tcopy);
                Q_lcfs.set(0, r, nr * Qcopy);
                Q_lcfs.set(1, r, njobs.get(r) - Q_lcfs.get(0, r));
            }
        }

        for (int r = 0; r < R; r++) {
            U_lcfs.set(0, r, T_lcfs.get(0, r) * alpha.get(r));
            U_lcfs.set(1, r, T_lcfs.get(1, r) * beta.get(r));
        }

        Matrix Q = new Matrix(M, R);
        Matrix U = new Matrix(M, R);
        Matrix T = new Matrix(M, R);
        Matrix R_resp = new Matrix(M, R);
        Matrix X = new Matrix(1, R);
        Matrix C = new Matrix(1, R);

        for (int r = 0; r < R; r++) {
            Q.set(lcfsStat, r, Q_lcfs.get(0, r));
            Q.set(lcfsprStat, r, Q_lcfs.get(1, r));
        }

        for (int r = 0; r < R; r++) {
            U.set(lcfsStat, r, U_lcfs.get(0, r));
            U.set(lcfsprStat, r, U_lcfs.get(1, r));
        }

        for (int r = 0; r < R; r++) {
            if (njobs.get(r) > 0) {
                X.set(r, T_lcfs.get(0, r));
                T.set(lcfsStat, r, T_lcfs.get(0, r));
                T.set(lcfsprStat, r, T_lcfs.get(1, r));
            }
        }

        List<Integer> stations = Arrays.asList(Integer.valueOf(lcfsStat), Integer.valueOf(lcfsprStat));
        for (Integer kBox : stations) {
            int k = kBox.intValue();
            for (int r = 0; r < R; r++) {
                if (T.get(k, r) > 0) {
                    R_resp.set(k, r, Q.get(k, r) / T.get(k, r));
                }
            }
        }

        for (int r = 0; r < R; r++) {
            if (njobs.get(r) > 0) {
                C.set(r, R_resp.get(lcfsStat, r) + R_resp.get(lcfsprStat, r));
            }
        }

        double runtime = (System.nanoTime() - startTime) / 1_000_000_000.0;

        return new SolverNC.SolverNCReturn(
                Q, U, R_resp, T, sn.nchains, X, lG, null, 1, runtime, "lcfsqn");
    }

    private static Matrix makeTx(Matrix alpha, Matrix beta, int xt, int K, int R, int r) {
        Matrix Tx = new Matrix(K - 1, K - 1);

        int idx = 0;
        for (int i = 0; i < R; i++) {
            if (i != r) {
                for (int j = 0; j < (xt - 1); j++) {
                    Tx.set(idx, j, Math.pow(alpha.get(i), j + 1));
                }
                for (int j = 0; j < (K - xt); j++) {
                    Tx.set(idx, xt - 1 + j, Math.pow(alpha.get(i), xt + j) * beta.get(i));
                }
                idx++;
            }
        }

        return Tx;
    }

    private static Matrix makeYx(Matrix alpha, Matrix beta, int xt, int K, int R, int r) {
        Matrix Y = new Matrix(R, K);

        for (int i = 0; i < R; i++) {
            for (int j = 0; j < xt; j++) {
                Y.set(i, j, Math.pow(alpha.get(i), j + 1));
            }
            for (int j = 0; j < (K - xt); j++) {
                if (i != r) {
                    Y.set(i, xt + j, Math.pow(alpha.get(i), xt + j) * beta.get(i));
                }
            }
        }

        return Y;
    }

}
