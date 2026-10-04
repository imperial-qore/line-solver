/**
 * @file Approximate MVA of closed networks with shortest-job-next stations
 *
 * Fixed-point (Bard-Schweitzer) counterpart of Pfqn_mvasjn for closed networks with
 * non-preemptive shortest-job-next (SJN/SJF) stations. Ported at parity from MATLAB
 * pfqn_amvasjn.m.
 *
 * Reference: K. Kant, "MVA approximations for SJN scheduling", Performance Evaluation 15(1):41-61,
 * 1992.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.io.InputOutput;
import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_amvasjn {
    private Pfqn_amvasjn() {}

    /**
     * Mean value analysis with shortest-job-next stations, through a Schweitzer fixed point.
     *
     * <p>{@link Pfqn_mvasjn} carries the conditional waiting time profile W(x,n) over the whole
     * population lattice, which costs prod(N+1) steps and rules the method out for large
     * populations. The closure used here rests on the observation that lam_k W_k(x,n) f_k(x) dx is
     * the mean number of queued class-k customers whose service requirement lies in (x, x+dx),
     * that is, the queue length resolved by job size. Schweitzer's assumption is applied to that
     * density rather than to its integral: removing one customer of class r scales the class-r
     * size-resolved queue length by (N_r-1)/N_r and leaves the other classes unchanged.
     * Integrating over x recovers the usual Schweitzer rule for the aggregate queue lengths, so
     * the closure is the exact analogue of the one applied at the ordinary stations.</p>
     *
     * <p>Cost per iteration is O(M R ns) against the prod(N+1) M R ns of the exact recursion. What
     * is given up is the population dependence of the SHAPE of W(x): the closure lets its level
     * scale but keeps its shape fixed, whereas the true profile stiffens with the load because the
     * denominator 1 - sum_k lam_k theta_k(x) sharpens. The error therefore concentrates at high
     * utilization, where the SJN approximation is already at its weakest.</p>
     *
     * <p>The iteration is started from the product-form Schweitzer solution, not from a light-load
     * guess: the latter puts the deflated utilization above one, where the response time equation
     * has no solution.</p>
     *
     * @param L service demand matrix (M x R) of the queueing stations
     * @param N population vector (1 x R)
     * @param Z think time vector (1 x R), may be null
     * @param scv squared coefficients of variation of the service times (M x R), may be null
     * @param sjnset zero-based indices of the stations scheduling by SJN, may be null
     * @param V visit ratios (M x R), so that the per-visit service time is L./V, may be null
     * @param options grid, priority, tolerance and cap options, may be null
     * @return mean performance measures and the converged conditional waiting time profiles
     */
    public static Pfqn_mvasjn.Result pfqn_amvasjn(Matrix L, Matrix N, Matrix Z, Matrix scv,
            int[] sjnset, Matrix V, SjnOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (options == null) {
            options = new SjnOptions();
        }
        options.validate(R);
        double[] Nv = SjnArgs.vector(N, R, "population vector");
        for (int r = 0; r < R; r++) {
            Nv[r] = Math.round(Nv[r]);
            if (Nv[r] < 0) {
                throw new IllegalArgumentException("negative class populations");
            }
        }
        double[] Zv = Z == null || Z.isEmpty() ? new double[R] : SjnArgs.vector(Z, R, "think times");
        double[][] Lm = SjnArgs.matrix(L, M, R);
        double[][] Vm = V == null || V.isEmpty() ? SjnArgs.ones(M, R) : SjnArgs.matrix(V, M, R);
        double[][] scvm = scv == null || scv.isEmpty() ? SjnArgs.ones(M, R)
                : SjnArgs.matrix(scv, M, R);
        double[][] S = SjnArgs.perVisit(Lm, Vm);
        int[] sjn = SjnArgs.stationSet(sjnset, M);
        boolean useprio = options.prio != null;

        int nsjn = sjn.length;
        int ngrid = options.ns + 1;
        SjnSupport.Grid[] G = new SjnSupport.Grid[nsjn];
        for (int q = 0; q < nsjn; q++) {
            G[q] = SjnSupport.setup(S[sjn[q]], scvm[sjn[q]], options.ns, options.Lfactor);
        }

        // start from the product-form Schweitzer solution: a light-load guess would put the
        // deflated utilization above one and the SJN denominator has no solution there
        Matrix Zmat = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            Zmat.set(0, r, Zv[r]);
        }
        Matrix Nmat = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            Nmat.set(0, r, Nv[r]);
        }
        Ret.pfqnAMVA bs = jline.api.pfqn.mva.Pfqn_bs.pfqn_bs(L, Nmat, Zmat, options.tol,
                options.iterMax);
        double[] X = new double[R];
        double[][] Q = new double[M][R];
        double[][] U = new double[M][R];
        double[][] C = new double[M][R];
        for (int r = 0; r < R; r++) {
            X[r] = bs.X.get(0, r);
            for (int m = 0; m < M; m++) {
                Q[m][r] = bs.Q.get(m, r);
                U[m][r] = bs.U.get(m, r);
                C[m][r] = bs.R.get(m, r);
            }
        }
        double[][][] W = new double[nsjn][ngrid][R];
        double[][][] P = new double[nsjn][ngrid][R];
        double[][] Iinf = new double[nsjn][R];
        double[][] T = new double[nsjn][R * 3];

        int it = 0;
        boolean capped = false;
        boolean converged = false;
        double delta = Double.POSITIVE_INFINITY;
        while (!converged && it < options.iterMax) {
            it++;
            double[][] Cit = new double[M][R];
            double[][][] Wit = SjnArgs.copy3(W);
            double[][][] Pit = SjnArgs.copy3(P);
            double[][] Iit = SjnArgs.copy2(Iinf);
            double[][] Tit = SjnArgs.copy2(T);
            for (int r = 0; r < R; r++) {
                if (Nv[r] == 0) {
                    continue;
                }
                double[] beta = SjnArgs.onesVector(R);
                beta[r] = (Nv[r] - 1) / Nv[r];
                for (int m = 0; m < M; m++) {
                    int q = SjnArgs.indexOf(sjn, m);
                    if (q < 0) {
                        double qsum = 0;
                        for (int k = 0; k < R; k++) {
                            qsum += beta[k] * Q[m][k];
                        }
                        Cit[m][r] = Lm[m][r] * (1 + qsum);
                        continue;
                    }
                    double[] lam = new double[R];
                    for (int k = 0; k < R; k++) {
                        lam[k] = X[k] * Vm[m][k];
                    }
                    SjnSupport.State st = new SjnSupport.State(lam, U[m], Q[m], W[q], P[q],
                            Iinf[q]);
                    SjnSupport.StationResult sr = SjnSupport.station(m, r, G[q], S[m], scvm[m],
                            Vm[m], st, beta, useprio, options.prio);
                    Cit[m][r] = sr.C;
                    for (int i = 0; i < ngrid; i++) {
                        Wit[q][i][r] = sr.W[i];
                        Pit[q][i][r] = sr.phi[i];
                    }
                    Iit[q][r] = sr.phiinf;
                    Tit[q][3 * r] = sr.tail[0];
                    Tit[q][3 * r + 1] = sr.tail[1];
                    Tit[q][3 * r + 2] = sr.tail[2];
                }
            }
            SjnSupport.CapResult cr = SjnSupport.cap(Cit, Lm, Nv, Zv, sjn, options.umax);
            if (cr.bound) {
                capped = true;
                for (int q = 0; q < nsjn; q++) {
                    for (int i = 0; i < ngrid; i++) {
                        for (int r = 0; r < R; r++) {
                            Wit[q][i][r] *= cr.kappa[q];
                            Pit[q][i][r] *= cr.kappa[q];
                        }
                    }
                    for (int r = 0; r < R; r++) {
                        Iit[q][r] *= cr.kappa[q];
                        Tit[q][3 * r] *= cr.kappa[q];
                        Tit[q][3 * r + 1] *= cr.kappa[q];
                    }
                }
            }
            delta = 0;
            for (int r = 0; r < R; r++) {
                for (int m = 0; m < M; m++) {
                    double qnew = cr.X[r] * cr.C[m][r];
                    delta = Math.max(delta, Math.abs(qnew - Q[m][r]));
                }
            }
            for (int q = 0; q < nsjn; q++) {
                for (int i = 0; i < ngrid; i++) {
                    for (int r = 0; r < R; r++) {
                        delta = Math.max(delta, Math.abs(Wit[q][i][r] - W[q][i][r]));
                    }
                }
            }
            for (int r = 0; r < R; r++) {
                X[r] = cr.X[r];
                for (int m = 0; m < M; m++) {
                    C[m][r] = cr.C[m][r];
                    Q[m][r] = cr.X[r] * cr.C[m][r];
                    U[m][r] = cr.X[r] * Lm[m][r];
                }
            }
            W = Wit;
            P = Pit;
            Iinf = Iit;
            T = Tit;
            converged = delta < options.tol;
        }
        if (!converged) {
            InputOutput.line_warning(Pfqn_amvasjn.class.getName(), "the SJN fixed point did not"
                    + " converge in " + options.iterMax + " iterations, residual " + delta);
        }
        if (capped) {
            SjnArgs.warnCapped(options.umax);
        }

        Pfqn_mvasjn.Profile[] profiles = new Pfqn_mvasjn.Profile[nsjn];
        for (int q = 0; q < nsjn; q++) {
            double[][] tail = new double[R][3];
            for (int r = 0; r < R; r++) {
                tail[r][0] = T[q][3 * r];
                tail[r][1] = T[q][3 * r + 1];
                tail[r][2] = T[q][3 * r + 2];
            }
            profiles[q] = new Pfqn_mvasjn.Profile(sjn[q], G[q].x, W[q], tail);
        }
        return new Pfqn_mvasjn.Result(SjnArgs.rowMatrix(X), SjnArgs.toMatrix(Q),
                SjnArgs.toMatrix(U), SjnArgs.toMatrix(C), profiles, it);
    }
}
