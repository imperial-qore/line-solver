/**
 * @file Traffic variability equations for the Robust Queueing Network Analyzer
 *
 * Assembles and solves the limiting variability equations and the
 * time-dependent index-of-dispersion-for-counts (IDC) equations of the RQNA of
 * W. Whitt and W. You (2018), "A Robust Queueing Network Analyzer Based on
 * Indices of Dispersion".
 *
 * Models a single-class open network of K single-server FCFS queues with
 * Markovian routing P (P[i,j]=p_{i,j}).
 *
 * @since LINE 3.0
 */
package jline.api.npfqn;

import jline.util.matrix.Matrix;

public final class Npfqn_traffic_idc {

    // ---- context fields (populated by the factory) ----
    public int K;
    public double[] lambda0;
    public double[] mu;
    public double[] cs2;
    public double[] c2a0;
    public Matrix P;
    public Matrix Xi;
    public double[] lambda;
    public double[] rho;
    public Matrix lam_ji;
    public Matrix c2alpha;
    public Matrix[] zetaAll;
    public double[] c2a;
    public Matrix c2aij;
    public double[] c2d;
    public double[] c2x;
    public IdcFunction a0IdcFun;
    public PerQueueIdc sIdcFun;

    private Npfqn_traffic_idc() {}

    /**
     * Traffic variability equations for the RQNA.
     *
     * @param lambda0  external arrival rate into each queue (length K)
     * @param P        KxK routing matrix among queues
     * @param c2a0     asymptotic IDC (SCV) of each external arrival process (length K)
     * @param a0IdcFun external arrival IDC handle, a0IdcFun(t) -&gt; length-K vector
     * @param mu       service rate at each queue (length K)
     * @param cs2      service SCV at each queue (length K)
     * @param sIdcFun  service IDC handle, sIdcFun(t) -&gt; length-K vector
     * @param useAlpha include the alpha_{i,j} correction terms (eq. 34)
     * @param useBeta  include the beta_i correction terms (eqs. 38-39)
     * @return populated context with a callable iaFun(t)
     */
    public static Npfqn_traffic_idc npfqn_traffic_idc(double[] lambda0, Matrix P, double[] c2a0,
            IdcFunction a0IdcFun, double[] mu, double[] cs2, PerQueueIdc sIdcFun,
            boolean useAlpha, boolean useBeta) {
        int K = mu.length;
        Npfqn_traffic_idc ctx = new Npfqn_traffic_idc();
        ctx.K = K;
        ctx.lambda0 = lambda0.clone();
        ctx.mu = mu.clone();
        ctx.cs2 = cs2.clone();
        ctx.c2a0 = c2a0.clone();
        ctx.P = P.copy();
        ctx.a0IdcFun = a0IdcFun;
        ctx.sIdcFun = sIdcFun;

        // ----- traffic rate equations (eq. 20-21) -----
        Matrix eyeK = Matrix.eye(K);
        Matrix Xi = eyeK.sub(P.transpose()).inv();  // (I - P')^{-1}
        ctx.Xi = Xi;
        double[] lambda = new double[K];
        double[] rho = new double[K];
        for (int i = 0; i < K; i++) {
            double s = 0;
            for (int j = 0; j < K; j++) {
                s += Xi.get(i, j) * lambda0[j];
            }
            lambda[i] = s;
            rho[i] = lambda[i] / mu[i];
        }
        ctx.lambda = lambda;
        ctx.rho = rho;
        Matrix lam_ji = new Matrix(K, K);   // lam_ji[j,i] = lambda_j p_{j,i}
        for (int j = 0; j < K; j++) {
            for (int i = 0; i < K; i++) {
                lam_ji.set(j, i, lambda[j] * P.get(j, i));
            }
        }
        ctx.lam_ji = lam_ji;

        // ----- alpha correction: c2alpha_{i,j} = 2 Xi_{i,j} p_{i,j} (1-p_{i,j}) -----
        Matrix c2alpha = new Matrix(K, K);
        if (useAlpha) {
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++) {
                    double p = P.get(i, j);
                    c2alpha.set(i, j, 2.0 * Xi.get(i, j) * p * (1.0 - p));
                }
            }
        }
        ctx.c2alpha = c2alpha;

        // ----- beta correction: zeta_{j,i;k,i} (eq. 39) -----
        Matrix[] Sigma = new Matrix[K];
        for (int l = 0; l < K; l++) {
            Matrix Sl = new Matrix(K, K);
            for (int a = 0; a < K; a++) {
                for (int b = 0; b < K; b++) {
                    Sl.set(a, b, -P.get(l, a) * P.get(l, b) * lambda[l]);
                }
            }
            for (int a = 0; a < K; a++) {
                Sl.set(a, a, P.get(l, a) * (1.0 - P.get(l, a)) * lambda[l]);
            }
            Sigma[l] = Sl;
        }
        Matrix Amat = new Matrix(K, K);
        for (int i = 0; i < K; i++) {
            Amat.set(i, i, c2a0[i] * lambda0[i]);
        }
        for (int l = 0; l < K; l++) {
            Amat = Amat.add(1.0, Sigma[l]);
        }
        Matrix[] zetaAll = new Matrix[K];
        double[] c2beta = new double[K];
        if (useBeta) {
            for (int i = 0; i < K; i++) {
                // nu[l,:] = p_{l,i} * Xi[l,:]
                Matrix nu = new Matrix(K, K);
                for (int l = 0; l < K; l++) {
                    for (int c = 0; c < K; c++) {
                        nu.set(l, c, P.get(l, i) * Xi.get(l, c));
                    }
                }
                Matrix Z = nu.mult(Amat).mult(nu.transpose());
                for (int j = 0; j < K; j++) {
                    for (int k = 0; k < K; k++) {
                        double add = 0;
                        for (int c = 0; c < K; c++) {
                            add += nu.get(k, c) * Sigma[j].get(c, i);
                            add += nu.get(j, c) * Sigma[k].get(c, i);
                        }
                        Z.set(j, k, Z.get(j, k) + add);
                    }
                }
                zetaAll[i] = Z;
                double s = 0;
                for (int j = 0; j < K; j++) {
                    for (int k = j + 1; k < K; k++) {
                        s += Z.get(j, k);
                    }
                }
                if (lambda[i] > 0) {
                    c2beta[i] = (2.0 / lambda[i]) * s;
                }
            }
        } else {
            for (int i = 0; i < K; i++) {
                zetaAll[i] = new Matrix(K, K);
            }
        }
        ctx.zetaAll = zetaAll;

        // ----- limiting variability equations (eq. 44): (E - Minf) c = binf -----
        int Na = K, Naij = K * K, N = Na + Naij + K;
        Matrix Minf = new Matrix(N, N);
        Matrix binf = new Matrix(N, 1);
        for (int i = 0; i < K; i++) {
            if (lambda[i] > 0) {
                for (int j = 0; j < K; j++) {
                    Minf.set(i, Na + j * K + i, lam_ji.get(j, i) / lambda[i]);
                }
                binf.set(i, 0, (lambda0[i] / lambda[i]) * c2a0[i] + c2beta[i]);
            }
            for (int j = 0; j < K; j++) {
                Minf.set(Na + i * K + j, Na + Naij + i, P.get(i, j));
                binf.set(Na + i * K + j, 0, (1.0 - P.get(i, j)) + c2alpha.get(i, j));
            }
            Minf.set(Na + Naij + i, i, 1.0);
        }
        Matrix Acoef = Matrix.eye(N).sub(Minf);
        Matrix csol = new Matrix(N, 1);
        Matrix.solve(Acoef, binf, csol);
        double[] c2a = new double[K];
        Matrix c2aij = new Matrix(K, K);
        double[] c2d = new double[K];
        double[] c2x = new double[K];
        for (int i = 0; i < K; i++) {
            c2a[i] = csol.get(i, 0);
        }
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                c2aij.set(i, j, csol.get(Na + i * K + j, 0));
            }
        }
        for (int i = 0; i < K; i++) {
            c2d[i] = csol.get(Na + Naij + i, 0);
            c2x[i] = c2a[i] + cs2[i];
        }
        ctx.c2a = c2a;
        ctx.c2aij = c2aij;
        ctx.c2d = c2d;
        ctx.c2x = c2x;
        return ctx;
    }

    /**
     * Solve the time-dependent IDC equations (eq. 43) at a single time t.
     *
     * @param t time argument
     * @return length-K vector of total arrival IDCs I_{a,i}(t)
     */
    public double[] iaFun(double t) {
        int K = this.K;
        double[] h = new double[K];
        for (int i = 0; i < K; i++) {
            h[i] = rho[i] * rho[i];   // tuning function h(rho)=rho^2
        }
        // departure weights w_i(t) = w*((1-rho_i)^2 lambda_i t /(h_i c2x_i))
        double[] warg = new double[K];
        for (int i = 0; i < K; i++) {
            if (h[i] > 0 && c2x[i] > 0) {
                warg[i] = (1.0 - rho[i]) * (1.0 - rho[i]) * lambda[i] * t / (h[i] * c2x[i]);
            } else {
                warg[i] = Double.POSITIVE_INFINITY;
            }
        }
        double[] w = Npfqn_rqna_weight.npfqn_rqna_weight(warg);

        double[] Ia0 = a0IdcFun.eval(t);
        double[] rhoT = new double[K];
        for (int i = 0; i < K; i++) {
            rhoT[i] = rho[i] * t;
        }
        double[] Is = sIdcFunEval(rhoT);   // service IDC at scaled time rho*t

        // see _kb/03-api-layer.md for rationale
        double[] beta_t = new double[K];
        for (int i = 0; i < K; i++) {
            Matrix Z = zetaAll[i];
            double[] wj = new double[K];
            for (int j = 0; j < K; j++) {
                if (h[j] > 0 && c2x[j] > 0 && P.get(j, i) > 0) {
                    double aj = (1.0 - rho[j]) * (1.0 - rho[j]) * P.get(j, i) * lambda[j] * t / (h[j] * c2x[j]);
                    wj[j] = Npfqn_rqna_weight.npfqn_rqna_weight(aj);
                }
            }
            double s = 0;
            for (int j = 0; j < K; j++) {
                for (int k = 0; k < K; k++) {
                    if (j != k) {
                        s += Z.get(j, k) * wj[j];
                    }
                }
            }
            if (lambda[i] > 0) {
                beta_t[i] = s / lambda[i];
            }
        }

        // assemble (E - M(t)) I = b(t)
        int Na = K, Naij = K * K, N = Na + Naij + K;
        Matrix M = new Matrix(N, N);
        Matrix b = new Matrix(N, 1);
        for (int i = 0; i < K; i++) {
            if (lambda[i] > 0) {
                for (int j = 0; j < K; j++) {
                    M.set(i, Na + j * K + i, lam_ji.get(j, i) / lambda[i]);
                }
                b.set(i, 0, (lambda0[i] / lambda[i]) * Ia0[i] + beta_t[i]);
            }
            for (int j = 0; j < K; j++) {
                M.set(Na + i * K + j, Na + Naij + i, P.get(i, j));
                b.set(Na + i * K + j, 0, (1.0 - P.get(i, j)) + c2alpha.get(i, j) * w[i]);
            }
            M.set(Na + Naij + i, i, w[i]);
            b.set(Na + Naij + i, 0, (1.0 - w[i]) * Is[i]);
        }
        Matrix Acoef = Matrix.eye(N).sub(M);
        Matrix sol = new Matrix(N, 1);
        Matrix.solve(Acoef, b, sol);
        double[] Ia = new double[K];
        for (int i = 0; i < K; i++) {
            Ia[i] = sol.get(i, 0);
        }
        return Ia;
    }

    // The service IDC handle is evaluated at the per-queue scaled time rho_i*t.
    private double[] sIdcFunEval(double[] scaledTimes) {
        return sIdcFun.evalVec(scaledTimes);
    }

    /**
     * Service-IDC handle whose per-queue evaluation time differs by queue. The
     * RQNA evaluates the service IDC of queue i at time rho_i*t.
     */
    public interface PerQueueIdc {
        double[] evalVec(double[] perQueueTimes);
    }
}
