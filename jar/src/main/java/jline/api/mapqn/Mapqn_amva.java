package jline.api.mapqn;

import jline.util.matrix.Matrix;

/**
 * Horizontal-cut mean value analysis for a MAP server (SolverMVA method 'amva.mapqn').
 *
 * Closed multiclass network of an exponential infinite-server station (think rate mu_r
 * for class r) and one FCFS single-server station whose class-r service is the MAP
 * (D0_r, D1_r); the MAP of class r moves only while a class-r job is in service and is
 * frozen otherwise, the convention of SolverCTMC. The recursion walks the population
 * lattice n <= N in lexicographic order and solves ONE linear R x R system per point.
 * Its unknowns are the per-phase means Q_r^k = E[n_r 1{k}] over the joint phase
 * k = (k_1..k_R), the busy laws U_r^k = P[serving r, k], the phase law pi_k and the
 * throughputs X_r.
 *
 * Exact relations: the joint phase balance, the class marginals U_r = X_r E[S_r] theta_r
 * and the per-class horizontal cut (generator balance of n_r 1{k}) of Casale-Smirni,
 * DSN 2009. Closures: the product busy law theta_r(k_r) prod_{s != r} phi_s(k_s), phi the
 * post-completion law of a frozen MAP, which solves the phase balance identically; the
 * service-age closure of the cross term E[n_r 1{serving s} 1{k}] (class r accumulates at
 * its throughput over the elapsed class-s service, whose mean given the phase is
 * theta_s (-D0_s)^{-1} / theta_s); Little's law resolved by arrival phase with the exact
 * FCFS response of the queue composition seen at n - e_r (the multiclass arrival
 * theorem). K_r = 1 for every class reproduces multiclass FCFS MVA on class means.
 * Mirrors matlab/src/api/mapqn/mapqn_amva.m and python api/mapqn/amva.py.
 */
public final class Mapqn_amva {
    private Mapqn_amva() {}

    /** X: class throughputs; Qq: mean queue lengths at the MAP station (job in service
     *  included); U: busy probability per class, X E[S]; ES: mean service times;
     *  pi: joint phase law at N (class R fastest). */
    public static final class Result {
        public double[] X;
        public double[] Qq;
        public double[] U;
        public double[] ES;
        public double[] pi;
    }

    public static Result solve(double[] mu, Matrix[] D0s, Matrix[] D1s, int[] N) {
        final int R = N.length;
        final int[] Ks = new int[R];
        final double[][][] D0 = new double[R][][];
        final double[][][] D1 = new double[R][][];
        for (int r = 0; r < R; r++) {
            Ks[r] = D0s[r].getNumRows();
            D0[r] = toArray(D0s[r]);
            D1[r] = toArray(D1s[r]);
        }
        int K = 1;
        for (int r = 0; r < R; r++) K *= Ks[r];
        final int[] stride = new int[R];
        for (int r = 0; r < R; r++) {
            int s = 1;
            for (int t = r + 1; t < R; t++) s *= Ks[t];
            stride[r] = s;                                     // class R-1 is the fastest index
        }
        final int[][] krOf = new int[K][R];
        for (int k = 0; k < K; k++)
            for (int r = 0; r < R; r++) krOf[k][r] = (k / stride[r]) % Ks[r];
        int Ntot = 0;
        for (int r = 0; r < R; r++) Ntot += N[r];

        final double[][] G = new double[R][];
        final double[][] th = new double[R][];
        final double[][] phi = new double[R][];
        final double[] ES = new double[R];
        final double[][][] T = new double[R][][];
        final double[][] age = new double[R][];
        final double[] abar = new double[R];
        final double[][][] Ainv = new double[R][][];
        for (int r = 0; r < R; r++) {
            final int Kr = Ks[r];
            double[][] Gr = new double[Kr][Kr];
            for (int i = 0; i < Kr; i++)
                for (int j = 0; j < Kr; j++) Gr[i][j] = D0[r][i][j] + D1[r][i][j];
            G[r] = flat(Gr, Kr);
            th[r] = stationary(Gr, Kr);
            double rate = 0.0;
            for (int i = 0; i < Kr; i++)
                for (int j = 0; j < Kr; j++) rate += th[r][i] * D1[r][i][j];
            ES[r] = 1.0 / rate;
            phi[r] = new double[Kr];
            for (int j = 0; j < Kr; j++) {
                double a = 0.0;
                for (int i = 0; i < Kr; i++) a += th[r][i] * D1[r][i][j];
                phi[r][j] = a * ES[r];                        // post-completion phase law
            }
            double[][] negD0inv = inverse(negate(D0[r], Kr), Kr);
            double[] s = new double[Kr];                      // mean time to the next completion
            for (int i = 0; i < Kr; i++) for (int j = 0; j < Kr; j++) s[i] += negD0inv[i][j];
            double[][] P = matmul(negD0inv, D1[r], Kr);        // embedded phase transition
            double[][] Tr = new double[Ntot + 2][Kr];          // Tr[j][k] = e_k'(I+P+..+P^(j-1)) s
            double[] acc = new double[Kr];
            double[] v = s.clone();
            for (int j = 1; j <= Ntot + 1; j++) {
                for (int i = 0; i < Kr; i++) acc[i] += v[i];
                Tr[j] = acc.clone();
                double[] nv = new double[Kr];
                for (int i = 0; i < Kr; i++) for (int m = 0; m < Kr; m++) nv[i] += P[i][m] * v[m];
                v = nv;
            }
            T[r] = Tr;
            double[] w = new double[Kr];                       // theta (-D0)^{-1}
            for (int j = 0; j < Kr; j++) for (int i = 0; i < Kr; i++) w[j] += th[r][i] * negD0inv[i][j];
            age[r] = new double[Kr];
            for (int j = 0; j < Kr; j++) { age[r][j] = w[j] / th[r][j]; abar[r] += w[j]; }
            double[][] A = new double[Kr][Kr];
            for (int i = 0; i < Kr; i++) for (int j = 0; j < Kr; j++) A[i][j] = Gr[i][j] - (i == j ? mu[r] : 0.0);
            Ainv[r] = inverse(A, Kr);
        }
        // joint-phase shapes: idle law F and class-r busy law u[r]
        final double[] F = new double[K];
        final double[][] u = new double[R][K];
        for (int k = 0; k < K; k++) {
            F[k] = 1.0;
            for (int r = 0; r < R; r++) F[k] *= phi[r][krOf[k][r]];
            for (int r = 0; r < R; r++) {
                u[r][k] = 1.0;
                for (int s = 0; s < R; s++) u[r][k] *= (s == r) ? th[s][krOf[k][s]] : phi[s][krOf[k][s]];
            }
        }
        // population lattice in lexicographic order: n - e_r always precedes n
        final int[] lstride = new int[R];
        int L = 1;
        for (int r = R - 1; r >= 0; r--) { lstride[r] = L; L *= (N[r] + 1); }
        final double[][][] Qs = new double[L][R][K];
        final double[][] pis = new double[L][K];
        final double[][] Xs = new double[L][R];
        pis[0] = F.clone();
        for (int l = 1; l < L; l++) {
            final int[] n = new int[R];
            for (int r = 0; r < R; r++) n[r] = (l / lstride[r]) % (N[r] + 1);
            final double[][][] b = new double[R][][];
            final double[][][] bN = new double[R][][];
            final double[][] Rk = new double[R][K];
            for (int r = 0; r < R; r++) {
                if (n[r] < 1) continue;
                final int lp = l - lstride[r];
                b[r] = new double[R][K];
                for (int t = 0; t < R; t++)
                    for (int k = 0; k < K; k++) b[r][t][k] = pis[lp][k] > 0 ? Qs[lp][t][k] / pis[lp][k] : 0.0;
                for (int k = 0; k < K; k++) {
                    double acc = 0.0;
                    for (int t = 0; t < R; t++) acc += tAt(T[t], b[r][t][k] + (t == r ? 1.0 : 0.0), krOf[k][t]);
                    Rk[r][k] = acc;
                }
                bN[r] = new double[R][K];
                for (int t = 0; t < R; t++) {
                    if (t == r || n[t] == 0) continue;
                    final double Xt = Xs[lp][t];
                    double Qt = 0.0;
                    for (int k = 0; k < K; k++) Qt += Qs[lp][t][k];
                    if (Xt > 0) {
                        final double W = Math.max(Qt / Xt - abar[r], 0.0);
                        for (int k = 0; k < K; k++)
                            bN[r][t][k] = Xt * Math.min(W + age[r][krOf[k][r]], n[t] / Xt);
                    }
                }
            }
            // the cut, linear in X: Q_r = c0[r] + sum_s X_s c1[r][s]
            final double[][] c0 = new double[R][];
            final double[][][] c1 = new double[R][R][];
            for (int r = 0; r < R; r++) {
                if (n[r] < 1) continue;
                double[] v0 = new double[K];
                for (int k = 0; k < K; k++) v0[k] = -mu[r] * n[r] * F[k];
                c0[r] = applyAxis(v0, Ainv[r], r, krOf, stride, Ks[r], K);
                for (int s = 0; s < R; s++) {
                    double[] term = new double[K];
                    for (int k = 0; k < K; k++) term[k] = -mu[r] * n[r] * ES[s] * (u[s][k] - F[k]);
                    if (s == r) {
                        double[] ud = applyAxis(u[r], D1[r], r, krOf, stride, Ks[r], K);
                        for (int k = 0; k < K; k++) term[k] += ES[r] * ud[k];
                    } else if (n[s] >= 1) {
                        double[] W = new double[K];
                        for (int k = 0; k < K; k++) W[k] = ES[s] * u[s][k] * bN[s][r][k];
                        double[] a1 = applyAxis(W, unflat(G[r], Ks[r]), r, krOf, stride, Ks[r], K);
                        double[] a2 = applyAxis(W, unflat(G[s], Ks[s]), s, krOf, stride, Ks[s], K);
                        for (int k = 0; k < K; k++) term[k] += a1[k] - a2[k];
                    }
                    c1[r][s] = applyAxis(term, Ainv[r], r, krOf, stride, Ks[r], K);
                }
            }
            // Little's law by arrival phase: one R x R solve
            final double[][] M = new double[R][R];
            final double[] v = new double[R];
            for (int r = 0; r < R; r++) M[r][r] = 1.0;
            for (int r = 0; r < R; r++) {
                if (n[r] < 1) continue;
                double acc = 0.0;
                for (int k = 0; k < K; k++) acc += (n[r] * F[k] - c0[r][k]) * Rk[r][k];
                v[r] = n[r] - mu[r] * acc;
                M[r][r] = 1.0 / mu[r];
                for (int s = 0; s < R; s++) {
                    if (n[s] < 1) continue;
                    double a = 0.0;
                    for (int k = 0; k < K; k++) a += (n[r] * ES[s] * (u[s][k] - F[k]) - c1[r][s][k]) * Rk[r][k];
                    M[r][s] += mu[r] * a;
                }
            }
            final double[] X = solveLinear(M, v, R);
            double[] pi = new double[K];
            double psum = 0.0;
            for (int k = 0; k < K; k++) {
                double p = F[k];
                for (int s = 0; s < R; s++) p += X[s] * ES[s] * (u[s][k] - F[k]);
                pi[k] = Math.max(p, 0.0);
                psum += pi[k];
            }
            for (int k = 0; k < K; k++) pi[k] /= psum;
            for (int r = 0; r < R; r++) {
                if (n[r] < 1) continue;
                double[] Ur = new double[K];
                double[] Qr = new double[K];
                double usum = 0.0, wsum = 0.0;
                double[] Wr = new double[K];
                for (int k = 0; k < K; k++) {
                    Ur[k] = X[r] * ES[r] * u[r][k];
                    Qr[k] = c0[r][k];
                    for (int s = 0; s < R; s++) Qr[k] += X[s] * c1[r][s][k];
                    Wr[k] = Math.max(Qr[k] - Ur[k], 0.0);
                    usum += Ur[k];
                    wsum += Wr[k];
                }
                // project onto Q >= U keeping the flow-balance total n_r - X_r/mu_r
                final double tot = Math.max(n[r] - X[r] / mu[r] - usum, 0.0);
                for (int k = 0; k < K; k++) Qs[l][r][k] = Ur[k] + (wsum > 0 ? Wr[k] * tot / wsum : 0.0);
            }
            pis[l] = pi;
            Xs[l] = X;
        }
        Result out = new Result();
        out.X = Xs[L - 1].clone();
        out.Qq = new double[R];
        out.U = new double[R];
        for (int r = 0; r < R; r++) {
            for (int k = 0; k < K; k++) out.Qq[r] += Qs[L - 1][r][k];
            out.U[r] = out.X[r] * ES[r];
        }
        out.ES = ES;
        out.pi = pis[L - 1].clone();
        return out;
    }

    /** out[k] = sum_h V[k with k_r -> h] M[h][k_r]: M applied along the class-r axis. */
    private static double[] applyAxis(double[] V, double[][] M, int r, int[][] krOf, int[] stride, int Kr, int K) {
        double[] out = new double[K];
        for (int k = 0; k < K; k++) {
            final int kr = krOf[k][r];
            final int base = k - kr * stride[r];
            double acc = 0.0;
            for (int h = 0; h < Kr; h++) acc += V[base + h * stride[r]] * M[h][kr];
            out[k] = acc;
        }
        return out;
    }

    private static double tAt(double[][] Tr, double b, int kr) {
        int j0 = (int) Math.floor(b);
        j0 = Math.max(0, Math.min(j0, Tr.length - 2));
        final double f = Math.min(Math.max(b - j0, 0.0), 1.0);
        return (1.0 - f) * Tr[j0][kr] + f * Tr[j0 + 1][kr];
    }

    private static double[] stationary(double[][] G, int K) {
        // th G = 0, th 1 = 1: transpose and replace the last equation by the normalization
        double[][] A = new double[K][K];
        double[] rhs = new double[K];
        for (int i = 0; i < K; i++)
            for (int j = 0; j < K; j++) A[i][j] = G[j][i];
        for (int j = 0; j < K; j++) A[K - 1][j] = 1.0;
        rhs[K - 1] = 1.0;
        double[] th = solveLinear(A, rhs, K);
        double s = 0.0;
        for (int i = 0; i < K; i++) { th[i] = Math.max(th[i], 0.0); s += th[i]; }
        for (int i = 0; i < K; i++) th[i] /= s;
        return th;
    }

    private static double[] solveLinear(double[][] Ain, double[] bin, int n) {
        double[][] A = new double[n][];
        for (int i = 0; i < n; i++) A[i] = Ain[i].clone();
        double[] b = bin.clone();
        for (int c = 0; c < n; c++) {
            int p = c;
            for (int i = c + 1; i < n; i++) if (Math.abs(A[i][c]) > Math.abs(A[p][c])) p = i;
            double[] tr = A[c]; A[c] = A[p]; A[p] = tr;
            double tb = b[c]; b[c] = b[p]; b[p] = tb;
            final double piv = A[c][c];
            for (int i = c + 1; i < n; i++) {
                final double f = A[i][c] / piv;
                if (f == 0.0) continue;
                for (int j = c; j < n; j++) A[i][j] -= f * A[c][j];
                b[i] -= f * b[c];
            }
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double s = b[i];
            for (int j = i + 1; j < n; j++) s -= A[i][j] * x[j];
            x[i] = s / A[i][i];
        }
        return x;
    }

    private static double[][] inverse(double[][] A, int n) {
        double[][] inv = new double[n][n];
        for (int c = 0; c < n; c++) {
            double[] e = new double[n];
            e[c] = 1.0;
            double[] col = solveLinear(A, e, n);
            for (int i = 0; i < n; i++) inv[i][c] = col[i];
        }
        return inv;
    }

    private static double[][] negate(double[][] A, int n) {
        double[][] B = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) B[i][j] = -A[i][j];
        return B;
    }

    private static double[][] matmul(double[][] A, double[][] B, int n) {
        double[][] C = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int m = 0; m < n; m++) s += A[i][m] * B[m][j];
            C[i][j] = s;
        }
        return C;
    }

    private static double[] flat(double[][] A, int n) {
        double[] f = new double[n * n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) f[i * n + j] = A[i][j];
        return f;
    }

    private static double[][] unflat(double[] f, int n) {
        double[][] A = new double[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = f[i * n + j];
        return A;
    }

    private static double[][] toArray(Matrix M) {
        final int n = M.getNumRows();
        final int m = M.getNumCols();
        double[][] A = new double[n][m];
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) A[i][j] = M.get(i, j);
        return A;
    }
}
