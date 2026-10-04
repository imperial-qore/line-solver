/**
 * @file Birman-Kogan asymptotic evaluation of closed networks with many stations
 *
 * Birman and Kogan (Communications in Statistics. Stochastic Models 8(3):543-563,
 * 1992) evaluate the multichain partition function by the saddle point method
 * applied to the Cauchy inversion of its generating function. This file carries
 * the three algorithms of that paper: the saddle point with bottleneck detection
 * (Propositions 1 and 3 with Algorithm 1), the van der Waerden uniform expansion
 * for a single chain, and the load concealment reduction of a multichain network to
 * single chain problems (Algorithm 2).
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.api.pfqn.mva.Pfqn_mva;
import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_bk {
    private Pfqn_bk() {}

    /**
     * Birman-Kogan saddle point normalizing constant with bottleneck detection.
     *
     * Stations that serve a single chain and appear only once (the paper's
     * dedicated single servers) stay outside the exponent as O(1) algebraic
     * factors, so their poles may be crossed by the saddle point; Algorithm 1
     * detects those chains and pins their coordinate on the pole. The remaining
     * stations are the paper's large groups of identical stations.
     *
     * @param L service demand matrix (stations x classes)
     * @param N population vector (1 x classes)
     * @param Z think time vector (1 x classes), may be empty
     * @return G, lG, the saddle point coordinates in X and the utilizations in Q
     */
    public static Ret.pfqnNc pfqn_bk(Matrix L, Matrix N, Matrix Z) {
        Ret.pfqnNc out = new Ret.pfqnNc(1.0, 0.0);
        out.method = "bk";
        if (L == null || L.isEmpty() || N == null || N.isEmpty()) {
            return out;
        }
        int M = L.getNumRows();
        int R = L.getNumCols();
        double[] Nv = new double[R];
        double[] Zv = new double[R];
        double Ntot = 0.0;
        for (int r = 0; r < R; r++) {
            Nv[r] = N.get(0, r);
            Ntot += Nv[r];
            if (Z != null && !Z.isEmpty()) {
                for (int i = 0; i < Z.getNumRows(); i++) Zv[r] += Z.get(i, r);
            }
        }
        out.X = new Matrix(1, R);
        out.Q = new Matrix(M, R);
        if (Ntot <= GlobalConstants.Zero) {
            return out;
        }
        // An empty class contributes a factor of 1 and has no saddle coordinate
        int nEmpty = 0;
        for (int r = 0; r < R; r++) if (Nv[r] <= GlobalConstants.Zero) nEmpty++;
        if (nEmpty > 0 && nEmpty < R) {
            int Rk = R - nEmpty;
            Matrix Lk = new Matrix(M, Rk);
            Matrix Nk = new Matrix(1, Rk);
            Matrix Zk = new Matrix(1, Rk);
            int[] map = new int[Rk];
            int c = 0;
            for (int r = 0; r < R; r++) {
                if (Nv[r] <= GlobalConstants.Zero) continue;
                for (int i = 0; i < M; i++) Lk.set(i, c, L.get(i, r));
                Nk.set(0, c, Nv[r]);
                Zk.set(0, c, Zv[r]);
                map[c] = r;
                c++;
            }
            Ret.pfqnNc red = pfqn_bk(Lk, Nk, Zk);
            out.G = red.G;
            out.lG = red.lG;
            for (int k = 0; k < Rk; k++) {
                out.X.set(0, map[k], red.X.get(0, k));
                for (int i = 0; i < M; i++) out.Q.set(i, map[k], red.Q.get(i, k));
            }
            return out;
        }

        // stations with no demand at all do not enter the generating function
        int Mq = 0;
        for (int i = 0; i < M; i++) {
            double s = 0.0;
            for (int r = 0; r < R; r++) s += L.get(i, r);
            if (s > 0) Mq++;
        }
        double[][] Lq = new double[Mq][R];
        int q = 0;
        for (int i = 0; i < M; i++) {
            double s = 0.0;
            for (int r = 0; r < R; r++) s += L.get(i, r);
            if (s <= 0) continue;
            for (int r = 0; r < R; r++) Lq[q][r] = L.get(i, r);
            q++;
        }

        // Dedicated station of each chain: single chain, no identical twin, no
        // think time, and only when the model holds a group of identical
        // stations, since it is against M_j >> 1 replicas that a lone station is
        // an O(1) factor rather than part of the exponent.
        double[] mu = new double[R];
        int[] poleRow = new int[R];
        for (int r = 0; r < R; r++) {
            mu[r] = Double.POSITIVE_INFINITY;
            poleRow[r] = -1;
        }
        boolean[] isPole = new boolean[Mq];
        int[] mult = multiplicity(Lq, R);
        boolean hasGroup = false;
        for (int i = 0; i < Mq; i++) if (mult[i] > 1) hasGroup = true;
        if (hasGroup) {
            for (int i = 0; i < Mq; i++) {
                if (mult[i] > 1) continue;
                int nz = -1;
                int cnt = 0;
                for (int r = 0; r < R; r++) {
                    if (Lq[i][r] > 0) { nz = r; cnt++; }
                }
                if (cnt != 1) continue;
                if (Zv[nz] > GlobalConstants.FineTol) continue;
                if (1.0 / Lq[i][nz] < mu[nz]) {
                    mu[nz] = 1.0 / Lq[i][nz];
                    poleRow[nz] = i;
                }
            }
            for (int r = 0; r < R; r++) if (poleRow[r] >= 0) isPole[poleRow[r]] = true;
        }
        int Mg = 0;
        for (int i = 0; i < Mq; i++) if (!isPole[i]) Mg++;
        double[][] Lg = new double[Mg][R];
        int g = 0;
        for (int i = 0; i < Mq; i++) {
            if (isPole[i]) continue;
            Lg[g++] = Lq[i];
        }

        boolean[] onBound = new boolean[R];
        double[] z = minimize(Lg, Nv, Zv, mu, onBound);
        int nFree = 0;
        for (int r = 0; r < R; r++) if (!onBound[r]) nFree++;
        for (int r = 0; r < R; r++) {
            out.X.set(0, r, z[r]);
            for (int i = 0; i < M; i++) {
                double u = L.get(i, r) * z[r];
                out.Q.set(i, r, onBound[r] ? FastMath.min(u, 1.0) : u);
            }
        }

        double psi0 = psi(z, Lg, Nv, Zv);
        double lG;
        if (nFree == 0) { // eq. (25): the residues carry everything
            lG = psi0;
        } else {
            double[][] H = hessian(z, Lg, Nv);
            double[][] Haa = new double[nFree][nFree];
            int[] free = new int[nFree];
            int a = 0;
            for (int r = 0; r < R; r++) if (!onBound[r]) free[a++] = r;
            for (int i = 0; i < nFree; i++) {
                for (int j = 0; j < nFree; j++) Haa[i][j] = H[free[i]][free[j]];
            }
            lG = psi0 - 0.5 * nFree * FastMath.log(2 * Math.PI) - 0.5 * logDet(Haa);
            for (int i = 0; i < nFree; i++) {
                int r = free[i];
                lG -= FastMath.log(z[r]);
                if (!Double.isInfinite(mu[r])) lG -= FastMath.log(1.0 - z[r] / mu[r]);
            }
        }
        out.lG = lG;
        out.G = FastMath.exp(lG);
        return out;
    }

    /**
     * Birman-Kogan uniform (van der Waerden) expansion for a single chain.
     *
     * The plain saddle point loses accuracy once the saddle approaches the
     * dominant pole of the integrand, which is the regime where the station
     * holding that pole saturates. The uniform expansion keeps the pole and the
     * saddle in one formula through the complementary error function.
     *
     * @param L service demand vector, single class
     * @param N population
     * @param Z think time
     * @return G and lG
     */
    public static Ret.pfqnNc pfqn_bkue(Matrix L, double N, double Z) {
        Ret.pfqnNc out = new Ret.pfqnNc(1.0, 0.0);
        out.method = "bkue";
        if (N <= GlobalConstants.Zero) {
            return out;
        }
        int n = (L == null || L.isEmpty()) ? 0 : L.getNumRows() * L.getNumCols();
        double[] all = new double[n];
        int cnt = 0;
        for (int i = 0; i < n; i++) {
            double v = L.get(i);
            if (v > 0) all[cnt++] = v;
        }
        if (cnt == 0) {
            out.lG = N * FastMath.log(Z) - Gamma.logGamma(N + 1.0);
            out.G = FastMath.exp(out.lG);
            return out;
        }
        double[] Lv = new double[cnt];
        System.arraycopy(all, 0, Lv, 0, cnt);
        int ipole = 0;
        for (int i = 1; i < cnt; i++) if (Lv[i] > Lv[ipole]) ipole = i;
        double dmax = Lv[ipole];
        double tolL = GlobalConstants.FineTol * FastMath.max(1.0, dmax);
        int ties = 0;
        for (int i = 0; i < cnt; i++) if (FastMath.abs(Lv[i] - dmax) <= tolL) ties++;
        double[] sorted = new double[cnt];
        System.arraycopy(Lv, 0, sorted, 0, cnt);
        java.util.Arrays.sort(sorted);
        boolean hasGroup = false;
        for (int i = 1; i < cnt; i++) if (sorted[i] - sorted[i - 1] <= tolL) hasGroup = true;
        boolean hasPole = ties == 1 && hasGroup;
        double[] D;
        double zp;
        if (hasPole) {
            D = new double[cnt - 1];
            int c = 0;
            for (int i = 0; i < cnt; i++) if (i != ipole) D[c++] = Lv[i];
            zp = 1.0 / dmax;
        } else {
            D = Lv;
            zp = Double.POSITIVE_INFINITY;
        }
        double z0 = saddle1(D, N, Z);
        double h2 = h1d2(z0, D, N);
        double h3 = h1d3(z0, D, N);
        double lG;
        if (Double.isInfinite(zp)) {
            // No pole to keep out of the exponent: the expansion degenerates to
            // the plain saddle point, and the third derivative term goes with the
            // pole it corrects.
            lG = h1(z0, D, N, Z) - FastMath.log(z0) - 0.5 * FastMath.log(2 * Math.PI * h2);
        } else {
            double t2 = (1.0 / z0 + h3 / (6 * h2)) / FastMath.sqrt(2 * Math.PI * h2);
            double b2 = FastMath.max(0.0, h1(zp, D, N, Z) - h1(z0, D, N, Z));
            if (zp >= z0) { // saddle before the pole
                lG = h1(z0, D, N, Z) + FastMath.log(0.5 * erfcx(FastMath.sqrt(b2)) + t2);
            } else { // the pole has been crossed and its residue leads
                lG = h1(zp, D, N, Z)
                        + FastMath.log(1.0 - 0.5 * Erf.erfc(FastMath.sqrt(b2)) + t2 * FastMath.exp(-b2));
            }
        }
        out.lG = lG;
        out.G = FastMath.exp(lG);
        return out;
    }

    /**
     * Birman-Kogan load concealment algorithm (Algorithm 2).
     *
     * Chain l is solved on its own with every station slowed by the residual
     * capacity the other chains leave it, A_i = 1 - sum_{k != l} L(i,k)*X_k.
     * Sweeping the chains in Gauss-Seidel order and iterating to a fixed point is
     * the load concealment algorithm.
     *
     * @param L service demand matrix (stations x classes)
     * @param N population vector (1 x classes)
     * @param Z think time vector (1 x classes), may be empty
     * @param method single chain solver, "mva" (default) or "ue"
     * @param tol convergence tolerance on the throughputs
     * @param maxiter maximum number of sweeps
     * @return throughputs, queue lengths, utilizations and the sweep count
     */
    public static Ret.pfqnBkLc pfqn_bklc(Matrix L, Matrix N, Matrix Z, String method,
                                             double tol, int maxiter) {
        int M = (L == null || L.isEmpty()) ? 0 : L.getNumRows();
        int R = (L == null || L.isEmpty()) ? 0 : L.getNumCols();
        Matrix X = new Matrix(1, Math.max(R, 1));
        Matrix Q = new Matrix(Math.max(M, 1), Math.max(R, 1));
        Matrix U = new Matrix(Math.max(M, 1), Math.max(R, 1));
        if (M == 0 || R == 0) {
            return new Ret.pfqnBkLc(X, Q, U, 0);
        }
        boolean ue = "ue".equalsIgnoreCase(method);
        if (tol <= 0) tol = 1e-10;
        if (maxiter <= 0) maxiter = 1000;
        double[] Nv = new double[R];
        double[] Zv = new double[R];
        for (int r = 0; r < R; r++) {
            Nv[r] = N.get(0, r);
            if (Z != null && !Z.isEmpty()) {
                for (int i = 0; i < Z.getNumRows(); i++) Zv[r] += Z.get(i, r);
            }
        }
        double Ntot = 0.0;
        for (int r = 0; r < R; r++) Ntot += Nv[r];
        if (Ntot <= GlobalConstants.Zero) {
            return new Ret.pfqnBkLc(X, Q, U, 0);
        }
        // Step 1: the saddle point utilizations of Corollary 1 seed the iteration
        double[] Xv = new double[R];
        Ret.pfqnNc seed = pfqn_bk(L, N, Z);
        for (int r = 0; r < R; r++) {
            double v = (seed.X == null || seed.X.isEmpty()) ? 0.0 : seed.X.get(0, r);
            if (!Double.isFinite(v) || v < 0) v = 0.0;
            Xv[r] = v;
        }
        for (int r = 0; r < R; r++) {
            double cap = 0.0;
            double sum = 0.0;
            for (int i = 0; i < M; i++) {
                cap = FastMath.max(cap, L.get(i, r));
                sum += L.get(i, r);
            }
            if (Xv[r] == 0 && Nv[r] > 0) Xv[r] = Nv[r] / (Zv[r] + sum);
            if (cap > 0) Xv[r] = FastMath.min(Xv[r], 1.0 / cap);
        }

        double[][] Qv = new double[M][R];
        int it = 0;
        for (it = 1; it <= maxiter; it++) {
            double[] Xold = new double[R];
            System.arraycopy(Xv, 0, Xold, 0, R);
            for (int l = 0; l < R; l++) {
                if (Nv[l] <= 0) {
                    Xv[l] = 0.0;
                    for (int i = 0; i < M; i++) Qv[i][l] = 0.0;
                    continue;
                }
                // Step 2a: residual capacity left to chain l at every station
                double[] D = new double[M];
                for (int i = 0; i < M; i++) {
                    double busy = 0.0;
                    for (int k = 0; k < R; k++) if (k != l) busy += L.get(i, k) * Xv[k];
                    double A = FastMath.max(1.0 - busy, GlobalConstants.FineTol);
                    D[i] = L.get(i, l) / A;
                }
                // Step 2b: solve the single chain network with the concealed rates
                double[] Qi = new double[M];
                double Xl;
                if (ue) {
                    Matrix Dm = new Matrix(M, 1);
                    for (int i = 0; i < M; i++) Dm.set(i, 0, D[i]);
                    double lgPrev = 0.0;
                    Xl = 0.0;
                    for (int nn = 1; nn <= (int) FastMath.round(Nv[l]); nn++) {
                        double lgn = pfqn_bkue(Dm, nn, Zv[l]).lG;
                        Xl = FastMath.exp(lgPrev - lgn);
                        for (int i = 0; i < M; i++) Qi[i] = D[i] * Xl * (1.0 + Qi[i]);
                        lgPrev = lgn;
                    }
                } else {
                    Matrix Dm = new Matrix(M, 1);
                    for (int i = 0; i < M; i++) Dm.set(i, 0, D[i]);
                    Matrix Nl = new Matrix(1, 1);
                    Nl.set(0, 0, Nv[l]);
                    Matrix Zl = new Matrix(1, 1);
                    Zl.set(0, 0, Zv[l]);
                    Ret.pfqnMVA mva = Pfqn_mva.pfqn_mva(Dm, Nl, Zl);
                    Xl = mva.X.get(0, 0);
                    for (int i = 0; i < M; i++) Qi[i] = mva.Q.get(i, 0);
                }
                Xv[l] = Xl;
                for (int i = 0; i < M; i++) Qv[i][l] = Qi[i];
            }
            double diff = 0.0;
            double xmax = 1.0;
            for (int r = 0; r < R; r++) {
                diff = FastMath.max(diff, FastMath.abs(Xv[r] - Xold[r]));
                xmax = FastMath.max(xmax, FastMath.abs(Xv[r]));
            }
            if (diff <= tol * xmax) break;
        }
        if (it > maxiter) it = maxiter;
        for (int r = 0; r < R; r++) {
            X.set(0, r, Xv[r]);
            for (int i = 0; i < M; i++) {
                Q.set(i, r, Qv[i][r]);
                U.set(i, r, L.get(i, r) * Xv[r]);
            }
        }
        return new Ret.pfqnBkLc(X, Q, U, it);
    }

    // ---- exponent of the integrand, groups only (eq. 11 and 23) ----

    private static double psi(double[] z, double[][] Lg, double[] N, double[] Z) {
        int R = z.length;
        double f = 0.0;
        for (int r = 0; r < R; r++) f += Z[r] * z[r] - N[r] * FastMath.log(z[r]);
        for (int i = 0; i < Lg.length; i++) {
            double u = 0.0;
            for (int r = 0; r < R; r++) u += Lg[i][r] * z[r];
            f -= FastMath.log(1.0 - u);
        }
        return f;
    }

    private static double[] grad(double[] z, double[][] Lg, double[] N, double[] Z) {
        int R = z.length;
        double[] g = new double[R];
        for (int r = 0; r < R; r++) g[r] = Z[r] - N[r] / z[r];
        for (int i = 0; i < Lg.length; i++) {
            double u = 0.0;
            for (int r = 0; r < R; r++) u += Lg[i][r] * z[r];
            double d = 1.0 / (1.0 - u);
            for (int r = 0; r < R; r++) g[r] += d * Lg[i][r];
        }
        return g;
    }

    private static double[][] hessian(double[] z, double[][] Lg, double[] N) {
        int R = z.length;
        double[][] H = new double[R][R];
        for (int r = 0; r < R; r++) H[r][r] = N[r] / (z[r] * z[r]);
        for (int i = 0; i < Lg.length; i++) {
            double u = 0.0;
            for (int r = 0; r < R; r++) u += Lg[i][r] * z[r];
            double d = 1.0 / (1.0 - u);
            double d2 = d * d;
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) H[r][s] += d2 * Lg[i][r] * Lg[i][s];
            }
        }
        return H;
    }

    /** Minimize psi over {z>0, Lg z<1, z<=mu}: Algorithm 1 as an active set method. */
    private static double[] minimize(double[][] Lg, double[] N, double[] Z, double[] mu,
                                     boolean[] onBound) {
        int R = N.length;
        double[] z = init(Lg, N, Z, mu);
        for (int outer = 0; outer <= R; outer++) {
            for (int r = 0; r < R; r++) if (onBound[r]) z[r] = mu[r];
            int nFree = 0;
            for (int r = 0; r < R; r++) if (!onBound[r]) nFree++;
            if (nFree == 0) break;
            int[] free = new int[nFree];
            int a = 0;
            for (int r = 0; r < R; r++) if (!onBound[r]) free[a++] = r;
            double Ntot = 0.0;
            for (int r = 0; r < R; r++) Ntot += N[r];
            for (int it = 0; it < 500; it++) {
                double[] g = grad(z, Lg, N, Z);
                double gnorm = 0.0;
                for (int i = 0; i < nFree; i++) gnorm += g[free[i]] * g[free[i]];
                gnorm = FastMath.sqrt(gnorm);
                if (gnorm <= 1e-12 * FastMath.max(1.0, Ntot)) break;
                double[][] H = hessian(z, Lg, N);
                Matrix Hf = new Matrix(nFree, nFree);
                Matrix rhs = new Matrix(nFree, 1);
                for (int i = 0; i < nFree; i++) {
                    for (int j = 0; j < nFree; j++) Hf.set(i, j, H[free[i]][free[j]]);
                    rhs.set(i, 0, -g[free[i]]);
                }
                Matrix sol = new Matrix(nFree, 1);
                if (!Matrix.solve(Hf, rhs, sol)) break;
                double alpha = 1.0;
                double[] zt = new double[R];
                boolean ok = false;
                while (alpha >= 1e-14) {
                    System.arraycopy(z, 0, zt, 0, R);
                    for (int i = 0; i < nFree; i++) zt[free[i]] = z[free[i]] + alpha * sol.get(i, 0);
                    ok = true;
                    for (int i = 0; i < nFree && ok; i++) {
                        int r = free[i];
                        if (!(zt[r] > 0) || zt[r] > mu[r]) ok = false;
                    }
                    for (int i = 0; i < Lg.length && ok; i++) {
                        double u = 0.0;
                        for (int r = 0; r < R; r++) u += Lg[i][r] * zt[r];
                        if (u >= 1.0) ok = false;
                    }
                    if (ok) break;
                    alpha /= 2;
                }
                if (!ok) break;
                System.arraycopy(zt, 0, z, 0, R);
            }
            // a chain whose descent direction still pushes past its pole is in B
            double[] g = grad(z, Lg, N, Z);
            boolean any = false;
            for (int i = 0; i < nFree; i++) {
                int r = free[i];
                if (z[r] >= mu[r] * (1 - 1e-9) && g[r] < 0) {
                    onBound[r] = true;
                    any = true;
                }
            }
            if (!any) break;
        }
        for (int r = 0; r < R; r++) if (onBound[r]) z[r] = mu[r];
        return z;
    }

    private static double[] init(double[][] Lg, double[] N, double[] Z, double[] mu) {
        int R = N.length;
        double[] z = new double[R];
        for (int r = 0; r < R; r++) {
            double den = Z[r];
            for (int i = 0; i < Lg.length; i++) den += Lg[i][r];
            z[r] = N[r] / FastMath.max(den, GlobalConstants.FineTol);
            z[r] = FastMath.min(z[r], 0.99 * mu[r]);
        }
        for (int it = 0; it < 200 && Lg.length > 0; it++) {
            double umax = 0.0;
            for (int i = 0; i < Lg.length; i++) {
                double u = 0.0;
                for (int r = 0; r < R; r++) u += Lg[i][r] * z[r];
                umax = FastMath.max(umax, u);
            }
            if (umax < 0.9) break;
            for (int r = 0; r < R; r++) z[r] *= 0.7;
        }
        return z;
    }

    /** Number of stations sharing each demand row, up to relative rounding. */
    private static int[] multiplicity(double[][] L, int R) {
        int M = L.length;
        int[] mult = new int[M];
        for (int i = 0; i < M; i++) mult[i] = 1;
        for (int i = 0; i < M; i++) {
            if (mult[i] > 1) continue;
            for (int j = i + 1; j < M; j++) {
                double scale = 1.0;
                double diff = 0.0;
                for (int r = 0; r < R; r++) {
                    scale = FastMath.max(scale, FastMath.max(FastMath.abs(L[i][r]), FastMath.abs(L[j][r])));
                    diff = FastMath.max(diff, FastMath.abs(L[i][r] - L[j][r]));
                }
                if (diff <= GlobalConstants.FineTol * scale) {
                    mult[i]++;
                    mult[j]++;
                }
            }
        }
        return mult;
    }

    /** log|det(A)| from an LU factorization with partial pivoting. */
    private static double logDet(double[][] A) {
        int n = A.length;
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++) System.arraycopy(A[i], 0, a[i], 0, n);
        double ld = 0.0;
        for (int k = 0; k < n; k++) {
            int p = k;
            for (int i = k + 1; i < n; i++) if (FastMath.abs(a[i][k]) > FastMath.abs(a[p][k])) p = i;
            if (p != k) {
                double[] t = a[p];
                a[p] = a[k];
                a[k] = t;
            }
            if (a[k][k] == 0.0) return Double.NEGATIVE_INFINITY;
            ld += FastMath.log(FastMath.abs(a[k][k]));
            for (int i = k + 1; i < n; i++) {
                double f = a[i][k] / a[k][k];
                for (int j = k; j < n; j++) a[i][j] -= f * a[k][j];
            }
        }
        return ld;
    }

    // ---- single chain exponent, excluding the pole station ----

    private static double h1(double z, double[] D, double N, double Z) {
        double f = Z * z - N * FastMath.log(z);
        for (int i = 0; i < D.length; i++) f -= FastMath.log(1.0 - D[i] * z);
        return f;
    }

    private static double h1d1(double z, double[] D, double N, double Z) {
        double g = Z - N / z;
        for (int i = 0; i < D.length; i++) g += D[i] / (1.0 - D[i] * z);
        return g;
    }

    private static double h1d2(double z, double[] D, double N) {
        double h = N / (z * z);
        for (int i = 0; i < D.length; i++) {
            double d = 1.0 - D[i] * z;
            h += D[i] * D[i] / (d * d);
        }
        return h;
    }

    private static double h1d3(double z, double[] D, double N) {
        double h = -2 * N / (z * z * z);
        for (int i = 0; i < D.length; i++) {
            double d = 1.0 - D[i] * z;
            h += 2 * D[i] * D[i] * D[i] / (d * d * d);
        }
        return h;
    }

    private static double saddle1(double[] D, double N, double Z) {
        if (D.length == 0) return N / Z;
        double dmax = 0.0;
        for (int i = 0; i < D.length; i++) dmax = FastMath.max(dmax, D[i]);
        double hi = 1.0 / dmax;
        double z = 0.5 * hi;
        for (int it = 0; it < 200; it++) {
            double g = h1d1(z, D, N, Z);
            if (FastMath.abs(g) <= 1e-14 * FastMath.max(1.0, N)) break;
            double dz = -g / h1d2(z, D, N);
            double alpha = 1.0;
            while (z + alpha * dz <= 0 || z + alpha * dz >= hi) {
                alpha /= 2;
                if (alpha < 1e-14) break;
            }
            if (alpha < 1e-14) break;
            z += alpha * dz;
        }
        return z;
    }

    /**
     * Scaled complementary error function exp(x^2)*erfc(x) for x >= 0. The direct
     * product overflows past x ~ 26, where the asymptotic series is already exact
     * to double precision.
     */
    static double erfcx(double x) {
        if (x < 25.0) {
            return FastMath.exp(x * x) * Erf.erfc(x);
        }
        double y = 1.0 / (2.0 * x * x);
        double term = 1.0;
        double sum = 1.0;
        for (int k = 1; k <= 12; k++) {
            term *= -(2 * k - 1) * y;
            sum += term;
        }
        return sum / (x * FastMath.sqrt(Math.PI));
    }
}
