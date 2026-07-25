/**
 * @file PANACEA asymptotic expansion for load-dependent closed networks.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

/**
 * Mitra-McKenna (JACM 33(3):568-592, 1986) load-dependent PANACEA: the
 * expansion coefficients A_n are linear combinations of partition functions of
 * a pseudonetwork whose load dependence is the phi(n) transform of the original
 * {f(n)}. See _kb/03-api-layer.md (pfqn/ family, panaceald).
 */
public final class Pfqn_panaceald {
    private Pfqn_panaceald() {}

    public static Ret.pfqnNc pfqn_panaceald(Matrix L, Matrix N, Matrix Z, Matrix mu) {
        return pfqn_panaceald(L, N, Z, mu, 3);
    }

    /**
     * Compute the load-dependent PANACEA approximation
     */
    public static Ret.pfqnNc pfqn_panaceald(Matrix L, Matrix N, Matrix Z, Matrix mu, int terms) {
        if (terms < 1 || terms > 3) {
            throw new RuntimeException("The terms parameter must be 1, 2, or 3 (higher-order coefficients are not implemented).");
        }
        int M = L.getNumRows();
        int R = L.getNumCols();
        String method = "panaceald";
        int Ntot = (int) FastMath.round(N.elementSum());
        if (Ntot == 0) {
            return new Ret.pfqnNc(1.0, 0.0, method);
        }

        double[] Ztot = new double[R];
        if (Z != null && !Z.isEmpty()) {
            Matrix Zs = Z.sumCols();
            for (int j = 0; j < R && j < Zs.length(); j++) {
                Ztot[j] = Zs.get(j);
            }
        }

        // Type-3 (infinite-server) rows are absent from the pseudonetwork and
        // enter only through rho_j0; Solver_ncld encodes them as mu(i,n)=n rows.
        boolean[] isIS = new boolean[M];
        for (int i = 0; i < M; i++) {
            isIS[i] = true;
            for (int n = 0; n < Ntot; n++) {
                if (FastMath.abs(muAt(mu, i, n) - (n + 1)) >= GlobalConstants.FineTol) {
                    isIS[i] = false;
                    break;
                }
            }
        }
        int Mq = 0;
        for (int i = 0; i < M; i++) {
            if (isIS[i]) {
                for (int j = 0; j < R; j++) {
                    Ztot[j] += L.get(i, j);
                }
            } else {
                Mq++;
            }
        }

        for (int j = 0; j < R; j++) {
            if (N.get(j) > 0 && Ztot[j] <= 0) {
                // no IS center on the route of a populated class: rho_j0 is
                // undefined and PANACEA does not apply
                return new Ret.pfqnNc(Double.NaN, Double.NaN, method);
            }
        }

        double lGdelay = -Matrix.factln(N).elementSum();
        for (int j = 0; j < R; j++) {
            if (N.get(j) != 0) {
                lGdelay += N.get(j) * FastMath.log(Ztot[j]);
            }
        }
        if (Mq == 0) {
            return new Ret.pfqnNc(FastMath.exp(lGdelay), lGdelay, method);
        }

        double[][] r = new double[Mq][R];
        double[][] muq = new double[Mq][Ntot];
        int qi = 0;
        for (int i = 0; i < M; i++) {
            if (isIS[i]) continue;
            for (int j = 0; j < R; j++) {
                r[qi][j] = Ztot[j] > 0 ? L.get(i, j) / Ztot[j] : 0.0;
            }
            for (int n = 0; n < Ntot; n++) {
                muq[qi][n] = muAt(mu, i, n);
                if (muq[qi][n] <= 0 || !Double.isFinite(muq[qi][n])) {
                    return new Ret.pfqnNc(Double.NaN, Double.NaN, method);
                }
            }
            qi++;
        }

        double[] lambda = new double[Mq];
        double[] alpha = new double[Mq];
        for (int i = 0; i < Mq; i++) {
            for (int j = 0; j < R; j++) {
                lambda[i] += N.get(j) * r[i][j];
            }
            alpha[i] = 1 - lambda[i] / muq[i][Ntot - 1];
            if (alpha[i] <= 0) {
                // model is not in normal usage: the {phi(n)} series diverges
                return new Ret.pfqnNc(Double.NaN, Double.NaN, method);
            }
        }

        // log-partial products log prod_{k=1}^{s} mu_i(k), s=0..Ntot
        double[][] lPi = new double[Mq][Ntot + 1];
        for (int i = 0; i < Mq; i++) {
            for (int s = 1; s <= Ntot; s++) {
                lPi[i][s] = lPi[i][s - 1] + FastMath.log(muq[i][s - 1]);
            }
        }

        int nmax = 2 * (terms - 1);
        double[][] lpsi = new double[Mq][nmax + 1];
        for (int i = 0; i < Mq; i++) {
            for (int n = 0; n <= nmax; n++) {
                lpsi[i][n] = logpsi(n, lambda[i], lPi[i], muq[i][Ntot - 1], alpha[i], Ntot);
            }
        }

        // load dependence of the pseudonetwork centers:
        // psi_i(n) = psi_i(0) n! / prod_{k=1}^{n} mups_i(k)
        double[][] mups = new double[Mq][FastMath.max(1, nmax)];
        for (int i = 0; i < Mq; i++) {
            for (int n = 1; n <= nmax; n++) {
                mups[i][n - 1] = FastMath.exp(FastMath.log(n) + lpsi[i][n - 1] - lpsi[i][n]);
            }
        }

        // Expansion coefficients (5.4). The large parameter N cancels
        // identically between beta_j=K_j/N, Gamma=N*r and the 1/N^n scaling, so
        // the demands are taken as r and beta as N.
        double[] A = new double[]{1.0, 0.0, 0.0};
        if (terms >= 2) {
            for (int j = 0; j < R; j++) {
                int[] k = new int[R];
                k[j] = 2;
                A[1] -= N.get(j) * pseudonet(r, k, mups);
            }
        }
        if (terms >= 3) {
            for (int j = 0; j < R; j++) {
                int[] k = new int[R];
                k[j] = 3;
                A[2] += 2 * N.get(j) * pseudonet(r, k, mups);
                k[j] = 4;
                A[2] += 3 * FastMath.pow(N.get(j), 2) * pseudonet(r, k, mups);
                for (int s = 0; s < R; s++) {
                    if (s == j) continue;
                    int[] k2 = new int[R];
                    k2[j] = 2;
                    k2[s] = 2;
                    A[2] += 0.5 * N.get(j) * N.get(s) * pseudonet(r, k2, mups);
                }
            }
        }
        double I = 0.0;
        for (int n = 0; n < terms; n++) I += A[n];
        if (I <= 0) {
            return new Ret.pfqnNc(Double.NaN, Double.NaN, method);
        }

        double lG = lGdelay + FastMath.log(I);
        for (int i = 0; i < Mq; i++) lG += lpsi[i][0];
        if (!Double.isFinite(lG)) {
            return new Ret.pfqnNc(Double.NaN, Double.NaN, method);
        }
        return new Ret.pfqnNc(FastMath.exp(lG), lG, method);
    }

    private static double muAt(Matrix mu, int i, int n) {
        if (mu == null || mu.isEmpty()) return 1.0;
        int cols = mu.getNumCols();
        return mu.get(i, FastMath.min(n, cols - 1));
    }

    /**
     * log of psi(n) = sum_{s&gt;=n} [s!/(s-n)!] lambda^(s-n) / prod_k mu(k), the
     * mu-free part of the phi(n) transform in eq. (3.7)-(3.8a). The series is
     * split into the exact head s&lt;=K and a geometric tail summed in closed form
     * via the Vandermonde identity, all terms positive.
     */
    private static double logpsi(int n, double lambda, double[] lPirow, double muK, double alpha, int K) {
        int head = FastMath.max(0, K - n + 1);
        double[] t = new double[head + n + 1];
        int c = 0;
        for (int s = n; s <= K; s++) {
            t[c++] = Maths.factln(s) - Maths.factln(s - n) + xlogy(s - n, lambda) - lPirow[s];
        }
        int T = FastMath.max(n, K + 1);
        for (int i = 0; i <= n; i++) {
            t[c++] = Maths.factln(n) + Maths.factln(T) - Maths.factln(n - i) - Maths.factln(T - n + i)
                    + xlogy(T + i - n, lambda) + (K - T - i) * FastMath.log(muK)
                    - (i + 1) * FastMath.log(alpha) - lPirow[K];
        }
        double tmax = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < c; i++) tmax = FastMath.max(tmax, t[i]);
        if (!Double.isFinite(tmax)) return tmax;
        double acc = 0.0;
        for (int i = 0; i < c; i++) acc += FastMath.exp(t[i] - tmax);
        return tmax + FastMath.log(acc);
    }

    /**
     * Partition function of the pseudonetwork at population k, normalized so
     * that G(0)=1. Populations are at most 2*(terms-1), so a direct
     * load-dependent convolution over the population lattice is used.
     */
    private static double pseudonet(double[][] gam, int[] k, double[][] mups) {
        int Mq = gam.length;
        int nz = 0;
        for (int j = 0; j < k.length; j++) if (k[j] > 0) nz++;
        int[] cls = new int[nz];
        int[] sizes = new int[nz];
        int c = 0;
        for (int j = 0; j < k.length; j++) {
            if (k[j] > 0) {
                cls[c] = j;
                sizes[c] = k[j] + 1;
                c++;
            }
        }
        int npop = 1;
        for (int j = 0; j < nz; j++) npop *= sizes[j];

        double[][] sterm = new double[Mq][npop];
        int[] m = new int[nz];
        for (int i = 0; i < Mq; i++) {
            for (int idx = 0; idx < npop; idx++) {
                idx2vec(idx, sizes, m);
                int sm = 0;
                for (int j = 0; j < nz; j++) sm += m[j];
                double v = Maths.factln(sm);
                boolean zero = false;
                for (int j = 0; j < nz; j++) {
                    if (m[j] > 0) {
                        double g = gam[i][cls[j]];
                        if (g <= 0) { zero = true; break; }
                        v += m[j] * FastMath.log(g) - Maths.factln(m[j]);
                    }
                }
                if (zero) {
                    sterm[i][idx] = 0.0;
                } else {
                    for (int l = 0; l < sm; l++) v -= FastMath.log(mups[i][l]);
                    sterm[i][idx] = FastMath.exp(v);
                }
            }
        }

        double[] gprev = new double[npop];
        double[] gcur = new double[npop];
        gprev[0] = 1.0;
        int[] nvec = new int[nz];
        for (int i = 0; i < Mq; i++) {
            for (int idx = 0; idx < npop; idx++) {
                idx2vec(idx, sizes, nvec);
                double acc = 0.0;
                for (int jdx = 0; jdx < npop; jdx++) {
                    idx2vec(jdx, sizes, m);
                    boolean fits = true;
                    for (int j = 0; j < nz; j++) {
                        if (m[j] > nvec[j]) { fits = false; break; }
                    }
                    if (!fits) continue;
                    int rem = 0;
                    int mult = 1;
                    for (int j = 0; j < nz; j++) {
                        rem += mult * (nvec[j] - m[j]);
                        mult *= sizes[j];
                    }
                    acc += gprev[rem] * sterm[i][jdx];
                }
                gcur[idx] = acc;
            }
            System.arraycopy(gcur, 0, gprev, 0, npop);
        }
        return gprev[npop - 1];
    }

    private static void idx2vec(int idx, int[] sizes, int[] out) {
        int t = idx;
        for (int j = 0; j < sizes.length; j++) {
            out[j] = t % sizes[j];
            t /= sizes[j];
        }
    }

    private static double xlogy(int e, double x) {
        return e == 0 ? 0.0 : e * FastMath.log(x);
    }
}
