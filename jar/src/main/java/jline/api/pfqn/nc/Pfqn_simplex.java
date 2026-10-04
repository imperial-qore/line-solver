/**
 * @file Shared machinery for closures of the simplex factor of the McKenna-Mitra integral
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Pieces backing the closure of the simplex factor, used by Pfqn_aghq.
 *
 * At Z = 0 the McKenna-Mitra integrand is homogeneous, so y = v*x separates and the radius
 * integrates exactly to gamma(N+M), leaving an integral over the unit simplex which is where
 * all of the error of the logistic expansion lives. With Z &gt; 0 that factorisation is gone
 * and the radius is integrated numerically here rather than closed, leaving the same M-1
 * simplex directions to a closure.
 *
 * See _kb/03-api-layer.md.
 */
final class Pfqn_simplex {
    private Pfqn_simplex() {}

    static final double TINY = Double.MIN_NORMAL;

    /** Log-integrand of the simplex factor, evaluated in logistic coordinates. */
    interface LogIntegrand {
        double at(double[] w);
    }

    /** log J(c) together with the moments of the tilted law of the radius. */
    static final class Radial {
        double lJ;
        double[] G;
        double vbar;
        double[][] Lam;
    }

    /** Mode, curvature and log-integrand at the mode of the simplex factor. */
    static final class Mode {
        double[] x;
        double[][] A;
        double ld;
        double h0;
    }

    /** Nodes and weights of the Gauss rule with the given Jacobi off-diagonal. */
    private static double[][] golubWelsch(double[] offDiag, double mu0) {
        int n = offDiag.length + 1;
        double[][] J = new double[n][n];
        for (int k = 0; k < n - 1; k++) {
            J[k][k + 1] = offDiag[k];
            J[k + 1][k] = offDiag[k];
        }
        EigenDecomposition ed = new EigenDecomposition(new Array2DRowRealMatrix(J, false));
        double[] x = new double[n];
        double[] w = new double[n];
        Integer[] order = new Integer[n];
        double[] ev = new double[n];
        for (int k = 0; k < n; k++) {
            ev[k] = ed.getRealEigenvalue(k);
            order[k] = Integer.valueOf(k);
        }
        java.util.Arrays.sort(order, new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Double.compare(ev[a.intValue()], ev[b.intValue()]);
            }
        });
        for (int k = 0; k < n; k++) {
            int j = order[k].intValue();
            x[k] = ev[j];
            double v0 = ed.getEigenvector(j).getEntry(0);
            w[k] = mu0 * v0 * v0;
        }
        return new double[][] {x, w};
    }

    /** N-point Gauss-Legendre rule on [-1,1]. */
    static double[][] gaussLegendre(int n) {
        double[] b = new double[n - 1];
        for (int k = 1; k < n; k++) {
            b[k - 1] = k / FastMath.sqrt(4.0 * k * k - 1.0);
        }
        return golubWelsch(b, 2.0);
    }

    /** Q-point Gauss-Hermite rule of the probabilists' weight exp(-z^2/2). */
    static double[][] gaussHermite(int q) {
        if (q == 1) {
            return new double[][] {{0.0}, {FastMath.sqrt(2 * FastMath.PI)}};
        }
        double[] b = new double[q - 1];
        for (int k = 1; k < q; k++) {
            b[k - 1] = FastMath.sqrt((double) k);
        }
        return golubWelsch(b, FastMath.sqrt(2 * FastMath.PI));
    }

    private static final double[][] GL64 = gaussLegendre(64);

    /** Log-determinant of a positive-definite matrix; 0 for the empty matrix. */
    static double logdet(double[][] A) {
        int n = A.length;
        if (n == 0) {
            return 0.0;
        }
        double[][] C = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double s = A[i][j];
                for (int k = 0; k < j; k++) {
                    s -= C[i][k] * C[j][k];
                }
                if (i == j) {
                    if (s <= 0.0) {
                        // Not positive definite: fall back to the eigenvalues, which the
                        // caller still needs a finite number from.
                        RealMatrix rm = new Array2DRowRealMatrix(A, true);
                        EigenDecomposition ed = new EigenDecomposition(rm);
                        double ld = 0.0;
                        for (int t = 0; t < n; t++) {
                            ld += FastMath.log(FastMath.abs(ed.getRealEigenvalue(t)));
                        }
                        return ld;
                    }
                    C[i][j] = FastMath.sqrt(s);
                } else {
                    C[i][j] = s / C[j][j];
                }
            }
        }
        double ld = 0.0;
        for (int i = 0; i < n; i++) {
            ld += 2.0 * FastMath.log(C[i][i]);
        }
        return ld;
    }

    /** Log-integrand of the radial integral in t = log v, Jacobian included. */
    static double radialLogf(double t, double[] c, double[] N, double[] Z, int M) {
        double v = FastMath.exp(t);
        double f = -v + M * t;
        for (int r = 0; r < c.length; r++) {
            f += N[r] * FastMath.log(FastMath.max(Z[r] + v * c[r], TINY));
        }
        return f;
    }

    /**
     * log J(c) = log int_0^inf exp(-v) v^(M-1) prod_r (Z_r + v*c_r)^N_r dv, plus the moments
     * of the tilted law of v that the simplex derivatives need: G_r = E[T_r], vbar = E[v]
     * and Lam = cov(T) - diag(E[T^2]/N) = grad^2_c log J, with T_r(v) = N_r*v/(Z_r + v*c_r).
     * Quadrature runs in t = log v, where the integrand is bounded at both ends, over two
     * Gauss-Legendre panels meeting at the mode.
     */
    static Radial radial(double[] c, double[] N, double[] Z, int M) {
        int R = c.length;
        double t = FastMath.log(sum(N) + M);
        for (int it = 0; it < 200; it++) {
            double v = FastMath.exp(t);
            double f1 = -v + M;
            double f2 = -v;
            for (int r = 0; r < R; r++) {
                double d = FastMath.max(Z[r] + v * c[r], TINY);
                f1 += N[r] * (v * c[r]) / d;
                f2 += N[r] * (v * c[r]) * Z[r] / (d * d);
            }
            if (f2 > -1e-300) {
                break;
            }
            double step = FastMath.max(FastMath.min(-f1 / f2, 2.0), -2.0);
            t += step;
            if (FastMath.abs(step) < 1e-13) {
                break;
            }
        }
        double v = FastMath.exp(t);
        double f2 = -v;
        for (int r = 0; r < R; r++) {
            double d = FastMath.max(Z[r] + v * c[r], TINY);
            f2 += N[r] * (v * c[r]) * Z[r] / (d * d);
        }
        double sig = f2 < -1e-300 ? 1.0 / FastMath.sqrt(-f2) : 1.0;
        double fm = radialLogf(t, c, N, Z, M);
        // Widen each half-window until the log-integrand has fallen by 60 nats, so the
        // discarded tails are below 1e-26 in relative terms.
        double a = FastMath.min(12.0 * sig, t + 745.0);
        for (int k = 0; k < 60; k++) {
            if (t - a <= -745.0 || radialLogf(t - a, c, N, Z, M) < fm - 60.0) {
                break;
            }
            a = FastMath.min(1.6 * a, t + 745.0);
        }
        double b = 12.0 * sig;
        for (int k = 0; k < 60; k++) {
            if (radialLogf(t + b, c, N, Z, M) < fm - 60.0) {
                break;
            }
            b *= 1.6;
        }
        double[] vg = GL64[0];
        double[] wg = GL64[1];
        int nq = 2 * vg.length;
        double[] tt = new double[nq];
        double[] W = new double[nq];
        for (int k = 0; k < vg.length; k++) {
            tt[k] = 0.5 * a * vg[k] + (t - 0.5 * a);
            W[k] = 0.5 * a * wg[k];
            tt[vg.length + k] = 0.5 * b * vg[k] + (t + 0.5 * b);
            W[vg.length + k] = 0.5 * b * wg[k];
        }
        double[] vv = new double[nq];
        double[][] D = new double[nq][R];
        double[] fv = new double[nq];
        double mx = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < nq; k++) {
            vv[k] = FastMath.exp(tt[k]);
            double f = -vv[k] + M * tt[k];
            for (int r = 0; r < R; r++) {
                D[k][r] = FastMath.max(Z[r] + vv[k] * c[r], TINY);
                f += N[r] * FastMath.log(D[k][r]);
            }
            fv[k] = f;
            if (f > mx) {
                mx = f;
            }
        }
        double se = 0.0;
        double[] e = new double[nq];
        for (int k = 0; k < nq; k++) {
            e[k] = W[k] * FastMath.exp(fv[k] - mx);
            se += e[k];
        }
        Radial out = new Radial();
        out.lJ = mx + FastMath.log(se);
        double[] p = new double[nq];
        for (int k = 0; k < nq; k++) {
            p[k] = e[k] / se;
        }
        double[][] T = new double[nq][R];
        out.G = new double[R];
        out.vbar = 0.0;
        for (int k = 0; k < nq; k++) {
            out.vbar += p[k] * vv[k];
            for (int r = 0; r < R; r++) {
                T[k][r] = N[r] * vv[k] / D[k][r];
                out.G[r] += p[k] * T[k][r];
            }
        }
        double[][] et2 = new double[R][R];
        for (int k = 0; k < nq; k++) {
            for (int r = 0; r < R; r++) {
                double pt = p[k] * T[k][r];
                for (int s = 0; s < R; s++) {
                    et2[r][s] += pt * T[k][s];
                }
            }
        }
        out.Lam = new double[R][R];
        for (int r = 0; r < R; r++) {
            for (int s = 0; s < R; s++) {
                out.Lam[r][s] = et2[r][s] - out.G[r] * out.G[s];
            }
            if (N[r] > 0) {
                out.Lam[r][r] -= et2[r][r] / N[r];
            }
        }
        symmetrize(out.Lam);
        return out;
    }

    /**
     * Mode and curvature of h(w) = log J(L'*x(w)) + sum_i log x_i with J the exact radial
     * integral. The fixed point x = (1 + x.*(L*G))/vbar is the Z &gt; 0 analogue of
     * pfqn_le_fpi: integrating by parts gives sum_i x_i*(L*G)_i = vbar - M, so the update is
     * normalised by construction, and at Z = 0 it reduces to pfqn_le_fpi. The term in the
     * second derivative of x(w) drops at the mode against sum_i x_i == 1.
     */
    static Mode simplexMode(double[][] L, double[] N, double[] Z) {
        int M = L.length;
        int R = N.length;
        Ret.pfqnLeFpiZ start = Pfqn_le_fpiZ.pfqn_le_fpiZ(toMatrix(L), toRow(N), toRow(Z));
        double[] x = new double[M];
        for (int i = 0; i < M; i++) {
            x[i] = start.u.get(i);
        }
        double[] x1 = new double[M];
        java.util.Arrays.fill(x1, Double.POSITIVE_INFINITY);
        for (int it = 0; it < 10000; it++) {
            double diff = 0.0;
            for (int i = 0; i < M; i++) {
                diff += FastMath.abs(x[i] - x1[i]);
            }
            if (diff <= 1e-11) {
                break;
            }
            x1 = x.clone();
            Radial rad = radial(matVec(x1, L, R), N, Z, M);
            double s = 0.0;
            for (int i = 0; i < M; i++) {
                double lg = 0.0;
                for (int r = 0; r < R; r++) {
                    lg += L[i][r] * rad.G[r];
                }
                x[i] = (1.0 + x1[i] * lg) / rad.vbar;
                s += x[i];
            }
            for (int i = 0; i < M; i++) {
                x[i] /= s;
            }
        }
        Radial rad = radial(matVec(x, L, R), N, Z, M);
        // P = L*Lam*L' - diag(1/x^2), then A = -Jm'*P*Jm with Jm = (diag(x)-x*x')(:,1:M-1).
        double[][] P = new double[M][M];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                double acc = 0.0;
                for (int r = 0; r < R; r++) {
                    double lr = 0.0;
                    for (int s = 0; s < R; s++) {
                        lr += rad.Lam[r][s] * L[j][s];
                    }
                    acc += L[i][r] * lr;
                }
                P[i][j] = acc;
            }
            P[i][i] -= 1.0 / (x[i] * x[i]);
        }
        int d = M - 1;
        double[][] Jm = new double[M][d];
        for (int i = 0; i < M; i++) {
            for (int aIdx = 0; aIdx < d; aIdx++) {
                Jm[i][aIdx] = (i == aIdx ? x[i] : 0.0) - x[i] * x[aIdx];
            }
        }
        double[][] A = new double[d][d];
        for (int aIdx = 0; aIdx < d; aIdx++) {
            for (int bIdx = 0; bIdx < d; bIdx++) {
                double acc = 0.0;
                for (int i = 0; i < M; i++) {
                    double pj = 0.0;
                    for (int j = 0; j < M; j++) {
                        pj += P[i][j] * Jm[j][bIdx];
                    }
                    acc += Jm[i][aIdx] * pj;
                }
                A[aIdx][bIdx] = -acc;
            }
        }
        symmetrize(A);
        Mode out = new Mode();
        out.x = x;
        out.A = A;
        out.ld = logdet(A);
        double sumLx = 0.0;
        for (int i = 0; i < M; i++) {
            sumLx += FastMath.log(x[i]);
        }
        out.h0 = rad.lJ + sumLx;
        return out;
    }

    /** softmax of [w; 0], the logistic parametrisation of the simplex with the gauge w_M = 0. */
    static double[] softmaxGauge(double[] w) {
        int M = w.length + 1;
        double[] a = new double[M];
        double mx = 0.0;
        for (int i = 0; i < w.length; i++) {
            a[i] = w[i];
            if (a[i] > mx) {
                mx = a[i];
            }
        }
        double s = 0.0;
        double[] x = new double[M];
        for (int i = 0; i < M; i++) {
            x[i] = FastMath.exp(a[i] - mx);
            s += x[i];
        }
        for (int i = 0; i < M; i++) {
            x[i] /= s;
        }
        return x;
    }

    /**
     * Log of the tensor Gauss-Hermite sum, accumulated with a running maximum. The
     * det(A)^(-1/2) of the rule is applied by the caller. A tensor rule is not invariant to
     * the choice of A^(-1/2); the principal-axis frame is used, as in the reference results.
     */
    static double aghqRule(LogIntegrand h, double[] w0, double h0, double[][] A, int q, int d) {
        if (d == 0) {
            return 0.0;
        }
        double nodesD = FastMath.pow((double) q, (double) d);
        if (nodesD > 1e7) {
            throw new RuntimeException("pfqn_aghq: the tensor rule needs q^(M-1)="
                    + (long) nodesD + " nodes; reduce q or use pfqn_le.");
        }
        int nodes = (int) nodesD;
        EigenDecomposition ed = new EigenDecomposition(new Array2DRowRealMatrix(A, true));
        double[][] B = new double[d][d];
        for (int j = 0; j < d; j++) {
            double lam = ed.getRealEigenvalue(j);
            if (lam <= 0.0) {
                return Double.NaN;
            }
            double sc = 1.0 / FastMath.sqrt(lam);
            for (int i = 0; i < d; i++) {
                B[i][j] = ed.getEigenvector(j).getEntry(i) * sc;
            }
        }
        double[][] gh = gaussHermite(q);
        double[] z = gh[0];
        double[] lwt = new double[q];
        for (int k = 0; k < q; k++) {
            lwt[k] = FastMath.log(gh[1][k]);
        }
        int[] idx = new int[d];
        double lmax = Double.NEGATIVE_INFINITY;
        double s = 0.0;
        double[] zz = new double[d];
        double[] w = new double[d];
        for (int k = 0; k < nodes; k++) {
            double lw = 0.0;
            double zsq = 0.0;
            for (int j = 0; j < d; j++) {
                zz[j] = z[idx[j]];
                lw += lwt[idx[j]];
                zsq += zz[j] * zz[j];
            }
            for (int i = 0; i < d; i++) {
                double acc = w0[i];
                for (int j = 0; j < d; j++) {
                    acc += B[i][j] * zz[j];
                }
                w[i] = acc;
            }
            double lt = lw + h.at(w) - h0 + 0.5 * zsq;
            if (lt > lmax) {
                s = s * FastMath.exp(lmax - lt) + 1.0;
                lmax = lt;
            } else {
                s += FastMath.exp(lt - lmax);
            }
            for (int j = d - 1; j >= 0; j--) {
                idx[j]++;
                if (idx[j] < q) {
                    break;
                }
                idx[j] = 0;
            }
        }
        return lmax + FastMath.log(s);
    }

    // ------------------------------------------------------------------ small utilities

    static double sum(double[] v) {
        double s = 0.0;
        for (int i = 0; i < v.length; i++) {
            s += v[i];
        }
        return s;
    }

    static void symmetrize(double[][] A) {
        for (int i = 0; i < A.length; i++) {
            for (int j = i + 1; j < A.length; j++) {
                double m = 0.5 * (A[i][j] + A[j][i]);
                A[i][j] = m;
                A[j][i] = m;
            }
        }
    }

    /** x' * L, i.e. the induced per-class demand at the simplex point x. */
    static double[] matVec(double[] x, double[][] L, int R) {
        double[] c = new double[R];
        for (int r = 0; r < R; r++) {
            double acc = 0.0;
            for (int i = 0; i < L.length; i++) {
                acc += x[i] * L[i][r];
            }
            c[r] = acc;
        }
        return c;
    }

    static double[][] toArray(Matrix L) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double[][] out = new double[M][R];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                out[i][r] = L.get(i, r);
            }
        }
        return out;
    }

    static double[] toVector(Matrix v, int n) {
        double[] out = new double[n];
        if (v != null && !v.isEmpty()) {
            for (int i = 0; i < FastMath.min(n, v.length()); i++) {
                out[i] = v.get(i);
            }
        }
        return out;
    }

    static Matrix toMatrix(double[][] A) {
        Matrix out = new Matrix(A.length, A.length == 0 ? 0 : A[0].length);
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[i].length; j++) {
                out.set(i, j, A[i][j]);
            }
        }
        return out;
    }

    static Matrix toRow(double[] v) {
        Matrix out = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            out.set(0, i, v[i]);
        }
        return out;
    }
}
