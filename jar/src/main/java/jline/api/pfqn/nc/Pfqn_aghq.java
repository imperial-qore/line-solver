/**
 * @file Adaptive Gauss-Hermite quadrature of the McKenna-Mitra integral
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

/**
 * Adaptive Gauss-Hermite quadrature of the simplex factor of the McKenna-Mitra integral.
 *
 * Rescaling by the logistic-expansion mode and curvature, w = w* + A^(-1/2)*z, and applying
 * the q-node probabilists' Gauss-Hermite rule in each of the M-1 directions gives a
 * convergent rule whose q = 1 member is Pfqn_le itself: a single node at the mode with
 * weight sqrt(2*pi). So LE is the first term of a convergent quadrature rather than an
 * approximation of unknown accuracy. The cost is q^(M-1) evaluations, which is what confines
 * the method to small M.
 *
 * A tensor rule is NOT invariant to the choice of square root of A: any B with B*B' = inv(A)
 * is admissible and they place the nodes differently. The principal-axis frame from the
 * eigendecomposition is used here, as in the reference results; where two curvatures are
 * close to equal the frame is close to arbitrary and two valid rules can part company well
 * above their own error, converging back together as q grows. Do not compare this across
 * codebases node by node.
 *
 * With Z &gt; 0 the radius is integrated numerically and the rule is applied
 * to the M-1 simplex directions, so every node costs one radial quadrature. Because the
 * radius is integrated rather than Laplaced, q = 1 there is the logistic expansion with an
 * exact radius, which is NOT Pfqn_le's own Z &gt; 0 branch.
 *
 * References:
 * J. McKenna, D. Mitra, "Integral Representations and Asymptotic Expansions for Closed
 * Markovian Queueing Networks: Normal Usage", Bell Syst. Tech. J. 61(5), 1982.
 * Q. Liu, D. A. Pierce, "A Note on Gauss-Hermite Quadrature", Biometrika 81(3), 1994.
 *
 * @since LINE 3.0
 */
public final class Pfqn_aghq {
    private Pfqn_aghq() {}

    public static Ret.pfqnNc pfqn_aghq(Matrix L, Matrix N, Matrix Z, int q) {
        final int M = L.getNumRows();
        int R = L.getNumCols();
        double lGn;

        if (L.isEmpty() || N.isEmpty() || N.elementSum() == 0.0
                || L.elementSum() < GlobalConstants.CoarseTol) {
            return Pfqn_le.pfqn_le(L, N, Z);
        }

        final double[] Nv = Pfqn_simplex.toVector(N, N.length());
        final double[] Zv = Pfqn_simplex.toVector(Z, N.length());
        final double[][] La = Pfqn_simplex.toArray(L);
        int d = M - 1;

        if (Z == null || Z.isEmpty() || Z.elementSum() < GlobalConstants.Zero) {
            Ret.pfqnLeFpi ret = Pfqn_le_fpi.pfqn_le_fpi(L, N);
            Matrix umax = ret.u;
            double[] u = Pfqn_simplex.toVector(umax, M);
            double[][] A =
                    Pfqn_simplex.toArray(Pfqn_le_hessian.pfqn_le_hessian(L, N, umax.transpose()));
            double ld = Pfqn_simplex.logdet(A);
            double S = 0.0;
            double sumLu = 0.0;
            for (int r = 0; r < R; r++) {
                double c = 0.0;
                for (int i = 0; i < M; i++) {
                    c += u[i] * L.get(i, r);
                }
                S += Nv[r] * FastMath.log(c);
            }
            for (int i = 0; i < M; i++) {
                sumLu += FastMath.log(u[i]);
            }
            double h0 = S + sumLu;
            double[] w0 = new double[d];
            for (int i = 0; i < d; i++) {
                w0[i] = FastMath.log(u[i] / u[M - 1]);
            }
            Pfqn_simplex.LogIntegrand h = new Pfqn_simplex.LogIntegrand() {
                public double at(double[] w) {
                    double[] x = Pfqn_simplex.softmaxGauge(w);
                    double[] c = Pfqn_simplex.matVec(x, La, Nv.length);
                    double acc = 0.0;
                    for (int r = 0; r < Nv.length; r++) {
                        acc += Nv[r] * FastMath.log(c[r]);
                    }
                    for (int i = 0; i < M; i++) {
                        acc += FastMath.log(x[i]);
                    }
                    return acc;
                }
            };
            double lacc = Pfqn_simplex.aghqRule(h, w0, h0, A, q, d);
            Matrix tmp = new Matrix(1, N.length() + 1);
            Matrix.extract(N, 0, 1, 0, N.length(), tmp, 0, 0);
            tmp.set(N.length(), (double) (M - 1));
            lGn = Maths.multinomialln(tmp) + Maths.factln(M - 1) + h0 + lacc - 0.5 * ld;
        } else {
            Pfqn_simplex.Mode mode = Pfqn_simplex.simplexMode(La, Nv, Zv);
            double[] w0 = new double[d];
            for (int i = 0; i < d; i++) {
                w0[i] = FastMath.log(mode.x[i] / mode.x[M - 1]);
            }
            Pfqn_simplex.LogIntegrand h = new Pfqn_simplex.LogIntegrand() {
                public double at(double[] w) {
                    double[] x = Pfqn_simplex.softmaxGauge(w);
                    double acc = Pfqn_simplex.radial(
                            Pfqn_simplex.matVec(x, La, Nv.length), Nv, Zv, M).lJ;
                    for (int i = 0; i < M; i++) {
                        acc += FastMath.log(x[i]);
                    }
                    return acc;
                }
            };
            double lacc = Pfqn_simplex.aghqRule(h, w0, mode.h0, mode.A, q, d);
            lGn = -Matrix.factln(N).elementSum() + mode.h0 + lacc - 0.5 * mode.ld;
        }
        return new Ret.pfqnNc(Double.valueOf(FastMath.exp(lGn)), Double.valueOf(lGn));
    }

    public static Ret.pfqnNc pfqn_aghq(Matrix L, Matrix N, Matrix Z) {
        return pfqn_aghq(L, N, Z, 3);
    }

    public static Ret.pfqnNc pfqn_aghq(Matrix L, Matrix N) {
        return pfqn_aghq(L, N, new Matrix(1, L.getNumCols()), 3);
    }
}
