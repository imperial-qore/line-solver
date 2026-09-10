/**
 * @file Norlund-Rice Probit (NRP) method for normalizing constants
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.util.FastMath;

import java.util.function.Function;

import jline.api.pfqn.ld.Pfqn_gld;
import jline.api.pfqn.ld.Pfqn_gld_complex;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.SerializableFunction;
import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;

public final class Pfqn_nrp {
    private Pfqn_nrp() {}

    public static double pfqn_nrp(Matrix Lin, Matrix N, Matrix Z, Matrix alphaIn, SolverOptions options) {
        Matrix L = Lin;
        Matrix alpha = alphaIn;
        double Nt = N.elementSum();
        if (Z.elementSum() > 0) {
            L = Matrix.concatRows(L, Z, null);
            // the delay is an infinite server: rates 1..Nt along row 0 of a
            // single-row block, which is the station concatRows appends
            Matrix alpha_tmp = new Matrix(1, (int) Nt);
            for (int k = 0; k < (int) Nt; k++) {
                alpha_tmp.set(0, k, k + 1.0);
            }
            alpha = Matrix.concatRows(alpha, alpha_tmp, null);
        }
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (M == 1 && Z.elementSum() == 0.0) {
            return Pfqn_gld.pfqn_gld(L, N, alpha, options).lG;
        }
        Matrix Lmax_mat = new Matrix(1, R);
        for (int col = 0; col < R; col++) {
            double colMax = Double.NEGATIVE_INFINITY;
            for (int row = 0; row < L.getNumRows(); row++) {
                if (L.get(row, col) > colMax) {
                    colMax = L.get(row, col);
                }
            }
            Lmax_mat.set(0, col, colMax == 0.0 ? 1.0 : colMax);
        }
        Matrix L_final = L.elementDivide(Lmax_mat.repmat(L.getNumRows(), 1));
        Matrix x0 = new Matrix(1, R);
        x0.zero();
        double lG = Maths.laplaceapprox_h_complex(x0, infradius_hnorm(L_final, N, alpha)).logI.getReal();
        Matrix Lmax_log_mat = Lmax_mat.copy();
        for (int i = 0; i < Lmax_mat.getNumElements(); i++) {
            Lmax_log_mat.set(i, FastMath.log(Lmax_mat.get(i)));
        }
        return lG + N.mult(Lmax_log_mat.transpose()).get(0);
    }

    public static SerializableFunction<Matrix, ComplexMatrix> infradius_hnorm(final Matrix L,
                                                                              final Matrix N,
                                                                              final Matrix alpha) {
        return new SerializableFunction<Matrix, ComplexMatrix>() {
            @Override
            public ComplexMatrix apply(Matrix x) {
                double MU = 0.0;
                double SIGMA = 1.0;
                double Nt = N.elementSum();
                Matrix beta = new Matrix(N.getNumRows(), N.getNumCols());
                N.divide(Nt, beta, true);
                NormalDistribution Z = new NormalDistribution(MU, SIGMA);

                Matrix t = x.copy();
                for (int i = 0; i < t.getNumRows(); i++) {
                    for (int j = 0; j < t.getNumCols(); j++) {
                        t.set(i, j, Z.cumulativeProbability(x.get(i, j)));
                    }
                }

                double tb = 0.0;
                for (int i = 0; i < t.getNumRows(); i++) {
                    for (int j = 0; j < t.getNumCols(); j++) {
                        tb += beta.get(i, j) * t.get(i, j);
                    }
                }

                Matrix tsubtb = t.copy();
                for (int i = 0; i < t.getNumElements(); i++) {
                    tsubtb.set(i, tsubtb.get(i) - tb);
                }

                ComplexMatrix y = new ComplexMatrix(new Matrix(x.getNumRows(), 1), new Matrix(x.getNumRows(), 1));
                Function<Matrix, Complex> h = nrp_h(L, tsubtb, Nt, alpha);
                for (int i = 0; i < x.getNumRows(); i++) {
                    Matrix xi = Matrix.extractRows(x, i, x.getNumRows(), null);
                    y.set(i, h.apply(xi));
                }

                return new ComplexMatrix(y.real);
            }
        };
    }

    public static Function<Matrix, Complex> nrp_h(final Matrix L,
                                                   final Matrix tsubtb,
                                                   final double Nt,
                                                   final Matrix alpha) {
        return new Function<Matrix, Complex>() {
            @Override
            public Complex apply(Matrix x) {
                int M = L.getNumRows();
                ComplexMatrix c = new ComplexMatrix(tsubtb.getNumRows(), tsubtb.getNumCols());

                for (int i = 0; i < c.getNumElements(); i++) {
                    c.set(i, new Complex(0.0, 2 * Math.PI * tsubtb.get(i)).exp());
                }

                ComplexMatrix LComplex = new ComplexMatrix(
                        L.elementMult(c.real.repmat(M, 1), null),
                        L.elementMult(c.im.repmat(M, 1), null));

                NormalDistribution Z = new NormalDistribution();
                Matrix normpdfX = x.copy();
                for (int i = 0; i < x.getNumElements(); i++) {
                    normpdfX.set(i, Z.density(x.get(i)));
                }

                return Pfqn_gld_complex.pfqn_gld_complex(
                        LComplex.sumRows(),
                        Matrix.singleton(Nt),
                        alpha,
                        null
                ).G.multiply(normpdfX.elementMult());
            }
        };
    }
}
