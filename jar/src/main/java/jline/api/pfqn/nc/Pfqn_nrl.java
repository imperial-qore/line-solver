/**
 * @file Norlund-Rice Logit (NRL) method for normalizing constants
 *
 * Implements the Norlund-Rice Logit approach for computing normalizing constants
 * using Laplace approximation with complex-valued service demands. Provides high-accuracy
 * approximation for closed networks through sophisticated lattice-based integration.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

import java.util.function.Function;

import jline.GlobalConstants;
import jline.api.pfqn.ld.Pfqn_gld;
import jline.api.pfqn.ld.Pfqn_gld_complex;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.SerializableFunction;
import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;

public final class Pfqn_nrl {
    private Pfqn_nrl() {}

    public static double pfqn_nrl(Matrix Lin, Matrix N, Matrix Z, Matrix alphaIn, SolverOptions options) {
        Matrix L = Lin;
        Matrix alpha = alphaIn;
        double Nt = N.elementSum();
        if (Nt < 0) {
            return GlobalConstants.NegInf;
        }
        if (Nt == 0.0) {
            return 0.0;
        }
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (Z.elementSum() > 0) {
            L = Matrix.concatRows(L, Z, null);
            java.util.List<Double> Ntrange = new java.util.ArrayList<Double>();
            double i = 1.0;
            while (i <= Nt) {
                Ntrange.add(i);
                i++;
            }
            alpha = Matrix.concatRows(alpha, new Matrix(Ntrange), null);
        }
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
        double lG = Maths.laplaceapprox_h_complex(x0, infradius_h(L_final, N, alpha)).logI.getReal();
        Matrix Lmax_log_mat = Lmax_mat.copy();
        for (int i = 0; i < Lmax_mat.getNumElements(); i++) {
            Lmax_log_mat.set(i, FastMath.log(Lmax_mat.get(i)));
        }
        return lG + N.mult(Lmax_log_mat.transpose()).get(0);
    }

    public static SerializableFunction<Matrix, ComplexMatrix> infradius_h(final Matrix L,
                                                                          final Matrix N,
                                                                          final Matrix alpha) {
        return new SerializableFunction<Matrix, ComplexMatrix>() {
            @Override
            public ComplexMatrix apply(Matrix x) {
                double Nt = N.elementSum();
                Matrix beta = new Matrix(N.getNumRows(), N.getNumCols());
                N.divide(Nt, beta, true);
                Matrix t = new Matrix(x.getNumRows(), x.getNumCols());
                Matrix tbMat = new Matrix(x.getNumRows(), 1);

                for (int i = 0; i < x.getNumRows(); i++) {
                    double tbRow = 0.0;
                    for (int j = 0; j < x.getNumCols(); j++) {
                        double expVal = Math.exp(x.get(i, j));
                        double tVal = expVal / (1.0 + expVal);
                        t.set(i, j, tVal);
                        tbRow += beta.get(i, j) * tVal;
                    }
                    tbMat.set(i, 0, tbRow);
                }

                double tb = tbMat.elementSum();
                Matrix tsubtb = t.copy();
                for (int i = 0; i < t.getNumElements(); i++) {
                    tsubtb.set(i, tsubtb.get(i) - tb);
                }

                Function<Matrix, Complex> h = nrl_h(L, tsubtb, Nt, alpha);
                ComplexMatrix y = new ComplexMatrix(new Matrix(x.getNumRows(), 1), new Matrix(x.getNumRows(), 1));
                for (int i = 0; i < x.getNumRows(); i++) {
                    Matrix xi = Matrix.extractRows(x, i, x.getNumRows(), null);
                    y.set(i, h.apply(xi));
                }

                return y;
            }
        };
    }

    public static Function<Matrix, Complex> nrl_h(final Matrix L,
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

                Matrix expX = x.copy();
                for (int i = 0; i < x.getNumElements(); i++) {
                    double expVal = Math.exp(x.get(i));
                    expX.set(i, expVal / ((1 + expVal) * (1 + expVal)));
                }

                return Pfqn_gld_complex.pfqn_gld_complex(
                        LComplex.sumRows(),
                        Matrix.singleton(Nt),
                        alpha,
                        null
                ).G.multiply(expX.elementMult());
            }
        };
    }
}
