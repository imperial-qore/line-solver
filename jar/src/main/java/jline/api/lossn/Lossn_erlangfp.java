/**
 * @file Loss Network Analysis via Erlang Fixed Point
 *
 * @since LINE 3.0
 */
package jline.api.lossn;

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;

import jline.GlobalConstants;
import jline.api.da.Da_fpi;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Lossn_erlangfp {
    private Lossn_erlangfp() {}

    /**
     * Erlang fixed point approximation for loss networks.
     */
    public static Ret.lossnErlangFP lossn_erlangfp(Matrix nuVec, Matrix Amat, Matrix cVec) {
        final double[] nu = nuVec.toArray1D();
        final double[][] A = Amat.toArray2D();
        final double[] c = cVec.toArray1D();
        final int R = nu.length;
        final int J = c.length;
        double[] E0 = new double[J];
        Arrays.fill(E0, 0.5);

        // Erlang fixed point on the link blocking probabilities, driven by
        // the generic DA successive-substitution driver
        Da_fpi.Options<double[]> fpopts = new Da_fpi.Options<double[]>(
                Integer.MAX_VALUE, 1e-8, new Da_fpi.Norm<double[]>() {
            @Override
            public double eval(double[] xnew, double[] xref) {
                return MathArrays.distance(xnew, xref);
            }
        });
        fpopts.nanstop = true; // legacy while-loop exited on NaN convergence measure
        Da_fpi.Result<double[]> fpres = Da_fpi.run(new Da_fpi.Sweep<double[]>() {
            @Override
            public Da_fpi.SweepResult<double[]> sweep(double[] E_1, int it) {
                double[] Enew = new double[J];
                for (int j = 0; j < J; j++) {
                    double rhoj_1 = 0.0;
                    for (int r = 0; r < R; r++) {
                        if (A[j][r] > 0) {
                            double termj = nu[r] * A[j][r];
                            for (int i = 0; i < J; i++) {
                                if (A[i][r] > GlobalConstants.Zero) {
                                    termj *= FastMath.pow(1 - E_1[i], A[i][r]);
                                }
                            }
                            rhoj_1 += termj;
                        }
                    }
                    rhoj_1 /= (1 - E_1[j]);
                    Enew[j] = ErlangB(rhoj_1, c[j]);
                }
                return new Da_fpi.SweepResult<double[]>(Enew, E_1);
            }
        }, E0, fpopts);
        double[] E = fpres.x;
        int niter = fpres.it;

        double[] QLen = nu.clone();
        for (int r = 0; r < R; r++) {
            for (int j = 0; j < J; j++) {
                QLen[r] *= FastMath.pow(1 - E[j], A[j][r]);
            }
        }
        for (int i = 0; i < QLen.length; i++) {
            QLen[i] = FastMath.max(QLen[i], 0.0);
        }

        double[] Loss = new double[R];
        for (int i = 0; i < R; i++) {
            Loss[i] = 1.0 - QLen[i] / nu[i];   // blocking = 1 - carried/offered
        }

        return new Ret.lossnErlangFP(new Matrix(QLen), new Matrix(Loss), new Matrix(E), niter);
    }

    private static double ErlangB(double nu, double C) {
        double den = 0.0;
        int i = 0;
        while (i <= C) {
            den += FastMath.exp(i * FastMath.log(nu) - Maths.factln(i));
            i++;
        }
        double blockProb = C * FastMath.log(nu) - Maths.factln(C) - FastMath.log(den);
        return FastMath.exp(blockProb);
    }
}
