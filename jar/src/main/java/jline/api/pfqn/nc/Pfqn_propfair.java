/**
 * @file Proportionally fair allocation approximation for normalizing constants
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import de.xypron.jcobyla.Calcfc;
import de.xypron.jcobyla.Cobyla;
import jline.io.Ret;
import jline.GlobalConstants;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

public final class Pfqn_propfair {
    private Pfqn_propfair() {}

    /**
     * Compute the proportionally fair allocation approximation.
     *
     * @param L demands at all stations
     * @param N number of jobs for each class
     * @param Z think time for each class
     * @return normalizing constant, its logarithm, and performance metrics
     */
    public static Ret.pfqnNcXQ pfqn_propfair(Matrix L, Matrix N, final Matrix Z) {
        final int M = L.getNumRows();
        final int R = L.getNumCols();

        final double[] Nvec = N.toArray1D();
        final double[] Zvec = Z.toArray1D();
        final double[][] Lmat = L.toArray2D();

        Calcfc objFun = new Calcfc() {
            @Override
            public double compute(int n, int m, double[] x, double[] con) {
                double obj = 0.0;
                for (int r = 0; r < R; r++) {
                    obj += ((Nvec[r] - x[r] * Zvec[r]) * FastMath.log(FastMath.abs(x[r]) + GlobalConstants.FineTol));
                }
                for (int i = 0; i < M; i++) {
                    con[i] = 1.0;
                    for (int r = 0; r < R; r++) {
                        con[i] -= Lmat[i][r] * x[r];
                    }
                }
                System.arraycopy(x, 0, con, M, R);
                return -obj;
            }
        };

        double[] Xasy = new double[R];
        Cobyla.findMinimum(objFun, R, M + R, Xasy, 1.0, 1.0e-8, 0, 10000); // iprint=0: silent (no COBYLA iteration trace)

        // Z.get(0,r) IS CORRECT HERE and a `Z.get(0,` sweep should not flag it:
        // the reference indexes Z as a VECTOR too (pfqn_propfair.m:44-48 write
        // Z(r)), so it is no more tolerant of a per-delay-node Z than this is.
        // Summing the rows would MANUFACTURE a divergence rather than remove
        // one. The defect this pattern signals is taking Z straight from
        // snGetProductFormParams with no sumCols() in between; propfair has no
        // such caller.
        double lG = 0.0;
        for (int r = 0; r < R; r++) {
            double x = Xasy[r];
            lG += ((N.get(0, r) - x * Z.get(0, r)) * FastMath.log(1.0 / (x + GlobalConstants.FineTol)));
        }
        for (int r = 0; r < R; r++) {
            lG -= FastMath.log(Maths.fact(Xasy[r] * Z.get(0, r)));
        }

        double G = FastMath.exp(lG);

        Matrix Xa = new Matrix(Xasy);
        Matrix Qa = Xa.copy();
        Qa.fill(Double.NaN);

        return new Ret.pfqnNcXQ(G, lG, Xa, Qa, "propfair");
    }
}
