/**
 * @file Grundmann-Möller cubature method for normalizing constant computation
 *
 * Implements the cubature (multi-dimensional integration) approach for computing normalizing
 * constants in closed product-form queueing networks. Uses Grundmann-Möller simplex quadrature
 * with configurable order for exact or approximate computation.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_cub {
    private Pfqn_cub() {}

    /** quadrature points of v in the think-time branch; pfqn_nc prices CUB against it */
    public static final int CUB_THINK_STEPS = 10000;

    /** integrand-evaluation budget above which pfqn_nc prefers le over cub */
    public static final double CUB_MAX_EVALS = 1e7;

    /**
     * Number of integrand evaluations pfqn_cub performs at this order. The
     * Grundmann-Moeller rule of degree order on the (M-1)-simplex evaluates
     * sum_{d=0..order} binom(M-1+2d, M-1) points, and a non-zero think time
     * repeats the whole rule at every v-quadrature step.
     *
     * @param M     - number of queueing stations
     * @param order - cubature order
     * @param Z     - think time per class, summed over delay stations
     * @return the integrand-evaluation count
     */
    public static double pfqn_cub_evals(int M, int order, Matrix Z) {
        int n = M - 1;
        double nodes = 0.0;
        for (int d = 0; d <= order; d++) {
            nodes += Maths.binomialCoeff(n + 2 * d, n);
        }
        boolean hasThink = Z != null && Z.elementSum() >= GlobalConstants.FineTol;
        return hasThink ? nodes * CUB_THINK_STEPS : nodes;
    }

    /**
     * Cubature method to compute the normalizing constant of a load-independent closed queueing network model
     */
    public static Ret.pfqnNc pfqn_cub(Matrix L, Matrix N, Matrix Z) {
        int order = (int) FastMath.ceil((N.elementSum() - 1) / 2.0);
        double atol = 1e-8;
        return pfqn_cub(L, N, Z, order, atol);
    }

    /**
     * Cubature method to compute the normalizing constant of a load-independent closed queueing network model
     *
     * @param L     - demands at all stations
     * @param N     - number of jobs for each class
     * @param Z     - think time for each class
     * @param order - cubature order ( ceil((|N|-1)/2) or above is exact)
     * @param atol  - numerical tolerance for simplex quadrature
     * @return normalizing constant, its logarithm, and mean performance metrics computed as a by-product
     */
    public static Ret.pfqnNc pfqn_cub(Matrix L, Matrix N, Matrix Z, int order, double atol) {
        int M = L.getNumRows();

        if (L.isEmpty() || N.isEmpty() || N.elementSum() == 0.0) {
            return new Ret.pfqnNc(1.0, 0.0);
        }

        if (order < 0) {
            order = (int) FastMath.ceil((N.elementSum() - 1) / 2.0);
        }

        if (atol <= 0) {
            atol = 1e-8;
        }

        if (Z.isEmpty() || Z.elementSum() < atol) {
            double Nt = N.elementSum();
            Matrix beta = N.scale(1.0 / Nt);

            Maths.simplexQuadResult simplexResult = Maths.simplexquad(new Ret.pfqnCUB(L, Nt, beta), M - 1, order, atol);
            double[] Q = simplexResult.Q;
            double Gn = Q[Q.length - 1] * FastMath.exp(Maths.factln(N.elementSum() + M - 1) - N.factln().elementSum());
            return new Ret.pfqnNc(Gn, FastMath.log(Gn));
        } else {
            int steps = CUB_THINK_STEPS;
            double Nt = N.elementSum();
            Matrix beta = N.scale(1 / Nt);
            double Gn = 0.0;
            double vmax = Nt * 10;
            double dv = vmax / steps;

            double v = 0.0;
            while (v <= vmax) {
                Matrix Lv = L.scale(v).add(Z.repmat(M, 1));
                Maths.simplexQuadResult simplexResult = Maths.simplexquad(new Ret.pfqnCUB(Lv, Nt, beta), M - 1, order, atol);
                double[] Q = simplexResult.Q;
                double dG = FastMath.exp(-v) * FastMath.pow(v, M - 1) * Q[Q.length - 1] * dv;
                Gn += dG;

                if (v > 0 && dG / Gn < atol) {
                    break;
                }
                v += dv;
            }

            Gn *= FastMath.exp(-N.factln().elementSum());
            return new Ret.pfqnNc(Gn, FastMath.log(Gn));
        }
    }
}
