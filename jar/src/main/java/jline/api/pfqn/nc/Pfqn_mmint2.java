/**
 * @file McKenna-Mitra (MMINT2) two-station integral via Simpson's rule
 *
 * Implements numerical integration for computing normalizing constants in multi-class
 * repairman models using Simpson's rule integration.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.util.FastMath;

public final class Pfqn_mmint2 {
    private Pfqn_mmint2() {}

    public static Ret.pfqnNc pfqn_mmint2(final Matrix L, final Matrix N, final Matrix Z) {
        final int[] nnzClasses = N.findNonNegative().toIntArray1D();
        int order = 12;

        UnivariateFunction func = new UnivariateFunction() {
            @Override
            public double value(double u) {
                double expTerm = FastMath.exp(-u);
                double prodTerm = 1.0;
                for (int j : nnzClasses) {
                    double term = Z.get(j) + L.get(j) * u;
                    prodTerm *= FastMath.pow(term, N.get(j));
                }
                return expTerm * prodTerm;
            }
        };

        double p = 1 - FastMath.pow(10.0, -(double) order);
        double exp1prctile = -FastMath.log(1 - p);

        BaseAbstractUnivariateIntegrator integrator = new SimpsonIntegrator(1e-12, 1e-8, 3, 64);

        double integralValue;
        try {
            integralValue = integrator.integrate(Integer.MAX_VALUE, func, 0.0, exp1prctile);
        } catch (Exception e) {
            throw new RuntimeException("Integration failed", e);
        }

        Matrix Nmat = new Matrix(N);
        double lG = FastMath.log(integralValue) - Nmat.factln().elementSum();
        double G = FastMath.exp(lG);

        return new Ret.pfqnNc(G, lG);
    }
}
