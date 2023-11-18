package jline.lib.ltinversion;
import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;

import java.math.BigDecimal;
import java.math.RoundingMode;

public class laguerre {
    // for performance optimisation
    static double[] coeff;
    static double[] lroots;
    static double[] coeffnplus1;
    static BigDecimal[] weights;

    public static BigDecimal Laguerre(UnivariateFunction fun, int n) {
        // get rid of the exponent since Laguerre already takes that into account
        UnivariateFunction temp = Math::exp;
        fun = FunctionUtils.multiply(fun, temp);
        // get Laguerre coefficients
        if (coeff == null)
            coeff = getLaguerreCoefficients(n);
        // get the roots
        if (lroots == null) {
            lroots = getLaguerreRoots(coeff);
        }
        // and weights
        if (weights == null)
            weights = getweight(lroots, n);
        BigDecimal toadd = BigDecimal.valueOf(0);
        for (int i = 0; i < n; i++) {
            toadd = toadd.add(weights[i].multiply(BigDecimal.valueOf(fun.value(lroots[i]))));
        }
        return toadd;
    }

    public static double[] getLaguerreRoots(double[] coeff) {
        LaguerreSolver L = new LaguerreSolver();
        Complex[] croots = L.solveAllComplex(coeff, 0);
        double[] roots = new double[croots.length];
        for (int i = 0; i < roots.length; i++)
            roots[i] = croots[i].getReal();
        return roots;
    }

    public static BigDecimal[] getweight(double[] roots, int n) {
        BigDecimal[] res = new BigDecimal[roots.length];
        if (coeffnplus1 == null)
            coeffnplus1 = getLaguerreCoefficients(n + 1);
        PolynomialFunction pf = new PolynomialFunction(coeffnplus1);
        //   System.out.println("pf = " + pf);
        for (int i = 0; i < roots.length; i++) {
            BigDecimal temp = BigDecimal.valueOf(pf.value(roots[i])).pow(2);
            BigDecimal temp2 = BigDecimal.valueOf(n + 1).pow(2);
            temp = temp.multiply(temp2);
            // res[i] = BigDecimal.valueOf(roots[i]/(Math.pow(pf.value(roots[i])*(n + 1), 2)));
            res[i] = BigDecimal.valueOf(roots[i]).divide(temp, 2500, RoundingMode.HALF_EVEN);
            //   res[i] = roots[i]/Math.pow((n + 1)*pf.value(roots[i]), 2);
        }
        return res;
    }

    public static double[] getLaguerreCoefficients(int n) {
        // generalised closed-form
        double[] res = new double[n + 1];
        for (int i = 0; i <= n; i++) {
            BigDecimal temp = new BigDecimal("-1");
            temp = temp.pow(i);
            temp = temp.multiply(new BigDecimal(abatewhitt.binom(n, i)));
            temp = temp.divide(new BigDecimal(abatewhitt.factorial(i)), 2500, RoundingMode.HALF_EVEN);
            res[i] = temp.doubleValue();
        }
        return res;
    }
}
