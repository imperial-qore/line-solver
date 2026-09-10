/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti;

import java.math.BigDecimal;
import java.math.RoundingMode;

import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;

public final class laguerre {
    private laguerre() {}

    // for performance optimisation
    public static double[] coeff = null;
    public static double[] lroots = null;
    public static double[] coeffnplus1 = null;
    public static BigDecimal[] weights = null;

    public static BigDecimal Laguerre(UnivariateFunction f, int n) {
        // get rid of the exponent since Laguerre already takes that into account
        UnivariateFunction temp = new UnivariateFunction() {
            @Override
            public double value(double a) {
                return Math.exp(a);
            }
        };
        f = FunctionUtils.multiply(f, temp);
        if (coeff == null) coeff = getLaguerreCoefficients(n);
        if (lroots == null) {
            lroots = getLaguerreRoots(coeff);
        }
        if (weights == null) weights = getweight(lroots, n);
        BigDecimal toadd = BigDecimal.valueOf(0);
        for (int i = 0; i < n; i++) {
            toadd = toadd.add(weights[i].multiply(BigDecimal.valueOf(f.value(lroots[i]))));
        }
        return toadd;
    }

    public static double[] getLaguerreRoots(double[] coeff) {
        LaguerreSolver L = new LaguerreSolver();
        Complex[] croots = L.solveAllComplex(coeff, 0.0);
        double[] roots = new double[croots.length];
        for (int i = 0; i < roots.length; i++) roots[i] = croots[i].getReal();
        return roots;
    }

    public static BigDecimal[] getweight(double[] roots, int n) {
        BigDecimal[] res = new BigDecimal[roots.length];
        if (coeffnplus1 == null) coeffnplus1 = getLaguerreCoefficients(n + 1);
        PolynomialFunction pf = new PolynomialFunction(coeffnplus1);
        for (int i = 0; i < roots.length; i++) {
            BigDecimal temp = BigDecimal.valueOf(pf.value(roots[i])).pow(2);
            BigDecimal temp2 = BigDecimal.valueOf((long) (n + 1)).pow(2);
            temp = temp.multiply(temp2);
            res[i] = BigDecimal.valueOf(roots[i]).divide(temp, 2500, RoundingMode.HALF_EVEN);
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
