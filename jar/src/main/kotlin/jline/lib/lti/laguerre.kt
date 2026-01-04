/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti

import jline.lib.lti.abatewhitt.binom
import jline.lib.lti.abatewhitt.factorial
import org.apache.commons.math3.analysis.FunctionUtils
import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction
import org.apache.commons.math3.analysis.solvers.LaguerreSolver
import java.math.BigDecimal
import java.math.RoundingMode
import kotlin.math.exp

object laguerre {
    // for performance optimisation
    var coeff: DoubleArray? = null
    var lroots: DoubleArray? = null
    var coeffnplus1: DoubleArray? = null
    var weights: Array<BigDecimal?>? = null

    fun Laguerre(f: UnivariateFunction, n: Int): BigDecimal {
        // get rid of the exponent since Laguerre already takes that into account
        var f = f
        val temp = UnivariateFunction { a: Double -> exp(a) }
        f = FunctionUtils.multiply(f, temp)
        // get Laguerre coefficients
        if (coeff == null) coeff = getLaguerreCoefficients(n)
        // get the roots
        if (lroots == null) {
            lroots = getLaguerreRoots(coeff)
        }
        // and weights
        if (weights == null) weights = getweight(lroots!!, n)
        var toadd = BigDecimal.valueOf(0)
        for (i in 0..<n) {
            toadd = toadd.add(weights!![i]!!.multiply(BigDecimal.valueOf(f.value(lroots!![i]))))
        }
        return toadd
    }

    fun getLaguerreRoots(coeff: DoubleArray?): DoubleArray {
        val L = LaguerreSolver()
        val croots = L.solveAllComplex(coeff, 0.0)
        val roots = DoubleArray(croots.size)
        for (i in roots.indices) roots[i] = croots[i].real
        return roots
    }

    fun getweight(roots: DoubleArray, n: Int): Array<BigDecimal?> {
        val res = arrayOfNulls<BigDecimal>(roots.size)
        if (coeffnplus1 == null) coeffnplus1 = getLaguerreCoefficients(n + 1)
        val pf = PolynomialFunction(coeffnplus1)
        //   System.out.println("pf = " + pf);
        for (i in roots.indices) {
            var temp = BigDecimal.valueOf(pf.value(roots[i])).pow(2)
            val temp2 = BigDecimal.valueOf((n + 1).toLong()).pow(2)
            temp = temp.multiply(temp2)
            // res[i] = BigDecimal.valueOf(roots[i]/(Math.pow(pf.value(roots[i])*(n + 1), 2)));
            res[i] = BigDecimal.valueOf(roots[i]).divide(temp, 2500, RoundingMode.HALF_EVEN)
            //   res[i] = roots[i]/Math.pow((n + 1)*pf.value(roots[i]), 2);
        }
        return res
    }

    fun getLaguerreCoefficients(n: Int): DoubleArray {
        // generalised closed-form
        val res = DoubleArray(n + 1)
        for (i in 0..n) {
            var temp = BigDecimal("-1")
            temp = temp.pow(i)
            temp = temp.multiply(BigDecimal(binom(n, i)))
            temp = temp.divide(BigDecimal(factorial(i)), 2500, RoundingMode.HALF_EVEN)
            res[i] = temp.toDouble()
        }
        return res
    }
}
