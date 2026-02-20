/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti

import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.util.FastMath
import java.math.BigDecimal
import java.math.RoundingMode

object customromberg {
    lateinit var dp: Array<Array<BigDecimal?>>
    var precision: Int = 35

    fun hn(a: Double, b: Double, n: Int): BigDecimal {
        val bg = BigDecimal(b - a)
        //      return bg.divide(BigDecimal.valueOf(Math.pow(2, n)), 999, RoundingMode.HALF_EVEN);
        return bg.divide(BigDecimal.valueOf(n.toLong()), precision, RoundingMode.HALF_EVEN)
    }

    fun solver(m: Int, n: Int, f: UnivariateFunction, lb: Double, ub: Double): BigDecimal? {
        if (dp[m][n] != null) return dp[m][n]
        if (m == 0 && n == 0) {
            dp[m][n] = hn(lb, ub, 1).multiply(BigDecimal.valueOf(f.value(lb) + f.value(ub)))
            return dp[m][n]
        } else if (n == 0) {
//            dp[m][n] = solver(m - 1, n, f, lb, ub).multiply(BigDecimal.valueOf(0.5));
//            BigDecimal tosum = new BigDecimal(0);
//            for (int k = 1; k <= FastMath.pow(2, n - 1); k++) {
//                tosum = tosum.add(BigDecimal.valueOf(f.value(lb + hn(lb, ub, m).multiply(BigDecimal.valueOf((2L * k - 1))).doubleValue())));
//            }
//            tosum = tosum.multiply(hn(lb, ub, m));
//            dp[m][n] = dp[m][n].add(tosum);
            // https://planetmath.org/compositetrapezoidalrule
            // https://www.math.usm.edu/lambers/mat460/fall09/lecture29.pdf
            var temp = hn(lb, ub, m).divide(BigDecimal.valueOf(2), precision, RoundingMode.HALF_EVEN)
            var toadd = BigDecimal.valueOf(f.value(lb))
            for (j in 1..m - 1) {
                toadd = toadd.add(BigDecimal.valueOf(2)
                    .multiply(BigDecimal.valueOf(f.value(lb + j * hn(lb, ub, m).toDouble()))))
            }
            toadd = toadd.add(BigDecimal.valueOf(f.value(ub)))
            temp = temp.multiply(toadd)
            dp[m][n] = temp
            return dp[m][n]
        } else {
            var temp = BigDecimal(1)
            var temp2 = BigDecimal(4).pow(n)
            temp2 = temp2.subtract(BigDecimal(1))
            temp = temp.divide(temp2, precision, RoundingMode.HALF_EVEN)
            temp = temp.multiply(solver(m, n - 1, f, lb, ub)!!.subtract(solver(m - 1, n - 1, f, lb, ub)))
            val res = solver(m, n - 1, f, lb, ub)!!.add(temp)
            dp[m][n] = res
            return dp[m][n]
        }
    }

    @JvmStatic
    fun main(args: Array<String>) {
        // for testing only
        // UnivariateFunction f = Math::exp;
        //   UnivariateFunction f = (t -> FastMath.pow(t, 3)*7 - 8*Math.pow(t, 2) - 3*t + 3);
        val f = (UnivariateFunction { t: Double -> FastMath.exp(-t) })
        //  UnivariateFunction f = (t -> FastMath.exp(-t*t));
        println(starter(f, 1.0, 2.0))
        for (i in dp.indices) {
            for (j in 0..i) {
                print(dp[i][j]!!.toDouble().toString() + " ")
            }
            println()
        }
    }

    fun starter(f: UnivariateFunction, lb: Double, ub: Double): BigDecimal? {
        val max_steps = 20
        var number_of_iterations = 0 // force a lower bound
        dp = Array(max_steps) { arrayOfNulls(max_steps) }
        val result = solver(max_steps - 1, max_steps - 1, f, lb, ub)
        // take a walk
        var prev_value = dp[0][0]
        for (i in 1..<max_steps) {
            for (j in 0..i) {
                if (dp[i][j]!!.subtract(prev_value).abs()
                        .compareTo(BigDecimal.valueOf(0.00000000000001)) < 0 && number_of_iterations >= 20) {
                    return dp[i][j]
                } else prev_value = dp[i][j]
                number_of_iterations++
            }
        }
        return result
    }
}
