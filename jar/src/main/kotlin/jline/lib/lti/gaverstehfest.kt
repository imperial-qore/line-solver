/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti

import jline.lib.lti.abatewhitt.binom
import jline.lib.lti.abatewhitt.factorial
import org.apache.commons.math3.util.FastMath
import java.math.BigDecimal
import java.math.RoundingMode
import java.util.*
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.pow

object gaverstehfest {
    fun getomega(n: Int): Array<BigDecimal?> {
        var n = n
        if (n % 2 == 1) {
            n-- // for gaver-stehfest, as it only supports even n
        }
        val res = arrayOfNulls<BigDecimal>(n)
        for (k in 1..n) {
            var `val` = BigDecimal("1")
            `val` = `val`.multiply(BigDecimal.valueOf(((-1.0).pow((n / 2.0) + k))))
            `val` = `val`.multiply(BigDecimal.valueOf(ln(2.0)))
            // summation
            var sum_val = BigDecimal("0")
            var j = FastMath.floor((k + 1) / 2.0).toInt()
            while (j <= FastMath.min(k.toDouble(), n / 2.0)) {
                var val2 = BigDecimal.valueOf(j.toDouble().pow(n / 2.0 + 1))
                val2 = val2.divide(BigDecimal(factorial(n / 2)),
                    9999,
                    RoundingMode.HALF_EVEN) // https://stackoverflow.com/questions/4591206/arithmeticexception-non-terminating-decimal-expansion-no-exact-representable
                val2 = val2.multiply(BigDecimal(binom(n / 2, j)))
                val2 = val2.multiply(BigDecimal(binom(2 * j, j)))
                val2 = val2.multiply(BigDecimal(binom(j, k - j)))
                sum_val = sum_val.add(val2)
                j++
            }
            `val` = `val`.multiply(sum_val)
            res[k - 1] = `val`
        }
        return res
    }

    fun getalpha(n: Int): Array<BigDecimal?> {
        var n = n
        if (n % 2 == 1) n-- // gaver-schefest only supports even n; same workaround done in inverselaplace.org [5].

        val res = arrayOfNulls<BigDecimal>(n)
        for (k in 1..n) {
            res[k - 1] = BigDecimal.valueOf(k * FastMath.log(2.0))
        }
        return res
    }

    fun laplace_result(d: Double): Double {
        // function = 1/(1 + s^2)
        //    return 1d / (1d + d * d);
        return (1.0 / d) / (exp(d) - 1)
    }

    @JvmStatic
    fun main(args: Array<String>) {
        val alpha = getalpha(12)
        val omega = getomega(12)
        val res = TreeMap<Double, BigDecimal>()
        var d = 0.1
        while (d <= 15) {
            var sum_val = BigDecimal.valueOf(0)
            for (i in alpha.indices) {
                sum_val = sum_val.add(omega[i]!!.multiply(BigDecimal.valueOf(laplace_result(alpha[i].toString()
                    .toDouble() / d))))
            }
            val r = sum_val.divide(BigDecimal.valueOf(d), 9999, RoundingMode.HALF_EVEN)
            res[d] = r
            d += 0.05
        }
        // just for verification
        for (d in res.keys) {
            println(d.toString() + " " + res[d]!!.toDouble())
        }
    }
}
