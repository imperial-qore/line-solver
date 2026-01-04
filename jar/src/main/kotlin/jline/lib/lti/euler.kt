/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti

import jline.lib.lti.abatewhitt.binom
import org.apache.commons.math3.complex.Complex
import org.apache.commons.math3.util.FastMath
import kotlin.math.pow

object euler {
    fun getalpha(n: Int): ArrayList<Complex> {
        val result = ArrayList<Complex>()
        for (i in 0..<n) {
            result.add(Complex((n - 1) * FastMath.log(10.0) / 6, FastMath.PI * i))
        }
        return result
    }

    fun geteta(n: Int): DoubleArray {
        val res = DoubleArray(n)
        res[0] = 0.5
        // euler defined only for odd n!
        for (i in 2..(n + 1) / 2) {
            res[i - 1] = 1.0
        }
        res[n - 1] = 1.0 / FastMath.pow(2.0, (n - 1) / 2.0)
        for (i in 1..<(n - 1) / 2) {
            res[n - i - 1] = res[n - i] + FastMath.pow(2.0, (1 - n) / 2.0) * binom((n - 1) / 2, i).toString().toDouble()
        }
        return res
    }

    fun getomega(n: Int): ArrayList<Complex> {/*
         * Page 6 of Numerical inverse Laplace transformation using concentrated matrix exponential distributions contains an error
         * specifically, that eta starts off negative, when it should be positive as k starts from 1 there
         * This can be validated by going to equation 36 of http://www.columbia.edu/~ww2040/AbateUnified2006.pdf - there k starts from 0 and hence (-1)^k is positive
         * the difference is that otherwise we would get -f(x) instead of f(x) as the result
         */
        val eta = geteta(n)
        val res = ArrayList<Complex>()
        for (i in 1..n) {
            res.add(Complex(10.0.pow((n - 1) / 6.0) * FastMath.pow(-1.0, i - 1) * eta[i - 1]))
        }
        return res
    }

    @JvmStatic
    fun main(args: Array<String>) {
        getalpha(99)
        getomega(99)
        //   abatewhitt.getResult(alpha, omega);
    }
}