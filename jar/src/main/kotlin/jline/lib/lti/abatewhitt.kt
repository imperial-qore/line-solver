/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti

import org.apache.commons.math3.util.FastMath
import org.apfloat.Apcomplex
import org.apfloat.Apfloat
import org.apfloat.ApfloatMath
import org.apfloat.Apint
import java.math.BigInteger
import java.util.*
import java.util.function.UnaryOperator

object abatewhitt {
    const val precision: Int = 35

    @JvmField
    var to_int: UnaryOperator<Apcomplex>? = null
    var factorial_memo: HashMap<Int, BigInteger> = HashMap()

    fun laplace_result(d: Apcomplex): Apcomplex {
        /**
         * If you want to evaluate a specific function, override to_int. For example,
         * to_int = (t -> (Apcomplex.ONE.divide(t)).divide(ApcomplexMath.exp(t).add(Apcomplex.ONE))); is (1/t)/(exp(t) + 1)
         * to_int = (t -> Apcomplex.ONE.divide(Apcomplex.ONE.add(d.multiply(d)))); is 1/(1 + s^2) (sine function)
         */
        return to_int!!.apply(d)
    }

    fun optimiser_laplace(d: Apcomplex): Apcomplex {
        /**
         * This is used for optimising lambda
         */
        return Apcomplex.ONE.divide(Apcomplex.ONE.add(d.multiply(d)))
    }

    fun optimise_lambda(alpha: ArrayList<Apcomplex>, omega: ArrayList<Apcomplex>): Apfloat {
        var lambda: Apfloat = Apfloat.ONE
        var bound_l = Apfloat(0.78, precision.toLong())
        var bound_r = Apfloat(1, precision.toLong())
        var bound = (bound_l.add(bound_r)).divide(Apfloat(2))
        var max_lambda: Apfloat = Apfloat.ONE
        val lambdavals = HashMap<Apfloat, Apfloat>()
        while (true) {
            bound = (bound_l.add(bound_r)).divide(Apint(2))
            val d = ApfloatMath.pi(precision.toLong()).divide(Apint(2)) // pi/2
            var ans = Apfloat(1, precision.toLong()).divide(d)
            var partialsum: Apcomplex = Apfloat(0, precision.toLong())
            for (i in alpha.indices) {
                val scaled_alpha = alpha[i].multiply(lambda)
                partialsum = partialsum.add(omega[i].multiply(optimiser_laplace(scaled_alpha.divide(d))))
            }
            ans = ans.multiply(partialsum.real())
            if ((ans.multiply(lambda)).compareTo(bound) >= 0) {
                //       System.out.println(bound + " " + bound_l + " " + bound_r + " ");
                bound_l = bound
                max_lambda = lambda
                lambda = Apfloat.ONE
                if ((bound_r.subtract(bound_l)).compareTo(Apfloat(0.0001, precision.toLong())) < 0) {
                    break
                }
            } else if (lambda.compareTo(Apfloat(9)) > 0) {
                lambda = Apfloat.ONE
                //     System.out.println(bound + " " + bound_l + " " + bound_r + " ");
                bound_r = bound
                if ((bound_r.subtract(bound_l)).compareTo(Apfloat(0.0001, precision.toLong())) < 0) {
                    println(lambdavals)
                }
            } else lambda = lambda.add(Apfloat(0.01, precision.toLong()))
            if (!lambdavals.containsKey(lambda)) lambdavals[lambda] = ans.multiply(lambda)
        }
        //   System.out.println(lambdavals);
        return max_lambda
    }

    fun getResult(alpha: ArrayList<Apcomplex>, omega: ArrayList<Apcomplex>) {
        val res = TreeMap<Apfloat, Apfloat>()
        val to_shift = 1 // determine if shifting is necessary
        var max_val = 0.0
        var unshifted_max = 0.0
        println("Entered getResult")
        println("Size of alpha " + alpha.size)
        val lambda = optimise_lambda(alpha, omega)
        println("optimised lambda $lambda")
        val shifter = Apcomplex(Apfloat(-1, precision.toLong()).divide(Apfloat(1, precision.toLong())), Apfloat.ZERO)
        //       System.out.println("Printing results");
        var d = Apfloat(0.01, precision.toLong())
        while (d.compareTo(Apfloat(70.01, precision.toLong())) <= 0) {
            //    BigDecimal ans = BigMath.ONE.divide(d, precision, RoundingMode.HALF_EVEN);
            var ans = Apfloat.ONE.divide(d)
            var ans_unshifted = Apfloat.ONE.divide(d)
            var partialsum: Apcomplex = Apcomplex.ZERO
            var partialsum_unshifted: Apcomplex = Apcomplex.ZERO
            for (i in alpha.indices) {
                if (to_shift == 0) {
                    val scaled_alpha = alpha[i].multiply(lambda)
                    partialsum = partialsum.add(omega[i].multiply(laplace_result(scaled_alpha.divide(d))))
                } else {
                    partialsum_unshifted =
                        partialsum_unshifted.add(omega[i].multiply(laplace_result(alpha[i].divide(d))))
                    partialsum = partialsum.add(omega[i].multiply(laplace_result((alpha[i].divide(d)).add(shifter))))
                }
            }
            if (to_shift == 0) {
                ans = ans.multiply((partialsum.real())).multiply(lambda)
                res[d] = ans
                println(d.toString() + " " + res[d]!!.toDouble())
            } else {
                ans = ans.multiply((partialsum.real()))
                ans_unshifted = ans_unshifted.multiply((partialsum_unshifted.real()))
                //     ans = ApfloatMath.exp(ans.multiply(shifter).real());
                ans = ans.multiply(ApfloatMath.exp(shifter.real().multiply(d)))
                max_val = FastMath.max(ans.toDouble(), max_val)
                unshifted_max = FastMath.max(ans_unshifted.toDouble(), unshifted_max)
                res[d] = ans
            }
            d = d.add(Apfloat(0.05, precision.toLong()))
        }
        if (to_shift == 1) {
            println("Shifted max $max_val")
            println("Unshifted max $unshifted_max")
            for (d in res.keys) {
                println(d.toString() + " " + res[d]!!.toDouble() / (max_val / unshifted_max))
            }
        }
    }

    @JvmStatic
    fun factorial(n: Int): BigInteger? {
        if (factorial_memo.containsKey(n)) return factorial_memo[n]
        var `val` = BigInteger(BigInteger.ONE.toString())
        if (n == 1) return `val`
        for (i in 1..n) `val` = `val`.multiply(BigInteger(i.toString()))
        factorial_memo[n] = `val`
        return `val`
    }

    @JvmStatic
    fun binom(n: Int, k: Int): BigInteger {
        // binomial formula: n!/(n - k!)k!
        return factorial(n)!!.divide(factorial(n - k)).divide(factorial(k))
    }
}
