/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti

import org.apache.commons.math3.analysis.FunctionUtils
import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.util.Pair
import org.apfloat.*
import java.util.*
import kotlin.math.cos
import kotlin.math.pow

object cme {
    var ul: IterativeLegendreGaussIntegrator = IterativeLegendreGaussIntegrator(50, 1e-9, 1e-9)

    fun moments(`fun`: UnivariateFunction?, k: Int): Apfloat {
        //     UnivariateIntegrator integrator = new RombergIntegrator(1e-9, 1e-9, 3, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        var `fun` = `fun`
        val tomultiply = (UnivariateFunction { t: Double -> FastMath.pow(t, k) })
        `fun` = FunctionUtils.multiply(`fun`, tomultiply)

        //    System.out.println(fun);
        //  return (laguerre.Laguerre(fun, 5000));
        //    return (customromberg.starter(fun, 0, 99));

        // return BigDecimal.valueOf(integrator.integrate(900000000, fun, 0, 99));
        // Legendre-Gaussian only
        // return BigDecimal.valueOf(ul.integrate(90000000, gaussiancombine(fun, 0, 99), -1, 1));
        return Apfloat((ul.integrate(900000000, `fun`, 0.0, 99.0)), abatewhitt.precision.toLong())
    }

    //    public static UnivariateFunction gaussiancombine(UnivariateFunction parent, double lb, double ub) {
    //        // parent(child))
    //        double eta = (ub - lb) / 2;
    //        double mu = (ub + lb) / 2;
    //        UnivariateFunction child = (t -> eta * (eta * t + mu));
    //        return FunctionUtils.compose(parent, child);
    //    }
    fun scv(`fun`: UnivariateFunction?): Apfloat {
        val res = moments(`fun`, 2).multiply(moments(`fun`, 0)).divide(ApfloatMath.pow(moments(`fun`, 1), 2))
        return res.subtract(Apfloat.ONE)
    }

    fun <T> deepcopy(tocopy: ArrayList<T>): ArrayList<T> {
        val result = ArrayList<T>()
        for (d in tocopy) result.add(d)
        return result
    }

    fun getnormalrandom(n: Int): ArrayList<Double> {
        // returns n randomly distributed numbers
        val rand = Random()
        val res = ArrayList<Double>()
        for (i in 0..<n) res.add(rand.nextGaussian())
        return res
    }

    fun getmefunction(params: ArrayList<Double>): UnivariateFunction {
        var temp =
            UnivariateFunction { t: Double -> FastMath.exp(-t) * FastMath.pow(cos(params[0] * t - params[1]), 2) }
        for (i in 0..<params.size - 2) {
            val finalI = i
            val temp2 = UnivariateFunction { t: Double -> FastMath.pow(cos(params[0] * t - params[finalI + 2]), 2) }
            temp = FunctionUtils.multiply(temp, temp2)
        }
        return temp
    }

    // only portion not under high precision
    // code directly adapted from http://webspn.hit.bme.hu/~illes/mincvnum.zip (as part of the "Concentrated matrix exponential distributions" paper)
    fun rechenberg(x: ArrayList<Double>): ArrayList<Double> {
        ArrayList<Double>()
        var opt = deepcopy(x)
        var sigma = 0.1
        val c = 0.9
        var siker = 0
        val n = x.size
        var best = deepcopy(opt)
        for (i in 1..1000) {
            val saveme = opt[0]
            val rands = getnormalrandom(n)
            val y = ArrayList<Double>()
            for (j in 0..<n) y.add(opt[j] + sigma * rands[j])
            if (y[0] <= 0) y[0] = saveme
            val mef1 = getmefunction(y) // in preparation for f(y) in rechenberg.m
            var mef2 = getmefunction(opt) // as above
            if (scv(mef1).compareTo(scv(mef2)) < 0) { // to prevent cases of negative SCV escaping
                siker++
                opt = deepcopy(y)
            }
            if (i % 20 == 0) {
                if (siker < 4) sigma *= c
                else if (siker > 4) sigma /= c
                siker = 0
            }
            mef2 = getmefunction(opt)
            val mef3 = getmefunction(best)
            if (scv(mef2).compareTo(scv(mef3)) < 0) {
                best = deepcopy(opt)
            }
            //   scv_opts.add(scv(getmefunction(opt)));
        }
        //     System.out.println("SCV opts");
        //    System.out.println(scv_opts);
        //    double arr = scv(getmefunction(opt));
        return opt
    }

    fun binaryCheck(num: Long, polysize: Int): BooleanArray {
        val res = BooleanArray(polysize)
        res[0] = true
        val bin = StringBuilder(java.lang.Long.toBinaryString(num))
        while (bin.length < polysize - 1) {
            bin.insert(0, "0")
        }
        bin.insert(0, "1")
        for (i in res.indices) {
            if (bin[i] == '1') res[i] = true
        }
        return res
    }

    fun <T> convert_to_arraylist(input: Array<T>): ArrayList<T> {
        return ArrayList(Arrays.asList(*input))
    }

    fun convert_to_laplace(res: ArrayList<Double>): Pair<ArrayList<Apcomplex?>, ArrayList<Apcomplex>> {
        // res[0] = omega, the rest are phi variables
        val omega = Apfloat(res[0], abatewhitt.precision.toLong())/*
        Given problem is
        (pi_{i = 1}^n cos (omega - phi_i))^2
        Then, use fact that (cos A cos B) = (1/2 (cos (A + B) + cos (A - B))) so that we get something of the form
        cos (A +- B +- C +- ...) + summed all over
        After that, the goal is to square this result. This will get a O((2^(n - 1))^2) sum.
        Then, use the fact that 2 cos(P) cos(Q) = (cos(P + Q) + cos(P - Q)). Where P = Q, this will naturally resolve to 1/2 (cos 2P + 1)
        At this point, if the result is real, we can mark it as such and sum them all together; if not, that's Apcomplex and should be resolved together.
        For the Apcomplex variables, we now have the form W = e^-t cos R, where R has both a real and an Apcomplex aspect. Note that e^(x + iy) = e^x (cos y + i sin y).
        Then, W can be split as W = 1/2 (U + V), where U = e^-t (cos R + i sin R) and V = e^-t (cos R - i sin R). Note that this is the same as -R in U
        After that, U and V can be written as e^(-t + R) and e^(-t - R), respectively. This is the same as e^-t(e^R + e^(-R)). Finally, take the t part outside.
         */
        // STEP 1: Get the polynomial cos(P) where P = omega*t - phi_i
        val polycoeff = ArrayList<Array<Apfloat>>()
        for (i in 1..<res.size) {
            polycoeff.add(arrayOf(Apfloat(-1 * res[i], abatewhitt.precision.toLong()), omega)) // -phi_i + omega_t
        }
        // cos (B + AX) where A = omega and B = phi
        // get part 2: this would be 2^{n - 1} polynomials of the form A +- B +- ...
        val polypart2 = ArrayList<Array<Apfloat?>>()
        val to_divide =
            ApfloatMath.pow(Apfloat(4, abatewhitt.precision.toLong()), Apfloat((polycoeff.size - 1).toLong()))
        var constant_term = Apfloat.ONE.divide(to_divide)
        // remember that the PDFs are unnormalised, so we will need to normalise
        // c int(f(t)) = 1
        val int_value = moments(getmefunction(res), 0)
        println(int_value)
        constant_term = constant_term.divide(int_value)
        //  constant_term = constant_term * (1.0 / int_value);
        for (l in 0..<(2.0.pow((polycoeff.size - 1).toDouble())).toLong()) {
            val b = binaryCheck(l, polycoeff.size)
            // loop and add as appropriate
            val newpoly = arrayOfNulls<Apfloat>(2)
            newpoly[0] = Apfloat(0, abatewhitt.precision.toLong())
            newpoly[1] = Apfloat(0, abatewhitt.precision.toLong())
            for (i in b.indices) {
                if (b[i]) {
                    newpoly[0] = newpoly[0]!!.add(polycoeff[i][0])
                    newpoly[1] = newpoly[1]!!.add(polycoeff[i][1])
                } else {
                    newpoly[0] = newpoly[0]!!.subtract(polycoeff[i][0])
                    newpoly[1] = newpoly[1]!!.subtract(polycoeff[i][1])
                }
            }
            polypart2.add(newpoly)
        }
        // part 3: do the squaring
        // also part 4: cos^2 x = (1/2 (cos 2x + 1)) -> constant term is 1/2
        val polypart3 = ArrayList<Array<Apfloat?>>(polypart2.size * polypart2.size + polypart2.size)
        val hs = HashSet<ArrayList<Apfloat?>>()
        for (i in polypart2.indices) {
            for (j in i..<polypart2.size) {
                // if (i != j), then 2 cos P cos Q = cos(P + Q) + cos (P - Q)
                if (i != j) {
                    var res2 = arrayOfNulls<Apfloat>(2)
                    res2[0] = polypart2[i][0]!!.add(polypart2[j][0])
                    res2[1] = polypart2[i][1]!!.add(polypart2[j][1])
                    polypart3.add(res2)
                    // (cos (P - Q))
                    res2 = arrayOfNulls(2)
                    res2[0] = polypart2[i][0]!!.subtract(polypart2[j][0])
                    res2[1] = polypart2[i][1]!!.subtract(polypart2[j][1])
                    polypart3.add(res2)
                } else {
                    // cos^2 P = 1/2 (cos 2P + 1).
                    // handle 1/2 term, looks like we will need to do it via a hashset
                    // this is because the entire structure is based on polynomials till date
                    var res2 = arrayOfNulls<Apfloat>(2)
                    res2[0] = Apfloat(2).multiply(polypart2[i][0])
                    res2[1] = Apfloat(2).multiply(polypart2[i][1])
                    polypart3.add(res2)
                    // add to hashset so that we can handle this specifically in the end
                    hs.add(convert_to_arraylist(res2))
                    // also handle the cos (0) part
                    res2 = arrayOfNulls(2)
                    res2[0] = Apfloat.ZERO
                    res2[1] = Apfloat.ZERO
                    polypart3.add(res2)
                    hs.add(convert_to_arraylist(res2))
                }
            }
        }
        //   constant_term /= 2;
        /* now, for all the polypart3 cases, check whether it is dependent on $t$.
        if it is, get a conjugate and finalise in that way [real, Apcomplex] + t*[real, Apcomplex]
        hint: e^-t (cos P) = 1/2 e^-t (cos P + i sin P) + 1/2 e^-t (cos P - i sin P)
        = 1/2 e^(-t + iP) + 1/2 e^(-t - iP)
        */
        var ptr = 0
        val eta_alphacombo = HashMap<Apcomplex, Apcomplex>()
        while (ptr < polypart3.size) {
            // check if it even depends on t
            val doubles = polypart3[ptr++]
            if (doubles[1] == Apfloat.ZERO) {
                // only real! alpha_i is easily 1
                // these terms can be directly evaluated
                //    alpha.add(Apcomplex.ONE);
                var C = Apcomplex((ApfloatMath.cos(doubles[0])), Apfloat.ZERO)
                if (hs.contains(convert_to_arraylist<Apfloat?>(doubles))) {
                    // needs to be divided by 2 as before
                    C = C.divide(Apint(2))
                }
                C = C.multiply(constant_term)
                //       eta.add(C);
                if (eta_alphacombo.containsKey(Apcomplex.ONE)) {
                    eta_alphacombo[Apcomplex.ONE] = eta_alphacombo[Apcomplex.ONE]!!.add(C)
                } else {
                    eta_alphacombo[Apcomplex.ONE] = C
                }
            } else {
                // computing eta1
                var C1 = Apcomplex(Apfloat(0, abatewhitt.precision.toLong()), doubles[0]) //
                C1 = ApcomplexMath.exp(C1)
                C1 = C1.divide(Apint(2))
                // computing alpha1
                val A1 = Apcomplex(Apfloat.ONE, Apfloat(-1).multiply(doubles[1])) // {1, -i omega}
                // eta2
                var C2 = Apcomplex(Apfloat(0, abatewhitt.precision.toLong()), Apfloat(-1).multiply(doubles[0])) //
                C2 = ApcomplexMath.exp(C2)
                C2 = C2.divide(Apint(2))
                // alpha2
                val A2 = Apcomplex(Apfloat.ONE, doubles[1]) // {1, i omega}
                if (hs.contains(convert_to_arraylist<Apfloat?>(doubles))) {
                    // the result MUST BE DIVIDED BY 2
                    // (only for etas though!)
                    C1 = C1.divide(Apint(2))
                    C2 = C2.divide(Apint(2))
                }
                // divide C1 and C2 with old normalising constant
                C1 = C1.multiply(constant_term)
                C2 = C2.multiply(constant_term)
                // and add it to alpha and eta
                //   eta.add(C1);
                //  alpha.add(A1);
                if (eta_alphacombo.containsKey(A1)) {
                    eta_alphacombo[A1] = eta_alphacombo[A1]!!.add(C1)
                } else {
                    eta_alphacombo[A1] = C1
                }
                //  eta.add(C2);
                //   alpha.add(A2);
                if (eta_alphacombo.containsKey(A2)) {
                    eta_alphacombo[A2] = eta_alphacombo[A2]!!.add(C2)
                } else {
                    eta_alphacombo[A2] = C2
                }
            }
        }
        val alpha = ArrayList<Apcomplex>()
        val eta = ArrayList<Apcomplex?>()
        for (e in eta_alphacombo.keys) {
            alpha.add(e)
            eta.add(eta_alphacombo[e])
        }
        return Pair(eta, alpha)
    }

}