/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti

import org.apache.commons.math3.complex.Complex
import org.apache.commons.math3.util.FastMath
import kotlin.math.tan

object talbot {
    // represent number as Pair<Double, Double>(real, complex)
    // pi = FastMath.PI
    fun getalpha(n: Int): ArrayList<Complex> {
        val arr = ArrayList<Complex>()
        // for k = 1
        arr.add(Complex(2.0 * n / 5.0, 0.0))
        // for k = 2 onwards
        for (i in 2..n) {
            // get real and imaginary part
            val real_part = 2 * (i - 1) * FastMath.PI / 5 * (1.0 / FastMath.tan((i - 1) * FastMath.PI / n))
            val imaginary_part = 2 * (i - 1) * FastMath.PI / 5
            arr.add(Complex(real_part, imaginary_part))
        }
        return arr
    }

    fun getomega(n: Int, alpha: ArrayList<Complex>): ArrayList<Complex> {
        val arr = ArrayList<Complex>()
        // for k = 1
        arr.add((alpha[0].exp()).divide(5.0))
        //  arr.add(new Complex(Math.exp(alpha.get(0).getReal())/5, 0d));
        // for k = 2 onwards
        for (i in 2..n) {
            val current_alpha = alpha[i - 1].exp()
            val temp_var = (i - 1) * FastMath.PI / n
            val result = current_alpha.multiply(2).divide(5.0).multiply(Complex(1.0,
                temp_var * (1 + 1.0 / FastMath.pow(tan(temp_var), 2)) - 1.0 / FastMath.tan(temp_var)))
            arr.add(result)
        }
        return arr
    } //    public static void main(String[] args) {
    //        ArrayList<Complex> alpha = getalpha(92);
    //        ArrayList<Complex> omega = getomega(92, alpha);
    //        abatewhitt.getResult(alpha, omega);
    // //       System.out.println(res);
    //    }
}
