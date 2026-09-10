/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;

public final class talbot {
    private talbot() {}

    // represent number as Pair<Double, Double>(real, complex)
    // pi = FastMath.PI
    public static ArrayList<Complex> getalpha(int n) {
        ArrayList<Complex> arr = new ArrayList<Complex>();
        // for k = 1
        arr.add(new Complex(2.0 * n / 5.0, 0.0));
        // for k = 2 onwards
        for (int i = 2; i <= n; i++) {
            // get real and imaginary part
            double real_part = 2 * (i - 1) * FastMath.PI / 5 * (1.0 / FastMath.tan((i - 1) * FastMath.PI / n));
            double imaginary_part = 2 * (i - 1) * FastMath.PI / 5;
            arr.add(new Complex(real_part, imaginary_part));
        }
        return arr;
    }

    public static ArrayList<Complex> getomega(int n, ArrayList<Complex> alpha) {
        ArrayList<Complex> arr = new ArrayList<Complex>();
        // for k = 1
        arr.add(alpha.get(0).exp().divide(5.0));
        //  arr.add(new Complex(Math.exp(alpha.get(0).getReal())/5, 0d));
        // for k = 2 onwards
        for (int i = 2; i <= n; i++) {
            Complex current_alpha = alpha.get(i - 1).exp();
            double temp_var = (i - 1) * FastMath.PI / n;
            Complex result = current_alpha.multiply(2).divide(5.0).multiply(
                    new Complex(1.0,
                            temp_var * (1 + 1.0 / FastMath.pow(Math.tan(temp_var), 2)) - 1.0 / FastMath.tan(temp_var)));
            arr.add(result);
        }
        return arr;
    }
}
