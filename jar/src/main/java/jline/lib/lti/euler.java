/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;

public final class euler {
    private euler() {}

    public static ArrayList<Complex> getalpha(int n) {
        ArrayList<Complex> result = new ArrayList<Complex>();
        for (int i = 0; i < n; i++) {
            result.add(new Complex((n - 1) * FastMath.log(10.0) / 6, FastMath.PI * i));
        }
        return result;
    }

    public static double[] geteta(int n) {
        double[] res = new double[n];
        res[0] = 0.5;
        // euler defined only for odd n!
        for (int i = 2; i <= (n + 1) / 2; i++) {
            res[i - 1] = 1.0;
        }
        res[n - 1] = 1.0 / FastMath.pow(2.0, (n - 1) / 2.0);
        for (int i = 1; i < (n - 1) / 2; i++) {
            res[n - i - 1] = res[n - i] + FastMath.pow(2.0, (1 - n) / 2.0)
                    * Double.parseDouble(abatewhitt.binom((n - 1) / 2, i).toString());
        }
        return res;
    }

    public static ArrayList<Complex> getomega(int n) {
        // eta sign erratum vs the reference paper: see _kb/03-api-layer.md ("jline.lib.lti")
        double[] eta = geteta(n);
        ArrayList<Complex> res = new ArrayList<Complex>();
        for (int i = 1; i <= n; i++) {
            res.add(new Complex(FastMath.pow(10.0, (n - 1) / 6.0)
                    * FastMath.pow(-1.0, i - 1) * eta[i - 1]));
        }
        return res;
    }

    public static void main(String[] args) {
        getalpha(99);
        getomega(99);
        // abatewhitt.getResult(alpha, omega);
    }
}
