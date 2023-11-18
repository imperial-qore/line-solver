package jline.lib.ltinversion;

import org.apache.commons.math3.complex.Complex;

import java.util.ArrayList;
import java.util.Arrays;

public class euler {
    public static ArrayList<Complex> getalpha(int n) {
        ArrayList<Complex> result = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            result.add(new Complex((n - 1) * Math.log(10) / 6, Math.PI * i));
        }
        return result;
    }

    public static double[] geteta(int n) {
        double[] res = new double[n];
        res[0] = 0.5;
        // euler defined only for odd n!
        for (int i = 2; i <= (n + 1) / 2; i++) {
            res[i - 1] = 1;
        }
        res[n - 1] = 1 / Math.pow(2, (n - 1) / 2d);
        for (int i = 1; i < (n - 1) / 2; i++) {
            res[n - i - 1] = res[n - i] + Math.pow(2, (1 - n) / 2d) * Double.parseDouble(abatewhitt.binom((n - 1) / 2, i).toString());
        }
        return res;
    }

    public static ArrayList<Complex> getomega(int n) {
        /*
         * Page 6 of Numerical inverse Laplace transformation using concentrated matrix exponential distributions contains an error
         * specifically, that eta starts off negative, when it should be positive as k starts from 1 there
         * This can be validated by going to equation 36 of http://www.columbia.edu/~ww2040/AbateUnified2006.pdf - there k starts from 0 and hence (-1)^k is positive
         * the difference is that otherwise we would get -f(x) instead of f(x) as the result
         */
        double[] eta = geteta(n);
        ArrayList<Complex> res = new ArrayList<>();
        for (int i = 1; i <= n; i++) {
            res.add(new Complex(Math.pow(10, (n - 1) / 6d) * Math.pow(-1, i - 1) * eta[i - 1]));
        }
        return res;
    }

    public static void main(String[] args) {
        ArrayList<Complex> alpha = getalpha(99);
        ArrayList<Complex> omega = getomega(99);
     //   abatewhitt.getResult(alpha, omega);
    }
}