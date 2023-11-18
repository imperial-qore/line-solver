package jline.lib.ltinversion;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.TreeMap;

import org.apache.commons.math3.complex.Complex;

public class talbot {
    // represent number as Pair<Double, Double>(real, complex)
    // pi = Math.PI
    public static ArrayList<Complex> getalpha(int n) {
        ArrayList<Complex> arr = new ArrayList<>();
        // for k = 1
        arr.add(new Complex(2d * n / 5d, 0));
        // for k = 2 onwards
        for (int i = 2; i <= n; i++) {
            // get real and imaginary part
            double real_part = 2 * (i - 1) * Math.PI / 5 * (1 / Math.tan((i - 1) * Math.PI / n));
            double imaginary_part = 2 * (i - 1) * Math.PI / 5;
            arr.add(new Complex(real_part, imaginary_part));
        }
        return arr;
    }

    public static ArrayList<Complex> getomega(int n, ArrayList<Complex> alpha) {
        ArrayList<Complex> arr = new ArrayList<>();
        // for k = 1
        arr.add((alpha.get(0).exp()).divide(5));
        //  arr.add(new Complex(Math.exp(alpha.get(0).getReal())/5, 0d));
        // for k = 2 onwards
        for (int i = 2; i <= n; i++) {
            Complex current_alpha = alpha.get(i - 1).exp();
            double temp_var = (i - 1) * Math.PI / n;
            Complex result = current_alpha.multiply(2).divide(5).multiply(new Complex(1, temp_var * (1 + 1 / Math.pow(Math.tan(temp_var), 2)) - 1 / Math.tan(temp_var)));
            arr.add(result);
        }
        return arr;
    }
//    public static void main(String[] args) {
//        ArrayList<Complex> alpha = getalpha(92);
//        ArrayList<Complex> omega = getomega(92, alpha);
//        abatewhitt.getResult(alpha, omega);
// //       System.out.println(res);
//    }
}
