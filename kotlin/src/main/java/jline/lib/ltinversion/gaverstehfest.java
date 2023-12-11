package jline.lib.ltinversion;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.TreeMap;

public class gaverstehfest {

    public static BigDecimal[] getomega(int n) {
        BigDecimal[] res = new BigDecimal[n];
        if (n % 2 == 1) {
            n--; // for gaver, as it does not support even n
        }
        for (int k = 1; k <= n; k++) {
            BigDecimal val = new BigDecimal("1");
            val = val.multiply(BigDecimal.valueOf((Math.pow(-1, (n / 2d) + k))));
            val = val.multiply(BigDecimal.valueOf(Math.log(2)));
            // summation
            BigDecimal sum_val = new BigDecimal("0");
            for (int j = (int) Math.floor((k + 1) / 2d); j <= Math.min(k, n / 2d); j++) {
                BigDecimal val2 = BigDecimal.valueOf(Math.pow(j, n / 2d + 1));
                val2 = val2.divide(new BigDecimal(abatewhitt.factorial(n / 2)), 9999, RoundingMode.HALF_EVEN); // https://stackoverflow.com/questions/4591206/arithmeticexception-non-terminating-decimal-expansion-no-exact-representable
                val2 = val2.multiply(new BigDecimal(abatewhitt.binom(n / 2, j)));
                val2 = val2.multiply(new BigDecimal(abatewhitt.binom(2 * j, j)));
                val2 = val2.multiply(new BigDecimal(abatewhitt.binom(j, k - j)));
                sum_val = sum_val.add(val2);
            }
            val = val.multiply(sum_val);
            res[k - 1] = val;
        }
        return res;
    }

    public static BigDecimal[] getalpha(int n) {
        if (n % 2 == 1)
            n--; // gaver-schefest only supports even n; same workaround done in inverselaplace.org [5].
        BigDecimal[] res = new BigDecimal[n];
        for (int k = 1; k <= n; k++) {
            res[k - 1] = BigDecimal.valueOf(k * Math.log(2));
        }
        return res;
    }

    public static double laplace_result(double d) {
        // function = 1/(1 + s^2)
    //    return 1d / (1d + d * d);
        return (1 / d) / (Math.exp(d) - 1);
    }

    public static void main(String[] args) {
        BigDecimal[] alpha = getalpha(12);
        BigDecimal[] omega = getomega(12);
        TreeMap<Double, BigDecimal> res = new TreeMap<>();
        for (double d = 0.1; d <= 15; d += 0.05) {
            BigDecimal sum_val = BigDecimal.valueOf(0);
            for (int i = 0; i < alpha.length; i++) {
                sum_val = sum_val.add(omega[i].multiply(BigDecimal.valueOf(laplace_result(Double.parseDouble(alpha[i].toString()) / d))));
            }
            BigDecimal r = sum_val.divide(BigDecimal.valueOf(d), 9999, RoundingMode.HALF_EVEN);
            res.put(d, r);
        }
        // just for verification
        for (Double d : res.keySet()) {
            System.out.println(d + " " + res.get(d).doubleValue());
        }

    }
}
