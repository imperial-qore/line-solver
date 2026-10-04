/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.TreeMap;

import org.apache.commons.math3.util.FastMath;

public final class gaverstehfest {
    private gaverstehfest() {}

    public static BigDecimal[] getomega(int n) {
        if (n % 2 == 1) {
            n--;
        }
        BigDecimal[] res = new BigDecimal[n];
        for (int k = 1; k <= n; k++) {
            BigDecimal val = new BigDecimal("1");
            val = val.multiply(BigDecimal.valueOf(Math.pow(-1.0, (n / 2.0) + k)));
            val = val.multiply(BigDecimal.valueOf(Math.log(2.0)));
            BigDecimal sum_val = new BigDecimal("0");
            int j = (int) FastMath.floor((k + 1) / 2.0);
            while (j <= FastMath.min((double) k, n / 2.0)) {
                BigDecimal val2 = BigDecimal.valueOf(Math.pow((double) j, n / 2.0 + 1));
                val2 = val2.divide(new BigDecimal(abatewhitt.factorial(n / 2)),
                        9999, RoundingMode.HALF_EVEN);
                val2 = val2.multiply(new BigDecimal(abatewhitt.binom(n / 2, j)));
                val2 = val2.multiply(new BigDecimal(abatewhitt.binom(2 * j, j)));
                val2 = val2.multiply(new BigDecimal(abatewhitt.binom(j, k - j)));
                sum_val = sum_val.add(val2);
                j++;
            }
            val = val.multiply(sum_val);
            res[k - 1] = val;
        }
        return res;
    }

    public static BigDecimal[] getalpha(int n) {
        if (n % 2 == 1) n--;

        BigDecimal[] res = new BigDecimal[n];
        for (int k = 1; k <= n; k++) {
            res[k - 1] = BigDecimal.valueOf(k * FastMath.log(2.0));
        }
        return res;
    }

    public static double laplace_result(double d) {
        return (1.0 / d) / (Math.exp(d) - 1);
    }

    public static void main(String[] args) {
        BigDecimal[] alpha = getalpha(12);
        BigDecimal[] omega = getomega(12);
        TreeMap<Double, BigDecimal> res = new TreeMap<Double, BigDecimal>();
        double d = 0.1;
        while (d <= 15) {
            BigDecimal sum_val = BigDecimal.valueOf(0);
            for (int i = 0; i < alpha.length; i++) {
                sum_val = sum_val.add(omega[i].multiply(BigDecimal.valueOf(
                        laplace_result(Double.parseDouble(alpha[i].toString()) / d))));
            }
            BigDecimal r = sum_val.divide(BigDecimal.valueOf(d), 9999, RoundingMode.HALF_EVEN);
            res.put(d, r);
            d += 0.05;
        }
        for (Double dk : res.keySet()) {
            System.out.println(dk + " " + res.get(dk).doubleValue());
        }
    }
}
