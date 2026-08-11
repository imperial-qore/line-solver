/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.lti;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.util.FastMath;

import java.math.BigDecimal;
import java.math.RoundingMode;

public final class customromberg {
    private customromberg() {}

    public static BigDecimal[][] dp;
    public static int precision = 35;

    public static BigDecimal hn(double a, double b, int n) {
        BigDecimal bg = new BigDecimal(b - a);
        return bg.divide(BigDecimal.valueOf((long) n), precision, RoundingMode.HALF_EVEN);
    }

    public static BigDecimal solver(int m, int n, UnivariateFunction f, double lb, double ub) {
        if (dp[m][n] != null) return dp[m][n];
        if (m == 0 && n == 0) {
            dp[m][n] = hn(lb, ub, 1).multiply(BigDecimal.valueOf(f.value(lb) + f.value(ub)));
            return dp[m][n];
        } else if (n == 0) {
            BigDecimal temp = hn(lb, ub, m).divide(BigDecimal.valueOf(2), precision, RoundingMode.HALF_EVEN);
            BigDecimal toadd = BigDecimal.valueOf(f.value(lb));
            for (int j = 1; j <= m - 1; j++) {
                toadd = toadd.add(BigDecimal.valueOf(2)
                        .multiply(BigDecimal.valueOf(f.value(lb + j * hn(lb, ub, m).doubleValue()))));
            }
            toadd = toadd.add(BigDecimal.valueOf(f.value(ub)));
            temp = temp.multiply(toadd);
            dp[m][n] = temp;
            return dp[m][n];
        } else {
            BigDecimal temp = new BigDecimal(1);
            BigDecimal temp2 = new BigDecimal(4).pow(n);
            temp2 = temp2.subtract(new BigDecimal(1));
            temp = temp.divide(temp2, precision, RoundingMode.HALF_EVEN);
            BigDecimal smn1 = solver(m, n - 1, f, lb, ub);
            BigDecimal sm1n1 = solver(m - 1, n - 1, f, lb, ub);
            temp = temp.multiply(smn1.subtract(sm1n1));
            BigDecimal res = solver(m, n - 1, f, lb, ub).add(temp);
            dp[m][n] = res;
            return dp[m][n];
        }
    }

    public static void main(String[] args) {
        UnivariateFunction f = new UnivariateFunction() {
            @Override
            public double value(double t) { return FastMath.exp(-t); }
        };
        System.out.println(starter(f, 1.0, 2.0));
        for (int i = 0; i < dp.length; i++) {
            for (int j = 0; j <= i; j++) {
                System.out.print(dp[i][j].doubleValue() + " ");
            }
            System.out.println();
        }
    }

    public static BigDecimal starter(UnivariateFunction f, double lb, double ub) {
        int max_steps = 20;
        int number_of_iterations = 0;
        dp = new BigDecimal[max_steps][max_steps];
        BigDecimal result = solver(max_steps - 1, max_steps - 1, f, lb, ub);
        BigDecimal prev_value = dp[0][0];
        for (int i = 1; i < max_steps; i++) {
            for (int j = 0; j <= i; j++) {
                if (dp[i][j].subtract(prev_value).abs()
                        .compareTo(BigDecimal.valueOf(0.00000000000001)) < 0 && number_of_iterations >= 20) {
                    return dp[i][j];
                } else {
                    prev_value = dp[i][j];
                }
                number_of_iterations++;
            }
        }
        return result;
    }
}
