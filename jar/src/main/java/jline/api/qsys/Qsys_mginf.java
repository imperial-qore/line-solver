package jline.api.qsys;

import java.util.HashMap;

import jline.GlobalConstants;

public final class Qsys_mginf {
    private Qsys_mginf() {}

    /**
     * M/G/inf queue analysis (infinite servers).
     */
    public static HashMap<String, Object> qsys_mginf(double lambda, double mu, double cv2) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        double rho = lambda / mu;

        double L = rho;
        double Lq = 0.0;
        double W = 1.0 / mu;
        double Wq = 0.0;
        double p0 = Math.exp(-rho);

        result.put("L", L);
        result.put("Lq", Lq);
        result.put("W", W);
        result.put("Wq", Wq);
        result.put("p0", p0);

        return result;
    }

    /**
     * M/G/inf queue analysis with state probability.
     */
    public static HashMap<String, Object> qsys_mginf(double lambda, double mu, double cv2, int k) {
        HashMap<String, Object> result = qsys_mginf(lambda, mu, cv2);

        double rho = lambda / mu;
        double pk = Math.exp(-rho) * Math.pow(rho, (double) k) / factorial(k);
        result.put("pk", pk);

        return result;
    }

    private static double factorial(int n) {
        if (n == 0 || n == 1) return 1.0;
        double result = 1.0;
        for (int i = 2; i <= n; i++) {
            result *= i;
            if (Double.isInfinite(result)) {
                return GlobalConstants.Inf;
            }
        }
        return result;
    }
}
