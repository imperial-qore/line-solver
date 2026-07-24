package jline.lib.butools;

import jline.GlobalConstants;

public final class APH3rdMomentUpperBound {
    private APH3rdMomentUpperBound() {}

    public static double APH3rdMomentUpperBound(double m1, double m2, int n) {
        double n2 = m2 / m1 / m1;
        if (n2 < (n + 1.0) / n) {
            return GlobalConstants.NegInf;
        } else if (n2 <= n / (n - 1.0)) {
            return m1 * m2 * (2.0 * (n - 2.0) * (n * n2 - n - 1.0) * Math.sqrt(1.0 + (n * (n2 - 2.0)) / (n - 1.0)) + (n + 2.0) * (3.0 * n * n2 - 2.0 * n - 2.0)) / (n * n * n2);
        } else {
            return GlobalConstants.Inf;
        }
    }
}
