package jline.lib.butools;

public final class APH2ndMomentLowerBound {
    private APH2ndMomentLowerBound() {}

    public static double APH2ndMomentLowerBound(double m1, int n) {
        return m1 * m1 * (n + 1) / n;
    }
}
