package jline.lib.butools;

public final class NormMomsFromMoms {
    private NormMomsFromMoms() {}

    public static double[] NormMomsFromMoms(double[] m) {
        int length = m.length;
        double[] nm = new double[length];
        nm[0] = m[0];

        for (int i = 1; i < length; i++) {
            nm[i] = m[i] / (m[i - 1] * m[0]);
        }

        return nm;
    }
}
