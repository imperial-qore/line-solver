package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_avgl1linfty {
    private Ms_avgl1linfty() {}

    /**
     * Average L1 L-infinity distance between two probability distributions.
     * Part of the combination family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Average L1 L-infinity distance
     */
    public static double ms_avgl1linfty(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sumL1 = 0.0;
        double maxLinfty = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double diff = Math.abs(P.get(i) - Q.get(i));
            sumL1 += diff;
            if (diff > maxLinfty) {
                maxLinfty = diff;
            }
        }

        return (sumL1 + maxLinfty) / 2.0;
    }
}
