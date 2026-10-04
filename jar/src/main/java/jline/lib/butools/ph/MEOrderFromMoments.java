/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.lib.butools.ReducedMomsFromMoms;
import jline.util.matrix.Matrix;

public final class MEOrderFromMoments {
    private MEOrderFromMoments() {}

    /**
     * Returns the order of ME distribution that can realize the given moments.
     *
     * <p>The order is determined using the Hankel determinant approach:
     * we build Hankel matrices of increasing size from the reduced moments (prepended with 1)
     * and check when the determinant becomes zero.
     *
     * <p>References:
     * L. Bodrog, A. Horvath, M. Telek, "Moment characterization of matrix exponential and Markovian
     * arrival processes," Annals of Operations Research, vol. 160, pp. 51-68, 2008.
     *
     * @param moms The list of moments
     * @param prec Precision used to detect if the determinant of the Hankel matrix is zero
     * @return The order of ME distribution that can realize the given moments
     */
    public static int meOrderFromMoments(double[] moms, double prec) {
        int sizem = (moms.length + 1) / 2;
        double[] rmoms = ReducedMomsFromMoms.ReducedMomsFromMoms(moms);
        // Prepend 1 to reduced moments: rmoms_full = [1, rmoms[0], rmoms[1], ...]
        double[] rmomsAll = new double[rmoms.length + 1];
        rmomsAll[0] = 1.0;
        for (int i = 0; i < rmoms.length; i++) {
            rmomsAll[i + 1] = rmoms[i];
        }

        for (int k = 1; k <= sizem; k++) {
            Matrix hankel = new Matrix(k, k);
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    int idx = i + j;
                    if (idx < rmomsAll.length) {
                        hankel.set(i, j, rmomsAll[idx]);
                    }
                }
            }
            if (Math.abs(hankel.det()) < prec) {
                return k - 1;
            }
        }
        return sizem;
    }

    public static int meOrderFromMoments(double[] moms) {
        return meOrderFromMoments(moms, 1e-12);
    }
}
