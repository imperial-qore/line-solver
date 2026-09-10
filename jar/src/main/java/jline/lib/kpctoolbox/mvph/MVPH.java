/**
 * @file Multivariate Phase-Type Distribution (MVPH) functions
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/mvph/
 *
 * @since LINE 3.0
 */
package jline.lib.kpctoolbox.mvph;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import jline.util.matrix.Matrix;

public final class MVPH {
    private MVPH() {}

    /**
     * Computes the joint moment E[X^n1 * Y^n2] of a bivariate phase-type distribution.
     */
    public static double mvph_joint(double[] alpha, Matrix S, Matrix T, Matrix D, int n1, int n2) {
        int sSize = S.getNumRows();
        int tSize = T.getNumRows();

        RealMatrix negS = MatrixUtils.createRealMatrix(sSize, sSize);
        for (int i = 0; i < sSize; i++) {
            for (int j = 0; j < sSize; j++) {
                negS.setEntry(i, j, -S.get(i, j));
            }
        }

        RealMatrix invNegS = new LUDecomposition(negS).getSolver().getInverse();

        RealMatrix invNegSPow = MatrixUtils.createRealIdentityMatrix(sSize);
        for (int p = 0; p < n1 + 1; p++) {
            invNegSPow = invNegSPow.multiply(invNegS);
        }

        RealMatrix negT = MatrixUtils.createRealMatrix(tSize, tSize);
        for (int i = 0; i < tSize; i++) {
            for (int j = 0; j < tSize; j++) {
                negT.setEntry(i, j, -T.get(i, j));
            }
        }

        RealMatrix invNegT = new LUDecomposition(negT).getSolver().getInverse();

        RealMatrix invNegTPow = MatrixUtils.createRealIdentityMatrix(tSize);
        for (int p = 0; p < n2 + 1; p++) {
            invNegTPow = invNegTPow.multiply(invNegT);
        }

        RealMatrix dMatrix = MatrixUtils.createRealMatrix(sSize, tSize);
        for (int i = 0; i < sSize; i++) {
            for (int j = 0; j < tSize; j++) {
                dMatrix.setEntry(i, j, D.get(i, j));
            }
        }

        RealMatrix minusT = MatrixUtils.createRealMatrix(tSize, tSize);
        for (int i = 0; i < tSize; i++) {
            for (int j = 0; j < tSize; j++) {
                minusT.setEntry(i, j, -T.get(i, j));
            }
        }

        RealMatrix alphaVec = MatrixUtils.createRowRealMatrix(alpha);
        RealMatrix step1 = alphaVec.multiply(invNegSPow);
        RealMatrix step2 = step1.multiply(dMatrix);
        RealMatrix step3 = step2.multiply(invNegTPow);
        RealMatrix step4 = step3.multiply(minusT);

        double[] onesArr = new double[tSize];
        for (int i = 0; i < tSize; i++) onesArr[i] = 1.0;
        RealMatrix ones = MatrixUtils.createColumnRealMatrix(onesArr);
        RealMatrix step5 = step4.multiply(ones);

        double result = step5.getEntry(0, 0);

        return factorial(n1) * factorial(n2) * result;
    }

    public static double mvph_mean_x(double[] alpha, Matrix S, Matrix T, Matrix D) {
        return mvph_joint(alpha, S, T, D, 1, 0);
    }

    public static double mvph_mean_y(double[] alpha, Matrix S, Matrix T, Matrix D) {
        return mvph_joint(alpha, S, T, D, 0, 1);
    }

    public static double mvph_cov(double[] alpha, Matrix S, Matrix T, Matrix D) {
        double exy = mvph_joint(alpha, S, T, D, 1, 1);
        double ex = mvph_joint(alpha, S, T, D, 1, 0);
        double ey = mvph_joint(alpha, S, T, D, 0, 1);
        return exy - ex * ey;
    }

    public static double mvph_corr(double[] alpha, Matrix S, Matrix T, Matrix D) {
        double ex = mvph_joint(alpha, S, T, D, 1, 0);
        double ey = mvph_joint(alpha, S, T, D, 0, 1);
        double ex2 = mvph_joint(alpha, S, T, D, 2, 0);
        double ey2 = mvph_joint(alpha, S, T, D, 0, 2);
        double exy = mvph_joint(alpha, S, T, D, 1, 1);

        double varX = ex2 - ex * ex;
        double varY = ey2 - ey * ey;
        double covXY = exy - ex * ey;

        if (varX <= 0 || varY <= 0) {
            return 0.0;
        }

        return covXY / (Math.sqrt(varX) * Math.sqrt(varY));
    }

    private static double factorial(int n) {
        if (n <= 1) return 1.0;
        double result = 1.0;
        for (int i = 2; i <= n; i++) {
            result *= (double) i;
        }
        return result;
    }
}
