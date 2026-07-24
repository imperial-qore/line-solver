/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * G. Horvath and M. Telek, "On the canonical representation of phase type
 * distributions," Performance Evaluation, vol. 66, no. 8, pp. 396-409, 2009.
 */
package jline.lib.butools.ph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.lib.butools.ReducedMomsFromMoms;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class PH3From5Moments {
    private PH3From5Moments() {}

    /**
     * Returns a PH(3) which has the same 5 moments as given.
     */
    public static PH3Representation ph3From5Moments(double[] moms, double prec) {
        double[] rmoms = ReducedMomsFromMoms.ReducedMomsFromMoms(moms);
        for (int i = 0; i < 5; i++) {
            rmoms[i] = rmoms[i] / Math.pow(moms[0], i + 1);
        }

        Matrix M = new Matrix(3, 3);
        M.set(0, 0, rmoms[2]); M.set(0, 1, -rmoms[1]); M.set(0, 2, rmoms[0]);
        M.set(1, 0, rmoms[3]); M.set(1, 1, -rmoms[2]); M.set(1, 2, rmoms[1]);
        M.set(2, 0, rmoms[4]); M.set(2, 1, -rmoms[3]); M.set(2, 2, rmoms[2]);

        Matrix b = new Matrix(3, 1);
        b.set(0, 0, 1.0);
        b.set(1, 0, rmoms[0]);
        b.set(2, 0, rmoms[1]);

        Matrix a = M.inv().mult(b);
        double a0 = a.get(0, 0);
        double a1 = a.get(1, 0);
        double a2 = a.get(2, 0);

        double discr = a2 * a2 - 3 * a1;
        if (discr < 0) {
            throw new IllegalArgumentException("Invalid characteristic polynomial!");
        }

        double gu = (a2 + 2 * Math.sqrt(discr)) / 3;
        double g0 = (a2 + Math.sqrt(discr)) / 3;

        double[] coeffs = new double[]{a0, a1, a2, 1.0};
        Complex[] roots = Maths.roots(coeffs);

        List<Complex> sortedRoots = new ArrayList<Complex>();
        for (int i = 0; i < roots.length; i++) {
            sortedRoots.add(roots[i]);
        }
        Collections.sort(sortedRoots, new Comparator<Complex>() {
            @Override
            public int compare(Complex x, Complex y) {
                return Double.compare(x.getReal(), y.getReal());
            }
        });
        List<Complex> lambda = new ArrayList<Complex>();
        for (int i = 0; i < sortedRoots.size(); i++) {
            lambda.add(sortedRoots.get(i).negate());
        }

        double d1 = a1 - a2 - a0 * rmoms[1];
        double d2 = a0 - a1 - a2 * d1;
        double d3 = -a0 - a1 * d1 - a2 * d2;

        if (d1 > prec || (Math.abs(d1) < prec && d2 > 0)) {
            throw new IllegalArgumentException("Negative density around 0!");
        }

        if (lambda.get(2).getReal() < 0) {
            throw new IllegalArgumentException("Invalid eigenvalues!");
        }

        double gl;
        if (Math.abs(lambda.get(0).getImaginary()) < prec) {
            gl = lambda.get(0).getReal();
        } else {
            gl = g0;
        }

        if (gl > gu + prec) {
            throw new IllegalArgumentException("Invalid eigenvalues (gl>gu detected)!");
        }
        double glAdj = (gl > gu) ? gu : gl;

        double g2;
        if (Math.abs(d1) < prec) {
            g2 = 0.0;
        } else {
            g2 = -d2 / d1;
        }

        if (g2 > gu + prec) {
            throw new IllegalArgumentException("alpha_2 is negative!");
        }
        double g2Adj = (g2 > gu) ? gu : g2;

        double x1 = Math.max(g2Adj, glAdj);

        double x13;
        if (lambda.get(0).getReal() == lambda.get(0).getReal()
                && Math.abs(lambda.get(0).getImaginary()) < prec
                && g2Adj < glAdj) {
            x13 = 0.0;
        } else {
            x13 = x1 - a0 / (x1 * x1 - a2 * x1 + a1);
        }

        double bels = Math.pow(a2 - x1, 2) - 4 * (x1 * x1 - a2 * x1 + a1);
        if (bels < 0 && bels > -prec) {
            bels = 0.0;
        }

        double x2 = (a2 - x1 + Math.sqrt(bels)) / 2;
        double x3 = (a2 - x1 - Math.sqrt(bels)) / 2;
        double p1 = d1 / (x13 - x1);
        double p2 = (x1 * d1 + d2) / (x13 - x1) / x2;
        double p3 = (x1 * x2 * d1 + x2 * d2 + x1 * d2 + d3) / (x13 - x1) / x2 / x3;

        Matrix A = new Matrix(3, 3);
        A.set(0, 0, -x1 / moms[0]); A.set(0, 1, 0.0); A.set(0, 2, x13 / moms[0]);
        A.set(1, 0, x2 / moms[0]); A.set(1, 1, -x2 / moms[0]); A.set(1, 2, 0.0);
        A.set(2, 0, 0.0); A.set(2, 1, x3 / moms[0]); A.set(2, 2, -x3 / moms[0]);

        Matrix alpha = new Matrix(1, 3);
        alpha.set(0, 0, p1);
        alpha.set(0, 1, p2);
        alpha.set(0, 2, p3);

        if (x13 < -prec || x13 > x1) {
            throw new IllegalArgumentException("Invalid generator!");
        }

        double[] alphaArr = alpha.toArray1D();
        double minAlpha = Double.POSITIVE_INFINITY;
        double maxAlpha = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < alphaArr.length; i++) {
            if (alphaArr[i] < minAlpha) minAlpha = alphaArr[i];
            if (alphaArr[i] > maxAlpha) maxAlpha = alphaArr[i];
        }

        if (minAlpha < -prec) {
            throw new IllegalArgumentException("Initial vector has negative entries!");
        }

        if (maxAlpha > 1 + prec) {
            throw new IllegalArgumentException("Initial vector has entries that are greater than 1!");
        }

        return new PH3Representation(alpha, A);
    }

    public static PH3Representation ph3From5Moments(double[] moms) {
        return ph3From5Moments(moms, 1e-10);
    }
}
