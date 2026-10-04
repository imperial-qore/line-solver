/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.lib.butools;

import org.apache.commons.math3.complex.Complex;

import jline.GlobalConstants;
import jline.lang.processes.APH;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class APHFrom3Moments {
    private APHFrom3Moments() {}

    public static APH APHFrom3Moments(double[] moms) {
        return APHFrom3Moments(moms, 100, 1e-14);
    }

    public static APH APHFrom3Moments(double[] moms, int maxSize) {
        return APHFrom3Moments(moms, maxSize, 1e-14);
    }

    public static APH APHFrom3Moments(double[] moms, int maxSize, double prec) {
        if (maxSize <= 0) {
            maxSize = 100;
        }
        if (prec <= 0) {
            prec = 1e-14;
        }

        double m1 = moms[0];
        double m2 = moms[1];
        double m3 = moms[2];

        // Detect number of phases needed
        int n = 2;
        while (n < maxSize && (APH2ndMomentLowerBound.APH2ndMomentLowerBound(m1, n) > m2
                || APH3rdMomentLowerBound.APH3rdMomentLowerBound(m1, m2, n) >= m3
                || APH3rdMomentUpperBound.APH3rdMomentUpperBound(m1, m2, n) <= m3)) {
            n++;
        }

        if (APH2ndMomentLowerBound.APH2ndMomentLowerBound(m1, n) > m2) {
            m2 = APH2ndMomentLowerBound.APH2ndMomentLowerBound(m1, n);
        }

        if (APH3rdMomentLowerBound.APH3rdMomentLowerBound(m1, m2, n) > m3) {
            m3 = APH3rdMomentLowerBound.APH3rdMomentLowerBound(m1, m2, n);
        }

        if (APH3rdMomentUpperBound.APH3rdMomentUpperBound(m1, m2, n) < m3) {
            m3 = APH3rdMomentUpperBound.APH3rdMomentUpperBound(m1, m2, n);
        }

        // Compute normalized moments
        double[] nmoms = NormMomsFromMoms.NormMomsFromMoms(new double[]{m1, m2, m3});
        double n1 = nmoms[0];
        double n2 = nmoms[1];
        double n3 = nmoms[2];

        if (n2 > 2 || n3 < 2 * n2 - 1) {
            Complex nComplex = new Complex((double) n, 0.0);
            Complex n2Complex = new Complex(n2, 0.0);
            Complex n3Complex = new Complex(n3, 0.0);
            Complex four = new Complex(4.0, 0.0);
            Complex two = new Complex(2.0, 0.0);
            Complex three = new Complex(3.0, 0.0);
            Complex twelve = new Complex(12.0, 0.0);
            Complex sixteen = new Complex(16.0, 0.0);
            Complex fifteen = new Complex(15.0, 0.0);
            Complex eight = new Complex(8.0, 0.0);
            Complex one = Complex.ONE;

            Complex numeratorB = two.multiply(four.subtract(nComplex.multiply(three.multiply(n2Complex).subtract(four))));
            Complex denominatorB = n2Complex.multiply(four.add(nComplex).subtract(nComplex.multiply(n3Complex)))
                    .add((nComplex.multiply(n2Complex)).sqrt()
                            .multiply((twelve.multiply(n2Complex.multiply(n2Complex)).multiply(nComplex.add(one))
                                    .add(sixteen.multiply(n3Complex).multiply(nComplex.add(one)))
                                    .add(n2Complex.multiply(nComplex.multiply(n3Complex.subtract(fifteen)).multiply(n3Complex.add(one))
                                            .subtract(eight.multiply(n3Complex.add(new Complex(3.0, 0.0)))))).sqrt())));

            Complex b = numeratorB.divide(denominatorB);

            Complex numeratorA = (b.multiply(n2Complex).subtract(two)).multiply(nComplex.subtract(one)).multiply(b);
            Complex denominatorA = b.subtract(one).multiply(nComplex);
            Complex a = numeratorA.divide(denominatorA);

            Complex p = (b.subtract(one)).divide(a);

            Complex lambda = (p.multiply(a).add(one)).divide(n1);

            Complex mu = (nComplex.subtract(one)).multiply(lambda).divide(a);

            // Construct representation
            double[] alpha = new double[n];
            alpha[0] = p.getReal();
            alpha[n - 1] = 1.0 - p.getReal();

            Matrix A = new Matrix(n, n);
            A.set(n - 1, n - 1, -lambda.getReal());
            for (int i = 0; i < n - 1; i++) {
                A.set(i, i, -mu.getReal());
                A.set(i, i + 1, mu.getReal());
            }
            return new APH(new Matrix(alpha).transpose(), A);
        } else {
            double c4 = n2 * (3.0 * n2 - 2.0 * n3) * Math.pow((double) (n - 1), 2.0);
            double c3 = 2.0 * n2 * (n3 - 3.0) * Math.pow((double) (n - 1), 2.0);
            double c2 = 6.0 * (n - 1.0) * (n - n2);
            double c1 = 4.0 * n * (2.0 - n);
            double c0 = n * (n - 2.0);

            double[] coefficients = new double[]{c0, c1, c2, c3, c4};
            Complex[] fs = Maths.roots(coefficients);

            for (Complex f : fs) {
                if (f.isNaN() || f.isInfinite()) {
                    continue;
                }

                Complex temp1 = f.multiply(f).multiply(n2).subtract(f.multiply(2.0)).add(2.0);
                Complex temp2 = temp1.multiply(n - 1);
                Complex temp3 = temp2.subtract((double) n);

                if (temp3.abs() < prec) {
                    continue;
                }

                Complex a = f.subtract(1.0).multiply(2.0).multiply(n - 1).divide(temp3);
                Complex p = f.subtract(1.0).multiply(a);
                Complex lambda = a.add(p).divide(n1);
                Complex mu = new Complex((double) (n - 1)).divide(n1 - p.divide(lambda).getReal());

                if (Math.abs(p.getImaginary()) < GlobalConstants.Zero) {
                    if (Math.abs(lambda.getImaginary()) < GlobalConstants.Zero) {
                        if (Math.abs(mu.getImaginary()) < GlobalConstants.Zero) {
                            if (p.getReal() >= 0 && p.getReal() <= 1 && lambda.getReal() > 0 && mu.getReal() > 0) {
                                double[] alpha = new double[n];
                                alpha[0] = p.getReal();
                                alpha[1] = 1.0 - p.getReal();

                                Matrix A = new Matrix(n, n);
                                A.set(0, 0, -lambda.getReal());
                                A.set(0, 1, lambda.getReal());
                                for (int j = 1; j < n; j++) {
                                    A.set(j, j, -mu.getReal());
                                    if (j < n - 1) {
                                        A.set(j, j + 1, mu.getReal());
                                    }
                                }
                                return new APH(new Matrix(alpha).transpose(), A);
                            }
                        }
                    }
                }
            }
        }

        throw new IllegalArgumentException("No APH found for the given 3 moments!");
    }
}
