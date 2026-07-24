/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.reptrans;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class TransformToMonocyclic {
    private TransformToMonocyclic() {}

    /**
     * Data class representing a Feedback Erlang Block (FEB).
     */
    private static final class FEB {
        final double lambda;
        final double z;
        final int n;
        final int multip;
        final List<Complex> evals;
        final List<Integer> emuls;

        FEB(double lambda, double z, int n, int multip, List<Complex> evals, List<Integer> emuls) {
            this.lambda = lambda;
            this.z = z;
            this.n = n;
            this.multip = multip;
            this.evals = evals;
            this.emuls = emuls;
        }
    }

    /**
     * Transforms an arbitrary matrix to a Markovian monocyclic matrix.
     */
    public static Matrix transformToMonocyclic(Matrix A, int maxSize, double precision) {
        Pair<List<Complex>, List<Integer>> evalsResult = eigvalc(A, precision);
        List<Complex> evalues = evalsResult.getLeft();
        List<Integer> repeats = evalsResult.getRight();

        List<FEB> febs = generateFEBs(evalues, repeats, maxSize, precision);

        int totalSize = 0;
        for (FEB feb : febs) {
            totalSize += feb.n * feb.multip;
        }

        Matrix B = Matrix.zeros(totalSize, totalSize);
        int pos = 0;
        for (int i = 0; i < febs.size(); i++) {
            int Ni = febs.get(i).n * febs.get(i).multip;
            Matrix febGen = febGenerator(febs.get(i).lambda, febs.get(i).z, febs.get(i).n, febs.get(i).multip);

            for (int r = 0; r < Ni; r++) {
                for (int c = 0; c < Ni; c++) {
                    B.set(pos + r, pos + c, febGen.get(r, c));
                }
            }

            if (i < febs.size() - 1) {
                double rowSum = 0.0;
                for (int c = 0; c < totalSize; c++) {
                    rowSum += B.get(pos + Ni - 1, c);
                }
                B.set(pos + Ni - 1, pos + Ni, -rowSum);
            }
            pos += Ni;
        }

        return B;
    }

    public static Matrix transformToMonocyclic(Matrix A, int maxSize) {
        return transformToMonocyclic(A, maxSize, 1e-14);
    }

    public static Matrix transformToMonocyclic(Matrix A) {
        return transformToMonocyclic(A, 100, 1e-14);
    }

    /**
     * Computes eigenvalues and their algebraic multiplicity.
     */
    private static Pair<List<Complex>, List<Integer>> eigvalc(Matrix A, double precision) {
        double tol = Math.sqrt(Math.ulp(1.0));

        List<Complex> eigArr = A.eig();
        List<Complex> eigenvalues = new ArrayList<Complex>(eigArr);
        Collections.sort(eigenvalues, new Comparator<Complex>() {
            @Override
            public int compare(Complex a, Complex b) {
                return Double.compare(a.getReal(), b.getReal());
            }
        });

        for (int i = 0; i < eigenvalues.size(); i++) {
            Complex rounded = new Complex(
                    Math.round(eigenvalues.get(i).getReal() / tol) * tol,
                    Math.round(eigenvalues.get(i).getImaginary() / tol) * tol);
            eigenvalues.set(i, rounded);
        }

        List<Complex> uniqueEvals = new ArrayList<Complex>();
        List<Integer> repeats = new ArrayList<Integer>();

        for (Complex ev : eigenvalues) {
            boolean found = false;
            for (int j = 0; j < uniqueEvals.size(); j++) {
                if (ev.subtract(uniqueEvals.get(j)).abs() <= tol) {
                    repeats.set(j, repeats.get(j) + 1);
                    found = true;
                    break;
                }
            }
            if (!found) {
                uniqueEvals.add(ev);
                repeats.add(1);
            }
        }

        return new Pair<List<Complex>, List<Integer>>(uniqueEvals, repeats);
    }

    /**
     * Generates Feedback Erlang Blocks from eigenvalues.
     */
    private static List<FEB> generateFEBs(List<Complex> evalues, List<Integer> repeats, int maxSize, double precision) {
        List<FEB> febs = new ArrayList<FEB>();
        int i = 0;
        int size = 0;

        while (i < evalues.size()) {
            int multip = repeats.get(i);
            double evalimag = -Math.abs(evalues.get(i).getImaginary());

            if (-evalimag < precision) {
                int n = 1;
                double sigma = -evalues.get(i).getReal();
                double z = 0.0;
                List<Complex> ev = new ArrayList<Complex>();
                ev.add(new Complex(evalues.get(i).getReal(), 0.0));
                List<Integer> em = new ArrayList<Integer>();
                em.add(multip);
                size += 1;

                febs.add(new FEB(sigma, z, n, multip, ev, em));
                i++;
            } else {
                int n = 3;
                size += 3;
                // The feedback probability z is strictly less than 1 only if the
                // ratio stays strictly below cot(pi/n). At the exact boundary
                // (ratio == cot(pi/n)) rounding can let the loop exit with z == 1,
                // which yields a block with zero exit rate, hence a degenerate
                // generator with a zero eigenvalue. Keep a margin so that the
                // boundary is always resolved upwards, as in MATLAB
                // TransformToMonocyclic.m.
                double boundaryMargin = Math.sqrt(Math.ulp(1.0));
                while (evalimag / evalues.get(i).getReal() >= 1.0 / Math.tan(Math.PI / n) - boundaryMargin) {
                    n++;
                    size++;
                    if (size > maxSize) {
                        throw new IllegalArgumentException(
                                "The representation is too large (>maxSize). No result returned.");
                    }
                }

                double sigma = -(2 * evalues.get(i).getReal()
                        + evalimag * (1.0 / Math.tan(Math.PI / n) - Math.tan(Math.PI / n))) / 2;
                double z = Math.pow(
                        -evalimag * (1.0 / Math.tan(Math.PI / n) + Math.tan(Math.PI / n)) / (2 * sigma),
                        n);

                List<Complex> ev = new ArrayList<Complex>();
                List<Integer> em = new ArrayList<Integer>();
                for (int k = 0; k < n; k++) {
                    double zRoot = Math.pow(z, 1.0 / n);
                    double angle = 2 * k * Math.PI / n;
                    ev.add(new Complex(
                            -(1 - zRoot * Math.cos(angle)) * sigma,
                            zRoot * Math.sin(angle) * sigma));
                    em.add(multip);
                }

                febs.add(new FEB(sigma, z, n, multip, ev, em));
                i += 2;
            }
        }

        Collections.sort(febs, new Comparator<FEB>() {
            @Override
            public int compare(FEB a, FEB b) {
                Complex aMax = null;
                for (Complex c : a.evals) {
                    if (aMax == null || c.abs() > aMax.abs()) aMax = c;
                }
                Complex bMax = null;
                for (Complex c : b.evals) {
                    if (bMax == null || c.abs() > bMax.abs()) bMax = c;
                }
                double aReal = (aMax != null) ? aMax.getReal() : Double.MIN_VALUE;
                double bReal = (bMax != null) ? bMax.getReal() : Double.MIN_VALUE;
                return Double.compare(bReal, aReal);
            }
        });

        return febs;
    }

    /**
     * Generates the generator matrix for a single FEB.
     */
    private static Matrix febGenerator(double lambda, double z, int n, int multip) {
        int size = multip * n;
        Matrix A = Matrix.zeros(size, size);

        for (int ii = 0; ii < multip; ii++) {
            for (int j = 0; j < n; j++) {
                A.set(ii * n + j, ii * n + j, -lambda);
                if (j < n - 1) {
                    A.set(ii * n + j, ii * n + j + 1, lambda);
                }
            }
            A.set(ii * n + n - 1, ii * n,
                    A.get(ii * n + n - 1, ii * n) + z * lambda);

            if (ii < multip - 1) {
                A.set(ii * n + n - 1, (ii + 1) * n, (1 - z) * lambda);
            }
        }

        return A;
    }
}
