/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.reptrans;

import jline.util.matrix.Matrix;

public final class ExtendToMarkovian {
    private ExtendToMarkovian() {}

    /**
     * Extends a non-Markovian initial vector to a Markovian one by appending an Erlang tail.
     */
    public static MarkovianRepresentation extendToMarkovian(Matrix alpha, Matrix A, int maxSize, double precision) {
        int N = A.getNumRows();

        double t0lower = 0.0;
        double t0upper = 1.0;
        Matrix beta = alpha.mult(A.scale(t0upper).expm());

        while (minOf(beta.toArray1D()) < -precision) {
            t0upper *= 2.0;
            beta = alpha.mult(A.scale(t0upper).expm());
        }

        while ((t0upper - t0lower) / (t0upper + t0lower) > precision) {
            double t0Try = (t0upper + t0lower) / 2.0;
            beta = alpha.mult(A.scale(t0Try).expm());
            if (minOf(beta.toArray1D()) < -precision) {
                t0lower = t0Try;
            } else {
                t0upper = t0Try;
            }
        }
        double t0 = t0upper;

        double increment = 1.1;
        double bestT0 = -1.0;
        int bestLupper = -1;

        for (int iter = 0; iter < 100; iter++) {
            int Llower = 1;
            int Lupper = 1;
            beta = inivecWithTail(alpha, A, Lupper, (double) Lupper / t0);

            while (minOf(beta.toArray1D()) < -precision && Lupper < maxSize) {
                Lupper *= 2;
                beta = inivecWithTail(alpha, A, Lupper, (double) Lupper / t0);
            }

            boolean success = minOf(beta.toArray1D()) >= -precision;

            if (success) {
                while (Lupper - Llower > 1) {
                    int L = (Lupper + Llower) / 2;
                    beta = inivecWithTail(alpha, A, L, (double) L / t0);
                    if (minOf(beta.toArray1D()) < -precision) {
                        Llower = L;
                    } else {
                        Lupper = L;
                    }
                }
            }

            if (success) {
                if (bestLupper >= 0 && Lupper > bestLupper) {
                    break;
                } else {
                    bestLupper = Lupper;
                    bestT0 = t0;
                    t0 *= increment;
                }
            } else {
                if (bestLupper >= 0) {
                    break;
                } else {
                    t0 *= increment;
                }
            }
        }

        if (bestLupper < 0) {
            throw new IllegalArgumentException("No positive representation found up to the given size!");
        }

        t0 = bestT0;
        int Lupper = bestLupper;

        beta = inivecWithTail(alpha, A, Lupper, (double) Lupper / t0);
        Matrix B = addErlangTail(A, Lupper, (double) Lupper / t0);

        return new MarkovianRepresentation(beta, B);
    }

    public static MarkovianRepresentation extendToMarkovian(Matrix alpha, Matrix A, int maxSize) {
        return extendToMarkovian(alpha, A, maxSize, 1e-14);
    }

    public static MarkovianRepresentation extendToMarkovian(Matrix alpha, Matrix A) {
        return extendToMarkovian(alpha, A, 100, 1e-14);
    }

    private static double minOf(double[] arr) {
        if (arr.length == 0) throw new IllegalArgumentException("empty array");
        double m = arr[0];
        for (int i = 1; i < arr.length; i++) {
            if (arr[i] < m) m = arr[i];
        }
        return m;
    }

    /**
     * Computes the initial vector with Erlang tail.
     */
    private static Matrix inivecWithTail(Matrix gamma, Matrix G, int tailLength, double mu) {
        int vlen = G.getNumRows() + tailLength;
        Matrix beta = Matrix.zeros(1, vlen);

        Matrix WG = Matrix.eye(G.getNumRows()).add(G.scale(1.0 / mu));
        Matrix opv = gamma.copy();

        Matrix clv = new Matrix(G.getNumRows(), 1);
        for (int i = 0; i < G.getNumRows(); i++) {
            double sum = 0.0;
            for (int j = 0; j < G.getNumCols(); j++) {
                sum += G.get(i, j) / mu;
            }
            clv.set(i, 0, -sum);
        }

        for (int k = vlen - 1; k >= G.getNumRows(); k--) {
            beta.set(0, k, opv.mult(clv).get(0, 0));
            opv = opv.mult(WG);
        }

        for (int i = 0; i < G.getNumRows(); i++) {
            beta.set(0, i, opv.get(0, i));
        }

        return beta;
    }

    /**
     * Adds an Erlang tail to a generator matrix.
     */
    private static Matrix addErlangTail(Matrix D, int len, double mu) {
        int DN = D.getNumRows();
        Matrix E = Matrix.zeros(DN + len, DN + len);

        for (int i = 0; i < DN; i++) {
            for (int j = 0; j < DN; j++) {
                E.set(i, j, D.get(i, j));
            }
        }

        double rowSum = 0.0;
        for (int j = 0; j < DN; j++) {
            rowSum += D.get(DN - 1, j);
        }
        E.set(DN - 1, DN, -rowSum);

        for (int ei = 0; ei < len; ei++) {
            E.set(DN + ei, DN + ei, -mu);
            if (ei < len - 1) {
                E.set(DN + ei, DN + ei + 1, mu);
            }
        }

        return E;
    }
}
