package jline.api.mam;

import jline.lang.constant.ProcessType;
import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Exact discrete phase-type representation of a lattice-valued law.
 *
 * <p>Returns (alpha, A) with P[X=k] = alpha*A^(k-1)*a, a = e - A*e, k = 1,2,...,
 * measuring X in slots. Geometric, Det and DiscreteUniform are represented
 * EXACTLY, not moment-matched: a fitted surrogate would leave the lattice the
 * caller relies on, so any other family is an error here.
 *
 * <p>MATLAB twin: dph_from_dist.m
 */
public final class Dph_from_dist {
    private Dph_from_dist() {}

    /**
     * @param procType process type, one of GEOMETRIC, DET, DUNIFORM
     * @param meanSlots mean of the law expressed in slots
     * @param scv squared coefficient of variation, used by DUNIFORM only
     * @return (alpha, A), the initial vector and the transient matrix
     */
    public static Pair<Matrix, Matrix> dph_from_dist(ProcessType procType, double meanSlots, double scv) {
        double tol = 1e-8;

        if (procType == ProcessType.GEOMETRIC) {
            double p = 1.0 / meanSlots;
            if (p > 1 + tol || p <= 0) {
                throw new RuntimeException("Geometric with mean " + meanSlots
                        + " slots is outside the support {1,2,...}.");
            }
            p = Math.min(1.0, p);
            Matrix alpha = new Matrix(1, 1);
            alpha.set(0, 0, 1.0);
            Matrix A = new Matrix(1, 1);
            A.set(0, 0, 1.0 - p);
            return new Pair<Matrix, Matrix>(alpha, A);
        }

        if (procType == ProcessType.DET) {
            long k = Math.round(meanSlots);
            if (Math.abs(meanSlots - k) > tol * Math.max(1.0, meanSlots) || k < 1) {
                throw new RuntimeException("Det of " + meanSlots
                        + " slots is not a positive integral number of slots.");
            }
            int n = (int) k;
            Matrix alpha = new Matrix(1, n);
            alpha.set(0, 0, 1.0);
            Matrix A = new Matrix(n, n);
            for (int i = 0; i < n - 1; i++) {
                A.set(i, i + 1, 1.0);
            }
            return new Pair<Matrix, Matrix>(alpha, A);
        }

        if (procType == ProcessType.DUNIFORM) {
            double varSlots = scv * meanSlots * meanSlots;
            double width = Math.sqrt(Math.max(0.0, 12 * varSlots + 1)) - 1;
            int lo = (int) Math.round(meanSlots - width / 2);
            int hi = (int) Math.round(meanSlots + width / 2);
            if (lo < 1 || hi < lo) {
                throw new RuntimeException("DiscreteUniform spanning [" + lo + "," + hi
                        + "] slots is outside the support {1,2,...}.");
            }
            Matrix alpha = new Matrix(1, hi);
            alpha.set(0, 0, 1.0);
            Matrix A = new Matrix(hi, hi);
            for (int j = 1; j < hi; j++) {
                // hazard of absorbing at step j, zero below the lower bound
                double h = j < lo ? 0.0 : 1.0 / (hi - j + 1);
                A.set(j - 1, j, 1.0 - h);
            }
            return new Pair<Matrix, Matrix>(alpha, A);
        }

        throw new RuntimeException("ProcessType " + procType
                + " has no exact discrete phase-type representation. The discrete-time path "
                + "accepts Geometric, Det on the slot lattice, DiscreteUniform and DMAP.");
    }
}
