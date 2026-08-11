/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import java.util.List;

import jline.lib.lti.iltcme;
import jline.util.matrix.Matrix;

/**
 * A Concentrated Matrix Exponential (CME) distribution.
 *
 * A CME is the matrix-exponential distribution of odd order 2*n+1 whose squared
 * coefficient of variation is (numerically) minimal for that order, from the tables of
 * Horvath, Horvath and Telek. Its SCV decays as O(1/n^2) and therefore goes far below
 * the Erlang bound 1/order attainable by a phase-type distribution of the same order:
 * order 101 gives SCV 3.9e-4, where Erlang-101 gives 9.9e-3.
 *
 * The density of the unit-mean CME with n harmonic terms is
 *
 * f(x) = mu1*exp(-mu1*x)*(c + sum_k [a_k*cos(k*w*mu1*x) + b_k*sin(k*w*mu1*x)])
 *
 * with w = omega, which is exactly alpha*expm(A*x)*(-A*e) for the block-diagonal
 * A = blkdiag(-mu1, mu1*[-1 -k*w; k*w -1], k=1..n). The parameters a, b, c, omega and
 * mu1 are read from the same iltcme.json table used by the CME inverse Laplace
 * transform, see {@link jline.lib.lti.iltcme}.
 *
 * The process type stays ME: a CME is a matrix-exponential representation, so every
 * solver gate, sn.procid entry and JSON key that accepts ME accepts it unchanged.
 */
public class CME extends ME {

    private final double mean;
    private final int order;

    /**
     * Creates a concentrated matrix-exponential distribution.
     *
     * @param mean  the mean of the distribution (positive and finite)
     * @param order the number of phases, an odd integer 2*n+1 with n present in the table
     * @throws IllegalArgumentException if the mean is not positive or the order is not tabulated
     */
    public CME(double mean, int order) {
        this(representation(order), mean, order);
    }

    /**
     * Delegating constructor: the representation must be built before the superclass
     * constructor runs, and Java allows only a single expression there, so the pair is
     * computed once by {@link #representation(int)} and threaded through.
     */
    private CME(Matrix[] rep, double mean, int order) {
        // see _kb/01-model-classes.md (CME sections) for rationale
        super(rep[0], scaleA(rep[1], mean), false);
        this.mean = mean;
        this.order = order;
    }

    private static Matrix scaleA(Matrix A, double mean) {
        if (Double.isNaN(mean) || Double.isInfinite(mean) || mean <= 0) {
            throw new IllegalArgumentException("CME mean must be a positive finite number, got " + mean);
        }
        return A.scale(1.0 / mean);
    }

    /**
     * Gets the CME order, i.e. the number of phases.
     *
     * @return the number of phases
     */
    public int getOrder() {
        return order;
    }

    /**
     * Gets the mean the distribution was constructed with.
     *
     * @return the mean
     */
    public double getConfiguredMean() {
        return mean;
    }

    /**
     * Selects the CME table entry realizing the given number of phases.
     *
     * The table is keyed by the number of harmonic terms n, so an order of 2*n+1 phases
     * maps to the entries with that n. Several entries can share an n (the "full" and
     * "approx" optimizations), and the most concentrated one is taken, matching the
     * selection rule of the CME inverse Laplace transform.
     *
     * @param order the number of phases
     * @return the selected table entry
     * @throws IllegalArgumentException if the order is not an odd integer or is not tabulated
     */
    public static iltcme.CmeEntry tableEntry(int order) {
        if (order < 3 || order % 2 == 0) {
            throw new IllegalArgumentException(
                    "CME order must be an odd integer of the form 2*n+1 with n >= 1, got " + order);
        }
        int n = (order - 1) / 2;
        iltcme.CmeEntry best = null;
        for (iltcme.CmeEntry entry : iltcme.parameters()) {
            if (entry.n == n && (best == null || entry.cv2 < best.cv2)) {
                best = entry;
            }
        }
        if (best == null) {
            int nearest = -1;
            for (int candidate : getSupportedOrders()) {
                if (nearest < 0 || Math.abs(candidate - order) < Math.abs(nearest - order)) {
                    nearest = candidate;
                }
            }
            throw new IllegalArgumentException("No tabulated CME of order " + order
                    + "; the nearest available order is " + nearest
                    + ". Use CME.getSupportedOrders() for the full list.");
        }
        return best;
    }

    /**
     * Gets the sorted list of CME orders (phase counts) present in the table.
     *
     * @return the available orders in increasing order
     */
    public static int[] getSupportedOrders() {
        List<iltcme.CmeEntry> params = iltcme.parameters();
        java.util.TreeSet<Integer> orders = new java.util.TreeSet<Integer>();
        for (iltcme.CmeEntry entry : params) {
            orders.add(2 * entry.n + 1);
        }
        int[] out = new int[orders.size()];
        int i = 0;
        for (Integer o : orders) {
            out[i++] = o;
        }
        return out;
    }

    /**
     * Gets the tabulated minimal SCV attained by a CME of the given order.
     *
     * @param order the number of phases
     * @return the minimal squared coefficient of variation
     */
    public static double getMinSCV(int order) {
        return tableEntry(order).cv2;
    }

    /**
     * Creates the lowest-order CME with the given mean whose SCV is at most the target.
     *
     * @param mean the target mean
     * @param scv  the target squared coefficient of variation, an upper bound: the lowest
     *             tabulated order whose minimal SCV does not exceed it is selected, so the
     *             result is at least as concentrated as requested
     * @return the CME of that order rescaled to the requested mean
     * @throws IllegalArgumentException if no tabulated order reaches the requested SCV
     */
    public static CME fitMeanAndSCV(double mean, double scv) {
        int bestOrder = -1;
        double minCv2 = Double.POSITIVE_INFINITY;
        int maxN = 0;
        for (iltcme.CmeEntry entry : iltcme.parameters()) {
            minCv2 = Math.min(minCv2, entry.cv2);
            maxN = Math.max(maxN, entry.n);
            if (entry.cv2 <= scv) {
                int order = 2 * entry.n + 1;
                if (bestOrder < 0 || order < bestOrder) {
                    bestOrder = order;
                }
            }
        }
        if (bestOrder < 0) {
            throw new IllegalArgumentException("No tabulated CME reaches SCV " + scv
                    + "; the most concentrated entry has SCV " + minCv2
                    + " at order " + (2 * maxN + 1));
        }
        return new CME(mean, bestOrder);
    }

    /**
     * Builds the unit-mean (alpha, A) matrix-exponential form of a CME.
     *
     * A = blkdiag(-mu1, mu1*[-1 -k*w; k*w -1], k=1..n) reproduces the exponential
     * envelope in its first phase and the k-th harmonic in its k-th 2x2 rotation block,
     * since expm(mu1*x*[-1 -kw; kw -1]) is exp(-mu1*x) times the rotation by k*w*mu1*x.
     * The entries of alpha follow by matching -alpha*expm(A*x)*A*e term by term: with
     * wk = k*w and d = 2*(1+wk^2),
     *
     * alpha_0 = c, alpha_{2k-1} = ((1+wk)*a_k - (1-wk)*b_k)/d,
     * alpha_{2k} = ((1-wk)*a_k + (1+wk)*b_k)/d.
     *
     * The result has unit mean and sums to one, as any (alpha, A) whose density
     * integrates to one must.
     *
     * @param order the number of phases
     * @return a two-element array holding the row vector alpha and the matrix A
     */
    public static Matrix[] representation(int order) {
        iltcme.CmeEntry entry = tableEntry(order);
        double c = entry.c;
        double mu1 = entry.mu1;
        double w = entry.omega;
        int n = entry.n;

        int size = 2 * n + 1;
        Matrix alpha = new Matrix(1, size);
        Matrix A = new Matrix(size, size);
        A.set(0, 0, -mu1);
        alpha.set(0, 0, c);
        for (int k = 1; k <= n; k++) {
            int i = 2 * k - 1;
            double wk = k * w;
            double ak = entry.a.get(k - 1);
            double bk = entry.b.get(k - 1);
            A.set(i, i, -mu1);
            A.set(i, i + 1, -wk * mu1);
            A.set(i + 1, i, wk * mu1);
            A.set(i + 1, i + 1, -mu1);
            double d = 2.0 * (1.0 + wk * wk);
            alpha.set(0, i, ((1.0 + wk) * ak - (1.0 - wk) * bk) / d);
            alpha.set(0, i + 1, ((1.0 - wk) * ak + (1.0 + wk) * bk) / d);
        }

        // see _kb/01-model-classes.md (CME sections) for rationale
        double sum = 0.0;
        for (int i = 0; i < size; i++) {
            sum += alpha.get(0, i);
        }
        alpha = alpha.scale(1.0 / sum);

        return new Matrix[]{alpha, A};
    }
}
