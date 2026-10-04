/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

import jline.util.matrix.Matrix;

/**
 * MGF arrival envelope of a MAP/MMPP flow with unit-size jobs.
 *
 * <p>For a Markovian arrival process (D0,D1) counting N(0,t) unit-work jobs,
 * {@code E[exp(theta*N(0,t))] = pi*expm((D0+D1*exp(theta))*t)*1}. With lstar
 * the eigenvalue of maximal real part of {@code A(theta)=D0+D1*e^theta} and
 * v &gt; 0 its right Perron eigenvector, bounding {@code 1 <= v/min(v)}
 * entrywise gives</p>
 *
 * <pre>
 *   rho(theta)   = lstar/theta,
 *   sigma(theta) = log(max(v)/min(v))/theta,
 * </pre>
 *
 * <p>the standard exponential-form envelope of a Markov-modulated source. The
 * burst term is what the modulating chain contributes: 0 for a one-phase MAP,
 * where this reproduces {@link Snc_env_poisson} exactly, and positive for an
 * MMPP.</p>
 *
 * <p>THE PERRON PAIR IS COMPUTED BY POWER ITERATION ON THE SHIFTED MATRIX
 * {@code A+cI}, not by a general eigensolver. A(theta) is essentially
 * nonnegative, so the shift makes it nonnegative with a positive diagonal,
 * hence primitive whenever the MAP is irreducible, and the iteration converges
 * to the pair the bound needs without asking a general solver which of its
 * eigenvectors is the positive one. The MATLAB reference uses {@code eig} and
 * agrees to machine precision.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_env_map.m. Reference: C.-S. Chang,
 * Performance Guarantees in Communication Networks, Springer 2000, Ch. 7.</p>
 */
public final class Snc_env_map {
    private Snc_env_map() {}

    private static final int MAX_ITER = 100000;
    private static final double TOL = 1e-14;

    /**
     * @param D0    hidden-transition generator block of the MAP
     * @param D1    arrival-transition block of the MAP
     * @param theta Chernoff parameter, theta &gt; 0
     * @return {sigma, rho}
     */
    public static double[] snc_env_map(Matrix D0, Matrix D1, double theta) {
        if (theta <= 0) {
            throw new IllegalArgumentException("snc_env_map: theta must be positive, got " + theta);
        }
        final int n = D0.getNumRows();
        if (D0.getNumCols() != n || D1.getNumRows() != n || D1.getNumCols() != n) {
            throw new IllegalArgumentException("snc_env_map: D0 and D1 must be square and of equal size.");
        }
        double[][] A = new double[n][n];
        double etheta = Math.exp(theta);
        double shift = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = D0.get(i, j) + D1.get(i, j) * etheta;
            }
            shift = Math.max(shift, -A[i][i]);
        }
        // A is essentially nonnegative: shifting by its most negative diagonal
        // entry makes it nonnegative with a positive diagonal, hence primitive.
        shift += 1.0;
        for (int i = 0; i < n; i++) {
            A[i][i] += shift;
        }

        double[] v = new double[n];
        for (int i = 0; i < n; i++) {
            v[i] = 1.0;
        }
        double lambdaShift = 0.0;
        for (int it = 0; it < MAX_ITER; it++) {
            double[] w = new double[n];
            for (int i = 0; i < n; i++) {
                double acc = 0.0;
                for (int j = 0; j < n; j++) {
                    acc += A[i][j] * v[j];
                }
                w[i] = acc;
            }
            double norm = 0.0;
            for (int i = 0; i < n; i++) {
                norm = Math.max(norm, Math.abs(w[i]));
            }
            if (!(norm > 0)) {
                throw new RuntimeException("snc_env_map: MAP is not irreducible: the Perron "
                        + "eigenvector is not positive.");
            }
            double delta = 0.0;
            for (int i = 0; i < n; i++) {
                w[i] /= norm;
                delta = Math.max(delta, Math.abs(w[i] - v[i]));
            }
            v = w;
            lambdaShift = norm;
            if (delta < TOL) {
                break;
            }
        }
        double vmin = Double.POSITIVE_INFINITY;
        double vmax = 0.0;
        for (int i = 0; i < n; i++) {
            vmin = Math.min(vmin, v[i]);
            vmax = Math.max(vmax, v[i]);
        }
        if (!(vmin > 0)) {
            throw new RuntimeException("snc_env_map: MAP is not irreducible: the Perron "
                    + "eigenvector is not positive.");
        }
        double lstar = lambdaShift - shift;
        return new double[] {Math.log(vmax / vmin) / theta, lstar / theta};
    }

    /**
     * @param D0 hidden-transition generator block
     * @param D1 arrival-transition block
     * @return the envelope as a function of theta
     */
    public static SncEnvelope of(final Matrix D0, final Matrix D1) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                return snc_env_map(D0, D1, theta);
            }
        };
    }
}
