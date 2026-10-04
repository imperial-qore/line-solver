/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import odesolver.LSODA;

import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;

/**
 * Sensitivity of the transient distribution of a CTMC to a scalar parameter.
 *
 * <p>Port of {@code matlab/src/api/mc/ctmc_transient_sens.m}, twin of
 * {@code cpp/include/line/api/mc/ctmc_transient_sens.h} and of the native
 * Python {@code api.mc.ctmc_transient_sens}. Differentiating the forward
 * equations d pi(t)/dt = pi(t) Q with respect to theta, with an initial vector
 * that does not depend on theta, gives Trivedi and Bobbio (2017), Eq. (9.82),
 *
 * <pre>
 *   d/dt (dpi/dtheta) = (dpi/dtheta) Q + pi (dQ/dtheta),   dpi(0)/dtheta = 0.
 * </pre>
 *
 * <p>The sensitivity equation is DRIVEN by pi(t), so the two cannot be advanced
 * separately: state and sensitivity are integrated as ONE augmented system of
 * size 2n, which is also what keeps them consistent at every returned time
 * point. The integrator is the same LSODA that {@link Ctmc_transient} uses, so
 * the accepted grid is the JAR's own and need not coincide with the reference's
 * ode23 grid; the trajectories agree, the abscissae need not.
 */
public final class Ctmc_transient_sens {

    private Ctmc_transient_sens() {
    }

    /** Distribution and its sensitivity on the integrator's accepted grid. */
    public static final class Result {
        /** Accepted time points, the first being t0. */
        public double[] t;
        /** length(t) x n, the distribution. */
        public Matrix pi;
        /** length(t) x n, its derivative with respect to theta. */
        public Matrix dpi;
    }

    /** From the uniform initial distribution over [0, t1]. */
    public static Result ctmc_transient_sens(Matrix Q, Matrix dQ, double t1) {
        int n = Q.getNumRows();
        Matrix pi0 = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            pi0.set(0, i, 1.0 / n);
        }
        return ctmc_transient_sens(Q, dQ, pi0, 0.0, t1);
    }

    /** From PI0 over [0, t1]. */
    public static Result ctmc_transient_sens(Matrix Q, Matrix dQ, Matrix pi0, double t1) {
        return ctmc_transient_sens(Q, dQ, pi0, 0.0, t1);
    }

    /**
     * @param Q   generator, n x n
     * @param dQ  derivative of the generator with respect to theta, same size
     * @param pi0 initial distribution, 1 x n
     * @param t0  initial time
     * @param t1  final time
     */
    public static Result ctmc_transient_sens(final Matrix Q, final Matrix dQ, Matrix pi0,
                                             double t0, double t1) {
        final int n = Q.getNumRows();
        if (Q.getNumCols() != n) {
            line_error("ctmc_transient_sens", "the generator is not square.");
        }
        if (dQ.getNumRows() != n || dQ.getNumCols() != n) {
            line_error("ctmc_transient_sens", "dQ must have the same size as Q.");
        }
        if (pi0.length() != n) {
            line_error("ctmc_transient_sens", "pi0 has the wrong length.");
        }

        // Augmented state v = [pi, dpi], with dpi(0) = 0 since pi(0) does not
        // depend on theta
        double[] v0 = new double[2 * n];
        for (int i = 0; i < n; i++) {
            v0[i] = pi0.get(i);
        }

        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 2 * n;
            }

            @Override
            public void computeDerivatives(double t, double[] v, double[] dv)
                    throws MaxCountExceededException, DimensionMismatchException {
                for (int j = 0; j < n; j++) {
                    double a = 0.0;
                    double b = 0.0;
                    for (int i = 0; i < n; i++) {
                        a += v[i] * Q.get(i, j);
                        b += v[n + i] * Q.get(i, j) + v[i] * dQ.get(i, j);
                    }
                    dv[j] = a;
                    dv[n + j] = b;
                }
            }
        };

        LSODA lsoda = new LSODA(0.0, 0.0, 1.0e-8, 1.0e-8, 12, 5);
        double[] yend = new double[2 * n];
        lsoda.integrate(ode, t0, v0, t1, yend);

        List<double[]> traj = new ArrayList<double[]>();
        for (Double[] row : lsoda.getYvec()) {
            double[] arr = new double[row.length];
            for (int k = 0; k < row.length; k++) {
                arr[k] = row[k].doubleValue();
            }
            traj.add(arr);
        }
        int m = traj.size();
        Result r = new Result();
        r.t = new double[m];
        for (int k = 0; k < m; k++) {
            r.t[k] = lsoda.getTvec().get(k).doubleValue();
        }
        r.pi = new Matrix(m, n);
        r.dpi = new Matrix(m, n);
        for (int k = 0; k < m; k++) {
            double[] row = traj.get(k);
            for (int j = 0; j < n; j++) {
                r.pi.set(k, j, row[j]);
                r.dpi.set(k, j, row[n + j]);
            }
        }
        return r;
    }
}
