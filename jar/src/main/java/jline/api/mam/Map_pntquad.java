/**
 * @file Markovian Arrival Process point process probability computation using quadrature
 *
 * Computes MAP point process probabilities using ODE quadrature methods with Runge-Kutta integration.
 * Provides alternative numerical approach for arrival probability calculations with controlled accuracy.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Map_pntquad {
    private Map_pntquad() {}

    /**
     * Compute MAP point process probabilities using ODE quadrature method.
     * This is a simplified implementation using Runge-Kutta 4th order method.
     *
     * @param MAP MAP process as Matrix[] where MAP[0] = D0 and MAP[1] = D1
     * @param na  Maximum number of arrivals to consider
     * @param t   Time interval
     * @return Pair&lt;Matrix, Matrix&gt; containing (Pnt, S) where:
     *         Pnt: Probability matrix for exactly na arrivals
     *         S: Verification sum matrix (should sum to identity)
     */
    public static Pair<Matrix, Matrix> map_pntquad(Matrix[] MAP, int na, double t) {
        if (MAP.length != 2) {
            throw new IllegalArgumentException("MAP must contain exactly 2 matrices [D0, D1]");
        }

        Matrix D0 = MAP[0];
        Matrix D1 = MAP[1];
        int Ki = D0.getNumRows();

        if (D0.getNumCols() != Ki || D1.getNumRows() != Ki || D1.getNumCols() != Ki) {
            throw new IllegalArgumentException("D0 and D1 must be square matrices of the same size");
        }

        // Initialize state vector P0
        int stateSize = Ki * Ki * (na + 1);
        double[] P0 = new double[stateSize];

        // Set initial condition: P0(1:Ki^2) = reshape(eye(Ki), 1, Ki^2)
        Matrix eye = Matrix.eye(Ki);
        for (int i = 0; i < Ki; i++) {
            for (int j = 0; j < Ki; j++) {
                P0[i * Ki + j] = eye.get(i, j);
            }
        }

        // All other components start at zero (already initialized)

        // Solve ODE using 4th-order Runge-Kutta method
        double[] result = rungeKutta4(P0, t, Ki, na, D0, D1);

        // Extract results
        Matrix Pnt = new Matrix(Ki, Ki);
        Matrix S = new Matrix(Ki, Ki);

        // Extract all probability matrices and compute sum
        for (int n = 0; n <= na; n++) {
            int startIdx = n * Ki * Ki;
            Matrix Pn = new Matrix(Ki, Ki);

            for (int i = 0; i < Ki; i++) {
                for (int j = 0; j < Ki; j++) {
                    int idx = startIdx + i * Ki + j;
                    Pn.set(i, j, result[idx]);
                }
            }

            S.addEq(Pn);

            // Keep the last one (na arrivals) as Pnt
            if (n == na) {
                for (int i = 0; i < Ki; i++) {
                    for (int j = 0; j < Ki; j++) {
                        Pnt.set(i, j, Pn.get(i, j));
                    }
                }
            }
        }

        return new Pair<Matrix, Matrix>(Pnt, S);
    }

    /**
     * 4th-order Runge-Kutta method for solving the ODE system.
     */
    private static double[] rungeKutta4(double[] P0, double t, int Ki, int na, Matrix D0, Matrix D1) {
        double h = t / 1000.0; // Step size
        int steps = (int) (t / h);

        double[] P = P0.clone();
        double currentT = 0.0;

        for (int step = 0; step < steps; step++) {
            double[] k1 = pntOde(currentT, P, Ki, na, D0, D1);

            double[] P_k1 = new double[P.length];
            for (int i = 0; i < P.length; i++) P_k1[i] = P[i] + h * k1[i] / 2.0;
            double[] k2 = pntOde(currentT + h / 2.0, P_k1, Ki, na, D0, D1);

            double[] P_k2 = new double[P.length];
            for (int i = 0; i < P.length; i++) P_k2[i] = P[i] + h * k2[i] / 2.0;
            double[] k3 = pntOde(currentT + h / 2.0, P_k2, Ki, na, D0, D1);

            double[] P_k3 = new double[P.length];
            for (int i = 0; i < P.length; i++) P_k3[i] = P[i] + h * k3[i];
            double[] k4 = pntOde(currentT + h, P_k3, Ki, na, D0, D1);

            // Update P
            for (int i = 0; i < P.length; i++) {
                P[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
                // Ensure non-negativity
                if (P[i] < 0.0) P[i] = 0.0;
            }

            currentT += h;
        }

        return P;
    }

    /**
     * ODE right-hand side function.
     * Implements the differential equation system for the probability evolution.
     */
    private static double[] pntOde(double t, double[] P, int Ki, int na, Matrix D0, Matrix D1) {
        double[] dP = new double[P.length];

        // For n=0: dP(1:Ki^2) = reshape(reshape(P(1:Ki^2), Ki, Ki) * D0, Ki^2, 1)
        Matrix P0_matrix = new Matrix(Ki, Ki);
        for (int i = 0; i < Ki; i++) {
            for (int j = 0; j < Ki; j++) {
                P0_matrix.set(i, j, P[i * Ki + j]);
            }
        }

        Matrix dP0_matrix = P0_matrix.mult(D0);
        for (int i = 0; i < Ki; i++) {
            for (int j = 0; j < Ki; j++) {
                dP[i * Ki + j] = dP0_matrix.get(i, j);
            }
        }

        // For n=1 to na
        for (int n = 1; n <= na; n++) {
            int startIdx = n * Ki * Ki;
            int prevStartIdx = (n - 1) * Ki * Ki;

            // Extract Pn-1(t) and Pn(t)
            Matrix Pn_1t = new Matrix(Ki, Ki);
            Matrix Pnt = new Matrix(Ki, Ki);

            for (int i = 0; i < Ki; i++) {
                for (int j = 0; j < Ki; j++) {
                    Pn_1t.set(i, j, P[prevStartIdx + i * Ki + j]);
                    Pnt.set(i, j, P[startIdx + i * Ki + j]);
                }
            }

            // Compute derivative: dPn/dt = Pn * D0 + Pn-1 * D1
            Matrix dPn_matrix = Pnt.mult(D0).add(1.0, Pn_1t.mult(D1));

            // Store back in dP array
            for (int i = 0; i < Ki; i++) {
                for (int j = 0; j < Ki; j++) {
                    dP[startIdx + i * Ki + j] = dPn_matrix.get(i, j);
                }
            }
        }

        return dP;
    }
}
