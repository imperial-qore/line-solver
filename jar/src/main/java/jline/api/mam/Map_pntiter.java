/**
 * @file Markovian Arrival Process point process probability computation using iteration
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.util.Pair;

import jline.util.matrix.Matrix;

public final class Map_pntiter {
    private Map_pntiter() {}

    /**
     * Probability of having exactly na arrivals within time interval t using iterative method.
     */
    public static Matrix map_pntiter(Matrix[] MAP, int na, double t, Integer M) {
        if (MAP.length != 2) {
            throw new IllegalArgumentException("MAP must contain exactly 2 matrices [D0, D1]");
        }

        double mapMean = Map_mean.map_mean(MAP[0], MAP[1]);
        int actualM;
        if (M != null) {
            actualM = M.intValue();
        } else {
            actualM = (int) Math.ceil(Math.log(t * 100 / mapMean) / Math.log(2.0));
        }

        if (actualM < 0) {
            Pair<Matrix, List<Matrix>> result = map_pntbisect(MAP, na, t);
            return result.getFirst();
        } else {
            Pair<Matrix, List<Matrix>> initialResult = map_pntbisect(MAP, na, t / Math.pow(2.0, actualM));
            List<Matrix> P = new ArrayList<Matrix>(initialResult.getSecond());

            for (int i = 1; i <= actualM; i++) {
                List<Matrix> Pold = new ArrayList<Matrix>(P.size());
                for (int x = 0; x < P.size(); x++) {
                    Pold.add(P.get(x).copy());
                }
                for (int n = 0; n <= na; n++) {
                    Matrix newP = P.get(n).copy();
                    newP.fill(0.0);

                    for (int j = 0; j <= n; j++) {
                        newP = newP.add(Pold.get(j).mult(Pold.get(n - j)));
                    }
                    P.set(n, newP);
                }
            }

            return P.get(na);
        }
    }

    public static Matrix map_pntiter(Matrix[] MAP, int na, double t) {
        return map_pntiter(MAP, na, t, null);
    }

    /**
     * Helper function implementing the bisection algorithm.
     */
    private static Pair<Matrix, List<Matrix>> map_pntbisect(Matrix[] MAP, int na, double t) {
        Matrix D0 = MAP[0];
        Matrix D1 = MAP[1];

        double tau = 0.0;
        for (int i = 0; i < D0.getNumRows(); i++) {
            double diagVal = -D0.get(i, i);
            if (diagVal > tau) tau = diagVal;
        }

        int N = findN(tau, t);
        Matrix[][] V = new Matrix[na + 1][N + 1];
        Matrix[] P = new Matrix[na + 1];

        Matrix I = Matrix.eye(D0.getNumRows());
        Matrix K = D0.scale(1.0 / tau).add(I);
        Matrix K1 = D1.scale(1.0 / tau);

        for (int n = 0; n <= na; n++) {
            P[n] = new Matrix(D0.getNumRows(), D0.getNumCols());
            for (int k = 0; k <= N; k++) {
                V[n][k] = new Matrix(D0.getNumRows(), D0.getNumCols());
            }
        }

        // Uniformization: V(n,k) holds the paths making exactly n arrivals in k
        // STEPS of the uniformized chain, so V(0,k) = V(0,k-1)*K is the
        // no-arrival path of length k and carries real mass for every k. It was
        // previously left at zero past k = 0, which truncated P_0 to its k = 0
        // term. A POISSON PROCESS CANNOT SEE THAT: there K = D0/tau + I = 0, so
        // V(0,k) really is zero for k >= 1 and the omission is invisible. The
        // identity that does see it is P_0(t) = expm(D0*t).
        V[0][0] = I.copy();
        P[0] = V[0][0].scale(br(tau, t, 0));
        for (int k = 1; k <= N; k++) {
            V[0][k] = V[0][k - 1].mult(K);
            P[0] = P[0].add(V[0][k].scale(br(tau, t, k)));
        }

        for (int n = 1; n <= na; n++) {
            V[n][0] = new Matrix(D0.getNumRows(), D0.getNumCols());
            P[n] = new Matrix(D0.getNumRows(), D0.getNumCols());

            for (int k = 1; k <= N; k++) {
                V[n][k] = V[n][k - 1].mult(K).add(V[n - 1][k - 1].mult(K1));
                P[n] = P[n].add(V[n][k].scale(br(tau, t, k)));
            }
        }

        return new Pair<Matrix, List<Matrix>>(P[na], Arrays.asList(P));
    }

    /**
     * Find optimal N parameter for the bisection algorithm.
     */
    private static int findN(double tau, double t) {
        double epsilon = Double.MIN_VALUE;
        int Nmax = 100;

        for (int N = 1; N <= Nmax; N++) {
            double S = 0.0;
            for (int n = N + 1; n <= Nmax; n++) {
                S += br(tau, t, n);
            }
            if (S < epsilon) {
                return N;
            }
        }

        return Nmax;
    }

    /**
     * Bernoulli coefficient: exp(-tau*t) * (tau*t)^r / r!
     */
    private static double br(double tau, double t, int r) {
        double tauT = tau * t;
        return Math.exp(-tauT) * Math.pow(tauT, r) / factorial(r);
    }

    /**
     * Factorial helper function.
     */
    private static double factorial(int n) {
        if (n < 0) throw new IllegalArgumentException("Factorial of negative number");
        if (n == 0 || n == 1) return 1.0;

        double result = 1.0;
        for (int i = 2; i <= n; i++) {
            result *= i;
        }
        return result;
    }
}
