/**
 * Harel et al. throughput bounds for closed queueing networks.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_harel_bounds {
    private Pfqn_harel_bounds() {}

    public static Ret.pfqnHarelBounds pfqn_harel_bounds(Matrix rho, int N) {
        return pfqn_harel_bounds(rho, N, 0.0, null);
    }

    public static Ret.pfqnHarelBounds pfqn_harel_bounds(Matrix rho, int N, double Z) {
        return pfqn_harel_bounds(rho, N, Z, null);
    }

    public static Ret.pfqnHarelBounds pfqn_harel_bounds(Matrix rho, int N, double Z, Integer maxUB) {
        if (Z != 0.0) {
            throw new IllegalArgumentException(
                    "Harel bounds are only valid for networks with zero think time. The provided think time Z=" + Z + " is nonzero.");
        }
        if (N < 1) throw new IllegalArgumentException("Population N must be at least 1, got N=" + N);
        int k = rho.getNumRows();
        if (k < 1) throw new IllegalArgumentException("Loading vector rho must have at least 1 element");
        for (int i = 0; i < k; i++) {
            if (rho.get(i, 0) <= 0) {
                throw new IllegalArgumentException("All loading factors must be positive, rho[" + i + "]=" + rho.get(i, 0));
            }
        }
        int effectiveMaxUB = (maxUB != null) ? maxUB : Math.min(N, 7);
        if (effectiveMaxUB > 7) {
            throw new IllegalArgumentException(
                    "Upper bounds UB(n) can only be computed for n <= 7 due to G(n) polynomial limitations. Requested maxUB=" + effectiveMaxUB);
        }
        double[] A = computeAValues(rho, N);
        double LB = computeLowerBound(A, N);
        double[] TH = computeThroughputValues(A, effectiveMaxUB);
        double[] UB = computeUpperBounds(A, N, TH, effectiveMaxUB);
        return new Ret.pfqnHarelBounds(LB, UB, TH, N, k, effectiveMaxUB);
    }

    public static double pfqn_harel_lb(Matrix rho, int N) {
        return pfqn_harel_lb(rho, N, 0.0);
    }

    public static double pfqn_harel_lb(Matrix rho, int N, double Z) {
        if (Z != 0.0) {
            throw new IllegalArgumentException("Harel LB is only valid for networks with zero think time. The provided think time Z=" + Z + " is nonzero.");
        }
        if (N < 1) throw new IllegalArgumentException("Population N must be at least 1, got N=" + N);
        double[] A = computeAValues(rho, N);
        return computeLowerBound(A, N);
    }

    public static double pfqn_harel_ub(Matrix rho, int N, int n) {
        return pfqn_harel_ub(rho, N, n, 0.0);
    }

    public static double pfqn_harel_ub(Matrix rho, int N, int n, double Z) {
        if (Z != 0.0) {
            throw new IllegalArgumentException("Harel UB is only valid for networks with zero think time. The provided think time Z=" + Z + " is nonzero.");
        }
        if (N < 1) throw new IllegalArgumentException("Population N must be at least 1, got N=" + N);
        if (n < 2) throw new IllegalArgumentException("n must be at least 2 for upper bound UB(n), got n=" + n);
        if (n > N) throw new IllegalArgumentException("n cannot exceed N for UB(n), got n=" + n + ", N=" + N);
        if (n > 7) throw new IllegalArgumentException("UB(n) can only be computed for n <= 7 due to G(n) polynomial limitations, got n=" + n);

        double[] A = computeAValues(rho, n);
        double Gn = computeG(A, n);
        double Gn_1 = computeG(A, n - 1);
        double THn = Gn_1 / Gn;
        double A1 = A[1];
        double nOverTH = (double) n / THn;
        double denominator = A1 + ((N - 1.0) / (n - 1.0)) * (nOverTH - A1);
        return (double) N / denominator;
    }

    private static double[] computeAValues(Matrix rho, int maxPower) {
        int k = rho.getNumRows();
        double[] A = new double[maxPower + 1];
        for (int i = 1; i <= maxPower; i++) {
            double sum = 0.0;
            for (int j = 0; j < k; j++) sum += Math.pow(rho.get(j, 0), i);
            A[i] = sum;
        }
        return A;
    }

    private static double computeLowerBound(double[] A, int N) {
        if (N == 1) return 1.0 / A[1];
        double A1 = A[1];
        double AN = A[N];
        double exponent = 1.0 / (N - 1);
        return (double) N / (A1 + (N - 1) * Math.pow(AN / A1, exponent));
    }

    private static double computeG(double[] A, int n) {
        switch (n) {
            case 0: return 1.0;
            case 1: return A[1];
            case 2: return (A[1] * A[1] + A[2]) / 2.0;
            case 3: {
                double A1 = A[1], A2 = A[2], A3 = A[3];
                return (A1 * A1 * A1 + 3 * A2 * A1 + 2 * A3) / 6.0;
            }
            case 4: {
                double A1 = A[1], A2 = A[2], A3 = A[3], A4 = A[4];
                return (Math.pow(A1, 4.0) + 6 * A2 * A1 * A1 + 8 * A3 * A1 + 3 * A2 * A2 + 6 * A4) / 24.0;
            }
            case 5: {
                double A1 = A[1], A2 = A[2], A3 = A[3], A4 = A[4], A5 = A[5];
                return (Math.pow(A1, 5.0) + 10 * A2 * Math.pow(A1, 3.0) + 20 * A3 * A1 * A1
                        + 15 * A1 * A2 * A2 + 30 * A4 * A1 + 20 * A2 * A3 + 24 * A5) / 120.0;
            }
            case 6: {
                double A1 = A[1], A2 = A[2], A3 = A[3], A4 = A[4], A5 = A[5], A6 = A[6];
                return (Math.pow(A1, 6.0) + 15 * A2 * Math.pow(A1, 4.0) + 40 * A3 * Math.pow(A1, 3.0)
                        + 45 * A1 * A1 * A2 * A2 + 90 * A4 * A1 * A1 + 120 * A1 * A2 * A3
                        + 144 * A5 * A1 + 15 * Math.pow(A2, 3.0) + 90 * A2 * A4 + 40 * A3 * A3 + 120 * A6) / 720.0;
            }
            case 7: {
                double A1 = A[1], A2 = A[2], A3 = A[3], A4 = A[4], A5 = A[5], A6 = A[6], A7 = A[7];
                return (Math.pow(A1, 7.0) + 21 * A2 * Math.pow(A1, 5.0) + 70 * A3 * Math.pow(A1, 4.0)
                        + 105 * Math.pow(A1, 3.0) * A2 * A2 + 210 * A4 * Math.pow(A1, 3.0) + 504 * A5 * A1 * A1
                        + 105 * A1 * Math.pow(A2, 3.0) + 280 * A1 * A3 * A3 + 210 * A2 * A2 * A3
                        + 504 * A2 * A5 + 420 * A3 * A4 + 840 * A6 * A1 + 630 * A1 * A2 * A4
                        + 420 * A1 * A1 * A2 * A3 + 720 * A7) / 5040.0;
            }
            default:
                throw new IllegalArgumentException("G(n) polynomial not available for n=" + n + " (max supported: 7)");
        }
    }

    private static double[] computeThroughputValues(double[] A, int maxN) {
        double[] TH = new double[maxN + 1];
        for (int n = 1; n <= maxN; n++) {
            double Gn = computeG(A, n);
            double Gn_1 = computeG(A, n - 1);
            TH[n] = Gn_1 / Gn;
        }
        return TH;
    }

    private static double[] computeUpperBounds(double[] A, int N, double[] TH, int maxN) {
        double[] UB = new double[maxN + 1];
        double A1 = A[1];
        for (int n = 2; n <= maxN; n++) {
            if (n > TH.length - 1) break;
            double nOverTH = (double) n / TH[n];
            double denominator = A1 + ((N - 1.0) / (n - 1.0)) * (nOverTH - A1);
            UB[n] = (double) N / denominator;
        }
        return UB;
    }
}
