/**
 * @file Upper and lower bounds for K-way Fork-Join response time
 *
 * Computes pessimistic (upper) and optimistic (lower) bounds for the mean response
 * time of a K-way Fork-Join queueing system with Poisson arrivals and exponential
 * service times.
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import java.util.Objects;

public final class FJ_bounds {
    private FJ_bounds() {}

    /**
     * Result container for Fork-Join bounds.
     */
    public static final class FJBoundsResult {
        public final double Rmax;
        public final double Rmin;

        public FJBoundsResult(double Rmax, double Rmin) {
            this.Rmax = Rmax;
            this.Rmin = Rmin;
        }

        public double getRmax() { return Rmax; }
        public double getRmin() { return Rmin; }

        public double component1() { return Rmax; }
        public double component2() { return Rmin; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof FJBoundsResult)) return false;
            FJBoundsResult that = (FJBoundsResult) o;
            return Double.compare(that.Rmax, Rmax) == 0 && Double.compare(that.Rmin, Rmin) == 0;
        }

        @Override
        public int hashCode() { return Objects.hash(Rmax, Rmin); }

        @Override
        public String toString() {
            return "FJBoundsResult(Rmax=" + Rmax + ", Rmin=" + Rmin + ")";
        }
    }

    /**
     * Compute upper and lower bounds for K-way F/J response time.
     *
     * @param K      Number of parallel servers (positive integer)
     * @param lambda Arrival rate
     * @param mu     Service rate (mu &gt; lambda for stability)
     * @return FJBoundsResult containing upper and lower bounds
     * @throws IllegalArgumentException if K &lt; 1 or system is unstable
     */
    public static FJBoundsResult fj_bounds(int K, double lambda, double mu) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }

        double rho = lambda / mu;
        if (rho >= 1.0) {
            throw new IllegalArgumentException(String.format(
                    "System is unstable: rho = lambda/mu = %.4f >= 1. Require lambda < mu.", rho));
        }

        double H_K = FJ_harmonic.fj_harmonic(K);

        // Upper bound (Eq. 1)
        double Rmax = H_K / (mu * (1 - rho));

        // Lower bound (Eq. 2)
        double S_K = 0.0;
        for (int j = 1; j <= K; j++) {
            S_K += (1.0 / j) * (rho / (j - rho));
        }
        double Rmin = (1.0 / mu) * (H_K + S_K);

        return new FJBoundsResult(Rmax, Rmin);
    }
}
