package jline.api.qsys;

import java.util.HashMap;
import java.util.Map;
import java.util.function.DoubleUnaryOperator;

import jline.api.mc.Dtmc_makestochastic;
import jline.api.mc.Dtmc_solve;
import jline.util.matrix.Matrix;

public final class Qsys_mg1k_loss {
    private Qsys_mg1k_loss() {}

    private static final Map<Integer, Double> factorialCache = new HashMap<Integer, Double>();

    /**
     * Exact M/G/1/K loss probability via the Markov chain embedded at
     * service-start epochs (transform-free analysis in the spirit of
     * Niu-Cooper).
     *
     * <p>State: number of customers waiting in the queue immediately after a
     * service start, q in {0,...,K-2}. With a_j = P(j Poisson arrivals during
     * a service time): from q=0, if no arrival occurs during the service the
     * system empties and the next service starts with the next arrival
     * (q'=0), so both a_0 and a_1 lead to q'=0 and j&gt;=2 arrivals lead to
     * q'=j-1; from q&gt;=1, q'=q-1+j with arrivals beyond the free capacity
     * lost (aggregated in the last column). The loss probability follows from
     * the renewal-reward argument P_loss = 1 - 1/(rho + sigma_0*a_0), with
     * sigma the stationary distribution at service-start epochs.
     *
     * <p>Reference: Niu, Cooper. Transform-Free Analysis of M/G/1/K and
     * Related Queues. Mathematics of Operations Research 18(2), 1993, 486-510.
     */
    public static HashMap<String, Object> qsys_mg1k_loss(double lambda, DoubleUnaryOperator svc_density, int K) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        // Choose the integration horizon so that it covers the service-time
        // distribution mass, then compute the mean service time
        final DoubleUnaryOperator svcDensRef = svc_density;
        DoubleUnaryOperator tTimesDensity = new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double t) {
                return t * svcDensRef.applyAsDouble(t);
            }
        };
        double T = 1.0 / lambda;
        for (int it = 0; it < 60; it++) {
            double mass = simpson(svcDensRef, 0.0, T, 65536);
            if (mass >= 1.0 - 1e-10) {
                break;
            }
            T = 2.0 * T;
        }
        double meanS = simpson(tTimesDensity, 0.0, T, 200000);
        double mu = 1.0 / meanS;

        // Integration horizon covering both the service-time mass and the
        // Poisson weight window
        double tmax = 100.0 * Math.max(meanS, 1.0 / lambda);

        // Arrival probabilities a_j = P(j arrivals during a service time)
        double[] a = new double[Math.max(K, 2)];
        for (int j = 0; j < a.length; j++) {
            double factj = factorial(j);
            final int jFinal = j;
            DoubleUnaryOperator integrand = new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double t) {
                    return Math.exp(-lambda * t) * Math.pow(lambda * t, jFinal) * svcDensRef.applyAsDouble(t);
                }
            };
            a[j] = simpson(integrand, 0.0, tmax, 200000) / factj;

            if (a[j] < 1e-12) {
                break;
            }
        }

        // Embedded chain at service-start epochs, states q=0..K-2
        int n = K - 1;
        Matrix P = new Matrix(n, n);

        // row 0 (q=0): idle period after an empty departure epoch
        P.set(0, 0, a[0] + (1 < a.length ? a[1] : 0.0));
        for (int i = 1; i <= K - 3; i++) {
            P.set(0, i, a[i + 1]);
        }
        P.set(0, n - 1, 1.0 - rowSumFirst(P, 0, n - 1));

        // row 1 (q=1): q' = number of arrivals during the service (capped)
        if (n >= 2) {
            for (int i = 0; i <= K - 3; i++) {
                P.set(1, i, a[i]);
            }
            P.set(1, n - 1, 1.0 - rowSumFirst(P, 1, n - 1));
        }

        // rows j>=2 (q=j): q' = q-1+arrivals (capped)
        for (int j = 2; j < n; j++) {
            for (int i = j - 1; i <= K - 3; i++) {
                P.set(j, i, a[i - j + 1]);
            }
            P.set(j, n - 1, 1.0 - rowSumFirst(P, j, n - 1));
        }

        // Make stochastic and solve
        Matrix Pstoch = Dtmc_makestochastic.dtmc_makestochastic(P);
        Matrix sigma = Dtmc_solve.dtmc_solve(Pstoch);

        double rho = lambda / mu;
        double lossprob = 1.0 - 1.0 / (sigma.get(0) * a[0] + rho);

        result.put("lossprob", lossprob);
        result.put("rho", rho);

        return result;
    }

    /** Sum of the first {@code count} entries of row {@code row}. */
    private static double rowSumFirst(Matrix P, int row, int count) {
        double s = 0.0;
        for (int i = 0; i < count; i++) {
            s += P.get(row, i);
        }
        return s;
    }

    /**
     * Composite Simpson quadrature with n subintervals (n rounded up to even).
     */
    private static double simpson(DoubleUnaryOperator f, double a, double b, int n) {
        if (n % 2 == 1) {
            n++;
        }
        double h = (b - a) / n;
        double sum = f.applyAsDouble(a) + f.applyAsDouble(b);
        for (int i = 1; i < n; i++) {
            double x = a + i * h;
            sum += f.applyAsDouble(x) * (i % 2 == 1 ? 4.0 : 2.0);
        }
        return sum * h / 3.0;
    }

    /**
     * Factorial function with memoization for efficiency
     */
    private static double factorial(int n) {
        if (n == 0 || n == 1) return 1.0;
        Double cached = factorialCache.get(n);
        if (cached != null) return cached;
        double result = n * factorial(n - 1);
        factorialCache.put(n, result);
        return result;
    }
}
