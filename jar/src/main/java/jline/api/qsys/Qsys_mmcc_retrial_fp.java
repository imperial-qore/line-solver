/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

/**
 * Fixed-point approximation for the M/M/c/c retrial queue.
 *
 * <p>Port of {@code matlab/src/api/qsys/qsys_mmcc_retrial_fp.m}. Customers
 * arrive at rate lambda to a system of c servers each of rate mu, with no
 * waiting room; a blocked customer joins an orbit and retries. Under the
 * assumption that the retrial rate is small relative to the service rate, the
 * total arrival flow (fresh plus retrial) is approximated by a Poisson process
 * of rate lambda + r, with r the solution of
 *
 * <pre>
 *   r = (lambda + r) B(lambda/mu + r/mu, c),
 * </pre>
 *
 * B(a, c) being the Erlang-B blocking probability at offered load a.
 *
 * <p>References: Cohen (1957); Phung-Duc, "Retrial Queueing Models: A Survey on
 * Theory and Applications", 2019, Eq. (1).
 */
public final class Qsys_mmcc_retrial_fp {

    private Qsys_mmcc_retrial_fp() {
    }

    /** Blocking probability, orbit-induced rate and iteration count. */
    public static final class Result {
        /** Fraction of arrivals blocked, hence retried. */
        public final double blocProb;
        /** Additional arrival rate contributed by retrials. */
        public final double r;
        /** Iterations taken to converge. */
        public final int niter;

        Result(double blocProb, double r, int niter) {
            this.blocProb = blocProb;
            this.r = r;
            this.niter = niter;
        }
    }

    /** The fixed point at the reference defaults tol = 1e-10, maxiter = 10000. */
    public static Result qsys_mmcc_retrial_fp(double lambda, double mu, int c) {
        return qsys_mmcc_retrial_fp(lambda, mu, c, 1e-10, 10000);
    }

    /**
     * The retrial fixed point.
     *
     * @param lambda  fresh arrival rate
     * @param mu      service rate per server
     * @param c       number of servers, equal to the capacity
     * @param tol     convergence tolerance on r
     * @param maxiter iteration cap
     */
    public static Result qsys_mmcc_retrial_fp(double lambda, double mu, int c, double tol,
                                              int maxiter) {
        double r = 0.0;
        int niter = 0;
        for (int iter = 1; iter <= maxiter; iter++) {
            niter = iter;
            double a = (lambda + r) / mu;
            double b = erlangB(a, c);
            double rNew = (lambda + r) * b;
            if (Math.abs(rNew - r) < tol) {
                r = rNew;
                break;
            }
            r = rNew;
        }
        return new Result(erlangB((lambda + r) / mu, c), r, niter);
    }

    /** Erlang-B by the recursion B_i = a B_{i-1}/(i + a B_{i-1}), numerically stable. */
    public static double erlangB(double a, int c) {
        double B = 1.0;
        for (int i = 1; i <= c; i++) {
            B = a * B / (i + a * B);
        }
        return B;
    }
}
