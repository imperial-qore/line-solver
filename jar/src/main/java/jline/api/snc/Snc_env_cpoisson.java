/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * MGF arrival envelope of a compound Poisson flow with Exp job sizes.
 *
 * <p>Jobs arrive Poisson at rate lambda and each carries an Exp(mu) amount of
 * work, so {@code rho(theta) = lambda/(mu-theta)} for {@code 0 &lt; theta &lt; mu}
 * and the burst is zero. Fed to a constant-rate server of rate mu
 * ({@link Snc_srv_rate}) this is the network calculus model of the M/M/1 queue
 * in units of WORK; the delay bound then decays at rate mu-lambda, the exact
 * asymptotic decay rate of the M/M/1 waiting time.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_env_cpoisson.m.</p>
 */
public final class Snc_env_cpoisson {
    private Snc_env_cpoisson() {}

    /**
     * @param lambda job arrival rate, jobs per slot
     * @param mu     rate of the Exp job size, so the mean work per job is 1/mu
     * @param theta  Chernoff parameter, theta &gt; 0
     * @return {sigma, rho}, with rho = Infinity when theta &gt;= mu
     */
    public static double[] snc_env_cpoisson(double lambda, double mu, double theta) {
        if (lambda < 0 || mu <= 0) {
            throw new IllegalArgumentException(
                    "snc_env_cpoisson: lambda must be nonnegative and mu positive, got "
                            + lambda + ", " + mu);
        }
        if (theta <= 0) {
            throw new IllegalArgumentException("snc_env_cpoisson: theta must be positive, got " + theta);
        }
        if (theta >= mu) {
            // the job-size MGF diverges, no envelope at this theta
            return new double[] {0.0, Double.POSITIVE_INFINITY};
        }
        return new double[] {0.0, lambda / (mu - theta)};
    }

    /**
     * @param lambda job arrival rate
     * @param mu     rate of the Exp job size
     * @return the envelope as a function of theta
     */
    public static SncEnvelope of(final double lambda, final double mu) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                return snc_env_cpoisson(lambda, mu, theta);
            }
        };
    }
}
