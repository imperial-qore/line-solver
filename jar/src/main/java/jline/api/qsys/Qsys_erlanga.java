/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

/**
 * Exact analysis of the Erlang A model M/M/s/r+M.
 *
 * <p>Poisson arrivals at rate lambda, exponential service of rate mu at each of
 * s servers and exponential patience of rate theta, so a waiting customer
 * abandons after an exponential time of mean 1/theta. The number in system is
 * the birth-and-death process with birth rate lambda and death rate
 * min(k,s)*mu + (k-s)^+ *theta, so every measure is EXACT: this is the special
 * case in which the state-dependent Markovian approximation of
 * {@link Qsys_mgisrgi_whitt} reproduces the model rather than approximating it
 * (eq. 7.12 of the reference). theta &gt; 0 makes the model ergodic at every load,
 * including lambda above s*mu.
 *
 * <p>Port of MATLAB qsys_erlanga.m.
 *
 * <p>Reference: W. Whitt (2005). Engineering solution of a basic call-center
 * model. Management Science 51(2), 221-235, Section 7 and eq. (7.12). The model
 * itself is due to C. Palm (1937, 1957).
 *
 * @since LINE 3.1.0
 */
public final class Qsys_erlanga {

    private Qsys_erlanga() {
    }

    /**
     * Erlang A with an infinite waiting room.
     *
     * @param lambda arrival rate
     * @param mu     service rate of one server
     * @param theta  abandonment rate of a waiting customer
     * @param s      number of servers
     * @return the steady-state measures
     */
    public static QsysAbandonResult qsys_erlanga(double lambda, double mu, double theta, int s) {
        return qsys_erlanga(lambda, mu, theta, s, Double.POSITIVE_INFINITY);
    }

    /**
     * Erlang A with r extra waiting spaces, so an arrival finding s+r customers
     * is blocked and lost.
     *
     * @param lambda arrival rate
     * @param mu     service rate of one server
     * @param theta  abandonment rate of a waiting customer
     * @param s      number of servers
     * @param r      extra waiting spaces, {@code Double.POSITIVE_INFINITY} if unbounded
     * @return the steady-state measures
     */
    public static QsysAbandonResult qsys_erlanga(double lambda, double mu, double theta, int s,
                                                 double r) {
        return qsys_erlanga(lambda, mu, theta, s, r, null);
    }

    /**
     * Erlang A, with the waiting-time cdfs when times are supplied.
     *
     * @param lambda  arrival rate
     * @param mu      service rate of one server
     * @param theta   abandonment rate of a waiting customer
     * @param s       number of servers
     * @param r       extra waiting spaces, {@code Double.POSITIVE_INFINITY} if unbounded
     * @param wPoints times at which to evaluate the waiting-time cdfs, or null
     * @return the steady-state measures
     */
    public static QsysAbandonResult qsys_erlanga(double lambda, double mu, double theta, int s,
                                                 double r, double[] wPoints) {
        if (theta <= 0 && Double.isInfinite(r) && lambda >= s * mu) {
            throw new RuntimeException("qsys_erlanga: without abandonment (theta = 0) and with an "
                    + "infinite waiting room the queue is unstable at lambda >= s*mu; give a "
                    + "finite r or a positive theta");
        }
        return Qsys_mgisrgi_whitt.qsys_mgisrgi_whitt(lambda, mu, s, r, Patience.exponential(theta),
                wPoints, Qsys_mgisrgi_whitt.DEFAULT_MAX_QUEUE, Qsys_mgisrgi_whitt.DEFAULT_TOL,
                "euler", 41);
    }
}
