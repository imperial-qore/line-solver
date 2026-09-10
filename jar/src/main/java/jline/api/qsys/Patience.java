/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.function.DoubleUnaryOperator;

/**
 * The patience (time-to-abandon) law of a queue with customer abandonment,
 * in the three forms accepted by {@link Qsys_mgisrgi_whitt}.
 *
 * <p>The engineering solution of W. Whitt (2005) needs the patience law only
 * through its HAZARD RATE near the origin, so the three factories below carry
 * exactly that: a constant hazard (exponential patience, for which the analysis
 * becomes exact), an explicit hazard function, or the complementary cdf from
 * which the hazard is recovered by integration.
 *
 * @since LINE 3.1.0
 */
public final class Patience {

    /** How the abandonment rates are read off this law. */
    public enum Form {
        /** Constant hazard theta: exponential patience, the Erlang A case. */
        EXPONENTIAL,
        /** Explicit hazard rate h(t), used pointwise as in eq. (3.3). */
        HAZARD,
        /** Complementary cdf G(t) = 1-F(t), integrated as in eq. (3.6). */
        CCDF
    }

    private final Form form;
    private final double theta;
    private final DoubleUnaryOperator fun;

    private Patience(Form form, double theta, DoubleUnaryOperator fun) {
        this.form = form;
        this.theta = theta;
        this.fun = fun;
    }

    /**
     * Exponential patience of rate theta, so h(t) = theta for every t.
     *
     * @param theta abandonment rate of a waiting customer, non-negative
     * @return the patience law
     */
    public static Patience exponential(double theta) {
        if (theta < 0) {
            throw new RuntimeException("Patience.exponential: the rate theta must be non-negative");
        }
        return new Patience(Form.EXPONENTIAL, theta, null);
    }

    /**
     * Patience given by its hazard rate h = f/(1-F).
     *
     * @param hazard the hazard rate function
     * @return the patience law
     */
    public static Patience hazard(DoubleUnaryOperator hazard) {
        if (hazard == null) {
            throw new RuntimeException("Patience.hazard: the hazard function must not be null");
        }
        return new Patience(Form.HAZARD, Double.NaN, hazard);
    }

    /**
     * Patience given by its complementary cdf G(t) = 1-F(t).
     *
     * @param ccdf the complementary cdf, positive on the range of interest
     * @return the patience law
     */
    public static Patience ccdf(DoubleUnaryOperator ccdf) {
        if (ccdf == null) {
            throw new RuntimeException("Patience.ccdf: the ccdf must not be null");
        }
        return new Patience(Form.CCDF, Double.NaN, ccdf);
    }

    /** @return which of the three forms this law carries. */
    public Form getForm() {
        return form;
    }

    /** @return the exponential rate, NaN unless the form is EXPONENTIAL. */
    public double getTheta() {
        return theta;
    }

    /** @return the hazard rate at t; only meaningful for EXPONENTIAL and HAZARD. */
    public double hazardAt(double t) {
        if (form == Form.EXPONENTIAL) {
            return theta;
        }
        if (form == Form.HAZARD) {
            return fun.applyAsDouble(t);
        }
        throw new RuntimeException("Patience.hazardAt: this law is given by its ccdf");
    }

    /** @return the complementary cdf at t; only meaningful for CCDF. */
    public double ccdfAt(double t) {
        if (form != Form.CCDF) {
            throw new RuntimeException("Patience.ccdfAt: this law is not given by its ccdf");
        }
        return fun.applyAsDouble(t);
    }
}
