/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.io.Serializable;

import org.apache.commons.math3.analysis.UnivariateFunction;

import jline.util.matrix.Matrix;

/**
 * Service-time descriptor of the MAP/G/1/K family, the JAR form of the MATLAB
 * {@code svc} struct of {@code qsys_mapg1k.m} and of the C++
 * {@code line::qsys::ServiceLaw}.
 *
 * <p>The law is NOT fitted to a phase-type distribution before use: it enters
 * the analysis exactly, through the uniformization coefficients
 * c_n = E[e^{-theta S} (theta S)^n / n!]. That is why the descriptor keeps the
 * family rather than a moment vector.
 *
 * <p>Build one with the static factories; the fields a given kind does not use
 * are never read.
 */
public final class QsysServiceLaw implements Serializable {

    private static final long serialVersionUID = 1L;

    /** Which family the service law belongs to. */
    public enum Kind {
        /** Gamma(shape, scale); shape 1 is exponential, integer shape Erlang. */
        GAMMA,
        /** Constant service time. */
        DETERMINISTIC,
        /** Phase type (alpha, T). */
        PHASE_TYPE,
        /** Arbitrary density on (0, inf), or on (0, tmax] when tmax is given. */
        DENSITY
    }

    public final Kind kind;
    /** Gamma shape alpha. */
    public final double shape;
    /** Gamma scale theta. */
    public final double scale;
    /** Deterministic service time d. */
    public final double det;
    /** PH initial probability row, 1 x p. */
    public final Matrix phAlpha;
    /** PH subgenerator, p x p. */
    public final Matrix phT;
    /** Density of the service time. */
    public final UnivariateFunction pdf;
    /** Upper support limit of the density, read only when tmaxFinite. */
    public final double tmax;
    /** Whether the density has bounded support. */
    public final boolean tmaxFinite;

    private QsysServiceLaw(Kind kind, double shape, double scale, double det, Matrix phAlpha,
                           Matrix phT, UnivariateFunction pdf, double tmax, boolean tmaxFinite) {
        this.kind = kind;
        this.shape = shape;
        this.scale = scale;
        this.det = det;
        this.phAlpha = phAlpha;
        this.phT = phT;
        this.pdf = pdf;
        this.tmax = tmax;
        this.tmaxFinite = tmaxFinite;
    }

    /** Gamma(shape, scale). Mean is shape*scale, SCV is 1/shape. */
    public static QsysServiceLaw gamma(double shape, double scale) {
        return new QsysServiceLaw(Kind.GAMMA, shape, scale, 0.0, null, null, null, 0.0, false);
    }

    /** Exponential with the given RATE, i.e. Gamma(1, 1/rate). */
    public static QsysServiceLaw exponential(double rate) {
        return gamma(1.0, 1.0 / rate);
    }

    /** Constant service time d. */
    public static QsysServiceLaw deterministic(double d) {
        return new QsysServiceLaw(Kind.DETERMINISTIC, 0.0, 0.0, d, null, null, null, 0.0, false);
    }

    /** Phase type with initial row ALPHA and subgenerator T. */
    public static QsysServiceLaw phaseType(Matrix alpha, Matrix T) {
        return new QsysServiceLaw(Kind.PHASE_TYPE, 0.0, 0.0, 0.0, alpha, T, null, 0.0, false);
    }

    /** Arbitrary density on (0, inf). */
    public static QsysServiceLaw density(UnivariateFunction f) {
        return new QsysServiceLaw(Kind.DENSITY, 0.0, 0.0, 0.0, null, null, f, 0.0, false);
    }

    /** Arbitrary density supported on (0, tmax]. */
    public static QsysServiceLaw density(UnivariateFunction f, double tmax) {
        return new QsysServiceLaw(Kind.DENSITY, 0.0, 0.0, 0.0, null, null, f, tmax, true);
    }
}
