/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Rate-scaled copies of a distribution, preserving its shape.
 *
 * <p>Java counterpart of MATLAB's {@code dist_scale_rate.m}.</p>
 *
 * @see jline.solvers.NetworkSolver#getSensitivityTable(String, double, String)
 */
public class DistributionScaling {

    private DistributionScaling() {
    }

    /**
     * Returns a new distribution whose rate is {@code factor} times the rate of
     * {@code distrib}, i.e. the time-scaled variable X/factor. The scaling is exact:
     * every moment of order n is divided by factor^n, so the mean is divided by
     * factor while the SCV, the skewness and the whole shape of the distribution are
     * preserved.
     *
     * <p>The scaled object is rebuilt from the parameters of the original, rather
     * than by rescaling its Markovian representation, so that the parameter list
     * stays coherent with the distribution family. Solvers that serialize the model
     * (JMT, LDES) read the parameters, and would otherwise export the unscaled
     * process.</p>
     *
     * <p>This is the perturbation primitive of the finite-difference branch of
     * {@code getSensitivityTable}: scaling the rate at a station-class by (1+h) is
     * exactly the perturbation the derivative d(.)/d(rate) is taken along.</p>
     *
     * @param distrib a distribution
     * @param factor  positive scaling factor for the rate
     * @return a new distribution of the same family with rate*factor
     */
    @SuppressWarnings("unchecked")
    public static Distribution scaleRate(Distribution distrib, double factor) {
        if (distrib == null) {
            line_error(mfilename(new Object() {
            }), "dist_scale_rate expects a Distribution object.");
        }
        if (Double.isNaN(factor) || Double.isInfinite(factor) || factor <= 0) {
            line_error(mfilename(new Object() {
            }), "The scaling factor must be a positive finite scalar.");
        }
        Class<?> cls = distrib.getClass();
        if (cls == Exp.class) {
            return new Exp(((Double) distrib.getParam(1).getValue()).doubleValue() * factor);
        } else if (cls == Erlang.class) {
            return new Erlang(((Double) distrib.getParam(1).getValue()).doubleValue() * factor,
                    ((Number) distrib.getParam(2).getValue()).intValue());
        } else if (cls == HyperExp.class) {
            // The n-phase form carries the whole rate vector; every rate is scaled.
            HyperExp he = (HyperExp) distrib;
            double[] lambda = he.getLambda();
            for (int i = 0; i < lambda.length; i++) {
                lambda[i] = lambda[i] * factor;
            }
            return new HyperExp(he.getP(), lambda);
        } else if (cls == Cox2.class) {
            // Completion probabilities are dimensionless and are left untouched;
            // only the phase rates carry the time scale.
            Coxian cx = (Coxian) distrib;
            Matrix mu = cx.getMu();
            Matrix phi = cx.getPhi();
            return new Cox2(mu.get(0) * factor, mu.get(1) * factor, phi.get(0));
        } else if (cls == Coxian.class) {
            Coxian cx = (Coxian) distrib;
            return new Coxian(cx.getMu().scale(factor), cx.getPhi());
        } else if (cls == APH.class) {
            return new APH((Matrix) distrib.getParam(2).getValue(),
                    ((Matrix) distrib.getParam(3).getValue()).scale(factor));
        } else if (cls == PH.class) {
            return new PH((Matrix) distrib.getParam(2).getValue(),
                    ((Matrix) distrib.getParam(3).getValue()).scale(factor));
        } else if (cls == MAP.class) {
            return new MAP(((Matrix) distrib.getParam(1).getValue()).scale(factor),
                    ((Matrix) distrib.getParam(2).getValue()).scale(factor));
        } else if (cls == MMPP2.class) {
            // Every rate of the modulating chain and of the arrival process is
            // scaled, which time-scales the whole process.
            return new MMPP2(((Double) distrib.getParam(1).getValue()).doubleValue() * factor,
                    ((Double) distrib.getParam(2).getValue()).doubleValue() * factor,
                    ((Double) distrib.getParam(3).getValue()).doubleValue() * factor,
                    ((Double) distrib.getParam(4).getValue()).doubleValue() * factor);
        } else if (cls == Det.class) {
            return new Det(((Double) distrib.getParam(1).getValue()).doubleValue() / factor);
        } else if (cls == Uniform.class) {
            return new Uniform(((Double) distrib.getParam(1).getValue()).doubleValue() / factor,
                    ((Double) distrib.getParam(2).getValue()).doubleValue() / factor);
        } else if (cls == Gamma.class) {
            // Gamma(shape, scale): the shape is dimensionless.
            return new Gamma(((Double) distrib.getParam(1).getValue()).doubleValue(),
                    ((Double) distrib.getParam(2).getValue()).doubleValue() / factor);
        } else if (cls == Pareto.class) {
            // Pareto(shape, scale): the scale is the minimum of the support.
            return new Pareto(((Double) distrib.getParam(1).getValue()).doubleValue(),
                    ((Double) distrib.getParam(2).getValue()).doubleValue() / factor);
        } else if (cls == Weibull.class) {
            // Weibull params are (1) the scale alpha and (2) the shape r, while the
            // constructor takes (shape, scale).
            return new Weibull(((Double) distrib.getParam(2).getValue()).doubleValue(),
                    ((Double) distrib.getParam(1).getValue()).doubleValue() / factor);
        } else if (cls == Lognormal.class) {
            // X/factor is lognormal with mu - log(factor) and the same sigma.
            return new Lognormal(((Double) distrib.getParam(1).getValue()).doubleValue() - Math.log(factor),
                    ((Double) distrib.getParam(2).getValue()).doubleValue());
        } else if (cls == NHPP.class) {
            // see _kb/01-model-classes.md (Java process-construction notes) for rationale
            NHPP nhpp = (NHPP) distrib;
            double[] breakpoints = nhpp.getBreakpoints();
            double[] rates = nhpp.getRates();
            double[] scaledBreakpoints = new double[breakpoints.length];
            for (int i = 0; i < breakpoints.length; i++) {
                scaledBreakpoints[i] = breakpoints[i] / factor;
            }
            double[] scaledRates = new double[rates.length];
            for (int i = 0; i < rates.length; i++) {
                scaledRates[i] = rates[i] * factor;
            }
            return new NHPP(scaledBreakpoints, scaledRates, nhpp.isCyclic());
        } else if (cls == Replayer.class) {
            // see _kb/01-model-classes.md (Java process-construction notes) for rationale
            Replayer replayer = (Replayer) distrib;
            if (replayer.getFileName() != null) {
                line_error(mfilename(new Object() {
                }), "Rate scaling is not defined for a file-backed Replayer: the trace "
                        + "is the parameter. Pass the samples as an array, "
                        + "new Replayer(double[]), to differentiate it.");
            }
            double[] data = replayer.getData();
            double[] scaledData = new double[data.length];
            for (int i = 0; i < data.length; i++) {
                scaledData[i] = data[i] / factor;
            }
            return new Replayer(scaledData);
        } else if (cls == Immediate.class) {
            return new Immediate();
        }
        line_error(mfilename(new Object() {
        }), "Rate scaling is not defined for a " + cls.getSimpleName() + " process. "
                + "Supported: Exp, Erlang, HyperExp, Coxian, Cox2, APH, PH, MAP, MMPP2, "
                + "Det, Uniform, Gamma, Pareto, Weibull, Lognormal, NHPP, Replayer, Immediate.");
        return null;
    }
}
