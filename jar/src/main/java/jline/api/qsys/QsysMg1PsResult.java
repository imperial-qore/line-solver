package jline.api.qsys;

import org.apache.commons.math3.complex.Complex;

/**
 * Result of {@link Qsys_mg1_ps#qsys_mg1_ps}: the sojourn time distribution of the
 * M/G/1 processor-sharing queue.
 *
 * <p>Besides the tabulated arrays, the object stays callable: {@link #lstCond},
 * {@link #lstExcess}, {@link #lstUncond} and {@link #dominantRoot} evaluate the
 * transforms at any argument, which is what the MATLAB twin returns as function
 * handles.</p>
 *
 * <p>Port of MATLAB qsys_mg1_ps.m.</p>
 */
public class QsysMg1PsResult {

    /** Utilization lambda*m1. */
    public double rho;
    /** Mean service requirement. */
    public double m1;
    /** Second moment of the service requirement, NaN when only an LST is given. */
    public double m2;
    /** Exact unconditional mean sojourn time m1/(1-rho). */
    public double meanUncond;
    /** Unconditional second moment of the sojourn time. */
    public double m2Uncond = Double.NaN;
    /** Unconditional variance of the sojourn time. */
    public double varUncond = Double.NaN;
    /** Mass (1-rho)*bhat(lambda) of the jobs that share the server with nobody. */
    public double atomUncond;

    /** Requested service requirements. */
    public double[] x;
    /** Requested transform arguments. */
    public double[] s;
    /** Requested times. */
    public double[] t;

    /** Atom (1-rho)*exp(-lambda*x) of V(x) at t = x. */
    public double[] atomCond;
    /** Exact conditional mean x/(1-rho). */
    public double[] meanCond;
    /** Conditional second moment of V(x). */
    public double[] m2Cond;
    /** Conditional variance of V(x). */
    public double[] varCond;
    /** Values of the conditional LST, indexed [x][s]. */
    public double[][] lstCondVal;
    /** Values of the unconditional LST. */
    public double[] lstUncondVal;
    /** Density of V(x), indexed [x][t], NaN on the lattice (k+1)*x. */
    public double[][] pdfCond;
    /** P(V(x) &lt;= t), indexed [x][t]. */
    public double[][] cdfCond;
    /** Density of the unconditional sojourn time. */
    public double[] pdfUncond;
    /** P(V &lt;= t). */
    public double[] cdfUncond;

    // state carried so that the transforms stay callable after the call returns
    double lambda;
    boolean isPH;
    double[] db;
    double[] nb;
    Qsys_mg1_ps.ServiceLST bhat;
    Qsys_mg1_ps.ServicePDF bpdf;
    double shift;
    int nterms;
    double[] yq;
    double[] wq;
    double[] bq;

    /** Conditional transform E[exp(-s V(x))]. */
    public Complex lstCond(Complex s, double x) {
        return Qsys_mg1_ps.lstCond(this, s, x, false);
    }

    /** Transform E[exp(-s (V(x)-x))] of the excess, bounded as s grows. */
    public Complex lstExcess(Complex s, double x) {
        return Qsys_mg1_ps.lstCond(this, s, x, true);
    }

    /** Unconditional transform E[exp(-s V)]. */
    public Complex lstUncond(Complex s) {
        return Qsys_mg1_ps.lstUncond(this, s);
    }

    /** Dominant singularity tau*(s), the right half plane root of tau = s + lambda*(1-bhat(tau)). */
    public Complex dominantRoot(Complex s) {
        return Qsys_mg1_ps.dominantRoot(this, s);
    }
}
