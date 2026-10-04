/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

/**
 * Outcome of {@link Infer_variational}.
 */
public class VariationalResult {
    /** (P) posterior Gamma shapes. */
    public double[] alpha;
    /** (P) posterior Gamma rates. */
    public double[] beta;
    /** (P) posterior mean rates, alpha/beta. */
    public double[] rates;
    /** (P) mean service times, beta/alpha. */
    public double[] meanServiceTime;
    /** evidence lower bound per iteration. */
    public double[] bound;
    /** (P x iter) posterior shape per iteration. */
    public double[][] alphaTrace;
    /** (P x iter) posterior rate per iteration. */
    public double[][] betaTrace;
    /** (T x G x ymax+1) transition-count marginals. */
    public double[][][] Y;
    /** (T x G x ymax+1) variational rates. */
    public double[][][] nu;
    /** (G) time grid. */
    public double[] tgrid;
    /** (G x M*R) expected queue lengths. */
    public double[][] qlen;
    public int iter;
    public boolean converged;
    /** largest probability mass sitting on the truncated top count. */
    public double tailmass;
}
