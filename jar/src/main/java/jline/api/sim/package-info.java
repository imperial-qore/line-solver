/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Simulation output analysis.
 *
 * <p>This package turns a simulation sample path into an interval estimate. Its
 * subject is the statistics of the output process, not the queueing model that
 * produced it, so nothing here takes a {@code NetworkStruct}: the input is a
 * sequence of observations such as successive waiting times exported from a
 * simulation run.
 *
 * <p>The steady-state <em>mean</em> is already covered elsewhere: the LDES
 * engine forms batch means and reports confidence-interval half-widths for
 * every average metric, with MSER-5 warmup detection. What this package adds is
 * the steady-state <em>quantile</em>, which is the quantity a tail
 * service-level agreement is written against and for which no interval was
 * previously available anywhere in LINE.
 *
 * <p>The two entry points are {@link jline.api.sim.Sim_fquest} for a single
 * sample path and {@link jline.api.sim.Sim_firquest} for independent
 * replications. Both rest on the standardized time series of the batched
 * quantile process built by {@link jline.api.sim.Sim_sts_quantile_areas}, and on
 * two hypothesis tests, {@link jline.api.sim.Sim_vonneumann} for randomness and
 * {@link jline.api.sim.Sim_shapirowilk} for normality, which are general-purpose
 * and usable on their own.
 *
 * @since LINE 3.1.0
 */
package jline.api.sim;
