/**
 * Package containing trace fitting utility functions.
 * Ported from BUTools V2.0.
 *
 * Functions include:
 * - SquaredDifference: Compute squared difference between vectors
 * - RelativeEntropy: Compute Kullback-Leibler divergence
 * - LikelihoodFromTracePH: Evaluate log-likelihood with PH distribution
 * - LikelihoodFromTraceMAP: Evaluate log-likelihood with MAP
 *
 * Note: Full trace fitting algorithms (PHFromTrace, MAPFromTrace) use
 * the EMpht library which is in jline.lib.empht package.
 */
package jline.lib.butools.fitting;
