/**
 * Stochastic network calculus (SNC).
 *
 * <p>The only api domain in LINE whose deliverable is a TAIL QUANTILE with a
 * certified violation probability rather than a mean. An envelope is the pair
 * {@code (sigma(theta), rho(theta))} in the exponential form</p>
 *
 * <pre>
 *   E[exp(theta*A(s,t))] &lt;= exp(theta*(rho*(t-s) + sigma)),   theta &gt; 0,
 * </pre>
 *
 * <p>and its service counterpart with the sign of theta reversed. The min-plus
 * operations ({@link jline.api.snc.Snc_leftover}, {@link jline.api.snc.Snc_conv},
 * {@link jline.api.snc.Snc_output}) compose envelopes across a feed-forward
 * network, and every bound is a Chernoff infimum over theta computed by
 * {@link jline.api.snc.Snc_thetaopt}.</p>
 *
 * <p>TIME IS SLOTTED with unit slot length, which is what makes the geometric
 * sum over the start of the backlogged period converge to
 * {@code 1/(1-exp(-theta*(rhoS-rhoA)))}; the continuous-time formulation would
 * give {@code 1/(theta*(rhoS-rhoA))} instead. INFEASIBILITY IS SIGNALLED BY
 * {@code Infinity}, never by an exception, so the theta search can discard the
 * point.</p>
 *
 * <p>Port of matlab/src/api/snc. Wired into {@link jline.solvers.ba.SolverBA} as
 * the method {@code snc.upper}.</p>
 *
 * <p>Reference: M. Fidler, A. Rizk, "A Guide to the Stochastic Network
 * Calculus", IEEE Communications Surveys and Tutorials 17(1), 92-105, 2015.</p>
 */
package jline.api.snc;
