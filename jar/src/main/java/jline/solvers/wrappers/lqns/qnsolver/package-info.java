/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The qnsolver engine behind the qns methods of
 * {@link jline.solvers.wrappers.lqns.SolverLQNS} on a flat Network.
 *
 * <p>qnsolver is the product-form MVA tool of the LQNS distribution. The model
 * is written as a JMVA document, qnsolver is run on it, and its results are read
 * back. A closed non-product-form Network does not come here: SolverLQNS
 * converts it with QN2LQN and solves it with lqns.
 *
 * Key classes:
 * - Solver_qns_analyzer: analyzer, called from SolverLQNS
 * - Solver_qns: core handler that writes the JMVA document and runs qnsolver
 *
 * Multiserver approximations, selected by the method qns.NAME:
 * - conway: Conway (1989), extending the multinomial all-servers-busy probability of de Souza e Silva and Muntz (Perform. Eval. 7(3), 1987)
 * - rolia: Rolia (PhD thesis, Toronto, 1992) as used in the method of layers (Rolia and Sevcik, IEEE TSE 21(8), 1995), in the per-class Rolia-Franks form of Franks (PhD thesis, Carleton, 1999)
 * - zhou: arrival-theorem binomial (AB) approximation, S. Zhou (M.A.Sc. thesis, Carleton, 2021) and Zhou and Woodside (ICPE Companion 2022)
 * - suri: Suri, Sahu and Vernon (IERC 2007)
 * - reiser: Reiser and Lavenberg (J. ACM 27(2), 1980) load-dependent MVA, see also Reiser (Perform. Eval. 1, 1981)
 * - schmidt: Schmidt (Perform. Eval. 29(4), 1997)
 *
 * The qnsolver executable must be available in the system PATH.
 */
package jline.solvers.wrappers.lqns.qnsolver;
