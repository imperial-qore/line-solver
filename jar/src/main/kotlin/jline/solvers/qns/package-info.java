/**
 * This package provides an implementation of SolverQNS.
 * Wrapper for the QNS utility part of the LQNS package.
 * 
 * This package provides a Java/Kotlin implementation of the QNS solver,
 * which is a wrapper around the external qnsolver command-line tool.
 * The solver supports various multiserver approximation methods for
 * analyzing queueing networks.
 * 
 * Key classes:
 * - SolverQNS: Main solver class extending NetworkSolver
 * - QNSResult: Result container for QNS solver output
 * - Solver_qns_analyzer: Analyzer implementation in Kotlin
 * - Solver_qns: Core solver handler in Kotlin
 * 
 * Supported multiserver approximation methods:
 * - conway: Conway's approximation
 * - rolia: Rolia's approximation  
 * - zhou: Zhou's approximation
 * - suri: Suri's approximation
 * - reiser: Reiser's approximation
 * - schmidt: Schmidt's approximation
 * 
 * The solver requires the external qnsolver executable to be available
 * in the system PATH.
 */
package jline.solvers.qns;