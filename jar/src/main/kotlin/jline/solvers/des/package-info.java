/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Discrete Event Simulation (DES) solver using SSJ library.
 *
 * <p>This package contains the DES solver implementation that uses the
 * SSJ (Stochastic Simulation in Java) library for discrete-event simulation
 * of queueing networks.</p>
 *
 * <p>Key components:
 * <ul>
 *   <li>{@link jline.solvers.des.SolverDES} - Main solver class</li>
 *   <li>{@link jline.solvers.des.DESResult} - Result container</li>
 *   <li>{@link jline.solvers.des.DESOptions} - Configuration options</li>
 * </ul>
 * </p>
 *
 * <p>Currently supports:
 * <ul>
 *   <li>M/M/1 queues (single server FCFS)</li>
 *   <li>Multiclass workloads</li>
 *   <li>Open queueing networks</li>
 * </ul>
 * </p>
 *
 * @see jline.solvers.des.SolverDES
 * @since 1.0
 */
package jline.solvers.des;
