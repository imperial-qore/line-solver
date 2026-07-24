/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * LINE Discrete Event Simulator (LDES) solver using SSJ library.
 *
 * <p>This package contains the LDES solver implementation that uses the
 * SSJ (Stochastic Simulation in Java) library for discrete-event simulation
 * of queueing networks.</p>
 *
 * <p>Key components:
 * <ul>
 *   <li>{@link jline.solvers.ldes.SolverLDES} - Main solver class</li>
 *   <li>{@link jline.solvers.ldes.LDESResult} - Result container</li>
 *   <li>{@link jline.solvers.ldes.LDESOptions} - Configuration options</li>
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
 * @see jline.solvers.ldes.SolverLDES
 * @since 1.0
 */
package jline.solvers.ldes;
