/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Closed Queueing Network (CQN) benchmarks for solver performance evaluation.
 *
 * <p>This package contains benchmark models and utilities specifically designed
 * for testing the performance and accuracy of solvers on closed queueing networks.
 * These benchmarks cover various network topologies, service disciplines, and
 * parameter ranges commonly encountered in practice.
 *
 * <h2>Benchmark Models</h2>
 * <ul>
 * <li>Single-class and multi-class closed networks</li>
 * <li>Load-dependent and load-independent services</li>
 * <li>Various routing topologies (cyclic, tree, mesh)</li>
 * <li>Different service time distributions</li>
 * <li>Networks with blocking and non-blocking queues</li>
 * </ul>
 *
 * <h2>Performance Metrics</h2>
 * <ul>
 * <li>Solution time for different solver algorithms</li>
 * <li>Memory usage during computation</li>
 * <li>Numerical accuracy compared to reference solutions</li>
 * <li>Scalability with respect to network size and population</li>
 * </ul>
 *
 * @see jline.bench
 * @see jline.solvers.mva
 * @see jline.solvers.nc
 * @since LINE 2.0
 */
package jline.bench.cqn;