/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Benchmarking utilities and performance test models for LINE solvers.
 *
 * <p>This package provides a comprehensive benchmarking framework for evaluating
 * the performance and accuracy of LINE's various solvers across different model
 * types and problem sizes. It includes predefined benchmark models, performance
 * measurement utilities, and comparison tools.
 *
 * <h2>Benchmark Categories</h2>
 * <ul>
 * <li><strong>{@link jline.bench.cqn}</strong> - Closed Queueing Network benchmarks</li>
 * <li><strong>{@link jline.bench.fj}</strong> - Fork-Join model benchmarks</li>
 * <li><strong>{@link jline.bench.lqn}</strong> - Layered Queueing Network benchmarks</li>
 * </ul>
 *
 * <h2>Features</h2>
 * <ul>
 * <li>Standardized benchmark models for solver comparison</li>
 * <li>Performance timing and memory usage measurement</li>
 * <li>Accuracy validation against reference solutions</li>
 * <li>Scalability testing with varying model sizes</li>
 * <li>Statistical analysis of solver performance</li>
 * </ul>
 *
 * <h2>Usage</h2>
 * <p>Benchmark classes typically provide static methods to generate test models
 * of various sizes and complexities. These can be used to evaluate solver
 * performance across different scenarios.
 *
 * <p>Example usage:
 * <pre>{@code
 * // Run CQN benchmark
 * Network model = CQNBench.generate(10, 3); // 10 stations, 3 classes
 * long startTime = System.nanoTime();
 * NetworkAvgTable result = new SolverMVA(model).getAvgTable();
 * long duration = System.nanoTime() - startTime;
 * }</pre>
 *
 * @see jline.solvers
 * @since LINE 2.0
 */
package jline.bench;