/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Fork-Join model benchmarks for parallel processing system evaluation.
 *
 * <p>This package provides benchmark models and test cases for fork-join
 * queueing systems, which are fundamental in modeling parallel processing
 * architectures, distributed systems, and multi-threaded applications.
 *
 * <h2>Fork-Join Architectures</h2>
 * <ul>
 * <li>Synchronous fork-join (all tasks must complete)</li>
 * <li>Asynchronous fork-join (subset completion)</li>
 * <li>Nested fork-join hierarchies</li>
 * <li>Load balancing in parallel queues</li>
 * <li>Barrier synchronization models</li>
 * </ul>
 *
 * <h2>Benchmark Scenarios</h2>
 * <ul>
 * <li>Varying degrees of parallelism</li>
 * <li>Heterogeneous server capabilities</li>
 * <li>Different task size distributions</li>
 * <li>Network communication delays</li>
 * <li>Failure and recovery modeling</li>
 * </ul>
 *
 * <h2>Applications</h2>
 * <ul>
 * <li>MapReduce and distributed computing frameworks</li>
 * <li>Multi-core processor performance analysis</li>
 * <li>Parallel database query processing</li>
 * <li>Scientific computing workflows</li>
 * </ul>
 *
 * @see jline.bench
 * @see jline.lang.nodes.Join
 * @see jline.lang.nodes.Fork
 * @since LINE 2.0
 */
package jline.bench.fj;