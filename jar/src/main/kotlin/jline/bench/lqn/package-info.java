/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Layered Queueing Network (LQN) benchmarks for hierarchical system evaluation.
 *
 * <p>This package contains benchmark models for layered queueing networks,
 * which are used to model multi-tier software architectures, distributed
 * systems, and service-oriented architectures where requests flow through
 * multiple layers of processing.
 *
 * <h2>LQN Model Types</h2>
 * <ul>
 * <li>Client-server architectures with multiple tiers</li>
 * <li>Web application models (presentation, business, data layers)</li>
 * <li>Microservices architectures</li>
 * <li>Database systems with caching layers</li>
 * <li>Cloud computing resource hierarchies</li>
 * </ul>
 *
 * <h2>Benchmark Features</h2>
 * <ul>
 * <li>Various layer interaction patterns</li>
 * <li>Different task and activity structures</li>
 * <li>Resource contention across layers</li>
 * <li>Synchronous and asynchronous communication</li>
 * <li>Load balancing and replication scenarios</li>
 * </ul>
 *
 * <h2>Performance Analysis</h2>
 * <ul>
 * <li>End-to-end response time analysis</li>
 * <li>Throughput bottleneck identification</li>
 * <li>Resource utilization across layers</li>
 * <li>Scalability analysis for tier expansion</li>
 * </ul>
 *
 * @see jline.bench
 * @see jline.lang.layered
 * @see jline.solvers.lqns
 * @since LINE 2.0
 */
package jline.bench.lqn;