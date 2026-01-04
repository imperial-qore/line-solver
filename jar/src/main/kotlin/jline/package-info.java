/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * LINE - Queueing theory algorithms.
 *
 * <p>LINE is an open-source package for analyzing queueing models via analytical methods
 * and simulation. It provides a complete modeling environment supporting open queueing
 * systems, closed and open queueing networks, and layered queueing networks.
 *
 * <h2>Core Packages</h2>
 * <ul>
 * <li><strong>{@link jline.lang}</strong> - Core modeling language and network definitions</li>
 * <li><strong>{@link jline.solvers}</strong> - Solution algorithms (MVA, JMT, SSA, CTMC, etc.)</li>
 * <li><strong>{@link jline.api}</strong> - Low-level algorithmic implementations</li>
 * <li><strong>{@link jline.examples}</strong> - Comprehensive example models and tutorials</li>
 * <li><strong>{@link jline.util}</strong> - Utility classes and helper functions</li>
 * <li><strong>{@link jline.lib}</strong> - External library dependencies</li>
 * <li><strong>{@link jline.bench}</strong> - Performance benchmarking utilities</li>
 * <li><strong>{@link jline.io}</strong> - Input/output and model serialization</li>
 * <li><strong>{@link jline.cli}</strong> - Command-line interface</li>
 * </ul>
 *
 * <h2>Supported Model Types</h2>
 * <ul>
 * <li><strong>Open Networks</strong>: Jobs arrive from external sources (e.g., Poisson arrivals)</li>
 * <li><strong>Closed Networks</strong>: Fixed population circulating within the system</li>
 * <li><strong>Mixed Networks</strong>: Combination of open and closed job classes</li>
 * <li><strong>Layered Networks</strong>: Hierarchical models with tasks and activities</li>
 * <li><strong>Petri Nets</strong>: Token-based models with places and transitions</li>
 * </ul>
 *
 * <h2>Solution Methods</h2>
 * <ul>
 * <li><strong>Exact Methods</strong>: MVA, Convolution, Matrix Analytic Methods</li>
 * <li><strong>Simulation</strong>: Discrete-event simulation via JMT integration</li>
 * <li><strong>Numerical</strong>: CTMC, Fluid analysis, SSA</li>
 * <li><strong>Approximation</strong>: Heavy-traffic, diffusion approximations</li>
 * </ul>
 *
 * <h2>Key Features</h2>
 * <ul>
 * <li>Support for complex routing patterns and class switching</li>
 * <li>Multiple scheduling disciplines (FCFS, PS, LCFS, HOL, DPS, etc.)</li>
 * <li>Various service time distributions (Exponential, Erlang, HyperExp, APH, MAP)</li>
 * <li>Load-dependent services and finite capacity queues</li>
 * <li>Fork-join synchronization and cache modeling</li>
 * <li>Random environments</li>
 * </ul>
 *
 * <p><strong>Getting Started:</strong>
 * <pre>{@code
 * // Create a simple M/M/1 queue
 * Network model = new Network("MM1");
 * Source source = new Source(model, "Source");
 * Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
 * Sink sink = new Sink(model, "Sink");
 * 
 * // Set parameters
 * OpenClass openClass = new OpenClass(model, "Class1", 0);
 * source.setArrival(openClass, Exp.fitMean(1.0));
 * queue.setService(openClass, Exp.fitMean(2.0));
 * 
 * // Connect nodes using routing matrix
 * RoutingMatrix routingMatrix = model.initRoutingMatrix();
 * routingMatrix.set(openClass, openClass, source, queue, 1.0);
 * routingMatrix.set(openClass, openClass, queue, sink, 1.0);
 * model.link(routingMatrix);
 * 
 * // Solve
 * NetworkAvgTable avgTable = new SolverJMT(model).getAvgTable();
 * avgTable.print();
 * }</pre>
 *
 * @see jline.lang.Network
 * @see jline.solvers
 * @see jline.examples
 * @since LINE 2.0
 * @version 3.0
 */
package jline;