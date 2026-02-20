/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Comprehensive examples and gallery of queueing network models for the LINE solver.
 *
 * <p>This package provides a extensive collection of example models demonstrating various
 * queueing network concepts, modeling techniques, and solution methods. The examples range
 * from classical single-queue systems to complex multi-class networks with advanced features.
 *
 * <h2>Example Categories</h2>
 * <ul>
 * <li><strong>Basic Networks</strong>: {@link OpenModel}, {@link ClosedModel}, {@link MixedModel}</li>
 * <li><strong>Advanced Topologies</strong>: {@link ForkJoinModel}, {@link CacheModel}, {@link LayeredModel}</li>
 * <li><strong>Specialized Models</strong>: {@link StochPetriNet}, {@link RandomEnvironment}, {@link LoadDependent}</li>
 * <li><strong>Routing and Control</strong>: {@link ClassSwitching}, {@link StateDepRouting}, {@link Prio}</li>
 * <li><strong>Classical Systems</strong>: {@link Gallery} (M/M/1, M/M/k, M/G/1, etc.)</li>
 * <li><strong>Analysis Tools</strong>: {@link StateProbabilities}, {@link CDFRespT}, {@link InitState}</li>
 * <li><strong>Getting Started</strong>: {@link GettingStarted}, {@link Misc}</li>
 * </ul>
 *
 * <h2>Model Types</h2>
 * <ul>
 * <li><strong>Open Networks</strong>: Jobs arrive from external sources and depart to sinks</li>
 * <li><strong>Closed Networks</strong>: Fixed population of jobs circulating within the system</li>
 * <li><strong>Mixed Networks</strong>: Combination of open and closed classes</li>
 * <li><strong>Layered Networks</strong>: Hierarchical models with tasks, processors, and activities</li>
 * <li><strong>Petri Nets</strong>: Token-based models with places, transitions, and firing rules</li>
 * </ul>
 *
 * <h2>Key Features Demonstrated</h2>
 * <ul>
 * <li>Multiple scheduling strategies (FCFS, PS, LCFS, HOL, DPS, etc.)</li>
 * <li>Various service time distributions (Exponential, Erlang, HyperExp, APH, MAP, etc.)</li>
 * <li>Complex routing patterns (probabilistic, state-dependent, class switching)</li>
 * <li>Load balancing and priority mechanisms</li>
 * <li>Multi-server stations and load-dependent services</li>
 * <li>Cache modeling and replacement strategies</li>
 * <li>Fork-join synchronization patterns</li>
 * <li>Random environment and Markovian arrival processes</li>
 * </ul>
 *
 * <h2>Usage</h2>
 * <p>Each example class contains static methods that return configured {@link jline.lang.Network}
 * objects ready for analysis. These can be solved using various solvers from the
 * {@link jline.solvers} package such as MVA, JMT, SSA, CTMC, etc.
 *
 * <p>Example usage:
 * <pre>{@code
 * Network model = GettingStarted.tut01_mm1_basics();
 * NetworkAvgTable results = new MVA(model).getAvgTable();
 * results.print();
 * }</pre>
 *
 * @see jline.lang.Network
 * @see jline.solvers
 * @since LINE 2.0
 */
package jline.examples.java;