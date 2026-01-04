/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Comprehensive examples and tutorials for LINE queueing network modeling.
 *
 * <p>This package serves as the main entry point for LINE examples, providing
 * a comprehensive gallery of queueing network models from basic single-queue
 * systems to complex multi-class networks with advanced features.
 *
 * <h2>Sub-packages</h2>
 * <ul>
 * <li><strong>{@link jline.examples.java}</strong> - Java implementation examples</li>
 * <li><strong>kotlin</strong> - Kotlin implementation examples (development versions)</li>
 * </ul>
 *
 * <h2>Example Categories</h2>
 * <ul>
 * <li><strong>Getting Started</strong> - Basic tutorials and first examples</li>
 * <li><strong>Classical Systems</strong> - M/M/1, M/M/k, M/G/1 queues</li>
 * <li><strong>Queueing Networks</strong> - Open, closed, and mixed networks</li>
 * <li><strong>Advanced Models</strong> - Fork-join, cache models, layered systems</li>
 * <li><strong>Specialized Features</strong> - Priority queues, class switching, state-dependent routing</li>
 * </ul>
 *
 * <h2>Model Complexity Levels</h2>
 * <ul>
 * <li><strong>Basic</strong> - Single-class networks with standard features</li>
 * <li><strong>Intermediate</strong> - Multi-class networks with routing</li>
 * <li><strong>Advanced</strong> - Complex topologies and specialized disciplines</li>
 * <li><strong>Research</strong> - Cutting-edge features and experimental models</li>
 * </ul>
 *
 * <h2>Usage Patterns</h2>
 * <p>Examples typically follow this structure:
 * <pre>{@code
 * public class ExampleClass {
 *     public static Network createModel() {
 *         // Model construction
 *         Network model = new Network("ExampleModel");
 *         // ... add nodes, classes, parameters
 *         return model;
 *     }
 *     
 *     public static void main(String[] args) {
 *         Network model = createModel();
 *         NetworkAvgTable results = new MVA(model).getAvgTable();
 *         results.print();
 *     }
 * }
 * }</pre>
 *
 * @see jline.lang.Network
 * @see jline.solvers
 * @since LINE 2.0
 */
package jline.examples;