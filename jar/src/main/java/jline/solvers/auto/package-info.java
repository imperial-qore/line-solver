/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Automatic solver selection and configuration for queueing networks.
 *
 * <p>This package provides an implementation of SolverAUTO.
 * It implements intelligent solver selection algorithms that
 * automatically choose the most appropriate solution method based on
 * network characteristics, model size, and accuracy requirements.
 *
 * <h2>Selection Criteria</h2>
 * <ul>
 * <li>Network topology and structure analysis</li>
 * <li>Problem size and complexity assessment</li>
 * <li>Required accuracy and performance trade-offs</li>
 * <li>Available computational resources</li>
 * <li>Model-specific feature requirements</li>
 * </ul>
 *
 * <h2>Supported Solvers</h2>
 * <ul>
 * <li><strong>MVA</strong> - For Product Form closed networks</li>
 * <li><strong>JMT</strong> - For complex networks requiring simulation</li>
 * <li><strong>CTMC</strong> - For small state-space models</li>
 * <li><strong>SSA</strong> - For steady-state analysis</li>
 * <li><strong>Fluid</strong> - For large population models</li>
 * </ul>
 *
 * <h2>Decision Algorithms</h2>
 * <ul>
 * <li>Rule-based selection for common patterns</li>
 * <li>Performance prediction models</li>
 * <li>Machine learning-based recommendations</li>
 * <li>Hybrid solver cascading strategies</li>
 * </ul>
 *
 * <h2>Configuration Features</h2>
 * <ul>
 * <li>Automatic parameter tuning</li>
 * <li>Convergence criteria optimization</li>
 * <li>Resource allocation management</li>
 * <li>Fallback strategy implementation</li>
 * </ul>
 *
 * @see jline.solvers
 * @since LINE 2.0
 */
package jline.solvers.auto;