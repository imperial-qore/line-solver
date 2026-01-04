/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Matrix operations and linear algebra utilities for queueing analysis.
 *
 * <p>This package provides comprehensive matrix operations, linear algebra
 * algorithms, and numerical methods essential for solving queueing models
 * and stochastic processes. It includes specialized routines for sparse
 * matrices and structured linear systems.
 *
 * <h2>Matrix Operations</h2>
 * <ul>
 * <li>Basic arithmetic (addition, multiplication, transpose)</li>
 * <li>Matrix decompositions (LU, QR, SVD, Cholesky)</li>
 * <li>Eigenvalue and eigenvector computation</li>
 * <li>Matrix exponential calculation</li>
 * <li>Kronecker products and operations</li>
 * </ul>
 *
 * <h2>Linear System Solvers</h2>
 * <ul>
 * <li>Direct methods for dense systems</li>
 * <li>Iterative methods for sparse systems</li>
 * <li>Specialized solvers for structured matrices</li>
 * <li>Condition number estimation</li>
 * <li>Error analysis and numerical stability</li>
 * </ul>
 *
 * <h2>Specialized Algorithms</h2>
 * <ul>
 * <li>Quasi-Birth-Death matrix operations</li>
 * <li>Stochastic matrix analysis</li>
 * <li>Doubly stochastic matrix generation</li>
 * <li>Matrix polynomial evaluation</li>
 * <li>Uniformization matrix computations</li>
 * </ul>
 *
 * <h2>Performance Features</h2>
 * <ul>
 * <li>Sparse matrix support for large systems</li>
 * <li>Memory-efficient storage schemes</li>
 * <li>Parallel computation capabilities</li>
 * <li>Numerical optimization for stability</li>
 * </ul>
 *
 * @see jline.util
 * @since LINE 2.0
 */
package jline.util.matrix;