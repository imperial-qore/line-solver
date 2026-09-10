/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Numerical inversion of Laplace transforms.
 *
 * <p>This package presents the same surface as the native Python
 * {@code api/lti}, the MATLAB {@code matlab/src/api/lti} and the C++
 * {@code line/api/lti}: Euler, Talbot, Gaver-Stehfest, CME and Weeks. The
 * coefficient generators of the first three, and the CME tables, already live
 * in {@code jline.lib.lti}; this package delegates to them rather than
 * duplicating them, and adds the Weeks (Laguerre) method, which had no
 * implementation in any codebase.
 *
 * @since LINE 3.0
 */
package jline.api.lti;
