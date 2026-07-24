/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Age of Information (AoI) analysis algorithms.
 *
 * <p>This package provides analytical formulas and Laplace-Stieltjes Transform
 * (LST) functions for computing Age of Information metrics in single-server
 * queueing systems.</p>
 *
 * <p><b>Scheduling Disciplines:</b>
 * <ul>
 *   <li>FCFS: First-Come First-Served</li>
 *   <li>LCFS-PR: Last-Come First-Served with Preemption</li>
 *   <li>LCFS-S: LCFS with Service discarding (set-aside)</li>
 *   <li>LCFS-D: LCFS with Departure discarding</li>
 * </ul>
 * </p>
 *
 * <p><b>Supported Queue Types:</b>
 * <ul>
 *   <li>M/M/1, M/D/1, D/M/1 (closed-form results)</li>
 *   <li>M/GI/1, GI/M/1 (semi-analytical with LST)</li>
 * </ul>
 * </p>
 *
 * <p><b>Reference:</b>
 * Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
 * the Stationary Distribution of the Age of Information and Its
 * Application to Single-Server Queues," IEEE Trans. Information Theory,
 * vol. 65, no. 12, pp. 8305-8324, 2019.</p>
 *
 * @since LINE 3.0
 */
package jline.api.aoi;
